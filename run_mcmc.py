#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 2016

@author: cbarbosa
"""

import os
import shutil

import pymc
import numpy as np
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import fmin, fminbound
from scipy.interpolate import LinearNDInterpolator
from sklearn.mixture import GMM
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.backends.backend_pdf import PdfPages

import cap_mpfit as mpfit
from config import *

class SSP:
    """ Interface to interpolate and call different models. """
    def __init__(self, modelname, idx=None):
        self.modelname = modelname
        self.idx = idx if idx is not None else np.arange(25)
        self.load_model()
        self.model = LinearNDInterpolator(self.pars, self.data)
        self.calc_lims()

    def load_model(self):
        """ Load data for model and set appropriate indices. """
        if self.modelname in ["TMJ10ext", "TMJ10", "TMJ10Padova", "TMJ10ext"]:
            tablename = dict([("TMJ10ext", "tmj_metal_extrapolated_ews.dat"),
                              ("TMJ10", "tmj.dat"),
                              ("TMJ10Padova", "tmj_padova.dat"), ])
            self.table = os.path.join(tables_dir, tablename[self.modelname])
            self.pars = np.loadtxt(self.table, usecols=range(3))
            self.data = np.loadtxt(self.table, usecols=np.arange(25)+3)
        else:
            raise NotImplementedError("Model {0} is not available.".format(
                self.modelname))
        return

    def calc_lims(self):
        """ Get limits of model. """
        self.lims = np.array([self.pars.min(axis=0), self.pars.max(axis=0)]).T
        self.ranges = np.array([self.data.min(axis=0),
                                self.data.max(axis=0)]).T
        return

    def __call__(self, *args):
        return self.model(*args)[self.idx]

class Dist():
    """ Simple class to handle the distribution data of MCMC. """
    def __init__(self, data, lims):
        self.data = data
        self.lims = lims
        self.genextreme = genextreme(self.data)
        self.norm = norm(self.data)
        dists = [self.genextreme, self.norm]
        idx = np.argmin([x.ks for x in dists])
        self.best = dists[idx]
        return


class genextreme():
    def __init__(self, data):
        self.dist = stats.genextreme
        self.distname = "genextreme"
        self.data = data
        self.p = self.dist.fit(self.data)
        self.frozen = self.dist(self.p[0], loc=self.p[1], scale=self.p[2])
        self.pdf = lambda x : self.frozen.pdf(x)
        self.sample = self.frozen.rvs(len(self.data))
        self.sample2 = self.frozen.rvs(100000)
        self.moments = self.frozen.stats(moments="mvsk")
        self.MAPP = fmin(lambda x: -self.pdf(x),
                        self.moments[0], disp=0)[0]
        self.mean = self.frozen.mean()
        self.median = self.frozen.median()
        self.sigmin, self.sigmax = self.frozen.interval(0.68)
        self.lerr = self.MAPP - self.sigmin
        self.uerr = self.sigmax - self.MAPP
        self.ks, self.pvalue = stats.kstest(self.data, self.distname,
                                            args=self.p)

class norm():
    def __init__(self, data):
        self.dist = stats.norm
        self.distname = "norm"
        self.data = data
        self.p = self.dist.fit(self.data)
        self.moments = self.p
        self.frozen = self.dist(*self.p)
        self.pdf = lambda x : self.frozen.pdf(x)
        self.MAPP = fmin(lambda x: -self.pdf(x),
                        self.moments[0], disp=0)[0]
        self.mean = self.frozen.mean()
        self.median = self.frozen.median()
        self.sigmin, self.sigmax = self.frozen.interval(0.68)
        self.lerr = self.MAPP - self.sigmin
        self.uerr = self.sigmax - self.MAPP
        self.ks, self.pvalue = stats.kstest(self.data, self.distname,
                                            args=self.p)
        return

class statdist():
    def __init__(self, data, dist, distname):
        self.dist = dist
        self.distname = distname
        self.data = data
        self.p = self.dist.fit(self.data)
        self.pdf = lambda x : self.dist.pdf(x, *self.p[:-2], loc=self.p[-2],
                                            scale=self.p[-1])
        self.sample = stats.norm.rvs(self.p[0], size=len(self.data),
                                           scale=self.p[-1])
        self.moments = self.dist.stats(*self.p, moments="mvsk")
        self.MAPP = fmin(lambda x: -self.pdf(x),
                        self.moments[0], disp=0)[0]

class gmm():
    def __init__(self, data):
        self.distname = "gmm"
        self.data = data
        self.n_components = np.arange(1,11)
        self.models = []
        self.X = np.reshape(self.data, (len(self.data),1))
        for i in self.n_components:
            self.models.append(GMM(i, covariance_type='full').fit(self.X))
        self.AIC = np.array([m.aic(self.X) for m in self.models])
        self.BIC = np.array([m.bic(self.X) for m in self.models])
        self.k = 2 * np.arange(1,11)
        self.n = len(self.data)
        self.AICc = self.AIC + 2*self.k * (self.k + 1) / (self.n - self.k - 1)
        self.imin = np.minimum(np.argmin(self.AIC), np.argmin(self.BIC))
        self.best = self.models[self.imin]

def run_mcmc(lick, error, modelname, idx, dbname):
    """ Run the MCMC routine. """
    model = SSP(modelname, idx)
    ##########################################################################
    # Setting the priors
    age_dist = pymc.Uniform(name="age_dist", lower=model.lims[0,0],
                            upper=model.lims[0,1])
    metal_dist = pymc.Uniform(name="metal_dist", lower=model.lims[1,0],
                              upper=model.lims[1,1])
    alpha_dist = pymc.Uniform(name="alpha_dist", lower=model.lims[2,0],
                              upper=model.lims[2,1])
    ##########################################################################
    taus = 1 / error[model.idx]**2
    @pymc.deterministic()
    def ssp(age=age_dist, metal=metal_dist, alpha=alpha_dist):
        return model(age, metal, alpha)
    y = pymc.Normal(name="y", mu=ssp, tau=taus,
                    value=lick[idx], observed=True)
    M = pymc.Model([y, age_dist, metal_dist, alpha_dist])
    mcmc = pymc.MCMC(M, db="txt", dbname=dbname)
    mcmc.sample(20000, 1000, 4)
    mcmc.db.close()
    return


def run_candidates():
    """ Calculate stellar populations in candidates. """
    os.chdir(data_dir)
    modelname="TMJ10ext"
    bands = os.path.join(tables_dir, "bands.txt")
    ltype = np.loadtxt(bands, usecols=(8,))
    w = np.diff(np.loadtxt(bands, usecols=(4,5))).T[0]
    idxs = np.where(ltype==1)[0]
    filename = "results.tab"
    specs = np.loadtxt(filename, usecols=(0,), dtype=str)
    lick = np.loadtxt(filename, usecols=np.arange(13,62,2))
    error = np.loadtxt(filename, usecols=np.arange(14,63,2))
    ##########################################################################
    # Convert indices to EWs
    for idx in idxs:
        lick[:,idx] = w[idx] * (1 - np.power(10, -0.4 * lick[:,idx]))
        error[:,idx] = np.abs(0.4 * np.log(10) * (w[idx] - lick[:,idx]) * \
                               error[:,idx])
    ##########################################################################
    for i, spec in enumerate(specs):
        print spec
        dbname = "mcmc2_{0}_{1}".format(spec.replace(".fits", ""), modelname)
        if os.path.exists(dbname):
            continue
        idx = np.array([0,1,8,12,16,17,18,19])
        # idx = np.arange(25)
        run_mcmc(lick[i], error[i], modelname, idx, dbname)

def convert_tmj_to_ews():
    """ Convert tables from TMJ models to EWs."""
    os.chdir(os.path.join(tables_dir, "tmj"))
    outdir = os.path.join(tables_dir, "tmjew")
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    tables = ["tmj_Ca.dat", "tmj_C.dat", "tmj_Cr.dat", "tmj.dat",
                "tmj_Mg.dat", "tmj_Na.dat", "tmj_N.dat", "tmj_Si.dat",
                "tmj_Ti.dat"]
    bands = os.path.join(tables_dir, "bands.txt")
    ltype = np.loadtxt(bands, usecols=(8,))
    idxs = np.where(ltype==1)[0] + 3
    w = np.diff(np.loadtxt(bands, usecols=(4,5))).T[0]
    for intable in tables:
        with open(intable, "r") as f:
            lines = f.readlines()
        comments = [x for x in lines if x.startswith("#")]
        lick = np.loadtxt(intable)
        for idx in idxs:
            lick[:,idx] = w[idx-3] * (1 - np.power(10, -0.4 * lick[:,idx]))
        outtable = os.path.join(outdir, intable)
        with open(outtable, "w") as f:
            f.write("".join(comments))
            np.savetxt(f, lick)
    return

def summary_table(specs, modelname, db):
    """ Make final table."""
    lines = []
    for spec in specs:
        folder = spec.replace(".fits", "_db{0}".format(db))
        logfile = os.path.join(working_dir, folder,
                               "summary.txt")
        if not os.path.exists(logfile):
            continue
        with open(logfile, "r") as f:
            header = f.readline()
            lines.append(f.readline())
    table = os.path.join(working_dir, "populations_{0}.txt".format(modelname))
    with open(table, "w") as f:
        f.write(header)
        f.write("\n".join(lines))

def run_analysis(chain=1):
    os.chdir(data_dir)
    modelname="TMJ10ext"
    dirs = sorted([x for x in os.listdir(".") if x.startswith("mcmc") and
                   os.path.isdir(x)])
    ssp = SSP(modelname)
    # pp = PdfPages(os.path.join(data_dir,
    #                            "mcmc_results_{0}.pdf".format(modelname)))
    # plt.figure(1, figsize=(9,6.5))
    # plt.minorticks_on()
    outtable = []
    for i, direc in enumerate(dirs):
        print "{0} / {1}".format(i+1, len(dirs))
        obj = "_".join(direc.split("_")[1:3])
        night = direc.split("_")[3]
        spec = "_".join(direc.split("_")[1:4]) + ".fits"
        ######################################################################
        # Reading traces
        trage = np.loadtxt(os.path.join(data_dir, direc,
                                    "Chain_{0}/age_dist.txt".format(chain)))
        trmetal = np.loadtxt(os.path.join(data_dir, direc,
                                    "Chain_{0}/metal_dist.txt".format(chain)))
        tralpha = np.loadtxt(os.path.join(data_dir, direc,
                                    "Chain_{0}/alpha_dist.txt".format(chain)))
        #######################################################################
        Age = Dist(trage, ssp.lims[0])
        Metal = Dist(trmetal, ssp.lims[1])
        Alpha = Dist(tralpha, ssp.lims[2])
        results = []
        for d in [Age, Metal, Alpha]:
            results += [d.best.MAPP, d.best.lerr, d.best.uerr]
        results = [spec] + ["{0:.5g}".format(x) for x in results]
        outtable.append(results)
        # results.append[[Age.best.mean, Age.lerr, Age.uerr,
        #                 Metal.best.mean, Metal.lerr, Metal.uerr,
        #                 Alpha.best.mean, Alpha.lerr, Alpha.uerr]]
        # print results
        #######################################################################
        # plot_results([Age, Metal, Alpha])
    #     plt.annotate("{0} {1}".format(obj.replace("_", " "), night), xy=(.7,.84),
    #                  xycoords="figure fraction", ha="center", size=20)
    #     pp.savefig()
    #     plt.show()
        # plt.clf()
    # pp.close()
    outtable = np.array(outtable)
    with open("populations_chain{0}.txt".format(chain), "w") as f:
        f.write("# Spec Age LERR UERR [Z/H] LERR UERR [alpha/Fe] LERR UERR\n")
        np.savetxt(f, outtable, fmt="%s")

def plot_results(dists):
    for i, d in enumerate(dists):
        ax = plt.subplot(3,3,(4*i)+1)
        N, bins, patches = plt.hist(d.data, color="b",ec="k", bins=30, \
                                    range=tuple(d.lims), normed=True, \
                                    edgecolor="k", histtype='bar',linewidth=1.)
        fracs = N.astype(float)/N.max()
        norm = Normalize(-.2* fracs.max(), 1.5 * fracs.max())
        for thisfrac, thispatch in zip(fracs, patches):
            color = cm.gray_r(norm(thisfrac))
            thispatch.set_facecolor(color)
            thispatch.set_edgecolor("w")
        x = np.linspace(d.data.min(), d.data.max(), 100)
        ylim = ax.get_ylim()
        plt.plot(x, d.best.pdf(x), "-r", lw=1.5, alpha=0.7)
        ax.set_ylim(ylim)
        plt.axvline(d.best.MAPP, c="r", ls="--", lw=1.5)
        plt.tick_params(labelright=True, labelleft=False, labelsize=10)
        plt.xlim(d.lims)
        plt.locator_params(axis='x',nbins=10)
        if i < 2:
            plt.setp(ax.get_xticklabels(), visible=False)
        else:
            plt.xlabel(r"[$\mathregular{\alpha}$ / Fe]")
        plt.minorticks_on()
    def hist2D(dist1, dist2):
        """ Plot distribution and confidence contours. """
        X, Y = np.mgrid[dist1.lims[0] : dist1.lims[1] : 20j,
                        dist2.lims[0] : dist2.lims[1] : 20j]
        extent = [dist1.lims[0], dist1.lims[1], dist2.lims[0], dist2.lims[1]]
        positions = np.vstack([X.ravel(), Y.ravel()])
        values = np.vstack([dist1.data, dist2.data])
        kernel = stats.gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)
        ax.imshow(np.rot90(Z), cmap="gray_r", extent=extent, aspect="auto",
                  interpolation="spline16")
        plt.axvline(dist1.best.MAPP, c="r", ls="--", lw=1.5)
        plt.axhline(dist2.best.MAPP, c="r", ls="--", lw=1.5)
        plt.tick_params(labelsize=10)
        ax.minorticks_on()
        plt.locator_params(axis='x',nbins=10)
        return
    ax = plt.subplot(3,3,4)
    hist2D(dists[0], dists[1])
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.ylabel("[Z/H]")
    plt.xlim(dists[0].lims)
    plt.ylim(dists[1].lims)
    ax = plt.subplot(3,3,7)
    hist2D(dists[0], dists[2])
    plt.ylabel(r"[$\mathregular{\alpha}$ / Fe]")
    plt.xlabel("log Age (yr)")
    plt.xlim(dists[0].lims)
    plt.ylim(dists[2].lims)
    ax = plt.subplot(3,3,8)
    plt.xlabel("[Z/H]")
    hist2D(dists[1], dists[2])
    plt.xlim(dists[1].lims)
    plt.ylim(dists[2].lims)
    return
    # Annotations
#     xys = [(.7,.84), (.7,.77), (.7,.70)]
#     line = r"{0:28s}".format(spec)
#     for j, par in enumerate([r"Log Age", r"[Z/H]", r"[$\alpha$/Fe]"]):
#         text = r"{0}={1[0]:.2f}$^{{+{1[1]:.2f}}}_"" \
#                ""{{-{1[2]:.2f}}}$ dex".format(par, summary[j])
#         plt.annotate(text, xy=xys[j], xycoords="figure fraction",
#                      ha="center", size=20)
#         line += "{0[1]:.5f}"
#     plt.tight_layout(pad=0.2)


def fit_interactive():
    """ Calculate stellar populations in candidates. """
    os.chdir(data_dir)
    modelname="TMJ10ext"
    bands = os.path.join(tables_dir, "bands.txt")
    ltype = np.loadtxt(bands, usecols=(8,))
    indices = np.loadtxt(bands, usecols=(0,), dtype=str)
    i0 = np.array([0,1,8,12,16,17,18,19])
    w = np.diff(np.loadtxt(bands, usecols=(4,5))).T[0]
    idxs = np.where(ltype==1)[0]
    filename = "results.tab"
    specs = np.loadtxt(filename, usecols=(0,), dtype=str)
    lick = np.loadtxt(filename, usecols=np.arange(13,62,2))
    error = np.loadtxt(filename, usecols=np.arange(14,63,2))
    ##########################################################################
    # Convert indices to EWs
    for idx in idxs:
        lick[:,idx] = w[idx] * (1 - np.power(10, -0.4 * lick[:,idx]))
        error[:,idx] = np.abs(0.4 * np.log(10) * (w[idx] - lick[:,idx]) * \
                               error[:,idx])
    ##########################################################################
    ssp = SSP(modelname, idx=i0)
    def fitfunc(p, fjac=None, x=None, y=None, err=None, model=None):
        status = 0
        return([status, (y-model(p)[0])/err])

    for i,spec in enumerate(specs):
        print spec
        p0 = np.array([5., 0., 0.])
        fa = {'y':lick[i][i0], 'err':error[i][i0], "model":ssp}
        m = mpfit.mpfit(fitfunc, p0, functkw=fa)
        print 'parameters = ', m.params
        raw_input()


if __name__ == "__main__":
    # convert_tmj_to_ews()
    # run_candidates()
    run_analysis()
    # fit_interactive()