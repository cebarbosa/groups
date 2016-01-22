# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 12:10:54 2014
Adapted to run in data frou groups of galaxies in Dec 22, 2015

@author: kadu

Run pPXF in data
"""
import os
import pickle

import numpy as np
import pyfits as pf
from scipy.signal import medfilt
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import gridspec

from ppxf import ppxf
import ppxf_util as util
from config import *
 
def run_ppxf(spectra, velscale, ncomp=None, has_emission=True, mdegree=-1,
             degree=20, pkls=None, plot=False, sky=None, start=None):
    """ Run pPXF in a list of spectra"""
    if isinstance(spectra, str):
        spectra = [spectra]
    if isinstance(pkls, str):
        pkls = [pkls]
    if pkls == None:
        pkls = [x.replace(".fits", ".pkl") for x in spectra]
    ##########################################################################
    # Load templates for both stars and gas
    star_templates = pf.getdata(os.path.join(templates_dir,
                                             'miles_FWHM_3.6.fits'), 0)
    logLam2 = pf.getdata(os.path.join(templates_dir, 'miles_FWHM_3.6.fits'), 1)
    miles = np.loadtxt(os.path.join(templates_dir, 'miles_FWHM_3.6.txt'),
                       dtype=str).tolist()
    gas_templates = pf.getdata(os.path.join(templates_dir,
                                             'emission_FWHM_3.6.fits'), 0)
    logLam_gas = pf.getdata(os.path.join(templates_dir, 'emission_FWHM_3.6.fits'),
                            1)
    gas_files = np.loadtxt(os.path.join(templates_dir, 'emission_FWHM_3.6.txt'),
                       dtype=str).tolist()

    ngas = len(gas_files)
    ntemplates = len(miles)
    ##########################################################################
    # Join templates in case emission lines are used.
    if has_emission:
        templates = np.column_stack((star_templates, gas_templates))
        templates_names = np.hstack((miles, gas_files))
    else:
        templates = star_templates
        templates_names = miles
        ngas = 0
    ##########################################################################
    if sky == None:
        nsky = 0
    else:
        nsky = sky.shape[1]
    if ncomp == 1:
        components = 0
        moments = [4]
    elif ncomp == 2:
        components = np.hstack((np.zeros(len(star_templates[0])),
                                np.ones(len(gas_templates[0]))))
        moments = [4,2]
    for i, spec in enumerate(spectra):
        print "pPXF run of spectrum {0} ({1} of {2})".format(spec, i+1,
              len(spectra))
        pkl = pkls[i]
        plt.clf()
        ######################################################################
        # Read galaxy spectrum and define the wavelength range
        hdu = pf.open(spec)
        spec_lin = hdu[0].data
        h1 = pf.getheader(spec)
        lamRange1 = h1['CRVAL1'] + np.array([0.,h1['CDELT1']*(h1['NAXIS1']-1)])
        ######################################################################
        # Rebin to log scale
        galaxy, logLam1, velscale = util.log_rebin(lamRange1, spec_lin, 
                                                   velscale=velscale)
        ######################################################################
        # First guess for the noise
        noise = np.ones_like(galaxy) * np.std(galaxy - medfilt(galaxy, 5))
        ######################################################################
        # Calculate difference of velocity between spectrum and templates
        # due to different initial wavelength
        dv = (logLam2[0]-logLam1[0])*c
        ######################################################################
        # Set first guess from setup files
        if start == None:
            start = [v0s[spec.split("_")[0]], 100]
        goodPixels = None
        ######################################################################
        # Expand start variable to include multiple components
        if ncomp > 1:
            start = [start, [start[0], 30]]
        ######################################################################
        # First pPXF interaction
        pp0 = ppxf(templates, galaxy, noise, velscale, start,
                   goodpixels=goodPixels, plot=False, moments=moments,
                   degree=12, mdegree=-1, vsyst=dv, component=components,
                   sky=sky)
        rms0 = galaxy[goodPixels] - pp0.bestfit[goodPixels]
        noise0 = 1.4826 * np.median(np.abs(rms0 - np.median(rms0)))
        noise0 = np.zeros_like(galaxy) + noise0
        # Second pPXF interaction, realistic noise estimation
        pp = ppxf(templates, galaxy, noise0, velscale, start,
                  goodpixels=goodPixels, plot=False, moments=moments,
                  degree=degree, mdegree=mdegree, vsyst=dv,
                  component=components, sky=sky)
        plt.title(spec.replace("_", "-"))
        plt.show(block=False)
        plt.savefig("logs/{0}".format(spec.replace(".fits", ".png")))
        ######################################################################
        # Adding other things to the pp object
        pp.template_files = templates_names
        pp.has_emission = has_emission
        pp.dv = dv
        pp.w = np.exp(logLam1)
        pp.velscale = velscale
        pp.spec = spec
        pp.ngas = ngas
        pp.ntemplates = ntemplates
        pp.nsky = nsky
        pp.templates = 0
        ######################################################################
        # Save fit to pickles file to keep session
        # print "Saving to file {0}".format(pkl)
        # with open(pkl, "w") as f:
        #     pickle.dump(pp, f)
        ######################################################################
        ppf = pPXF(spec, velscale, pp)
        logdir = os.path.join(os.getcwd(), "logs")
        if not os.path.exists(logdir):
            os.mkdir(logdir)
        ppf.plot("logs/{0}".format(spec.replace(".fits", ".png")))
        ######################################################################
        # # Save to output text file
        sol = [val for pair in zip(pp.sol, pp.error) for val in pair]
        sol = ["{0:12s}".format("{0:.3g}".format(x)) for x in sol]
        sol.append("{0:12s}".format("{0:.3g}".format(pp0.chi2)))
        sol = ["{0:30s}".format(spec)] + sol
        with open(pkl.replace(".pkl", ".txt"), "w") as f:
            f.write("".join(sol))
    return

def read_setup_file(gal, logw, mask_emline=True):
    """ Read setup file to set first guess and regions to be avoided. """
    w = np.exp(logw)
    filename = gal + ".setup"
    with open(filename) as f:
        f.readline()
        start = f.readline().split()
    start = np.array(start, dtype=float)
    ranges = np.loadtxt(filename, skiprows=5)
    ##########################################################################
    # Mask all the marked regions in the setup file
    if mask_emline:
        for i, (w1, w2) in enumerate(ranges.reshape((len(ranges)/2, 2))):
            if i == 0:
                good = np.where(np.logical_and(w > w1, w < w2))[0]
            else:
                good = np.hstack((good, np.where(np.logical_and(w > w1,
                                                                w < w2))[0]))
        good = np.array(good)
        good.sort()

    ##########################################################################
    # Mask only regions in the beginning and in the end of the spectra plus
    # the residuals in the emission line at 5577 Angstroms
    else:
        ranges = [[np.min(ranges), 5577. - 15], [5577. + 15, np.max(ranges)]]
        for i, (w1, w2) in enumerate(ranges):
            if w1 >= w2:
                continue
            if i == 0:
                good = np.where(np.logical_and(w > w1, w < w2))[0]
            else:
                good = np.hstack((good, np.where(np.logical_and(w > w1,
                                                                w < w2))[0]))
        good = np.array(good)
        good.sort()
        return start, good

def make_table(spectra, output, mc=False, nsim=200, clean=True, pkls=None):
    """ Make table with results.

    ===================
    Input Parameters
    ===================
    spectra : list
        Names of the spectra to be processed. Should end in "fits".
    mc : bool
        Calculate the errors using a Monte Carlo method.
    nsim : int
        Number of simulations in case mc keyword is True.
    clean : bool
        Remove lines for which the velocity dispersion is 1000 km/s.
    pkls : list
        Specify list of pkl files to be used. Default value replaces fits for
        pkl
    ==================
    Output file
    ==================
    In case mc is False, the function produces a file called ppxf_results.dat.
    Otherwise, the name of the file is named ppfx_results_mc_nsim.dat.

    """
    print "Producing summary table..."
    head = ("{0:<30}{1:<14}{2:<14}{3:<14}{4:<14}{5:<14}{6:<14}{7:<14}"
             "{8:<14}{9:<14}{10:<14}{11:<14}{12:<14}{13:<14}\n".format("# FILE",
             "V", "dV", "S", "dS", "h3", "dh3", "h4", "dh4", "chi/DOF",
             "S/N (/ pixel)", "ADEGREE", "MDEGREE", "100*S/N/sigma"))
    results = []
    ##########################################################################
    # Initiate the output file
    ##########################################################################
    with open(output, "w") as f:
        f.write(head)
    ##########################################################################
    if pkls== None:
        pkls = [x.replace(".fits", ".pkl") for x in spectra]
    for i, (spec, pkl) in enumerate(zip(spectra, pkls)):
        print "Working on spectra {0} ({1} / {2})".format(spec, i+1,
                                                          len(spectra))
        if not os.path.exists(spec.replace(".fits", ".pkl")):
            continue
        pp = pPXF(spec, velscale, pklfile=pkl)
        pp.calc_sn()
        sn = pp.sn
        if mc:
            pp.mc_errors(nsim=nsim)
        if pp.ncomp > 1:
            pp.sol = pp.sol[0]
            pp.error = pp.error[0]
        data = [pp.sol[0], pp.error[0],
                pp.sol[1], pp.error[1], pp.sol[2], pp.error[2], pp.sol[3],
                pp.error[3], pp.chi2, sn]
        data = ["{0:12.3f} ".format(x) for x in data]
        if clean and pp.sol[1] == 1000.:
            comment = "#"
        else:
            comment = ""
        line = ["{0}{1}".format(comment, spec)] + data + \
               ["{0:12}".format(pp.degree), "{0:12}".format(pp.mdegree),
                "{0:12.3f}".format(100 * sn / pp.sol[1])]
        # Append results to outfile
        with open(output, "a") as f:
            f.write("".join(line) + "\n")

class pPXF():
    """ Class to read pPXF pkl files """
    def __init__(self, spec, velscale, pp):
        self.__dict__ = pp.__dict__.copy()
        self.calc_arrays()
        self.calc_sn()
        return

    def calc_arrays(self):
        """ Calculate the different useful arrays."""
        # Slice matrix into components
        self.m_poly = self.matrix[:,:self.degree + 1]
        self.matrix = self.matrix[:,self.degree + 1:]
        self.m_ssps = self.matrix[:,:self.ntemplates]
        self.matrix = self.matrix[:,self.ntemplates:]
        self.m_gas = self.matrix[:,:self.ngas]
        self.matrix = self.matrix[:,self.ngas:]
        self.m_sky = self.matrix
        # Slice weights
        self.w_poly = self.polyweights
        self.w_ssps = self.weights[:self.ntemplates]
        self.weights = self.weights[self.ntemplates:]
        self.w_gas = self.weights[:self.ngas]
        self.weights = self.weights[self.ngas:]
        self.w_sky = self.weights
        # Calculating components
        self.poly = self.m_poly.dot(self.w_poly)
        self.ssps = self.m_ssps.dot(self.w_ssps)
        self.gas = self.m_gas.dot(self.w_gas)
        self.bestsky = self.m_sky.dot(self.w_sky)
        return

    def calc_sn(self, w1=4200., w2=6300.):
        idx = np.logical_and(self.w >=w1, self.w <=w2)
        self.res = self.galaxy - self.bestfit
        # Using robust method to calculate noise using median deviation
        self.noise = 1.4826 * np.median(np.abs(self.res[idx] -
                                               np.median(self.res[idx])))
        self.signal = np.sum(self.ssps[idx]) / len(self.ssps[idx])
        self.sn = self.signal / self.noise
        return

    def mc_errors(self, nsim=200):
        """ Calculate the errors using MC simulations"""
        errs = np.zeros((nsim, len(self.error)))
        for i in range(nsim):
            y = self.bestfit + np.random.normal(0, self.noise,
                                                len(self.galaxy))

            noise = np.ones_like(self.galaxy) * self.noise
            sim = ppxf(self.bestfit_unbroad, y, noise, velscale,
                       [0, self.sol[1]],
                       goodpixels=self.goodpixels, plot=False, moments=4,
                       degree=-1, mdegree=-1,
                       vsyst=self.vsyst, lam=self.lam, quiet=True, bias=0.)
            errs[i] = sim.sol
        median = np.ma.median(errs, axis=0)
        error = 1.4826 * np.ma.median(np.ma.abs(errs - median), axis=0)
        # Here I am using always the maximum error between the simulated
        # and the values given by pPXF.
        self.error = np.maximum(error, self.error)
        return

    def plot(self, output, xlims = [4000, 6500], fignumber=1, textsize = 16):
        """ Plot pPXF run in a output file"""
        plt.figure(fignumber)
        plt.clf()
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        ax = plt.subplot(gs[0])
        ax.minorticks_on()
        if self.ncomp > 1:
            sol = self.sol[0]
            error = self.error[0]
            sol2 = self.sol[1]
            error2 = self.error[1]
        else:
            sol = self.sol
            error = self.error
        ax.plot(self.w, self.galaxy, "-k", lw=2., label=self.spec.replace("_", \
                                                    "-").replace(".fits", ""))
        ax.plot(self.w[self.goodpixels], self.ssps[self.goodpixels] + self.poly,
                "-", lw=1, c="r", label="SSPs")
        if self.has_emission:
            ax.plot(self.w[self.goodpixels], self.gas[self.goodpixels], "-b",
                    lw=1., label="Emission Lines")
        if self.sky != None:
            ax.plot(self.w[self.goodpixels], self.bestsky[self.goodpixels], \
                    "-", lw=1, c="g", label="Sky")
        plt.axhline(y=0, ls="--", c="k")
        plt.ylabel(r"Flux (Counts)", size=18)
        plt.tight_layout()
        plt.annotate(r"$\chi^2=${0:.2f}".format(self.chi2),
                     xycoords='axes fraction',
                    xy=(0.05,0.9), size=textsize)
        plt.annotate(r"S/N={0}".format(np.around(self.sn,1)),
                     xycoords='axes fraction', xy=(0.25,0.9), size=textsize)
        plt.annotate(r"V={0} km/s".format(np.around(sol[0])),
                     xycoords='axes fraction', xy=(0.45,0.9), size=textsize,
                     color="r")
        plt.annotate(r"$\sigma$={0} km/s".format(np.around(sol[1])),
                     xycoords='axes fraction', xy=(0.75,0.9), size=textsize,
                     color="r")
        if self.ncomp > 1:
            plt.annotate(r"V={0} km/s".format(np.around(sol2[0])),
                         xycoords='axes fraction', xy=(0.45,0.88),
                         size=textsize, color="b")
            plt.annotate(r"$\sigma$={0} km/s".format(np.around(sol2[1])),
                         xycoords='axes fraction', xy=(0.75,0.88),
                         size=textsize, color="b")
        leg = plt.legend(loc=6, prop={"size":10})
        leg.get_frame().set_linewidth(0.0)
        ax1 = plt.subplot(gs[1])
        ax1.minorticks_on()
        ax.set_xlim(self.w[0], self.w[-1])
        ax1.set_xlim(self.w[0], self.w[-1])
        ax1.plot(self.w[self.goodpixels], self.galaxy[self.goodpixels] - \
                 self.bestfit[self.goodpixels], "-k")
        ax1.axhline(y=0, ls="--", c="k")
        ax1.axhline(y=self.noise, ls="--", c="r")
        ax1.axhline(y=-self.noise, ls="--", c="r")
        ax1.set_xlabel(r"$\lambda$ ($\AA$)", size=18)
        ax.xaxis.set_ticklabels([])
        ax1.set_ylabel(r"$\Delta$Flux", size=18)
        gs.update(hspace=0.075, left=0.125, bottom=0.1, top=0.98, right=0.97)
        plt.savefig(output)




def fits_to_matrix(filenames):
    """ Load several fits of the same dimension into a matrix.  """
    f = pf.getdata(filenames[0])
    a = np.zeros((f.shape[0], len(filenames)))
    for i in range(len(filenames)):
        try:
            a[:,i] = pf.getdata(filenames[i])
        except:
            print filenames[i]
    return a

def load_sky(filenames, velscale):
    """ Load and rebin sky files. """
    skydata = fits_to_matrix(filenames)
    h1 = pf.getheader(filenames[0])
    lamRange1 = h1['CRVAL1'] + np.array([0.,h1['CDELT1']*(h1['NAXIS1']-1)])
    sky1, logLam1, velscale = util.log_rebin(lamRange1, skydata[:,0],
                                                   velscale=velscale)

    skylog = np.zeros((sky1.shape[0], len(filenames)))
    for i in range(len(filenames)):
        skylog[:,i], logLam1, velscale = util.log_rebin(lamRange1,
                                        skydata[:,i], velscale=velscale)
    return skylog

def make_table_from_txt():
    """ Make a summary table using the txt outputs. """
    nights = sorted(os.listdir(data_dir))
    for night in nights:
        wdir = os.path.join(data_dir, night)
        os.chdir(wdir)
        specs = [x for x in os.listdir(".") if x.endswith(".fits")]
        txts = [x.replace(".fits", ".txt") for x in specs]
        txts = [x for x in txts if  os.path.exists(os.path.join(wdir, x))]
        txts.sort()
        with open("ppxf_results.txt", "w") as fout:
            for txt in txts:
                with open(txt) as fin:
                    fout.write(fin.read() + "\n")
    return

def select_specs():
    """ Select which spectra can be used in the analysis. """
    nights = sorted(os.listdir(data_dir))
    for night in nights:
        wdir = os.path.join(data_dir, night, "logs")
        os.chdir(wdir)
        pngs = sorted([x for x in os.listdir(".") if x.endswith(".png")])
        comments = []
        for image in pngs:
            img = mpimg.imread(image)
            plt.imshow(img)
            plt.axis("off")
            plt.pause(0.001)
            plt.show()
            comm = raw_input("Skip this spectrum in the anaylis? (y/N) ")
            if comm.lower().strip() in ["y", "ye", "yes"]:
                comments.append("#")
            else:
                comments.append("")
            plt.clf()
        ignorelist = ["{0}{1}".format(x,y.replace(".png", ".fits")) for \
                      x,y in zip(comments, pngs)]
        output = os.path.join(data_dir, night, "ignore.txt")
        with open(output, "w") as f:
            f.write("\n".join(ignorelist))

def run_over_all():
    """ Run pPXF in a generic way over all data. """
    nights = sorted(os.listdir(data_dir))
    for night in nights:
        print "Working in run ", night
        wdir = os.path.join(data_dir, night)
        os.chdir(wdir)
        log_dir = os.path.join(wdir, "logs")
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)
        fits = [x for x in os.listdir(".") if x.endswith(".fits")]
        skies =  [x for x in fits if x.startswith("sky")]
        specs = [x for x in fits if x not in skies]
        specs.sort()
        skies.sort()
        sky = load_sky(skies, velscale)
        # #################################################################
        # # Go to the main routine of fitting
        run_ppxf(specs, velscale, ncomp=1, has_emission=1, mdegree=-1,
                 degree=12, plot=True, sky=sky)
    return

def run_list(night, specs):
    """ Run pPXF on a given list of spectra of the same night. """
    wdir = os.path.join(data_dir, night)
    os.chdir(wdir)
    fits = [x for x in os.listdir(".") if x.endswith(".fits")]
    skies =  sorted([x for x in fits if x.startswith("sky")])
    sky = load_sky(skies, velscale)
    run_ppxf(specs, velscale, ncomp=1, has_emission=1, mdegree=-1,
                 degree=12, plot=True, sky=sky)

if __name__ == '__main__':
    # run_list("blanco10n1", ["hcg62_4.fits"])
    run_over_all()



