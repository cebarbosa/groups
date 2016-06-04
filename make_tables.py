# -*- coding: utf-8 -*-
"""

Created on 11/05/16

@author: Carlos Eduardo Barbosa

Merge results from tables in single ASCII and Tex tables.

"""
import os

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import Angle
from astropy import units

from config import *
from loess.cap_loess_2d import loess_2d

def intersect_n_arrays(arrays, assume_unique=False):
    """Find the intersection of any number of 1D arrays.
    Return the sorted, unique values that are in all of the input arrays.
    Adapted from numpy.lib.arraysetops.intersect1d"""
    N = len(arrays)
    if N == 0:
        return np.asarray(arrays)
    arrays = list(arrays) # allow assignment
    if not assume_unique:
        for i, arr in enumerate(arrays):
            arrays[i] = np.unique(arr)
    aux = np.concatenate(arrays) # one long 1D array
    aux.sort() # sorted
    if N == 1:
        return aux
    shift = N-1
    return aux[aux[shift:] == aux[:-shift]]

def merge_tables():
    """ Make single table containing kinematics and stellar populations. """
    kintable = os.path.join(data_dir, "ppxf_results.dat")
    licktab = os.path.join(data_dir, "lick.txt")
    lickerrtab = os.path.join(data_dir, "lickerr_mc50.txt")
    poptable = os.path.join(data_dir, "populations_chain1.txt")
    tables = [kintable, licktab, lickerrtab, poptable]
    ##########################################################################
    specs = [np.loadtxt(x, usecols=(0,), dtype=str) for x in tables]
    ref = intersect_n_arrays(specs)
    data = [np.loadtxt(x, dtype=str) for x in tables]
    idx = [np.searchsorted(x[:,0], ref) for x in data]
    data = [d[i] for (i,d) in zip(idx,data)]
    ###########################################################################
    # Interwine Lick indices and errors
    lick = np.insert(data[2], np.arange(len(data[1][0])), data[1], axis=1)
    data[1] = lick[:,1:]
    del data[2]
    ###########################################################################
    combined = []
    for i in range(len(data)):
        if i == 0:
            combined = data[i]
        else:
            combined = np.column_stack((combined, data[i][:,1:]))
    ##########################################################################
    # Obtaining coordinates of each galaxy
    radec = np.loadtxt(os.path.join(tables_dir, "candidates_radec.dat"),
                       dtype=str, skiprows=1)
    radec = dict([(x,[y,z]) for x,y,z in radec])
    gal = ["_".join(x.split("_")[:2]) for x in combined[:,0]]
    ra, dec = np.array([radec[x] for x in gal]).T
    ra =  np.array([str(x) for x in Angle(ra, unit=units.hourangle)])
    dec = np.array([str(x) for x in Angle(dec, unit=units.degree)])
    combined = np.column_stack((combined, ra, dec))
    ###########################################################################
    # Making headers
    header = ["File", 'V', 'dV', 'S', 'dS', 'h3', 'dh3', 'h4', 'dh4',
              'chi/DOF', 'S/N/Angstrom', 'ADEGREE', 'MDEGREE']
    indices = np.loadtxt(os.path.join(tables_dir, "bands.txt"), usecols=(0,),
                         dtype=str)
    err = np.array(25 * ["err    "], dtype=str)
    indices = np.insert(err, np.arange(25), indices)
    header += indices.tolist()
    header += ["Age", "LERR", "UERR", "[Z/H]", "LERR", "UERR", "[alpha]/Fe",
              "LERR", "UERR"]
    header += ["RAJ2000", "DECJ2000"]
    header = ["# ({0}) {1}".format(i,h) for i,h in enumerate(header)]
    ###########################################################################
    with open(os.path.join(data_dir, "results.tab"), "w") as f:
        f.write("\n".join(header) + "\n")
        np.savetxt(f, combined, fmt="%s")
    return

def comment_table():
    table = os.path.join(data_dir, "results.tab")
    with open(table) as f:
        header = [x for x in f.readlines() if x.startswith("#")]
    data = np.loadtxt(table, dtype=str)
    comment = ["ngc193_119_blanco11bn3.fits", "ngc7619_166_blanco11bn3.fits"]
    specs = data[:,0]
    comments = np.array(["" if x not in comment else "#" for x in specs])
    data = np.column_stack((comments, data))
    with open(table, "w") as f:
        f.write("".join(header))
        np.savetxt(f, data, fmt="%s")


def print_stats():
    table = os.path.join(data_dir, "results.tab")
    specs = np.loadtxt(table, usecols=(0,), dtype=str)
    objs = ["_".join(x.split("_")[:2]) for x in specs]
    galaxies = np.unique(objs)
    print "total of analysed spectra: ", len(objs)
    print "Total of galaxies:", len(galaxies)
    for group in groups:
        members = [x for x in galaxies if x.startswith(group)]
        print "Number of members in {0}: {1}".format(group, len(members))

def sigma_mgb():
    sigma, sigerr, mgb, mgberr = np.loadtxt(os.path.join(data_dir,
                                 "results.tab"), usecols=(3,4,45,46)).T
    fig = plt.figure(1)
    ax = plt.subplot(111)
    ax.errorbar(np.log10(sigma), mgb, xerr=sigerr/(sigma * np.log(10.)) ,
                yerr=mgberr, fmt="o", ecolor="0.7",
                capsize=0)
    plt.show()

def make_latex():
    os.chdir(data_dir)
    table = "results.tab"
    specs = np.loadtxt(table, usecols=(0,), dtype=str)
    objs = ["_".join(x.split("_")[:2]) for x in specs]
    objs = [x.replace("_", "\_").upper() for x in objs]
    runs = ["({0})".format(obsrun[x.replace(".fits", "").split("_")[2]]) for x in specs]
    ra, dec = np.loadtxt(table, usecols=(72,73), dtype=str).T
    ra = np.array([str(x) for x in Angle(ra, unit=units.hour)])
    dec = np.array([str(x) for x in Angle(dec, unit=units.degree)])
    out = np.column_stack((objs, runs, ra, dec))
    # Kinematics table
    kcols = (1,3,5,7)
    for col in kcols:
        kdata = np.loadtxt(table, usecols=(col,))
        kerr = np.loadtxt(table, usecols=(col+1,))
        out = np.column_stack((out, make_str_err(kdata, kerr)))
    sn = np.loadtxt(table, usecols=(10,))
    sn = ["{0:.1f}".format(x) for x in sn]
    out = np.column_stack((out, sn))
    out = [" & ".join(x) + "\\\\" for x in out]
    with open("groups_kinematics.tex", "w") as f:
        f.write("\n".join(out))
    # icols = np.arange(13, 62, 2)
    # print len(cols)
    # for col in icols:
    #     idata = np.loadtxt(table, usecols=(col,))
    #     ierr = np.loadtxt(table, usecols=(col+1,))
    #     out = np.column_stack((out, make_str_err(idata, ierr)))



def make_str_err(vals, errs):
    strs = []
    for val, err in zip(vals, errs):
        if np.isnan(val) or np.isnan(err):
            strs.append("--")
        elif err >= 1:
            strs.append("${0}\pm{1}$".format(int(round(val)), int(round(err))))
        else:
            ndig = int(np.ceil(np.abs(np.log10(err))))
            ndig = ".{0}".format(ndig)
            strs.append("${0:{2}f}\pm{1:{2}f}$".format(val, err,ndig))
    return strs

    # val = float(d[7+3*i])
    # lerr =  np.abs(float(d[8+3*i]))
    # uerr = np.abs(float(d[9+3*i]))
    # ndig = np.minimum(np.log10(lerr), np.log10(uerr))
    # ndig = int(np.ceil(np.abs(ndig)))
    # ndig = ".{0}".format(ndig)
    # valstr = "{4}=${0:{3}f}_{{-{1:{3}f}}}^{{+{2:{3}f}}}$".format(val, lerr,
    #                                             uerr, ndig, prop[i])
    # stpop.append(valstr)

def internal_check(verbose=False):
    table = os.path.join(data_dir, "results.tab")
    specs = np.loadtxt(table, usecols=(0,), dtype=str)
    data = np.loadtxt(table, usecols=(1,3))
    err = np.loadtxt(table, usecols=(2,4))
    objs = np.array(["_".join(x.split("_")[:2]) for x in specs])
    vgroups = np.array([v0s[x.split("_")[0]] for x in specs])
    data[:,0] -= vgroups
    galaxies = np.unique(objs)
    repdata, reperr = [], []
    for gal in galaxies:
        idx = np.where(objs == gal)[0]
        if len(idx) > 1:
            if verbose:
                print gal, data[idx].T.ravel()
            if type(repdata) is list:
                repdata = np.array(data[idx].T.ravel())
                reperr = np.array(err[idx].T.ravel())
            else:
                repdata = np.column_stack((repdata, data[idx].T.ravel()))
                reperr = np.column_stack((reperr, err[idx].T.ravel()))
    fig = plt.figure(figsize=(8,3.5))
    xlims = [[-1000, 1000], [0,350]]
    ylims=[[-50,50],[-80,80]]
    xlabels = ["$\mathbf{V-V_{group}}$ (km/s)", "$\mathbf{\sigma}$ (km/s)"]
    ylabels = ["$\mathbf{\Delta V}$ (km/s)", "$\mathbf{\Delta \sigma}$ (km/s)"]
    for i,(i1,i2) in enumerate(np.arange(4).reshape(2,2)):
        ax = plt.subplot(1,2,i+1)
        ax.minorticks_on()
        ax.errorbar(repdata[i1], repdata[i2] - repdata[i1], xerr=None,
                    yerr=np.sqrt(reperr[i2]**2 + reperr[i1]**2),
                    ecolor="0.7", color="r", mec="r", fmt="o")
        ax.set_xlim(xlims[i][0], xlims[i][1])
        ax.set_ylim(ylims[i][0], ylims[i][1])
        # ax.plot(np.linspace(lims[i][0], lims[i][1], 100),
        #         np.linspace(lims[i][0], lims[i][1], 100), "--k")
        ax.axhline(ls="--", c="k")
        ax.set_xlabel(xlabels[i])
        ax.set_ylabel(ylabels[i])
    fig.subplots_adjust(bottom=0.16, hspace=0.3, right=0.98, wspace=0.35,
                        top=0.95, left=0.1)
    plt.savefig(os.path.join(plots_dir, "internal_check.png"))
    # plt.show()

def mass_pop():
    table = os.path.join(data_dir, "results.tab")
    sigma, sigerr = np.loadtxt(table, usecols=(3,4)).T
    logsigma = np.log10(sigma)
    logsigerr = np.abs(sigerr / sigma / np.log(10))
    stpops = np.loadtxt(table, usecols=(63,66,69)).T
    lerrs = np.abs(np.loadtxt(table, usecols=(64,67,70)).T)
    uerrs = np.abs(np.loadtxt(table, usecols=(65,68,71)).T)
    lerrs[0] = np.abs(lerrs[0]/stpops[0]/np.log(10))
    uerrs[0] = np.abs(uerrs[0]/stpops[0]/np.log(10))
    stpops[0] = np.log10(stpops[0])
    fig = plt.figure(figsize=(11,3.5))
    deg = [3, 2, 2]
    ylabels = ["$\log$ Age (Gyr)", "[Z/H]", r"[$\alpha$/Fe]"]
    ylim = [[None, 1.2], [None, None], [None, None]]
    for i, (stpop, lerr, uerr) in enumerate(zip(stpops, lerrs, uerrs)):
        ax = plt.subplot(1,3,i+1)
        ax.minorticks_on()
        ax.errorbar(logsigma, stpop, xerr=logsigerr,
                    yerr=[lerr, uerr],
                    ecolor="0.7", color="r", mec="r", fmt="o", ms=6)
        # zout, wout = loess_2d(logsigma, stpop, np.power(10,stpops[0]), degree=deg[i])
        # ax.scatter(logsigma, stpop, c=zout, cmap="rainbow", s=40,
        #            zorder=100)
        ax.set_xlim(0., 3)
        ax.set_ylim(ylim[i])
        ax.set_xlabel("$\mathbf{\log\sigma}$ (km/s)")
        ax.set_ylabel(ylabels[i])
        # ax.set_ylim(stpop.min(), stpop.max())
    fig.subplots_adjust(bottom=0.16, hspace=0.3, right=0.98, wspace=0.35,
                        top=0.95, left=0.1)
    plt.savefig(os.path.join(plots_dir, "sigma_stpop.png"))

def make_unique():
    """ Make table averaging repeated galaxies."""
    table = os.path.join(data_dir, "results.tab")
    with open(table) as f:
        header = [x for x in f.readlines() if x.startswith("#")]
    specs =  np.loadtxt(table, usecols=(0,), dtype=str)
    specs = np.array(["_".join(x.split("_")[:2]) for x in specs])
    data = np.loadtxt(table, dtype=str)
    data[:,0] = specs
    # Defining index type: 0 = strings, 1 = mean, 2 = error
    i0 = np.array([0,72,73])
    i1 = np.hstack(([1,3,5,7,9,10,11,12], np.arange(13,62,2)))
    i2 = np.hstack(([2,4,6,8,14,16], np.arange(14,63,2)))
    i3 = np.arange(63,72)
    new = []
    unique = np.unique(specs)
    results = np.empty((len(unique), len(data[0])), dtype="S32")
    for i,gal in enumerate(unique):
        idx = np.where(gal == specs)[0]
        if len(idx) == 1:
            results[i] = data[idx][0]
        else:
            subdata = data[idx]
            newdata = subdata[0]
            newdata[i1] = subdata[:,i1].astype(float).mean(axis=0).astype(str)
            error = np.sqrt(np.sum(subdata[:,i2].astype(float)**2, axis=0))/ \
                    np.sqrt(len(idx))
            newdata[i2] = error.astype(str)
            results[i] = newdata
    output = os.path.join(data_dir, "results_unique.tab")
    with open(output, "w") as f:
        f.write("".join(header))
        np.savetxt(f, results, fmt="%s")
    return



if __name__ == "__main__":
    merge_tables()
    comment_table()
    make_latex()
    make_unique()
    # print_stats()
    # sigma_mgb()
    # internal_check()
    # mass_pop()