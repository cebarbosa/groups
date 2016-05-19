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
    coords = [radec[x] for x in gal]
    combined = np.column_stack((combined, coords))
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


if __name__ == "__main__":
    # merge_tables()
    comment = ["hcg42_011_blanco11an2", "hcg42_021_blanco11an2",
               "hcg62_66_blanco10n1", "hcg62_67_blanco10n2",
               "hcg62_n0127_blanco10n2", "hcg90_eso466d_blanco11bn3",
               "ngc193_119_blanco11bn2"]
    make_latex()
    # sigma_mgb()
