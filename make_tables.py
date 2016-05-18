# -*- coding: utf-8 -*-
"""

Created on 11/05/16

@author: Carlos Eduardo Barbosa

Merge results from tables in single ASCII and Tex tables.

"""
import os

import numpy as np
import matplotlib.pyplot as plt

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
    tables = [kintable, licktab, lickerrtab]
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


if __name__ == "__main__":
    merge_tables()
    # comment = ["hcg42_011_blanco11an2", "hcg42_021_blanco11an2",
    #
    #            "hcg62_66_blanco10n1",
    #            "hcg62_67_blanco10n2", "hcg62_n0127_blanco10n2",
    #            "hcg90_eso466d_blanco11bn3", "ngc193_119_blanco11bn2"]
    # sigma_mgb()
