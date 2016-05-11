# -*- coding: utf-8 -*-
"""

Created on 11/05/16

@author: Carlos Eduardo Barbosa

Merge results from tables in single ASCII and Tex tables.

"""
import os

import numpy as np

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
    kintable = os.path.join(data_dir, "ppxf_results.dat")
    licktab = os.path.join(data_dir, "lick.txt")
    lickerrtab = os.path.join(data_dir, "lickerr_mc50.txt")
    tables = [kintable, licktab, lickerrtab]
    specs = [np.loadtxt(x, usecols=(0,), dtype=str) for x in tables]
    ref = intersect_n_arrays(specs)
    data = [np.loadtxt(x, dtype=str) for x in tables]
    idx = [np.searchsorted(ref, x[:,0]) for x in data]
    data = [d[i] for (i,d) in zip(idx,data)]
    print data[0][:,0]
    print data[1][:,0]




#     s1 = np.genfromtxt(files[0], usecols=(0,), dtype=None).tolist()
#     s2 = np.genfromtxt(files[1], usecols=(0,), dtype=None).tolist()
#     s3 = np.genfromtxt(files[2], usecols=(0,), dtype=None).tolist()
#     s4 = np.genfromtxt(files[3], usecols=(0,), dtype=None).tolist()
#     s5 = np.genfromtxt(files[4], usecols=(0,), dtype=None).tolist()
#     s6 = np.genfromtxt(files[5], usecols=(0,), dtype=None).tolist()
#     sref = list(set(s1) & set(s2) & set(s3) & set(s4) & set(s5) & set(s6))
#     ignore = ["fin1_n3311{0}.fits".format(x) for x in ignore_slits]
#     sref = [x for x in sref if x not in ignore]
#     if os.path.exists("template_mismatches.dat"):
#         temp_mismatches = np.loadtxt("template_mismatches.dat",
#                                      dtype=str).tolist()
#         sref = [x for x in sref if x not in temp_mismatches]
#     sref.sort()
#     # if workdir in [data_dir, minus_1pc_dir, plus_1pc_dir, best_dir,
#     #                rerun_dir]:
#     x,y = get_positions(sref).T
#     coords = get_coords(sref)
#     # elif workdir == binning_dir:
#     #     sref = ["s{0}.fits".format(i) for i in range(500)
#     #              if "s{0}.fits".format(i) in sref]
# #        sref = ["s{0}_v3800.fits".format(i) for i in range(500)
# #                 if "s{0}_v3800.fits".format(i) in sref]
# #         coords = np.loadtxt("spec_xy.txt", usecols=(3,4))
# #         x, y = coords.T
#     r = np.sqrt(x*x + y*y)
#     pa = np.rad2deg(np.arctan2(x, y))
#     data1 = np.loadtxt(files[0], usecols=np.arange(1,11))
#     c = 299792.458
#     ##########################################################################
#     # Account for difference in resolution
#     # Not used anymore because the resolution is now matched in pPXF
#     # fwhm_dif = (2.5 - 2.1) * c / 5500. / 2.3548
#     # data1[:,2] = np.sqrt(data1[:,2]**2 - fwhm_dif**2)
#     ##########################################################################
#     # Loading files
#     data2 = np.loadtxt(files[1], usecols=np.arange(1,26))
#     data3 = np.loadtxt(files[2], usecols=(1,2,3,5,6,7,9,10,11))
#     data4 = np.loadtxt(files[3], usecols=np.arange(1,26))
#     data5 = np.loadtxt(files[4], usecols=(1,))
#     data6 = np.loadtxt(files[5], usecols=(1,))
#     ##########################################################################
#     # Homogenization of the data
#     data1 = match_data(s1, sref, data1)
#     data2 = match_data(s2, sref, data2)
#     data3 = match_data(s3, sref, data3)
#     data4 = match_data(s4, sref, data4)
#     data5 = match_data(s5, sref, data5)
#     data6 = match_data(s6, sref, data6)
#     ##########################################################################
#     # Calculating composite indices: <Fe>, [MgFe]' and Mg b / <Fe>
#     data24 = np.zeros_like(np.column_stack((data2, data4)))
#     for i in np.arange(25):
#         data24[:, 2*i] = data2[:,i]
#         data24[:, 2*i+1] = data4[:,i]
#     fe5270 = data2[:,17]
#     fe5270_e = data4[:,17]
#     fe5335 = data2[:,18]
#     fe5335_e = data4[:,18]
#     mgb = data2[:,16]
#     mgb_e = data4[:,16]
#     meanfe = 0.5 * (fe5270 + fe5335)
#     meanfeerr = 0.5 * np.sqrt(fe5270_e**2 + fe5335_e**2)
#     term = (0.72 * fe5270 + 0.28 * fe5335)
#     mgfeprime = np.sqrt(mgb *  term)
#     mgfeprimeerr = 0.5 * np.sqrt(term /  mgb * (mgb_e**2) +
#     mgb / term * ((0.72 * fe5270_e)**2 + (0.28 * fe5335_e)**2))
#     sref = np.array(sref)
#     mgb_meanfe = mgb/meanfe
#     mgb_meanfe[mgb_meanfe>10] = np.nan
#     mgb_meanfe_err = np.sqrt((mgb * meanfeerr / meanfe**2)**2 +
#                               (mgb_e/meanfe)**2)
#     ##########################################################################
#     # Calculating [Fe / H]
#     ##########################################################################
#     feh = np.zeros((data3.shape[0], 3))
#     feh[:,0] = data3[:,3] - 0.94 * data3[:,6]
#     feh[:,1] = data3[:,4] -0.94 * data3[:,8]
#     feh[:,2] = data3[:,5] - 0.94 * data3[:,7]
#     ##########################################################################
#     # Saving results
#     results = np.column_stack((sref, x, y, r, pa, data1, data24, meanfe,
#               meanfeerr, mgfeprime, mgfeprimeerr, data3, coords,
#               mgb_meanfe, mgb_meanfe_err, data5, data6, feh))
#     header = ['FILE', "X[kpc]", "Y[kpc]", "R[kpc]", "PA", 'V', 'dV', 'S', 'dS',
#               'h3', 'dh3', 'h4', 'dh4',  'chi/DOF', 'S/N', 'Hd_A', 'dHd_A',
#               'Hd_F', 'dHd_F', 'CN_1', 'dCN_1', 'CN_2', 'dCN_2', 'Ca4227',
#               'dCa4227', 'G4300', 'dG4300', 'Hg_A', 'dHg_A', 'Hg_F', 'dHg_F',
#               'Fe4383', 'dFe4383', 'Ca4455', 'dCa4455', 'Fe4531',  'dFe4531',
#               'C4668', 'dC4668', 'H_beta', 'dH_beta', 'Fe5015', 'dFe5015',
#               'Mg_1', 'dMg_1', 'Mg_2',  'dMg_2', 'Mg_b',  'dMg_b', 'Fe5270',
#               'dFe5270', 'Fe5335', 'dFe5335', 'Fe5406',  'dFe5406', 'Fe5709',
#               'dFe5709', 'Fe5782', 'dFe5782', 'Na_D', 'dNa_D', 'TiO_1',
#               'dTiO_1', 'TiO_2',  'dTiO_2', "<Fe>", "d<Fe>", "[MgFe]'",
#               "d[MgFe]'", 'Age(Gyr)', 'Age-', 'Age+', '[Z/H]', '[Z/H]-',
#               '[Z/H]+', '[alpha/Fe]', '[alpha/Fe]-', '[alpha/Fe]+', "RA",
#               "DEC", "Mg b / <Fe>", "d Mg b / <Fe>",
#               "V-band surface brightness (mag arcsec-2)",
#               "Residual V-band surface brightness (mag arcsec-2)",
#               "[Fe / H]", "[Fe / H] lower limit", "[Fe / H] upper limit"]
#     with open(outtable, "w") as f:
#         for i,field in enumerate(header):
#             print "# {0} : {1}\n".format(i, field)
#             f.write("# {0} : {1}\n".format(i, field))
#         np.savetxt(f, results, fmt="%s")
    return

if __name__ == "__main__":
    merge_tables()
