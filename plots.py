# -*- coding: utf-8 -*-
"""

Created on Aug 29, 206

@author: Carlos Eduardo Barbosa

Creates various plots for the project.

"""
import os

import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
from scipy.stats import spearmanr

from config import *

def sigma_lick():
    """ Plot correlations between Lick indices and the central velocity
    dispersion """
    table = os.path.join(data_dir, "results_unique.tab")
    licktable = os.path.join(tables_dir, "bands.txt")
    names, units = np.loadtxt(licktable, usecols=(0,9), dtype=str).T
    w = np.loadtxt(licktable, usecols=(4,5)).T
    deltaw = np.diff(w, axis=0)[0]
    names = [x.replace("_", "").replace("beta", "$\\beta$") for x in names]
    cols = np.arange(14,63,2)
    emission = np.loadtxt(table, usecols=(13,), dtype=str)
    sigma = np.loadtxt(table, usecols=(3,4)).T
    lick = np.loadtxt(table, usecols=cols).T
    error = np.loadtxt(table, usecols=cols+1).T
    for j,unit in enumerate(units):
        if unit == "Ang":
            error[j] = np.abs(2.5 / np.log(10) * error[j] / (deltaw[j] - lick[j]))
            lick[j] = -2.5 * np.log10(1 - lick[j] / deltaw[j])
    idxem = emission == "yes"
    idxpas = ~idxem
    fig = plt.figure(1, figsize=(14,8))
    k00 = k00_scaling_relations()
    corr_none, corr_weak, corr_mod, corr_strg = [], [], [], []
    for i in range(25):
        std = np.std(sigma_clip(lick[i]))
        loc = np.nanmedian(lick[i])
        ax = plt.subplot(5,5,i+1)
        ax.minorticks_on()
        ax.errorbar(np.log10(sigma[0]), lick[i],
                    xerr=sigma[1] / sigma[0] / np.log(10),
                    yerr = error[i],
                    marker=None, ecolor="0.8", fmt="none")
        ax.plot(np.log10(sigma[0])[idxem], lick[i][idxem], "^b", ms=8)
        ax.plot(np.log10(sigma[0])[idxpas], lick[i][idxpas], "or", ms=8)
        ax.set_ylabel(names[i] + " (mag)")
        plt.ylim(loc - 3.5 * std,  loc + 5 * std)
        s, p = spearmanr(np.log10(sigma[0]), lick[i])
        if np.abs(s) <= 0.1:
            corr_none.append(names[i])
        elif s <= 0.3:
            corr_weak.append(names[i])
        elif s<= 0.5:
            corr_mod.append(names[i])
        else:
            corr_strg.append(names[i])
        if i < 20:
            ax.get_xaxis().set_ticklabels([])
        else:
            ax.set_xlabel("$\log \sigma_0$ (km/s)")
        ax.get_yaxis().set_label_coords(-0.22,0.5)
        if k00[i] is not None:
            ls0 = np.linspace(1.6, 2.6, 100)
            ax.plot(ls0, k00[i](ls0), "-g", lw=3)
        plt.xlim(1, 2.6)
        plt.legend([], [], title="r={0:.2f} p={1:.2f}".format(s,p), loc=2,
                   prop={'size':8})
    plt.subplots_adjust(left=0.05, right=0.99, top=0.98, bottom=0.08,
                        hspace=0.05, wspace=0.38)
    plt.savefig(os.path.join(plots_dir, "sigma_lick.png"))
    for c in [corr_none, corr_weak, corr_mod, corr_strg]:
        print len(c), c
    # plt.plot(np.log10(data[:,0])[idxpas], np.log10(data[:,1][idxpas]), "or")

    # plt.ylim(0,0.8)
    # plt.show()

def k00_scaling_relations():
    """ Scaling relations for Lick indices of galaxies in the Fornax cluster
    by Kuntschner 2000"""
    k00 = 25 * [None]
    k00[15] = lambda x: 0.191 * x -0.127
    k00[14] = lambda x: 0.136 * x -0.158
    k00[16] = lambda x: 0.102 * x -0.056
    k00[11] = lambda x: 0.090 * x - 0.11
    k00[8] = lambda x: 0.043 * x + 0.37
    k00[10] = lambda x: 0.036 * x + 0.009
    k00[13] = lambda x: 0.036 * x + 0.002
    k00[17] = lambda x: 0.029 * x + 0.024
    k00[18] = lambda x: 0.043 * x - 0.020
    k00[19] = lambda x: 0.023 * x + 0.023
    k00[9] = lambda x: 0.035 * x + 0.014
    k00[12] = lambda x: -0.020 * x + 0.106
    k00[6] = lambda x: -0.045 * x -0.038
    k00[7] = lambda x: -0.049 * x - 0.018
    return k00

def sigma_stpop():
    """ Plot correlations between Lick indices and the central velocity
    dispersion """
    table = os.path.join(data_dir, "results_unique.tab")
    cols = np.array((64,67,70,67))
    emission = np.loadtxt(table, usecols=(13,), dtype=str)
    sigma = np.loadtxt(table, usecols=(3,4)).T
    pop = np.loadtxt(table, usecols=cols).T
    error1 = np.loadtxt(table, usecols=cols+1).T
    error2 =  np.loadtxt(table, usecols=cols+2).T
    # Calculating Fe/H
    pop[3] -= 0.94 * pop[2]
    error1[3] = np.sqrt(error1[1]**2 + error1[2]**2)
    error2[3] = np.sqrt(error2[1]**2 + error2[2]**2)
    ###########################################################################
    idxem = emission == "yes"
    idxpas = ~idxem
    fig = plt.figure(1, figsize=(7,6))
    ylims = [[-1.5, 1], [-0.3, 0.5], [-1.5,1]]
    names = ["$\log $(Age) (Gyr)", "[Z/H]", r"[$\alpha$/Fe]", "[Fe/H]"]
    P04 = [None, [(1.7, 2.6), (-0.2,0.67)], None, None]
    for i in [1,2,3]:
        ax = plt.subplot(3,2,2*i-1)
        ax.minorticks_on()
        ax.errorbar(np.log10(sigma[0]), pop[i],
                    xerr=sigma[1] / sigma[0] / np.log(10),
                    yerr=[error1[i], error2[i]],
                    marker=None, ecolor="0.8", fmt="none")
        ax.plot(np.log10(sigma[0])[idxem], pop[i][idxem], "^b", ms=8)
        ax.plot(np.log10(sigma[0])[idxpas], pop[i][idxpas], "or", ms=8)
        ax.set_xlim(1., 2.6)
        ax.set_ylim(ylims[i-1])
        ax.set_ylabel(names[i], size=14)
        if i < 3:
            ax.get_xaxis().set_ticklabels([])
        else:
            ax.set_xlabel("$\log \sigma_0$ (km/s)", size=14)
        if P04[i] is not None:
            ax.plot(P04[i][0],P04[i][1], "-g", lw=3)
        # Correlations
        s, p = spearmanr(np.log10(sigma[0]), pop[i])
        plt.legend([], [], title="r={0:.2f} p={1:.2f}".format(s,p), loc=2,
                   prop={'size':8})
    P04 = [None, None, [(0.1,1.45),(-0.16,0.38)],  [(0.1,1.45),(-0.4,-0.35)]]
    for i in [1,2,3]:
        ax = plt.subplot(3,2,2*i)
        ax.minorticks_on()
        ax.errorbar(np.log10(pop[0]), pop[i],
                    xerr=[error1[0]/pop[0]/np.log(10), error2[0]]/pop[0]
                         /np.log(10),
                    yerr=[error1[i], error2[i]],
                    marker=None, ecolor="0.8", fmt="none")
        ax.plot(np.log10(pop[0])[idxem], pop[i][idxem], "^b", ms=8)
        ax.plot(np.log10(pop[0])[idxpas], pop[i][idxpas], "or", ms=8)
        ax.get_yaxis().set_ticklabels([])
        ax.set_xlim(0.5, 1.2)
        ax.set_ylim(ylims[i-1])
        if i < 3:
            ax.get_xaxis().set_ticklabels([])
        else:
            ax.set_xlabel(names[0], size=14)
        if P04[i] is not None:
            ax.plot(P04[i][0], P04[i][1], "-g", lw=3)
        # Correlations
        s, p = spearmanr(np.log10(sigma[0]), pop[i])
        plt.legend([], [], title="r={0:.2f} p={1:.2f}".format(s, p), loc=2,
                   prop={'size': 8})
    plt.subplots_adjust(left=0.11, right=0.97, top=0.98, bottom=0.10,
                        hspace=0.07, wspace=0.07)
    plt.savefig(os.path.join(plots_dir, "sigma_age_metallicities.png"))
    return

if __name__ == "__main__":
    # sigma_lick()
    sigma_stpop()