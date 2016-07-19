# -*- coding: utf-8 -*-
"""

Created on 25/05/16

@author: Carlos Eduardo Barbosa

Make calculations of the dynamics of groups

"""

import os

import numpy as np
from astropy.io import ascii
import astropy.coordinates as coord
import astropy.units as units
from astropy import constants
import astropy.stats.funcs as funcs
from astropy.stats import sigma_clip
from astropy.cosmology import Planck15 as cosmo
import matplotlib.pyplot as plt

from config import *
from make_tables import make_str_err

def get_zm98():
    """ Add galaxies from the catalog of Zabludoff and Mulchaey 1998."""
    path = os.path.join(tables_dir, "zm98")
    table = ascii.read(os.path.join(path, 'table1.dat'),
                       readme=os.path.join(path, "ReadMe"))
    ra = coord.Angle(table["RAh"] + table["RAm"]/60. + table["RAs"]/3600.,
                     unit=units.hour)
    signal = np.ones_like(ra.degree)
    signal[table["DE-"]=="-"] = -1
    dec = coord.Angle(signal * (table["DEd"] + table["DEm"]/60. +
                                table["DEs"]/3600.), unit=units.degree)
    return (table["Name"], ra, dec, table["HRV"],
            table["e_HRV"])

def group_clip(data, v0):
    """ Remove galaxies that are not members. """
    idx = np.where(np.abs(data[3] - v0) < 1500)[0]
    d = tuple([x[idx] for x in data])
    newdata = sigma_clip(d[3], sigma=3, cenfunc=np.median,
                               stdfunc=funcs.biweight_midvariance,
                               copy=True, iters=10)
    mask = ~newdata._mask
    return tuple([x[mask] for x in d])

def biweight_clipped(v):
    return  funcs.biweight_midvariance(sigma_clip(v, sigma=4, cenfunc=np.mean,
                               stdfunc=funcs.biweight_midvariance,
                               copy=False, iters=5))

def gapper(v):
    n = len(v)
    gi = np.diff(np.sort(v))
    i = np.arange(1, n)
    wi = i * (n - i)
    return np.sqrt(np.pi) / (n * (n-1)) * np.sum(wi * gi)

def virial_mass(v, r):
    """ Calculate the virial mass according to Heisler, Tremaine and Bahcall
    (1985)"""
    vcentroid = np.mean(v)
    N = len(r)
    term1 = np.power(v - vcentroid, 2).sum()
    ii, jj = np.triu_indices(N, k=1)
    term2 = 0
    for i,j in zip(ii,jj):
        term2 += 1 / np.abs(r[i]-r[j])
    return 3. * np.pi * N / (2. * G) * (term1 * units.km**2 / units.s**2)/ \
           (term2 * units.kpc**-1)

def harmonic_projected_radius(ra, dec, V):
    """ Virial and projected radius according to Ramella+ 1989"""
    N = len(ra)
    ra1, ra2 = np.meshgrid(ra.degree, ra.degree)
    dec1, dec2 =  np.meshgrid(dec.degree, dec.degree)
    c1 = coord.SkyCoord(ra=ra1.flatten()*units.degree,
                        dec=dec1.flatten()*units.degree)
    c2 = coord.SkyCoord(ra=ra2.flatten()*units.degree,
                        dec=dec2.flatten()*units.degree)
    theta = np.reshape(c1.separation(c2).radian, (N,N))
    Rh = np.pi * V / cosmo.H0 * \
         np.sin(0.5 * (0.5 * N * (N-1)) * np.sum(np.triu(1/theta, k=1))**-1)
    Rp = 8. * V / np.pi / cosmo.H0 * \
        np.sin(np.sum(np.triu(theta, k=1))/(N * (N-1)))
    return Rh, Rp

def crossing_time(Rv, sigma):
    """ Crossing time according to Ramella+ 1989"""
    return 3. / np.power(5., 1.5) * Rv * cosmo.H0 / sigma

def virial_radius_mass(ra, dec, D, sigma, alpha=2.6):
    """ Virial radius and mass according to definitions in Tully 2015."""
    N = len(ra)
    ra1, ra2 = np.meshgrid(ra, ra)
    dec1, dec2 =  np.meshgrid(dec, dec)
    c2 = coord.SkyCoord(ra=ra1.flatten(), dec=dec1.flatten(),
                        unit=units.degree)
    c3 = coord.SkyCoord(ra=ra2.flatten(), dec=dec2.flatten(),
                        unit=units.degree)
    sep = D * np.tan(c2.separation(c3).radian.reshape((N,N)))
    ij = np.triu_indices(N,k=1)
    Rg = N * N / np.sum(1. / sep[ij])
    Mv = alpha * np.pi / 2. * sigma**2 * Rg / G
    return Rg, Mv


def calc_zm98():
    """ Calculate dynamics for data in Z&M98. """
    print "Results from ZM98."
    data = get_zm98()
    objs = np.array([x.split("_")[0] for x in data[0]])
    dist = {"H42":56, "H62":63.6, "H90":39.5, "N4325":112.9, 'N491':56.5,
                 'N5129':101.5, 'N533':82.9, 'N5846':28.0, 'N664':75.1,
            "N2563":64.5, "N741":83., "N7582":22.9}
    table = os.path.join(tables_dir, "zm98/table2.dat")
    readme = os.path.join(tables_dir, "zm98/ReadMe")
    table2 = ascii.read(table, readme=readme)
    sample = [x.replace(" ", "").replace("GC","").replace("CG","") for
              x in table2["Group"]]
    results = []
    for j, group in enumerate(sample):
        print group,
        idx = np.array([i for i,x in enumerate(objs) if x.startswith(group)])
        zmdata = []
        for d in data:
            zmdata.append(d[idx])
        D = dist[group] * units.Mpc
        v0 = cosmo.H0 * D
        zmgals = group_clip(zmdata, v0)
        N = len(zmgals[0])
        ra = zmgals[1]
        dec = zmgals[2]
        print biweight_clipped(zmgals[3]), np.std(zmgals[3])
        continue
        line = dynamics(group, ra, dec, zmgals[3], zmgals[4], D)
        results.append(line)
    results = np.array(results)
    # with open(os.path.join(tables_dir, "zm98_biweight.txt"), "w") as f:
    #     np.savetxt(f, results, fmt="%s")
    return

def calc_our():
    """ Calculate dynamics for our groups. """
    table = os.path.join(data_dir, "results_unique.tab")
    specs = np.loadtxt(table, usecols=(0,), dtype=str)
    vs, verrs = np.loadtxt(table, usecols=(1,2)).T
    ras, decs = np.loadtxt(table, usecols=(72,73), dtype=str).T
    print "Results from our work."
    results = []
    for group in groups[1:]:
        print group
        D = distances[group]* units.Mpc

        idx = np.array([i for i,spec in enumerate(specs) if
                        spec.startswith(group)])
        ra = coord.Angle(ras[idx], unit=units.hourangle)
        dec = coord.Angle(decs[idx], unit=units.degree)
        line = dynamics(group, ra, dec, vs[idx], verrs[idx], D)
        results.append(line)
    results = np.array(results)
    with open(os.path.join(data_dir, "virial.txt"), "w") as f:
        np.savetxt(f, results, fmt="%s")
    return

def bootstrap(v, verr, f, nsim=1000):
    """ Make bootstrapping calculations."""
    results = np.zeros(nsim)
    for i in range(nsim):
        vsim = v + verr * np.random.normal(size=len(verr))
        vboot = np.random.choice(vsim, size=vsim.size, replace=True)
        results[i] = f(vboot)
    return np.nanstd(results)

def plot_zm98():
    """ Comparison with ZM98. """
    table = os.path.join(tables_dir, "zm98/table2.dat")
    readme = os.path.join(tables_dir, "zm98/ReadMe")
    table2 = ascii.read(table, readme=readme)
    table = os.path.join(tables_dir, "zm98_test.txt")
    data = np.loadtxt(table, usecols=(0,2,4,6,8,10)).T
    err =  np.loadtxt(table, usecols=(1,3,5,7,9,11)).T
    zmfields = ["HRV", "sigma", "rp", "rh", "Mvir", "tc/tH"]
    zmferrs = ['e_HRV', 'e_sigma', None, None, None, None]
    lims = [[1000, 8000], [0, 500], [0.7, 1.2], [0.4,  0.9], [0, 2.5],
            [0, 0.1]]
    div = [1., 1., 0.65 , 0.65, 0.65, 1.]
    labels = ["$v$ (km/s)", "$\sigma_v$ (km/s)",
              "$R_p$ (Mpc)", "$R_h$ (Mpc)",
              "$M_v$ ($10^{14}$M$_\odot$)", "$t_c/t_H$"]
    fig = plt.figure(1, figsize=(9,5.6))
    scale = [1., 1., 1., 1., 1e-14, 1.]
    panels = ["A", "B", "C", "D", "E", "F"]
    for i in range(6):
        ax = plt.subplot(2,3,i+1)
        ax.minorticks_on()
        err2 = None if zmferrs[i] is None else table2[zmferrs[i]]/div[i]
        ax.errorbar(table2[zmfields[i]]/div[i]*scale[i], data[i] * scale[i],
                    yerr=err[i] * scale[i], xerr=err2, ecolor="0.8",
                    fmt="o", color="r",
                    mec="r", ms=8)
        ax.set_xlim(lims[i])
        ax.set_ylim(lims[i])
        plt.plot(lims[i], lims[i], "--k")
        ax.set_ylabel(labels[i] + " - this work")
        ax.set_xlabel(labels[i] + " - ZM98")
        ax.locator_params(nbins=5)
        ax.legend([], [], title="({0})".format(panels[i]), loc=2, frameon=False)
    plt.subplots_adjust(left=0.09, right=0.96, top=0.97, hspace=0.3,
                        wspace=0.35)
    plt.savefig(os.path.join(plots_dir, "zmcomp.png"))
    return

def make_latex():
    intables = [os.path.join(tables_dir, "zm98_biweight.txt"),
              os.path.join(data_dir, "virial.txt")]
    outtables = [x.replace(".txt", ".tex") for x in intables]
    for intable, outtable in zip(intables, outtables):
        out = np.loadtxt(intable, usecols=(0,1), dtype=str)
        # Kinematics table
        kcols = (2,4,6,8,10,12)
        div = [1, 1, 1, 1, 1e12,1]
        for i, col in enumerate(kcols):
            kdata = np.loadtxt(intable, usecols=(col,)) / div[i]
            kerr = np.loadtxt(intable, usecols=(col+1,)) / div[i]
            out = np.column_stack((out, make_str_err(kdata, kerr)))
        out = [" & ".join(x) + "\\\\" for x in out]
        with open(outtable, "w") as f:
            f.write("\n".join(out))

def dynamics(group, ra, dec, v, verr, D):
    sigma = biweight_clipped(v) * units.km / units.s
    sigmaerr = bootstrap(v, verr, biweight_clipped) * units.km / units.s
    N = len(ra)
    V = np.mean(v) * units.km / units.s
    Verr = bootstrap(v, verr, np.mean) * units.km / units.s
    Rg, Mv = virial_radius_mass(ra.degree, dec.degree, D, sigma, alpha=2.6)
    Rg_err = Rg / cosmo.H0.value * 5
    Mverr = np.abs(2 * Mv / sigma * sigmaerr)
    Rh, Rp = harmonic_projected_radius(ra, dec, V)
    Rh_err = Rh / V * Verr
    tc = crossing_time(Rh, sigma)
    tc_err = np.sqrt((tc/cosmo.H0*0.9*units.km / units.s / units.Mpc)**2 +
                     (tc*sigmaerr/sigma)**2 + (tc/Rg*Rg_err)**2)
    line = np.array([group, N, V, Verr, sigma.value, sigmaerr.value,
                     Rh.value, Rh_err.value,
                     Rg.value, Rg_err.value,
                     Mv.value, Mverr.value, tc.value, tc_err.value])
    return line


if __name__ == "__main__":
    G = constants.G.to("(km2 kpc)/(M_sun s2)")
    # calc_zm98()
    plot_zm98()
    # calc_our()
    # make_latex()
