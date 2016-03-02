# -*- coding: utf-8 -*-
"""
Created on Mon March 2 2016

@author: cbarbosa

Calculate the Lick indices in the groups data
"""

import os

import numpy as np
import pyfits as pf
from scipy.interpolate import NearestNDInterpolator as interpolator
from scipy.interpolate import interp1d
from scipy.ndimage.filters import convolve1d
import matplotlib.pyplot as plt

from config import *
import lector as lector
from run_ppxf import pPXF, ppload


def wavelength_array(spec, axis=3, extension=0):
    """ Produces array for wavelenght of a given array. """
    w0 = pf.getval(spec, "CRVAL{0}".format(axis), extension)
    deltaw = pf.getval(spec, "CD{0}_{0}".format(axis), extension)
    pix0 = pf.getval(spec, "CRPIX{0}".format(axis), extension)
    npix = pf.getval(spec, "NAXIS{0}".format(axis), extension)
    return w0 + deltaw * (np.arange(npix) + 1 - pix0)

def correct_indices(indices, inderr, indtempl, indtempl_b):
    """ Make corrections for the broadening in the spectra."""
    bands = os.path.join(tables_dir, "bands.txt")
    types = np.loadtxt(bands, usecols=(8,))
    corrected = np.zeros_like(indices)
    errors = np.zeros_like(indices)
    for i,t in enumerate(types):
        if t == 0:
            C = indtempl[i] / indtempl_b[i]
            if C >= 1:
                corrected[i] = indices[i] * C
                errors[i] = C * inderr[i]
            else:
                corrected[i] = indices[i] 
                errors[i] = inderr[i]
        elif t in [1, 2]:
            C = indtempl[i] - indtempl_b[i]
            corrected[i] = indices[i] + C
            errors[i] = inderr[i]
    return corrected, errors

def check_intervals(setupfile, bands, vel):
    """ Check which indices are defined in the spectrum. """
    c = 299792.458 # speed of light in km/s
    with open(setupfile) as f:
        lines = [x for x in f.readlines()]
    lines = [x for x in lines if x.strip()]
    intervals = np.array(lines[5:]).astype(float)
    intervals = intervals.reshape((len(intervals)/2, 2))
    bands = np.loadtxt(bands, usecols=(2,7))
    bands *= np.sqrt((1 + vel/c)/(1 - vel/c))
    goodbands = np.zeros(len(bands))
    for i, (b1, b2) in enumerate(bands):
        for (i1, i2) in intervals:
            if i1 < b1 and b2 < i2:
                goodbands[i] = 1
    return np.where(goodbands == 1, 1, np.nan)

class BroadCorr:
    """ Wrapper for the interpolated model."""
    def __init__(self, table): 
        self.interpolate(table)
    
    def interpolate(self, table):
        sigmas = np.loadtxt(table, dtype=np.double, usecols=(0,))
        inds = np.loadtxt(table, dtype=np.double, usecols=np.arange(1,51, 2)).T
        corrs = np.loadtxt(table, dtype=np.double,usecols=np.arange(2,51, 2)).T
        self.fs = []
        self.indlims = []
        for i, (idx, corr) in enumerate(zip(inds, corrs)):
            f = interpolator(np.column_stack((sigmas, idx)), corr)
            self.fs.append(f)
            self.indlims.append([idx.min(), idx.max()])
        self.siglims = [sigmas.min(), sigmas.max()]
        return
        
    def __call__(self, sigma, lick): 
        b = np.zeros(len(lick))
        for i,l in enumerate(lick):
            if np.isnan(l):
                b[i] = 0.
            else:
                b[i] = self.fs[i](sigma, l)
        return b

class Vdisp_corr_k04():
    """ Correction for LOSVD only for multiplicative indices from
        Kuntschner 2004."""
    def __init__(self):
        """ Load tables. """
        table = os.path.join(tables_dir, "kuntschner2004.tab")
        bands = os.path.join(tables_dir, "BANDS")
        self.indices_k04 = np.loadtxt(table, usecols=(1,), dtype=str).tolist()
        self.type_k04 = np.loadtxt(table, usecols=(2,),
                                            dtype=str).tolist()
        self.coeff_k04 = np.loadtxt(table, usecols=np.arange(3,10))
        self.lick_indices = np.loadtxt(bands, usecols=(0,), dtype=str)
        self.lick_types = np.loadtxt(bands, usecols=(8,))

    def __call__(self, lick, sigma, h3=0., h4=0.):
        newlick = np.zeros(25)
        for i,index in enumerate(lick_indices):
            if index in self.indices_k04:
                idx = self.indices_k04.index(index)
                a1, a2, a3, b1, b2, c1, c2 = self.coeff_k04[idx]
                if self.type_k04[idx] == "m":
                    C_k04 = 1. + a1 * sigma + a2 * sigma**2 + \
                            a3 * sigma**3 + b1 * sigma * h3 + \
                            b2 * sigma**2 * h3 + \
                            c1 * sigma * h4 + c2 * sigma**2 * h4
                    newlick[i] = C_k04 * lick[i]
                else:
                    C_k04 = a1 * sigma    + a2 * sigma**2 + \
                            a3 * sigma**3 + b1 * sigma * h3 + \
                            b2 * sigma**2 * h3 + \
                            c1 * sigma * h4 + c2 * sigma**2 * h4
                    newlick[i] = lick[i] + C_k04
            else:
                newlick[i] = lick[i]
        return newlick


def losvd_convolve(spec, losvd):
    """ Apply LOSVD to a given spectra given that both wavelength and spec
     arrays are log-binned. """
    # Convert to pixel scale
    pars = np.copy(losvd)
    pars[:2] /= velscale
    dx = int(np.ceil(np.max(abs(pars[0]) + 5*pars[1])))
    nl = 2*dx + 1
    x = np.linspace(-dx, dx, nl)   # Evaluate the Gaussian using steps of 1/factor pixel
    vel = pars[0]
    w = (x - vel)/pars[1]
    w2 = w**2
    gauss = np.exp(-0.5*w2)
    profile = gauss/gauss.sum()
    # Hermite polynomials normalized as in Appendix A of van der Marel & Franx (1993).
    # Coefficients for h5, h6 are given e.g. in Appendix C of Cappellari et al. (2002)
    if losvd.size > 2:        # h_3 h_4
        poly = 1 + pars[2]/np.sqrt(3)*(w*(2*w2-3)) \
                 + pars[3]/np.sqrt(24)*(w2*(4*w2-12)+3)
        if len(losvd) == 6:  # h_5 h_6
            poly += pars[4]/np.sqrt(60)*(w*(w2*(4*w2-20)+15)) \
                  + pars[5]/np.sqrt(720)*(w2*(w2*(8*w2-60)+90)-15)
        profile *= poly
    return convolve1d(spec, profile)

def run(nights, specnames):
    for night, spectrum in zip(nights, specnames):
        print "Working with object {0}, {1}".format(spectrum, night)
        os.chdir(os.path.join(data_dir, night))
        spec = pf.getdata(spectrum)
        w = wavelength_array(spectrum, axis=1, extension=0)
        ssps_file = os.path.join(templates_dir, 'miles_FWHM_3.6.fits')
        ssps = pf.getdata(ssps_file, 0)
        wssps = np.power(np.e, pf.getdata(ssps_file, 1))
        pp = ppload("logs/{0}".format(spectrum.replace(".fits", "")))
        pp = pPXF(spectrum, velscale, pp)
        ######################################################################
        # Make best sky interpolation for subtraction
        sky_interp = interp1d(pp.w, pp.bestsky, kind="linear",
                              bounds_error=False, fill_value=0.)
        sky = sky_interp(w)
        spec -= sky
        ######################################################################
        # Make interpolation of polynomial
        poly = np.interp(wssps, pp.w, pp.poly)
        ######################################################################
        # Make the combination
        ssps_unbroad_v0 = ssps.dot(pp.w_ssps)
        ssps_broad = losvd_convolve(ssps_unbroad_v0, pp.sol) + poly
        sol_unb = np.zeros_like(pp.sol)
        sol_unb[0] = pp.sol[0]
        sol_unb[1] = 0.1 * velscale
        ssps_unbroad = losvd_convolve(ssps_unbroad_v0, np.array(sol_unb)) + poly
        ######################################################################
        # Correct for emission line if present
        em = interp1d(pp.w, pp.gas, kind="linear", bounds_error=False,
                      fill_value=0.)
        emission = em(w)
        spec -= emission
        ######################################################################
        noise =  np.ones_like(spec) * np.nanstd(pp.galaxy - pp.bestfit)
        ######################################################################
        # Broaden to Lick/IDS system
        spec = lector.broad2lick(w, spec, 3.6, vel=pp.sol[0])
        ssps_broad = lector.broad2lick(wssps, ssps_broad, 3.6, vel=pp.sol[0])
        ssps_unbroad = lector.broad2lick(wssps, ssps_unbroad, 3.6, vel=pp.sol[0])
        # plt.plot(w, spec, "-k")
        # plt.plot(wssps, ssps_broad , "-r")
        # plt.plot(wssps, ssps_unbroad , "-b")
        # plt.show()
        lick, lickerr = lector.lector(w, spec, noise, bands, vel=pp.sol[0])
        noise_temp = np.ones_like(wssps) * noise[0]
        lick_unb, tmp = lector.lector(wssps, ssps_unbroad, noise_temp, bands,
                                 vel=pp.sol[0])
        lick_br, tmp = lector.lector(wssps, ssps_broad, noise_temp, bands,
                     vel=pp.sol[0])
        ##################################################################
        # LOSVD correction using best fit templates
        ##################################################################
        lickc, lick_cerr = correct_indices(lick, lickerr, lick_unb,
                                           lick_br)
    #     print w[0], w[-1]
    #     print lickc
    #     raw_input()
    # #     ################################################################
    # #     # Convert to string
    # #     ################################################################
    #     lick = "".join(["{0:14}".format("{0:.5f}".format(x)) for x in lick])
    #     lickc = "".join(["{0:14}".format("{0:.5f}".format(x)) for x in lickc])
    # #     # Append to output
    #     logfile1 = "logs_sn{0}/lick_{1}_bin{2:04d}_raw.txt".format(targetSN,
    #                                                            field, bin)
    #     logfile2 = "logs_sn{0}/lick_{1}_bin{2:04d}_corr.txt".format(targetSN,
    #                                                             field, bin)
    #     with open(logfile1, "w") as f:
    #         f.write("{0:16s}".format(name) + lick + "\n")
    #     with open(logfile2, "w") as f:
    #         f.write("{0:16s}".format(name) + lickc + "\n")
    #     results.append("{0:16s}".format(name) + lick)
    #     resultsc.append("{0:16s}".format(name) + lickc)

def make_table(fields, targetSN, ltype="corr"):
    """ Gather information of Lick indices of a given target SN in a single
    table. """
    for field in fields:
        os.chdir(os.path.join(data_dir, "combined_{0}".format(field),
                              "logs_sn{0}".format(targetSN)))
        logfiles = sorted([x for x in os.listdir(".") if x.startswith("lick")
                           and x.endswith("{0}.txt".format(ltype))])
        table = []
        for logfile in logfiles:
            with open(logfile) as f:
                line = f.readline()
            table.append(line)
        output = os.path.join(data_dir, "combined_{0}".format(field),
                              "lick_{0}_sn{1}.txt".format(ltype, targetSN))
        with open(output, "w") as f:
            f.write("".join(table))

def lick_standards_table():
    """ Prepare a table with the Lick indices of the standard stars"""
    folder = os.path.join(home, "data/standards")
    fits = [x for x in os.listdir(folder) if x.endswith(".fits")]
    stars = [x.split(".")[0].lower() for x in fits]
    stars = [x.upper() for x in stars if x.startswith("hr") or
             x.startswith("hd")]
    stars = list(set([x.replace("HD0", "HD") for x in stars]))
    table = os.path.join(tables_dir, "lick_standards.dat")
    with open(table) as f:
        header = f.readline()
        data = f.readlines()
    cols = np.arange(39,240,8)
    lick = np.zeros((len(data), 25))
    for i in range(25):
        d = [x[cols[i]:cols[i+1]].strip() for x in data]
        d = np.array([x if x else np.nan for x in d])
        lick[:, i] = d
    hd = ["HD" + x[:9].replace(" ", "") for x in data]
    hd = [x.replace("HD0", "HD") for x in hd]
    hr = [x[9:19].replace(" ", "") for x in data]
    table = []
    for star in stars:
        if star in hr:
            idx = hr.index(star)
        elif star in hd:
            idx = hd.index(star)
        else:
            continue
        l = ["{0:.5f}".format(x) for x in lick[idx]]
        l = ["{0:10s}".format(x) for x in l]
        table.append("{0:12s}".format(star) + "".join(l))
    with open(os.path.join(tables_dir, "lick_standards.txt"), "w") as f:
        f.write("\n".join(table))

if __name__ == "__main__":
    bands = os.path.join(tables_dir, "bands.txt")
    lick_standards_table()
    # specs, nights = np.loadtxt(os.path.join(tables_dir, "spec_location.txt"),
    #                            dtype=str, usecols=(0,1)).T
    # run(nights[3:], specs[3:])