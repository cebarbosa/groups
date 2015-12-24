#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on 24/09/14 13:38

@author: Carlos Eduardo Barbosa

Program to load templates from the MILES library for pPXF use.

"""
import os

import numpy as np
import pyfits as pf

import ppxf_util as util

from config import *

def load_templates_regul(velscale):
    """ Load templates into 2D array for regularization"""
    # cd into template folder to make calculations
    current_dir = os.getcwd()
    os.chdir(os.path.join(home, "miles_models"))
    # Extract the wavelength range and logarithmically rebin one spectrum
    # to the same velocity scale of the SAURON galaxy spectrum, to determine
    # the size needed for the array which will contain the template spectra.
    miles = [x for x in os.listdir(".") if x.endswith(".fits")]
    hdu = pf.open(miles[0])
    ssp = hdu[0].data
    h2 = hdu[0].header
    lamRange2 = h2['CRVAL1'] + np.array([0.,h2['CDELT1']*(h2['NAXIS1']-1)])

    sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp,
                                               velscale=velscale)
    # Ordered array of metallicities
    Zs = set([x.split("Z")[1].split("T")[0] for x in miles])
    Zs = [float(x.replace("m", "-").replace("p", "")) for x in Zs]
    Zs.sort()
    Z2 = Zs[:]
    for i in range(len(Zs)):
        if Zs[i] < 0:
            Zs[i] = "{0:.2f}".format(Zs[i]).replace("-", "m")
        else:
            Zs[i] = "p{0:.2f}".format(Zs[i])
    # Ordered array of ages
    Ts = list(set([x.split("T")[1].split(".fits")[0] for x in miles]))
    Ts.sort()
    # Create a three dimensional array to store the
    # two dimensional grid of model spectra
    #
    nAges = len(Ts)
    nMetal = len(Zs)
    templates = np.empty((sspNew.size,nAges,nMetal))

    # Here we make sure the spectra are sorted in both [M/H]
    # and Age along the two axes of the rectangular grid of templates.
    # A simple alphabetical ordering of Vazdekis's naming convention
    # does not sort the files by [M/H], so we do it explicitly below
    miles = []
    for k in range(nMetal):
        for j in range(nAges):
            filename = "Mun1.30Z{0}T{1}.fits".format(Zs[k], Ts[j])
            ssp = pf.getdata(filename)
            sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp,
                                                       velscale=velscale)
            templates[:,j,k] = sspNew # Templates are *not* normalized here
            miles.append(filename)
    templates /= np.median(templates) # Normalizes templates by a scalar
    os.chdir(current_dir)
    return templates, logLam2, Ts, Z2, miles, h2['CDELT1']

def stellar_templates(velscale):
    """ Load files with stellar library used as templates. """
    current_dir = os.getcwd()
    # Template directory is also set in config.py
    os.chdir(templates_dir)
    miles = [x for x in os.listdir(".") if x.startswith("Mun") and
             x.endswith(".fits")]
    # Ordered array of metallicities
    Zs = set([x.split("Z")[1].split("T")[0] for x in miles])
    Zs = [float(x.replace("m", "-").replace("p", "").replace("_", "."))
          for x in Zs]
    Zs.sort()
    for i in range(len(Zs)):
        if Zs[i] < 0:
            Zs[i] = "{0:.2f}".format(Zs[i]).replace("-", "m")
        else:
            Zs[i] = "p{0:.2f}".format(Zs[i])
    # Ordered array of ages
    Ts = list(set([x.split("T")[1][:7]
                   for x in miles]))
    Ts.sort()
    miles = []
    metal_ages = []
    for m in Zs:
        for t in Ts:
            filename = "Mun1.30Z{0}T{1}_linear_FWHM_3.6.fits".format(m, t)
            if os.path.exists(filename):
                miles.append(filename)
                metal_ages.append([m.replace("_", ".").replace("p",
                       "+").replace("m", "-"),t.replace("_", ".")])
    hdu = pf.open(miles[0])
    ssp = hdu[0].data
    h2 = hdu[0].header
    lamRange2 = h2['CRVAL1'] + np.array([0.,h2['CDELT1']*(h2['NAXIS1']-1)])
    sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp,
                                               velscale=velscale)
    templates = np.empty((sspNew.size,len(miles)))
    for j in range(len(miles)):
        hdu = pf.open(miles[j])
        ssp = hdu[0].data
        sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp,
                                                   velscale=velscale)
        templates[:,j] = sspNew
    os.chdir(current_dir)
    return templates, logLam2, h2['CDELT1'], miles

def emission_templates(velscale):
    """ Load files with stellar library used as templates. """
    current_dir = os.getcwd()
    # Template directory is also set in setyp.py
    os.chdir(templates_dir)
    emission = [x for x in os.listdir(".") if x.startswith("emission") and
             x.endswith(".fits")]
    emission.sort()
    c = 299792.458
    FWHM_tem = 2.54 # MILES library spectra have a resolution FWHM of 2.54A.
    # Extract the wavelength range and logarithmically rebin one spectrum
    # to the same velocity scale of the SAURON galaxy spectrum, to determine
    # the size needed for the array which will contain the template spectra.
    #
    hdu = pf.open(emission[0])
    ssp = hdu[0].data
    h2 = hdu[0].header
    lamRange2 = h2['CRVAL1'] + np.array([0.,h2['CDELT1']*(h2['NAXIS1']-1)])
    sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp,
                                               velscale=velscale)
    templates = np.empty((sspNew.size,len(emission)))
    for j in range(len(emission)):
        hdu = pf.open(emission[j])
        ssp = hdu[0].data
        sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp,
                                                   velscale=velscale)
        templates[:,j] = sspNew
    # templates *= 1e5 # Normalize templates
    os.chdir(current_dir)
    return templates, logLam2, h2['CDELT1'], emission

def emission_line_template(lines, velscale, res=3.6, intens=None, resamp=15,
                           return_log=True):
    lines = np.atleast_1d(lines)
    if intens == None:
        intens = np.ones_like(lines) * 1e-5
    current_dir = os.getcwd()
    os.chdir(templates_dir)
    refspec = [x for x in os.listdir(".") if x.endswith(".fits") and
               x.startswith("Mun")][0]
    lamb = wavelength_array(refspec)
    delta = lamb[1] - lamb[0]
    lamb2 = np.linspace(lamb[0]-delta/2., lamb[-1] + delta/2., len(lamb+1)*resamp)
    sigma = res / (2. * np.sqrt(2. * np.log(2.)))
    spec = np.zeros_like(lamb2)
    for line, intensity in zip(lines, intens):
        spec += intensity * np.exp(- (lamb2 - line)**2 / (2 * sigma * sigma))
    spec = np.sum(spec.reshape(len(lamb), resamp), axis=1)
    if not return_log:
        return spec
    specNew, logLam2, velscale = util.log_rebin([lamb[0], lamb[-1]], spec,
                                                   velscale=velscale)
    os.chdir(current_dir)
    return specNew

def make_fits(spec, outfile):
    hdu = pf.PrimaryHDU(spec)
    miles = [x for x in os.listdir(".") if x.startswith("Mun") and
             x.endswith(".fits")][0]
    w0 = pf.getval(miles, "CRVAL1")
    deltaw = pf.getval(miles, "CDELT1")
    pix0 = pf.getval(miles, "CRPIX1")
    hdu.header["CRVAL1"] = w0
    hdu.header["CDELT1"] = deltaw
    hdu.header["CRPIX1"] = pix0
    pf.writeto(outfile, hdu.data, hdu.header, clobber=True)
    return

def wavelength_array(spec):
    """ Produces array for wavelenght of a given array. """
    w0 = pf.getval(spec, "CRVAL1")
    deltaw = pf.getval(spec, "CDELT1")
    pix0 = pf.getval(spec, "CRPIX1")
    npix = pf.getval(spec, "NAXIS1")
    return w0 + deltaw * (np.arange(npix) + 1 - pix0)

if __name__ == "__main__":
    os.chdir(templates_dir)
    templates, logLam2, h, miles = stellar_templates(velscale)
    # em_OIII = emission_line_template([5006.84, 4958.91], velscale,
    #                                   intens=[1,0.33], return_log=0)
    # em_hbeta = emission_line_template(4861.333, velscale, return_log=0)
    # em_NI = emission_line_template(5200.257, velscale, return_log=0)
    # make_fits(em_OIII, "emission_OIII_fwhm3.6.fits")
    # make_fits(em_hbeta, "emission_hbeta_fwhm3.6.fits")
    # make_fits(em_NI, "emission_NI_fwhm3.6.fits")

