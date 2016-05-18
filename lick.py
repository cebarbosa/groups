# -*- coding: utf-8 -*-
"""

Created on 16/05/16

@author: Carlos Eduardo Barbosa

Program to calculate lick indices

"""

import numpy as np
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from scipy.ndimage.filters import gaussian_filter1d
from scipy.constants import c

from lector import lector

ckms = c / 1000. # Convert speed of light to km / s

class Lick():
    def __init__(self, wave, galaxy, bands, vel=0):
        self.galaxy = galaxy
        self.wave = wave
        self.dw = np.diff(self.wave)
        if not np.all(self.dw == self.dw[0]):
            raise "Dispersion not linear"
        self.dw = self.dw[0]
        self.bands = bands * np.sqrt((1 + vel/ckms)/(1 - vel/ckms))
        self.classic_integration()


    def classic_integration(self):
        k_order=10.
        # Define number of points for integration
        npoints = 2**k_order + 1
        self.R = np.zeros(self.bands.shape[0])
        self.Ia = np.zeros_like(self.R)
        self.Im = np.zeros_like(self.R)
        for i, w in enumerate(self.bands):
            if (w[0] - self.dw < self.wave[0]) or \
               (w[-1] + self.dw > self.wave[-1]):
                self.R[i] = np.nan
            # Defining indices for each section
            idxb = np.where(((self.wave > w[0] - 2 * self.dw) &
                                 (self.wave < w[1] + 2 * self.dw)))
            idxr = np.where(((self.wave > w[4] - 2 * self.dw) &
                                (self.wave < w[5] + 2 * self.dw)))
            idxcen = np.where(((self.wave > w[2] - 2 * self.dw) &
                                (self.wave < w[3] + 2 * self.dw)))
            # Defining wavelenght samples
            wb = self.wave[idxb]
            wr = self.wave[idxr]
            wcen = self.wave[idxcen]
            # Defining intensity samples
            fb = self.galaxy[idxb]
            fr = self.galaxy[idxr]
            fcen = self.galaxy[idxcen]
            # Making interpolation functions
            sb = InterpolatedUnivariateSpline(wb, fb)
            sr = InterpolatedUnivariateSpline(wr, fr)
            # Make oversampled arrays of wavelenghts
            xb = np.linspace(w[0], w[1], npoints)
            xr = np.linspace(w[4], w[5], npoints)
            xcen = np.linspace(w[2], w[3], npoints)
            # Calculating the mean fluxes for the pseudocontinuum
            fp1 = sb.integral(w[0], w[1]) / (w[1] - w[0])
            fp2 = sr.integral(w[4], w[5]) / (w[5] - w[4])
            # Making pseudocontinuum vector
            x0 = (w[2] + w[3])/2.
            x1 = (w[0] + w[1])/2.
            x2 = (w[4] + w[5])/2.
            fc = fp1 + (fp2 - fp1)/ (x2 - x1) * (wcen - x1)
            # Calculating indices
            ffc = InterpolatedUnivariateSpline(wcen, fcen/fc/(w[3]-w[2]))
            self.R[i] =  ffc.integral(w[2], w[3])
            self.Ia[i] = (1 - self.R[i]) * (w[3]-w[2])
            self.Im[i] = -2.5 * np.log10(self.R[i])
        self.classic = np.copy(self.Ia)
        idx = np.array([2,3,14,15,23,24])
        self.classic[idx] = self.Im[idx]
        return





def broad2lick(wl, intens, obsres, vel=0):
    """ Convolve spectra to match the Lick resolution.

    Broad a given spectra to the Lick system resolution. As the resolution
    in the Lick system varies as function of the wavelength, we use the
    interpolated values from Worthey and Ottaviani 1997.

    ================
    Input parameters
    ================
    wl: array_like
        Wavelenght 1-D array in Angstroms.

    intens: array_like
        Intensity 1-D array of Intensity, in arbitrary units. The lenght has
        to be the same as wl.

    obsres: float or array
        Value of the observed resolution Full Width at Half Maximum (FWHM) in
        Angstroms.

    vel: float
        Recession velocity of the measured spectrum.

    =================
    Output parameters
    =================
    array_like
        The convolved intensity 1-D array.

    """
    dw = wl[1] - wl[0]
    if not isinstance(obsres, np.ndarray):
        obsres = np.ones_like(wl) * obsres
    wlick = np.array([2000., 4000., 4400., 4900., 5400., 6000., 8000.]) * \
            np.sqrt((1 + vel/ckms)/(1 - vel/ckms))
    lickres = np.array([11.5, 11.5, 9.2, 8.4, 8.4, 9.8, 9.8])
    flick = interp1d(wlick, lickres, kind="linear", bounds_error=False,
                         fill_value="extrapolate")
    fwhm_lick = flick(wl)
    fwhm_broad = np.sqrt(fwhm_lick**2 - obsres**2)
    sigma_b = fwhm_broad/ 2.3548 / dw
    intens2D = np.diag(intens)
    for i in range(len(sigma_b)):
        intens2D[i] = gaussian_filter1d(intens2D[i], sigma_b[i],
                      mode="constant", cval=0.0)
    return intens2D.sum(axis=0)

def bands_shift(bands, vel):
    return  bands * np.sqrt((1 + vel/ckms)/(1 - vel/ckms))

def example():
    import os
    import pyfits as pf
    from config import data_dir, tables_dir
    from run_ppxf import wavelength_array
    os.chdir(data_dir)
    filename = "hcg22_023_blanco11an2.fits"
    v = 2.74319141e+03 # km/s

    bandsfile = os.path.join(tables_dir, "bands.txt")
    bands = np.loadtxt(bandsfile, usecols=np.arange(2,8))
    spec = pf.getdata(filename)
    noise = np.ones_like(spec)
    wave = wavelength_array(filename)
    lick, err = lector(wave, spec, noise, bandsfile, vel=v)
    print lick
    ll = Lick(wave, spec, bands, vel=v)
    print ll.classic
    return

if __name__ == "__main__":
    example()
