# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 14:41:39 2013

@author: cbarbosa

Lick indices measurement program. 

Last update: August 6, 2013
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from scipy.integrate import romb
from scipy.ndimage.filters import gaussian_filter1d
from scipy.constants import c

c /= 1000. # Convert to km / s

def lector(wl, intens, noise, infile, vel=0, cols=(0,8,2,3,4,5,6,7)):
    """ Make the measurement of Lick indices from file infile in spectrum. 
    
    This function make a direct measurement of the Lick (or any other 
    absorption line feature defines as the Lick indices). The 
    input spectra has to be already broadened to the correct values for 
    the measurement previous to this task. The errors are calculated using 
    the analytic expressions of Cardiel et al. 1998. 
    
    ================
    Input parameters 
    ================    
    wl : array
        Wavelenght numpy array. Units have to be consistent with all other 
        wavelenght parameters.    
              
    intens : array
        Intensity numpy array. Can be used in arbitrary units (au).
        
    noise : array
        Noise vector in the same units as the intensity vector.  
    
    infile : string
        Filename with definition of indices to be used. See cols.
    
    vel : float
        Relative velocity of the object, in km/s. 
    
    cols : array
        Order of columns in the infile. The required fields are:
        0. Index name
        1. Index type number: 0 and 2 for EW, 1 for magnitude 
        2. 3. Blue continuum wavelengths (blue and red)
        4. 5. Indice wavelengths (blue and red)
        6. 7. Red continuum (blue and red)
    
    =================
    Output parameters
    =================
    results : array
        Numpy array with the measured indices. 
    
    errors : array
        Numpy array with indice errors according to analytical expressions of 
        Cardiel et al. 1998.
        
    """
    k_order=10.
    # Define number of points for integration
    npoints = 2**k_order + 1
    # Calculate dispersion
    disp = wl[1] - wl[0]
    # Read file for indices definitions
    indnames = np.loadtxt(infile, usecols = (cols[0],), dtype='|S16')
    indtype = np.loadtxt(infile, usecols = (cols[1],))
    indices = np.loadtxt(infile, usecols = cols[2:])
    # Correct for velocity recession
    indices *= np.sqrt((1 + vel/c)/(1 - vel/c))
    # Initiate array of index
    results = np.empty((len(indnames)))
    errors = np.empty((len(indnames)))
    # Measure the indices
    for i, w in enumerate(indices): 
        # Test whether indice in defined in the wavelenght domain
        if wl[0] > w[0] or wl[-1] < w[5]:
            results[i] = np.nan
            errors[i] = np.nan
            continue
        itype = indtype[i]
        # Defining indices for each section
        sec_blue = np.where(((wl > w[0] - 2 * disp) & (wl < w[1] + 2 * disp)))
        sec_red = np.where(((wl > w[4] - 2 * disp) & (wl < w[5] + 2 * disp)))
        sec_ind = np.where(((wl > w[2] - 2 * disp) & (wl < w[3] + 2 * disp)))
        # Defining wavelenght samples
        w_blue = wl[sec_blue]
        w_red = wl[sec_red]
        w_ind = wl[sec_ind]
        # Defining intensity samples
        int_blue = intens[sec_blue]
        int_red = intens[sec_red]
        int_ind = intens[sec_ind]
        # Making interpolation functions
        s_blue = interp1d(w_blue, int_blue, kind="linear")
        s_red = interp1d(w_red, int_red, kind="linear")
        # Make oversampled arrays of wavelenghts
        xblue = np.linspace(w[0], w[1], npoints)
        xred = np.linspace(w[4], w[5], npoints)
        xind = np.linspace(w[2], w[3], npoints)
        # Calculating the mean fluxes for the pseudocontinuum
        fp1 = romb(s_blue(xblue), dx=xblue[1] - xblue[0])/ (w[1] - w[0])
        fp2 = romb(s_red(xred), dx=xred[1] - xred[0]) / (w[5] - w[4])
        # Making pseudocontinuum vector
        x0 = (w[2] + w[3])/2.
        x1 = (w[0] + w[1])/2.
        x2 = (w[4] + w[5])/2.
        fc = fp1 + (fp2 - fp1)/ (x2 - x1) * (w_ind - x1)
        # Calculating term C2 of Cardiel et al. 1998 of formula 44 for errors
        c2 = np.sqrt( 1 / (w[3]- w[2]) +
                     np.power((x1 - x0) / (x1 - x2), 2.) / (w[5] - w[4]) + 
                     np.power((x0 - x2) / (x1 - x2), 2.) / (w[1] - w[0]))
        # Calculating S/N using Cardiel et al. 1998 formula. 
        sec = np.concatenate((sec_blue[0], sec_ind[0], sec_red[0]))
        sec = np.unique(sec)
        SN = np.sum(intens[sec]) / np.std(noise[sec]) / (len(sec) * disp)
        # Calculating index according to type: 0 and 2 in angstroms and 1 in 
        # mags, and respective errors
        if itype in [0,2]:
            si = interp1d(w_ind, 1 - int_ind / fc, kind="linear")
            c1 = (w[3] - w[2]) * c2
            results[i] = romb(si(xind), dx = xind[1] - xind[0])
            errors[i] = (c1 - c2 * results[i]) / SN
        elif itype == 1:
            si = interp1d(w_ind, int_ind / fc, kind="linear")

            results[i] = -2.5 * np.log10(romb(si(xind), 
                        dx = xind[1] - xind[0]) / (w[3] - w[2]))
            errors[i] = 2.5 * c2 * np.log10(np.e) / SN
    return results, errors
    
def broad2lick(wl, intens, obsres, vel=0):
    """ Convolve spectra to Lick resolution. 
        
    Broad a given spectra to the Lick system resolution. As the resolution 
    in the Lick system varies as function of the wavelenght, we use the 
    interpolated values from Worthey and Ottaviani 1997.
    
    ================
    Input parameters
    ================
    wl: array_like
        Wavelenght 1-D array in Angstroms. 
    
    intens: array_like
        Intensity 1-D array of Intensity, in arbitrary units. The lenght has 
        to be the same as wl. 
        
    obsres: float
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
    wlick = np.array([2000., 4000., 4400., 4900., 5400., 6000., 8000.]) * \
            np.sqrt((1 + vel/c)/(1 - vel/c))
    lickres = np.array([11.5, 11.5, 9.2, 8.4, 8.4, 9.8, 9.8])
    sigma_b = np.sqrt(lickres * lickres - obsres * obsres) / 2.3548
    sigmas = interp1d(wlick, sigma_b, kind="linear")(wl)
    intens2D = np.diag(intens)
    for i in range(len(sigmas)):
        intens2D[i] = gaussian_filter1d(intens2D[i], sigmas[i], 
                      mode="constant", cval=0.0)
    return intens2D.sum(axis=0)
    
if __name__ == "__main__":
    pass
    








