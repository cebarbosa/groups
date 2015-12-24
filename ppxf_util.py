#######################################################################
#
# Copyright (C) 2001-2014, Michele Cappellari
# E-mail: cappellari_at_astro.ox.ac.uk
#
# This software is provided as is without any warranty whatsoever.
# Permission to use, for non-commercial purposes is granted.
# Permission to modify for personal or internal use is granted,
# provided this copyright and disclaimer are included unchanged
# at the beginning of the file. All other rights are reserved.
#
#######################################################################
#
# NAME:
#   LOG_REBIN
#
# PURPOSE:
#   Logarithmically rebin a spectrum, while rigorously conserving the flux. 
#   Basically the photons in the spectrum are simply ridistributed according 
#   to a new grid of pixels, with non-uniform size in the spectral direction.
#
#   This routine makes the `standard' zero-order assumption that the spectrum
#   is *constant* within each pixels. It is possible to perform log-rebinning
#   by assuming the spectrum is represented by a piece-wise polynomial of
#   higer degree, while still obtaining a uniquely defined linear problem,
#   but this reduces to a deconvolution and amplifies noise.
#
#   This same routine can be used to compute approximate errors
#   of the log-rebinned spectrum. To do this type the command
#
#       LOG_REBIN, lamRange, err^2, err2New
#
#   and the desired errors will be given by SQRT(err2New).
#   NB: This rebinning of the error-spectrum is very *approximate* as 
#   it does not consider the correlation introduced by the rebinning!
#
# CALLING SEQUENCE:
#   LOG_REBIN, lamRange, spec, specNew, logLam, $
#       OVERSAMPLE=oversample, VELSCALE=velScale, /FLUX
#
# INPUTS:
#   LAMRANGE: two elements vector containing the central wavelength
#       of the first and last pixels in the spectrum, which is assumed
#       to have constant wavelength scale! E.g. from the values in the
#       standard FITS keywords: LAMRANGE = CRVAL1 + [0,CDELT1*(NAXIS1-1)].
#       It must be LAMRANGE[0] < LAMRANGE[1].
#   SPEC: input spectrum.
#
# OUTPUTS:
#   SPECNEW: logarithmically rebinned spectrum.
#   LOGLAM: log(lambda) (*natural* logarithm: ALOG) of the central 
#       wavelength of each pixel. This is the log of the geometric 
#       mean of the borders of each pixel.
#
# KEYWORDS:
#   FLUX: Set this keyword to preserve total flux. In this case the 
#       log rebinning changes the pixels flux in proportion to their 
#       dLam so the following command will show large differences 
#       beween the spectral shape before and after LOG_REBIN:
#       
#           plot, exp(logLam), specNew  # Plot log-rebinned spectrum
#           oplot, range(lamRange[0],lamRange[1],n_elements(spec)), spec
#           
#       By defaul, when this keyword is *not* set, the above two lines 
#       produce two spectra that almost perfectly overlap each other.
#   OVERSAMPLE: Oversampling can be done, not to loose spectral resolution, 
#       especally for extended wavelength ranges and to avoid aliasing.
#       Default: OVERSAMPLE=1 ==> Same number of output pixels as input.
#   VELSCALE: velocity scale in km/s per pixels. If this variable is
#       not defined, then it will contain in output the velocity scale.
#       If this variable is defined by the user it will be used
#       to set the output number of pixels and wavelength scale.
#
# MODIFICATION HISTORY:
#   V1.0: Using interpolation. Michele Cappellari, Leiden, 22 October 2001
#   V2.0: Analytic flux conservation. MC, Potsdam, 15 June 2003
#   V2.1: Allow a velocity scale to be specified by the user.
#       MC, Leiden, 2 August 2003
#   V2.2: Output the optional logarithmically spaced wavelength at the
#       geometric mean of the wavelength at the border of each pixel.
#       Thanks to Jesus Falcon-Barroso. MC, Leiden, 5 November 2003
#   V2.21: Verify that lamRange[0] < lamRange[1]. 
#       MC, Vicenza, 29 December 2004
#   V2.22: Modified the documentation after feedback from James Price.
#       MC, Oxford, 21 October 2010
#   V2.3: By default now preserve the shape of the spectrum, not the 
#       total flux. This seems what most users expect from the procedure.
#       Set the keyword /FLUX to preserve flux like in previous version.
#       MC, Oxford, 30 November 2011
#   V2.4: Translated from IDL into Python. MC, Santiago, 23 November 2013
#   V3.0: Fully vectorized. Typical speed up by two orders of magnitude.
#       MC, Oxford, 4 March 2014
#
#----------------------------------------------------------------------

import numpy as np    

def log_rebin(lamRange, spec, oversample=False, velscale=None, flux=False):
    """
    Logarithmically rebin a spectrum, while rigorously conserving the flux. 
    Basically the photons in the spectrum are simply ridistributed according 
    to a new grid of pixels, with non-uniform size in the spectral direction.
    
    """
    lamRange = np.asarray(lamRange)
    if len(lamRange) != 2:
        raise ValueError('lamRange must contain two elements')
    if lamRange[0] >= lamRange[1]:
        raise ValueError('It must be lamRange[0] < lamRange[1]')
    s = spec.shape
    if len(s) != 1:
        raise ValueError('input spectrum must be a vector')
    n = s[0]
    if oversample:
        m = int(n*oversample)
    else:
        m = int(n)
    
    dLam = np.diff(lamRange)/(n - 1.)        # Assume constant dLam
    lim = lamRange/dLam + [-0.5, 0.5]        # All in units of dLam
    borders = np.linspace(*lim, num=n+1)     # Linearly
    logLim = np.log(lim)
    
    c = 299792.458                           # Speed of light in km/s
    if velscale is None:                     # Velocity scale is set by user
        velscale = np.diff(logLim)/m*c       # Only for output
    else:
        logScale = velscale/c
        m = int(np.diff(logLim)/logScale)    # Number of output pixels
        logLim[1] = logLim[0] + m*logScale
    
    newBorders = np.exp(np.linspace(*logLim, num=m+1)) # Logarithmically
    k = (newBorders - lim[0]).clip(0, n-1).astype(np.int)
     
    specNew = np.add.reduceat(spec, k)[:-1]  # Do analytic integral
    specNew *= np.diff(k) > 0    # fix for design flaw of reduceat()
    specNew += np.diff((newBorders - borders[k])*spec[k])

    if not flux:
        specNew /= np.diff(newBorders)

    # Output log(wavelength): log of geometric mean
    logLam = np.log(np.sqrt(newBorders[1:]*newBorders[:-1])*dLam)

    return specNew, logLam, velscale

#----------------------------------------------------------------------
#
# PPXF_DETERMINE_GOODPIXELS: Example routine to generate the vector of goodPixels 
#     to be used as input keyword for the routine PPXF. This is useful to mask 
#     gas emission lines or atmospheric absorptions. 
#     It can be trivially adapted to mask different lines.
# 
# INPUT PARAMETERS:
# - LOGLAM: Natural logarithm ALOG(wave) of the wavelength in Angstrom 
#     of each pixel of the log rebinned *galaxy* spectrum.
# - LAMRANGETEMP: Two elements vectors [lamMin2,lamMax2] with the minimum and
#     maximum wavelength in Angstrom in the stellar *template* used in PPXF.
# - VEL: Estimate of the galaxy velocity in km/s.
# 
# V1.0: Michele Cappellari, Leiden, 9 September 2005
# V1.01: Made a separate routine and included additional common emission lines. 
#   MC, Oxford 12 January 2012
# V2.0: Translated from IDL into Python. MC, Oxford, 10 December 2013
# V2.01: Updated line list. MC, Oxford, 8 January 2014

def determine_goodpixels(logLam, lamRangeTemp, vel):
    """
    Generates a list of goodpixels to mask a given set of gas emission
    lines. This is meant to be used as input for PPXF.
    
    """
#                     -----[OII]-----    Hdelta   Hgamma   Hbeta   -----[OIII]-----   [OI]    -----[NII]-----   Halpha   -----[SII]-----  
    lines = np.array([3726.03, 3728.82, 4101.76, 4340.47, 4861.33, 4958.92, 5006.84, 6300.30, 6548.03, 6583.41, 6562.80, 6716.47, 6730.85])
    dv = lines*0 + 800 # width/2 of masked gas emission region in km/s
    c = 299792.458 # speed of light in km/s
    
    flag = logLam < 0 # empy mask    
    for j in range(lines.size):
        flag |= (logLam > np.log(lines[j]) + (vel - dv[j])/c) \
              & (logLam < np.log(lines[j]) + (vel + dv[j])/c)
    
    flag |= logLam < np.log(lamRangeTemp[0]) + (vel + 900.)/c # Mask edges of
    flag |= logLam > np.log(lamRangeTemp[1]) + (vel - 900.)/c # stellar library
    
    return np.where(flag == 0)[0]

#------------------------------------------------------------------------------
# V1.0: Michele Cappellari, Oxford, 7 January 2014

def emission_lines(logLam_temp, FWHM_gal):
    """
    Generates an array of Gaussian emission lines to be used as templates in PPXF.
    logLam is the natural log of the wavelength of the templates in Angstrom.
    logLam should be the same as that of the stellar templates.
    FWHM_gal is the FWHM of the galaxy spectrum under study in Angstrom. 
    
    """
#   In this routine all lines are free to have independent intensities
#   One can fix the intensity ratio of different lines (e.g. the [OIII] doublet) 
#   by placing them in the same emission template
    
#                     -----[OII]-----    Hdelta   Hgamma   Hbeta   -----[OIII]-----   [OI]    -----[NII]-----   Halpha   -----[SII]-----  
    lines = np.array([3726.03, 3728.82, 4101.76, 4340.47, 4861.33, 4958.92, 5006.84, 6300.30, 6548.03, 6583.41, 6562.80, 6716.47, 6730.85])
    lam = np.exp(logLam_temp)
    lines = lines[(lines > np.min(lam)) & (lines < np.max(lam))]
    sigma = FWHM_gal/2.355 # Assumes instrumental sigma is constant in Angstrom
    emission_lines = np.exp(-0.5*((lam[:,np.newaxis] - lines)/sigma)**2)
            
    return emission_lines
    
#------------------------------------------------------------------------------
