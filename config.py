# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 17:24:25 2014

@author: kadu

Configuration file for the work in the groups of galaxies for my PhD work.
"""

import os 

home = "/data/kadu/Dropbox/groups"
templates_dir = os.path.join(home, "miles_models")


# Constants
c = 299792.458 # Speed of light
FWHM_tem = 2.54 # MILES library spectra have a resolution FWHM of 2.54A.
FWHM_spec = 3.6 # FWHM of data

# Resolution for binning with pPXF
velscale = 30. # km/s
