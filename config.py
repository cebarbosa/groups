# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 17:24:25 2014

@author: kadu

Configuration file for the work in the groups of galaxies for my PhD work.
"""

home = "/home/kadu/Dropbox/groups"
templates_dir = home + "/miles_models"
data_dir = home + "/data/combined"
tables_dir = home + "/tables"


# Constants
c = 299792.458 # Speed of light
FWHM_tem = 3.6 # MILES library spectra have a resolution FWHM of 2.54A.
FWHM_spec = 3.6 # FWHM of data

# Resolution for binning with pPXF
velscale = 30. # km/s

# Velocities of the groups according to NED
v0s = dict([("hcg22", 2698.), ("hcg62", 4413.), ("hcg90", 2638.),
            ("hcg42", 3987.), ("ngc193", 4400.), ("ngc7619", 3439.)])

# Defining systems
groups = ["hcg22", "hcg42", "hcg62", "hcg90", "ngc193", "ngc7619"]

# Defining observing runs
nights = ["blanco10n1", "blanco10n2", "blanco11an1", "blanco11an2",
          "blanco11bn1", "blanco11bn2", "blanco11bn3"]
