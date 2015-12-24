# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 12:56:26 2013

@author: kadu

This program splits spectra from multifiber files into individual spectra. 

It is possible to set the spectra type using the spectype parameter. In case 
of 0, sky spectra are extracted, in case of 1 object spectra is extracted. 

Notice that, in case of object spectra, it is necessary to use the routine 
renaming.py to stardartize the name of the files before stacking. 

"""
import os

import numpy as np
import pyfits as pf
from pyraf import iraf, iraffunctions

from multispec import MultiSpec

iraf.onedspec(_doprint=0)

def chdir(path):
    """General chdir updating also pyraf's path """
    os.chdir(path)
    iraffunctions.chdir(path)
    return

if __name__ == "__main__":
    home = '/home/kadu/Dropbox/groups/Blanco/dohydra/'
    outhome = '/home/kadu/Dropbox/groups/data/reduced'
    dirs = [x for x in os.listdir(home) if
            os.path.isdir(os.path.join(home, x))]
    spectype = '1'
    for d in dirs:
        chdir(os.path.join(home, d))
        outdir = os.path.join(outhome, d)
        print outdir
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        sample = np.loadtxt("objs.txt", dtype=str).tolist()
        for im in sample:
            spectrum = im.replace('.fits', '.ms.fits')
            if not os.path.exists(spectrum):
                continue
            print "working spectrum: ", spectrum
            mspec = MultiSpec(im, spectrum)
            objfibers = mspec.hydra_fibselect(ftype= spectype)
            print objfibers
            for i,f in enumerate(objfibers):
                print f,
                finfo = pf.getval(mspec.image, 'SLFIB%s' % f)
                name = finfo.split()[4].lower()
                outspec_name = '%s_%s' % (mspec.image.replace('.fits', ''), name)
                print outspec_name
                outspec = os.path.join(outdir, outspec_name)
                iraf.scopy(input=spectrum, output=outspec, w1="INDEF", w2="INDEF",
                           apertures=objfibers[i], clobber = 'yes')
        
    
    