# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 15:37:08 2012

@author: kadu

Files with classes for reduction of spectra.
"""

import pyfits as pf
import numpy as np

class MultiSpec:
    """ This classes is intended to deal with Multi Spectra of instruments 
        with fibers, i.e., Blanco's Hydra. This class requires one image:
        containing header information about the position of the fibers and
        one spectrum containing the data.
        """
    def __init__(self, im, spec, verbose=False, instrument="hydra"):
        self.image = im
        self.id = im.replace('.fits', '')
        self.spec = spec
        self.verbose = verbose        
        self.data = pf.getdata(spec)
        self.header = pf.getheader(spec)
        self.imheader = pf.getheader(im)
        self.shape_info()
        self.check_fields()
        if instrument == "hydra": 
            self.reverse_spec()
        return

    def check_fields(self):
        """Check mandatory fields in the header, and update to default values
        in case they don't exist. '"""
        if self.refpix() == None: 
            if self.verbose: print " Warning: CRPIX1 param does not exist."
            self.header.update("CRPIX1", 1.)
        if self.refpix_value() == None: 
            if self.verbose: print " Warning: CRVAL1 param does not exist."            
            self.header.update("CRVAL1", 1.)
        if self.refpix_delta() == None: 
            if self.verbose: print " Warning: CDELT1 param does not exist."            
            self.header.update("CDELT1", 1.)
        return
        
    def center_coords(self):
        """ Returns the central coordinates of the field using 
            header information."""
        s = self.imheader.get('FLDCNTR')
        ra = self.decimal_ra(s.split()[0])
        dec = self.decimal_dec(s.split()[1])
        return [ra, dec]
        
    def decimal_ra(self, obj):
        """ Convert a string or a list of values of right ascencion 
            from HH:MM:SS to decimal coordinates in degrees. """
        if isinstance(obj, str): 
            l = obj.split(':')
            ra_hours = (float(l[0]) + float(l[1]) / 60 + float(l[2]) / 3600)
            return 15. * ra_hours
        elif isinstance(obj, list):
            data = [[float(y) for y in x.split(':')] for x in obj]
            ra = [(x[0] + x[1]/60 + x[2]/3600) for x in data]
            return 15. * np.array(ra)
    
    def decimal_dec(self, obj):
        """ Convert a string or a list of values of declination 
            from DD:MM:SS to decimal coordinates in degrees. """
        if isinstance(obj, str):
            sign = -1. if obj.startswith('-') else 1.
            l = obj.replace('-', '').split(':')
            return (float(l[0]) + float(l[1]) / 60. + float(l[2]) /3600) *sign
        elif isinstance(obj, list):
            data = [[float(y) for y in x.split(':')] for x in obj]
            dec = [(x[0] + x[1]/60 + x[2]/3600) for x in data]
            return np.array(dec)
    
    def fiball(self):
        """ This task defines the default set of fibers to be considered,
            which are all the fibers available. """
        return np.arange(self.nfibers)+1
            
    def fcoords(self, fibers=False):
        """ Give the coordinates for one fiber if fiber is a float or
            a list of coordinates if a list is given. """
        if isinstance(fibers, np.ndarray):
            ras = np.empty_like(fibers).astype(np.float64)
            decs = np.empty_like(fibers).astype(np.float64)
            for i,f in enumerate(fibers):
                ras[i], decs[i] = self.fcoords(f)
            return [ras, decs]
        if fibers == False: fibers = self.fiball()
        if isinstance(fibers, int):
            s = self.imheader.get("SLFIB%i" % fibers).split()
            if len(s) >= 4:
                ra = self.decimal_ra(s[2])
                dec = self.decimal_dec(s[3])
                return [ra, dec]
            else:
                return [0., 0.]
            
    def fradius(self, fibers=False):
        """ Returns the angular distance to the center of the field for a 
            given fiber or array of fibers. """
        ra, dec = self.fcoords(fibers) 
        rac, decc = self.center_coords()
        return np.sqrt(np.power(ra - rac, 2) + np.power(dec - decc, 2))
            
    def fspec(self, fiber):
        """ Returns the spectrum of a given fiber. """
        return self.xarray(), self.data[fiber - 1]
        
    def hydra_fibselect(self, ftype = '1'):
        """ Define the fibers that are useful using the information in the
            image header. The number trailing the fiber number gives 
            this piece of information: 0 for sky fibers, 1 for observed 
            objects, -1 for broken fibers. """
        fibers = [pf.getval(self.image, "SLFIB%i" % (i+1)) for i in
                                                       range(self.nfibers)]
        fibers = [x.split() for x in fibers]
        return np.array([int(x[0]) for x in fibers if x[1] == ftype])

    def refpix(self):
        """ Get the reference pixel for the x coordinate. """
        return self.header.get("CRPIX1")
        
    def refpix_delta(self):
        """Get the displacement between pixel for the x axis. """
        return self.header.get("CDELT1")
    
    def refpix_value(self):
        """Get the value of the reference pixel. """
        return self.header.get("CRVAL1")

    def shape_info(self):
        """ Get information about the dimensions of the spectrum. """
        self.nfibers, self.xsize = np.shape(self.data)
        if self.verbose:
            print " File with %s fibers and %s pixels long." % (self.nfibers, 
                                                                self.xsize) 
    def reverse_spec(self):
        """ For Hydra data, the number of the fiber and the aperture number 
            in the spectra are in inverse order. To deal with this, this task
            simply reverse the data order, in a way that the number of the 
            fiber (nf) is related to the array index i simply by 
            nf = i+1, as usual. """
        self.data = self.data[::-1]
        return
                                                                
    def xarray(self):
        """ Make an array for the x axis using information about the 
            reference pixel in the header. """
        xmin = self.refpix_value()
        xmax = self.refpix_value() + self.xsize * self.refpix_delta()-1
        return np.arange(xmin, xmax, self.refpix_delta())     
        