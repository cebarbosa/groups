# -*- coding: utf-8 -*-
"""

Created on 25/03/16

@author: Carlos Eduardo Barbosa

Using pyraf to calculate VHELIO correction

"""
import os

import pyfits as pf
from pyraf import iraf

from config import *

iraf.rv(_doprint=0)

if __name__ == "__main__":
    os.chdir(data_dir)
    fits = sorted([x for x in os.listdir(".") if x.endswith(".fits")])
    table = "VHELIO"
    for img in fits:
        h = pf.getheader(img)
        print img,
        if "VHELIO" not in h.keys():
            print
            continue
            # try:
            #     iraf.hedit (images=img, fields="EPOCH", value=2000., add=True,
            #                 delete=False, verify=False, show=False,
            #                 update=True)
            #     iraf.rvcorrect(images=img, imupdat=True, epoch=2000.,
            #                    ut=h["UT"], ra=h["RA"], dec=h["DEC"])
            # except:
            #     continue
        print h["VHELIO"]


