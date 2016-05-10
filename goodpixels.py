# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 16:48:17 2013

@author: cbarbosa

Program so define regions to be used by PPXF 
"""

import os

import numpy as np
import pyfits as pf
import matplotlib.pyplot as plt

from config import *
from run_ppxf import wavelength_array

def set_coords(spec):
    fig = plt.figure(1)
    ax = plt.subplot(111)
    data = pf.getdata(spec)
    wave = wavelength_array(spec)
    coords = [wave[0], wave[-1]]
    plt.cla()
    ax.plot(wave, data, "-k")
    for x in coords:
        ax.axvline(x, ls="--", c="r")
    plt.pause(0.001)
    plt.show(   )
    accept = raw_input("Accept x0={0[0]}, x1={0[1]}? (Y/n)".format(coords))
    if accept.lower().strip() in ["", "y", "ye", "yes"]:
        return coords
    x0 = raw_input("New x0: ")
    x1 = raw_input("New x1: ")
    coords = [float(x0), float(x1)]
    return coords

if __name__ == "__main__":
    os.chdir(data_dir)
    specs = sorted([x for x in os.listdir(".") if x.endswith(".fits")])
    filename = "xlims.txt"
    specs_done = []
    if not os.path.exists(filename):
        with open("xlims.txt", "w") as f:
            f.write("# Spec x0 x1 \n")
    else:
        specs_done = np.loadtxt(filename, dtype=str, usecols=(0,))
    specs = [x for x in specs if x not in specs_done]
    for i, spec in enumerate(specs):
        print "Working with spectrum {0} ({1} of {2})".format(spec, i+1,
                                                              len(specs))
        coords = set_coords(spec)
        with open(filename, "a") as f:
            f.write("{0} {1[0]:.1f} {1[1]:.1f}\n".format(spec, coords))