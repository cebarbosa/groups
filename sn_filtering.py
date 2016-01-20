# -*- coding: utf-8 -*-
"""

Created on 28/12/15

@author: Carlos Eduardo Barbosa

"""
import os

import numpy as np
import pyfits as pf
import matplotlib.pyplot as plt

from config import *

def wavelength_array(spec):
    """ Produces array for wavelenght of a given array. """
    w0 = pf.getval(spec, "CRVAL1")
    deltaw = pf.getval(spec, "CD1_1")
    pix0 = pf.getval(spec, "CRPIX1")
    npix = pf.getval(spec, "NAXIS1")
    return w0 + deltaw * (np.arange(npix) + 1 - pix0)

def signal_analysis(filenames, w1=4800, w2=5400):
    signal = np.zeros(len(filenames))
    for i, f in enumerate(filenames):
        print i, f
        data = pf.getdata(f)
        w = wavelength_array(f)
        goodpixels = np.logical_and(w > w1, w < w2)
        signal[i] =  np.median(data[goodpixels])
    return signal
        # ax = plt.subplot(111)
        # ax.minorticks_on()
        # ax.plot(w[goodpixels], data[goodpixels], "-k")
        # ax.axhline(s, ls="-", c="r")
        # ax.set_xlim(w1, w2)
        # ax.set_ylim(0, data[goodpixels].max())
        # ax.set_xlabel("Wavelength (Angstrom)")
        # ax.set_ylabel("Counts")
        # plt.title("{1}\tS = {0:.1f}".format(s, f.replace("_", "")))

if __name__ == "__main__":
    reduced_path = os.path.join(home, "data/reduced")
    nights = os.listdir(reduced_path)
    fig = plt.figure(1)
    for night in nights:
        os.chdir(os.path.join(reduced_path, night))
        objs = [x for x in os.listdir(".") if x.startswith("crobj") and
                x.endswith(".fits")]
        skies =  [x for x in os.listdir(".") if x.startswith("crobj") and
                 "sky" in x]
        objs = [x for x in objs if x not in skies]
        objs.sort()
        skies.sort()
        multispecs = list(set([x.split("_")[0] for x in objs]))
        multispecs.sort()
        for ms in multispecs:
            specs = [x for x in objs if x.startswith(ms)]
            specs.sort()
            sky = [x for x in skies if x.startswith(ms)]
            sky.sort()
            print specs
            print sky
            raw_input()
            signal = signal_analysis(specs)
            signal_sky = signal_analysis(sky)
            median = np.median(signal)
            mad = np.median(np.abs(signal - median))
            plt.clf()
            ax = plt.subplot(111)
            ax.hist(signal, range=(0, median + 3 * mad), bins=30)
            ax.hist(signal_sky, range=(0, median + 3 * mad), bins=30)
            plt.pause(2)
