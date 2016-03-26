# -*- coding: utf-8 -*-
"""

Created on 04/01/16

@author: Carlos Eduardo Barbosa

Combine spectra night by night

"""

import os
import shutil

import numpy as np
import pyfits as pf
from pyraf import iraf, iraffunctions
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from config import *
from run_ppxf import wavelength_array

iraf.onedspec(_doprint=0, Stdout=1)

def select_specs(specs):
    specs = np.array(specs, dtype=str)
    exptimes = np.zeros(len(specs))
    for i, spec in enumerate(specs):
        exptimes[i] = pf.getval(spec, "exptime")
    index = exptimes > 100.
    return specs[index].tolist()

if __name__ == "__main__":
    wdir = os.path.join(home, "data/reduced")
    outroot = wdir.replace("reduced", "combined")
    outroot2 = wdir.replace("reduced", "single")
    if not os.path.exists(outroot):
        os.mkdir(outroot)
    os.chdir(wdir)
    for night in nights:
        print "Working in night ", night
        outdir = os.path.join(outroot, night)
        outdir2 = os.path.join(outroot2, night)
        os.chdir(os.path.join(wdir, night))
        iraffunctions.chdir(os.path.join(wdir, night))
        outfile = PdfPages("log_scombine.pdf")
        fig = plt.figure(1)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if not os.path.exists(outdir2):
            os.mkdir(outdir2)
        fits = [x for x in os.listdir(".") if x.endswith("fits")]
        objs = set(["_".join(x.split("_")[1:]) for x in fits])
        for i, obj in enumerate(objs):
            print "File: {0} ({1}/{2})".format(obj, i+1, len(objs))
            lis = [x for x in fits if x.endswith(obj)]
            goodlis = select_specs(lis)
            filenames = ", ".join(goodlis)
            output = os.path.join(outdir, obj)
            if  len(goodlis) == 1:
                shutil.copy(goodlis[0], os.path.join(outdir2,
                    goodlis[0].replace("_hcg_", "_").replace("_h62_", "_")))


            continue
            if os.path.exists(output):
                continue
            iraf.scombine(input = filenames, output = output, group = 'all',
                          combine = 'sum', reject="none",
                          weight="none", w1 = 4000., w2=6500.,dw=1.25)
            ax = plt.subplot(111)
            for l in lis:
                w = wavelength_array(l)
                intens = pf.getdata(l)
                c = "k" if l in goodlis else "0.5"
                ax.plot(w, intens, "-", color=c, lw=0.4)
            w = wavelength_array(output)
            intens = pf.getdata(output)
            ax.plot(w, intens, "-r", lw=2)
            ax.set_xlabel("Wavelength (Angstrom)")
            ax.set_ylabel("Counts")
            plt.title(obj.replace("_", "-"))
            plt.pause(0.1)
            plt.show(block=0)
            outfile.savefig()
            plt.clf()
        outfile.close()



