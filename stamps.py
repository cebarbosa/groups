# -*- coding: utf-8 -*-
"""

Created on 13/05/16

@author: Carlos Eduardo Barbosa

Produce stamps for candidates.

"""
import os
from subprocess import call

import pywcs
import numpy as np
import pyfits as pf
from astropy import units
import matplotlib as mpl
import astropy.coordinates as coord
import matplotlib.pyplot as plt

from config import *
from run_ppxf import wavelength_array

def download_images(data):
    """ Use dss server to retrieve images for candidates. """
    lines = []
    for d in data:
        args = [d["specs"]]
        args += d["RA"].split(":")
        args += d["DEC"].split(":")
        args += ["5", "5"]
        lines.append(" ".join(args))
    with open("input.lis", "w") as f:
        f.write("\n".join(lines))
    for band in ["red", "blue", "IR"]:
        os.chdir(os.path.join(wdir, band))
        call(["dss2", band, "-i", "../input.lis"])
    return

def plot():
    table = os.path.join(data_dir, "results.tab")
    data = np.loadtxt(table, usecols=(0,72,73,1,2,3,4,63,64,65,66,67,68,69,70,
                                      71,72), dtype=str)
    specs = data[:,0]
    ras = coord.Angle(data[:,1], unit=units.hour)
    decs = coord.Angle(data[:,2], unit=units.degree)
    # download_images()
    reddir = os.path.join(wdir, "red")
    os.chdir(reddir)
    filenames = os.listdir(".")
    for d, ra, dec in zip(data, ras, decs):
        spec = d[0]
        print spec
        ######################################################################
        galaxy = "_".join(spec.split("_")[:2])
        run  = obsrun[spec.replace(".fits", "").split("_")[2]]
        imgfile = [x for x in filenames if x.startswith(galaxy)][0]
        img = pf.getdata(imgfile)
        raim = wavelength_array(imgfile, axis=1)
        decim = wavelength_array(imgfile, axis=2)
        wcs = pywcs.WCS(pf.getheader(imgfile))
        xsize = np.abs(raim[-1] - raim[0]) * 60
        fig = plt.figure(1, figsize=(5,4.6))
        ax = plt.subplot(111)
        ax.minorticks_on()
        plt.imshow(img, origin="bottom", cmap="hot",
                   extent=[raim[0], raim[-1], decim[0], decim[-1]],
                   aspect="equal")
        # labels = [str(float(x)/15) for x in ax.get_xticks().tolist()]
        # xlabels = coord.Angle(labels, unit=units.hour)
        # ax.set_xticklabels(xlabels)
        # plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
        ax.tick_params(axis='both', which='major', labelsize=9)
        ax.set_xlabel("RA J2000 (degree)")
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        plt.ylabel("DEC J2000 (degree)")
        xy = lambda r, phi: (r * np.cos(phi), r * np.sin(phi))
        xcirc, ycirc = xy(2./60 / 2, np.arange(0,6.28,0.1))
        plt.plot(ra.degree + xcirc, dec.degree + ycirc, "--y")
        plt.xlim(raim[0], raim[-1])
        plt.ylim(decim[0], decim[-1])
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        #######################################################################
        # Include text
        v = round(float(d[3]))
        dv = round(float(d[4]))
        vstr = "$V={0:.0f}\pm{1:.0f}$ km/s".format(v,dv)
        sig = round(float(d[5]))
        dsig = round(float(d[6]))
        sigstr = "$\sigma={0:.0f}\pm{1:.0f}$ km/s".format(sig,dsig)
        prop = ["Age", "[Z/H]", r"[$\alpha$/Fe]"]
        stpop = []
        for i in range(3):
            val = float(d[7+3*i])
            lerr =  np.abs(float(d[8+3*i]))
            uerr = np.abs(float(d[9+3*i]))
            ndig = np.minimum(np.log10(lerr), np.log10(uerr))
            ndig = int(np.ceil(np.abs(ndig)))
            ndig = ".{0}".format(ndig)
            valstr = "{4}=${0:{3}f}_{{-{1:{3}f}}}^{{+{2:{3}f}}}$".format(val, lerr,
                                                        uerr, ndig, prop[i])
            stpop.append(valstr)
        line1 = "{0} ({1})".format(galaxy.replace("_", "\_").upper(), run)
        line2 = "{0:18s}, {1}".format(vstr, sigstr)
        plt.figtext(0.20, 0.9, line1, color="w")
        plt.figtext(0.20, 0.84, line2, color="w", fontsize=14)
        plt.figtext(0.20, 0.25, stpop[0]+"Gyr", color="w", fontsize=14)
        plt.figtext(0.60, 0.18, stpop[1], color="w", fontsize=14)
        plt.figtext(0.20, 0.18, stpop[2], color="w", fontsize=14)
        #######################################################################
        output = os.path.join(plots_dir, "stamps/{0}.png".format(
                 spec.replace(".fits", "")))
        plt.subplots_adjust(bottom=0.09, right=0.96, top=0.99, left=0.19)
        plt.savefig(output)
        # plt.show()
        plt.close(1)

def prepare_figures():
    os.chdir(os.path.join(plots_dir, "stamps"))
    figures = sorted(os.listdir("."))
    objs = sorted(list(set(["_".join(x.split("_")[:2]) for x in figures])))
    for obj in objs:
        filenames = [x for x in figures if x.startswith(obj)]
        print "\includegraphics[width=0.32\linewidth]{{figs/stamps/{0}}}".format(filenames[0])


if __name__ == "__main__":
    mpl.rc('text', usetex=True)
    mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
    wdir = os.path.join(home, "images/dss")
    os.chdir(wdir)
    # plot()
    prepare_figures()



