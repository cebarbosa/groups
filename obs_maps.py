# -*- coding: utf-8 -*-
"""

Created on 26/03/16

@author: Carlos Eduardo Barbosa

Make images of groups and members

"""
import os

import numpy as np
import pyfits as pf
from astropy.io import ascii
import astropy.coordinates as coord
import astropy.units as units
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

from config import *
from run_ppxf import wavelength_array

def make_table_coordinates(redo=False):
    """ Make table with coordinates of all spectra. """
    output = os.path.join(tables_dir, "spec_radec.dat")
    if not redo:
        return
    coords =ascii.read(os.path.join(tables_dir, "coords.dat"))
    ra = coord.Angle(coords["ra"], unit=units.hour)
    dec = coord.Angle(coords["dec"], unit=units.degree)
    from config import rac, decc
    rac = coord.Angle(rac, unit=units.hour)
    decc = coord.Angle(decc, unit=units.degree)
    table = []
    for c,r,d in zip(coords, ra.degree, dec.degree):
         radius = np.sqrt((r - rac.degree)**2 + (d - decc.degree)**2)
         idx = np.argmin(radius)
         group = groups[idx]
         line = "{0:30s}{1:18s}{2:15s}{3:15s}{4:15s}{5:>15.9f}{6:>15.9f}" \
                "{7:>15.9f}".format(c["spec"], c["obsrun"], c["ra"], c["dec"],
                                    group, r, d, radius[idx])
         table.append(line)
    with open(output, "w") as f:
        f.write("#Spec Obsrun RAstr DECstr Group RA DEC RADIUS(deg)\n")
        f.write("\n".join(table ))
    return

def get_observational_data():
    """ Read files with observational data. """
    obs = ascii.read(os.path.join(tables_dir, "spec_radec.dat"))
    obs = [obs["Spec"], obs["Group"], obs["RA"], obs["DEC"], obs["RADIUS(deg)"]]
    candidates = ascii.read(os.path.join(tables_dir, "candidates_radec.dat"))
    cand_groups = np.array([x.split("_")[0] for x in candidates["specs"]])
    cand_ra = coord.Angle(candidates["RA"], unit=units.hour)
    cand_ra = cand_ra.degree
    cand_dec = coord.Angle(candidates["DEC"], unit=units.degree)
    cand_dec = cand_dec.degree
    cand = [candidates["specs"], cand_groups, cand_ra, cand_dec]
    return obs, cand

if __name__ == "__main__":
    make_table_coordinates()
    obs, candidates = get_observational_data()
    for i,group in enumerate(groups):
        print rac[i], decc[i]
        continue
        cat = ascii.read(os.path.join(home, "images/{0}_2mass.tbl".format(group)))
        ra2mass = coord.Angle(cat["ra"], unit=units.degree)
        dec2mass = coord.Angle(cat["dec"], unit=units.degree)
        idx = np.where(obs[1] == group)
        idx2 = np.where(candidates[1] == group)
        imgfile = os.path.join(home, "images/{0}_dss.fits".format(group))
        im = pf.getdata(imgfile)
        ra_im = wavelength_array(imgfile, axis=1)
        dec_im = wavelength_array(imgfile, axis=2)
        fig, ax = plt.subplots()
        ax.minorticks_on()
        plt.imshow(im, origin="bottom",
                   extent=[ra_im[0], ra_im[-1], dec_im[0], dec_im[-1]])
        for r,d in zip(obs[2][idx], obs[3][idx]):
            circle = plt.Circle((r,d), 30/3600., color="r", fill=False,)
            ax.add_artist(circle)
        for r,d in zip(candidates[2][idx2], candidates[3][idx2]):
            circle = plt.Circle((r,d), 60/3600., color="y", fill=False,)
            ax.add_artist(circle)
        for r,d in zip(cat["ra"], cat["dec"]):
            circle = plt.Circle((r,d), 40/3600., color="k", fill=False)
            ax.add_artist(circle)


        # plt.plot(ra[idx], dec[idx], "o", mec="r", markerfacecolor="none",
        #          mew=0.5)
        plt.show()



