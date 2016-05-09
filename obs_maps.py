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
from matplotlib.collections import PatchCollection
import pywcs

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
         if str(c["mask"]).lower().startswith(("hr", "std", "hd")):
             continue
         mask =  c["mask"].lower().replace("h22", "hcg22").replace("h42",
                 "hcg42").replace("h62", "hcg62").replace("mask_1", "mask1").upper()
         line = "{0:30s}{1:18s}{2:25s}{3:15s}{4:15s}{5:15s}{6:>15.9f}{7:>15.9f}" \
                "{8:>15.9f}".format(c["spec"], c["obsrun"], mask,
                c["ra"], c["dec"], group, r, d, radius[idx])
         table.append(line)
    with open(output, "w") as f:
        f.write("#Spec Obsrun Mask RAstr DECstr Group RA DEC RADIUS(deg)\n")
        f.write("\n".join(table))
    return

def get_observational_data():
    """ Read files with observational data. """
    obs = ascii.read(os.path.join(tables_dir, "spec_radec.dat"))
    obs = [obs["Spec"], obs["Group"], obs["RA"], obs["DEC"],
           obs["RADIUS(deg)"], obs["Mask"]]
    candidates = ascii.read(os.path.join(tables_dir, "candidates_radec.dat"))
    cand_groups = np.array([x.split("_")[0] for x in candidates["specs"]])
    cand_ra = coord.Angle(candidates["RA"], unit=units.hour)
    cand_ra = cand_ra.degree
    cand_dec = coord.Angle(candidates["DEC"], unit=units.degree)
    cand_dec = cand_dec.degree
    cand = [candidates["specs"], cand_groups, cand_ra, cand_dec]
    return obs, cand

def print_dss_coords(ra1, ra2, dec1, dec2):
    print 0.5 * (ra1 + ra2), 0.5 * (dec1 + dec2), 60 * (ra2 - ra1),
    print 60 * (dec2 - dec1)
    return

def get_zm98(group):
    """ Add galaxies from the catalog of Zabludoff and Mulchaey 1998."""
    gname = group.replace("hcg", "h").upper()
    path = os.path.join(tables_dir, "zm98")
    table = ascii.read(os.path.join(path, 'table1.dat'),
                       readme=os.path.join(path, "ReadMe"))
    ra = coord.Angle(table["RAh"] + table["RAm"]/60. + table["RAs"]/3600.,
                     unit=units.hour)
    signal = np.ones_like(ra.degree)
    signal[table["DE-"]=="-"] = -1
    dec = coord.Angle(signal * (table["DEd"] + table["DEm"]/60. +
                                table["DEs"]/3600.), unit=units.degree)
    objs = np.array([x.split("_")[0] for x in table["Name"]])
    idx = np.where(objs==gname)
    return table["Name"][idx], ra[idx], dec[idx], table["HRV"][idx]

def make_map(group):
    """ Make the plot"""
    obs, candidates = get_observational_data()
    masks_labels = {"hcg22" : ["HCG 22 Mask 1"], "hcg42" : ["HCG 42 Mask 1"],
             "hcg62" : ["HCG 62 Centre", "HCG 62 Outskirt 2",
                        "HCG 62 Outskirt 1"],
             "hcg90" : ["HCG 90 Mask 1", "HCG 90 Mask 2"],
             "ngc193" : ["NGC 193 Mask 1", "NGC 193 Mask 2"],
             "ngc7619" : ["NGC 7619 Mask 1", "NGC 7619 Mask 2"]}
    #######################################################################
    # Display image
    imgfile = os.path.join(home, "images/{0}/{0}_dss.fits".format(group))
    im = pf.getdata(imgfile)
    ydim, xdim = im.shape
    wcs = pywcs.WCS(pf.getheader(imgfile))
    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.subplots_adjust(left=0.15, right=0.98, top=0.98)
    ax.minorticks_on()
    ra_im = wavelength_array(imgfile, axis=1)
    dec_im = wavelength_array(imgfile, axis=2)
    plt.imshow(im, origin="bottom", cmap="cubehelix_r", vmin=6000, vmax=8000,
               extent=[ra_im[0], ra_im[-1], dec_im[0], dec_im[-1]])
    #######################################################################
    # Show all observed fibers
    idx = np.where(obs[1] == group)
    data = list(set(zip(obs[2][idx], obs[3][idx], obs[5][idx])))
    ras, decs, masks = np.array(list(data)).T
    ras = ras.astype(float)
    decs = decs.astype(float)
    colors = ["r", "b", "y"]
    for j, mask in enumerate(set(masks)):
        idx = np.where(masks == mask)
        x, y = wcs.wcs_sky2pix(np.column_stack((ras[idx], decs[idx])), 1).T
        r = x/xdim*(ra_im[-1]-ra_im[0])+ra_im[0]
        d = y/ydim*(dec_im[-1]-dec_im[0])+dec_im[0]
        ax.plot(r, d, "o", c="none", ms=10, mec=colors[j],
                label=masks_labels[group][j])
    ax.set_xlabel("RA J2000 (deg)")
    ax.set_ylabel("DEC J2000 (deg)")
    #########################################################################
    # Indicate candidates
    idx = np.where(candidates[1] == group)
    x, y = wcs.wcs_sky2pix(np.column_stack((candidates[2][idx], \
                                            candidates[3][idx])), 1).T
    r = x/xdim*(ra_im[-1]-ra_im[0])+ra_im[0]
    d = y/ydim*(dec_im[-1]-dec_im[0])+dec_im[0]
    ax.plot(r, d, "s", c="none", ms=10, mec="k",
                label="This work members")
    ##########################################################################
    # Plot galaxies from Z&M 98
    if group in ["hcg62", "hcg42", "hcg90"]:
        objs, ra, dec, vel = get_zm98(group)
        idx = np.where((vel < v0s[group] + 2000.) & (vel > v0s[group] - 2000))
        x, y = wcs.wcs_sky2pix(np.column_stack((ra.degree[idx], \
                                                dec.degree[idx])), 1).T
        r = x/xdim*(ra_im[-1]-ra_im[0])+ra_im[0]
        d = y/ydim*(dec_im[-1]-dec_im[0])+dec_im[0]
        ax.plot(r, d, "^", c="none", ms=10, mec="g",
                label="ZM98 members")
    ##########################################################################
    outdir = os.path.join(home, "plots")
    output = os.path.join(outdir, "obs_{0}.png".format(group))
    ax.set_xlim(ra_im[0], ra_im[-1])
    ax.set_ylim(dec_im[0], dec_im[-1])
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    plt.legend(loc=2,prop={'size':10}, frameon=True)
    plt.savefig(output)
    # plt.show()
    # plt.close()
    return

if __name__ == "__main__":
    make_table_coordinates(redo=False)
    # make_map("hcg22")
    # make_map("hcg42")
    # make_map("hcg62")
    # make_map("hcg90")
    # make_map("ngc7619")
    make_map("ngc193")

