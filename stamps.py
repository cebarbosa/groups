# -*- coding: utf-8 -*-
"""

Created on 13/05/16

@author: Carlos Eduardo Barbosa

Produce stamps for candidates.

"""
import os
from subprocess import call

from astropy.io import ascii

from config import *

def download_images():
    """ Use dss server to retrieve images for candidates. """
    os.chdir(wdir)
    data = ascii.read(os.path.join(tables_dir, "candidates_radec.dat"))
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

if __name__ == "__main__":
    wdir = os.path.join(home, "images/dss")
    download_images()


