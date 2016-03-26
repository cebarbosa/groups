# -*- coding: utf-8 -*-
"""

Created on 22/03/16

@author: Carlos Eduardo Barbosa

Get coordinates for different spectra and fibers.

"""
import os

import numpy as np
import pyfits as pf

from config import *

if __name__ == "__main__":
    wdir = os.path.join(home, "data/reduced")
    candidates = np.loadtxt(os.path.join(tables_dir, "candidates.txt"),
                            dtype=str)
    specs = []
    for night in nights:
        specs.append([x for x in os.listdir(os.path.join(wdir, night)) if
                      x.endswith(".fits")])
    radecs = []
    for obj in candidates:
        matches = [] # Check if matches are unique
        coords = []
        for i, night in enumerate(nights):
            match = [s for s in specs[i] if "{0}.fits".format(obj) in s]
            for m in match:
                os.chdir(os.path.join(wdir, night))
                finfo = pf.getval(m, "FINFO")
                radec = " ".join(finfo.split()[2:4])
                coords.append(radec)
                matches.append(m)
        radecs.append("{0:15s}{1}".format(obj, list(set(coords))[0]))
    with open(os.path.join(tables_dir, "candidates_radec.dat"), "w") as f:
        f.write("\n".join(radecs))



