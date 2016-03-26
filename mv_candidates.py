# -*- coding: utf-8 -*-
"""

Created on 23/03/16

@author: Carlos Eduardo Barbosa

Move candidates to a new working directory

"""
import os
import shutil

import numpy as np

from config import *

if __name__ == "__main__":
    wdir = os.path.join(home, "data/combined")
    outdir = os.path.join(home, "data/candidates")
    if not os.path.exists(outdir):
        os.mkdir(outdir)
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
                shutil.move(m, os.path.join(outdir,
                                m.replace(".fits", "_{0}.fits".format(night))))