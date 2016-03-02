# -*- coding: utf-8 -*-
"""

Created on 02/03/16

@author: Carlos Eduardo Barbosa

Make a list containing the directory of the candidate spectra

"""
import os

import numpy as np

from config import *

if __name__ == "__main__":
    candidates = np.loadtxt(os.path.join(tables_dir, "candidates.txt"),
                            dtype=str)
    table = []
    for night in nights:
        os.chdir(os.path.join(data_dir, night))
        fits = sorted([x for x in os.listdir(".") if x.endswith(".fits")])
        good = [x for x in fits if x.replace(".fits", "") in candidates]
        n = [night] * len(good)
        table += zip(good, n)
    table = ["{0:20s}{1:20s}".format(x,y) for (x,y) in table]
    with open(os.path.join(tables_dir, "spec_location.txt"), "w") as f:
        f.write("{0:20s}{1:20s}\n".format("#Spectrum", "Directory"))
        f.write("\n".join(table))