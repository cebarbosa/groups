# -*- coding: utf-8 -*-
"""

Created on 02/03/16

@author: Carlos Eduardo Barbosa

Check the candidates for membership of groups of galaxies using results from
FXCOR

"""

import os

import numpy as np

from config import *

if __name__ == "__main__":
    wdir = os.path.join(home, "Blanco/members")
    candidates = np.loadtxt(os.path.join(home, "tables/candidates.txt"),
                            dtype=str)
    for group in groups:
        print group.upper(), v0s[group]
        table = os.path.join(wdir, group, "recession_velocities.txt")
        specs = np.loadtxt(table, usecols=(0,), dtype=str)
        specs = [x.replace("-", "_") for x in specs]
        specs = np.array(["_".join([group] + x.split("_")[1:]) for x in specs])
        vs = np.loadtxt(table, usecols=(1,))
        idx = np.abs(vs - v0s[group]) < 2000.
        for spec, v in zip(specs[idx], vs[idx]):
            if spec.replace(".fits", "") not in candidates:
                print spec, v
        print "\n"