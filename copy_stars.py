# -*- coding: utf-8 -*-
"""

Created on 14/03/16

@author: Carlos Eduardo Barbosa

Copy standard stars according to observing night.

"""
import os
import shutil

from config import *

if __name__ == "__main__":
    wdir = os.path.join(home, "Blanco/dohydra")
    for night in nights:
        os.chdir(os.path.join(wdir, night))
        fits = [x for x in os.listdir(".") if x.endswith("fits")]
        others = ("obj", "comp", "crobj", "flat", "sky", "sflat", "dflat",
                  "zero", "temp", "test", "crsky", "h90", "std")
        fits = [x for x in fits if not x.lower().startswith(others)]
        outdir = os.path.join(home, "data/standards", night)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        for star in fits:
            shutil.copy(star, outdir)

