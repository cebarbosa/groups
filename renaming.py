# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 19:38:17 2013

@author: kadu

This program is used for the renaming of spectra extracted with extracting.py
"""

import os

import pyfits as pf

from config import *

def rename(l1, l2):
    for n1, n2 in zip(l1, l2):
        if not os.path.exists(n1):
            Warning("File {0} does not exists!".format(n1))
            continue
        print "%s --> %s" % (n1, n2)
        os.rename(n1, n2)
        
if __name__ == "__main__":
    wdir = os.path.join(home, "data/reduced")
    for night in nights:
        os.chdir(os.path.join(wdir, night))
        objs = [x for x in os.listdir(".") if x.startswith("crobj") and
                x.endswith(".fits")]
        skies =  [x for x in os.listdir(".") if x.startswith("crobj") and
                 "sky" in x]
        objs = [x for x in objs if x not in skies]
        objs.sort()
        skies.sort()
        multispecs = list(set([x.split("_")[0] for x in objs]))
        multispecs.sort()
        for ms in multispecs:
            specs = [x for x in objs if x.startswith(ms)]
            split_specs = [x.replace('-', '_').split('_') for x in specs]
            groups = list(set([x[1] for x in split_specs]))
            if groups[0] == 'n19':
                for n in split_specs:
                    n[1] = 'ngc193'
                nspecs = ['_'.join(x) for x in split_specs]
                rename(specs, nspecs)
            elif groups[0] == 'h42':
                for n in split_specs:
                    n[1] = 'hcg42'
                nspecs = ['_'.join(x) for x in split_specs]
                rename(specs, nspecs)
            elif groups[0] == 'h22':
                for n in split_specs:
                    n[1] = 'hcg22'
                nspecs = ['_'.join(x) for x in split_specs]
                rename(specs, nspecs)
            elif groups[0] in ['h62']:
                for n in split_specs:
                    n[1] = 'hcg62'
                nspecs = ['_'.join(x) for x in split_specs]
                rename(specs, nspecs)
            elif groups[0] == 'n76':
                for n in split_specs:
                    n[1] = 'ngc7619'
                nspecs = ['_'.join(x) for x in split_specs]
                rename(specs, nspecs)
            if groups[0] in ["hcg22", "hcg42", "ngc193", "ngc7619", "hcg62",
                             "hcg90"]:
                continue
            else:
                nspecs = []
                for spec in specs:
                    new = spec.replace("hcg", "_hcg90")
                    new = new.replace("_eso", "_hcg90_eso")
                    new = new.replace("_466-064", "_hcg90_466-064")
                    new = new.replace("_ngc", "_hcg90_ngc").replace("__", "_")
                    nspecs.append(new)
                rename(specs, nspecs)
    print "Done!"

            
