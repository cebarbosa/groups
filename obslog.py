# -*- coding: utf-8 -*-
"""

Created on 02/05/16

@author: Carlos Eduardo Barbosa

Read fits files to create and observational log.

"""
import os
import datetime

import numpy as np
import pyfits as pf

from config import *

if __name__ == "__main__":
    wdir = os.path.join(home, "Blanco/dohydra")
    objs, dates, exptimes = [], [], []
    for night in nights:
        os.chdir(os.path.join(wdir, night))
        filenames = sorted([x for x in os.listdir(".") if
                         x.startswith("crobj") and not x.endswith("ms.fits")])
        for fits in filenames:
            obj = pf.getval(fits, "object")
            if obj.lower().startswith(("hr", "hd", "std")):

                continue
            date = pf.getval(fits, "date-obs").split("T")[0]
            exptime = pf.getval(fits, "exptime")
            if exptime < 100:
                continue
            objs.append(obj)
            dates.append(date)
            exptimes.append(exptime)
    objs = np.array(objs)
    dates = np.array(dates)
    exptimes = np.array(exptimes)
    times, masks, obsdate = [], [], []
    for obj in set(objs):
        masks.append(obj)
        idx = np.where(objs == obj)
        expt = []
        for exptime in set(exptimes[idx]):
            t = len(np.where(exptimes[idx] == exptime)[0])
            expt.append("{0} x {1:.0f}".format(t, exptime))
        times.append(", ".join(expt))
        obsdate.append(list(set(dates[idx]))[0])
    data = zip(obsdate, masks, times)
    data.sort(key = lambda t: t[0])
    data = np.array(data)
    lines = []
    for i, dd in enumerate(sorted(set(obsdate))):
        dt = dd.split("-")
        dt = datetime.datetime(*[int(x) for x in dd.split("-")])
        date = dt.strftime('%b %d %Y')
        idx = np.where(dd == data[:,0])[0]
        # date = "\multirow{{{1}}}{{*}}{{{0}}}".format(date, len(idx))
        line = ["({0})".format(i+1)]
        for k,j in enumerate(idx):
            if k == 0:
                line.append(date)
            else:
                line.append("")
            line.append(data[j,1])
            line.append("{0}s\\\\\n".format(data[j,2]))

        lines.append(" & ".join(line))
    print "".join(lines)

