# -*- coding: utf-8 -*-
"""

Created on 29/03/16

@author: Carlos Eduardo Barbosa

Match resolution of standard stars

"""
from run_ppxf import *

def match_resolution(velscale):
    """ Run over stellar templates. """
    temp_dir = os.path.join(home, "stellar_templates/MILES_FWHM_2.5")
    standards_dir = os.path.join(home, "data/standards")
    table = os.path.join(tables_dir, "lick_standards.txt")
    ids = np.loadtxt(table, usecols=(0,), dtype=str).tolist()
    star_pars = np.loadtxt(table, usecols=(26,27,28,))
    for night in nights:
        cdir = os.path.join(standards_dir, night)
        os.chdir(cdir)
        standards = sorted([x for x in os.listdir(".") if x.endswith(".fits")])
        for standard in standards:
            name = standard.split(".")[0].upper()
            if name not in ids:
                continue
            idx = ids.index(name)
            T, logg, FeH = star_pars[idx]
            tempfile= "MILES_Teff{0:.2f}_Logg{1:.2f}_MH{2:.2f}" \
                       "_linear_FWHM_2.50.fits".format(T, logg, FeH )
            os.chdir(temp_dir)
            template = pf.getdata(tempfile)
            htemp = pf.getheader(tempfile)
            wtemp = htemp["CRVAL1"] + htemp["CDELT1"] * \
                            (np.arange(htemp["NAXIS1"]) + 1 - htemp["CRPIX1"])
            os.chdir(cdir)
            data = pf.getdata(standard)
            w = wavelength_array(standard)
            lamRange1 = np.array([w[0], w[-1]])
            lamRange2 = np.array([wtemp[0], wtemp[-1]])
            # Rebin to log scale
            star, logLam1, velscale = util.log_rebin(lamRange1, data,
                                                       velscale=velscale)
            temp, logLam2, velscale = util.log_rebin(lamRange2, template,
                                                       velscale=velscale)
            noise = np.ones_like(star)
            dv = (logLam2[0]-logLam1[0])*c
            pp0 = ppxf(temp, star, noise, velscale, [300.,5], plot=False,
                       moments=2, degree=20, mdegree=-1, vsyst=dv, quiet=True)
            noise = np.ones_like(noise) * np.nanstd(star - pp0.bestfit)
            pp0 = ppxf(temp, star, noise, velscale, [0.,5], plot=False,
                       moments=2, degree=20, mdegree=-1, vsyst=dv)
            pp0.w = np.exp(logLam1)
            pp0.wtemp = np.exp(logLam2)
            pp0.template_linear = [wtemp, template]
            pp0.temp = temp
            pp0.ntemplates = 1
            pp0.ngas = 0
            pp0.has_emission = False
            pp0.dv = dv
            pp0.velscale = velscale
            pp0.ngas = 0
            pp0.nsky = 0
            if not os.path.exists("logs"):
                os.mkdir("logs")
            ppsave(pp0, "logs/{0}".format(standard.replace(".fits", "")))
            pp = ppload("logs/{0}".format(standard.replace(".fits", "")))
            pp = pPXF(standard, velscale, pp)
            pp.plot("logs/{0}".format(standard.replace(".fits", ".png")))
    return

if __name__ == "__main__":
    match_resolution(velscale)
