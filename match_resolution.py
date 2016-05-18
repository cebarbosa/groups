# -*- coding: utf-8 -*-
"""

Created on 29/03/16

@author: Carlos Eduardo Barbosa

Match resolution of standard stars

"""
from run_ppxf import *

def match_resolution(velscale):
    """ Run over stellar templates. """
    temp_dir = os.path.join(home, "miles_models")
    standards_dir = os.path.join(home, "data/standards")
    table = os.path.join(tables_dir, "lick_standards.txt")
    ids = np.loadtxt(table, usecols=(0,), dtype=str).tolist()
    star_pars = np.loadtxt(table, usecols=(26,27,28,))
    results = []
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
            if not os.path.exists(tempfile):
                os.chdir(cdir)
                continue
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
            ##################################################################
            # First run to set errors
            idx0 = np.where((4000.<np.exp(logLam1)) & (np.exp(logLam1)<5500))
            star0 = star[idx0]
            logLam0 = logLam1[idx0]
            dv0 = (logLam2[0]-logLam0[0])*c
            noise0 = np.ones_like(star0)
            pp0 = ppxf(temp, star0, noise0, velscale, [0.,5], plot=False,
                       moments=2, degree=20, mdegree=-1, vsyst=dv0, quiet=True)
            noise = np.std(pp0.bestfit - pp0.galaxy)
            ##################################################################
            deltaw = 300
            wtest = np.arange(4000., 6300, 100)
            for w1 in wtest:
                w2 = w1 + deltaw
                idx1 = np.where((w1<np.exp(logLam1)) & (np.exp(logLam1)<w2))
                idx2 = np.where((w1-50<np.exp(logLam2)) & (np.exp(logLam2)<w2+50))
                star1 = star[idx1]
                logLam3 = logLam1[idx1]
                logLam4 = logLam1[idx2]
                temp1 = temp[idx2]
                dv1 = (logLam2[idx2][0]-logLam3[0])*c
                noise1 = np.ones_like(star1) * noise
                pp1= ppxf(temp1, star1, noise1, velscale, [0.,5], plot=False,
                          moments=2, degree=20, mdegree=-1, vsyst=dv1,
                          quiet=True)
                print w1 + 0.5 * deltaw, pp1.sol[1], pp1.error[1]
                results.append([w1 + 0.5 * deltaw, pp1.sol[1], pp1.error[1]])
    results = np.array(results)
    with open(os.path.join(tables_dir, "w_sig_standard.dat"),"w") as f:
        np.savetxt(f, results)
    return

def plot():
    filename = os.path.join(tables_dir, "w_sig_standard.dat")
    sig2fwhm = 2.335
    wave, sigma, err = np.loadtxt(filename).T
    fig, ax = plt.subplots()
    ax.minorticks_on()
    ws, ms = [], []
    for w in np.unique(wave):
        label = "Standard Stars" if w==wave[0] else None
        idx = np.where(wave==w)
        ax.plot(wave[idx], sigma[idx], ".", color="0.8", label=label,
                marker=(5,1,0), ms=10, mec="0.6")
        ws.append(w)
        ms.append(np.median(sigma[idx]))
    ax.plot(ws, ms, "sr", mec="r", label="Median Velocity\nDispersion", ms=8)
    obsres = np.array(ws) * np.array(ms) / c * sig2fwhm
    res = np.sqrt(obsres**2 + 2.5**2)
    resolution = ws / res
    print np.median(resolution)
    print resolution.min(), resolution.max()
    z = np.polyfit(ws, res, 5)
    p = np.poly1d(z)
    w = np.linspace(4000., 6500, 2500)
    # ax.plot(w, p(w), "-r", label="Best fit")
    ax.set_ylim(0,220)
    ax.set_xlabel("Wavelength (\AA)")
    ax.set_ylabel("$\sigma$ (km/s)", color="r")
    plt.legend(loc=0, frameon=False,prop={'size':16})
    ax2 = ax.twinx()
    ax2.minorticks_on()
    ax2.plot(ws, res, "ob", ms=8, mec="none", label="Hydra-CTIO\n resolution")
    ax2.plot(w,p(w),"-b", label="Polynomial Fit")
    ax2.set_ylabel("FWHM (\AA)", color="b")
    for tl in ax2.get_yticklabels():
        tl.set_color('b')
    for tl in ax.get_yticklabels():
        tl.set_color('r')
    plt.legend(loc=0, frameon=False,prop={'size':16})
    plt.subplots_adjust(left=0.08, right=0.92, top=0.98)
    plt.savefig("/home/kadu/projects/thesis/figs/hydra_resolution.png",
                dpi=200)
    np.savetxt(os.path.join(tables_dir, "wave_fwhm_standards.dat"),
               np.column_stack((w, p(w))))
    return

if __name__ == "__main__":
    velscale = 3.
    # match_resolution(velscale)
    plot()
