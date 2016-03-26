# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 12:10:54 2014
Adapted to run in data frou groups of galaxies in Dec 22, 2015

@author: kadu

Run pPXF in data
"""
import os
import pickle

import numpy as np
import pyfits as pf
from scipy.signal import medfilt
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import gridspec
import matplotlib.cm as cm
from scipy.ndimage.filters import convolve1d

from ppxf import ppxf
import ppxf_util as util
from config import *

def wavelength_array(spec, axis=1, extension=0):
    """ Produces array for wavelenght of a given array. """
    w0 = pf.getval(spec, "CRVAL{0}".format(axis), extension)
    deltaw = pf.getval(spec, "CD{0}_{0}".format(axis), extension)
    pix0 = pf.getval(spec, "CRPIX{0}".format(axis), extension)
    npix = pf.getval(spec, "NAXIS{0}".format(axis), extension)
    return w0 + deltaw * (np.arange(npix) + 1 - pix0)

def losvd_convolve(spec, losvd, velscale):
    """ Apply LOSVD to a given spectra given that both wavelength and spec
     arrays are log-binned. """
    # Convert to pixel scale
    pars = np.copy(losvd)
    pars[:2] /= velscale
    dx = int(np.ceil(np.max(abs(pars[0]) + 5*pars[1])))
    nl = 2*dx + 1
    x = np.linspace(-dx, dx, nl)   # Evaluate the Gaussian using steps of 1/factor pixel
    vel = pars[0]
    w = (x - vel)/pars[1]
    w2 = w**2
    gauss = np.exp(-0.5*w2)
    profile = gauss/gauss.sum()
    # Hermite polynomials normalized as in Appendix A of van der Marel & Franx (1993).
    # Coefficients for h5, h6 are given e.g. in Appendix C of Cappellari et al. (2002)
    if losvd.size > 2:        # h_3 h_4
        poly = 1 + pars[2]/np.sqrt(3)*(w*(2*w2-3)) \
                 + pars[3]/np.sqrt(24)*(w2*(4*w2-12)+3)
        if len(losvd) == 6:  # h_5 h_6
            poly += pars[4]/np.sqrt(60)*(w*(w2*(4*w2-20)+15)) \
                  + pars[5]/np.sqrt(720)*(w2*(w2*(8*w2-60)+90)-15)
        profile *= poly
    return convolve1d(spec, profile)
 
def run_ppxf(spectra, velscale, ncomp=None, has_emission=True, mdegree=-1,
             degree=20, plot=False, sky=None, start=None, moments=None):
    """ Run pPXF in a list of spectra"""
    if isinstance(spectra, str):
        spectra = [spectra]
    ##########################################################################
    # Load templates for both stars and gas
    star_templates = pf.getdata(os.path.join(templates_dir,
                                             'miles_FWHM_3.6.fits'), 0)
    logLam2 = pf.getdata(os.path.join(templates_dir, 'miles_FWHM_3.6.fits'), 1)
    miles = np.loadtxt(os.path.join(templates_dir, 'miles_FWHM_3.6.txt'),
                       dtype=str).tolist()
    gas_templates = pf.getdata(os.path.join(templates_dir,
                                             'emission_FWHM_3.6.fits'), 0)
    logLam_gas = pf.getdata(os.path.join(templates_dir, 'emission_FWHM_3.6.fits'),
                            1)
    gas_files = np.loadtxt(os.path.join(templates_dir, 'emission_FWHM_3.6.txt'),
                       dtype=str).tolist()

    ngas = len(gas_files)
    ntemplates = len(miles)
    ##########################################################################
    # Join templates in case emission lines are used.
    if has_emission:
        templates = np.column_stack((star_templates, gas_templates))
        templates_names = np.hstack((miles, gas_files))
    else:
        templates = star_templates
        templates_names = miles
        ngas = 0
    ##########################################################################
    if sky == None:
        nsky = 0
    else:
        nsky = sky.shape[1]
    if ncomp == 1:
        components = 0
    elif ncomp == 2:
        components = np.hstack((np.zeros(len(star_templates[0])),
                                np.ones(len(gas_templates[0]))))
    if moments == None:
        moments = [4] if ncomp == 1 else [4,2]
    for i, spec in enumerate(spectra):
        print "pPXF run of spectrum {0} ({1} of {2})".format(spec, i+1,
              len(spectra))
        plt.clf()
        ######################################################################
        # Read galaxy spectrum and define the wavelength range
        hdu = pf.open(spec)
        spec_lin = hdu[0].data
        h1 = pf.getheader(spec)
        lamRange1 = h1['CRVAL1'] + np.array([0.,h1['CDELT1']*(h1['NAXIS1']-1)])
        ######################################################################
        # Rebin to log scale
        galaxy, logLam1, velscale = util.log_rebin(lamRange1, spec_lin, 
                                                   velscale=velscale)
        ######################################################################
        # First guess for the noise
        noise = np.ones_like(galaxy) * np.std(galaxy - medfilt(galaxy, 5))
        ######################################################################
        # Calculate difference of velocity between spectrum and templates
        # due to different initial wavelength
        dv = (logLam2[0]-logLam1[0])*c
        ######################################################################
        # Set first guess from setup files
        if start == None:
            start = [v0s[spec.split("_")[0]], 100]
        goodPixels = None
        ######################################################################
        # Expand start variable to include multiple components
        if ncomp > 1:
            start = [start, [start[0], 30]]
        ######################################################################
        # First pPXF interaction
        pp0 = ppxf(templates, galaxy, noise, velscale, start,
                   goodpixels=goodPixels, plot=False, moments=moments,
                   degree=12, mdegree=-1, vsyst=dv, component=components,
                   sky=sky)
        rms0 = galaxy[goodPixels] - pp0.bestfit[goodPixels]
        noise0 = 1.4826 * np.median(np.abs(rms0 - np.median(rms0)))
        noise0 = np.zeros_like(galaxy) + noise0
        # Second pPXF interaction, realistic noise estimation
        pp = ppxf(templates, galaxy, noise0, velscale, start,
                  goodpixels=goodPixels, plot=False, moments=moments,
                  degree=degree, mdegree=mdegree, vsyst=dv,
                  component=components, sky=sky)
        plt.title(spec.replace("_", "-"))
        plt.show(block=False)
        plt.savefig("logs/{0}".format(spec.replace(".fits", ".png")))
        ######################################################################
        # Adding other things to the pp object
        pp.template_files = templates_names
        pp.has_emission = has_emission
        pp.dv = dv
        pp.w = np.exp(logLam1)
        pp.velscale = velscale
        pp.spec = spec
        pp.ngas = ngas
        pp.ntemplates = ntemplates
        pp.nsky = nsky
        pp.templates = 0
        ######################################################################
        # Save fit to pickles file to keep session
        ppsave(pp, "logs/{0}".format(spec.replace(".fits", "")))
        pp = ppload("logs/{0}".format(spec.replace(".fits", "")))
        ######################################################################
        ppf = pPXF(spec, velscale, pp)
        logdir = os.path.join(os.getcwd(), "logs")
        if not os.path.exists(logdir):
            os.mkdir(logdir)
        ppf.plot("logs/{0}".format(spec.replace(".fits", ".png")))
        ######################################################################
        # # Save to output text file
        # if ncomp > 1:
        #     pp.sol = pp.sol[0]
        #     pp.error = pp.error[0]
        # sol = [val for pair in zip(pp.sol, pp.error) for val in pair]
        # sol = ["{0:12s}".format("{0:.3g}".format(x)) for x in sol]
        # sol.append("{0:12s}".format("{0:.3g}".format(pp0.chi2)))
        # sol = ["{0:30s}".format(spec)] + sol
    return

def read_setup_file(gal, logw, mask_emline=True):
    """ Read setup file to set first guess and regions to be avoided. """
    w = np.exp(logw)
    filename = gal + ".setup"
    with open(filename) as f:
        f.readline()
        start = f.readline().split()
    start = np.array(start, dtype=float)
    ranges = np.loadtxt(filename, skiprows=5)
    ##########################################################################
    # Mask all the marked regions in the setup file
    if mask_emline:
        for i, (w1, w2) in enumerate(ranges.reshape((len(ranges)/2, 2))):
            if i == 0:
                good = np.where(np.logical_and(w > w1, w < w2))[0]
            else:
                good = np.hstack((good, np.where(np.logical_and(w > w1,
                                                                w < w2))[0]))
        good = np.array(good)
        good.sort()

    ##########################################################################
    # Mask only regions in the beginning and in the end of the spectra plus
    # the residuals in the emission line at 5577 Angstroms
    else:
        ranges = [[np.min(ranges), 5577. - 15], [5577. + 15, np.max(ranges)]]
        for i, (w1, w2) in enumerate(ranges):
            if w1 >= w2:
                continue
            if i == 0:
                good = np.where(np.logical_and(w > w1, w < w2))[0]
            else:
                good = np.hstack((good, np.where(np.logical_and(w > w1,
                                                                w < w2))[0]))
        good = np.array(good)
        good.sort()
        return start, good

def make_table(spectra, output, mc=False, nsim=200, clean=True, pkls=None):
    """ Make table with results.

    ===================
    Input Parameters
    ===================
    spectra : list
        Names of the spectra to be processed. Should end in "fits".
    mc : bool
        Calculate the errors using a Monte Carlo method.
    nsim : int
        Number of simulations in case mc keyword is True.
    clean : bool
        Remove lines for which the velocity dispersion is 1000 km/s.
    pkls : list
        Specify list of pkl files to be used. Default value replaces fits for
        pkl
    ==================
    Output file
    ==================
    In case mc is False, the function produces a file called ppxf_results.dat.
    Otherwise, the name of the file is named ppfx_results_mc_nsim.dat.

    """
    print "Producing summary table..."
    head = ("{0:<30}{1:<14}{2:<14}{3:<14}{4:<14}{5:<14}{6:<14}{7:<14}"
             "{8:<14}{9:<14}{10:<14}{11:<14}{12:<14}{13:<14}\n".format("# FILE",
             "V", "dV", "S", "dS", "h3", "dh3", "h4", "dh4", "chi/DOF",
             "S/N (/ pixel)", "ADEGREE", "MDEGREE", "100*S/N/sigma"))
    results = []
    ##########################################################################
    # Initiate the output file
    ##########################################################################
    with open(output, "w") as f:
        f.write(head)
    ##########################################################################
    if pkls== None:
        pkls = [x.replace(".fits", ".pkl") for x in spectra]
    for i, (spec, pkl) in enumerate(zip(spectra, pkls)):
        print "Working on spectra {0} ({1} / {2})".format(spec, i+1,
                                                          len(spectra))
        if not os.path.exists(spec.replace(".fits", ".pkl")):
            continue
        pp = pPXF(spec, velscale, pklfile=pkl)
        sn = pp.sn
        if mc:
            pp.mc_errors(nsim=nsim)
        if pp.ncomp > 1:
            pp.sol = pp.sol[0]
            pp.error = pp.error[0]
        data = [pp.sol[0], pp.error[0],
                pp.sol[1], pp.error[1], pp.sol[2], pp.error[2], pp.sol[3],
                pp.error[3], pp.chi2, sn]
        data = ["{0:12.3f} ".format(x) for x in data]
        if clean and pp.sol[1] == 1000.:
            comment = "#"
        else:
            comment = ""
        line = ["{0}{1}".format(comment, spec)] + data + \
               ["{0:12}".format(pp.degree), "{0:12}".format(pp.mdegree),
                "{0:12.3f}".format(100 * sn / pp.sol[1])]
        # Append results to outfile
        with open(output, "a") as f:
            f.write("".join(line) + "\n")

class pPXF():
    """ Class to read pPXF pkl files """
    def __init__(self, spec, velscale, pp):
        self.__dict__ = pp.__dict__.copy()
        self.dw = 0.7 # Angstrom / pixel
        self.calc_arrays()
        self.calc_sn()
        return

    def calc_arrays(self):
        """ Calculate the different useful arrays."""
        # Slice matrix into components
        self.m_poly = self.matrix[:,:self.degree + 1]
        self.matrix = self.matrix[:,self.degree + 1:]
        self.m_ssps = self.matrix[:,:self.ntemplates]
        self.matrix = self.matrix[:,self.ntemplates:]
        self.m_gas = self.matrix[:,:self.ngas]
        self.matrix = self.matrix[:,self.ngas:]
        self.m_sky = self.matrix
        # Slice weights
        if hasattr(self, "polyweights"):
            self.w_poly = self.polyweights
            self.poly = self.m_poly.dot(self.w_poly)
        else:
            self.poly = np.zeros_like(self.galaxy)
        if hasattr(self, "mpolyweights"):
            x = np.linspace(-1, 1, len(self.galaxy))
            self.mpoly = np.polynomial.legendre.legval(x, np.append(1, self.mpolyweights))
        else:
            self.mpoly = np.ones_like(self.galaxy)
        self.w_ssps = self.weights[:self.ntemplates]
        self.weights = self.weights[self.ntemplates:]
        self.w_gas = self.weights[:self.ngas]
        self.weights = self.weights[self.ngas:]
        self.w_sky = self.weights
        # Calculating components
        self.ssps = self.m_ssps.dot(self.w_ssps)
        self.gas = self.m_gas.dot(self.w_gas)
        self.bestsky = self.m_sky.dot(self.w_sky)
        return

    def calc_sn(self, w1=4200., w2=6000.):
        idx = np.logical_and(self.w >=w1, self.w <=w2)
        self.res = self.galaxy - self.bestfit
        # Using robust method to calculate noise using median deviation
        self.noise = np.nanstd(self.res[idx])
        self.signal = np.sum(self.mpoly[idx] * (self.ssps[idx] + \
                      self.poly[idx])) / len(self.ssps[idx]) / np.sqrt(self.dw)
        self.sn = self.signal / self.noise
        return

    def mc_errors(self, nsim=200):
        """ Calculate the errors using MC simulations"""
        errs = np.zeros((nsim, len(self.error)))
        for i in range(nsim):
            y = self.bestfit + np.random.normal(0, self.noise,
                                                len(self.galaxy))

            noise = np.ones_like(self.galaxy) * self.noise
            sim = ppxf(self.bestfit_unbroad, y, noise, velscale,
                       [0, self.sol[1]],
                       goodpixels=self.goodpixels, plot=False, moments=4,
                       degree=-1, mdegree=-1,
                       vsyst=self.vsyst, lam=self.lam, quiet=True, bias=0.)
            errs[i] = sim.sol
        median = np.ma.median(errs, axis=0)
        error = 1.4826 * np.ma.median(np.ma.abs(errs - median), axis=0)
        # Here I am using always the maximum error between the simulated
        # and the values given by pPXF.
        self.error = np.maximum(error, self.error)
        return

    def plot(self, output, xlims = [4000, 6500], fignumber=1, textsize = 16):
        """ Plot pPXF run in a output file"""
        if self.ncomp > 1:
            sol = self.sol[0]
            error = self.error[0]
            sol2 = self.sol[1]
            error2 = self.error[1]
        else:
            sol = self.sol
            error = self.error
        plt.figure(fignumber)
        plt.clf()
        gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
        ax = plt.subplot(gs[0])
        ax.minorticks_on()
        ax.plot(self.w, self.galaxy - self.bestsky, "-k", lw=2., label=self.spec.replace("_", \
                                                    "-").replace(".fits", ""))
        ax.plot(self.w, self.bestfit - self.bestsky, "-", lw=2., c="r", label="Bestfit")
        ax.xaxis.set_ticklabels([])
        if self.has_emission:
            ax.plot(self.w[self.goodpixels], self.gas[self.goodpixels], "-b",
                    lw=1., label="Emission Lines")
        # if self.sky != None:
        #     ax.plot(self.w[self.goodpixels], self.bestsky[self.goodpixels], \
        #             "-", lw=1, c="g", label="Sky")
        ax.set_xlim(self.w[0], self.w[-1])
        ax.set_ylim(-5 * self.noise, None)
        leg = plt.legend(loc=6, prop={"size":10})
        leg.get_frame().set_linewidth(0.0)
        plt.axhline(y=0, ls="--", c="k")
        plt.ylabel(r"Flux (Counts)", size=18)
        plt.annotate(r"$\chi^2=${0:.2f}".format(self.chi2),
                     xycoords='axes fraction',
                    xy=(0.05,0.9), size=textsize)
        plt.annotate(r"S/N={0}".format(np.around(self.sn,1)),
                     xycoords='axes fraction', xy=(0.25,0.9), size=textsize)
        plt.annotate(r"V={0} km/s".format(np.around(sol[0])),
                     xycoords='axes fraction', xy=(0.45,0.9), size=textsize,
                     color="r")
        plt.annotate(r"$\sigma$={0} km/s".format(np.around(sol[1])),
                     xycoords='axes fraction', xy=(0.75,0.9), size=textsize,
                     color="r")
        if self.ncomp > 1:
            plt.annotate(r"V={0} km/s".format(np.around(sol2[0])),
                         xycoords='axes fraction', xy=(0.45,0.84),
                         size=textsize, color="b")
            plt.annotate(r"$\sigma$={0} km/s".format(np.around(sol2[1])),
                         xycoords='axes fraction', xy=(0.75,0.84),
                         size=textsize, color="b")
        ax1 = plt.subplot(gs[1])
        ax1.minorticks_on()
        ax1.set_xlim(self.w[0], self.w[-1])
        ax1.plot(self.w[self.goodpixels], (self.galaxy[self.goodpixels] - \
                 self.bestfit[self.goodpixels]), "-k")
        ax1.axhline(y=0, ls="--", c="k")
        ax1.set_ylim(-5 * self.noise, 5 * self.noise)
        ax1.set_xlabel(r"$\lambda$ ($\AA$)", size=18)
        ax1.set_ylabel(r"$\Delta$Flux", size=18)
        gs.update(hspace=0.075, left=0.15, bottom=0.1, top=0.98, right=0.97)
        plt.savefig(output)

    def plot2(self, output, xlims = [4000, 6500], fignumber=1, textsize = 16):
        """ Plot pPXF run in a output file"""
        if self.ncomp > 1:
            sol = self.sol[0]
            error = self.error[0]
            sol2 = self.sol[1]
            error2 = self.error[1]
        else:
            sol = self.sol
            error = self.error
        plt.figure(fignumber)
        plt.clf()
        gs = gridspec.GridSpec(5, 1)
        ax = plt.subplot(gs[0])
        ax.minorticks_on()
        lab0 = "{0}, S/N={1:.1f}".format(self.spec.replace("_","-" \
                                              ).replace(".fits", ""), self.sn)
        ax.plot(self.w, self.galaxy, "-k", lw=2., label=lab0)
        leg = plt.legend(loc=0, prop={"size":10})
        leg.get_frame().set_linewidth(0.0)
        ax1 = plt.subplot(gs[1])
        ax1.minorticks_on()
        ax1.plot(self.w[self.goodpixels], self.ssps[self.goodpixels]+self.poly,
                "-", lw=1, c="r", label="SSPs")
        ax1.xaxis.set_ticklabels([])
        leg = plt.legend(loc=0, prop={"size":10})
        leg.get_frame().set_linewidth(0.0)
        if self.has_emission:
            ax2 = plt.subplot(gs[2])
            ax2.plot(self.w[self.goodpixels], self.gas[self.goodpixels], "-b",
                    lw=1., label="Emission Lines")
            ax2.xaxis.set_ticklabels([])
            ax2.set_xlim(self.w[0], self.w[-1])
            leg = plt.legend(loc=0, prop={"size":10})
            leg.get_frame().set_linewidth(0.0)
        if self.sky != None:
            ax3 = plt.subplot(gs[3])
            ax3.plot(self.w[self.goodpixels], self.bestsky[self.goodpixels], \
                    "-", lw=1, c="g", label="Sky")
            ax3.xaxis.set_ticklabels([])
            ax3.set_xlim(self.w[0], self.w[-1])
            leg = plt.legend(loc=0, prop={"size":10})
            leg.get_frame().set_linewidth(0.0)
        ax4 = plt.subplot(gs[4])
        ax4.minorticks_on()
        ax.set_xlim(self.w[0], self.w[-1])
        ax1.set_xlim(self.w[0], self.w[-1])
        ax4.plot(self.w[self.goodpixels], self.galaxy[self.goodpixels] - \
                 self.bestfit[self.goodpixels], "-k")
        ax4.axhline(y=0, ls="--", c="k")
        ax4.axhline(y=self.noise, ls="--", c="r")
        ax4.axhline(y=-self.noise, ls="--", c="r")
        ax4.set_xlabel(r"$\lambda$ ($\AA$)", size=18)
        ax.xaxis.set_ticklabels([])
        ax4.set_ylabel(r"$\Delta$Flux", size=18)
        plt.axhline(y=0, ls="--", c="k")
        plt.ylabel(r"Flux (Counts)", size=18)
        plt.tight_layout()
        plt.annotate(r"V={0} km/s".format(np.around(sol[0])),
                     xycoords='axes fraction', xy=(0.45,0.9), size=textsize,
                     color="r")
        plt.annotate(r"$\sigma$={0} km/s".format(np.around(sol[1])),
                     xycoords='axes fraction', xy=(0.75,0.9), size=textsize,
                     color="r")
        if self.ncomp > 1:
            plt.annotate(r"V={0} km/s".format(np.around(sol2[0])),
                         xycoords='axes fraction', xy=(0.45,0.88),
                         size=textsize, color="b")
            plt.annotate(r"$\sigma$={0} km/s".format(np.around(sol2[1])),
                         xycoords='axes fraction', xy=(0.75,0.88),
                         size=textsize, color="b")
        gs.update(hspace=0.2, left=0.15, bottom=0.1, top=0.98, right=0.97)
        plt.savefig(output)

def fits_to_matrix(filenames):
    """ Load several fits of the same dimension into a matrix.  """
    f = pf.getdata(filenames[0])
    a = np.zeros((f.shape[0], len(filenames)))
    for i in range(len(filenames)):
        try:
            a[:,i] = pf.getdata(filenames[i])
        except:
            print filenames[i]
    return a

def load_sky(filenames, velscale, full_output=False):
    """ Load and rebin sky files. """
    skydata = fits_to_matrix(filenames)
    h1 = pf.getheader(filenames[0])
    lamRange1 = h1['CRVAL1'] + np.array([0.,h1['CDELT1']*(h1['NAXIS1']-1)])
    sky1, logLam1, velscale = util.log_rebin(lamRange1, skydata[:,0],
                                                   velscale=velscale)

    skylog = np.zeros((sky1.shape[0], len(filenames)))
    for i in range(len(filenames)):
        skylog[:,i], logLam1, velscale = util.log_rebin(lamRange1,
                                        skydata[:,i], velscale=velscale)
    if full_output:
        return skylog, logLam1
    return skylog

def make_table_from_txt():
    """ Make a summary table using the txt outputs. """
    nights = sorted(os.listdir(data_dir))
    for night in nights:
        wdir = os.path.join(data_dir, night)
        os.chdir(wdir)
        specs = [x for x in os.listdir(".") if x.endswith(".fits")]
        txts = [x.replace(".fits", ".txt") for x in specs]
        txts = [x for x in txts if  os.path.exists(os.path.join(wdir, x))]
        txts.sort()
        with open("ppxf_results.txt", "w") as fout:
            for txt in txts:
                with open(txt) as fin:
                    fout.write(fin.read() + "\n")
    return

def select_specs():
    """ Select which spectra can be used in the analysis. """
    nights = sorted(os.listdir(data_dir))
    for night in nights:
        wdir = os.path.join(data_dir, night, "logs")
        os.chdir(wdir)
        pngs = sorted([x for x in os.listdir(".") if x.endswith(".png")])
        comments = []
        for image in pngs:
            img = mpimg.imread(image)
            plt.imshow(img)
            plt.axis("off")
            plt.pause(0.001)
            plt.show()
            comm = raw_input("Skip this spectrum in the anaylis? (y/N) ")
            if comm.lower().strip() in ["y", "ye", "yes"]:
                comments.append("#")
            else:
                comments.append("")
            plt.clf()
        ignorelist = ["{0}{1}".format(x,y.replace(".png", ".fits")) for \
                      x,y in zip(comments, pngs)]
        output = os.path.join(data_dir, night, "ignore.txt")
        with open(output, "w") as f:
            f.write("\n".join(ignorelist))

def run_over_all():
    """ Run pPXF in a generic way over all data. """
    nights = sorted(os.listdir(data_dir))
    for night in nights:
        print "Working in run ", night
        wdir = os.path.join(data_dir, night)
        os.chdir(wdir)
        log_dir = os.path.join(wdir, "logs")
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)
        fits = [x for x in os.listdir(".") if x.endswith(".fits")]
        skies =  [x for x in fits if x.startswith("sky")]
        specs = [x for x in fits if x not in skies]
        specs.sort()
        skies.sort()
        sky = load_sky(skies, velscale)
        # #################################################################
        # # Go to the main routine of fitting
        run_ppxf(specs, velscale, ncomp=1, has_emission=1, mdegree=-1,
                 degree=12, plot=True, sky=sky)
    return

def ppsave(pp, outroot="logs/out"):
    """ Produces output files for a ppxf object. """
    arrays = ["matrix", "w", "bestfit", "goodpixels", "galaxy", "noise"]
    # delattr(pp, "star_rfft")
    # delattr(pp, "star")
    hdus = []
    for i,att in enumerate(arrays):
        if i == 0:
            hdus.append(pf.PrimaryHDU(getattr(pp, att)))
        else:
            hdus.append(pf.ImageHDU(getattr(pp, att)))
        delattr(pp, att)
    hdulist = pf.HDUList(hdus)
    hdulist.writeto(outroot + ".fits", clobber=True)
    with open(outroot + ".pkl", "w") as f:
        pickle.dump(pp, f)

def ppload(inroot="logs/out"):
    """ Read ppxf arrays. """
    with open(inroot + ".pkl") as f:
        pp = pickle.load(f)
    arrays = ["matrix", "w", "bestfit", "goodpixels", "galaxy", "noise"]
    for i, item in enumerate(arrays):
        setattr(pp, item, pf.getdata(inroot + ".fits", i))
    return pp

def plot_all():
    """ Make plot of all fits. """
    os.chdir(data_dir)
    specs = sorted([x for x in os.listdir(".") if x.endswith(".fits")])
    for i,spec in enumerate(specs):
        print "Working on spec {0} ({1}/{2})".format(spec, i+1, len(specs))
        pp = ppload("logs/{0}".format(spec.replace(".fits", "")))
        pp = pPXF(spec, velscale, pp)
        pp.plot("logs/{0}".format(spec.replace(".fits", ".png")))

def run_list(night, specs, start=None):
    """ Run pPXF on a given list of spectra of the same night. """
    wdir = os.path.join(data_dir, night)
    if start is None:
        start = [2000., 50.]
    os.chdir(wdir)
    fits = [x for x in os.listdir(".") if x.endswith(".fits")]
    skies =  sorted([x for x in fits if x.startswith("sky")])
    sky, loglam = load_sky(skies, velscale, full_output=True)
    # make_sky_fig(sky, loglam, skies)]
    run_ppxf(specs, velscale, ncomp=2, has_emission=1, mdegree=-1,
                 degree=12, plot=True, sky=sky, start=start,
             moments=[4, 2])

def make_sky_fig(skies, loglam, filenames):
    w = np.exp(loglam)
    fig = plt.figure(1)
    ax = plt.subplot(111)
    for i, sky in enumerate(skies.T):
        ax.plot(w, sky / np.median(sky), "-", label=filenames[i])
    ax.set_ylim(0,20)
    plt.legend(loc=0, prop={"size":6})
    plt.savefig("logs/sky.png")
    return

def run_stellar_templates(velscale):
    """ Run over stellar templates. """
    temp_dir = os.path.join(home, "stellar_templates")
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
            print standard
            idx = ids.index(name)
            T, logg, FeH = star_pars[idx]
            tempfile= "MILES_Teff{0:.2f}_Logg{1:.2f}_MH{2:.2f}" \
                       "_linear_FWHM_3.6.fits".format(T, logg, FeH )
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
            pp0 = ppxf(temp, star, noise, velscale, [300.,20], plot=False,
                       moments=2, degree=-1, mdegree=25, vsyst=dv)
            plt.title("{0} {1}".format(night, standard))
            # plt.show()
            plt.clf()
            ######################################################################
            pp0.w = np.exp(logLam1)
            pp0.wtemp = np.exp(logLam2)
            pp0.template_linear = [wtemp, template]
            pp0.temp = temp
            pp0.ntemplates = 1
            pp0.ngas = 0
            if not os.path.exists("logs"):
                os.mkdir("logs")
            ppsave(pp0, "logs/{0}".format(standard.replace(".fits", "")))
    return

def flux_calibration_test(velscale):
    standards_dir = os.path.join(home, "data/standards")
    for night in nights:
        os.chdir(os.path.join(standards_dir, night))
        standards = sorted([x for x in os.listdir(".") if
                            x.endswith(".fits")])
        standards = [x for x in standards if
                     os.path.exists("logs/{0}".format(x))]
        fibers = np.array([int(x.split(".")[1]) for x in standards])
        cmap = cm.get_cmap("rainbow")
        color = np.linspace(0,1,len(nights))
        for i,standard in enumerate(standards):
            pp = ppload("logs/{0}".format(standard.replace(".fits", "")))
            pp = pPXF(standard, velscale, pp)
            h = pf.getheader(standard)
            # plt.plot(pp.w, pp.galaxy, "-k")
            # plt.plot(pp.w, pp.bestfit, "-r")
            lab = night if i == 0 else None
            plt.plot(pp.w, pp.mpoly, "-", c=cmap(color[nights.index(night)]),
                         label=lab)
    plt.legend(loc=0)
    plt.show()
    return

def run_candidates(velscale, filenames=None):
    """ Run pPXF over candidates. """
    os.chdir(data_dir)
    log_dir = os.path.join(data_dir, "logs")
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    if filenames is None:
        filenames = [x for x in os.listdir(".") if x.endswith(".fits")]
    for night in nights:
        print night
        nspecs = sorted([x for x in os.listdir(".") if \
                x.endswith("{0}.fits".format(night))])
        specs = [x for x in filenames if x in nspecs]
        if len(specs) == 0:
            continue
        os.chdir(os.path.join(home, "data/combined", night))
        skies =  sorted([x for x in os.listdir(".") if x.startswith("sky") and
                         x.endswith(".fits")])
        specs.sort()
        skies.sort()
        sky = load_sky(skies, velscale)
        os.chdir(data_dir)
        # #################################################################
        # # Go to the main routine of fitting
        run_ppxf(specs, velscale, ncomp=1, has_emission=False, mdegree=12,
                 degree=-1, plot=True, sky=sky, start=[4334., 30])
    return

def make_table():
    """ Make table with results. """
    head = ("{0:<30}{1:<14}{2:<14}{3:<14}{4:<14}{5:<14}{6:<14}{7:<14}"
             "{8:<14}{9:<14}{10:<14}{11:<14}{12:<14}\n".format("# FILE",
             "V", "dV", "S", "dS", "h3", "dh3", "h4", "dh4", "chi/DOF",
             "S/N (/ pixel)", "ADEGREE", "MDEGREE"))
    os.chdir(data_dir)
    specs = sorted([x for x in os.listdir(".") if x.endswith(".fits")])
    results = []
    for spec in specs:
        print spec
        vhelio = pf.getval(spec, "VHELIO")
        output = os.path.join(data_dir, "ppxf_results.dat")
        pp = ppload("logs/" + spec.replace(".fits", ""))
        pp = pPXF(spec, velscale, pp)
        sol = pp.sol if pp.ncomp == 1 else pp.sol[0]
        sol[0] += vhelio
        error = pp.error if pp.ncomp == 1 else pp.error[0]
        cond = error[1] > 300. and pp.sn > 10.
        name = spec if cond else "#{0}".format(spec)
        line = np.zeros((sol.size + error.size,))
        line[0::2] = sol
        line[1::2] = error
        line = np.append(line, [pp.chi2, pp.sn])
        line = ["{0:12.3f}".format(x) for x in line]
        line = ["{0:30s}".format(name)] + line + \
               ["{0:12}".format(pp.degree), "{0:12}".format(pp.mdegree)]
        results.append("".join(line))
    # Append results to outfile
    with open(output, "w") as f:
        f.write(head)
        f.write("\n".join(results))

if __name__ == '__main__':
    # run_stellar_templates(velscale)
    # flux_calibration_test(velscale)
    # run_list("blanco10n1", ["hcg62_14.fits"])
    # run_over_all()
    run_candidates(velscale, filenames=["hcg62_7_blanco10n2.fits"])
    # plot_all()
    # make_table()