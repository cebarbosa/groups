# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 12:10:54 2014
Adapted to run in data frou groups of galaxies in Dec 22, 2015

@author: kadu

Run pPXF in data
"""
import os
import pickle
import fileinput

import numpy as np
import pyfits as pf
from scipy import ndimage
from scipy.signal import medfilt
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from ppxf import ppxf
import ppxf_util as util
from config import *
from load_templates import stellar_templates, emission_templates, \
                            wavelength_array
 
def run_ppxf(spectra, velscale, ncomp=None, has_emission=True, mdegree=-1,
             degree=20, pkls=None, plot=False, sky=None, lamRange1=None):
    """ Run pPXF in a list of spectra"""
    if isinstance(spectra, str):
        spectra = [spectra]
    if isinstance(pkls, str):
        pkls = [pkls]
    if pkls == None:
        pkls = [x.replace(".fits", ".pkl") for x in spectra]
    ##########################################################################
    # Load templates for both stars and gas
    star_templates, logLam2, delta, miles= stellar_templates(velscale)
    gas_templates,logLam_gas, delta_gas, gas_files=emission_templates(velscale)
    ##########################################################################
    # Join templates in case emission lines are used.
    if has_emission:
        templates = np.column_stack((star_templates, gas_templates))
        templates_names = np.hstack((miles, gas_files))
    else:
        templates = star_templates
        templates_names = miles
    ##########################################################################
    if ncomp == 1:
        components = 0
        moments = [4]
    elif ncomp == 2:
        components = np.hstack((np.zeros(len(star_templates[0])),
                                np.ones(len(gas_templates[0]))))
        moments = [4,2]
    for i, spec in enumerate(spectra):
        print "pPXF run of spectrum {0} ({1} of {2})".format(spec, i+1,
              len(spectra))
        pkl = pkls[i]
        plt.clf()
        ######################################################################
        # Read one galaxy spectrum and define the wavelength range
        specfile = os.path.join(wdir, spec)
        hdu = pf.open(specfile)
        spec_lin = hdu[0].data
        h1 = pf.getheader(specfile)
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
        start = [v0s[spec.split("_")[0]], 100]
        goodPixels = None
        ######################################################################
        # Expand start variable to include multiple components
        if ncomp > 1:
            start = [start, [start[0], 30]]
        ######################################################################
        # First pPXF interaction
        if os.path.exists(spec.replace(".fits", ".pkl")):
            pp0 = pPXF(spec, velscale, pklfile=spec.replace(".fits", ".pkl"))
            noise0 = pp0.noise
        else:
            pp0 = ppxf(templates, galaxy, noise, velscale, start,
                       goodpixels=goodPixels, plot=False, moments=moments,
                       degree=12, mdegree=-1, vsyst=dv, component=components,
                       sky=sky)
            rms0 = galaxy[goodPixels] - pp0.bestfit[goodPixels]
            noise0 = 1.4826 * np.median(np.abs(rms0 - np.median(rms0)))
            noise0 = np.zeros_like(galaxy) + noise0
        # Second pPXF interaction, realistic noise estimation
        pp = ppxf(templates, galaxy, noise0, velscale, start,
                  goodpixels=goodPixels, plot=plot, moments=moments,
                  degree=degree, mdegree=mdegree, vsyst=dv,
                  component=components, sky=sky)
        plt.title(spec.replace("_", "-"))
        plt.show(block=False)
        plt.savefig("logs/{0}".format(spec.replace(".fits", ".png")))
        pp.template_files = templates_names
        pp.has_emission = has_emission
        ######################################################################
        ppf = pPXF(spec, velscale, pp)
        logdir = os.path.join(os.getcwd(), "logs")
        if not os.path.exists(logdir):
            os.mkdir(logdir)
        ppf.plot("logs/{0}".format(spec.replace(".fits", ".png")))
        ######################################################################
        # # Save to output file to keep session
        sol = [val for pair in zip(pp.sol, pp.error) for val in pair]
        sol = ["{0:12s}".format("{0:.3g}".format(x)) for x in sol]
        sol.append("{0:12s}".format("{0:.3g}".format(pp0.chi2)))
        sol = ["{0:30s}".format(spec)] + sol
        with open(pkl.replace(".pkl", ".txt"), "w") as f:
            f.write("".join(sol))
        ######################################################################
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
        pp.calc_sn()
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
        self.spec = spec
        self.velscale = velscale
        self.w = wavelength_array(spec)
        self.flux = pf.getdata(self.spec)
        self.flux_log, self.logw, velscale = util.log_rebin(
                        [self.w[0], self.w[-1]], self.flux, velscale=velscale)
        self.w_log = np.exp(self.logw)
        self.header = pf.getheader(self.spec)
        self.lam = self.header['CRVAL1'] + np.array([0.,
                              self.header['CDELT1']*(self.header['NAXIS1']-1)])
        ######################################################################
        # # Read templates
        star_templates, self.logLam2, self.delta, miles= stellar_templates(velscale)
        ######################################################################
        # Convolve our spectra to match MILES resolution
        FWHM_dif = np.sqrt(FWHM_tem**2 - FWHM_spec**2)
        sigma = FWHM_dif/2.355/self.delta # Sigma difference in pixels
        ######################################################################
        spec_lin = ndimage.gaussian_filter1d(self.flux,sigma)
        # Rebin to logarithm scale
        galaxy, self.logLam1, velscale = util.log_rebin(self.lam, spec_lin,
                                                   velscale=velscale)
        self.dv = (self.logLam2[0]-self.logLam1[0])*c
        return

    def calc_sn(self, w1=5200., w2=5500.):
        idx = np.logical_and(self.w >=w1, self.w <=w2)
        self.res = self.galaxy[idx] - self.bestfit[idx]
        # Using robust method to calculate noise using median deviation
        self.noise = 1.4826 * np.median(np.abs(self.res - np.median(self.res)))
        self.signal = np.sum(self.galaxy[idx]) / len(self.galaxy[idx])
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

    def calc_arrays_emission(self):
        """ Calculate arrays correcting for emission lines. """
        if self.has_emission:
            if self.sky == None:
                em_weights = self.weights[-3:]
                em_matrix = self.matrix[:,-3:]
            else:
                em_weights = self.weights[-4:-1]
                em_matrix = self.matrix[:,-4:-1]
            self.em = em_matrix.dot(em_weights)
            f = interp1d(self.w_log, self.em, kind="linear", bounds_error=False,
                         fill_value=0. )
            self.em_linear = f(self.w)
        else:
            self.em_linear = np.zeros_like(self.flux)
            self.em = np.zeros_like(self.bestfit)
        return

    def sky_sub(self):
        """ Make sky subtraction in case it was not done before. """
        if self.sky != None:
            ntemplates = self.star.shape[1]
            self.bestsky = self.sky.dot(self.weights[ntemplates:])
            self.galaxy-= self.bestsky
            self.bestfit -= self.bestsky
            # self.bestfit_unbroad -= self.bestsky
            f = interp1d(self.w_log, self.sky.T[0] * self.weights[-1],
                         kind="linear",
                         bounds_error=False, fill_value=0. )
            self.sky_lin = f(self.w)
            self.flux -= self.sky_lin
        else:
            print "Warning: No sky templates for this run."

    def plot(self, output, xlims = [4000, 6500], fignumber=1, textsize = 16):
        """ Plot pPXF run in a output file"""
        plt.figure(fignumber)
        plt.clf()
        ax = plt.subplot(111)
        ax.minorticks_on()
        self.calc_sn()
        self.calc_arrays_emission()
        if self.ncomp > 1:
            sol = self.sol[0]
            error = self.error[0]
            sol2 = self.sol[1]
            error2 = self.error[1]
        else:
            sol = self.sol
            error = self.error
        self.sky_sub()
        ax.plot(self.w_log, self.galaxy, "-k")
        ax.plot(self.w_log[self.goodpixels], self.bestfit[self.goodpixels],
                "-r", lw=1.5)
        if self.has_emission:
            ax.plot(self.w_log[self.goodpixels],
                     self.bestfit[self.goodpixels] - self.em[self.goodpixels],
                    "--y")
            ax.plot(self.w_log[self.goodpixels], self.em[self.goodpixels], "-b",
                    lw=1.5)
            # plt.plot(self.w, self.flux - self.em_linear, "--y")
        diff = self.galaxy[self.goodpixels] - self.bestfit[self.goodpixels]
        plt.plot(self.w_log[self.goodpixels], diff, ".g", ms=0.5)
        badpixels = np.setdiff1d(np.arange(len((self.w_log))), self.goodpixels)
        badpixels.sort()
        ax.set_xlim(xlims[0], xlims[1])
        plt.plot(self.w_log[badpixels], 
                 self.flux_log[badpixels] - self.bestfit[badpixels], 
                 ".k", ms=0.5)
        ax.set_ylim(-3 * self.noise, 4 * np.median(self.galaxy))
        plt.axhline(y=0, ls="--", c="k")
        plt.xlabel(r"$\lambda$ ($\AA$)", size=18)
        plt.ylabel(r"Flux (Counts)", size=18)
        plt.tight_layout()
        plt.annotate(r"$\chi^2=${0:.2f}".format(self.chi2),
                     xycoords='axes fraction',
                    xy=(0.05,0.88), size=textsize)
        plt.annotate(r"S/N={0}".format(np.around(self.sn,1)),
                     xycoords='axes fraction', xy=(0.25,0.95), size=textsize)
        plt.annotate(r"V={0} km/s".format(np.around(sol[0])),
                     xycoords='axes fraction', xy=(0.45,0.95), size=textsize,
                     color="r")
        plt.annotate(r"$\sigma$={0} km/s".format(np.around(sol[1])),
                     xycoords='axes fraction', xy=(0.75,0.95), size=textsize,
                     color="r")
        if self.ncomp > 1:
            plt.annotate(r"V={0} km/s".format(np.around(sol2[0])),
                         xycoords='axes fraction', xy=(0.45,0.88),
                         size=textsize, color="b")
            plt.annotate(r"$\sigma$={0} km/s".format(np.around(sol2[1])),
                         xycoords='axes fraction', xy=(0.75,0.88),
                         size=textsize, color="b")
        plt.savefig(output)




def read_sky(filenames):
    """ Prepare matrix with linear sky files. """
    f = pf.getdata(filenames[0])
    a = np.zeros((f.shape[0], len(filenames)))
    for i in range(len(filenames)):
        try:
            a[:,i] = pf.getdata(filenames[i])
        except:
            print filenames[i]
    return a

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
        spec1 = pf.getdata(specs[0])
        h1 = pf.getheader(specs[0])
        lamRange1 = h1['CRVAL1'] + np.array([0.,h1['CDELT1']*(h1['NAXIS1']-1)])
        galaxy1, logLam1, velscale = util.log_rebin(lamRange1, spec1,
                                                   velscale=velscale)
        skydata = read_sky(skies)
        skylog = np.zeros((galaxy1.shape[0], len(skies)))
        for i in range(len(skies)):
            skylog[:,i], logLam1, velscale = util.log_rebin(lamRange1,
                                            skydata[:,i], velscale=velscale)
        # #################################################################
        # # Go to the main routine of fitting
        run_ppxf(specs, velscale, ncomp=1, has_emission=1, mdegree=-1,
                 degree=12, plot=True, sky=skylog)
    return

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

if __name__ == '__main__':
    select_specs()

