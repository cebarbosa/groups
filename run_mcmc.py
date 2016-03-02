#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 2016

@author: cbarbosa
"""

import os

import numpy as np
from scipy.interpolate import LinearNDInterpolator
import pymc
import multiprocessing as mp

from config import *

class SSP:
    """ Interface to interpolate and call different models. """
    def __init__(self, modelname, idx=None):
        self.modelname = modelname
        self.idx = idx if idx != None else np.arange(25)
        self.load_model()
        self.calc_lims()

    def load_model(self):
        """ Load data for model and set appropriate indices. """
        if self.modelname in ["TMJ10ext", "TMJ10", "TMJ10Padova"]:
            tablename = dict([("TMJ10ext", "tmj_metal_extrapolated.dat"),
                              ("TMJ10", "tmj.dat"),
                              ("TMJ10Padova", "tmj_padova.dat")])
            self.idx = np.intersect1d(self.idx,np.arange(25) )
            self.table = os.path.join(tables_dir, tablename[self.modelname])
            self.pars = np.loadtxt(self.table, usecols=range(3))
            self.data = np.loadtxt(self.table, usecols=self.idx+3)
        else:
            raise NotImplementedError("Model {0} is not available.".format(
                self.modelname))
        self.model = LinearNDInterpolator(self.pars, self.data)
        return

    def calc_lims(self):
        """ Get limits of model. """
        self.lims = np.array([self.pars.min(axis=0), self.pars.max(axis=0)]).T
        self.ranges = np.array([self.data.min(axis=0),
                                self.data.max(axis=0)]).T
        return

    def __call__(self, *args): 
        return self.model(*args)

def run_mcmc(lick, error, modelname, idx, dbname):
    """ Run the MCMC routine. """
    model = SSP(modelname, idx)
    ##########################################################################
    # Setting the priors
    age_dist = pymc.Uniform(name="age_dist", lower=model.lims[0,0],
                            upper=model.ranges[0,1])
    metal_dist = pymc.Uniform(name="metal_dist", lower=model.lims[1,0],
                              upper=model.lims[1,1])
    alpha_dist = pymc.Uniform(name="alpha_dist", lower=model.lims[2,0],
                              upper=model.lims[2,1])
    ##########################################################################
    taus = 1 / error[model.idx]**2
    @pymc.deterministic()
    def ssp(age=age_dist, metal=metal_dist, alpha=alpha_dist):
        return model(age, metal, alpha)
    y = pymc.Normal(name="y", mu=ssp, tau=taus,
                    value=lick[model.idx], observed=True)
    M = pymc.Model([y, age_dist, metal_dist, alpha_dist])
    mcmc = pymc.MCMC(M, db="txt", dbname=dbname)
    mcmc.sample(20000, 1000, 4)
    mcmc.db.close()
    return


def run_fields(fields, targetSN, nsim=30, modelname="TMJ10ext"):
    """ Interface to run MCMC in a list of fields including all bins"""
    for field in fields:
        os.chdir(os.path.join(data_dir, "combined_{0}".format(field),
                              "logs_sn{0}".format(targetSN)))
        lickfiles = sorted([x for x in os.listdir(".") if x.startswith("lick")
                            and x.endswith("corr.txt")])
        for flick in lickfiles:
            errorfile = "mc" + flick.replace("corr.txt", "nsim{0}.txt".format(
                        nsim))
            if not os.path.exists(errorfile):
                continue
            dbname = "mcmc_{0}_{1}_{2}".format(field, flick.split("_")[2],
                                               modelname)
            if os.path.exists(dbname):
                continue
            lick = np.loadtxt(flick, usecols=np.arange(1,26))
            error = np.nanstd(np.loadtxt(errorfile), axis=0)
            idx = np.arange(25)[np.logical_and(~np.isnan(error),
                                               ~np.isnan(lick))]
            # Remove Mg1 and Mg2
            idx = idx[~np.in1d(idx, [14,15])]
            run_mcmc(lick, error, modelname, idx, dbname)

if __name__ == "__main__":
    targetSN = 100
    run_fields(fields[:1], targetSN, nsim=30)