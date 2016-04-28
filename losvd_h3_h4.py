# -*- coding: utf-8 -*-
"""

Created on 04/04/16

@author: Carlos Eduardo Barbosa

Make plot of line of sight profiles using h3 and h4.

"""

import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import romb

from config import *

def gauss_hermite(v, mu, sigma, h3, h4):
    gauss = np.exp(-(v-mu)**2/(2 * sigma**2)) / \
                                  (sigma * np.sqrt(2 * np.pi))
    w = (v - mu)/sigma
    w2 = w**2
    poly = 1 + h3/np.sqrt(3)*(w*(2*w2-3)) \
             + h4/np.sqrt(24)*(w2*(4*w2-12)+3)
    return gauss * poly


if __name__ == "__main__":
    os.chdir(os.path.join(home, "plots"))
    v = np.linspace(-5,5,1001)
    h3s = np.linspace(0.2,-0.2,5)
    h4s = np.linspace(-0.2,0.2,5)
    moments = np.array(np.meshgrid(h3s,h4s)).T.reshape((-1,2))
    fig = plt.figure(figsize=(11,11))
    gs = gridspec.GridSpec(5,5)
    gs.update(left=0.08,right=0.99, top=0.99, bottom=0.06, wspace=0.08,
              hspace=0.08)
    for i, (h3,h4) in enumerate(moments):
        ax = plt.subplot(gs[i])
        ax.minorticks_on()
        ax.plot(v, gauss_hermite(v,0,1,h3,h4), "-r", lw=2.)
        # ax.plot(v, gauss_hermite(v,0,1,h3,h4) / np.trapz(gauss_hermite(v,0,1,h3,h4),v), "-b")
        # print h3, h4, np.trapz(gauss_hermite(v,0,1,h3,h4), v)
        ax.set_xlim(-5,5)
        ax.set_ylim(-0.1, 0.6)
        ax.axhline(y=0, ls="--", c="k")
        ax.text(-3.5, 0.47, "$h_3={0}$, $h_4={1}$".format(h3,h4), size=12)
        if h3==-0.2:
            ax.set_xlabel("$y$")
        else:
            ax.xaxis.set_ticklabels([])
        if h4==-.2:
            ax.set_ylabel("$L(v)$")
        else:
            ax.yaxis.set_ticklabels([])
        plt.locator_params(nbins=7)
    plt.savefig("h3h4.png")
