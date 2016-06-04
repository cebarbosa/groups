"""

V1.0.0: Michele Cappellari, Oxford, 04 April 2014
V1.0.1: Included angle keyword and kwargs.
    MC, Campos do Jordao, Brazil, 23 November 2015
V1.0.2: Updated documentation and included usage warning, MC, Oxford, 11 January 2015

"""

import warnings
import numpy as np

from cap_display_pixels import display_pixels

def display_bins_generators(xBin, yBin, velBin, x, y, angle=None, **kwargs):
    """
    Displays a Voronoi binned map starting from the original coordinates of the pixels
    and the coordinates of the *generators* (not the centroids!) of the Voronoi
    tessellation, as provided in output e.g. by my voronoi_2d_binning routine.

    NB: When possible, insted of this routine, one should use the more general display_bins
    routine which uses the binNumber of every spaxel insted of the Voronoi generators.

    :param xBin: coordinates of the *generators* of the Voronoi tessellation
    :param yBin:
    :param velBin:
    :param x: coordinates of the original spaxels
    :param y:
    :param angle:
    :param kwargs:
    :return: image
    """

    warnings.warn('When possible, usage of the routine display_bins is preferred to display_bins_generators')

    if not (xBin.size == yBin.size == velBin.size):
        raise ValueError('The vectors (XBIN, YBIN, VEL) must have the same size')
    if x.size != y.size:
        raise ValueError('The vectors (X, Y) must have the same size')
    if x.size < xBin.size:
        raise ValueError('The vectors (X, Y) cannot be smaller than (XBIN, YBIN)')
    
    # Perform a Voronoi tessellation starting from the coordinates
    # of the generators and the coordinates of the original pixels
    #
    binNum = np.argmin((x[:, np.newaxis] - xBin)**2 + (y[:, np.newaxis] - yBin)**2, axis=1)
    f = display_pixels(x, y, velBin[binNum], angle=angle, **kwargs)
    
    return f
    

