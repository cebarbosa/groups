"""
    Copyright (C) 2016, Michele Cappellari
    E-mail: michele.cappellari_at_physics.ox.ac.uk

    Updated versions of the software are available from my web page
    http://purl.org/cappellari/software

    See example at the bottom for usage instructions.

    V1.0.0: Created to emulate my IDL procedure with the same name.
        Michele Cappellari, Oxford, 28 March 2014
    V1.0.1: Fixed treatment of optional parameters. MC, Oxford, 6 June 2014
    V1.0.2: Avoid potential runtime warning. MC, Oxford, 2 October 2014
    V1.0.3: Return axis. MC, Oxford, 26 March 2015
    V1.0.4: Return image instead of axis. MC, Oxford, 15 July 2015
    V1.0.5: Removes white gaps from rotated images using edgecolors.
        MC, Oxford, 5 October 2015
    V1.0.6: Pass kwargs to graphics functions.
        MC, Campos do Jordao, Brazil, 23 November 2015
    V1.0.7: Check that input (x,y) come from an axis-aligned image.
        MC, Oxford, 28 January 2016
    V1.0.8: Fixed deprecation warning in Numpy 1.11. MC, Oxford, 22 April 2016
    V1.1.0: Fixed program stop with kwargs. Included `colorbar` keyword.
        MC, Oxford, 18 May 2016

"""
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from scipy.spatial import distance
import numpy as np

from sauron_colormap import sauron

##############################################################################

def display_pixels(x, y, val, pixelsize=None, vmin=None, vmax=None,
                   angle=None, colorbar=False, label=None, nticks=7, **kwargs):
    """
    Display vectors of square pixels at coordinates (x,y) coloured with "val".
    An optional rotation around the origin can be applied to the whole image.
    
    The pixels are assumed to be taken from a regular cartesian grid with 
    constant spacing (like an axis-aligned image), but not all elements of
    the grid are required (missing data are OK).

    This routine is designed to be fast even with large images and to produce
    minimal file sizes when the output is saved in a vector format like PDF.

    """
    x, y, val = map(np.ravel, [x, y, val])

    if not (x.size == y.size == val.size):
        raise ValueError('The vectors (x, y, val) must have the same size')

    if vmin is None:
        vmin = np.min(val)

    if vmax is None:
        vmax = np.max(val)

    if pixelsize is None:
        pixelsize = np.min(distance.pdist(np.column_stack([x, y])))

    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    x1 = (x - xmin)/pixelsize
    y1 = (y - ymin)/pixelsize
    nx = int(round((xmax - xmin)/pixelsize) + 1)
    ny = int(round((ymax - ymin)/pixelsize) + 1)
    mask = np.ones((nx, ny), dtype=bool)
    img = np.empty((nx, ny))
    j = np.round(x1).astype(int)
    k = np.round(y1).astype(int)
    mask[j, k] = 0
    img[j, k] = val
    img = np.ma.masked_array(img, mask)

    if np.any(np.abs(np.append(j - x1, k - y1)) > 0.1):
        raise ValueError('The coordinates (x, y) must come from an axis-aligned image')

    ax = plt.gca()

    if (angle is None) or (angle == 0):

        img = ax.imshow(np.rot90(img), interpolation='none',
                        cmap=kwargs.get("cmap", sauron), vmin=vmin, vmax=vmax,
                        extent=[xmin-pixelsize/2, xmax+pixelsize/2,
                                ymin-pixelsize/2, ymax+pixelsize/2])

    else:

        x, y = np.ogrid[xmin-pixelsize/2 : xmax+pixelsize/2 : (nx+1)*1j,
                        ymin-pixelsize/2 : ymax+pixelsize/2 : (ny+1)*1j]
        ang = np.radians(angle)
        x, y = x*np.cos(ang) - y*np.sin(ang), x*np.sin(ang) + y*np.cos(ang)

        mask1 = np.ones_like(x, dtype=bool)
        mask1[:-1, :-1] *= mask  # Flag the four corners of the mesh
        mask1[:-1, 1:] *= mask
        mask1[1:, :-1] *= mask
        mask1[1:, 1:] *= mask
        x = np.ma.masked_array(x, mask1)  # Mask is used for proper plot range
        y = np.ma.masked_array(y, mask1)

        img = ax.pcolormesh(x, y, img, cmap=kwargs.get("cmap", sauron),
                            vmin=vmin, vmax=vmax, edgecolors="face")
        ax.axis('image')

    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        ticks = MaxNLocator(nticks).tick_values(vmin, vmax)
        cbar = plt.colorbar(img, cax=cax, ticks=ticks)
        cbar.solids.set_edgecolor("face")  # Remove gaps in PDF http://stackoverflow.com/a/15021541
        if label:
            cbar.set_label(label)
        plt.sca(ax)  # Activate main plot before returning


    ax.minorticks_on()
    ax.tick_params(length=10, width=1, which='major')
    ax.tick_params(length=5, width=1, which='minor')

    return img

##############################################################################

# Usage example for display_pixels()

if __name__ == '__main__':

    n = 41  # 1 arcsec pixels
    x = np.linspace(-20, 20, n)
    y = np.linspace(-20, 20, n)
    xx, yy = np.meshgrid(x,y)
    counts = xx**2 - 2*yy**2
    w = xx**2 + 2*yy**2 < 20.1**2

    plt.clf()
    ax = display_pixels(xx[w], yy[w], counts[w], pixelsize=x[1]-x[0],
                        angle=20, colorbar=1)
    plt.pause(0.01)
