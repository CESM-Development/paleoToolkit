from mpl_toolkits.basemap import cm as bcm
import numpy as np
import matplotlib as mpl
import math, sys


def topography_cmap(n, start=0., end=1.0):
    """
    Constructs and returns a new GMT_globe colormap.
    ARGUMENTS:
        n   - number of levels in the colormap
        start - the normalized start fraction in the colormap (default 0.0, i.e start of original colormap)
        end - the normalized end fraction in the colormap (default 1.0, i.e end of original colormap)
    RETURNS:
        instance of mpl.colors.ListedColormap
    """
    # Creating a new colormap
    lvTmp  = np.linspace(start, end, n-1)
    cmTmp  = bcm.GMT_globe(lvTmp)
    newCmap= mpl.colors.ListedColormap(cmTmp)
    # Setting the colours of the out of range values
    newCmap.set_over(bcm.GMT_globe(end), alpha=1.0)
    newCmap.set_under(bcm.GMT_globe(start), alpha=1.0)
    return newCmap



def make_balanced(ll=None, ul=None):
    if (ll is not None) and (ul is not None):
        print("make_balanced: Please only specify lower limit, or upper limit")
        sys.exit()
    
    if ll is not None:
        ul = abs(round(ll/1.47, 2))
    else:
        ll = -1.*round(ul*1.481, 2)
    
    return ll, ul
