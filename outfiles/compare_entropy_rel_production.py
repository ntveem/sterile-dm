#! /usr/bin/env python
#-*- coding: utf-8 -*-
#==============================================================================
#
#          FILE: compare_entropy_rel_production.py
#
#   DESCRIPTION: This file plots n_s/entopy density vs temperature for 
#                m_s = 7.15 keV, for a number of runs with different 
#                \sin^2 2\theta. It takes in data from parameter and data files 
#                in the appropriate sub-directories. Run it from the root 
#                directory w/ all sub-dirs.
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tejaswi N
#       CREATED: 11-03-2014
# LAST MODIFIED: Wed May 27 10:47:32 2015
#      REVISION: ---
#==============================================================================

import numpy as np
import sys, os
import re, fnmatch
import matplotlib.pyplot as plt
import matplotlib.cm as cm # colormap
import matplotlib.colorbar as cb # colobar
import matplotlib.colors as cs # colors to normalize

from scipy import interpolate as ip

def main():
    """docstring for main"""

    paramfile = 'params.dat'
    statefile = 'state.dat'

    results = {}
    params = {}
    states = {}
    #num_dens = {}

    numf2 = 0.182684 # (3\zeta(3)/(2\pi^2))
    nudens = 1.06556e4 # nu density*\rho_crit/h^2
    gstarfac = 2.0*np.pi**2/45.0

    # read T vs g_{*,s}
    try:
        with open('../data/tables/SMgstar.dat', 'r') as g:
            gstart = np.loadtxt(g)
    except IOError:
        sys.exit('File not readable')
    except IndexError:
        sys.exit('File not passed')

    # define splines for g_*,s vs log(T)
    gstars = ip.splrep(np.log(gstart[:,0]), gstart[:,2])

    for dir in os.walk('.').next()[1]:
        # each entry is a subdirectory, possibly with data
        # read in parameter and state files
        try:
            with open(os.path.join('.',dir,paramfile), 'r') as g:
                params[dir] = np.loadtxt(g)
            with open(os.path.join('.',dir,statefile), 'r') as f:
                states[dir] = np.loadtxt(f, skiprows=1)
        except IOError:
            # This isn't an output directory, skip.
            continue
        except IndexError:
            # This isn't an output directory, skip.
            continue
        
    if (len(os.walk('.').next()[1])==0):
        sys.exit('Give me something to work with!')

    # Now create a sorted set of tuples for params based on s2
    sorted_params = params.items()
    sorted_params.sort(key=lambda tup: tup[1][1])

    #ms = params[0]
    #s2 = params[1]
    #l = params[2]
    #omegas = params[3]

    # to calculate colors, define function that equals 0 at min,
    # 1 at max, and logarithmically varies with s2\theta in between
    s2min = sorted_params[0][1][1]
    s2max = sorted_params[-1][1][1]
    s2array = np.array(
	[sorted_params[i][1][1] for i in range(len(sorted_params))])
    # set colors
    colors = cm.rainbow(np.log(s2array/s2min)/np.log(s2max/s2min))

    fig = plt.figure()
    fig.suptitle('$n_s/s$ vs temperature for $m_s = {0:.2g} \ \\rm MeV$'.
            format(sorted_params[0][1][0])) # Plot title

    ax = fig.add_axes([0.125, 0.1, 0.725, 0.8])
    ax.set_xscale('log')
    #ax.set_yscale('log')

    for (root, pars), c in zip(sorted_params, colors):
        ax.plot(states[root][:,0], states[root][:,4]*numf2/(params[root][0]*nudens*gstarfac*ip.splev(np.log(states[root][:,0]),gstars)), 
                lw=1.5, color=c)

    # draw colorbar
    ax2 = fig.add_axes([0.85, 0.1, 0.025, 0.8])
    norm = cs.LogNorm(vmin=s2min, vmax=s2max)
    cb1 = cb.ColorbarBase(ax2, cmap=cm.rainbow, norm=norm, 
            ticks=s2array, orientation='vertical')
    cb1.ax.set_yticklabels(["{0:.2g}".format(sorted_params[i][1][1]) for i
                           in range(len(sorted_params))])
    cb1.ax.set_title('$\sin^2 2\\theta$')


    ax.set_xlabel('$T (\\rm MeV)$')
    ax.set_ylabel('$n_s/s$')

    plt.savefig('n_s_all_v2.pdf')

if __name__ == "__main__":
     main()

pass
