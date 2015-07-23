#! /usr/bin/env python
#-*- coding: utf-8 -*-
#==============================================================================
#
#          FILE: compare_PSDs.py
#
#   DESCRIPTION: This file plots output PSDs for m_s = 7.15 keV, for a number
#                of runs with omega_s h^2 \simeq 0.12. It takes in the
#                data from parameter and data files in the appropriate sub-dire
#                ctories. Run it from the root directory w/ all sub-dirs.
#                It is set up to plot e^2 f(e) as in arxiv:1403.0954.
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tejaswi N
#       CREATED: 11-03-2014
# LAST MODIFIED: Wed Mar 25 16:01:30 2015
#      REVISION: ---
#==============================================================================

import numpy as np
import sys, os
import re, fnmatch
import matplotlib.pyplot as plt
import matplotlib.cm as cm # colormap
import matplotlib.colorbar as cb # colobar
import matplotlib.colors as cs # colors to normalize

def main():
    """docstring for main"""

    paramfile = 'params.dat'
    statefile = 'state.dat'

    results = {}
    params = {}
    states = {}
    temps = {}

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
        # read in names of data files, and save indices to sort
        datafiles = fnmatch.filter(os.listdir(dir), 'Snapshot*.dat')
        index = {}
        for filename in datafiles:
            # read off index
            index_reg = re.search(r"Snapshot([0-9]*?)\.dat", filename)
            index[filename] = int(index_reg.group(1))
        # Now sort data files based on index
        datafiles.sort(key=index.get)
        # now read in file with least temp == greatest index
        try:
            with open(os.path.join('.',dir,datafiles[-1]), 'r') as f:
                results[dir] = np.loadtxt(f, skiprows=1)
        except IOError:
            sys.exit('Data file is unreadable')
        except IndexError:
            sys.exit('No data file passed')
        # save the temperature for later use
        temps[dir] = states[dir][-1,0]
        
    if (len(os.walk('.').next()[1])==0):
        sys.exit('Give me something to work with!')

    # Now create a sorted set of tuples for params based on s2
    sorted_params = params.items()
    sorted_params.sort(key=lambda tup: tup[1][1])

    #ms = params[0]
    #s2 = params[1]
    #l = params[2]
    #omegas = params[3]
    
    # set colors
    colors = cm.rainbow(np.linspace(0, 1, len(sorted_params)))

    fig = plt.figure()
    ax = fig.add_subplot('111')
    fig.suptitle('PSDs for $m_s = {0:.2g} \ \\rm MeV$, $(\Omega_s h^2)_f = {1:.2g}$'.
            format(sorted_params[0][1][0], sorted_params[0][1][3])) # Plot title

    for (root, pars), c in zip(sorted_params, colors):
        ax.plot((1.0/temps[root])*results[root][:,0],
                100.0*((1.0/temps[root])*results[root][:,0])**2*results[root][:,1], 
                lw=1.5, color=c, label='$\sin^2 2\\theta = {0:.2g}$,' \
                '$(L/n_\gamma)_i = {1:.2g}$'.format(pars[1], pars[2]))

    ax.legend(frameon=False, numpoints=1)
    ##  plot upto p/T=8
    ax.set_xlim([0.0,8.0])

    ax.set_xlabel('$\epsilon = p/T$')
    ax.set_ylabel('$\epsilon^2 f(\epsilon) (\\times 100)$')

    plt.savefig('e2fe_all.pdf')

if __name__ == "__main__":
     main()

pass
