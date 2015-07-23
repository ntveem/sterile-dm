#! /usr/bin/env python
#-*- coding: utf-8 -*-
#==============================================================================
#
#          FILE: plot_production.py
#
#   DESCRIPTION: This file plots L/n_gamma and \Omega_s h^2 vs temperature, for 
#		 a single run from a state file in the same directory. Copy it
#		 into a run's output directory.
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tejaswi N
#       CREATED: 19-03-2015
# LAST MODIFIED: Tue Mar 24 16:26:15 2015
#      REVISION: ---
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt

def main():
    """docstring for main"""

    paramfile = 'params.dat'
    statefile = 'state.dat'

    try:    # Load param file
            with open(paramfile, 'r') as f:
                    params = np.loadtxt(f)
    except IOError:
            sys.exit('Param file is unreadable')
    except IndexError:
            sys.exit('No param file present')

    ms = params[0]
    s2 = params[1]
    l = params[2]
    omegas = params[3]

    try:    # Load state file
        with open(statefile, 'r') as f:
            statearray = np.loadtxt(f, skiprows=1)
    except IOError:
        sys.exit('State file is unreadable')
    except IndexError:
        sys.exit('No state file present')

    colors = ['k', 'b']

    fig = plt.figure()
    fig.suptitle('$m_s = {0:.2g} \ \\rm MeV$,'\
                 ' $\sin^2 2\\theta = {1:.2g}$,'\
                 ' $(L/n_\gamma)_i = {2:.2g}$,'\
                 ' $(\Omega_s h^2)_f = {3:.2g}$'.
                 format(ms, s2, l, omegas)) # Plot title

    ax = fig.add_subplot('111')
    ax.set_xscale('log')

    ax.plot( statearray[:,0], statearray[:,3], '-o',
            color=colors[0], markersize=3.0, lw=1.2)
    ax.set_xlabel('$T (\\rm MeV)$')
    ax.set_ylabel('$L/n_\gamma$', color=colors[0])
    for t1 in ax.get_yticklabels():
        t1.set_color(colors[0])

    ax2 = ax.twinx()
    ax2.plot( statearray[:,0], statearray[:,4], '-^',
            color=colors[1], markersize=3.0, lw=1.2)
    ax2.set_ylabel('$\Omega_s h^2$', color=colors[1])
    for t1 in ax2.get_yticklabels():
        t1.set_color(colors[1])

    plt.savefig('LvsT.pdf')

if __name__ == "__main__":
     main()

pass
