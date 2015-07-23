#! /usr/bin/env python
#-*- coding: utf-8 -*-
#==============================================================================
#
#          FILE: plot_PSDs.py
#
#   DESCRIPTION: This file plots input PSDs over temperature, for a single run
#                It takes in the values of m_s, \sin^2 \theta and L, the 
#                temperatures and the data from parameter, state and data files 
#                in the same directory. Copy it into a run's output directory.
#                It is set up to plot e^2 f(e) as in arxiv:1403.0954.
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tejaswi N
#       CREATED: 11-03-2014
# LAST MODIFIED: Tue Mar 24 16:36:01 2015
#      REVISION: ---
#==============================================================================

import numpy as np
import sys
import glob, re
import matplotlib.pyplot as plt
import matplotlib.cm as cm # colormap
import matplotlib.colorbar as cb # colobar
import matplotlib.colors as cs # colors to normalize

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

     # all the names of the data files
    fnames = glob.glob('./Snapshot*.dat')

    if (len(fnames) != len(statearray)): # Sanity check
        sys.exit('There should be one output file for each line'\
                 ' in the state file')
   
    results = {}    # dictionary with all the results
    temp = {}       # dictionary with all the temperatures
    for name in fnames:
        # temperature index in state file
        temp_reg = re.search(r"Snapshot([0-9]*?)\.dat", name)
        idx = int(temp_reg.group(1))
        temp[name] = statearray[idx-1,0] # Fortran index starts at 1

        try:    # Load file with P.S.D
                with open(name, 'r') as f:
                        results[name] = np.loadtxt(f, skiprows=1)
        except IOError:
                print name
                sys.exit('Result file is unreadable')
        except IndexError:
                sys.exit('No result file passed')

    # Now sort filenames based on temperatures
    fnames.sort(key=temp.get)

    # set colors
    # note that log temperatures are equally spaced
    colors = cm.rainbow(np.linspace(0, 1, len(fnames)))

    fig = plt.figure()
    fig.suptitle('$m_s = {0:.2g} \ \\rm MeV$,'\
                 ' $\sin^2 2\\theta = {1:.2g}$,'\
                 ' $(L/n_\gamma)_i = {2:.2g}$,'\
                 ' $(\Omega_s h^2)_f = {3:.2g}$'.
                 format(ms, s2, l, omegas)) # Plot title

    ax = fig.add_axes([0.125, 0.1, 0.725, 0.8])
        
    #ax.set_xscale('log')
    #ax.set_yscale('log')

    for name, c in zip(fnames, colors):
        ax.plot((1.0/temp[name])*results[name][:,0],
                100.0*((1.0/temp[name])*results[name][:,0])**2*
                results[name][:,1], lw=1.5, color=c)

    #ax.legend(frameon=False, numpoints=1)
    #  plot upto p/T=8
    ax.set_xlim([0.0,8.0])

    # draw colorbar
    ax2 = fig.add_axes([0.85, 0.1, 0.025, 0.8])
    norm = cs.LogNorm(vmin=temp[fnames[0]], vmax=temp[fnames[-1]])
    cb1 = cb.ColorbarBase(ax2, cmap=cm.rainbow, norm=norm, 
            ticks=[temp[fnames[4*i]] for i in range(len(fnames)/4)],
            orientation='vertical')
    cb1.ax.set_yticklabels(["{0:.1f}".format(temp[fnames[4*i]]) for i
                           in range(len(fnames)/4)])
    cb1.ax.set_title('$T ({\\rm MeV})$')

    ax.set_xlabel('$\epsilon = p/T$')
    ax.set_ylabel('$\epsilon^2 f(\epsilon) (\\times 100)$')

    plt.savefig('e2fe.pdf')
    #plt.savefig('loge2fe.pdf')

if __name__ == "__main__":
     main()

pass
