#! /usr/bin/env python
#-*- coding: utf-8 -*-
#==============================================================================
#
#          FILE: plot_functions.py
#
#   DESCRIPTION: This file contains a few plotting functions to visualize the 
#		 output
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tejaswi N
#       CREATED: 11-03-2014
# LAST MODIFIED: Wed Mar  4 14:52:59 2015
#      REVISION: ---
#==============================================================================

import matplotlib.pyplot as plt
import numpy as np
import sys

def main():
    """docstring for main"""

    try:
	    with open('hiresgstarKTstat.dat', 'r') as f:
		    resultsnew = np.loadtxt(f)
	    with open('hiresgstarKT.dat', 'r') as g:
		    resultsold = np.loadtxt(g)
	    with open('SMgstar.dat', 'r') as h:
		    resultspaper = np.loadtxt(h, skiprows=2)
    except IOError:
	    sys.exit('Result file is unreadable')
    except IndexError:
	    sys.exit('No result file passed')

    # to plot the gstar
    fig = plt.figure()
    
    ax = fig.add_subplot(111)
    
    ax.set_xscale('log')
    #ax.set_yscale('log')

    ax.scatter(resultsold[:,0], resultsold[:,1], s=5, c=u'k', edgecolors='none', label='$g_* \ {\\rm KT}$')
    ax.scatter(resultsold[:,0], resultsold[:,2], s=5, c=u'k', marker=u'*', edgecolors='none', label='$g_*,{\\rm s} \ {\\rm KT}$')
    ax.scatter(resultsnew[:,0], resultsnew[:,1], s=5, c=u'b', edgecolors='none', label='$g_* \ \\rm {KT, stat}$')
    ax.scatter(resultsnew[:,0], resultsnew[:,2], s=5, c=u'b', marker=u'*', edgecolors='none', label='$g_*,{\\rm s} \ \\rm{KT, stat}$')
    ax.scatter(resultspaper[:,0], resultspaper[:,1], s=5, c=u'r', edgecolors='none', label='$g_* \ \\rm{LS (2006)}$')
    ax.scatter(resultspaper[:,0], resultspaper[:,2], s=5, c=u'r', marker=u'*', edgecolors='none', label='$g_*,{\\rm s} \ \\rm{LS (2006)}$')

    ax.legend(loc='upper left', frameon=False, scatterpoints=1)

    ax.set_xlabel('$T \ ({\\rm MeV})$')
    ax.set_ylabel('$g_*$')
    
    plt.savefig('pltgstar.pdf')


if __name__ == "__main__":
     main()

pass
