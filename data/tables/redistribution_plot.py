#! /usr/bin/env python
#-*- coding: utf-8 -*-
#==============================================================================
#
#          FILE: redistribution_plot.py
#
#   DESCRIPTION: This file plots the asymmetry potential
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tejaswi N
#       CREATED: 19-02-2015
# LAST MODIFIED: Mon Jul  6 17:46:51 2015
#      REVISION: ---
#==============================================================================

import sys
import numpy as np
from scipy import integrate as intg
from scipy import interpolate as ip
import matplotlib.pyplot as plt
from matplotlib import rcParams # Equate size to that of companion plot

rcParams.update({'figure.autolayout': True})
plt.rc('text', usetex=True)
rcParams['text.latex.preamble'] = [r'\boldmath']

me = 0.511          # Mass of electron in MeV
mmu = 105.7         # Mass of muon in MeV
mtau = 1777         # Mass of tau lepton in MeV
mnue = 0.0          # Mass of electron neutrino in MeV
mnumu = 0.0         # Mass of mu neutrino in MeV
mnutau = 0.0        # Mass of tau neutrino in MeV
gf = 1.16637e-11    # Fermi constant in MeV^{-2}
s2thetaw = 0.23     # sin^2(weak mixing angle)
mz = 91.1876e+3     # Mass of the Z boson in MeV
mw = 80.385e+3      # Mass of the W boson in MeV

def main():
    """docstring for main"""

    redistarray = np.loadtxt('dmudLmu_new.dat', skiprows=1)
    redistarray_sb = np.loadtxt('dmudLmu_new_sb.dat', skiprows=1)

    chiarray = np.loadtxt('ChiTable_alltemp_new.dat', skiprows=1)
    chiarray_sb = np.loadtxt('ChiTable_Stefan_Boltz.dat', skiprows=1)

    # \partial \hat{V}^L/\partial \hat{\mathcal{L}}_\mu
    potentials_l_mu = np.array( [Vl(
	    redistarray[i,1:], 
	    chiarray[i,0], 
	    chiarray[i,1], 
	    chiarray[i,2], 
	    chiarray[i,3])[1] for i in range(len(chiarray))] )

    potentials_l_mu_sb = np.array( [Vl(
	    redistarray_sb[i,1:], 
	    chiarray_sb[i,0], 
	    chiarray_sb[i,1], 
	    chiarray_sb[i,2], 
	    chiarray_sb[i,3])[1] for i in range(len(chiarray_sb))] )

    spl_redist = ip.splrep( np.log10(chiarray[:,0]), 
	    potentials_l_mu/(gf*chiarray[:,0]**2) )

    spl_redist_sb = ip.splrep( np.log10(chiarray_sb[:,0]), 
	    potentials_l_mu_sb/(gf*chiarray_sb[:,0]**2) )

    # set of temperatures for plotting
    plot_temps = np.logspace(np.log10(chiarray[0,0]), 
	    np.log10(chiarray[-1,0]), num=200)

    no_redist = 2.0*np.sqrt(2.0)*np.ones(len(plot_temps))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xscale('Log')

    ax.set_xlabel(r"$T \ [\rm MeV]$", fontsize=18)
    ax.set_ylabel(r"$\partial V^{\rm L}_{\nu_\mu}/G_{\rm F} \partial \mathcal{L}_\mu$", fontsize=18)
    ax.plot( plot_temps, ip.splev(np.log10(plot_temps), spl_redist), c='k', 
		    lw=1.5, ls='-', label=r"${\rm Fit \ (this \ work)}$")
    ax.plot( plot_temps, ip.splev(np.log10(plot_temps), spl_redist_sb), c='k', 
		    lw=1.5, ls='--', label=r"${\rm Stefan \ Boltzmann}$")
    ax.plot( plot_temps, no_redist, c='k', 
		    lw=1.5, ls=':', label=r"${\rm No \ redist}$")

    ax.legend(frameon=False, loc=0)
    ax.tick_params(axis='both', which='major', labelsize=18)

    plt.savefig('potentials.pdf')
    #plt.show()
  
    #potentials = V(soln.x, pars.T, pars.p, chiq2, chib2, chibq11)

    #return potentials

    #mt = 0.1*np.arange(0,2001)

    #eft = np.array( [ef(0.0,m,1.0) for m in mt] )

    #np.savetxt('Efermi.dat', np.c_[ mt, np.log(eft) ], fmt='%15.4e', header='! m/T  log_e(E_f in units of T^4)')

    
def pchi(m):
    """This calculates the particle number susceptibility (in units of T^2)
    for a fermi dirac distribution with m (units of T), and degeneracy factor one"""

    result = intg.quad(lambda x: ((x**2)/(np.pi**2))*(np.exp(-np.sqrt(x**2 + m**2))/((np.exp(-np.sqrt(x**2 + m**2)) + 1.0)**2)), 0.0, np.inf)

    return result[0]

def ef(mu,m,g):
    """This calculates the energy density (in units of T^4) in a fermi-dirac 
    distribution with rest-mass m (in units of T), chemical potential mu (in 
    units of T) and degeneracy factor g"""

    result = intg.quad(lambda x: g*((x**2)/(2.0*np.pi**2))*np.sqrt(x**2 + m**2)*np.exp(-np.sqrt(x**2 + m**2))/(np.exp(-mu) + np.exp(-np.sqrt(x**2 + m**2))), 0, np.inf)

    return result[0]

def Vl(mu_array, T, chib2, chiq2, chibq11):
    """This function calculates the asymmetry potentials/T for neutrinos
    and anti-neutrinos of all generations given the chemical potentials 
    (in units of T), and temperature"""
    
    # [n_e- - n_e+, n_mu- - n_mu+, n_tau- - n_tau+, n_nu_e - n_nu_e_bar, 
    # n_nu_mu - n_nu_mu_bar, n_nu_tau - n_nu_tau_bar, n_q, n_b] in units
    # of T^3
    delta_n = np.zeros(8) 

    delta_n[0] = 2.0*pchi(me/T)*mu_array[0]
    delta_n[1] = 2.0*pchi(mmu/T)*mu_array[1]
    delta_n[2] = 2.0*pchi(mtau/T)*mu_array[2]
    
    delta_n[3] = pchi(mnue/T)*mu_array[3]
    delta_n[4] = pchi(mnumu/T)*mu_array[4]
    delta_n[5] = pchi(mnutau/T)*mu_array[5]

    delta_n[6] = chiq2*mu_array[6] + chibq11*mu_array[7]
    delta_n[7] = chibq11*mu_array[6] + chib2*mu_array[7]

    b0_array = np.zeros(3) # [b0_nu_e, b0_nu_mu, b0_nu_tau], Eq. (4.29) of the notes
 
    b0_array[0] = (np.sqrt(2)*gf*T**2*(
            (0.5 + 2.0*s2thetaw)*delta_n[0] + 
            (-0.5 + 2.0*s2thetaw)*delta_n[1] + 
            (-0.5 + 2.0*s2thetaw)*delta_n[2] + 
            2.0*delta_n[3] + 
            delta_n[4] + 
            delta_n[5] - 
            0.5*delta_n[7] + 
            (1.0 - 2.0*s2thetaw)*delta_n[6]))
    b0_array[1] = (np.sqrt(2)*gf*T**2*(
            (-0.5 + 2.0*s2thetaw)*delta_n[0] + 
            (0.5 + 2.0*s2thetaw)*delta_n[1] + 
            (-0.5 + 2.0*s2thetaw)*delta_n[2] + 
            delta_n[3] + 
            2.0*delta_n[4] + 
            delta_n[5] - 
            0.5*delta_n[7] + 
            (1.0 - 2.0*s2thetaw)*delta_n[6]))
    b0_array[2] = (np.sqrt(2)*gf*T**2*(
            (-0.5 + 2.0*s2thetaw)*delta_n[0] + 
            (-0.5 + 2.0*s2thetaw)*delta_n[1] + 
            (0.5 + 2.0*s2thetaw)*delta_n[2] + 
            delta_n[3] + 
            delta_n[4] + 
            2.0*delta_n[5] - 
            0.5*delta_n[7] + 
            (1.0 - 2.0*s2thetaw)*delta_n[6]))

    return b0_array

def Vth(T):
    """This function calculates the thermal potential per unit energy 
    for neutrinos and anti-neutrinos of all generations"""
    
    # [rho_e, rho_mu, rho_tau, rho_nu_e, rho_nu_mu, rho_nu_tau]
    rho_array = np.zeros(6)

    rho_array[0] = 2.0*ef(0,me/T,2)
    rho_array[1] = 2.0*ef(0,mmu/T,2)
    rho_array[2] = 2.0*ef(0,mtau/T,2)
    rho_array[3] = 2.0*ef(0,mnue/T,1)
    rho_array[4] = 2.0*ef(0,mnumu/T,1)
    rho_array[5] = 2.0*ef(0,mnutau/T,1)

    b1_array = np.zeros(3) # [b1_nu_e, b1_nu_mu, b1_nu_tau], Eq. (4.28) of the notes
 
    b1_array[0] = -((8.0*np.sqrt(2)/3.0)*gf*T**4*(
	    (1.0/mz)**2*rho_array[3] +
	    (1.0/mw)**2*rho_array[0]))
    b1_array[1] = -((8.0*np.sqrt(2)/3.0)*gf*T**4*(
	    (1.0/mz)**2*rho_array[4] +
	    (1.0/mw)**2*rho_array[1]))
    b1_array[2] = -((8.0*np.sqrt(2)/3.0)*gf*T**4*(
	    (1.0/mz)**2*rho_array[5] +
	    (1.0/mw)**2*rho_array[2]))

    return b1_array

if __name__ == "__main__":
    main()

pass

#def V(mu_array, T, p, chiq2, chib2, chibq11):
    #"""This function calculates the asymmetry potentials for the neutrinos and 
    #anti-neutrinos of all generations given the chemical potentials, 
    #temperature and | p | in MeV"""

    #gf = 1.16637e-11    # Fermi constant in MeV^{-2}
    #mz = 91.1876e+3     # Mass of the Z boson in MeV
    #mw = 80.385e+3      # Mass of the W boson in MeV
    #s2thetaw = 0.23     # sin^2(weak mixing angle)
    
    ## [n_e- - n_e+, n_mu- - n_mu+, n_tau- - n_tau+, n_nu_e - n_nu_e_bar, 
    ## n_nu_mu - n_nu_mu_bar, n_nu_tau - n_nu_tau_bar, n_q, n_b]
    #delta_n = np.zeros(8)  
    
    #delta_n[0] = nf(mu_array[0],me,T,2)[0] - nf(-mu_array[0],me,T,2)[0]
    #delta_n[1] = nf(mu_array[1],mmu,T,2)[0] - nf(-mu_array[1],mmu,T,2)[0]
    #delta_n[2] = nf(mu_array[2],mtau,T,2)[0] - nf(-mu_array[2],mtau,T,2)[0]

    #delta_n[3] = nf(mu_array[3],mnue,T,1)[0] - nf(-mu_array[3],mnue,T,1)[0]
    #delta_n[4] = nf(mu_array[4],mnumu,T,1)[0] - nf(-mu_array[4],mnumu,T,1)[0]
    #delta_n[5] = nf(mu_array[5],mnutau,T,1)[0] - nf(-mu_array[5],mnutau,T,1)[0]

    #delta_n[6] = chiq2*mu_array[6] + chibq11*mu_array[7]
    #delta_n[7] = chibq11*mu_array[6] + chib2*mu_array[7]

    ## [rho_e, rho_mu, rho_tau, rho_nu_e, rho_nu_mu, rho_nu_tau]
    #rho_array = np.zeros(6)

    #rho_array[0] = ef(mu_array[0],me,T,2)[0] + ef(-mu_array[0],me,T,2)[0]
    #rho_array[1] = ef(mu_array[1],mmu,T,2)[0] + ef(-mu_array[1],mmu,T,2)[0]
    #rho_array[2] = ef(mu_array[2],mtau,T,2)[0] + ef(-mu_array[2],mtau,T,2)[0]
    #rho_array[3] = ef(mu_array[3],mnue,T,1)[0] + ef(-mu_array[3],mnue,T,1)[0]
    #rho_array[4] = ef(mu_array[4],mnumu,T,1)[0] + ef(-mu_array[4],mnumu,T,1)[0]
    #rho_array[5] = ef(mu_array[5],mnutau,T,1)[0] + ef(-mu_array[5],mnutau,T,1)[0]

    #b0_array = np.zeros(3) # [b0_nu_e, b0_nu_mu, b0_nu_tau], Eq. (4.29) of the notes
    #b1_array = np.zeros(3) # [b1_nu_e, b1_nu_mu, b1_nu_tau], Eq. (4.28) of the notes
 
    #b0_array[0] = (np.sqrt(2)*gf*T**3*(
            #(0.5 + 2.0*s2thetaw)*delta_n[0] + 
            #(-0.5 + 2.0*s2thetaw)*delta_n[1] + 
            #(-0.5 + 2.0*s2thetaw)*delta_n[2] + 
            #2.0*delta_n[3] + 
            #delta_n[4] + 
            #delta_n[5] - 
            #0.5*delta_n[7] + 
            #(1.0 - 2.0*s2thetaw)*delta_n[6]))
    #b0_array[1] = (np.sqrt(2)*gf*T**3*(
            #(-0.5 + 2.0*s2thetaw)*delta_n[0] + 
            #(0.5 + 2.0*s2thetaw)*delta_n[1] + 
            #(-0.5 + 2.0*s2thetaw)*delta_n[2] + 
            #delta_n[3] + 
            #2.0*delta_n[4] + 
            #delta_n[5] - 
            #0.5*delta_n[7] + 
            #(1.0 - 2.0*s2thetaw)*delta_n[6]))
    #b0_array[2] = (np.sqrt(2)*gf*T**3*(
            #(-0.5 + 2.0*s2thetaw)*delta_n[0] + 
            #(-0.5 + 2.0*s2thetaw)*delta_n[1] + 
            #(0.5 + 2.0*s2thetaw)*delta_n[2] + 
            #delta_n[3] + 
            #delta_n[4] + 
            #2.0*delta_n[5] - 
            #0.5*delta_n[7] + 
            #(1.0 - 2.0*s2thetaw)*delta_n[6]))

    #b1_array[0] = -((8.0*np.sqrt(2)/3.0)*gf*T**4*(
            #(1.0/mz)**2*rho_array[3] +
            #(1.0/mw)**2*rho_array[0]))
    #b1_array[1] = -((8.0*np.sqrt(2)/3.0)*gf*T**4*(
            #(1.0/mz)**2*rho_array[4] +
            #(1.0/mw)**2*rho_array[1]))
    #b1_array[2] = -((8.0*np.sqrt(2)/3.0)*gf*T**4*(
            #(1.0/mz)**2*rho_array[5] +
            #(1.0/mw)**2*rho_array[2]))

    ## [V_nu_e, V_nu_mu, V_nu_tau, V_nu_e_bar, V_nu_mu_bar, V_nu_tau_bar ]
    #V_array = np.zeros(6) 

    ##[b0_nu_e, b0_nu_mu, b0_nu_tau, b1_nu_e, b1_nu_mu, b1_nu_tau]
    ##V_array = np.zeros(6)
    ##V_array[0] = b0_array[0]
    ##V_array[1] = b0_array[1]
    ##V_array[2] = b0_array[2]
    ##V_array[3] = b1_array[0]
    ##V_array[4] = b1_array[1]
    ##V_array[5] = b1_array[2]

    #V_array[0] = b0_array[0] + b1_array[0]*p
    #V_array[1] = b0_array[1] + b1_array[1]*p
    #V_array[2] = b0_array[2] + b1_array[2]*p

    #V_array[3] = -b0_array[0] + b1_array[0]*p
    #V_array[4] = -b0_array[1] + b1_array[1]*p
    #V_array[5] = -b0_array[2] + b1_array[2]*p

    #return V_array

    #ax.plot( chiarray[:,0], 
	     #np.array( [ 2.0*pchi(me/x) for x in chiarray[:,0] ] )*muderivative[:,0],
	     #label=r'$e^-$' )
    #ax.plot( chiarray[:,0], 
	     #np.array( [ 2.0*pchi(mmu/x) for x in chiarray[:,0] ] )*muderivative[:,1],
	     #label=r'$\mu^-$' )
    #ax.plot( chiarray[:,0], 
	     #np.array( [ 2.0*pchi(mtau/x) for x in chiarray[:,0] ] )*muderivative[:,2],
	     #label=r'$\tau^-$' )
    #ax.plot( chiarray[:,0], 
	     #np.array( [ pchi(mnue/x) for x in chiarray[:,0] ] )*muderivative[:,3],
	     #label=r'$\nu_e$' )
    #ax.plot( chiarray[:,0], 
	     #np.array( [ pchi(mnumu/x) for x in chiarray[:,0] ] )*muderivative[:,4],
	     #label=r'$\nu_\mu$' )
    #ax.plot( chiarray[:,0], 
	     #np.array( [ pchi(mnutau/x) for x in chiarray[:,0] ] )*muderivative[:,5],
	     #label=r'$\nu_\tau$' )
    #ax.plot( chiarray[:,0], 
	     #chiarray[:,2]*muderivative[:,6] + chiarray[:,3]*muderivative[:,7],
	     #label=r'$Q$' )


