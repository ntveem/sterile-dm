#! /usr/bin/env python
#-*- coding: utf-8 -*-
#==============================================================================
#
#          FILE: redistribution.py
#
#   DESCRIPTION: This file calculates the effect of asymmetry redistribution,
#                and outputs a table of d\mu_X/dL_Y and dV_A/dL_Y
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tejaswi N
#       CREATED: 19-02-2015
# LAST MODIFIED: Mon 27 Jul 2015 04:46:15 AM PDT
#      REVISION: ---
#==============================================================================

import sys
import numpy as np
from scipy import integrate as intg
import matplotlib.pyplot as plt

me = 0.511                          # Mass of electron in MeV
mmu = 105.7                         # Mass of muon in MeV
mtau = 1777                         # Mass of tau lepton in MeV
mnue = 0.0                          # Mass of electron neutrino in MeV
mnumu = 0.0                         # Mass of mu neutrino in MeV
mnutau = 0.0                        # Mass of tau neutrino in MeV

def main():
    """docstring for main"""

    chiarray = np.loadtxt('ChiTable_alltemp_new.dat', skiprows=1)

    # nx8x8 array, where n is the number of temperatures tabulated
    # array relates chemical potentials to conserved quantities, 8 eqs
    coeffarray = np.dstack(
            (
                np.c_[
                    np.array( [ 2.0*pchi(me/x) for x in chiarray[:,0] ] ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    -np.ones( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    -np.array( [ 2.0*pchi(me/x) for x in chiarray[:,0] ] )
                    ],
                np.c_[
                    np.zeros( len(chiarray[:,0]) ),
                    np.array( [ 2.0*pchi(mmu/x) for x in chiarray[:,0] ] ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    -np.ones( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    -np.array( [ 2.0*pchi(mmu/x) for x in chiarray[:,0] ] )
                    ],
                np.c_[
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.array( [ 2.0*pchi(mtau/x) for x in chiarray[:,0] ] ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    -np.ones( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    -np.array( [ 2.0*pchi(mtau/x) for x in chiarray[:,0] ] )
                    ],
                np.c_[
                    np.array( [ pchi(mnue/x) for x in chiarray[:,0] ] ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.ones( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) )
                    ],
                np.c_[
                    np.zeros( len(chiarray[:,0]) ),
                    np.array( [ pchi(mnumu/x) for x in chiarray[:,0] ] ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.ones( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) )
                    ],
                np.c_[
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.array( [ pchi(mnutau/x) for x in chiarray[:,0] ] ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.ones( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) )
                    ],
                np.c_[
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    -np.ones( len(chiarray[:,0]) ),
                    -np.ones( len(chiarray[:,0]) ),
                    -np.ones( len(chiarray[:,0]) ),
                    chiarray[:,3],
                    chiarray[:,2]
                    ],
                np.c_[
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    np.zeros( len(chiarray[:,0]) ),
                    chiarray[:,1],
                    chiarray[:,3]
                    ]
                )
            )

    # result matrices for three flavors of lepton asymmetry
    Leresult = np.c_[
            np.ones( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) )
                    ]
    Lmuresult = np.c_[
            np.zeros( len(chiarray[:,0]) ),
            np.ones( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) )
                    ]
    Ltauresult = np.c_[
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.ones( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) ),
            np.zeros( len(chiarray[:,0]) )
                    ]

    Leresult_sb = np.c_[
            np.ones( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) )
                    ]
    Lmuresult_sb = np.c_[
            np.zeros( len(chiarray_sb[:,0]) ),
            np.ones( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) )
                    ]
    Ltauresult_sb = np.c_[
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.ones( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) ),
            np.zeros( len(chiarray_sb[:,0]) )
                    ]

    # derivatives w.r.t e, mu and tau asymmetries
    # each is an nx8 array, where each row is d\mu_\alpha/dL_x
    ederivative = np.linalg.solve(coeffarray, Leresult) 
    muderivative = np.linalg.solve(coeffarray, Lmuresult) 
    tauderivative = np.linalg.solve(coeffarray, Ltauresult)

    np.savetxt('dmudLe_new.dat', np.c_[ chiarray[:,0], ederivative  ], fmt='%15.4e', 
            header='!T (MeV) d\mu_e/dL_e d\mu_\mu/dL_e d\mu_\\tau/dL_e d\mu_{\\nu_e}/dL_e d\mu_{\\nu_\mu}/dL_e d\mu_{\\nu_\\tau}/dL_e d\mu_Q/dL_e d\mu_B/dL_e')
    np.savetxt('dmudLmu_new.dat', np.c_[ chiarray[:,0], muderivative  ], fmt='%15.4e',
            header='!T (MeV) d\mu_e/dL_\mu d\mu_\mu/dL_\mu d\mu_\\tau/dL_\mu d\mu_{\\nu_e}/dL_\mu d\mu_{\\nu_\mu}/dL_\mu d\mu_{\\nu_\\tau}/dL_\mu d\mu_Q/dL_\mu d\mu_B/dL_\mu')
    np.savetxt('dmudLtau_new.dat', np.c_[ chiarray[:,0], tauderivative  ], fmt='%15.4e',
            header='!T (MeV) d\mu_e/dL_\\tau d\mu_\mu/dL_\\tau d\mu_\\tau/dL_\\tau d\mu_{\\nu_e}/dL_\\tau d\mu_{\\nu_\mu}/dL_\\tau d\mu_{\\nu_\\tau}/dL_\\tau d\mu_Q/dL_\\tau d\mu_B/dL_\\tau')

       
def pchi(m):
    """This calculates the particle number susceptibility (in units of T^2)
    for a fermi dirac distribution with m (units of T), and degeneracy factor one"""

    result = intg.quad(lambda x: ((x**2)/(np.pi**2))*(np.exp(-np.sqrt(x**2 + m**2))/((np.exp(-np.sqrt(x**2 + m**2)) + 1.0)**2)), 0.0, np.inf)

    return result[0]

def V(mu_array, T, chib2, chiq2, chibq11):
    """This function calculates the asymmetry potentials/T for neutrinos
    and anti-neutrinos of all generations given the chemical potentials 
    (in units of T), and temperature"""

    gf = 1.16637e-11    # Fermi constant in MeV^{-2}
    s2thetaw = 0.23     # sin^2(weak mixing angle)
    
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

#def ef(mu,m,g):
    #"""This calculates the energy density (in units of T^4) in a fermi-dirac 
    #distribution with rest-mass m (in units of T), chemical potential mu (in 
    #units of T) and degeneracy factor g"""

    #result = intg.quad(lambda x: g*((x**2)/(2.0*np.pi**2))*np.sqrt(x**2 + m**2)*np.exp(-np.sqrt(x**2 + m**2))/(np.exp(-mu) + np.exp(-np.sqrt(x**2 + m**2))), 0, np.inf)

    #return result[0]
