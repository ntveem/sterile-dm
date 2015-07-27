! Some constants needed by the code. 
! This file contains the following modules:
! 1. Phconst -- Defines some mathematical and physical constants.
! 2. Cosmoparams -- Reads and stores relevant cosmological parameters.      
!---------------------------------------------------------------------
      module Phconst!{{{
      !---------------------------------------------------------------
      implicit none
      public

      ! Mathematical constants
      double precision :: pi, zeta3

      ! Physical constants
      double precision :: mp, Gf, hbar, kb
      double precision :: me, mmu, mtau
      double precision :: mnue, mnumu, mnutau
      double precision :: mz, mw, s2w

      parameter (pi=3.1415926535897931D0,zeta3=1.202056903159594D0)
      ! m_p = 1/sqrt{G} in MeV, G_f in MeV^-2
      parameter (mp=1.22089201D+22,Gf=1.1663787D-11)
      ! hbar in MeV-s, kb in MeV K^-1
      parameter (hbar=6.58211928D-22,kb=8.6173325D-11)
      ! Particle masses in MeV
      parameter (me=0.511D0,mmu=1.05658D+2,mtau=1.77682D+3)
      parameter (mnue=0.0D0,mnumu=0.0D0,mnutau=0.0D0)
      parameter (mz=91.1876D+3,mw=80.385D+3)
      ! Weak mxing angle, \sin^2(2\theta_W)
      parameter (s2w=0.23126D0)

      save
      !---------------------------------------------------------------
      end module ! Phconst!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      module CPars!{{{
      !---------------------------------------------------------------
      use IniFile
      use Phconst
      implicit None
      public

      double precision :: CPARS_omegadmh2
      double precision :: CPARS_TcmbMeV
      ! factor relating nu and photon temperatures after e^-e^+ epoch
      double precision, parameter :: CPARS_nfactor = (4.0D0/11.0D0)**   &
     &                                               (1.0D0/3.0D0)
      ! (critical density/h^2) today in MeV^-4
      double precision, parameter :: CPARS_rhocrit = 8.095921884D-35

      save
      contains
      !-------------------------------------------------------------
        subroutine CPars_Read!{{{
        ! Assigns cosmoparams from read param file

          double precision :: Tcmb

          CPARS_omegadmh2 = Ini_Read_Double('omegadmh2')
          Tcmb = Ini_Read_Double('Tcmb')
          CPARS_TcmbMeV = Tcmb*kb
        end subroutine ! CPars_Read!}}}
      !-------------------------------------------------------------
      !---------------------------------------------------------------
      end module ! CPars!}}}
!---------------------------------------------------------------------

