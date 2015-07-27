! Default parameters controlling the code's operation that can be 
! changed by the user if needed. This file contains the following
! modules:
! 1. Defaults
!---------------------------------------------------------------------
      module Defaults!{{{
     !---------------------------------------------------------------
      implicit none
      public
     
      !-------------------------------------------------------------
      ! Parameters for the param file parser!{{{
      ! Maximum size of name field in param file
      integer, parameter :: INI_MAX_NAME_LEN = 128
      ! Maximum size of row in param file
      integer, parameter :: INI_MAX_STRING_LEN = 1024
      ! What to do when parameter is not passed in. <t> todo: add 
      ! default params, then it can be set to false </t>
      logical :: INI_FAIL_ON_NOT_FOUND = .true.!}}}

      !-------------------------------------------------------------
      ! Parameters for the data file parser!{{{
      ! Maximum size of row in data file
      integer, parameter :: TAB_MAX_ROW_LEN = 2048
      ! Make IO verbose 
      logical :: IO_VERBOSE = .true.!}}}

      !-------------------------------------------------------------
      ! Parameters for functions and psd evolution !{{{
      ! High and low temperatures in MeV. Ensure that low is above
      ! weak decoupling at \simeq 2 MeV, and that the temperatures 
      ! in the data files cover the range
      double precision, parameter :: T_HIGH = 1.0D+4
      double precision, parameter :: T_LOW = 1.0D+1
      ! Number of momentum bins
      integer, parameter :: N_P_BINS = 1000
      ! Lowest & highest momentum bins at start, in temperature units
      double precision, parameter :: PTMIN_IN = 1.0D-3
      double precision, parameter :: PTMAX_IN = 20.0D0 - 1.0D-5
      ! Acceptable relative error in integrator
      double precision, parameter :: D_TOL = 1.0e-6!}}}

      !-------------------------------------------------------------
      ! Parameters for ODE integrator's storage!{{{
      integer, parameter :: MAXSTP = 1000000000
      ! Maximum number of internal variables stored
      integer, parameter :: NMAX = 4096
      ! Maximum number of intermediate steps stored
      integer, parameter :: KMAXX = 1024
      ! Flag to output intermediate temperatures in production
      logical, parameter :: INTFLAG = .true.
      ! Number of intermediate temperatures needed, if applicable (>1)
      integer, parameter :: NSAVE_MAX = 100!}}}

      !-------------------------------------------------------------
      ! Parameters for brent solver for lepton asymmetry!{{{
      ! limits for L_i/T^3 (hardcoded for now, change later)
      double precision, parameter :: L_LOW = 1.000D-3 
      double precision, parameter :: L_HIGH = 6.000D-3
      ! tolerance
      double precision, parameter :: L_TOL = 5.0D-6!}}}

      !-------------------------------------------------------------
      Save
      end module ! Defaults!}}}
!---------------------------------------------------------------------

