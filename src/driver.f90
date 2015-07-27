! Driver routine for sterile neutrino production
! Usage: ./<executable name> <parameter file name>
! This file also contains subroutines to read parameters annd free
! memory 
!---------------------------------------------------------------------
      program driver!{{{
      !---------------------------------------------------------------
      use Defaults
      use Sterile
      use ODE
      use DArrays
      implicit none

      character(LEN=Ini_max_name_len) :: parfile
      integer :: i, j

      if (iargc() < 1) then 
        stop "usage: ./<executable name> <parameter file name>"
      else
        call getarg(1, parfile)
      end if

      call Initialize(parfile)          ! Read in parameters and data,
                                        ! and setup functions

      do i=1, Sterile_n_models
        call ODE_Init_From_File(i)      ! Initialize model
        if (Sterile_leptasym_mode) then
          call ODE_Driver_Simple        ! Use given lepton asymmetry
        else
          call ODE_Driver_Brent         ! Find leptasym for req. closure
        end if
        do j=1, ODE_save_count
          call ODE_output(j)            ! Save results to run directory
        end do
      end do

      call Cleanup                      ! Free memory
      
      stop
      !---------------------------------------------------------------
      end program ! driver!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine Initialize(pfilename)!{{{
      !---------------------------------------------------------------
      ! Loads external files and initializes some modules that don't 
      ! change during execution
      use IniFile
      use Data_Repo
      use CPars
      use DArrays
      use Functions
      use Sterile
      implicit none

      character (LEN=*), intent(IN) :: pfilename

      call Ini_Open(pfilename, 80) ! Initialize parameter repository and
                                   ! read in parameter file
      call Data_Repo_Init          ! Initialize Data repository

      call CPars_Read              ! Read in cosmoparams
      call DArrays_Read_Files      ! Read in data into repository
      call DArrays_Init            ! Initialize data arrays to be used
                                   ! by functions
      call Functions_Splines_Init  ! Initialize Splines

      call Sterile_Params_Read     ! Read in sterile neutrino parameters
      return
      !---------------------------------------------------------------
      end subroutine ! Initialize!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine Cleanup!{{{
      !---------------------------------------------------------------
      ! Frees up memory
      use IniFile
      use Data_Repo
      use DArrays
      use Functions
      use Sterile

      call Ini_Close
      call Data_Repo_Close
      call DArrays_Close
      call Functions_Splines_Close
      call Sterile_Params_Close
      !---------------------------------------------------------------
      end subroutine ! Cleanup!}}}
!---------------------------------------------------------------------
