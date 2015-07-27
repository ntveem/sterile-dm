! This file reads in and saves the sterile neutrino parameters specified
! in the param file, for use in the rest of the code. This file contains
! the modules:
! 1. Sterile      
!---------------------------------------------------------------------
      module Sterile!{{{
      !---------------------------------------------------------------
      ! Reads in sterile neutrino parameters from paramfile and saves
      ! them for later use
      use IniFile
      use Defaults
      implicit None

      integer :: Sterile_n_models
      logical :: Sterile_leptasym_mode
      character (LEN=INI_MAX_NAME_LEN) :: Sterile_common_flavor
      double precision, pointer :: Sterile_ms_a(:), Sterile_s2_a(:),    &
     &                             Sterile_leptasymi_a(:)

      save
      contains
      !------------------------------------------------------------
        subroutine Sterile_Params_Read!{{{
        ! Reads in number and nature of sterile neutrino models from
        ! the parfile, and read in their parameters

          integer :: nms, ns2, nleptasymi, i

          ! Flavor of active neutrino mixing w/ steriles 
          ! (For now we only have mu opacities)
          Sterile_common_flavor = Ini_Read_String('flavor')

          ! Find number of models to run
          nms = Ini_NumberOf('ms') ! Number of Masses
          ns2 = Ini_NumberOf('s2') ! Number of Mixing angles
    
          if (nms /= ns2 .or. nms < 0) Then
            stop 'error in param file, each model should have ms and s2'
          else
            Sterile_n_models = nms
          end if
    
          allocate(Sterile_ms_a(Sterile_n_models))
          allocate(Sterile_s2_a(Sterile_n_models))
          allocate(Sterile_leptasymi_a(Sterile_n_models))

          do i=1, Sterile_n_models
            Sterile_ms_a(i) = Ini_Read_Double('ms',i)
            Sterile_s2_a(i) = Ini_Read_Double('s2',i)
          end do

          ! Find how we are fixing the initial asymmetry
          nleptasymi = Ini_NumberOf('leptasymi') ! # leptasyms
          if (nleptasymi < 1) Then
            ! We're fixing the lepton asymmetry by the closure relation
              Sterile_leptasym_mode = .false.
          else if (nleptasymi == Sterile_n_models) Then
            Sterile_leptasym_mode = .true.
            do i=1, Sterile_n_models
              Sterile_leptasymi_a(i) = Ini_Read_Double('leptasymi',i)
            end do
          else
                  stop "error in param file, if using given Ls, their   &
     & number should equal that of ms and s2"
          end if
        end subroutine ! Sterile_Params_Read!}}}
      !------------------------------------------------------------

      !------------------------------------------------------------
        subroutine Sterile_Params_Close!{{{
        ! Free up memory

          deallocate(Sterile_ms_a)
          deallocate(Sterile_s2_a)
          deallocate(Sterile_leptasymi_a)

        end subroutine ! Sterile_Params_Close!}}}
      !------------------------------------------------------------
      !---------------------------------------------------------------
      end module ! Sterile!}}}
!---------------------------------------------------------------------

