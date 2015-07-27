! ODE integrator and associated variables and functions. This file
! contains functions to initialize storage space used by driver routines
! to integrate the sterile neutrino p.s.d from T_HIGH to T_LOW using the
! lsode solver, available at:
!          https://computation.llnl.gov/casc/odepack/odepack_home.html
! and included as src/lsode.f. There are two kinds of routines,
! depending on whether the initial lepton asymmetry is input or fixed to
! produce the right DM closure fraction. This file contans the following
! modules:
! 1. ODE
!---------------------------------------------------------------------
      module ODE!{{{
      !---------------------------------------------------------------
      use Defaults
      use Sterile
      use Functions
      use Phconst
      use CPars
      implicit none
      public

      ! odeintegrator's internal storage
      integer            :: ODE_save_count
      double precision   :: ODE_dx_save, ODE_T_save(KMAXX),             &
     &                      ODE_y_save(NMAX,KMAXX)
      double precision   :: ODE_p_bin_save(NMAX,KMAXX)

      ! sterile parameters
      double precision   :: ODE_ms, ODE_s2, ODE_leptasymi
      character(len=INI_MAX_NAME_LEN) :: ODE_flavor
      ! solver's internal variables
      ! number of variables
      integer, parameter :: ODE_nvars = 2*N_P_BINS + 3
      ! time, efolds, leptasym, psds of nu and nubar
      double precision   :: ODE_vars(ODE_nvars)
      ! momentum bins
      double precision   :: ODE_p_bins_init(N_P_BINS),                  &
     &                        ODE_p_bins(N_P_BINS)

      save
      contains
      !-------------------------------------------------------------
        subroutine ODE_Init_From_File(i)!{{{
        ! Init sterile neutrino params from those in the param file

          integer, intent(IN) :: i

          if (ODE_nvars.gt.NMAX) Then
            print *,ODE_nvars,NMAX
            stop "Not enough memory allocated to store internal         &
     &            variables of the integrator. Increase NMAX."
          end if


          if ((i > Sterile_n_models) .or. (i < 1)) then
            print *,"i=",i
            print *,"number of models =",Sterile_n_models
            print *,"Too many models requested."
            stop
          end if

          ODE_ms = Sterile_ms_a(i)
          ODE_s2 = Sterile_s2_a(i)
          ODE_leptasymi = Sterile_leptasymi_a(i)
          ODE_flavor = Sterile_common_flavor

        end subroutine ! ODE_Init_From_File!}}}
      !-------------------------------------------------------------

      !-------------------------------------------------------------
        subroutine ODE_Init_From_Params(ms_in, s2_in, l_in, fl_in)!{{{
        ! Init or update sterile neutrino params individually
          double precision, optional, intent(IN) :: ms_in, s2_in, l_in
          character (LEN=*), optional, intent(IN):: fl_in

          if (present(ms_in)) ODE_ms = ms_in
          if (present(s2_in)) ODE_s2 = s2_in
          if (present(l_in)) ODE_leptasymi = l_in
          if (present(fl_in)) ODE_flavor = fl_in

        end subroutine ! ODE_Init_From_Params!}}}
      !-------------------------------------------------------------

      !-------------------------------------------------------------
        subroutine ODE_Vars_Init!{{{
        ! Initialize ode integrator's internal variables from 
        ! parameters in default file and ode integrator's stored 
        ! sterile neutrino parameters

        double precision :: gstar_init, t_init 
        double precision :: pmin_in, pmax_in, dp_in
        integer :: i

        gstar_init = gstar(T_HIGH,1,0)
        t_init = sqrt(45.0d0/(16.0d0*pi**3))*mp*hbar/                   &
     &             (sqrt(gstar_init)*(T_HIGH)**2)

        pmin_in = PTMIN_IN*T_HIGH
        pmax_in = PTMAX_IN*T_HIGH

        ! initialize time, efolds, leptasym and p.s.ds
        ODE_vars(1) = t_init
        ODE_vars(2) = 0.0D0
        ODE_vars(3) = ODE_leptasymi
        Do i=4,ODE_nvars
          ODE_vars(i)=0.0D0
        End do 

        dp_in = (-dlog10(pmin_in)+dlog10(pmax_in))/dble(N_P_BINS-1)

        Do i=1,N_P_BINS
          ODE_p_bins_init(i) = 10.D0**(dlog10(pmin_in)+dble(i-1)*dp_in)
          ODE_p_bins(i) = ODE_p_bins_init(i) 
        End do

        end subroutine ! ODE_Vars_Init!}}}
      !-------------------------------------------------------------

      !-------------------------------------------------------------
        subroutine ODE_Driver_Simple!{{{
        ! driver routine for integrating the p.s.ds using the saved
        ! values of mass, mixing angle and lepton asymmetry and from 
        ! parameters in default file
          external :: derivstemp

          call ODE_Vars_Init
          call lsodeint(ODE_vars, ODE_nvars, T_HIGH, T_LOW, D_TOL,      &
     &                  INTFLAG, derivstemp)
          return
        end subroutine ! ODE_Driver_Simple!}}}
      !-------------------------------------------------------------

      !-------------------------------------------------------------
        subroutine ODE_Driver_Brent!{{{

        external :: zeroin
        double precision :: zeroin
        double precision :: L_sol

        L_sol = zeroin(L_LOW, L_HIGH, ODE_closure, L_TOL)
        ODE_leptasymi = L_sol
        return
        end subroutine ! ODE_Driver_Brent!}}}
      !-------------------------------------------------------------

      !-------------------------------------------------------------
        function ODE_closure(L) result (Value)!{{{
        ! Function that equals zero when the lepton asymmetry is the 
        ! right value to produce the planck closure fraction
          double precision, intent(IN) :: L
          double precision :: Value

          double precision :: omegas, omegasbar

          call ODE_Init_From_Params(l_in = L) ! Set initial leptasym 
          call ODE_Driver_Simple ! Run with this leptasym
          call ODE_sterint(ODE_save_count, omegas, omegasbar)
          Value = omegas + omegasbar - CPARS_omegadmh2
          return
        end function ! closure!}}}
      !-------------------------------------------------------------

      !-------------------------------------------------------------
        subroutine ODE_sterint(ind, omegas, omegasbar)!{{{
        ! Computes integrals over the sterile nu p.s.ds from the
        ! ind^th stored state of the ODE integrator. In order to map 
        ! to z=0, it assumes that it is called above the epoch of weak
        ! decoupling i.e. above ~1 MeV, but below the quark hadron 
        ! transition, so that the neutrino temperature ~1/a. Stores 
        ! \omega_s h^2 in omegas, \omega_sbar h^2 in omegasbar,

          integer, intent(IN) :: ind
          double precision, intent(OUT) :: omegas, omegasbar

          ! Local variables
          ! h^2. Set to 1 for omegah^2
          double precision, parameter :: hub02 = 1.0D0

          double precision :: p_bin_z0(N_P_BINS), p2_bin_z0(N_P_BINS)
          double precision :: rhois(N_P_BINS), rhoisbar(N_P_BINS)
          double precision :: rho_s, rho_sbar, temp
          double precision, external :: NIntegrate

          if ((ind.le.0).or.(ind.gt.ODE_save_count)) Then
            stop "sterint called with invalid index."
          end if

          temp = ODE_T_save(ind)

          ! physical momenta today 
          ! p_today in MeV = (p/T_nu)*T_nu_today
          p_bin_z0 = (CPARS_TcmbMeV*CPARS_nfactor/temp)*                &
     &                ODE_p_bin_save(1:N_P_BINS, ind)
          p2_bin_z0 = p_bin_z0**2

          ! Integrands for energy densities today in units of MeV^4, if 
          ! production were to cease below this temperature
          rhois = (0.5d0/pi**2)*p2_bin_z0*sqrt(p2_bin_z0 + ODE_ms**2)*  &
     &             ODE_y_save(4:N_P_BINS+3, ind)
          rhoisbar = (0.5d0/pi**2)*p2_bin_z0*                           &
     &               sqrt(p2_bin_z0 + ODE_ms**2)*                       &
     &               ODE_y_save(N_P_BINS+4:2*N_P_BINS+3, ind)
          rho_s = NIntegrate(p_bin_z0, rhois, N_P_BINS, p_bin_z0(1),    &
     &                       p_bin_z0(N_P_BINS))
          rho_sbar = NIntegrate(p_bin_z0, rhoisbar, N_P_BINS,           &
     &                          p_bin_z0(1), p_bin_z0(N_P_BINS))

          ! \Omega_s h^2 = rho_s/(\rho_c/h^2)
          omegas = rho_s/CPARS_rhocrit 
          omegasbar = rho_sbar/CPARS_rhocrit
          return
        end subroutine ! ODE_sterint!}}}
      !-------------------------------------------------------------

      !-------------------------------------------------------------
        subroutine lsodeint(ystart, ny, T1, T2, eps, outflag, derivs)!{{{
        ! Subroutine for integration using lsode
          implicit none

          integer :: ny
          double precision :: ystart(ny), T1, T2, eps
          logical :: outflag
          external :: derivs

          ! Local variables
          double precision :: y(ny), T, TOUT, RTOL, ATOL 
          integer          :: ITOL, ITASK, ISTATE, IOPT, LRW, LIW, MF
          ! use non-stiff method
          integer          :: IWORK(20)
          double precision :: RWORK(20+ny*16)
          integer          :: i, nstp, kmax
          ! Should actually be external/interface for the Jacobian,
          ! but it fails during link-time unless I use this.
          double precision :: JDUM

          ! LSODE compulsory parameters and flags
          ITOL   = 1            ! Using scalar absolute tolerance
          RTOL   = eps          ! Relative tolerance
          ATOL   = 1.0D-30      ! Absolute tolerance = 0, 
                                ! set small for f=0
          ITASK  = 4            ! Don't overshoot
          ISTATE = 1            ! First call, so initialize
          IOPT   = 1            ! We will be using optional inputs
          LRW    = 20+ny*16     ! Length of work array
          LIW    = 20           ! Length of integer work array
          MF     = 10           ! use non-stiff method

          ! Run parameters
          ODE_save_count = 0
          T    = T1
          TOUT = T       ! Always begin by outputing the starting state

          ! LSODE optional inputs
          Do i=5, 10
            RWORK(i) = 0.0D0
            IWORK(i) = 0
          end do
          ! Don't overshoot T2, which is lesser than T1
          RWORK(1) = MAX(TOUT, T2)
          IWORK(6) = MAXSTP     ! Maximum number of steps

          if ((outflag) .and. (NSAVE_MAX > 1)) then 
            kmax = NSAVE_MAX - 1 ! Output at 100 log spaced temperatures.
            ODE_dx_save = abs(log(T2)-log(T1))/kmax
          else 
            kmax = 0
          end if 

          ! Copy into interface variables
          do i=1, ny 
            y(i) = ystart(i)
          end do

          !---------------------------------------------------------
          ! Boilerplate code to output variables, here initial state
          ! Perform time evolution
          Call DLSODE(derivs, ny, y, T, TOUT, ITOL, RTOL, ATOL, ITASK,  &
     &                ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JDUM, MF) 
          ! Check correct execution
          IF (ISTATE < 0) GO TO 80
          ! Copy into integrator's module variables
          ODE_save_count = ODE_save_count + 1
          ODE_T_save(ODE_save_count) = T
          do i=1, ODE_nvars
            ODE_y_save(i, ODE_save_count) = y(i)
          end do
          do i=1, N_P_BINS
            ODE_p_bin_save(i, ODE_save_count) = ODE_p_bins(i)
          end do
          !---------------------------------------------------------

          if (kmax.eq.0) then 
            TOUT = T2 ! Do not output intermediate states
            ! Don't overshoot end
            RWORK(1) = max(TOUT, T2)
            ! Perform time evolution
            call DLSODE(derivs, ny, y, T, TOUT, ITOL, RTOL, ATOL,       &
     &                  ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW,    &
     &                  JDUM, MF)
            ! Check correct execution
            if (ISTATE < 0) go to 80
            ! Copy into integrator's variables
            ODE_save_count = ODE_save_count + 1
            ODE_T_save(ODE_save_count) = T
            do i=1, ODE_nvars
               ODE_y_save(i, ODE_save_count) = y(i)
            end do
            do i=1, N_P_BINS
               ODE_p_bin_save(i, ODE_save_count) = ODE_p_bins(i)
            end do
          else ! Output intermediate states
            do nstp = 1, kmax
              TOUT = TOUT*exp(-ODE_dx_save)
              RWORK(1) = max(TOUT, T2)
              call DLSODE(derivs, ny, y, T, TOUT, ITOL, RTOL,           &
     &                    ATOL, ITASK, ISTATE, IOPT, RWORK, LRW,        &
     &                    IWORK, LIW, JDUM, MF)
              ! Check correct execution
              if (ISTATE .LT. 0) go to 80
              ! Copy into integrator's variables
              ODE_save_count = ODE_save_count + 1
              ODE_T_save(ODE_save_count) = T
              do i=1, ODE_nvars
                ODE_y_save(i, ODE_save_count) = y(i)
              end do
              do i=1, N_P_BINS
                ODE_p_bin_save(i, ODE_save_count) = ODE_p_bins(i)
              end do
              ISTATE = 2 ! It's no longer the first call
            end do
          end if

          ! Output inspection of LSODE run to screen if successful
          WRITE(6,60) IWORK(11),IWORK(12),IWORK(13),IWORK(19),RWORK(15)
 60       FORMAT(' No. steps =',I8,'  No. f-s =',I8,'  No. J-s =',I8,   &
     &           ' Method last used =',I2,'   Last switch was at T=',   &
     &           D12.4)
          return
      
          ! Error message if run fails
 80       WRITE(6,90) ISTATE
 90       FORMAT(///' Error halt.. ISTATE =',I3)
          return
        end subroutine ! lsodeint!}}}
      !-------------------------------------------------------------

      !-------------------------------------------------------------
        subroutine ODE_output(ind)!{{{
        ! (re)creates and populates an output directory with the 
        ! state of the ode integrator. (Over)writes a param file, and 
        ! a file with p, f_\nu_s, f_\nu_sbar. Input the required index.
        ! <t> For now, works only on a unix env </t>

          integer, intent(IN) :: ind ! Index in list of saved states
    
          ! Local variables
          double precision :: temp, lkev
          double precision :: omeganet, omegas, omegasbar
          Logical          :: state ! check if state file exists
          ! File units
          integer, parameter  :: fppar=25, fpint=26, fpstate=27
          character(LEN=9), parameter  :: dmake = 'mkdir -p'
          character(LEN=64), parameter :: parname = 'params.dat'
          character(LEN=64), parameter :: statename = 'state.dat'
          character*64     :: dname, fname

          if (ODE_save_count==0) then
            stop "ODE_output called without any saved values."
          end if

          if ((ind<=0).or.(ind>ODE_save_count)) then
            stop "ODE_output called with invalid index."
          end if

          temp = ODE_T_save(ind) ! Temperature
          ! Evaluate (\Omega_s, sbar) h^2 
          call ODE_sterint(ind, omegas, omegasbar)
          omeganet = omegas + omegasbar

          ! First (re)create output directory
          write(dname, "(A9,A2,1pE9.3E2,A2,1pE9.3E2,A1,1pE9.3E2)")      &
     &    'outfiles/','ms',ODE_ms,'s2',ODE_s2,'L',ODE_leptasymi
          call system(dmake // trim(dname))

          ! Then (over)write a param file
          ! Initial lepton asymm in Kev's units
          lkev = ODE_leptasymi*pi**2/(2.0D0*zeta3)
          Open(unit=fppar,file=trim(dname)//'/'//trim(parname),         &
     &         status='unknown')
          Write(fppar,"(1pE24.16,A6)")  ODE_ms, ' # m_s'
          Write(fppar,"(1pE24.16,A16)") ODE_s2, ' # \sin^2 \theta'
          Write(fppar,"(1pE24.16,A12)") lkev, ' # L/n_gamma'
          Write(fppar,"(1pE24.16,A17)") omeganet, ' # \Omega_wdm h^2'
          Write(fppar,"(1pE24.16,A15)") omegas, ' # \Omega_s h^2'
          Write(fppar,"(1pE24.16,A18)") omegasbar,' # \Omega_sbar h^2'
          Close(unit=fppar)

          ! Then output the PSDs
          Write(fname, "(A8,1I3.3,A4)") 'Snapshot',ind,'.dat'
          Open(unit=fpint,file=trim(dname)//'/'//trim(fname),           &
     &         status='unknown')
          call writearray(ODE_p_bin_save, ODE_y_save, N_P_BINS, 2,      &
     &                    (ind-1)*NMAX, (ind-1)*NMAX+3,                 &
     &            'p(MeV) \t \delta f_\\nu \t \delta f_\\nubar', fpint)
          close(unit=fpint)

          ! Then append/output auxiliary run variables into state file
          ! Current Lepton asymm in Kev's units
          lkev = ODE_y_save(3,ind)*pi**2/(2.0D0*zeta3) 
          ! Check if file already exists, as we are appending to it
          Inquire(file=trim(dname)//"/"//trim(statename), exist=state)
          if (state) then
            open(unit=fpstate, file=trim(dname)//"/"//trim(statename),  &
     &           status='old',position='append',action="write")
          else
            open(unit=fpstate,file=trim(dname)//"/"//trim(statename),   &
     &           status="new",action="write")
          ! Write header
            write(fpstate,"(A1, A23,$)") '!','T (MeV)'
            write(fpstate,"(4A24)") 't (sec)', 'ln(a/a_i)',             &
     &                              'L/n_gamma', '\Omega_wdm h^2'
          end if
          write(fpstate,"(1p 5E24.16)") temp, ODE_y_save(1,ind),        &
     &                         ODE_y_save(2,ind), lkev, omeganet
          close(unit=fpstate)
    
          ! Also write to stdout
          print *,'T=', temp, 't=', ODE_y_save(1,ind), 'ln(a/a_i)=',    &
     &             ODE_y_save(2,ind), 'L/T^3=', ODE_y_save(3,ind),      &
     &             'Omega_wdm h^2=',omeganet
          return
        end subroutine ! ODE_output!}}}
        !-------------------------------------------------------------
      end module ! ODE!}}}
