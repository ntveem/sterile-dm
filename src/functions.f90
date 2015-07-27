! functions needed by the code. 
! This file contains the following modules:
! 1. DArrays   -- Reads in data from the data files specified in the 
!                 parameter file. It initializes data structures that
!                 will be used by the splines in the functions.
! 2. functions -- Defines a number of functions required in the 
!                 production calculation. It takes in data from dataIO 
!                 and computes splines for 
!                  a) relativistic degrees of freedom
!                  b) Stefan-Boltzmann particle number suscpetibilities
!                  c) Fermi-Dirac energy densities
!                  d) strong fluid's suceptibilities
!                  e) lepton asymmetry redistribution functions
!                  f) neutrino opacities
!                 and includes a function for the asymmetry potential 
!---------------------------------------------------------------------
      module DArrays!{{{
      !---------------------------------------------------------------
      use Defaults
      use IniFile
      use Data_Repo
      implicit none
      public

      ! Filenames!{{{
      character(LEN=INI_MAX_STRING_LEN) :: gstarfile, pchifile, EFfile, &
     &                                     chifile, ratefile_mu,        &
     &                                     redistfile_e, redistfile_mu, &
     &                                     redistfile_tau

      ! parameters for data structures to be used by physics functions

      ! Relativistic d.o.f -- two for g_* and g_*,s
      integer :: ngstar, ngstar_f
      double precision, pointer :: ltempgstar_t(:), gstar_t(:,:)

      ! Particle number susceptibility (Stefan-Boltzmann)
      integer :: nlpchi, nlpchi_f
      double precision, pointer :: masslpchi_t(:), lpchi_t(:,:)

      ! Energy density for FD distribution
      integer :: nlEF, nlEF_f
      double precision, pointer :: masslEF_t(:), lEF_t(:,:)

      ! Strongly interacting fluid's susceptibilities
      integer :: nscps, nscps_f
      ! three for B, Q and QB
      double precision, pointer :: temps_t(:), scps_t(:,:)

      ! asymmetry redistribution functions
      integer :: nscpe, nscpmu, nscptau
      integer :: nscpe_f, nscpmu_f, nscptau_f
      double precision, pointer :: tempe_t(:), tempmu_t(:), temptau_t(:)
      double precision, pointer :: scpe_t(:,:), scpmu_t(:,:),           &
     &                             scptau_t(:,:)

      ! neutrino opacities
      integer :: nlpTrate, nlTrate 
      double precision, pointer :: lpTrate_t(:), lTrate_t(:),           &
     &                             lrate_t(:,:)

      ! parameters of the interpolation scheme
      double precision :: qcdtemp, Tc!}}}

      save
      contains
      !------------------------------------------------------------
        subroutine DArrays_Read_Files!{{{
        ! Reads in appropriate data file names from param file and 
        ! reads in their contents into the repo

          ! Read Filenames
          gstarfile = Ini_Read_String('gstarfile')
          pchifile = Ini_Read_String('pchifile')
          EFfile = Ini_Read_String('EFfile')
          chifile = Ini_Read_String('chifile')
          ratefile_mu = Ini_Read_String('ratefile_mu')
          redistfile_e = Ini_Read_String('redistfile_e')
          redistfile_mu = Ini_Read_String('redistfile_mu')
          redistfile_tau = Ini_Read_String('redistfile_tau')

          ! Read in data in files into Data repository
          call Data_Open(gstarfile, 80, 'gstar')
          call Data_Open(pchifile, 80, 'pchi')
          call Data_Open(EFfile, 80, 'EF')
          call Data_Open(chifile, 80, 'chi')
          ! Read asymmetry functions
          call Data_Open(redistfile_e, 80, 'dmudLe')
          call Data_Open(redistfile_mu, 80, 'dmudLmu')
          call Data_Open(redistfile_tau, 80, 'dmudLtau')
          ! Also read in temperatures from second line in rate file
          call Data_Open(ratefile_mu, 80, 'rates_mu', 2)
        end subroutine ! DArrays_Read_Files!}}}
      !------------------------------------------------------------

      !------------------------------------------------------------
        subroutine DArrays_Init!{{{
        ! reads in the data from the repository into appropriate
        ! arays for use by the physics functions

        integer :: i, j, k 

        ! Read the Relativistic d.o.f
        ngstar = Data_Book_NIndep('gstar')
        ngstar_f = Data_Book_NDep('gstar')
        allocate(ltempgstar_t(ngstar))
        allocate(gstar_t(ngstar,ngstar_f))
        call Data_Read('gstar', ltempgstar_t, ngstar, i,                &
                                gstar_t, ngstar, ngstar_f, j)

        ! Read the Stefan-Bolztmann susceptibilities
        nlpchi = Data_Book_NIndep('pchi')
        nlpchi_f = Data_Book_NDep('pchi')
        allocate(masslpchi_t(nlpchi))
        allocate(lpchi_t(nlpchi,nlpchi_f))
        call Data_Read('pchi', masslpchi_t, nlpchi, i,                  &
     &                         lpchi_t, nlpchi, nlpchi_f, j)

        ! Read the energy densities for FD distributions
        nlEF = Data_Book_NIndep('EF')
        nlEF_f = Data_Book_NDep('EF')
        allocate(masslEF_t(nlEF))
        allocate(lEF_t(nlEF,nlEF_f))
        call Data_Read('EF', masslEF_t, nlEF, i, lEF_t, nlEF, nlEF_f, j)

        ! Read the strong fluid's susceptibilities
        nscps = Data_Book_NIndep('chi')
        nscps_f = Data_Book_NDep('chi')
        allocate(temps_t(nscps))
        allocate(scps_t(nscps,nscps_f))
        call Data_Read('chi', temps_t, nscps, i,                        &
     &                        scps_t, nscps, nscps_f, j)

        ! Read in the redistribution functions
        nscpe = Data_Book_NIndep('dmudLe')
        nscpe_f = Data_Book_NDep('dmudLe')
        nscpmu = Data_Book_NIndep('dmudLmu')
        nscpmu_f = Data_Book_NDep('dmudLmu')
        nscptau = Data_Book_NIndep('dmudLtau')
        nscptau_f = Data_Book_NDep('dmudLtau')
        allocate(tempe_t(nscpe))
        allocate(scpe_t(nscpe,nscpe_f))
        allocate(tempmu_t(nscpmu))
        allocate(scpmu_t(nscpmu,nscpmu_f))
        allocate(temptau_t(nscptau))
        allocate(scptau_t(nscptau,nscptau_f))
        call Data_Read('dmudLe', tempe_t, nscpe, i,                     &
     &                           scpe_t, nscpe, nscpe_f, j)
        call Data_Read('dmudLmu', tempmu_t, nscpmu, i,                  &
     &                            scpmu_t, nscpmu, nscpmu_f, j)
        call Data_Read('dmudLtau', temptau_t, nscptau, i,               &
     &                             scptau_t, nscptau, nscptau_f, j)

        ! Read in the neutrino opacities, with extra storage for 
        ! temperatures
        nlpTrate = Data_Book_NIndep('rates_mu')
        nlTrate = Data_Book_NDep('rates_mu')
        allocate(lpTrate_t(nlpTrate))
        allocate(lrate_t(nlpTrate,nlTrate))
        allocate(lTrate_t(nlTrate))
        call Data_Read('rates_mu', lpTrate_t, nlpTrate, i,              &
     &                             lrate_t, nlpTrate, nlTrate, j,       &
     &                             lTrate_t, nlTrate, k)
        if (k /= nlTrate) then
           print *,"Given rates at ",nlTrate," temperatures, but",      &
     &             " given ",k," temperatures"
           stop
        end if
        ! Read in the parameters of the rate interpolation scheme
        qcdtemp = Ini_Read_double('qcdtemp')
        Tc = Ini_Read_double('Tc')
        end subroutine ! DArrays_Init!}}}
      !------------------------------------------------------------

      !------------------------------------------------------------
        subroutine DArrays_Close!{{{
          ! Frees up memory
          deallocate(ltempgstar_t)
          deallocate(gstar_t)

          deallocate(masslpchi_t)
          deallocate(lpchi_t)

          deallocate(masslEF_t)
          deallocate(lEF_t)

          deallocate(temps_t)
          deallocate(scps_t)

          deallocate(tempe_t)
          deallocate(tempmu_t)
          deallocate(temptau_t)
          deallocate(scpe_t)
          deallocate(scpmu_t)
          deallocate(scptau_t)

          deallocate(lpTrate_t)
          deallocate(lTrate_t)
          deallocate(lrate_t)
        end subroutine ! DArrays_Close!}}}
      !------------------------------------------------------------
      !---------------------------------------------------------------
      end module ! DArrays!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
     module Functions!{{{
      !---------------------------------------------------------------
      Use Phconst
      Use DArrays
      Use Defaults ! maximum p bins needed for working space for
                   ! spline interpolation
      Implicit none
      Public

      ! Degrees of splines!{{{
      integer, parameter :: kgstar = 3
      integer, parameter :: klpchi = 3
      integer, parameter :: klEF = 3
      integer, parameter :: kscps = 3
      integer, parameter :: kscpe = 3
      integer, parameter :: kscpmu = 3
      integer, parameter :: kscptau = 3
      integer, parameter :: klpTrate = 3
      integer, parameter :: klTrate = 3

      ! numbers of knots and arrays of spline knots and coefficients
      integer :: ngstarknots
      double precision, pointer :: gstarknots(:,:), gstarcoeffs(:,:)

      integer :: nlpchiknots
      double precision, pointer :: lpchiknots(:,:), lpchicoeffs(:,:)

      integer :: nlEFknots
      double precision, pointer :: lEFknots(:,:), lEFcoeffs(:,:)

      integer :: nscpsknots
      double precision, pointer :: scpsknots(:,:), scpscoeffs(:,:)

      integer :: nscpeknots, nscpmuknots, nscptauknots
      double precision, pointer :: scpeknots(:,:), scpecoeffs(:,:)
      double precision, pointer :: scpmuknots(:,:), scpmucoeffs(:,:)
      double precision, pointer :: scptauknots(:,:), scptaucoeffs(:,:)

      ! rates are special because we will possibly truncate the array
      double precision :: lpTmin, lpTmax, lTmin, lTmax
      ! parameters of the interpolation scheme
      double precision :: hightcutoff
      double precision, parameter :: delta = 1.0D-5
      ! parameters of the truncated array
      integer :: nlTrate_sub
      double precision, pointer :: lTrate_t_sub(:), lrate_t_sub(:)
      ! nlTrateknots and nlpTrateknots are the numbers of knots for a 
      ! 2D spline of log(scaled rate) vs log(T/MeV) and log(p/T)  
      ! (coordinate style, array style is reversed) and with degrees 
      integer :: nlTrateknots,nlpTrateknots
      double precision, pointer :: lTrateknots(:), lpTrateknots(:),     &
     &                             lratecoeffs(:)

      ! Working space for rate spline evaluation
      ! Maximum size of pT by T array we will be evaluating over
      integer, parameter :: nlpTrateevalmax = N_P_BINS
      integer, parameter :: nlTrateevalmax = 1

      integer lwrkrateeval, kwrkrateeval
      parameter (lwrkrateeval=nlpTrateevalmax*(klpTrate+1)+             &
     &           nlTrateevalmax*(klTrate+1))
      parameter (kwrkrateeval=nlpTrateevalmax+nlTrateevalmax)

      double precision wrkrateeval(lwrkrateeval)
      integer iwrkrateeval(kwrkrateeval)!}}}

      save
      contains
      !------------------------------------------------------------
        subroutine Functions_Splines_Init!{{{
        ! Initialize a number of splines that are used by the 
        ! functions

          integer :: i, j
          integer :: indlowT, indhighT

          ! Splines for g_* and g_*,s vs log(T)
          allocate(gstarknots(ngstar+kgstar+1,ngstar_f))
          allocate(gstarcoeffs(ngstar+kgstar+1,ngstar_f))
          ! Logscale temperatures for gstar
          do i=1, ngstar
            ltempgstar_t(i) = log(ltempgstar_t(i))
          end do
          do i=1, ngstar_f
            Call SplFit(ltempgstar_t, gstar_t(1,i), ngstar,             &
     &                  gstarknots(1,i), gstarcoeffs(1,i), ngstarknots)
          end do

          ! Spline for log(Stefan-Boltzmann susceptibilities/T^2) vs m/T
          allocate(lpchiknots(nlpchi+klpchi+1,nlpchi_f))
          allocate(lpchicoeffs(nlpchi+klpchi+1,nlpchi_f))
          ! Logscale susceptibilities for SB
          do i=1, nlpchi
            do j=1, nlpchi_f
              lpchi_t(i,j) = log(lpchi_t(i,j))
            end do
          end do
          do i=1, nlpchi_f
            Call SplFit(masslpchi_t, lpchi_t(1,i), nlpchi,              &
     &                  lpchiknots(1,i), lpchicoeffs(1,i), nlpchiknots)
          end do

          ! Spline for log(E_Fermi-Dirac/T^4) vs m/T
          allocate(lEFknots(nlEF+klEF+1,nlEF_f))
          allocate(lEFcoeffs(nlEF+klEF+1,nlEF_f))
          ! Logscale energy densities for FD
          do i=1, nlEF
            do j=1, nlEF_f
              lEF_t(i,j) = log(lEF_t(i,j))
            end do
          end do
          do i=1, nlEF_f
            Call SplFit(masslEF_t, lEF_t(1,i), nlEF,                    &
     &                  lEFknots(1,i), lEFcoeffs(1,i), nlEFknots)
          end do

          ! Splines for \chi^B, \chi^Q and \chi^QB vs T
          allocate(scpsknots(nscps+kscps+1,nscps_f))
          allocate(scpscoeffs(nscps+kscps+1,nscps_f))
          do i=1, nscps_f
            Call SplFit(temps_t, scps_t(1,i), nscps,                    &
     &                  scpsknots(1,i), scpscoeffs(1,i), nscpsknots)
          end do

          ! Splines for dmu/dL vs T
          allocate(scpeknots(nscpe+kscpe+1,nscpe_f))
          allocate(scpmuknots(nscpmu+kscpmu+1,nscpmu_f))
          allocate(scptauknots(nscptau+kscptau+1,nscptau_f))
          allocate(scpecoeffs(nscpe+kscpe+1,nscpe_f))
          allocate(scpmucoeffs(nscpmu+kscpmu+1,nscpmu_f))
          allocate(scptaucoeffs(nscptau+kscptau+1,nscptau_f))
          do i=1, nscpe_f
            Call SplFit(tempe_t, scpe_t(1,i), nscpe,                    &
     &                  scpeknots(1,i), scpecoeffs(1,i), nscpeknots)
          end do
          do i=1, nscpmu_f
            Call SplFit(tempmu_t, scpmu_t(1,i), nscpmu,                 &
     &                  scpmuknots(1,i), scpmucoeffs(1,i), nscpmuknots)
          end do
          do i=1, nscptau_f
            Call SplFit(temptau_t, scptau_t(1,i), nscptau,              &
     &               scptauknots(1,i), scptaucoeffs(1,i), nscptauknots)
          end do

          ! A table of log(scaled rate) vs log(p/T) and log(T/MeV). 
          ! z vs x and y is in coordinate style.
          ! First define lower temp cutoff for the interpolation scheme
          if ((qcdtemp < Tc) .and. (Tc <= (qcdtemp + delta))) then
            hightcutoff = qcdtemp + delta
          Else if (Tc > (qcdtemp + delta)) then
            hightcutoff = Tc - delta
          Else
            Print *,"Tc=",Tc
            Print *,"qcdtemp=",qcdtemp
            Print *,"Need Tc>qcdtemp"
            stop
          end if
          ! Next def. subset for cutoff rates
          ! Inefficient, but traverse T array and find temperature limits
          indlowT = 0
          indhighT = 1
          do i=1, nlTrate
             if (lTrate_t(i).le.qcdtemp) indlowT = indlowT + 1
             if (lTrate_t(i).lt.hightcutoff) indhighT = indhighT + 1
          end do
          nlTrate_sub = indlowT + nlTrate - indhighT + 1
          ! Allocate space for temperature and transposed array 
          ! (latter for compatibility reasons)
          allocate(lTrate_t_sub(nlTrate_sub))
          allocate(lrate_t_sub(nlTrate_sub*nlpTrate))
          ! Now copy subset of working array, traversing it row-wise
          ! and populating column-wise
          do i=1, nlpTrate
             ! First below the transition
             do j=1, indlowT
                lTrate_t_sub(j) = lTrate_t(j) ! many times...
                lrate_t_sub((i-1)*nlTrate_sub+j) = lrate_t(i,j)
             end do
             ! then above the cutoff
             do j=indhighT, nlTrate
                lTrate_t_sub(indlowT+j-indhighT+1) = lTrate_t(j)
                lrate_t_sub((i-1)*nlTrate_sub+indlowT+j-indhighT+1) =   &
     &                                                      lrate_t(i,j)
             end do
          end do
          ! Allocate space for spline knots and coefficients
          allocate(lpTrateknots(nlpTrate + klpTrate + 1))
          allocate(lTrateknots(nlTrate_sub + klTrate + 1))
          allocate(lratecoeffs(nlpTrate*nlTrate_sub))
          ! Scale rates w/ p. and take logs because the splines 
          ! are of log(scaled rate) vs log(p/T) and log(T)
          ! Traverse in column major 
          do i=1, nlpTrate
             do j=1, nlTrate_sub
                lrate_t_sub((i-1)*nlTrate_sub+j) =                      &
     &            log(lrate_t_sub((i-1)*nlTrate_sub+j)/lpTrate_t(i))
             end do
          end do
          do i=1, nlpTrate
             lpTrate_t(i) = log(lpTrate_t(i))
          end do
          do i=1, nlTrate_sub
             lTrate_t_sub(i) = log(lTrate_t_sub(i))
          end do
          ! Bounds, needed to check 2d spline evaluation
          lpTmin = minval(lpTrate_t)
          lpTmax = maxval(lpTrate_t)
          lTmin = minval(lTrate_t_sub)
          lTmax = maxval(lTrate_t_sub) 
          ! Find spline
          Call TwoDimSplFit(lpTrate_t, lTrate_t_sub, lrate_t_sub,       &
     &                      nlpTrate, nlTrate_sub, lpTrateknots,        &
     &                      lTrateknots, lratecoeffs, nlpTrateknots,    &
     &                      nlTrateknots)
        end subroutine ! functions_Splines_Init!}}}
      !------------------------------------------------------------

      !------------------------------------------------------------
        subroutine Functions_Splines_Close!{{{
        ! Free up memory

          deallocate(gstarknots)
          deallocate(gstarcoeffs)

          deallocate(lpchiknots)
          deallocate(lpchicoeffs)

          deallocate(lEFknots)
          deallocate(lEFcoeffs)

          deallocate(scpsknots)
          deallocate(scpscoeffs)

          deallocate(scpeknots)
          deallocate(scpmuknots)
          deallocate(scptauknots)
          deallocate(scpecoeffs)
          deallocate(scpmucoeffs)
          deallocate(scptaucoeffs)

          deallocate(lTrate_t_sub)
          deallocate(lrate_t_sub)

          deallocate(lpTrateknots)
          deallocate(lTrateknots)
          deallocate(lratecoeffs)
        end subroutine ! functions_Splines_Close!}}}
      !------------------------------------------------------------

      !------------------------------------------------------------
        function gstar(T, alpha, beta)!{{{
        ! Evaluates (d/dT)^\beta(g_*/g_*,s) for alpha=1/2, as a function
        ! of temperature (in MeV) from precalculated spline knots and 
        ! coefficients, for beta = 0 or 1. 
        ! The derivative is in units of MeV^{-1}.

          double precision :: T
          integer :: alpha, beta
          double precision :: gstar
          ! Internal variables
          double precision, pointer :: wrk(:)
          double precision :: logT, der, val
          integer :: ier
          integer, parameter :: nvals = 1

          allocate(wrk(ngstarknots))
          logT = log(T)

          if ((alpha > ngstar_f) .or. (alpha <= 0)) then
            Print *,"alpha = ",alpha
            Print *,"gstar(T,alpha,beta) defined only for alpha <= ",   &
     &            ngstar_f
            stop
          end if

          if (beta.eq.0) then
                  call splev(gstarknots(1,alpha), ngstarknots,          &
     &                       gstarcoeffs(1,alpha), kgstar, logT,        &
     &                       val, nvals, ier)
          Else if (beta.eq.1) then
                  call splder(gstarknots(1,alpha), ngstarknots,         &
     &                        gstarcoeffs(1,alpha), kgstar, beta,       &
     &                        logT, der, nvals, wrk, ier)
                  val = (1.0d0/T)*der
          Else
                  Print *,beta
          ! As a curiosity, the general formula for converting 
          ! derivatives in terms of Log(x) to those in x involves 
          ! Stirling Numbers of the First kind :)
                  Stop "function gstar only evaluates function or       &
     &                  its first derivative, i.e., beta = 0 or 1."
          end if

          if (ier.gt.0) then
                  Print *,"ier=",ier
                  Print *,"T=",T
                  Stop "Error in function gstar, check if evaluated     &
     &                  outside temperature range"
          end if

          gstar = val
          deallocate(wrk)
          Return
        end function ! gstar!}}}
      !------------------------------------------------------------

      !------------------------------------------------------------
        function pchi(m)!{{{
        ! Evaluate particle number susceptibility for a Fermi-Dirac
        ! distribution with degeneracy factor of one, from pre-computed
        ! spline knots and coefficients for log(chi/T^2). 
        ! Returns chi/T^2 given m/T

          double precision :: m
          double precision :: pchi
          ! Internal variables
          double precision :: val
          integer :: ier
          integer, parameter :: nvals = 1

          CALL splev(lpchiknots(1,1), nlpchiknots, lpchicoeffs(1,1),    &
     &               klpchi, m, val, nvals, ier)

          if (ier.gt.0) then
                  Print *,"ier=",ier
                  Print *,"m/T=",m
                  Stop "Error in function pchi, check if evaluated      &
     &                  outside temperature range."
          end if

          pchi = exp(val)
          Return
        end function ! pchi!}}}
      !------------------------------------------------------------

      !------------------------------------------------------------
        function EF(m)!{{{
        ! Evaluate energy density for a Fermi-Dirac distribution
        ! with degeneracy factor of one, from pre-computed spline knots 
        ! and coefficients for log(E_Fermi Dirac/T^4). 
        ! Returns E_F/T^4 given m/T

          double precision :: m
          double precision :: EF
          ! Internal variables
          double precision :: val
          integer :: ier
          integer, parameter :: nvals = 1
   
          CALL splev(lEFknots(1,1), nlEFknots, lEFcoeffs(1,1),          &
     &               klEF, m, val, nvals, ier)

          if (ier.gt.0) then
                  Print *,"ier=",ier
                  Print *,"m/T=",m
                  Stop "Error in function EF, check if evaluated        &
     &                  outside temperature range."
          end if

          EF = exp(val)
          Return
        end function ! EF!}}}
      !------------------------------------------------------------

      !------------------------------------------------------------
        function scps(T, x)!{{{
        ! Evaluate susceptibilities \chi^B_2, \chi^Q_2 and \chi^BQ_11
        ! of the strong fluid (in units of T^2) as a function of temp
        ! (in MeV) from pre-calculated spline knots and coefficients.

          double precision :: T
          integer :: x ! \Chi^B_2, \chi^Q_2, \chi^BQ_11 for x = 1, 2, 3
          double precision :: scps
          ! Internal variables
          double precision :: val
          integer :: ier ! number of values and error flag
          integer, parameter :: nvals = 1

          if ((x > nscps_f) .or. (x <= 0)) then
            Print *,"x = ",x
            Print *,"scps(T,x) defined only for x <= ",nscps_f
            stop
          end if
    
          CALL splev(scpsknots(1,x), nscpsknots, scpscoeffs(1,x),       &
     &               kscps, T, val, nvals, ier)

          if (ier.gt.0) then
                  Print *,"ier=",ier
                  Print *,"T=",T
                  Print *,"x=",x
                  Stop "Error in function scps, check if evaluated      &
     &                  outside temperature range."
          end if

          scps = val
          Return
        end function ! scps!}}}
      !------------------------------------------------------------

      !------------------------------------------------------------
        function dmudlx(T, alpha, X)!{{{
        ! Evaluates dmu_{\alpha}/dL_X in units of T^{-2} as a function
        ! of temperature T (in MeV), from precalculated spline knots &
        ! coefficients. \alpha = 1, ... 8 runs over e, mu, tau, nu_e, 
        ! nu_mu, nu_tau, Q, B. X = e, mu or tau

          double Precision :: T
          integer :: alpha 
          Character (LEN=*) :: X
          double precision :: dmudlx
          ! Internal variables
          double precision val
          integer :: ier ! flag
          integer, parameter :: nvals = 1 ! Number of values requested
   
          if (X(1:1) .eq. 'e') then
            if ((alpha > nscpe_f) .or. (alpha <= 0)) goto 40
            CALL splev(scpeknots(1,alpha), nscpeknots,                  &
     &                 scpecoeffs(1,alpha), kscpe, T, val, nvals, ier)
          Else if (X(1:1) .eq. 'm') then
            if ((alpha > nscpmu_f) .or. (alpha <= 0)) goto 40
            CALL splev(scpmuknots(1,alpha), nscpmuknots,                &
     &                 scpmucoeffs(1,alpha), kscpmu, T, val, nvals, ier)
          Else if (X(1:1) .eq. 't') then
            if ((alpha > nscptau_f) .or. (alpha <= 0)) goto 40
            CALL splev(scptauknots(1,alpha), nscptauknots,              &
     &                 scptaucoeffs(1,alpha), kscptau, T, val, nvals,   &
     &                 ier)
          Else
                  Print *,X
                  Stop "function dmudlx needs X = e, mu or tau."
          end if
          if (ier.gt.0) then
                  Print *,"ier=",ier
                  Print *,"T=",T
                  Stop "Error in function dmudlx, check if evaluated    &
     &                  outside temperature range."
          end if

          dmudlx = val
          return

 40       Print *,alpha
          stop "alpha out of range for dmudlx(T,alpha,X)"
        end function ! dmudlx!}}}
      !------------------------------------------------------------

      !------------------------------------------------------------
        function vasym(T, X)!{{{
        ! Evaluates dimensionless factor relating the asymmetry 
        ! potential V_A to G_F*(L_X in T^3), at temperature T (in MeV). 
          double precision :: T
          Character (LEN=*) :: X
          double precision :: vasym

          ! Internal variables
          ! contributions of charged, neutral and strong sectors
          double precision :: vlc, vln, vlq

          ! First the parts which don't depend on flavor
          vlc = sqrt(2.0D0)*(-0.5D0 + 2.0D0*s2w)*(                      &
     &           2.0D0*pchi(me/T)*dmudlx(T,1,X) +                       &
     &           2.0D0*pchi(mmu/T)*dmudlx(T,2,X) +                      &
     &           2.0D0*pchi(mtau/T)*dmudlx(T,3,X) )

          vln = sqrt(2.0D0)*(                                           &
     &           pchi(mnue/T)*dmudlx(T,4,X) +                           &
     &           pchi(mnumu/T)*dmudlx(T,5,X) +                          &
     &           pchi(mnutau/T)*dmudlx(T,6,X) )

          vlq = sqrt(2.0D0)*(1.0D0 - 2.0D0*s2w)*(                       &
     &           scps(T,2)*dmudlx(T,7,X) +                              &
     &           scps(T,3)*dmudlx(T,8,X))

          ! then parts which depend on flavor
          if (X(1:1).eq.'e') then
            vlc = vlc + sqrt(2.0D0)*2.0D0*pchi(me/T)*dmudlx(T,1,X)
            vln = vln + sqrt(2.0D0)*pchi(mnue/T)*dmudlx(T,4,X)
          Else if (X(1:1).eq.'m') then
            vlc = vlc + sqrt(2.0D0)*2.0D0*pchi(mmu/T)*dmudlx(T,2,X)
            vln = vln + sqrt(2.0D0)*pchi(mnumu/T)*dmudlx(T,5,X)
          Else if (X(1:1).eq.'t') then
            vlc = vlc + sqrt(2.0D0)*2.0D0*pchi(mtau/T)*dmudlx(T,3,X)
            vln = vln + sqrt(2.0D0)*pchi(mnutau/T)*dmudlx(T,6,X)
          Else
            Print *,X
            Stop "function vl takes X = e, m or t."
          end if
          vasym = vlc + vln + vlq
          RETURN
        end function !vasym!}}}
      !------------------------------------------------------------

      !------------------------------------------------------------
        subroutine GamFac(npT, p, T, flavor, scrates)!{{{
        ! subroutine returns the scaled scattering rates, 
        ! \Gamma(p)/G_F^2 p T^4 given a set of n input ps and T in MeV 
          integer npT
          double precision p(npT), T, scrates(npT)
          Character*8 flavor
          ! Internal variables
          double precision logpT(npT), logT, val(npT)
          integer i, ier, nlT ! ier=Error flag
          double precision lpTmin_in, lpTmax_in
          parameter (nlT=1)

          if (flavor(1:1).ne.'m') then
                  Print *,"flavor=",flavor
                  Stop "For now we only have the rates for mu."
          end if

          do i=1, npT
             logpT(i) = log(p(i)/T)
          end do
          logT = log(T)

          Call bispev(lpTrateknots,nlpTrateknots,lTrateknots,           &
     &                nlTrateknots,lratecoeffs,klpTrate,klTrate,        &
     &                logpT,npT,logT,nlT,val,wrkrateeval,lwrkrateeval,  &
     &                iwrkrateeval,kwrkrateeval,ier)

          ! Check bounds because bispev doesn't do that for you
          lpTmin_in = minval(logpT)
          lpTmax_in = maxval(logpT)
          if ((lpTmin_in < lpTmin) .or. (lpTmax_in > lpTmax) .or.       &
     &        (logT < lTmin) .or. (logT > lTmax)) ier = 1

          if (ier.gt.0) then
                  Print *,"ier=",ier
                  Print *,"log(p/T)(1)=",logpT(1)
                  Print *,"log(p/T)(npT)=",logpT(npT)
                  Print *,"T=",T
                  Print *,"min(log(p/T))=",lpTmin
                  Print *,"max(log(p/T))=",lpTmax
                  Stop "Error in subroutine gamfac, check if evaluated  &
     &                  outside p/T and T range covered by rate tables. &
     &                  if so change pT(min/max)_in in defaults."
          end if

          do i=1, npT
             scrates(i) = exp(val(i))
          end do
        end subroutine ! GamFac!}}}
      !------------------------------------------------------------
      !---------------------------------------------------------------
      end module ! Functions!}}}
!---------------------------------------------------------------------
