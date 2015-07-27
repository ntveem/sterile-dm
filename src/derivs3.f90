! Function defining the derivative w.r.t temperature of all the
! quantities being evolved. It follows the conventions of the lsode
! solver -- i.e., its arguments are of the form 
!                derivs(nvars, x, y(nvars), dydx(nvars)), 
! where x is the temperature, and y and dydx are arrays of dimension 
! )2*N_P_BINS + 3): 
!                y(1) = FRW coordinate time. 
!                y(2) = number of efolds since T_HIGH = \int H dt
!                y(3) = L/T^3
!                y(4:N_P_BINS+3) = p.s.ds of sterile neutrinos
!                y(N_P_BINS + 4:) = p.s.ds of sterile antineutrinos   
!---------------------------------------------------------------------
      subroutine derivstemp(nvr,temp,y,dydtemp)!{{{
     !---------------------------------------------------------------
      use Phconst
      use Functions
      use Ode
      use Defaults
      implicit None

      integer, intent(IN) :: nvr
      Double precision, intent(IN) :: temp, y(nvr) 
      Double precision, intent(OUT) :: dydtemp(nvr)
    
      ! LOCAL VARIABLES
      ! Parameter definitions
      Double precision :: scalef, leptasym
      Double precision :: mx, mnux, munux ! Flavor dependent params
      Double precision :: ex, enux
      ! potentials and scaled \Gamma
      Double precision :: vl, vtf, sctf(N_P_BINS) 

      ! Local versions of variables on the heap
      Double precision :: pi_l, Gf_l, hbar_l, mp_l
      Double precision :: ms_l, s2_l, dm2_l
      Integer :: n_p_bins_l

      ! Sterile nu integrands, calculated on the stack
      Double precision p_bins_l(N_P_BINS)
      Double precision p2_bins(N_P_BINS), estert(N_P_BINS) 
      Double precision rhosteri(N_P_BINS), rhoster
      ! Production of sterile neutrinos     
      ! Opacities in matter
      Double precision gams(N_P_BINS), gambs(N_P_BINS)
      ! SM neutrino PSDs
      Double precision fnux(N_P_BINS), fnuxbar(N_P_BINS)
      ! Time derivatives
      Double precision dydt(nvr)
      ! Integrands for tder
      Double precision dLdti(N_P_BINS), drhosdti(N_P_BINS) 
      ! Standard model integrands
      Double precision rhoSM, entropySM 
      ! Time and asymmetry evolution
      Double precision hubcst, drhosdt, drhoSMdT, dtdT

      ! V_asymmetry, scattering rate factor, spline intergral
      external :: NIntegrate
      double precision :: NIntegrate
      integer i

      if (nvr.lt.(2*N_P_BINS+3)) then
        stop "derivs called with too few variables."
      end if

      ! Parameter definitions -- populating the stack
      !--------------------------------------------------------------
      scalef = exp(-y(2)) ! (a_start/a_current)
      leptasym = y(3) ! Leptasym in units of T^3

      ! Flavor dependent parameters
      mx = 0.0D0
      mnux = 0.0D0
      munux = 0.0D0
      If (ODE_flavor(1:1).eq.'e') Then
        mx = me
        mnux = mnue
        munux = dmudlx(temp,4,ODE_flavor)*leptasym ! mu/T for \nu_X
      Else if (ODE_flavor(1:1).eq.'m') Then
        mx = mmu
        mnux = mnumu
        munux = dmudlx(temp,5,ODE_flavor)*leptasym ! mu/T for \nu_X
      Else if (ODE_flavor(1:1).eq.'t') Then
        mx = mtau
        mnux = mnutau
        munux = dmudlx(temp,6,ODE_flavor)*leptasym ! mu/T for \nu_X
      Else
        Print *,ODE_flavor
        Stop "Flavor needs to be one of e, m or t"
      End if
      ex = 4.0d0*EF(mx/temp) ! Net energy of charged leptons
      enux = 2.0d0*EF(mnux/temp) ! Net energy of neutral leptons

      ! Asymmetry potential in MeV
      vl = vasym(temp, ODE_flavor)*Gf*leptasym*temp**3
      ! Dimensionless factor relating thermal potential and p
      vtf = -(8.0*sqrt(2.0D0)/3.0)*Gf*temp**4*(enux/mz**2 + ex/mw**2)

      ! Read into stack
      pi_l = pi
      Gf_l = Gf
      mp_l = mp
      hbar_l = hbar

      ms_l = ODE_ms
      dm2_l = (ODE_ms - mnux)**2
      s2_l = ODE_s2
      n_p_bins_l = N_P_BINS

      ! Define Lagrangian momentum bins on stack, and update heap 
      Do i=1, N_P_BINS
        ODE_p_bins(i) = scalef*ODE_p_bins_init(i)
        p_bins_l(i) = ODE_p_bins(i) 
      End do

      ! Dimensionless factors relating scattering rates to G_F^2 p T^4
      Call GamFac(n_p_bins_l, p_bins_l, temp, ODE_flavor, sctf)
      ! Dimensionless factors relating the scattering rate to p 
      Do i=1, n_p_bins_l
        sctf(i) = sctf(i)*(Gf**2)*temp**4
      End do

      ! Calculating df/dt entirely on the stack
      !------------------------------------------------------------------
      !$OMP PARALLEL DO
      Do i=1, n_p_bins_l
        p2_bins(i) = p_bins_l(i)**2
        estert(i) = sqrt(p2_bins(i) + ms_l**2)
        ! Integrand for \rho_sterile
        rhosteri(i) = (0.5d0/pi_l**2)*p2_bins(i)*estert(i)*             &
     &                ( y(i+3) + y(i + 3 + n_p_bins_l) ) 

        ! Sterile production, code to evaluate d/dt(f_\nu, f_\bar{\nu})
        ! SM neutrino phase space densities
        fnux(i) = 1.0d0/                                                &
     &     ( 1.0d0 + exp((sqrt(p2_bins(i) + mnux**2)/temp) - munux))
        fnuxbar(i) = 1.0d0/                                             &
     &     ( 1.0d0 + exp((sqrt(p2_bins(i) + mnux**2)/temp) + munux))

         ! Interaction rates in matter, in MeV
        gams(i) = 0.25d0*sctf(i)*p_bins_l(i)*s2_l/                      &
     &            ( s2_l + (sctf(i)*p2_bins(i)/dm2_l)**2 +              &
     &            ( sqrt(1.0D0 - s2_l) - vtf*2.0d0*p2_bins(i)/dm2_l     &
     &            - vl*2.0d0*p_bins_l(i)/dm2_l)**2)
        gambs(i) = 0.25d0*sctf(i)*p_bins_l(i)*s2_l/                     &
     &             ( s2_l + (sctf(i)*p2_bins(i)/dm2_l)**2 +             &
     &             ( sqrt(1.0D0 - s2_l) - vtf*2.0d0*p2_bins(i)/dm2_l    &
     &             + vl*2.0d0*p_bins_l(i)/dm2_l)**2)

        ! Sterile PSD derivatives w.r.t time
        dydt(i+3) = (1.0d0/hbar_l)*gams(i)*( fnux(i) - y(i+3) )
        dydt(i+3+n_p_bins_l) = (1.0d0/hbar_l)*gambs(i)*                 &
     &                         ( fnuxbar(i) - y(i+3+n_p_bins_l) )

        ! Integrands for leptasym and sterile energy rates of change
        ! solely due to sterile neutrino production
        dLdti(i) = -(0.5d0/(pi_l**2*temp**3))*p2_bins(i)*               &
     &              ( dydt(i+3) - dydt(i+3+n_p_bins_l) )
        drhosdti(i) = (0.5d0/pi_l**2)*p2_bins(i)*estert(i)*             &
     &                ( dydt(i+3) + dydt(i+3+n_p_bins_l) )
      end do
      !$OMP END PARALLEL DO

      ! Time, efold and Leptasym derivatives w.r.t time
      !------------------------------------------------------------------
      ! Net sterile neutrino energy density
      rhoster = NIntegrate(p_bins_l, rhosteri, n_p_bins_l,              &
     &                     p_bins_l(1), p_bins_l(n_p_bins_l))
      ! Energy and Entropy densities in SM species
      rhoSM = (pi_l**2/30.0d0)*gstar(temp,1,0)*temp**4
      entropySM = (2.0d0*pi_l**2/45.0d0)*gstar(temp,2,0)*temp**3

      ! Hubble rate, in s^-1
      hubcst = sqrt(8.0D0*pi_l/3.0D0)*(1.0D0/(mp_l*hbar_l))*            &
     &         sqrt(rhoSM + rhoster)

      dydt(1) = 1.0D0 ! dtdt
      dydt(2) = hubcst ! d/dt(\int H dt)
      dydt(3) = NIntegrate(p_bins_l, dLdti, n_p_bins_l,                 &
     &                     p_bins_l(1), p_bins_l(n_p_bins_l))

      ! drho_{sterile neutrino}/dt
      drhosdt = NIntegrate(p_bins_l, drhosdti, n_p_bins_l,              &
     &                     p_bins_l(1), p_bins_l(n_p_bins_l))

      ! drho_(standard model)/dT
      drhoSMdT = (pi_l**2/30.0)*( 4.0d0*gstar(temp,1,0)*temp**3 +       &
     &                            gstar(temp,1,1)*temp**4 )

      ! dt/dT = -(d\rho_SM/dT)/(3 H (rho_SM + p_SM ) + drho_s/dt) 
      ! extra factor to account for sterile production
      dtdT = -drhoSMdT/( drhosdt + 3.0d0*hubcst*entropySM*temp )

      ! Convert time derivatives to temperature derivatives
      !------------------------------------------------------------------
      Do i=1, nvr
         dydtemp(i) = dtdT*dydt(i)
      End do
      ! Correction to dL/dT since T doesn't go as a^{-1}
      dydtemp(3) = dydtemp(3) -                                         &
     &             3.0D0*( hubcst*dtdT + (1.0D0/temp) )*leptasym
      !------------------------------------------------------------------
      RETURN
      END Subroutine ! derivstemp!}}}
!---------------------------------------------------------------------
