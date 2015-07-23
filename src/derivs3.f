!******************************************************************
      SUBROUTINE derivstemp(nvr,x,y,dydx) !RK derivative
!******************************************************************!{{{
!     version 5.0 - added a more realistic set of scattering rates !{{{
!     through the quark hadron transition.     

!     version 4.0 - added effect of asymmetry redistribution,
!     corrections for proper gstar, backreaction of production on
!     temp evolution, and correct evolution of L

!     version 3.0 - corrected for the different evolution of 
!     scale-factor and T with time, for varying no. of relativistic
!     d.o.f

!     version 2.0 - changed dL/dt to direct difference of sterile
!     neutrino distributions to be applicable with large lepton numbers

!     version 1.7 (version after 1.4) - with a linear scaling of the
!     number of scatterers in the softening phase of the QH transition
!     -- also includes some optimizations

!     version 1.3 - with an actual iterated form of the actual abundance
!     of scatterers, and an iteration for the thermal potential (vtch)
!     for the mu and tau lepton abundance evolution!}}}

      Use Phconst!{{{
      Use Varpr
      Use Pathvar
      Use FermiDirac
      Use Suscep
      Use Rates
      Use DOF
      Implicit None

      Integer nvr
      Double precision x, y(nvr), dydx(nvr) ! RK assignment
      
!---- LOCAL VARIABLES
      ! Parameter definitions
      Double precision temp, scalef, leptasym
      Double precision mx, mnux, munux ! Flavor dependent params
      Double precision ex, enux
      Double precision vl, vtf, sctf(NZS) ! potentials and scaled \Gamma

      ! Local versions of variables on the heap
      Double precision pil, Gfl, hbarl
      Double precision msl, dm2l, s2l, cl
      Integer NZSl

      ! Sterile nu integrands
      Double precision p(NZS), yy2(NZS), estert(NZS) 
      Double precision rhosteri(NZS), rhoster
      ! Production of sterile neutrinos      
      Double precision gams(NZS), gambs(NZS)! \Gamma in matter
      Double precision fnux(NZS), fnuxbar(NZS) ! SM neutrino PSDs
      Double precision dydt(nvr) ! Time derivatives
      Double precision dLdti(NZS), drhosdti(NZS) ! Integrands for tder
      ! Standard model integrands
      Double precision rhoSM, entropySM 
      ! Time and asymmetry evolution
      Double precision hubcst, drhosdt, drhoSMdT, dtdT

!---  V_asymmetry, scattering rate factor, Integral of a spline
      External vasym, NIntegrate ! gamfac
      Double precision vasym, NIntegrate ! gamfac

      Integer i!}}}

      If (nvr.lt.(2*NZS+3)) Then
              Stop "derivs called with too few variables."
      End if

C---  Parameter definitions -- populating the stack
!------------------------------------------------------------------
      temp = x
      scalef = exp(-y(2)) ! (a_start/a_current)
      leptasym = y(3) ! Leptasym in units of T^3

!---  Flavor dependent parameters
      mx = 0.0D0
      mnux = 0.0D0
      munux = 0.0D0
      If (flavor(1:1).eq.'e') Then
              mx = me
              mnux = mnue
              munux = dmudlx(temp,4,flavor)*leptasym ! mu/T for \nu_X
      Else if (flavor(1:1).eq.'m') Then
              mx = mmu
              mnux = mnumu
              munux = dmudlx(temp,5,flavor)*leptasym ! mu/T for \nu_X
      Else if (flavor(1:1).eq.'t') Then
              mx = mtau
              mnux = mnutau
              munux = dmudlx(temp,6,flavor)*leptasym ! mu/T for \nu_X
      Else
              Print *,flavor
              Stop "Flavor needs to be one of e, m or t"
      End if
      ex = 4.0d0*EF(mx/temp) ! Net energy of charged leptons
      enux = 2.0d0*EF(mnux/temp) ! Net energy of neutral leptons

      ! Asymmetry potential in MeV
      vl = vasym(temp, flavor)*Gf*leptasym*temp**3
      ! Dimensionless factor relating thermal potential and p
      vtf = -(8.0*sqrt(2.0D0)/3.0)*Gf*temp**4*( enux/mz**2 + ex/mw**2 )

      ! read into stack
      pil = pi
      Gfl = Gf
      hbarl = hbar
      msl = ms
      dm2l = dm2
      s2l = s2
      cl = c
      NZSl = NZS

!---  Define Lagrangian momentum bins on stack, and update heap 
      Do i=1, NZS
         p(i) = scalef*yyinit(i)
         yy(i) = p(i) 
      End do

!---  Dimensionless factors relating scattering rates to G_F^2 p T^4
      Call GamFac(NZS, p, temp, flavor, sctf)
      ! Dimensionless factors relating the scattering rate to p 
      Do i=1, NZS
         sctf(i) = sctf(i)*(Gf**2)*temp**4
      End do

!------------------------------------------------------------------

C---  Calculating df/dt entirely on the stack
!------------------------------------------------------------------
      !$OMP PARALLEL DO
      Do i=1, NZSl
         yy2(i) = p(i)**2 
         estert(i) = sqrt(yy2(i) + msl**2)
         ! Integrand for \rho_sterile
         rhosteri(i) = (0.5d0/pil**2)*yy2(i)*estert(i)*
     *                 ( y(i+3) + y(i + 3 + NZSl) )

         ! Sterile production, code to evaluate d/dt(f_\nu, f_\bar{\nu})
         ! SM neutrino phase space densities
         fnux(i) = 1.0d0/
     *             ( 1.0d0 + exp((sqrt(yy2(i) + mnux**2)/temp) - munux))
         fnuxbar(i) = 1.0d0/
     *             ( 1.0d0 + exp((sqrt(yy2(i) + mnux**2)/temp) + munux))

         ! Interaction rates in matter, in MeV
         gams(i) = 0.25d0*sctf(i)*p(i)*s2l/
     *             ( s2l + (sctf(i)*yy2(i)/dm2l)**2 + 
     *             (cl - vtf*2.0d0*yy2(i)/dm2l - vl*2.0d0*p(i)/dm2l)**2)
         gambs(i) = 0.25d0*sctf(i)*p(i)*s2l/
     *             ( s2l + (sctf(i)*yy2(i)/dm2l)**2 + 
     *             (cl - vtf*2.0d0*yy2(i)/dm2l + vl*2.0d0*p(i)/dm2l)**2)

         ! Sterile PSD derivatives w.r.t time
         dydt(i+3) = (1.0d0/hbarl)*gams(i)*( fnux(i) - y(i+3) )
         dydt(i+3+NZSl) = (1.0d0/hbarl)*gambs(i)*
     *                    ( fnuxbar(i) - y(i+3+NZSl) )

         ! Integrands for leptasym and sterile energy rates of change
         ! solely due to sterile neutrino production
         dLdti(i) = -(0.5d0/(pil**2*temp**3))*yy2(i)*
     *               ( dydt(i+3) - dydt(i+3+NZSl) )
         drhosdti(i) = (0.5d0/pil**2)*yy2(i)*estert(i)*
     *                 ( dydt(i+3) + dydt(i+3+NZSl) )
      End do
      !$OMP END PARALLEL DO
!------------------------------------------------------------------

C---  Time, efold and Leptasym derivatives w.r.t time
!------------------------------------------------------------------
      ! Net sterile neutrino energy density
      rhoster = NIntegrate(p, rhosteri, NZSl, p(1), p(NZSl))
      ! Energy and Entropy densities in SM species
      rhoSM = (pil**2/30.0d0)*gstar(temp,1,0)*temp**4
      entropySM = (2.0d0*pil**2/45.0d0)*gstar(temp,2,0)*temp**3

      ! Hubble rate, in s^-1
      hubcst = sqrt(8.0D0*pi/3.0D0)*(1.0D0/(mp*hbar))
     *         *sqrt(rhoSM + rhoster)
      !hubcst = 0.360112d0*sqrt(rhoSM + rhoster)

      dydt(1) = 1.0D0 ! dtdt
      dydt(2) = hubcst ! d/dt(\int H dt)
      dydt(3) = NIntegrate(p, dLdti, NZSl, p(1), p(NZSl))

      ! drho_{sterile neutrino}/dt
      drhosdt = NIntegrate(p, drhosdti, NZSl, p(1), p(NZSl))
      ! drho_(standard model)/dT
      drhoSMdT = (pil**2/30.0)*( 4.0d0*gstar(temp,1,0)*temp**3 + 
     *                        gstar(temp,1,1)*temp**4 )

      ! dt/dT = -(d\rho_SM/dT)/(3 H (rho_SM + p_SM ) + drho_s/dt) 
      ! extra factor to account for sterile production
      dtdT = -drhoSMdT/( drhosdt + 3.0d0*hubcst*entropySM*temp )
!------------------------------------------------------------------

C--- Convert time derivatives to temperature derivatives
!------------------------------------------------------------------
      Do i=1, nvr
         dydx(i) = dtdT*dydt(i)
      End do
      ! Correction to dL/dT since T doesn't go as a^{-1}
      dydx(3) = dydx(3) - 3.0D0*( hubcst*dtdT + (1.0D0/temp) )*leptasym
!------------------------------------------------------------------
      RETURN
!******************************************************************!}}}
      END Subroutine ! derivstemp
!******************************************************************
