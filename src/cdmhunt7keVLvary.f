!******************************************************************
      program cdmqcd
!******************************************************************!{{{
C Main program to run and test the code 

      Use Varpr!{{{
      Use Odeintmemory
      Use Rates
      Implicit None

      ! Root finder routine, and \Omega_s h^2(L_i)
      External zeroin, evalstd
      Double precision zeroin, evalstd

      ! Local variables
      Double precision Llow, Lhigh, tol, leptasymi
      Integer i, j, k!}}}

      ! Limits for initial lepton asymmetry in units of T^3.
      Llow = 1.000D-3
      Lhigh = 6.000D-3
      tol = 5.0D-6 ! Tolerance to which we want L_i
      ! Note Kev uses L/(2 \zeta(3)/pi^2) T^3. 
      ! Set background params that don't depend on L_i
      Call start

      !Print *,"nlTrate=",nlTrate
      !Print *,"nlTratesub=",nlTratesub
      !Print *,"nlpTrate=",nlpTrate
      !Call writearray(lTratetsub,lratetsub,nlTratesub,1,0,
      !*                (nlpTrate-1)*nlTratesub,'',6)

      Do i=0, 100, 10 ! index for \sin^2\theta
         Do j=0, 100, 10 ! index for m_s
            Print *,"i=",i,"j=",j
            ! Initialize sin^2(2\theta) and m_s
            Call assignpm(i, j)
            ! Find value of l that satisfies \omega_s h^2 = planck
            leptasymi = zeroin(Llow, Lhigh, evalstd, tol)
            !Print *,"leptasymi=",leptasymi

            ! Now output initial and final states
            !Call outputode(1)
            !Call outputode(kount)
            Do k=1, kount
               Call outputode(k)
            End do
         End do
      End do

      Stop
!******************************************************************!}}}
      END ! PROGRAM cdmqcd
!******************************************************************

!******************************************************************
      SUBROUTINE start
!******************************************************************!{{{
C Loads external files and initializes some modules that don't change 
C during execution. These set run parameters that don't depend on L_i

      Use Phconst!{{{
      Use FermiDirac
      Use Suscep
      Use Rates
      Use DOF
      Use Varpr
      Use Pathvar
      Implicit none

!---- Local Variables
      Double precision g1, g2!}}}

!---- Initialize all the tables and splines
      Call InitFD
      Call InitSuscep
      Call InitRates
      Call InitDOF
      ! Now the knots and coeffs are stored for future use by pchi, EF
      ! dmudlx, scps, GamFac and gstar

!---- Range of Sterile parameters being explored
      ! Masses
      minms = 6.75e-3  ! MeV
      maxms = 7.5e-3   ! MeV
      ! sin^2 \theta
      ! 2 sigma range for stat+sys quadratic combination
      minsin = 0.8e-11
      maxsin = 20.e-11
      ! Step sizes
      ddm2 = (dlog10(maxms)-dlog10(minms))/dble(ni2) 
      ds  = (dlog10(maxsin)-dlog10(minsin))/dble(ni1)

      !! Initialize ms etc. for testing purposes
      flavor = 'm'  ! Flavor of active neutrino that mixes with the
                    ! steriles     

!---- Initial lepton asymmetries and background parameters
      tmevi = 1.0D+4 ! <t> Start at lower values below M_Z </t>
      tmevf = 1.0D+1 ! <t> End above weak decoupling </t>

      g1=gstar(tmevi,1,0)
      g2=gstar(tmevf,1,0)

      ! Start and stop times, w/ approximate relation to temperatures
      ! as t \sim 1/2H. Exact ones not needed, as we only measure offset
      ti = sqrt(45.0d0/(16.0d0*pi**3))*mp*hbar/(sqrt(g1)*(tmevi)**2) 
      tf = sqrt(45.0d0/(16.0d0*pi**3))*mp*hbar/(sqrt(g2)*(tmevf)**2)

      ! Initialize bins in the calling function

      Return
!******************************************************************!}}}
      END !subroutine start
!******************************************************************

!******************************************************************
      Double precision Function evalstd(leptasymi) 
!******************************************************************!{{{
C Integrates from initial temperature down to final temperature. Calls 
C the ODE integrator, and then calls function to compute integrals over 
C the sterile distribution in order to compute \Omega_s,bar{s} h^2. 
C Finds \Omega_net h^2 - (planck value)      

      Use Pathvar!{{{
      Use Varpr
      Use Odeintmemory
      Use Phconst
      Implicit None
      Double precision leptasymi

      ! local variables
      Double precision dtiny ! dtiny ensures spline range
      Parameter (dtiny=1.0D-5)
      Double precision pmin, pmax 
      Double precision omegas, omegasbar, avgp, nnu, freest
      Double precision eps
      External derivstemp, rkqs!}}}

!---- The sterile distribution
      ! Lowest and highest initial momentum bins, in units of T
      pmin = 1.0D-3 
      pmax = 20.0d0 - dtiny
      ! Initialize the time, efolds, leptasym and sterile PSD
      Call Inity(ti, leptasymi)
      ! Initialize the sterile momentum bins
      Call Inityy(tmevi, pmin, pmax)
      ! Acceptable relative error
      eps = 1.0e-6      ! Acceptable relative error
      ! Call integrator      
      !CALL odepackint(yinit,nvar,tmevi,tmevf,eps,.false.,derivstemp)
      CALL odepackint(yinit,nvar,tmevi,tmevf,eps,.true.,derivstemp)
      ! Compute integrals 
      Call sterint(kount, omegas, omegasbar, nnu, avgp, freest)

      evalstd = omegas + omegasbar - omegacont
      
      Return
!******************************************************************!}}}
      End Function !evalstd
!******************************************************************

!******************************************************************
      SUBROUTINE sterint(ind, omegas, omegasbar, nnu, avgp, freest)
!******************************************************************!{{{
C Computes integrals over the sterile nu phase-space distribution from
C the i^th stored state of the ODE integrator. In order to map to z=0, 
C it assumes that it is called above the epoch of weak decoupling i.e. 
C above ~1 MeV, but below the quark hadron transition, so that the 
C neutrino temperature ~1/a. Stores \omega_s h^2 in omegas, 
C \omega_sbar h^2 in omegasbar, g_{*,eff} for one species in nnu (i.e. 
C for both \nu_s and \bar{\nu_s} together, g_{*,eff} = 2*nnu) <p/T_nu> 
C in avgp and free-streaming length (in Mpc) in freest

      Use Odeintmemory!{{{
      Use Phconst
      Use Params
      Use Varpr
      Implicit none

      Integer ind
      Double precision omegas,omegasbar,nnu,avgp,freest

!---  Local variables
      Double precision hub02, temp
      Parameter (hub02=1.d0) ! h^2. Set to 1 for calculating omegah^2
      Double precision yyl(NZS), yy2l(NZS)
      Double precision numi(NZS), rhois(NZS), rhoisbar(NZS)
      Double precision weightedpi(NZS)
      Double precision num, rhos, rhosbar, weightedp
      Double precision numf2, rhob2
      External NIntegrate
      Double precision NIntegrate!}}}

      If ((ind.le.0).or.(ind.gt.kount)) Then
              Stop "sterint called with invalid index."
      End if

      temp = xp(ind)

      !yyl = yyp(1:NZS,ind) ! Local copy of momentum bins
      !yy2l = yyl**2 ! p^2
      !! f_\nu + f_\bar{\nu}
      !psd = yp(4:NZS+3,ind) + yp(NZS+4:2*NZS+3,ind)

      !! Integrands for number density, energy and average momentum
      !numi = (0.5d0/pi**2)*yy2l*psd
      !rhoi = (0.5d0/pi**2)*yy2l*sqrt(yy2l + ms**2)*psd
      !weightedpi = (0.5d0/pi**2)*yyl**3*psd

      !num = NIntegrate(yyl, numi, NZS, yyl(1), yyl(NZS))
      !rho = NIntegrate(yyl, rhoi, NZS, yyl(1), yyl(NZS))
      !weightedp = NIntegrate(yyl, weightedpi, NZS, yyl(1), yyl(NZS))

      !! const = (3\zeta(3)/(2\pi^2)), for massless \nu, \nubar with g=1
      !numf2 = 0.182684d0*temp**3
      !! const = 2 \pi^2/30, for a massless boson + anti-boson
      !rhob2 = 0.657974d0*temp**4

      !! omegas = ms*112.31d0*(int/numf2)/1.054d-2
      !! 112.31 is the number density of neutrinos today in cm^{-3}
      !! 1.054d-2 is rho_crit/h^2 in units of MeV cm^{-3}      
      !omegas = ms*1.06556d4*(num/numf2)/hub02
      !nnu = rho/rhob2
      !avgp = weightedp/(num*temp)
      !! freest = 40.d0*3.d-5/ms*(avgp/3.151)  ! Mpc
      !freest = 3.80832d-4*avgp/ms  ! Mpc

      ! physical momenta today, p_today in MeV = (p/T_nu)*T_nu_today
      yyl = (Tcmb*nfactor/temp)*yyp(1:NZS,ind)
      yy2l = yyl**2 ! p_today^2

      ! Integrands for number density, energy and average momentum
      ! today in units of MeV^3, MeV^4, MeV^4, if production were to
      ! completely stop below this temperature
      numi = (0.5d0/pi**2)*yy2l*
     *       (yp(4:NZS+3,ind) + yp(NZS+4:2*NZS+3,ind))
      rhois = (0.5d0/pi**2)*yy2l*sqrt(yy2l + ms**2)*yp(4:NZS+3,ind)
      rhoisbar = (0.5d0/pi**2)*yy2l*sqrt(yy2l + ms**2)*
     *           yp(NZS+4:2*NZS+3,ind)
      weightedpi = (0.5d0/pi**2)*yyl**3*
     *             (yp(4:NZS+3,ind) + yp(NZS+4:2*NZS+3,ind))

      num = NIntegrate(yyl, numi, NZS, yyl(1), yyl(NZS))
      rhos = NIntegrate(yyl, rhois, NZS, yyl(1), yyl(NZS))
      rhosbar = NIntegrate(yyl, rhoisbar, NZS, yyl(1), yyl(NZS))
      weightedp = NIntegrate(yyl, weightedpi, NZS, yyl(1), yyl(NZS))

      omegas = rhos/rhocrit ! \Omega_s h^2 = rho_s/(\rho_c/h^2)
      omegasbar = rhosbar/rhocrit

      ! const = 2 \pi^2/30, for a massless boson + anti-boson
      rhob2 = 0.657974d0*(Tcmb*nfactor)**4
      nnu = (rhos+rhosbar)/rhob2
      avgp = weightedp/(num*Tcmb*nfactor)
      ! freest = 40.d0*3.d-5/ms*(avgp/3.151)  ! Mpc ! <t> Remove </t>
      freest = 3.80832d-4*avgp/ms  ! Mpc

      Return
!******************************************************************!}}}
      END !Subroutine sterint
!******************************************************************
    
!******************************************************************
      Double precision FUNCTION vasym(T, X)
!******************************************************************!{{{
C Evaluates dimensionless factor relating the asymmetry potential
C V_A to G_F*(L_X in T^3), at temperature T (in MeV). The net baryon 
C number is assumed to equal zero.

      Use Phconst!{{{
      Use FermiDirac
      Use Suscep
      Implicit None

      Double precision T
      Character*8 X

      ! Internal variables
      Double precision vlc, vln, vlq ! contributions of charged,
                                     ! neutral and strong sectors!}}}

      ! First the parts which don't depend on flavor
      vlc = sqrt(2.0D0)*(-0.5D0 + 2.0D0*s2w)*(
     *           2.0D0*pchi(me/T)*dmudlx(T,1,X) +
     *           2.0D0*pchi(mmu/T)*dmudlx(T,2,X) +
     *           2.0D0*pchi(mtau/T)*dmudlx(T,3,X) )

      vln = sqrt(2.0D0)*(
     *           pchi(mnue/T)*dmudlx(T,4,X) +
     *           pchi(mnumu/T)*dmudlx(T,5,X) +
     *           pchi(mnutau/T)*dmudlx(T,6,X) )

      vlq = sqrt(2.0D0)*(1.0D0 - 2.0D0*s2w)*(
     *           scps(T,1)*dmudlx(T,7,X) +
     *           scps(T,2)*dmudlx(T,8,X))

      ! Then parts which depend on flavor

      If (X(1:1).eq.'e') Then
             vlc = vlc + sqrt(2.0D0)*2.0D0*pchi(me/T)*dmudlx(T,1,X)
             vln = vln + sqrt(2.0D0)*pchi(mnue/T)*dmudlx(T,4,X)
      Else if (X(1:1).eq.'m') Then
             vlc = vlc + sqrt(2.0D0)*2.0D0*pchi(mmu/T)*dmudlx(T,2,X)
             vln = vln + sqrt(2.0D0)*pchi(mnumu/T)*dmudlx(T,5,X)
      Else if (X(1:1).eq.'t') Then
             vlc = vlc + sqrt(2.0D0)*2.0D0*pchi(mtau/T)*dmudlx(T,3,X)
             vln = vln + sqrt(2.0D0)*pchi(mnutau/T)*dmudlx(T,6,X)
      Else
             Print *,X
             Stop "Function vl takes X = e, m or t."
      End If
      vasym = vlc + vln + vlq

      RETURN
!******************************************************************!}}}
      END FUNCTION !vasym
!******************************************************************

!******************************************************************
      SUBROUTINE writearray(x, y, m, n, xstride, ystride, header, fp)
!******************************************************************!{{{
C Subroutine to write to a file indexed by fp the values 
C x(xstride+1:xstride+m) (abscissae), and y(ystride+1:ystride+m*n) 
C (n columns of ordinates) in double precision. Works for n-d arrays 
C stored in column major. This form makes it easy to output everything. 
C Typecast inputs to this function, if so needed!
      
      Integer m, n, xstride, ystride, fp!{{{
      Double precision x(xstride+m), y(ystride + m*n)
      Character*(*) header
      
      Integer i, j!}}}

      Write(fp,*) '!',header

      i=1
      j=1

      ! Gymnastics to avoid spurious newline at end.
      ! First write upto penultimate row
      Do while (i.le.m-1)
         Write(fp,100) x(xstride+i)
         ! Write upto penultimate element, array stored in column major.
         Do while (j.le.n-1)
            Write(fp,100) y(ystride + (j-1)*m + i)
            j = j + 1
         End do
         ! Write the last element
         Write(fp,110) y(ystride + (n-1)*m + i)
         j = 1
         i = i + 1
      End do
      ! Now write the last row
      Write(fp,100) x(xstride+m)
      Do j=1, n
         Write(fp,100) y(ystride + j*m)
      End do

 100  Format(1pE24.16, $)
 110  Format(1pE24.16)

      Return
!******************************************************************!}}}
      End Subroutine ! writearray
!******************************************************************

!******************************************************************
      SUBROUTINE outputode(ind)
!******************************************************************!{{{
C Subroutine to (re)create and populate an output directory with the 
C state of the ode integrator. (Over)writes a param file, and 
C a file with p, f_\nu_s, f_\nu_sbar. Input the required index.
C <t> For now, works only on a unix env </t>
C <pt> Modified to output p/T(min) and p/T(max) </pt>

      Use Odeintmemory!{{{
      Use Phconst
      Use Params
      Use Pathvar 
      Use Varpr
      Implicit None

      Integer ind ! Index in list of saved states

      ! Local variables
      Double precision temp, lkev
      Double precision omeganet, omegas, omegasbar, nnu, avgp, freest
      Double precision pTmin, pTmax
      Logical state ! To check if state file already exists
      Integer fppar, fpint, fpstate ! File units
      Character*9 dmake
      Character*64 dname, parname, fname, statename!}}}

      If (kount.eq.0) Then
              Stop "outputode called without any saved values."
      End if

      If ((ind.le.0).or.(ind.gt.kount)) Then
              Stop "outputode called with invalid index."
      End if

      temp = xp(ind) ! Temperature
      ! Evaluate (\Omega_s, sbar) h^2 
      call sterint(ind, omegas, omegasbar, nnu, avgp, freest)
      omeganet = omegas + omegasbar

      ! First (re)create output directory
      dmake = 'mkdir -p '
      Write(dname, "(A9,A2,1pE9.3E2,A2,1pE9.3E2,A1,1pE9.3E2)")
     *      'outfiles/','ms',ms,'s2',s2,'L',li
      Call system(dmake // trim(dname))

      ! Then (over)write a param file
      ! Initial lepton asymm in Kev's units
      lkev = li*pi**2/(2.0D0*zeta3)
      parname = 'params.dat'
      fppar = 25
      Open(unit=fppar,file=trim(dname)//'/'//trim(parname),
     *     status='unknown')
      Write(fppar,"(1pE24.16,A6)") ms, ' # m_s'
      Write(fppar,"(1pE24.16,A16)") s2, ' # \sin^2 \theta'
      Write(fppar,"(1pE24.16,A12)") lkev, ' # L/n_gamma'
      Write(fppar,"(1pE24.16,A17)") omeganet, ' # \Omega_wdm h^2'
      Write(fppar,"(1pE24.16,A15)") omegas, ' # \Omega_s h^2'
      Write(fppar,"(1pE24.16,A18)") omegasbar, ' # \Omega_sbar h^2'
      Write(fppar,"(1pE24.16,A8)") avgp, ' # <p/T>'
      Close(unit=fppar)

      ! Then output the PSDs
      Write(fname, "(A8,1I3.3,A4)") 'Snapshot',ind,'.dat'
      fpint = 26
      Open(unit=fpint,file=trim(dname)//'/'//trim(fname),
     *     status='unknown')
      Call writearray(yyp, yp, NZS, 2, (ind-1)*NMAX, (ind-1)*NMAX+3,
     *     'p(MeV) \t \delta f_\\nu \t \delta f_\\nubar', fpint)
      Close(unit=fpint)

      ! Then append/output auxiliary run variables into a state file
      statename = 'state.dat'
      fpstate = 27
      ! Current Lepton asymm in Kev's units
      lkev = yp(3,ind)*pi**2/(2.0D0*zeta3) 
      ! Minimum and maximum p/T
      pTmin = yyp(1,ind)/temp
      pTmax = yyp(NZS,ind)/temp
      ! Check if file already exists, as we are appending to it
      Inquire(file=trim(dname)//"/"//trim(statename), exist=state)
      If (state) Then
              Open(unit=fpstate,file=trim(dname)//"/"//trim(statename),
     *             status='old',position='append',action="write")
      Else
              Open(unit=fpstate,file=trim(dname)//"/"//trim(statename),
     *             status="new",action="write")
              ! Write header
              Write(fpstate,"(A1, A23,$)") '!','T (MeV)'
              Write(fpstate,"(6A24)") 't (sec)', 'ln(a/a_i)',
     *                                'L/n_gamma', '\Omega_wdm h^2',
     *                                'p/T(min)', 'p/T(max)'         
      End if
      Write(fpstate,"(1p 7E24.16)") temp, yp(1,ind), yp(2,ind), 
     *                              lkev, omeganet, pTmin, pTmax
      Close(unit=fpstate)

      ! Also write to stdout
      Print *,'T=',temp,'t=',yp(1,ind),'ln(a/a_i)=',yp(2,ind),
     *        'L/T^3=',yp(3,ind),'Omega_wdm h^2=',omeganet,
     *        'p/T(min)=',pTmin,'p/T(max)=',pTmax 

      Return
!******************************************************************!}}}
      End Subroutine ! outputode
!******************************************************************

!******************************************************************
      SUBROUTINE odeint(ystart,nvars,x1,x2,eps,h1,hmin,nok,nbad,outflag
     $     ,derivs,rkqs)
!******************************************************************!{{{
C Driver subroutine for ode integration. It calls subroutine outputode
C to save the state at some intermediate times      
C <temp> Changed to output error calculation </temp>      

      Use Odeintmemory!{{{
      Implicit None
      Integer nvars, nok, nbad
      Double precision ystart(nvars), x1, x2, eps, h1, hmin
      LOGICAL outflag
      External derivs, rkqs

      ! Local
      Double precision h, hdid, hnext, x, xsav,
     $                 y(nvars), dydx(nvars), yscal(nvars)
      Integer i,nstp!}}}

      If (nvars.gt.NMAX) Then
              Print *,nvars,NMAX
              Stop "Not enough memory allocated to store internal
     *              variables of the integrator. Increase NMAX."
      End if

      x    = x1
      h    = sign(h1,x2-x1)
      nok  = 0
      nbad = 0
      kount= 0
      xsav = x

      if(outflag)then 
         kmax = 100 ! Output at 40 log spaced temperatures.
         dxsav = abs(log(x2)-log(x1))/kmax
      else 
         kmax = 0
      endif
      do 11 i=1,nvars
         y(i)=ystart(i)
 11   continue
      if (kmax.gt.0) xsav=x*exp(-2.D0*dxsav) ! Assures storage of first step
      print*, MAXSTP
      do 16 nstp=1,MAXSTP
         if (MOD(nstp,10000) .EQ. 0) print*, nstp
         CALL derivs(nvars,x,y,dydx)
         do 12 i=1,nvars
            yscal(i)=dabs(y(i))+dabs(h*dydx(i))+fTINY
 12      continue
         if(kmax.gt.0)then
            if(abs(log(x)-log(xsav)).gt.abs(dxsav)) then
               if(kount.lt.kmax-1)then
                  kount=kount+1
                  xp(kount)=x
                  do 13 i=1,nvars
                     yp(i,kount)=y(i)
 13               continue
                  xsav=x
                  ! Output the state
                  Call outputode(kount)
               endif
            endif
         endif
         if((x+h-x2)*(x+h-x1).gt.0.D0) h=x2-x ! To prevent overshooting.
         Call rkqs(y,dydx,nvars,x,h,eps,yscal,hdid,hnext,derivs)
         if(hdid.eq.h)then
            nok=nok+1
         else
            nbad=nbad+1
         endif

         if((x-x2)*(x2-x1).ge.0.D0)then !if this was the last step
!            do 14 i=1,nvars
!               ystart(i)=y(i) ! Write back into input array
! 14         continue
            if(kmax.ne.0)then ! Save final step
               kount=kount+1
               xp(kount)=x
               do 15 i=1,nvars
                  yp(i,kount)=y(i)
 15            continue
               ! Output the state
               Call outputode(kount)
            endif
            Print *,'n_ok =',nok,'n_bad =',nbad ! nsteps
            Return ! Break out and end odeint
         endif
         if(dabs(hnext).lt.hmin) stop
     *        'stepsize smaller than minimum in odeint'
         h=hnext
 16   continue
      stop 'too many steps in odeint'

      Return
!******************************************************************!}}}
      END ! SUBROUTINE odeint
!     (C) Copr. 1986-92 Numerical Recipes Software V,3.
!******************************************************************

!******************************************************************
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
!******************************************************************!{{{
! Adaptive stepsize control for Runge Kutta stepper

      Use Odeintmemory!{{{
      Implicit None
      Integer n
      Double precision y(n), dydx(n), x, htry, 
     *                 eps, yscal(n), hdid, hnext
      External derivs
!U    USES derivs,rkck

      ! Local
      Double precision h, ytemp(NMAX), yerr(NMAX)
      Double precision errmax,xnew

      Double precision SAFETY, PGROW, PSHRNK, ERRCON
      Parameter (SAFETY=0.9D0,PGROW=-.2D0,PSHRNK=-.25D0,ERRCON=1.89D-4)
      Integer i!}}}

      If (n.gt.NMAX) Then
              Print *,n,NMAX
              Stop "Not enough memory allocated to store internal
     *              variables of the integrator. Increase NMAX."
      End if

      h=htry

 1    Call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)

      errmax=0.D0
      Do 11 i=1,n
         errmax=max(errmax,dabs(yerr(i)/yscal(i)))
 11   continue
      errmax=errmax/eps
      if(errmax.gt.1.D0)then
          h=SAFETY*h*(errmax**PSHRNK)
          if(h.lt.0.1D0*h)then
            h=.1D0*h
          endif
          xnew=x+h
          if(xnew.eq.x)stop 'stepsize underflow in rkqs'
          goto 1
      else
          if(errmax.gt.ERRCON)then
            hnext=SAFETY*h*(errmax**PGROW)
          else
            hnext=5.D0*h
          endif
          hdid=h
          x=x+h
          do 12 i=1,n
            y(i)=ytemp(i)
 12       continue
          return
      endif
 100  FORMAT('L Sign Change',//,4(1pd13.5))
 130  FORMAT(2(1pe13.5))
 150  FORMAT(i8)
!******************************************************************!}}}
      END ! SUBROUTINE rkqs
!     (C) Copr. 1986-92 Numerical Recipes Software V,3.
!******************************************************************

!******************************************************************
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
!******************************************************************!{{{
C Dumb Runge Kutta stepper routine

      Use Odeintmemory!{{{
      Implicit None

      Integer n
      Double precision y(n), dydx(n), x, h, yout(n), yerr(n)
      External derivs

!U    USES derivs
      ! Local
      Integer i
      Double precision ak2(NMAX),ak3(NMAX),ak4(NMAX),
     $ak5(NMAX),ak6(NMAX),
     *ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      Parameter (A2=.2D0,A3=.3D0,A4=.6D0,A5=1.D0,A6=.875D0,B21=.2D0,
     $B31=3.D0/40.D0,
     *B32=9.D0/40.D0,B41=.3D0,B42=-.9D0,B43=1.2D0,B51=-11.D0/54.D0,
     $B52=2.5D0,
     *B53=-70.D0/27.D0,B54=35.D0/27.D0,B61=1631.D0/55296.D0,
     $B62=175.D0/512.D0,
     *B63=575.D0/13824.D0,B64=44275.D0/110592.D0,B65=253.D0/4096.D0,
     $C1=37.D0/378.D0,
     *C3=250.D0/621.D0,C4=125.D0/594.D0,C6=512.D0/1771.D0,
     $DC1=C1-2825.D0/27648.D0,
     *DC3=C3-18575.D0/48384.D0,DC4=C4-13525.D0/55296.D0,
     $DC5=-277.D0/14336.D0,
     *DC6=C6-.25D0)!}}}
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
 11   continue
      Call derivs(n,x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
 12   continue
      Call derivs(n,x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
 13   continue
      Call derivs(n,x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
 14   continue
      Call derivs(n,x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
 15   continue
      Call derivs(n,x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
 16   continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
 17   continue
      return
!******************************************************************!}}}
      END ! SUBROUTINE rkck
!     (C) Copr. 1986-92 Numerical Recipes Software V,3.
!******************************************************************

!******************************************************************
      SUBROUTINE odepackint(ystart,nvars,x1,x2,eps,outflag,derivs)
!******************************************************************!{{{
C Driver subroutine for ode integration using LSODE. It also calls
C subroutine outputode to save the state at intermediate times. 

      Use Odeintmemory!{{{
      Use Params
      Use Pathvar ! use pathvar to save momenta at each step
      Implicit None
      Integer nvars
      Double precision ystart(nvars), x1, x2, eps
      LOGICAL outflag
      External derivs

      ! Local
      Double precision y(nvars), T, TOUT, RTOL, ATOL 
      Integer ITOL, ITASK, ISTATE, IOPT, 
     *        LRW, LIW, MF
      !! Use full Jacobian
      !Integer IWORK(20+nvars)
      !Double precision RWORK(22+nvars*MAX(16,nvars+9))
      ! Use non-stiff method
      Integer IWORK(20)
      Double precision RWORK(20+nvars*16)
      Integer i, nstp
      ! Should actually be external/interface for the Jacobian, but it 
      ! fails during link-time unless I use this. Mysterious.
      Double precision JDUM!}}}

      If (nvars.gt.NMAX) Then
              Print *,nvars,NMAX
              Stop "Not enough memory allocated to store internal
     *              variables of the integrator. Increase NMAX."
      End if

      ! LSODE compulsory parameters and flags
      ITOL = 1          ! Using scalar absolute tolerance
      RTOL = eps        ! Relative tolerance
      ATOL = 1.0D-30    ! Absolute tolerance = 0, set to small for f=0
      ITASK = 4         ! Don't overshoot
      ISTATE = 1        ! First call, so initialize
      IOPT = 1          ! We will be using optional inputs
      !LRW = 22+nvars*MAX(16,nvars+9) ! Length of work array
      !LIW = 20+nvars    ! Length of integer work array
      LRW = 20+nvars*16 ! Length of work array
      LIW = 20          ! Length of integer work array
      !MF = 2            ! Internally generate full Jacobian if needed
      MF = 10            ! Use non-stiff method

      ! Run parameters
      kount = 0
      T     = x1
      TOUT = T ! Always begin by outputing the starting state

      ! LSODE optional inputs
      Do i=5, 10
         RWORK(i) = 0.0D0
         IWORK(i) = 0
      End do
      RWORK(1) = max(TOUT, 10.0D0)  ! Limit not to shoot beyond. 
                                    ! Set to current TOUT/table limit
      IWORK(6) = 1000000  ! Maximum number of steps

      If (outflag) Then 
         kmax = 99 ! Output at 100 log spaced temperatures.
         dxsav = abs(log(x2)-log(x1))/kmax
      Else 
         kmax = 0
      End if

      If (kmax.gt.KMAXX) Then
              Print *,kmax,KMAXX
              Stop "Not enough memory allocated to store internal
     *              variables of the integrator. Increase KMAXX."
      End if

      Do 11 i=1,nvars ! Copy into interface variables
         y(i)=ystart(i)
 11   Continue

      !---------------------------------------------------------
      ! Boilerplate code to output variables, here initial state
      ! Perform time evolution
      Call DLSODE(derivs, nvars, y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     *            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JDUM, MF)
      ! Check correct execution
      IF (ISTATE .LT. 0) GO TO 80
      ! Copy into integrator's module variables
      kount = kount + 1
      xp(kount) = T
      Do i=1, nvars
         yp(i,kount) = y(i)
      End do
      Do i=1, NZS
         yyp(i,kount) = yy(i)
      End do
      ! Output state
      !Call outputode(kount)! <t> Converted function </t>
      !---------------------------------------------------------

      If (kmax.eq.0) Then 
              TOUT = x2 ! Do not output intermediate states
              RWORK(1) = max(TOUT, 10.0D0) ! Don't overshoot end
              ! Perform time evolution
              Call DLSODE(derivs, nvars, y, T, TOUT, ITOL, RTOL, ATOL,
     *                    ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW,
     *                    JDUM, MF)
              ! Check correct execution
              IF (ISTATE .LT. 0) GO TO 80
              ! Copy into integrator's variables
              kount = kount + 1
              xp(kount) = T
              Do i=1, nvars
                 yp(i,kount) = y(i)
              End do
              Do i=1, NZS
                 yyp(i,kount) = yy(i)
              End do
              ! Output final state
              !Call outputode(kount)! <t> Converted function </t>
      Else ! Output intermediate states
              Do nstp=1,kmax
                 TOUT = TOUT*exp(-dxsav)
                 RWORK(1) = max(TOUT, 10.0D0) ! Don't overshoot end
                 Call DLSODE(derivs, nvars, y, T, TOUT, ITOL, RTOL,
     *                       ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, 
     *                       IWORK, LIW, JDUM, MF)
                 ! Check correct execution
                 IF (ISTATE .LT. 0) GO TO 80
                 ! Copy into integrator's variables
                 kount = kount + 1
                 xp(kount) = T
                 Do i=1, nvars
                    yp(i,kount) = y(i)
                 End do
                 Do i=1, NZS
                    yyp(i,kount) = yy(i)
                 End do
                 ! Output intermediate/final state
                 !Call outputode(kount)! <t> Converted function </t>
                 ISTATE = 2 ! It's no longer the first call
              End do
      End if

      ! Output inspection of LSODE run to screen if successful
      WRITE(6,60) IWORK(11),IWORK(12),IWORK(13),IWORK(19),RWORK(15)
 60   FORMAT(/' No. steps =',I8,'  No. f-s =',I8,'  No. J-s =',I8/
     1   ' Method last used =',I2,'   Last switch was at T=',D12.4)
      Return

      ! Error message if run fails
 80   WRITE(6,90) ISTATE
 90   FORMAT(///' Error halt.. ISTATE =',I3)
      Return
!******************************************************************!}}}
      END ! SUBROUTINE odepackint
!******************************************************************

!******************************************************************
      !SUBROUTINE contour
!!******************************************************************!{{{
!C Calls the integrator and determines the contours in the \delta s^2 
!C vs sin^2 \theta plane. Older version.

      !Use Phconst!{{{
      !Use Varpr
      !Implicit none

!!---- Local variables
      !Double precision  omegas
      !Integer i1, i2!}}}
      
      !omegas = 0.0D0 ! Sterile closure fraction
!!---- Start at min sin^2 \theta and max m_s
      !i1 = 0
      !i2 = 0
      !! Find contour on the first row
      !Do while (omegas.lt.omegacont)
         !If (i1.gt.ni1) Then
                 !Print *,ms,s2,omegas
                 !Stop "No omega contour found on first row,
      !*                 consider changing the range of leptasymi."
         !End if
         !Call assignpm(i1, i2)
         !Call evalstd(omegas)
         !Print *,'initial run right', i1, i2, s2, ms, omegas
         !i1 = i1 + 1
      !End do
      !i1 = i1 - 1
      !Write(12,130) i1, i2, s2, ms, omegas
      !Call FLUSH(12) ! Non-standard
!!---- Now proceed with the other rows 
      !Do i2 = 1, ni2
         !Call searchrow(omegas, omegacont, i1, i2)
         !Write(12, 130) i1, i2, s2, ms, omegas
         !Call FLUSH(12) ! Non-standard
      !End do   

!130  Format(2(i8), 3(1pe13.5))

      !Return
!!******************************************************************!}}}
      !end !subroutine contour
!******************************************************************

!******************************************************************
      !SUBROUTINE searchrow(omegas, omegacont, i1, i2)
!!******************************************************************!{{{
!C Starts at \sin^2 \theta and m_s indexed by i1, i2, and searches
!C along i1, i.e. moves along constant m_s. It assumes that the sterile 
!C production, & \Omega_s h^2 increase with \sin^2\theta. Outputs 
!C maximum i1 s.t. omegas \leq omegacont.      

      !Use Varpr!{{{
      !Implicit None

      !Double precision omegas, omegacont
      !Integer i1, i2!}}}

      !Call assignpm(i1, i2)
      !Call evalstd(omegas)
      !If (omegas.lt.omegacont) Then
              !Do while (omegas.lt.omegacont)
                 !If (i1.gt.ni1) Then
                         !Print *, i1, i2, s2, ms, omegas
                         !Stop "No omega found on row, contour will be
      !*                       found at larger values of \sin^2 \\theta."
                 !End if
                 !Call assignpm(i1, i2)
                 !Call evalstd(omegas)
                 !i1 = i1 + 1
              !End do
              !i1 = i1 - 1
              !Call assignpm(i1, i2)
              !Return
      !Else if (omegas.gt.omegacont) Then
              !Do while (omegas.gt.omegacont)
                 !If (i1.lt.0) Then
                         !Print *, i1, i2, s2, ms, omegas
                         !Stop "No omega found on row, contour will be
      !*                       found at lower values of \sin^2 \\theta."
                 !End if
                 !Call assignpm(i1, i2)
                 !Call evalstd(omegas)
                 !i1 = i1 - 1
              !End do
              !If (omegas.eq.omegacont) Then
                      !i1 = i1 + 1
              !End if
              !Call assignpm(i1, i2)
              !Return
      !Else
              !Call assignpm(i1, i2)
              !Return
      !End if

      !Return
!!******************************************************************!}}}
      !END Subroutine !searchrow
!******************************************************************


