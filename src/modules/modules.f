!******************************************************************
      Module Phconst
!******************************************************************!{{{
! Some Physical and mathematical constants
      Implicit None

      Double precision pi, zeta3, nfactor
      Double precision mp, Gf, hbar
      ! \Omega_DM h^2, T_cmb, (critical density/h^2) today
      Double precision omegacont, Tcmb, rhocrit 
      Double precision me, mmu, mtau
      Double precision mnue, mnumu, mnutau
      Double precision mz, mw, s2w

      Parameter (pi=3.1415926535897931, zeta3=1.202056903159594)
      ! factor relating neutrino and photon temps after e^-e^+ epoch
      Parameter (nfactor=(4.0D0/11.0D0)**(1.0D0/3.0D0))
      ! MeV^-2, (mp = 1/sqrt{G}) MeV, MeV-s
      Parameter (Gf=1.1663787D-11,mp=1.22089201D+22,hbar=6.58211928D-22)
      ! Tcmb in MeV 
      Parameter (omegacont=0.1188d0, Tcmb=2.348653946D-10)
      ! rhocrit/h^2 in MeV^4
      Parameter (rhocrit=8.095921884D-35)

      Parameter (me=0.511D0,mmu=1.05658D+2,mtau=1.77682D+3)
      Parameter (mnue=0.0D0,mnumu=0.0D0,mnutau=0.0D0)
      Parameter (mz=91.1876D+3,mw=80.385D+3,s2w=0.23126)

      Save
!******************************************************************!}}}
      End Module ! Phconst
!******************************************************************

!******************************************************************
      Module Params
!******************************************************************!{{{
! Some parameters controlling the code
      Implicit None

      ! Temperature demarcating low and high T, and lower cutoff for 
      ! high temperatures (if so desired, else set to qcdtemp+delta)
      Double precision qcdtemp, hightcutoff, delta

      Parameter (delta=1.0D-5)

      Parameter (qcdtemp=178.0D0)
      Parameter (hightcutoff = 1.0D+3 - delta)

      ! Filenames
      Character*50 ratefile, gstarfile, pchifile, EFfile
      Character *50 Lefile, Lmufile, Ltaufile, chifile

      Parameter (ratefile='data/tables/rate_total.dat')
      Parameter (gstarfile='data/tables/SMgstar.dat')
      Parameter (pchifile='data/tables/ParticleChiTable.dat')
      Parameter (EFfile='data/tables/EFermi.dat')
      Parameter (Lefile='data/tables/dmudLe_new.dat')
      Parameter (Lmufile='data/tables/dmudLmu_new.dat')
      Parameter (Ltaufile='data/tables/dmudLtau_new.dat')
      Parameter (chifile='data/tables/ChiTable_alltemp_new.dat')

      ! nlpTrate and nlTrate are the numbers of scaled momenta (p/T) 
      ! and temperatures at which rates have been sampled. 
      Integer nlpTrate,nlTrate
      Parameter (nlpTrate=101,nlTrate=45)
      !Parameter (nlpTrate=200,nlTrate=200)

      ! ngstar is the number of g_* vs log(temperature)
      Integer ngstar
      Parameter (ngstar=461)

      Integer NZS,nvar ! Number of momentum bins, and ODE variables
      Parameter (NZS=1000,nvar=2*NZS+3)
      ! extra three parameters are time, efolds, and lepton asymmetry 
      ! in units of T^3

      ! nscpx is the number of tabulated dmu/dL_X vs temperature, 
      ! nscps the same for chi_Q_2, and BQ_11
      Integer nscpe, nscpmu, nscptau
      Integer nscps

      Parameter (nscpe=1999,nscpmu=1999,nscptau=1999)
      Parameter (nscps=1999)

      Save
!******************************************************************!}}}
      End Module ! Params
!******************************************************************

!******************************************************************
      Module Varpr
!******************************************************************!{{{
! Lepton parameters varied at every point on the grid
      Implicit None

      ! Description of parameter grid
      Integer ni1, ni2 ! N(sin^2 \theta) = ni1+1, N(ms^2) = ni2 + 1
      Parameter (ni1=100,ni2=100)

      ! Parameter range variables
      Double precision minms, maxms, minsin, maxsin, ddm2, ds
      !Double precision leptasymi

      ! Current lepton parameter values
      Double precision ms, dm2, s2, c
      Character*8 flavor

      Save
      Contains
      !------------------------------------------------------------
        Subroutine assignpm(i1, i2)
        !---- Assigns m_s, and \sin^2 \theta!{{{

          Integer i1,i2

          ms  = maxms*10.0d0**(-dble(i2)*ddm2)
          dm2 = ms**2
          s2  = minsin*10.0d0**(dble(i1)*ds)
          c   = sqrt(1.D0-s2)!}}}
        End Subroutine ! assignpm
      !------------------------------------------------------------
!******************************************************************!}}}
      End Module ! Varpr
!******************************************************************

!******************************************************************
      Module Pathvar
!******************************************************************!{{{
! Parameters describing the system's initial state, and internal state 
! at each time step
      Use Params
      Implicit None

      ! Initial and final states
      Double precision ti, tf, tmevi, tmevf, li
      Double precision yinit(nvar) ! Initial time, efolds, L, and PSD
      Double precision yyinit(NZS) ! Initial momentum bins
 
      ! Current momentum bins
      Double precision yy(NZS)

      Save
      Contains
      !------------------------------------------------------------
        Subroutine Inity(time, leptasym)
        ! Initializes time, leptasym and null-initializes efolds!{{{
        ! and sterile PSD

          Double precision time, leptasym
          Integer i

          li = leptasym

          yinit(1) = time
          yinit(2) = 0.0D0
          yinit(3) = leptasym
          Do i=4,nvar
             yinit(i)=0.0D0
          End do!}}}
        End Subroutine ! Inity
      !------------------------------------------------------------

      !------------------------------------------------------------
        Subroutine Inityy(T, pmin, pmax)
        ! Initializes NZS momentum bins at the starting temperature!{{{
        ! T (in MeV, say), given pmin and pmax in units of T

          Double precision T, pmin, pmax
          ! Internal variables
          Double precision pminmev, pmaxmev, dyy
          Integer i

          pminmev = pmin*T
          pmaxmev = pmax*T
          dyy = (-dlog10(pminmev)+dlog10(pmaxmev))/dble(NZS-1)

          Do i=1,NZS
             yyinit(i) = 10.D0**(dlog10(pminmev)+dble(i-1)*dyy)
             yy(i) = yyinit(i) ! Not really needed, derivs takes care
          End do!}}}
        End Subroutine ! Inityy
      !------------------------------------------------------------
!******************************************************************!}}}
      End Module ! Pathvar
!******************************************************************

!******************************************************************
      Module Odeintmemory
!******************************************************************!{{{
! Variables the integrator saves
      Implicit None

      Integer MAXSTP,NMAX,KMAXX
      Double precision fTINY
      Parameter (MAXSTP=1000000000,NMAX=3000,KMAXX=1000,fTINY=1.D-90)

      Integer kount, kmax
      Double precision dxsav, xp(KMAXX), yp(NMAX,KMAXX)
      Double precision yyp(NMAX,KMAXX)

!******************************************************************!}}}
      End Module ! Integrator
!******************************************************************

!******************************************************************
      Module FermiDirac
!******************************************************************!{{{
!Defines a number of splines and functions to do with FD distribs
      Use Params
      Implicit None

      Integer nlpchi,klpchi,nlpchiknots
      Integer nlEF,klEF,nlEFknots
      ! nlpchi is the number of tabulated log(particle number 
      ! susceptibilities) vs masses, nlpchiknots is the no. of knots for
      ! a spline of degree klpchi
      ! nlEF, klEF, and nlEFknots the same for log(E_Fermi Dirac)
      Parameter (nlpchi=2001,klpchi=3)
      Parameter (nlEF=2001,klEF=3)

!---  Data for various splines, and their knots and coefficients
!---  log(particle number susceptibility), log(E_Fermi Dirac), 
!---  strong susceptibility, redistribution vars, and the gstars
      Double precision masslpchit(nlpchi),lpchit(nlpchi)
      Double precision lpchiknots(nlpchi + klpchi + 1),
     *                 lpchicoeffs(nlpchi + klpchi + 1)

      Double precision masslEFt(nlEF),lEFt(nlEF)
      Double precision lEFknots(nlEF + klEF + 1),
     *                 lEFcoeffs(nlEF + klEF + 1)

      Save
      Contains
      !------------------------------------------------------------
        Subroutine InitFD
        ! Reads in the tables, and assigns the splines!{{{

          Integer i

          Open(unit=90,file=trim(pchifile),status='old')
          Open(unit=91,file=trim(EFfile),status='old')

          Read(90,*) ! Ignore header
          Do i=1, nlpchi
             Read(90,60) masslpchit(i),lpchit(i)
          End do

          Read(91,*) ! Ignore header
          Do i=1, nlEF
             Read(91,60) masslEFt(i), lEFt(i)
          End do

 60       Format(1E8.1, 1E19.16)

          Close(unit=90)
          Close(unit=91)

          ! Find knots and coefficients of splines for 
          ! a) log(particle susceptibility) vs mass
          ! b) log(E_fermi dirac) vs mass
          Call SplFit(masslpchit,lpchit,nlpchi,lpchiknots,
     *                lpchicoeffs,nlpchiknots)
          Call SplFit(masslEFt,lEFt,nlEF,lEFknots,lEFcoeffs,
     *                nlEFknots)!}}}
        End Subroutine ! Initfd
      !------------------------------------------------------------

      !------------------------------------------------------------
        Double precision Function pchi(m)
        ! Evaluate particle number susceptibility for a Fermi-Dirac!{{{
        ! distribution with degeneracy factor of one, from pre-computed
        ! spline knots and coefficients for log(chi/T^2). 
        ! Returns chi/T^2 given m/T

          Double precision m
          ! Internal variables
          Double precision val
          Integer nvals, ier
          Parameter (nvals=1)

          CALL splev(lpchiknots,nlpchiknots,lpchicoeffs,
     *               klpchi,m,val,nvals,ier)

          If (ier.gt.0) Then
                  Print *,"ier=",ier
                  Print *,"m/T=",m
                  Stop "Error in function pchi, check if evaluated
     *                  outside temperature range."
          End If

          pchi = exp(val)

          Return!}}}
        End Function ! pchi
      !------------------------------------------------------------

      !------------------------------------------------------------
        Double precision Function EF(m)
        ! Evaluate energy density for a Fermi-Dirac distribution!{{{
        ! with degeneracy factor of one, from pre-computed spline knots 
        ! and coefficients for log(E_Fermi Dirac/T^4). 
        ! Returns E_F/T^4 given m/T

          Double precision m      
          ! Internal variables
          Double precision val
          Integer nvals, ier
          Parameter (nvals=1)
   
          CALL splev(lEFknots,nlEFknots,lEFcoeffs,klEF,m,val,nvals,ier)

          If (ier.gt.0) Then
                  Print *,"ier=",ier
                  Print *,"m/T=",m
                  Stop "Error in function EF, check if evaluated
     *                  outside temperature range."
          End If

          EF = exp(val)

          Return!}}}
        End Function ! EF
      !------------------------------------------------------------
!******************************************************************!}}}
      End Module ! FermiDirac
!******************************************************************

!******************************************************************
      Module Suscep
!******************************************************************!{{{
!Defines a number of splines and functions dealing with redistribution
      Use Params
      Implicit None

      Integer kscpe,nscpeknots,
     *        kscpmu,nscpmuknots,
     *        kscptau,nscptauknots
      Integer kscps,nscpsknots
      ! nscpxknots is the no. of knots for a spline of degree kscpx.
      ! nscpsknots the same for chi_Q_2, and BQ_11
      Parameter (kscpe=3)
      Parameter (kscpmu=3)
      Parameter (kscptau=3)
      Parameter (kscps=3)

!---  Data for various splines, and their rates and coefficients
      ! 8 because the redistribution tables have 8 columns: \mu_X,
      ! \mu_{\nu_X}, \mu_Q, \mu_B
      Double precision tempet(nscpe),scpet(nscpe,8)
      Double precision scpeknots(nscpe + kscpe + 1, 8),
     *                 scpecoeffs(nscpe + kscpe + 1, 8)
      Double precision tempmut(nscpmu),scpmut(nscpmu,8)
      Double precision scpmuknots(nscpmu + kscpmu + 1, 8),
     *                 scpmucoeffs(nscpmu + kscpmu + 1, 8)
      Double precision temptaut(nscptau),scptaut(nscptau,8)
      Double precision scptauknots(nscptau + kscptau + 1, 8),
     *                 scptaucoeffs(nscptau + kscptau + 1, 8)

      ! 2 because all we need is chi^Q_2, and chi^BQ_11
      Double precision tempst(nscps),scpst(nscps,2)
      Double precision scpsknots(nscps + kscps + 1, 2),
     *                 scpscoeffs(nscps + kscps + 1, 2)

      Save
      Contains
      !------------------------------------------------------------
        Subroutine InitSuscep
        ! Reads in the tables and assigns the splines!{{{

          Integer i
          Double precision var

          Open(unit=90,file=trim(Lefile),status='old')
          Open(unit=91,file=trim(Lmufile),status='old')
          Open(unit=92,file=trim(Ltaufile),status='old')
          Open(unit=93,file=trim(chifile),status='old')

          Read(90,*) ! Ignore header
          Do i=1, nscpe
             Read(90,70) tempet(i), scpet(i,1), scpet(i,2), 
     *                   scpet(i,3), scpet(i,4), scpet(i,5), 
     *                   scpet(i,6), scpet(i,7), scpet(i,8)
          End do

          Read(91,*) ! Ignore header
          Do i=1, nscpmu
             Read(91,70) tempmut(i), scpmut(i,1), scpmut(i,2), 
     *                   scpmut(i,3), scpmut(i,4), scpmut(i,5), 
     *                   scpmut(i,6), scpmut(i,7), scpmut(i,8)
          End do

          Read(92,*) ! Ignore header
          Do i=1, nscptau
             Read(92,70) temptaut(i), scptaut(i,1), scptaut(i,2), 
     *                   scptaut(i,3), scptaut(i,4), scptaut(i,5), 
     *                   scptaut(i,6), scptaut(i,7), scptaut(i,8)
          End do

          Read(93,*) ! Ignore header
          Do i=1, nscps
             Read(93,80) tempst(i), var, scpst(i,1), scpst(i,2)
          End do

 70       Format(9E16.4)
 80       Format(1E8.1, 3E24.15)

          Close(unit=90)
          Close(unit=91)
          Close(unit=92)
          Close(unit=93)

          ! Find knots and coefficients of splines for
          ! a) dmu_\alpha/dL_X
          ! b) chi^Q_2 and chi^BQ_11
          ! Note arrays are stored in column major
          Do i=1, 8
             Call SplFit(tempet,scpet(1,i),nscpe,scpeknots(1,i),
     *                   scpecoeffs(1,i),nscpeknots)
             Call SplFit(tempmut,scpmut(1,i),nscpmu,scpmuknots(1,i),
     *                   scpmucoeffs(1,i),nscpmuknots)
             Call SplFit(temptaut,scptaut(1,i),nscptau,scptauknots(1,i),
     *                   scptaucoeffs(1,i),nscptauknots)
          End do
          Do i=1, 2
             Call SplFit(tempst,scpst(1,i),nscps,scpsknots(1,i),
     *                   scpscoeffs(1,i),nscpsknots)
          End do!}}}
        End Subroutine ! InitSuscep
      !------------------------------------------------------------

      !------------------------------------------------------------
        Double precision Function dmudlx(T, alpha, X)
        ! Evaluates dmu_{\alpha}/dL_X in units of T^{-2} as a function !{{{
        ! of temperature T (in MeV), from precalculated spline knots &
        ! coefficients. \alpha = 1, ... 8 runs over e, mu, tau, nu_e, 
        ! nu_mu, nu_tau, Q, B. X = e,mu or tau

          Double Precision T
          Integer alpha 
          Character*8 X
          ! Internal variables
          Double precision val, lowT, highT
          Integer nvals, ier ! number of values requested and flag
          Parameter (nvals=1)

          highT = 0.0D0
          lowT = 0.0D0

          If (X(1:1) .eq. 'e') Then
                  lowT = tempet(1)
                  highT = tempet(nscpe)
          Else If (X(1:1) .eq. 'm') Then
                  lowT = tempmut(1)
                  highT = tempmut(nscpmu)
          Else If (X(1:1) .eq. 't') Then
                  lowT = temptaut(1)
                  highT = temptaut(nscptau)
          Else
                  Print *,X
                  Stop "Function dmudlx needs X = e, mu or tau."
          End If

          If ((alpha .lt. 1) .or. (alpha .gt. 8)) Then
                  Print *,alpha
                  Stop "Function dmudlx needs alpha \\in [1,8]."
          End If
    
          If (X(1:1) .eq. 'e') Then
                  CALL splev(scpeknots(1,alpha),nscpeknots,
     *                       scpecoeffs(1,alpha),kscpe,T,val,nvals,ier)
          Else If (X(1:1) .eq. 'm') Then
                  CALL splev(scpmuknots(1,alpha),nscpmuknots,
     *                      scpmucoeffs(1,alpha),kscpmu,T,val,nvals,ier)
          Else If (X(1:1) .eq. 't') Then
                  CALL splev(scptauknots(1,alpha),nscptauknots,
     *                    scptaucoeffs(1,alpha),kscptau,T,val,nvals,ier)
          Else
                  Print *,X
                  Stop "Function dmudlx needs X = e, mu or tau."
          End If

          If (ier.gt.0) Then
                  Print *,"ier=",ier
                  Print *,"T=",T
                  Stop "Error in function dmudlx, check if evaluated
     *                  outside temperature range."
          End If
         
          dmudlx = val

          Return!}}}
        End Function ! dmudlx
      !------------------------------------------------------------

      !------------------------------------------------------------
        Double precision Function scps(T, x)
        ! Evaluate susceptibilities \chi^Q_2 and \chi^BQ_11 of the !{{{
        ! strong fluid (in units of T^2) as a function of temperature 
        ! (in MeV) from pre-calculated spline knots and coefficients.

          Double precision T
          Integer x ! \chi^Q_2 and \chi^BQ_11 for x = 1 and 2
          ! Internal variables
          Double precision val
          Integer nvals, ier ! number of values and error flag
          Parameter (nvals=1)

          If ((x.ne.1).and.(x.ne.2)) Then
                  Print *,x
                  Stop "Function scps defined only for x = 1 or 2."
          End if
    
          CALL splev(scpsknots(1,x),nscpsknots,scpscoeffs(1,x),
     *               kscps,T,val,nvals,ier)

          If (ier.gt.0) Then
                  Print *,"ier=",ier
                  Print *,"T=",T
                  Print *,"x=",x
                  Stop "Error in function scps, check if evaluated
     *                  outside temperature range."
          End If

          scps = val

          Return!}}}
        End Function ! scps
      !------------------------------------------------------------
!******************************************************************!}}}
      End Module ! Suscep
!******************************************************************

!******************************************************************
      Module Rates
!******************************************************************!{{{
! Defines splines for scattering rates
      Use Params
      Implicit None

      Integer klpTrate,klTrate
      Integer nlpTrateknots,nlTrateknots
      ! nlpTrateknots and nlTrateknots are the numbers of knots for a 
      ! 2D spline of log(scaled rate) vs log(p/T) and log(T/MeV) 
      ! (coordinate style, array style is reversed) and with degrees 
      ! klpTrate and klTrate
      Parameter (klpTrate=3,klTrate=3)

!---  Input data, subset for spline, and spline knots and coefficients
      Double precision lpTratet(nlpTrate),lTratet(nlTrate)
      Double precision lratet(nlpTrate*nlTrate)
!<a>---------------------------------------------------------------
      Integer nlTratesub
      Double precision lTratetsub(nlTrate) ! both larger than needed
      Double precision lratetsub(nlpTrate*nlTrate)
!--------------------------------------------------------------</a>
      Double precision lpTrateknots(nlpTrate + klpTrate + 1),
     *                 lTrateknots(nlTrate + klTrate + 1),
     *                 lratecoeffs(nlpTrate*nlTrate)

!---  Working space for spline evaluation, define once
      ! Maximum size of pT by T array we will be evaluating over
      Integer nlpTrateevalmax, nlTrateevalmax

      Parameter (nlpTrateevalmax=1000,nlTrateevalmax=1)

      Integer lwrkrateeval, kwrkrateeval
      Parameter (lwrkrateeval=nlpTrateevalmax*(klpTrate+1)+
     *           nlTrateevalmax*(klTrate+1))
      Parameter (kwrkrateeval=nlpTrateevalmax+nlTrateevalmax)

      Double precision wrkrateeval(lwrkrateeval)
      Integer iwrkrateeval(kwrkrateeval)

      Save
      Contains
      !------------------------------------------------------------
        Subroutine InitRates
        ! Read in the tables and define splines. When we say z vs !{{{
        ! x and y below, we mean in coordinate style.
        ! We finally make a table of scaled rate vs p/T and T, because:
        ! 1) the tables are stored as rates vs T and p/T 
        ! 2) we read in row by row 
        ! 3) Fortran stores arrays in column major.

          Integer i, j
          Character*15 buf ! Ignore
          Character*30 Trow ! Format string for temperature row
          Character*30 datarow ! Format string for data row

          !<a>-----------------------------------------------------------
          Integer indlowT
          Integer indhighT
          !----------------------------------------------------------</a>

          ! First read in full tables before choosing data and taking logs

          ! Format strings for reading input
          Write(Trow, "(A6,I4,A6)") "(A15, ",nlTrate,"E15.6)"
          Write(datarow, "(A8,I4,A6)") "(E15.6, ",nlTrate,"E15.6)"

          Open(unit=90,file=trim(ratefile),status='old')

          Read(90,*) ! Ignore explanatory header
          Read(90,FMT=trim(Trow)) buf, lTratet ! read in temperatures
          Do i=1, nlpTrate
             Read(90,FMT=trim(datarow)) lpTratet(i), 
     *        (lratet((i-1)*nlTrate+j),j=1,nlTrate)
          End do

          Close(unit=90)

          !<a>-----------------------------------------------------------
          ! Inefficient, but traverse T array and find temperature limits
          ! Need to do this to avoid yet another transpose.
          indlowT = 0
          indhighT = 1
          Do i=1, nlTrate
             If (lTratet(i).le.qcdtemp) indlowT = indlowT + 1
             If (lTratet(i).lt.hightcutoff) indhighT = indhighT + 1
          End do
          nlTratesub = indlowT + nlTrate - indhighT + 1
          ! Now copy subset of working array, traversing it column-wise
          Do i=1, nlpTrate
             ! First below the transition
             Do j=1, indlowT
                lTratetsub(j) = lTratet(j) ! many times...
                lratetsub((i-1)*nlTratesub+j) = lratet((i-1)*nlTrate+j)
             End do
             ! Then above the transition
             Do j=indhighT, nlTrate
                lTratetsub(indlowT+j-indhighT+1) = lTratet(j)
                lratetsub((i-1)*nlTratesub+indlowT+j-indhighT+1) = 
     *                   lratet((i-1)*nlTrate+j)
             End do
          End do
          !----------------------------------------------------------</a>
          !-- <a> Changed everything underneath to sub </a>
          ! Scale rates w/ p. and take logs because the splines 
          ! are of log(scaled rate) vs log(p/T) and log(T)
          ! Traverse in column major 
          Do i=1, nlpTrate
             Do j=1, nlTratesub
                lratetsub((i-1)*nlTratesub+j) = 
     *          log(lratetsub((i-1)*nlTratesub+j)/lpTratet(i))
             End do
          End do
          Do i=1, nlpTrate
             lpTratet(i) = log(lpTratet(i))
          End do
          Do i=1, nlTratesub
             lTratetsub(i) = log(lTratetsub(i))
          End do

          ! Find spline
          Call TwoDimSplFit(lpTratet,lTratetsub,lratetsub,
     *                      nlpTrate,nlTratesub,
     *                      lpTrateknots,lTrateknots,lratecoeffs,
     *                      nlpTrateknots,nlTrateknots)!}}}
        End Subroutine ! InitRates
      !------------------------------------------------------------

      !------------------------------------------------------------
        Subroutine GamFac(npT, p, T, flavor, scrates)
        ! Subroutine returns the scaled scattering rates, !{{{
        ! \Gamma(p)/G_F^2 p T^4 given a set of n input ps and T in MeV
          
          Integer npT
          Double precision p(npT), T, scrates(npT)
          Character*8 flavor
          ! Internal variables
          Double precision logpT(npT), logT, val(npT)
          Integer i, ier, nlT ! ier=Error flag

          Parameter (nlT=1)

          If (flavor(1:1).ne.'m') Then
                  Print *,"flavor=",flavor
                  Stop "For now we only have the rates for mu."
          End if

          Do i=1, npT
             logpT(i) = log(p(i)/T)
          End do
          logT = log(T)

          Call bispev(lpTrateknots,nlpTrateknots,lTrateknots,
     *                nlTrateknots,lratecoeffs,klpTrate,klTrate,
     *                logpT,npT,logT,nlT,val,wrkrateeval,lwrkrateeval,
     *                iwrkrateeval,kwrkrateeval,ier)

          If (ier.gt.0) Then
                  Print *,"ier=",ier
                  Print *,"npT=",npT
                  Print *,"nlpTrateevalmax=",nlpTrateevalmax
                  Print *,"p(1)=",p(1)
                  Print *,"p(npT)=",p(npT)
                  Print *,"T=",T
                  Stop "Error in subroutine gamfac, check if evaluated
     *                  outside temperature range, or if size of p
     *                  array is larger than nlpTrateevalmax."
          End If

          Do i=1, npT
             scrates(i) = exp(val(i))
          End do!}}}
        End Subroutine ! GamFac 
      !------------------------------------------------------------
!******************************************************************!}}}
      End Module ! Rates
!******************************************************************

!******************************************************************
      Module DOF
!******************************************************************!{{{
!Defines a number of splines for relativistic degrees of freedom
      Use Params
      Implicit None

      Integer kgstar,ngstarknots
      ! ngstarknots is the no. of knots for a spline of degree kgstar.
      Parameter (kgstar=3)

!---  Various splines
      ! 2 because we have g_* and g_*,s
      Double precision ltempgstart(ngstar),gstart(ngstar,2)
      Double precision gstarknots(ngstar + kgstar + 1, 2),
     *                 gstarcoeffs(ngstar + kgstar + 1, 2)

      Save
      Contains
      !------------------------------------------------------------
        Subroutine InitDOF
        ! Read in the tables and define splines!{{{

          Integer i

          Open(unit=90,file=trim(gstarfile),status='old')

          Read(90,*) ! Ignore header
          Read(90,*) ! Ignore header
          Do i=1, ngstar
             Read(90,90) ltempgstart(i), gstart(i,1), gstart(i,2)
          End do
      
 90       Format(3E24.16)

          Close(unit=90)

          ! Logscale temperatures
          Do i=1, ngstar
             ltempgstart(i) = log(ltempgstart(i))
          End do

          ! Find knots and coefficients of splines for 
          ! g_* and g_*,s vs log(T/MeV)      
          ! Note arrays are stored in column major
          Do i=1, 2
             Call SplFit(ltempgstart,gstart(1,i),ngstar,gstarknots(1,i),
     *                   gstarcoeffs(1,i),ngstarknots)
          End do!}}}
        End Subroutine ! InitDOF
      !------------------------------------------------------------

      !------------------------------------------------------------
        Double precision Function gstar(T, alpha, beta)
        ! Evaluates (d/dT)^\beta(g_*/g_*,s) for alpha=1/2, as a function!{{{
        ! of temperature (in MeV) from precalculated spline knots and 
        ! coefficients, for beta = 0 or 1. 
        ! The derivative is in units of MeV^{-1}.

          Double precision T
          Integer alpha, beta
          ! Internal variables
          Double precision logT, der, wrk(ngstarknots), val
          Integer nvals, ier
          Parameter (nvals=1)

          logT = log(T)
    
          If ((alpha.ne.1) .and. (alpha.ne.2)) Then
                  Print *,alpha
                  Stop "Function gstar needs alpha = 1 or 2."
          End if

          If (beta.eq.0) Then
                  call splev(gstarknots(1,alpha),ngstarknots,
     *                       gstarcoeffs(1,alpha),kgstar,logT,
     *                       val,nvals,ier)
          Else if (beta.eq.1) Then
                  call splder(gstarknots(1,alpha),ngstarknots,
     *                        gstarcoeffs(1,alpha),kgstar,beta,
     *                        logT,der,nvals,wrk,ier)
                  val = (1.0d0/T)*der
          Else
                  Print *,beta
          ! As a curiosity, the general formula for converting 
          ! derivatives in terms of Log(x) to those in x involves 
          ! Stirling Numbers of the First kind :)
                  Stop "Function gstar only evaluates function or
     *                  its first derivative, i.e., beta = 0 or 1."
          End if

          If (ier.gt.0) Then
                  Print *,"ier=",ier
                  Print *,"T=",T
                  Stop "Error in function gstar, check if evaluated 
     *                  outside temperature range"
          End If

          gstar = val

          Return!}}}
        End Function ! gstar
      !------------------------------------------------------------
!******************************************************************!}}}
      End Module ! DOF
!******************************************************************
