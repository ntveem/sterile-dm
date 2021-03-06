! This file contains utilities for finding the root of a one-dimensional
! function using Brent's method, and for one and two dimensional spline 
! fitting and evaluation. The first part is a set of interface 
! functions for convenience, and the second part is a set of routines      
! obtained from http://www.netlib.org/
!---------------------------------------------------------------------

! Begin interface functions

!---------------------------------------------------------------------
      SUBROUTINE SplFit(x,y,m,t,c,n)!{{{
C Subroutine SplFit takes in data (x(i),y(i)), i=1,...,m and finds
C spline of degree k=3 that interpolates between the data. The knots
C and control points of the spline are stored within t(j), c(j) 
C with j=1,...,n=m+k+1=m+4, the total no. of knots is stored in n.

      Implicit none!{{{
C k is the spline's degree, nest is an upper bound on the knot number,
C also the length of the knot and ctrl point vectors, xb and xe are 
C lower and upper limits on x coordinates      
      Integer m, k, nest ! nest = m+k+1
      Double precision x(m), y(m)
      Double precision xb, xe

      Parameter (k=3)

C Things which are modified by the code     
C n is the final no. of knots, ier is a success flag
C t, c are knots and control points, fp is sum of residuals
      Integer n, ier
      Double precision t(m+k+1), c(m+k+1), fp

C Flags controlling code's operation      
C iopt >= 0 makes it choose its own knots, s=0 makes it interpolating
C w(i) i=1,..,m is a set of weights, we make them all unity     
      Integer iopt, i
      Double precision s, w(m)

      Parameter (iopt=0,s=0.0E0)

C Internal working space
      Double precision wrk(m*(k+1)+(m+k+1)*(7+3*k))
      Integer lwrk, iwrk(m+k+1)!}}}
 
      do i=1, m
         w(i) = 1
      end do

      nest = m+k+1
      xb = minval(x)
      xe = maxval(x)

      lwrk=m*(k+1)+(m+k+1)*(7+3*k)

      CALL curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,
     * wrk,lwrk,iwrk,ier)

      If (ier.gt.0) Then
              Print *,"ier=",ier
              Print *,"iopt=",iopt
              Print *,"m=",m
              Print *,"xb=",xb
              Print *,"xe=",xe
              Print *,"k=",k
              Print *,"s=",s
              Print *,"nest=",nest
              Print *,"n=",n
              Print *,"fp=",fp
              Print *,"lwrk=",lwrk
              Stop "Error in SplFit."
      End If
      END SUBROUTINE ! SplFit!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      SUBROUTINE TwoDimSplFit(x,y,z,mx,my,tx,ty,c,nx,ny)!{{{
C Subroutine TwoDimSplFit takes in data z(i,j) on the rectangular grid 
C (x(i),y(j)), i=1,...,mx;j=1,...,my and finds a two dimensional 
C interpolating spline approximant of degree kx=ky=3 in both directions.
C The knots in the x and y directions are stored in tx(i),
C i=1,...,nx=mx+kx+1=mx+4 and ty(j),j=1,...,ny=my+ky+1=my+4. The 
C b-spline coefficients are stored in the matrix c((ny-ky-1)*(i-1)+j),
C i=1,...,nx-kx-1=mx;j=1,...,ny-ky-1=my, and the knot numbers are stored
C in nx and ny. 

      Implicit none!{{{
C k(x/y) are the spline's degrees, n(x/y)est are upper bounds on the 
C knot numbers. (x/y)b and (x/y)e are lower and upper limits on (x/y) 
C coordinates. 
      Integer mx, my, kx, ky, nxest, nyest ! n(x/y)est = m(x/y)+k(x/y)+1
      Double precision x(mx), y(my), z(mx*my)
      Double precision xb, xe, yb, ye

      Parameter (kx=3,ky=3)

C Things which are modified by the code     
C n(x/y) are final nos. of knots, ier is a success flag
C tx, ty are knots, c are spline coeffs, fp is sum of residuals
      Integer nx, ny, ier
      Double precision tx(mx+kx+1), ty(my+ky+1), c(mx*my), fp

C Flags controlling code's operation      
C iopt >= 0 makes it choose its own knots, s=0 makes it interpolating
      Integer iopt
      Double precision s

      Parameter (iopt=0,s=0.0D0)

C Internal working space
      Double precision wrk(4+(mx+kx+1)*(my+2*kx+5)+(my+ky+1)*(2*ky+5)+
     *                    mx*(kx+1)+my*(ky+1)+max(my,mx+kx+1))
      Integer lwrk, kwrk, iwrk(3+mx+my+(mx+kx+1)+(my+ky+1))!}}}

      nxest = mx+kx+1
      nyest = my+ky+1
      xb = minval(x)
      xe = maxval(x)
      yb = minval(y)
      ye = maxval(y)

      lwrk = 4+(mx+kx+1)*(my+2*kx+5)+(my+ky+1)*(2*ky+5)+
     *       mx*(kx+1)+my*(ky+1)+max(my,mx+kx+1)
      kwrk = 3+mx+my+(mx+kx+1)+(my+ky+1)

      CALL regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
     * nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)

      If (ier.gt.0) Then
              Print *,"ier=",ier
              Print *,"iopt=",iopt
              Print *,"mx=",mx
              Print *,"my=",my
              Print *,"xb=",xb
              Print *,"xe=",xe
              Print *,"yb=",yb
              Print *,"ye=",ye
              Print *,"kx=",kx
              Print *,"ky=",ky
              Print *,"s=",s
              Print *,"nxest=",nxest
              Print *,"nyest=",nyest
              Print *,"nx=",nx
              Print *,"ny=",ny
              Print *,"fp=",fp
              Print *,"lwrk=",lwrk
              Print *,"kwrk=",kwrk
              Stop "Error in TwoDimSplFit."
      End If

      END SUBROUTINE ! TwoDimSplFit!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION NIntegrate(x,y,m,a,b)!{{{
C Function NIntegrate takes in data (x(i),y(i)), i=1,...,m and finds
C the definite integral between limits a and b, of a spline of degree 
C three that interpolates between the data.

      Implicit none!{{{
C k is the spline's degree, n is an upper bound on the knot number,
C also the length of the knot and ctrl point vectors, t, and c
      Integer m, k, n ! n = m+k+1
      Double precision x(m), y(m), a, b

      Parameter (k=3)

C knots and control points
      Double precision t(m+k+1), c(m+k+1), wrk(m+k+1)

C External function
      External splint
      Double precision splint!}}}

      If ((a.lt.x(1)) .OR. (a.gt.x(m)) .OR.
     * (b.lt.x(1)) .OR. (b.gt.x(m))) Then
              Print *,a,b,x(1),x(m)
              Stop "Function NIntegrate evaluated outside range."
      End If

      CALL SplFit(x,y,m,t,c,n) ! Find spline

      NIntegrate = splint(t,n,c,k,a,b,wrk)

      RETURN
      END FUNCTION ! NIntegrate!}}}
!---------------------------------------------------------------------

! End interface functions

!---------------------------------------------------------------------
      double precision function zeroin(ax,bx,f,tol)!{{{
      double precision ax,bx,f,tol
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result (.ge.0.)
c
c  output..
c
c  zeroin abscissa approximating a zero of  f  in the interval ax,bx
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  this is checked, and an error message is printed if this is not
c  satisfied.   zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
c  the  relative machine precision defined as the smallest representable
c  number such that  1.+macheps .gt. 1.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice-hall, inc. (1973).
c
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision  dabs, d1mach
   10 eps = d1mach(4)
      tol1 = eps+1.0d0
c
      a=ax
      b=bx
      fa=f(a)
      fb=f(b)
c     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
      if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20
         write(6,2500)
2500     format(1x,'f(ax) and f(bx) do not have different signs,',
     1             ' zeroin is aborting')
         return
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if (dabs(fc).ge.dabs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0d0*eps*dabs(b)+0.5d0*tol
      xm = 0.5d0*(c-b)
      if ((dabs(xm).le.tol1).or.(fb.eq.0.0d0)) go to 150
c
c see if a bisection is forced
c
      if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
      if (a.ne.c) go to 60
c
c linear interpolation
c
      p=2.0d0*xm*s
      q=1.0d0-s
      go to 70
c
c inverse quadratic interpolation
c
   60 q=fa/fc
      r=fb/fc
      p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
      q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
   70 if (p.le.0.0d0) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.(p.ge.
     *dabs(0.5d0*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if (dabs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0d0) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
  140 fb=f(b)
      if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20
      go to 30
  150 zeroin=b
      return
      end!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION D1MACH(I)!{{{
      INTEGER I
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C  D1MACH( 5) = LOG10(B)
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, CRAY1(38), J
      COMMON /D9MACH/ CRAY1
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      DOUBLE PRECISION DMACH(5)
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
C  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
C  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
C  MANY MACHINES YET.
C  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
C  ON THE NEXT LINE
      DATA SC/0/
C  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
C  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
C          mail netlib@research.bell-labs.com
C          send old1mach from blas
C  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS.
C      DATA SMALL(1),SMALL(2) /    8388608,           0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
C      DATA DIVER(1),DIVER(2) /  620756992,           0 /
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
C
C     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532
     *       .AND. SMALL(1) .EQ. -448790528) THEN
*           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935
     *       .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943
     *       .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 1120022684
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** CONVEX C-1 ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 1019215872
            RIGHT(2) = 0
            DIVER(1) = 1020264448
            DIVER(2) = 0
            LOG10(1) = 1072907283
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 815547074
     *       .AND. SMALL(2) .EQ. 58688) THEN
*           *** VAX G-FLOATING ***
            SMALL(1) = 16
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 15552
            RIGHT(2) = 0
            DIVER(1) = 15568
            DIVER(2) = 0
            LOG10(1) = 1142112243
            LOG10(2) = 2046775455
         ELSE
            DMACH(2) = 1.D27 + 1
            DMACH(3) = 1.D27
            LARGE(2) = LARGE(2) - RIGHT(2)
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
               CRAY1(1) = 67291416
               DO 10 J = 1, 20
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 10               CONTINUE
               CRAY1(22) = CRAY1(21) + 321322
               DO 20 J = 22, 37
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 20               CONTINUE
               IF (CRAY1(38) .EQ. SMALL(1)) THEN
*                  *** CRAY ***
                  CALL I1MCRY(SMALL(1), J, 8285, 8388608, 0)
                  SMALL(2) = 0
                  CALL I1MCRY(LARGE(1), J, 24574, 16777215, 16777215)
                  CALL I1MCRY(LARGE(2), J, 0, 16777215, 16777214)
                  CALL I1MCRY(RIGHT(1), J, 16291, 8388608, 0)
                  RIGHT(2) = 0
                  CALL I1MCRY(DIVER(1), J, 16292, 8388608, 0)
                  DIVER(2) = 0
                  CALL I1MCRY(LOG10(1), J, 16383, 10100890, 8715215)
                  CALL I1MCRY(LOG10(2), J, 0, 16226447, 9001388)
               ELSE
                  WRITE(*,9000)
                  STOP 779
                  END IF
            ELSE
               WRITE(*,9000)
               STOP 779
               END IF
            END IF
         SC = 987
         END IF
*    SANITY CHECK
      IF (DMACH(4) .GE. 1.0D0) STOP 778
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      D1MACH = DMACH(I)
      RETURN
 9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/
     *' appropriate for your machine.')
* /* Standard C source for D1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*double d1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return DBL_MIN;
*	  case 2: return DBL_MAX;
*	  case 3: return DBL_EPSILON/FLT_RADIX;
*	  case 4: return DBL_EPSILON;
*	  case 5: return log10((double)FLT_RADIX);
*	  }
*	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
*	exit(1); return 0; /* some compilers demand return values */
*}
      END
      SUBROUTINE I1MCRY(A, A1, B, C, D)
**** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,!{{{
     * wrk,lwrk,iwrk,ier)
c  given the set of data points (x(i),y(i)) and the set of positive
c  numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
c  approximation of degree k on the interval xb <= x <= xe.
c  if iopt=-1 curfit calculates the weighted least-squares spline
c  according to a given set of knots.
c  if iopt>=0 the number of knots of the spline s(x) and the position
c  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
c  ness of s(x) is then achieved by minimalizing the discontinuity
c  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
c  n-k-1. the amount of smoothness is determined by the condition that
c  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
c  negative constant, called the smoothing factor.
c  the fit s(x) is given in the b-spline representation (b-spline coef-
c  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
c  subroutine splev.
c
c  calling sequence:
c     call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,
c    * lwrk,iwrk,ier)
c
c  parameters:
c   iopt  : integer flag. on entry iopt must specify whether a weighted
c           least-squares spline (iopt=-1) or a smoothing spline (iopt=
c           0 or 1) must be determined. if iopt=0 the routine will start
c           with an initial set of knots t(i)=xb, t(i+k+1)=xe, i=1,2,...
c           k+1. if iopt=1 the routine will continue with the knots
c           found at the last call of the routine.
c           attention: a call with iopt=1 must always be immediately
c           preceded by another call with iopt=1 or iopt=0.
c           unchanged on exit.
c   m     : integer. on entry m must specify the number of data points.
c           m > k. unchanged on exit.
c   x     : double precision array of dimension at least (m). before entry, x(i)
c           must be set to the i-th value of the independent variable x,
c           for i=1,2,...,m. these values must be supplied in strictly
c           ascending order. unchanged on exit.
c   y     : double precision array of dimension at least (m). before entry, y(i)
c           must be set to the i-th value of the dependent variable y,
c           for i=1,2,...,m. unchanged on exit.
c   w     : double precision array of dimension at least (m). before entry, w(i)
c           must be set to the i-th value in the set of weights. the
c           w(i) must be strictly positive. unchanged on exit.
c           see also further comments.
c   xb,xe : double precision values. on entry xb and xe must specify the boundaries
c           of the approximation interval. xb<=x(1), xe>=x(m).
c           unchanged on exit.
c   k     : integer. on entry k must specify the degree of the spline.
c           1<=k<=5. it is recommended to use cubic splines (k=3).
c           the user is strongly dissuaded from choosing k even,together
c           with a small s-value. unchanged on exit.
c   s     : double precision.on entry (in case iopt>=0) s must specify the smoothing
c           factor. s >=0. unchanged on exit.
c           for advice on the choice of s see further comments.
c   nest  : integer. on entry nest must contain an over-estimate of the
c           total number of knots of the spline returned, to indicate
c           the storage space available to the routine. nest >=2*k+2.
c           in most practical situation nest=m/2 will be sufficient.
c           always large enough is  nest=m+k+1, the number of knots
c           needed for interpolation (s=0). unchanged on exit.
c   n     : integer.
c           unless ier =10 (in case iopt >=0), n will contain the
c           total number of knots of the spline approximation returned.
c           if the computation mode iopt=1 is used this value of n
c           should be left unchanged between subsequent calls.
c           in case iopt=-1, the value of n must be specified on entry.
c   t     : double precision array of dimension at least (nest).
c           on succesful exit, this array will contain the knots of the
c           spline,i.e. the position of the interior knots t(k+2),t(k+3)
c           ...,t(n-k-1) as well as the position of the additional knots
c           t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=t(n)=xe needed for
c           the b-spline representation.
c           if the computation mode iopt=1 is used, the values of t(1),
c           t(2),...,t(n) should be left unchanged between subsequent
c           calls. if the computation mode iopt=-1 is used, the values
c           t(k+2),...,t(n-k-1) must be supplied by the user, before
c           entry. see also the restrictions (ier=10).
c   c     : double precision array of dimension at least (nest).
c           on succesful exit, this array will contain the coefficients
c           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
c   fp    : double precision. unless ier=10, fp contains the weighted sum of
c           squared residuals of the spline approximation returned.
c   wrk   : double precision array of dimension at least (m*(k+1)+nest*(7+3*k)).
c           used as working space. if the computation mode iopt=1 is
c           used, the values wrk(1),...,wrk(n) should be left unchanged
c           between subsequent calls.
c   lwrk  : integer. on entry,lwrk must specify the actual dimension of
c           the array wrk as declared in the calling (sub)program.lwrk
c           must not be too small (see wrk). unchanged on exit.
c   iwrk  : integer array of dimension at least (nest).
c           used as working space. if the computation mode iopt=1 is
c           used,the values iwrk(1),...,iwrk(n) should be left unchanged
c           between subsequent calls.
c   ier   : integer. unless the routine detects an error, ier contains a
c           non-positive value on exit, i.e.
c    ier=0  : normal return. the spline returned has a residual sum of
c             squares fp such that abs(fp-s)/s <= tol with tol a relat-
c             ive tolerance set to 0.001 by the program.
c    ier=-1 : normal return. the spline returned is an interpolating
c             spline (fp=0).
c    ier=-2 : normal return. the spline returned is the weighted least-
c             squares polynomial of degree k. in this extreme case fp
c             gives the upper bound fp0 for the smoothing factor s.
c    ier=1  : error. the required storage space exceeds the available
c             storage space, as specified by the parameter nest.
c             probably causes : nest too small. if nest is already
c             large (say nest > m/2), it may also indicate that s is
c             too small
c             the approximation returned is the weighted least-squares
c             spline according to the knots t(1),t(2),...,t(n). (n=nest)
c             the parameter fp gives the corresponding weighted sum of
c             squared residuals (fp>s).
c    ier=2  : error. a theoretically impossible result was found during
c             the iteration proces for finding a smoothing spline with
c             fp = s. probably causes : s too small.
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=3  : error. the maximal number of iterations maxit (set to 20
c             by the program) allowed for finding a smoothing spline
c             with fp=s has been reached. probably causes : s too small
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=10 : error. on entry, the input data are controlled on validity
c             the following restrictions must be satisfied.
c             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
c             xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)
c             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
c                         xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe
c                       the schoenberg-whitney conditions, i.e. there
c                       must be a subset of data points xx(j) such that
c                         t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1
c             if iopt>=0: s>=0
c                         if s=0 : nest >= m+k+1
c             if one of these conditions is found to be violated,control
c             is immediately repassed to the calling program. in that
c             case there is no approximation returned.
c
c  further comments:
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the weighted least-squares polynomial of degree k if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in y(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   polynomial and the corresponding upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximation shows more detail) to obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if curfit is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. but, if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   curfit once more with the selected value for s but now with iopt=0.
c   indeed, curfit may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c
c  other subroutines required:
c    fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati,fprota
c
c  references:
c   dierckx p. : an algorithm for smoothing, differentiation and integ-
c                ration of experimental data using spline functions,
c                j.comp.appl.maths 1 (1975) 165-184.
c   dierckx p. : a fast algorithm for smoothing data on a rectangular
c                grid while using spline functions, siam j.numer.anal.
c                19 (1982) 1286-1304.
c   dierckx p. : an improved algorithm for curve fitting with spline
c                functions, report tw54, dept. computer science,k.u.
c                leuven, 1981.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1987
c
c  ..
c  ..scalar arguments..
      double precision xb,xe,s,fp
      integer iopt,m,k,nest,n,lwrk,ier
c  ..array arguments..
      double precision x(m),y(m),w(m),t(nest),c(nest),wrk(lwrk)
      integer iwrk(nest)
c  ..local scalars..
      double precision tol
      integer i,ia,ib,ifp,ig,iq,iz,j,k1,k2,lwest,maxit,nmin
c  ..
c  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1D-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(k.le.0 .or. k.gt.5) go to 50
      k1 = k+1
      k2 = k1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 50
      nmin = 2*k1
      if(m.lt.k1 .or. nest.lt.nmin) go to 50
      lwest = m*k1+nest*(7+3*k)
      if(lwrk.lt.lwest) go to 50
      if(xb.gt.x(1) .or. xe.lt.x(m) .or. w(1).le.0.) go to 50
      do 10 i=2,m
         if(x(i-1).ge.x(i) .or. w(i).le.0.) go to 50
  10  continue
      if(iopt.ge.0) go to 30
      if(n.lt.nmin .or. n.gt.nest) go to 50
      j = n
      do 20 i=1,k1
         t(i) = xb
         t(j) = xe
         j = j-1
  20  continue
      call fpchec(x,m,t,n,k,ier)
      if(ier) 50,40,50
  30  if(s.lt.0.0D0) go to 50
      if(s.eq.0.0D0 .and. nest.lt.(m+k1)) go to 50
      ier = 0
c we partition the working space and determine the spline approximation.
  40  ifp = 1
      iz = ifp+nest
      ia = iz+nest
      ib = ia+nest*k1
      ig = ib+nest*k2
      iq = ig+nest*k2
      call fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,n,t,c,fp,
     * wrk(ifp),wrk(iz),wrk(ia),wrk(ib),wrk(ig),wrk(iq),iwrk,ier)
  50  return
      end!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine fpback(a,z,n,k,c,nest)!{{{
c  subroutine fpback calculates the solution of the system of
c  equations a*c = z with a a n x n upper triangular matrix
c  of bandwidth k.
c  ..
c  ..scalar arguments..
      integer n,k,nest
c  ..array arguments..
      double precision a(nest,k),z(n),c(n)
c  ..local scalars..
      double precision store
      integer i,i1,j,k1,l,m
c  ..
      k1 = k-1
      c(n) = z(n)/a(n,1)
      i = n-1
      if(i.eq.0) go to 30
      do 20 j=2,n
        store = z(i)
        i1 = k1
        if(j.le.k1) i1 = j-1
        m = i
        do 10 l=1,i1
          m = m+1
          store = store-c(m)*a(i,l+1)
  10    continue
        c(i) = store/a(i,1)
        i = i-1
  20  continue
  30  return
      end !}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine fpbspl(t,n,k,x,l,h)!{{{
c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
c  degree k at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  ..
c  ..scalar arguments..
      double precision x
      integer n,k,l
c  ..array arguments..
      double precision t(n),h(6)
c  ..local scalars..
      double precision f,one
      integer i,j,li,lj
c  ..local arrays..
      double precision hh(5)
c  ..
      one = 0.1D+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.0D0
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end !}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine fpchec(x,m,t,n,k,ier)!{{{
c  subroutine fpchec verifies the number and the position of the knots
c  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
c  and the position of the data points x(i),i=1,2,...,m. if all of the
c  following conditions are fulfilled, the error parameter ier is set
c  to zero. if one of the conditions is violated ier is set to ten.
c      1) k+1 <= n-k-1 <= m
c      2) t(1) <= t(2) <= ... <= t(k+1)
c         t(n-k) <= t(n-k+1) <= ... <= t(n)
c      3) t(k+1) < t(k+2) < ... < t(n-k)
c      4) t(k+1) <= x(i) <= t(n-k)
c      5) the conditions specified by schoenberg and whitney must hold
c         for at least one subset of data points, i.e. there must be a
c         subset of data points y(j) such that
c             t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
c  ..
c  ..scalar arguments..
      integer m,n,k,ier
c  ..array arguments..
      double precision x(m),t(n)
c  ..local scalars..
      integer i,j,k1,k2,l,nk1,nk2,nk3
      double precision tj,tl
c  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      ier = 10
c  check condition no 1
      if(nk1.lt.k1 .or. nk1.gt.m) go to 80
c  check condition no 2
      j = n
      do 20 i=1,k
        if(t(i).gt.t(i+1)) go to 80
        if(t(j).lt.t(j-1)) go to 80
        j = j-1
  20  continue
c  check condition no 3
      do 30 i=k2,nk2
        if(t(i).le.t(i-1)) go to 80
  30  continue
c  check condition no 4
      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 80
c  check condition no 5
      if(x(1).ge.t(k2) .or. x(m).le.t(nk1)) go to 80
      i = 1
      l = k2
      nk3 = nk1-1
      if(nk3.lt.2) go to 70
      do 60 j=2,nk3
        tj = t(j)
        l = l+1
        tl = t(l)
  40    i = i+1
        if(i.ge.m) go to 80
        if(x(i).le.tj) go to 40
        if(x(i).ge.tl) go to 80
  60  continue
  70  ier = 0
  80  return
      end !}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,!{{{
     * n,t,c,fp,fpint,z,a,b,g,q,nrdata,ier)
c  ..
c  ..scalar arguments..
      double precision xb,xe,s,tol,fp
      integer iopt,m,k,nest,maxit,k1,k2,n,ier
c  ..array arguments..
      double precision x(m),y(m),w(m),t(nest),c(nest),fpint(nest),
     * z(nest),a(nest,k1),b(nest,k2),g(nest,k2),q(m,k1)
      integer nrdata(nest)
c  ..local scalars..
      double precision acc,con1,con4,con9,cos,half,fpart,fpms,fpold,fp0,
     * f1,f2,f3,one,p,pinv,piv,p1,p2,p3,rn,sin,store,term,wi,xi,yi
      integer i,ich1,ich3,it,iter,i1,i2,i3,j,k3,l,l0,
     * mk1,new,nk1,nmax,nmin,nplus,npl1,nrint,n8
c  ..local arrays..
      double precision h(7)
c  ..function references
      double precision abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota
c  ..
c  set constants
      one = 0.1D+01
      con1 = 0.1D0
      con9 = 0.9D0
      con4 = 0.4D-01
      half = 0.5D0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1: determination of the number of knots and their position     c
c  **************************************************************      c
c  given a set of knots we compute the least-squares spline sinf(x),   c
c  and the corresponding sum of squared residuals fp=f(p=inf).         c
c  if iopt=-1 sinf(x) is the requested approximation.                  c
c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
c    if fp <=s we will continue with the current set of knots.         c
c    if fp > s we will increase the number of knots and compute the    c
c       corresponding least-squares spline until finally fp<=s.        c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmax = m+k+1.                                        c
c    if s > 0 and                                                      c
c      iopt=0 we first compute the least-squares polynomial of         c
c      degree k; n = nmin = 2*k+2                                      c
c      iopt=1 we start with the set of knots found at the last         c
c      call of the routine, except for the case that s > fp0; then     c
c      we compute directly the least-squares polynomial of degree k.   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  determine nmin, the number of knots for polynomial approximation.
      nmin = 2*k1
      if(iopt.lt.0) go to 60
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  determine nmax, the number of knots for spline interpolation.
      nmax = m+k1
      if(s.gt.0.0D0) go to 45
c  if s=0, s(x) is an interpolating spline.
c  test whether the required storage space exceeds the available one.
      n = nmax
      if(nmax.gt.nest) go to 420
c  find the position of the interior knots in case of interpolation.
  10  mk1 = m-k1
      if(mk1.eq.0) go to 60
      k3 = k/2
      i = k2
      j = k3+2
      if(k3*2.eq.k) go to 30
      do 20 l=1,mk1
        t(i) = x(j)
        i = i+1
        j = j+1
  20  continue
      go to 60
  30  do 40 l=1,mk1
        t(i) = (x(j)+x(j-1))*half
        i = i+1
        j = j+1
  40  continue
      go to 60
c  if s>0 our initial choice of knots depends on the value of iopt.
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  polynomial of degree k which is a spline without interior knots.
c  if iopt=1 and fp0>s we start computing the least squares spline
c  according to the set of knots found at the last call of the routine.
  45  if(iopt.eq.0) go to 50
      if(n.eq.nmin) go to 50
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0.gt.s) go to 60
  50  n = nmin
      fpold = 0.0D0
      nplus = 0
      nrdata(1) = m-2
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  60  do 200 iter = 1,m
        if(n.eq.nmin) ier = -2
c  find nrint, tne number of knot intervals.
        nrint = n-nmin+1
c  find the position of the additional knots which are needed for
c  the b-spline representation of s(x).
        nk1 = n-k1
        i = n
        do 70 j=1,k1
          t(j) = xb
          t(i) = xe
          i = i-1
  70    continue
c  compute the b-spline coefficients of the least-squares spline
c  sinf(x). the observation matrix a is built up row by row and
c  reduced to upper triangular form by givens transformations.
c  at the same time fp=f(p=inf) is computed.
        fp = 0.0D0
c  initialize the observation matrix a.
        do 80 i=1,nk1
          z(i) = 0.0D0
          do 80 j=1,k1
            a(i,j) = 0.0D0
  80    continue
        l = k1
        do 130 it=1,m
c  fetch the current data point x(it),y(it).
          xi = x(it)
          wi = w(it)
          yi = y(it)*wi
c  search for knot interval t(l) <= xi < t(l+1).
  85      if(xi.lt.t(l+1) .or. l.eq.nk1) go to 90
          l = l+1
          go to 85
c  evaluate the (k+1) non-zero b-splines at xi and store them in q.
  90      call fpbspl(t,n,k,xi,l,h)
          do 95 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  95      continue
c  rotate the new row of the observation matrix into triangle.
          j = l-k1
          do 110 i=1,k1
            j = j+1
            piv = h(i)
            if(piv.eq.0.0D0) go to 110
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(j,1),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,z(j))
            if(i.eq.k1) go to 120
            i2 = 1
            i3 = i+1
            do 100 i1 = i3,k1
              i2 = i2+1
c  transformations to left hand side.
              call fprota(cos,sin,h(i1),a(j,i2))
 100        continue
 110      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 120      fp = fp+yi**2
 130    continue
        if(ier.eq.(-2)) fp0 = fp
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
c  backward substitution to obtain the b-spline coefficients.
        call fpback(a,z,nk1,k1,c,nest)
c  test whether the approximation sinf(x) is an acceptable solution.
        if(iopt.lt.0) go to 440
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  if f(p=inf) < s accept the choice of knots.
        if(fpms.lt.0.0D0) go to 250
c  if n = nmax, sinf(x) is an interpolating spline.
        if(n.eq.nmax) go to 430
c  increase the number of knots.
c  if n=nest we cannot increase the number of knots because of
c  the storage capacity limitation.
        if(n.eq.nest) go to 420
c  determine the number of knots nplus we are going to add.
        if(ier.eq.0) go to 140
        nplus = 1
        ier = 0
        go to 150
 140    npl1 = nplus*2
        rn = nplus
        if(fpold-fp.gt.acc) npl1 = int(rn*fpms/(fpold-fp))
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
 150    fpold = fp
c  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
c  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.0D0
        i = 1
        l = k2
        new = 0
        do 180 it=1,m
          if(x(it).lt.t(l) .or. l.gt.nk1) go to 160
          new = 1
          l = l+1
 160      term = 0.0D0
          l0 = l-k2
          do 170 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 170      continue
          term = (w(it)*(term-y(it)))**2
          fpart = fpart+term
          if(new.eq.0) go to 180
          store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 180    continue
        fpint(nrint) = fpart
        do 190 l=1,nplus
c  add a new knot.
          call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
c  if n=nmax we locate the knots as for interpolation.
          if(n.eq.nmax) go to 10
c  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 200
 190    continue
c  restart the computations with the new set of knots.
 200  continue
c  test whether the least-squares kth degree polynomial is a solution
c  of our approximation problem.
 250  if(ier.eq.(-2)) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 2: determination of the smoothing spline sp(x).                c
c  ***************************************************                 c
c  we have determined the number of knots and their position.          c
c  we now compute the b-spline coefficients of the smoothing spline    c
c  sp(x). the observation matrix a is extended by the rows of matrix   c
c  b expressing that the kth derivative discontinuities of sp(x) at    c
c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
c  ponding weights of these additional rows are set to 1/p.            c
c  iteratively we then have to determine the value of p such that      c
c  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
c  the least-squares kth degree polynomial corresponds to p=0, and     c
c  that the least-squares spline corresponds to p=infinity. the        c
c  iteration process which is proposed here, makes use of rational     c
c  interpolation. since f(p) is a convex and strictly decreasing       c
c  function of p, it can be approximated by a rational function        c
c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
c  to calculate the new value of p such that r(p)=s. convergence is    c
c  guaranteed by taking f1>0 and f3<0.                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jump of the kth derivative of the
c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
      call fpdisc(t,n,k2,b,nest)
c  initial value for p.
      p1 = 0.0D0
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = 0.0D0
      do 255 i=1,nk1
         p = p+a(i,1)
 255  continue
      rn = nk1
      p = rn/p
      ich1 = 0
      ich3 = 0
      n8 = n-nmin
c  iteration process to find the root of f(p) = s.
      do 360 iter=1,maxit
c  the rows of matrix b with weight 1/p are rotated into the
c  triangularised observation matrix a which is stored in g.
        pinv = one/p
        do 260 i=1,nk1
          c(i) = z(i)
          g(i,k2) = 0.0D0
          do 260 j=1,k1
            g(i,j) = a(i,j)
 260    continue
        do 300 it=1,n8
c  the row of matrix b is rotated into triangle by givens transformation
          do 270 i=1,k2
            h(i) = b(it,i)*pinv
 270      continue
          yi = 0.0D0
          do 290 j=it,nk1
            piv = h(1)
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,g(j,1),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,c(j))
            if(j.eq.nk1) go to 300
            i2 = k1
            if(j.gt.n8) i2 = nk1-j
            do 280 i=1,i2
c  transformations to left hand side.
              i1 = i+1
              call fprota(cos,sin,h(i1),g(j,i1))
              h(i) = h(i1)
 280        continue
            h(i2+1) = 0.0D0
 290      continue
 300    continue
c  backward substitution to obtain the b-spline coefficients.
        call fpback(g,c,nk1,k2,c,nest)
c  computation of f(p).
        fp = 0.0D0
        l = k2
        do 330 it=1,m
          if(x(it).lt.t(l) .or. l.gt.nk1) go to 310
          l = l+1
 310      l0 = l-k2
          term = 0.0D0
          do 320 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 320      continue
          fp = fp+(w(it)*(term-y(it)))**2
 330    continue
c  test whether the approximation sp(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximal number of iterations is reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 340
        if((f2-f3).gt.acc) go to 335
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p=p1*con9 + p2*con1
        go to 360
 335    if(f2.lt.0.0D0) ich3=1
 340    if(ich1.ne.0) go to 350
        if((f1-f2).gt.acc) go to 345
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.0D0) go to 360
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 360
 345    if(f2.gt.0.0D0) ich1=1
c  test whether the iteration process proceeds as theoretically
c  expected.
 350    if(f2.ge.f1 .or. f2.le.f3) go to 410
c  find the new value for p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 360  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
 440  return
      end !}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine fpdisc(t,n,k2,b,nest)!{{{
c  subroutine fpdisc calculates the discontinuity jumps of the kth
c  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
c  ..scalar arguments..
      integer n,k2,nest
c  ..array arguments..
      double precision t(n),b(nest,k2)
c  ..local scalars..
      double precision an,fac,prod
      integer i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
c  ..local array..
      double precision h(12)
c  ..
      k1 = k2-1
      k = k1-1
      nk1 = n-k1
      nrint = nk1-k
      an = nrint
      fac = an/(t(nk1+1)-t(k1))
      do 40 l=k2,nk1
        lmk = l-k1
        do 10 j=1,k1
          ik = j+k1
          lj = l+j
          lk = lj-k2
          h(j) = t(l)-t(lk)
          h(ik) = t(l)-t(lj)
  10    continue
        lp = lmk
        do 30 j=1,k2
          jk = j
          prod = h(j)
          do 20 i=1,k
            jk = jk+1
            prod = prod*h(jk)*fac
  20      continue
          lk = lp+k1
          b(lmk,j) = (t(lk)-t(lp))/prod
          lp = lp+1
  30    continue
  40  continue
      return
      end !}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine fpgivs(piv,ww,cos,sin)!{{{
c  subroutine fpgivs calculates the parameters of a givens
c  transformation .
c  ..
c  ..scalar arguments..
      double precision piv,ww,cos,sin
c  ..local scalars..
      double precision dd,one,store
c  ..function references..
      double precision abs,sqrt
c  ..
      one = 0.1D+01
      store = abs(piv)
      if(store.ge.ww) dd = store*sqrt(one+(ww/piv)**2)
      if(store.lt.ww) dd = ww*sqrt(one+(piv/ww)**2)
      cos = ww/dd
      sin = piv/dd
      ww = dd
      return
      end !}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine fpintb(t,n,bint,nk1,x,y)!{{{
c  subroutine fpintb calculates integrals of the normalized b-splines
c  nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n.
c  it makes use of the formulae of gaffney for the calculation of
c  indefinite integrals of b-splines.
c
c  calling sequence:
c     call fpintb(t,n,bint,nk1,x,y)
c
c  input parameters:
c    t    : double precision array,length n, containing the position of the knots.
c    n    : integer value, giving the number of knots.
c    nk1  : integer value, giving the number of b-splines of degree k,
c           defined on the set of knots ,i.e. nk1 = n-k-1.
c    x,y  : double precision values, containing the end points of the integration
c           interval.
c  output parameter:
c    bint : array,length nk1, containing the integrals of the b-splines.
c  ..
c  ..scalars arguments..
      integer n,nk1
      double precision x,y
c  ..array arguments..
      double precision t(n),bint(nk1)
c  ..local scalars..
      integer i,ia,ib,it,j,j1,k,k1,l,li,lj,lk,l0,min
      double precision a,ak,arg,b,f,one
c  ..local arrays..
      double precision aint(6),h(6),h1(6)
c  initialization.
      one = 0.1D+01
      k1 = n-nk1
      ak = k1
      k = k1-1
      do 10 i=1,nk1
        bint(i) = 0.0D0
  10  continue
c  the integration limits are arranged in increasing order.
      a = x
      b = y
      min = 0
      if(a-b) 30,160,20
  20  a = y
      b = x
      min = 1
  30  if(a.lt.t(k1)) a = t(k1)
      if(b.gt.t(nk1+1)) b = t(nk1+1)
c  using the expression of gaffney for the indefinite integral of a
c  b-spline we find that
c  bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1)
c    where for t(l) <= x < t(l+1)
c    res(j,x) = 0, j=1,2,...,l-k-1
c             = 1, j=l+1,l+2,...,nk1
c             = aint(j+k-l+1), j=l-k,l-k+1,...,l
c               = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i)))
c                 i=0,1,...,k
      l = k1
      l0 = l+1
c  set arg = a.
      arg = a
      do 90 it=1,2
c  search for the knot interval t(l) <= arg < t(l+1).
  40    if(arg.lt.t(l0) .or. l.eq.nk1) go to 50
        l = l0
        l0 = l+1
        go to 40
c  calculation of aint(j), j=1,2,...,k+1.
c  initialization.
  50    do 55 j=1,k1
          aint(j) = 0.0D0
  55    continue
        aint(1) = (arg-t(l))/(t(l+1)-t(l))
        h1(1) = one
        do 70 j=1,k
c  evaluation of the non-zero b-splines of degree j at arg,i.e.
c    h(i+1) = nl-j+i,j(arg), i=0,1,...,j.
          h(1) = 0.0D0
          do 60 i=1,j
            li = l+i
            lj = li-j
            f = h1(i)/(t(li)-t(lj))
            h(i) = h(i)+f*(t(li)-arg)
            h(i+1) = f*(arg-t(lj))
  60      continue
c  updating of the integrals aint.
          j1 = j+1
          do 70 i=1,j1
            li = l+i
            lj = li-j1
            aint(i) = aint(i)+h(i)*(arg-t(lj))/(t(li)-t(lj))
            h1(i) = h(i)
  70    continue
        if(it.eq.2) go to 100
c  updating of the integrals bint
        lk = l-k
        ia = lk
        do 80 i=1,k1
          bint(lk) = -aint(i)
          lk = lk+1
  80    continue
c  set arg = b.
        arg = b
  90  continue
c  updating of the integrals bint.
 100  lk = l-k
      ib = lk-1
      do 110 i=1,k1
        bint(lk) = bint(lk)+aint(i)
        lk = lk+1
 110  continue
      if(ib.lt.ia) go to 130
      do 120 i=ia,ib
        bint(i) = bint(i)+one
 120  continue
c  the scaling factors are taken into account.
 130  f = one/ak
      do 140 i=1,nk1
        j = i+k1
        bint(i) = bint(i)*(t(j)-t(i))*f
 140  continue
c  the order of the integration limits is taken into account.
      if(min.eq.0) go to 160
      do 150 i=1,nk1
        bint(i) = -bint(i)
 150  continue
 160  return
      end !}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)!{{{
c  subroutine fpknot locates an additional knot for a spline of degree
c  k and adjusts the corresponding parameters,i.e.
c    t     : the position of the knots.
c    n     : the number of knots.
c    nrint : the number of knotintervals.
c    fpint : the sum of squares of residual right hand sides
c            for each knot interval.
c    nrdata: the number of data points inside each knot interval.
c  istart indicates that the smallest data point at which the new knot
c  may be added is x(istart+1)
c  ..
c  ..scalar arguments..
      integer m,n,nrint,nest,istart
c  ..array arguments..
      double precision x(m),t(nest),fpint(nest)
      integer nrdata(nest)
c  ..local scalars..
      double precision an,am,fpmax
      integer ihalf,j,jbegin,jj,jk,jpoint,k,maxbeg,maxpt,
     * next,nrx,number
c  ..
      k = (n-nrint-1)/2
c  search for knot interval t(number+k) <= x <= t(number+k+1) where
c  fpint(number) is maximal on the condition that nrdata(number)
c  not equals zero.
      fpmax = 0.0D0
      jbegin = istart
      do 20 j=1,nrint
        jpoint = nrdata(j)
        if(fpmax.ge.fpint(j) .or. jpoint.eq.0) go to 10
        fpmax = fpint(j)
        number = j
        maxpt = jpoint
        maxbeg = jbegin
  10    jbegin = jbegin+jpoint+1
  20  continue
c  let coincide the new knot t(number+k+1) with a data point x(nrx)
c  inside the old knot interval t(number+k) <= x <= t(number+k+1).
      ihalf = maxpt/2+1
      nrx = maxbeg+ihalf
      next = number+1
      if(next.gt.nrint) go to 40
c  adjust the different parameters.
      do 30 j=next,nrint
        jj = next+nrint-j
        fpint(jj+1) = fpint(jj)
        nrdata(jj+1) = nrdata(jj)
        jk = jj+k
        t(jk+1) = t(jk)
  30  continue
  40  nrdata(number) = ihalf-1
      nrdata(next) = maxpt-ihalf
      am = maxpt
      an = nrdata(number)
      fpint(number) = fpmax*an/am
      an = nrdata(next)
      fpint(next) = fpmax*an/am
      jk = next+k
      t(jk) = x(nrx)
      n = n+1
      nrint = nrint+1
      return
      end !}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      double precision function fprati(p1,f1,p2,f2,p3,f3)!{{{
c  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
c  gives the value of p such that the rational interpolating function
c  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
c  ..
c  ..scalar arguments..
      double precision p1,f1,p2,f2,p3,f3
c  ..local scalars..
      double precision h1,h2,h3,p
c  ..
      if(p3.gt.0.0D0) go to 10
c  value of p in case p3 = infinity.
      p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      go to 20
c  value of p in case p3 ^= infinity.
  10  h1 = f1*(f2-f3)
      h2 = f2*(f3-f1)
      h3 = f3*(f1-f2)
      p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
c  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  20  if(f2.lt.0.0D0) go to 30
      p1 = p2
      f1 = f2
      go to 40
  30  p3 = p2
      f3 = f2
  40  fprati = p
      return
      end !}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine fprota(cos,sin,a,b)!{{{
c  subroutine fprota applies a givens rotation to a and b.
c  ..
c  ..scalar arguments..
      double precision cos,sin,a,b
c ..local scalars..
      double precision stor1,stor2
c  ..
      stor1 = a
      stor2 = b
      b = cos*stor2+sin*stor1
      a = cos*stor1-sin*stor2
      return
      end !}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine splder(t,n,c,k,nu,x,y,m,wrk,ier)!{{{
c  subroutine splder evaluates in a number of points x(i),i=1,2,...,m
c  the derivative of order nu of a spline s(x) of degree k,given in
c  its b-spline representation.
c
c  calling sequence:
c     call splder(t,n,c,k,nu,x,y,m,wrk,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    nu   : integer, specifying the order of the derivative. 0<=nu<=k
c    x    : array,length m, which contains the points where the deriv-
c           ative of s(x) must be evaluated.
c    m    : integer, giving the number of points where the derivative
c           of s(x) must be evaluated
c    wrk  : double precision array of dimension n. used as working space.
c
c  output parameters:
c    y    : array,length m, giving the value of the derivative of s(x)
c           at the different points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    0 <= nu <= k
c    m >= 1
c    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer n,k,nu,m,ier
c  ..array arguments..
      double precision t(n),c(n),x(m),y(m),wrk(n)
c  ..local scalars..
      integer i,j,kk,k1,k2,l,ll,l1,l2,nk1,nk2,nn
      double precision ak,arg,fac,sp,tb,te
c  ..local arrays ..
      double precision h(6)
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nu.lt.0 .or. nu.gt.k) go to 200
      if(m-1) 200,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) go to 200
  20  continue
  30  ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
c  the derivative of order nu of a spline of degree k is a spline of
c  degree k-nu,the b-spline coefficients wrk(i) of which can be found
c  using the recurrence scheme of de boor.
      l = 1
      kk = k
      nn = n
      do 40 i=1,nk1
         wrk(i) = c(i)
  40  continue
      if(nu.eq.0) go to 100
      nk2 = nk1
      do 60 j=1,nu
         ak = kk
         nk2 = nk2-1
         l1 = l
         do 50 i=1,nk2
            l1 = l1+1
            l2 = l1+kk
            fac = t(l2)-t(l1)
            if(fac.le.0.0D0) go to 50
            wrk(i) = ak*(wrk(i+1)-wrk(i))/fac
  50     continue
         l = l+1
         kk = kk-1
  60  continue
      if(kk.ne.0) go to 100
c  if nu=k the derivative is a piecewise constant function
      j = 1
      do 90 i=1,m
         arg = x(i)
  70     if(arg.lt.t(l+1) .or. l.eq.nk1) go to 80
         l = l+1
         j = j+1
         go to 70
  80     y(i) = wrk(j)
  90  continue
      go to 200
 100  l = k1
      l1 = l+1
      k2 = k1-nu
c  main loop for the different points.
      do 180 i=1,m
c  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
c  search for knot interval t(l) <= arg < t(l+1)
 140    if(arg.lt.t(l1) .or. l.eq.nk1) go to 150
        l = l1
        l1 = l+1
        go to 140
c  evaluate the non-zero b-splines of degree k-nu at arg.
 150    call fpbspl(t,n,kk,arg,l,h)
c  find the value of the derivative at x=arg.
        sp = 0.0D0
        ll = l-k1
        do 160 j=1,k2
          ll = ll+1
          sp = sp+wrk(ll)*h(j)
 160    continue
        y(i) = sp
 180  continue
 200  return
      end !}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine splev(t,n,c,k,x,y,m,ier)!{{{
c  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
c  a spline s(x) of degree k, given in its b-spline representation.
c
c  calling sequence:
c     call splev(t,n,c,k,x,y,m,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    x    : array,length m, which contains the points where s(x) must
c           be evaluated.
c    m    : integer, giving the number of points where s(x) must be
c           evaluated.
c
c  output parameter:
c    y    : array,length m, giving the value of s(x) at the different
c           points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    m >= 1
c    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl.
c
c  references :
c    de boor c  : on calculating with b-splines, j. approximation theory
c                 6 (1972) 50-62.
c    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
c                 applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer n,k,m,ier
c  ..array arguments..
      double precision t(n),c(n),x(m),y(m)
c  ..local scalars..
      integer i,j,k1,l,ll,l1,nk1
      double precision arg,sp,tb,te
c  ..local array..
      double precision h(6)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(m-1) 100,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) go to 100
  20  continue
  30  ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
      l = k1
      l1 = l+1
c  main loop for the different points.
      do 80 i=1,m
c  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
c  search for knot interval t(l) <= arg < t(l+1)
  40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50
        l = l1
        l1 = l+1
        go to 40
c  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t,n,k,arg,l,h)
c  find the value of s(x) at x=arg.
        sp = 0.0D0
        ll = l-k1
        do 60 j=1,k1
          ll = ll+1
          sp = sp+c(ll)*h(j)
  60    continue
        y(i) = sp
  80  continue
 100  return
      end !}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      double precision function splint(t,n,c,k,a,b,wrk)!{{{
c  function splint calculates the integral of a spline function s(x)
c  of degree k, which is given in its normalized b-spline representation
c
c  calling sequence:
c     aint = splint(t,n,c,k,a,b,wrk)
c
c  input parameters:
c    t    : array,length n,which contains the position of the knots
c           of s(x).
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, containing the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    a,b  : double precision values, containing the end points of the integration
c           interval. s(x) is considered to be identically zero outside
c           the interval (t(k+1),t(n-k)).
c
c  output parameter:
c    aint : double precision, containing the integral of s(x) between a and b.
c    wrk  : double precision array, length n.  used as working space
c           on output, wrk will contain the integrals of the normalized
c           b-splines defined on the set of knots.
c
c  other subroutines required: fpintb.
c
c  references :
c    gaffney p.w. : the calculation of indefinite integrals of b-splines
c                   j. inst. maths applics 17 (1976) 37-41.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      double precision a,b
      integer n,k
c  ..array arguments..
      double precision t(n),c(n),wrk(n)
c  ..local scalars..
      integer i,nk1
c  ..
      nk1 = n-k-1
c  calculate the integrals wrk(i) of the normalized b-splines
c  ni,k+1(x), i=1,2,...nk1.
      call fpintb(t,n,wrk,nk1,a,b)
c  calculate the integral of s(x).
      splint = 0.0D0
      do 10 i=1,nk1
        splint = splint+c(i)*wrk(i)
  10  continue
      return
      end!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,!{{{
     * ty,ny,p,c,nc,fp,fpx,fpy,mm,mynx,kx1,kx2,ky1,ky2,spx,spy,right,q,
     * ax,ay,bx,by,nrx,nry)
c  ..
c  ..scalar arguments..
      double precision p,fp
      integer ifsx,ifsy,ifbx,ifby,mx,my,mz,kx,ky,nx,ny,nc,mm,mynx,
     * kx1,kx2,ky1,ky2
c  ..array arguments..
      double precision x(mx),y(my),z(mz),tx(nx),ty(ny),c(nc),spx(mx,kx1)
     * ,spy(my,ky1),right(mm),q(mynx),ax(nx,kx2),bx(nx,kx2),ay(ny,ky2),
     * by(ny,ky2),fpx(nx),fpy(ny)
      integer nrx(mx),nry(my)
c  ..local scalars..
      double precision arg,cos,fac,pinv,piv,sin,term,one,half
      integer i,ibandx,ibandy,ic,iq,irot,it,iz,i1,i2,i3,j,k,k1,k2,l,
     * l1,l2,ncof,nk1x,nk1y,nrold,nroldx,nroldy,number,numx,numx1,
     * numy,numy1,n1
c  ..local arrays..
      double precision h(7)
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fprota
c  ..
c  the b-spline coefficients of the smoothing spline are calculated as
c  the least-squares solution of the over-determined linear system of
c  equations  (ay) c (ax)' = q       where
c
c               |   (spx)    |            |   (spy)    |
c        (ax) = | ---------- |     (ay) = | ---------- |
c               | (1/p) (bx) |            | (1/p) (by) |
c
c                                | z  ' 0 |
c                            q = | ------ |
c                                | 0  ' 0 |
c
c  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the
c                b-spline coefficients.
c       z      : the my x mx matrix which contains the function values.
c       spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation
c                matrices according to the least-squares problems in
c                the x- and y-direction.
c       bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1)
c                matrices which contain the discontinuity jumps of the
c                derivatives of the b-splines in the x- and y-direction.
      one = 1
      half = 0.5
      nk1x = nx-kx1
      nk1y = ny-ky1
      if(p.gt.0.) pinv = one/p
c  it depends on the value of the flags ifsx,ifsy,ifbx and ifby and on
c  the value of p whether the matrices (spx),(spy),(bx) and (by) still
c  must be determined.
      if(ifsx.ne.0) go to 50
c  calculate the non-zero elements of the matrix (spx) which is the
c  observation matrix according to the least-squares spline approximat-
c  ion problem in the x-direction.
      l = kx1
      l1 = kx2
      number = 0
      do 40 it=1,mx
        arg = x(it)
  10    if(arg.lt.tx(l1) .or. l.eq.nk1x) go to 20
        l = l1
        l1 = l+1
        number = number+1
        go to 10
  20    call fpbspl(tx,nx,kx,arg,l,h)
        do 30 i=1,kx1
          spx(it,i) = h(i)
  30    continue
        nrx(it) = number
  40  continue
      ifsx = 1
  50  if(ifsy.ne.0) go to 100
c  calculate the non-zero elements of the matrix (spy) which is the
c  observation matrix according to the least-squares spline approximat-
c  ion problem in the y-direction.
      l = ky1
      l1 = ky2
      number = 0
      do 90 it=1,my
        arg = y(it)
  60    if(arg.lt.ty(l1) .or. l.eq.nk1y) go to 70
        l = l1
        l1 = l+1
        number = number+1
        go to 60
  70    call fpbspl(ty,ny,ky,arg,l,h)
        do 80 i=1,ky1
          spy(it,i) = h(i)
  80    continue
        nry(it) = number
  90  continue
      ifsy = 1
 100  if(p.le.0.) go to 120
c  calculate the non-zero elements of the matrix (bx).
      if(ifbx.ne.0 .or. nx.eq.2*kx1) go to 110
      call fpdisc(tx,nx,kx2,bx,nx)
      ifbx = 1
c  calculate the non-zero elements of the matrix (by).
 110  if(ifby.ne.0 .or. ny.eq.2*ky1) go to 120
      call fpdisc(ty,ny,ky2,by,ny)
      ifby = 1
c  reduce the matrix (ax) to upper triangular form (rx) using givens
c  rotations. apply the same transformations to the rows of matrix q
c  to obtain the my x (nx-kx-1) matrix g.
c  store matrix (rx) into (ax) and g into q.
 120  l = my*nk1x
c  initialization.
      do 130 i=1,l
        q(i) = 0.
 130  continue
      do 140 i=1,nk1x
        do 140 j=1,kx2
          ax(i,j) = 0.
 140  continue
      l = 0
      nrold = 0
c  ibandx denotes the bandwidth of the matrices (ax) and (rx).
      ibandx = kx1
      do 270 it=1,mx
        number = nrx(it)
 150    if(nrold.eq.number) go to 180
        if(p.le.0.) go to 260
        ibandx = kx2
c  fetch a new row of matrix (bx).
        n1 = nrold+1
        do 160 j=1,kx2
          h(j) = bx(n1,j)*pinv
 160    continue
c  find the appropriate column of q.
        do 170 j=1,my
          right(j) = 0.
 170    continue
        irot = nrold
        go to 210
c  fetch a new row of matrix (spx).
 180    h(ibandx) = 0.
        do 190 j=1,kx1
          h(j) = spx(it,j)
 190    continue
c  find the appropriate column of q.
        do 200 j=1,my
          l = l+1
          right(j) = z(l)
 200    continue
        irot = number
c  rotate the new row of matrix (ax) into triangle.
 210    do 240 i=1,ibandx
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 240
c  calculate the parameters of the givens transformation.
          call fpgivs(piv,ax(irot,1),cos,sin)
c  apply that transformation to the rows of matrix q.
          iq = (irot-1)*my
          do 220 j=1,my
            iq = iq+1
            call fprota(cos,sin,right(j),q(iq))
 220      continue
c  apply that transformation to the columns of (ax).
          if(i.eq.ibandx) go to 250
          i2 = 1
          i3 = i+1
          do 230 j=i3,ibandx
            i2 = i2+1
            call fprota(cos,sin,h(j),ax(irot,i2))
 230      continue
 240    continue
 250    if(nrold.eq.number) go to 270
 260    nrold = nrold+1
        go to 150
 270  continue
c  reduce the matrix (ay) to upper triangular form (ry) using givens
c  rotations. apply the same transformations to the columns of matrix g
c  to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
c  store matrix (ry) into (ay) and h into c.
      ncof = nk1x*nk1y
c  initialization.
      do 280 i=1,ncof
        c(i) = 0.
 280  continue
      do 290 i=1,nk1y
        do 290 j=1,ky2
          ay(i,j) = 0.
 290  continue
      nrold = 0
c  ibandy denotes the bandwidth of the matrices (ay) and (ry).
      ibandy = ky1
      do 420 it=1,my
        number = nry(it)
 300    if(nrold.eq.number) go to 330
        if(p.le.0.) go to 410
        ibandy = ky2
c  fetch a new row of matrix (by).
        n1 = nrold+1
        do 310 j=1,ky2
          h(j) = by(n1,j)*pinv
 310    continue
c  find the appropiate row of g.
        do 320 j=1,nk1x
          right(j) = 0.
 320    continue
        irot = nrold
        go to 360
c  fetch a new row of matrix (spy)
 330    h(ibandy) = 0.
        do 340 j=1,ky1
          h(j) = spy(it,j)
 340    continue
c  find the appropiate row of g.
        l = it
        do 350 j=1,nk1x
          right(j) = q(l)
          l = l+my
 350    continue
        irot = number
c  rotate the new row of matrix (ay) into triangle.
 360    do 390 i=1,ibandy
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 390
c  calculate the parameters of the givens transformation.
          call fpgivs(piv,ay(irot,1),cos,sin)
c  apply that transformation to the colums of matrix g.
          ic = irot
          do 370 j=1,nk1x
            call fprota(cos,sin,right(j),c(ic))
            ic = ic+nk1y
 370      continue
c  apply that transformation to the columns of matrix (ay).
          if(i.eq.ibandy) go to 400
          i2 = 1
          i3 = i+1
          do 380 j=i3,ibandy
            i2 = i2+1
            call fprota(cos,sin,h(j),ay(irot,i2))
 380      continue
 390    continue
 400    if(nrold.eq.number) go to 420
 410    nrold = nrold+1
        go to 300
 420  continue
c  backward substitution to obtain the b-spline coefficients as the
c  solution of the linear system    (ry) c (rx)' = h.
c  first step: solve the system  (ry) (c1) = h.
      k = 1
      do 450 i=1,nk1x
        call fpback(ay,c(k),nk1y,ibandy,c(k),ny)
        k = k+nk1y
 450  continue
c  second step: solve the system  c (rx)' = (c1).
      k = 0
      do 480 j=1,nk1y
        k = k+1
        l = k
        do 460 i=1,nk1x
          right(i) = c(l)
          l = l+nk1y
 460    continue
        call fpback(ax,right,nk1x,ibandx,right,nx)
        l = k
        do 470 i=1,nk1x
          c(l) = right(i)
          l = l+nk1y
 470    continue
 480  continue
c  calculate the quantities
c    res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my
c    fp = sumi=1,mx(sumj=1,my(res(i,j)))
c    fpx(r) = sum''i(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1
c                  tx(r+kx) <= x(i) <= tx(r+kx+1)
c    fpy(r) = sumi=1,mx(sum''j(res(i,j))) , r=1,2,...,ny-2*ky-1
c                  ty(r+ky) <= y(j) <= ty(r+ky+1)
      fp = 0.
      do 490 i=1,nx
        fpx(i) = 0.
 490  continue
      do 500 i=1,ny
        fpy(i) = 0.
 500  continue
      nk1y = ny-ky1
      iz = 0
      nroldx = 0
c  main loop for the different grid points.
      do 550 i1=1,mx
        numx = nrx(i1)
        numx1 = numx+1
        nroldy = 0
        do 540 i2=1,my
          numy = nry(i2)
          numy1 = numy+1
          iz = iz+1
c  evaluate s(x,y) at the current grid point by making the sum of the
c  cross products of the non-zero b-splines at (x,y), multiplied with
c  the appropiate b-spline coefficients.
          term = 0.
          k1 = numx*nk1y+numy
          do 520 l1=1,kx1
            k2 = k1
            fac = spx(i1,l1)
            do 510 l2=1,ky1
              k2 = k2+1
              term = term+fac*spy(i2,l2)*c(k2)
 510        continue
            k1 = k1+nk1y
 520      continue
c  calculate the squared residual at the current grid point.
          term = (z(iz)-term)**2
c  adjust the different parameters.
          fp = fp+term
          fpx(numx1) = fpx(numx1)+term
          fpy(numy1) = fpy(numy1)+term
          fac = term*half
          if(numy.eq.nroldy) go to 530
          fpy(numy1) = fpy(numy1)-fac
          fpy(numy) = fpy(numy)+fac
 530      nroldy = numy
          if(numx.eq.nroldx) go to 540
          fpx(numx1) = fpx(numx1)-fac
          fpx(numx) = fpx(numx)+fac
 540    continue
        nroldx = numx
 550  continue
      return
      end!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,!{{{
     * nxest,nyest,tol,maxit,nc,nx,tx,ny,ty,c,fp,fp0,fpold,reducx,
     * reducy,fpintx,fpinty,lastdi,nplusx,nplusy,nrx,nry,nrdatx,nrdaty,
     * wrk,lwrk,ier)
c  ..
c  ..scalar arguments..
      double precision xb,xe,yb,ye,s,tol,fp,fp0,fpold,reducx,reducy
      integer iopt,mx,my,mz,kx,ky,nxest,nyest,maxit,nc,nx,ny,lastdi,
     * nplusx,nplusy,lwrk,ier
c  ..array arguments..
      double precision x(mx),y(my),z(mz),tx(nxest),ty(nyest),c(nc),
     * fpintx(nxest),fpinty(nyest),wrk(lwrk)
      integer nrdatx(nxest),nrdaty(nyest),nrx(mx),nry(my)
c  ..local scalars
      double precision acc,fpms,f1,f2,f3,p,p1,p2,p3,rn,one,half,con1,
     * con9,con4
      integer i,ich1,ich3,ifbx,ifby,ifsx,ifsy,iter,j,kx1,kx2,ky1,ky2,
     * k3,l,lax,lay,lbx,lby,lq,lri,lsx,lsy,mk1,mm,mpm,mynx,ncof,
     * nk1x,nk1y,nmaxx,nmaxy,nminx,nminy,nplx,nply,npl1,nrintx,
     * nrinty,nxe,nxk,nye
c  ..function references..
      double precision abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpgrre,fpknot
c  ..
c   set constants
      one = 1
      half = 0.5e0
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
c  we partition the working space.
      kx1 = kx+1
      ky1 = ky+1
      kx2 = kx1+1
      ky2 = ky1+1
      lsx = 1
      lsy = lsx+mx*kx1
      lri = lsy+my*ky1
      mm = max0(nxest,my)
      lq = lri+mm
      mynx = nxest*my
      lax = lq+mynx
      nxk = nxest*kx2
      lbx = lax+nxk
      lay = lbx+nxk
      lby = lay+nyest*ky2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 1: determination of the number of knots and their position.     c
c ****************************************************************     c
c  given a set of knots we compute the least-squares spline sinf(x,y), c
c  and the corresponding sum of squared residuals fp=f(p=inf).         c
c  if iopt=-1  sinf(x,y) is the requested approximation.               c
c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
c    if fp <=s we will continue with the current set of knots.         c
c    if fp > s we will increase the number of knots and compute the    c
c       corresponding least-squares spline until finally fp<=s.        c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmaxx = mx+kx+1  and  nmaxy = my+ky+1.               c
c    if s>0 and                                                        c
c     *iopt=0 we first compute the least-squares polynomial of degree  c
c      kx in x and ky in y; nx=nminx=2*kx+2 and ny=nymin=2*ky+2.       c
c     *iopt=1 we start with the knots found at the last call of the    c
c      routine, except for the case that s > fp0; then we can compute  c
c      the least-squares polynomial directly.                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  determine the number of knots for polynomial approximation.
      nminx = 2*kx1
      nminy = 2*ky1
      if(iopt.lt.0) go to 120
c  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  find nmaxx and nmaxy which denote the number of knots in x- and y-
c  direction in case of spline interpolation.
      nmaxx = mx+kx1
      nmaxy = my+ky1
c  find nxe and nye which denote the maximum number of knots
c  allowed in each direction
      nxe = min0(nmaxx,nxest)
      nye = min0(nmaxy,nyest)
      if(s.gt.0.) go to 100
c  if s = 0, s(x,y) is an interpolating spline.
      nx = nmaxx
      ny = nmaxy
c  test whether the required storage space exceeds the available one.
      if(ny.gt.nyest .or. nx.gt.nxest) go to 420
c  find the position of the interior knots in case of interpolation.
c  the knots in the x-direction.
      mk1 = mx-kx1
      if(mk1.eq.0) go to 60
      k3 = kx/2
      i = kx1+1
      j = k3+2
      if(k3*2.eq.kx) go to 40
      do 30 l=1,mk1
        tx(i) = x(j)
        i = i+1
        j = j+1
  30  continue
      go to 60
  40  do 50 l=1,mk1
        tx(i) = (x(j)+x(j-1))*half
        i = i+1
        j = j+1
  50  continue
c  the knots in the y-direction.
  60  mk1 = my-ky1
      if(mk1.eq.0) go to 120
      k3 = ky/2
      i = ky1+1
      j = k3+2
      if(k3*2.eq.ky) go to 80
      do 70 l=1,mk1
        ty(i) = y(j)
        i = i+1
        j = j+1
  70  continue
      go to 120
  80  do 90 l=1,mk1
        ty(i) = (y(j)+y(j-1))*half
        i = i+1
        j = j+1
  90  continue
      go to 120
c  if s > 0 our initial choice of knots depends on the value of iopt.
 100  if(iopt.eq.0) go to 115
      if(fp0.le.s) go to 115
c  if iopt=1 and fp0 > s we start computing the least- squares spline
c  according to the set of knots found at the last call of the routine.
c  we determine the number of grid coordinates x(i) inside each knot
c  interval (tx(l),tx(l+1)).
      l = kx2
      j = 1
      nrdatx(1) = 0
      mpm = mx-1
      do 105 i=2,mpm
        nrdatx(j) = nrdatx(j)+1
        if(x(i).lt.tx(l)) go to 105
        nrdatx(j) = nrdatx(j)-1
        l = l+1
        j = j+1
        nrdatx(j) = 0
 105  continue
c  we determine the number of grid coordinates y(i) inside each knot
c  interval (ty(l),ty(l+1)).
      l = ky2
      j = 1
      nrdaty(1) = 0
      mpm = my-1
      do 110 i=2,mpm
        nrdaty(j) = nrdaty(j)+1
        if(y(i).lt.ty(l)) go to 110
        nrdaty(j) = nrdaty(j)-1
        l = l+1
        j = j+1
        nrdaty(j) = 0
 110  continue
      go to 120
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  polynomial of degree kx in x and ky in y (which is a spline without
c  interior knots).
 115  nx = nminx
      ny = nminy
      nrdatx(1) = mx-2
      nrdaty(1) = my-2
      lastdi = 0
      nplusx = 0
      nplusy = 0
      fp0 = 0.
      fpold = 0.
      reducx = 0.
      reducy = 0.
 120  mpm = mx+my
      ifsx = 0
      ifsy = 0
      ifbx = 0
      ifby = 0
      p = -one
c  main loop for the different sets of knots.mpm=mx+my is a save upper
c  bound for the number of trials.
      do 250 iter=1,mpm
        if(nx.eq.nminx .and. ny.eq.nminy) ier = -2
c  find nrintx (nrinty) which is the number of knot intervals in the
c  x-direction (y-direction).
        nrintx = nx-nminx+1
        nrinty = ny-nminy+1
c  find ncof, the number of b-spline coefficients for the current set
c  of knots.
        nk1x = nx-kx1
        nk1y = ny-ky1
        ncof = nk1x*nk1y
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(x,y).
        i = nx
        do 130 j=1,kx1
          tx(j) = xb
          tx(i) = xe
          i = i-1
 130    continue
        i = ny
        do 140 j=1,ky1
          ty(j) = yb
          ty(i) = ye
          i = i-1
 140    continue
c  find the least-squares spline sinf(x,y) and calculate for each knot
c  interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=y<=ty(j+ky+1)) the sum
c  of squared residuals fpintx(j),j=1,2,...,nx-2*kx-1 (fpinty(j),j=1,2,
c  ...,ny-2*ky-1) for the data points having their absciss (ordinate)-
c  value belonging to that interval.
c  fp gives the total sum of squared residuals.
        call fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty,
     *  ny,p,c,nc,fp,fpintx,fpinty,mm,mynx,kx1,kx2,ky1,ky2,wrk(lsx),
     *  wrk(lsy),wrk(lri),wrk(lq),wrk(lax),wrk(lay),wrk(lbx),wrk(lby),
     *  nrx,nry)
        if(ier.eq.(-2)) fp0 = fp
c  test whether the least-squares spline is an acceptable solution.
        if(iopt.lt.0) go to 440
        fpms = fp-s
        if(abs(fpms) .lt. acc) go to 440
c  if f(p=inf) < s, we accept the choice of knots.
        if(fpms.lt.0.) go to 300
c  if nx=nmaxx and ny=nmaxy, sinf(x,y) is an interpolating spline.
        if(nx.eq.nmaxx .and. ny.eq.nmaxy) go to 430
c  increase the number of knots.
c  if nx=nxe and ny=nye we cannot further increase the number of knots
c  because of the storage capacity limitation.
        if(nx.eq.nxe .and. ny.eq.nye) go to 420
        ier = 0
c  adjust the parameter reducx or reducy according to the direction
c  in which the last added knots were located.
        if(lastdi) 150,170,160
 150    reducx = fpold-fp
        go to 170
 160    reducy = fpold-fp
c  store the sum of squared residuals for the current set of knots.
 170    fpold = fp
c  find nplx, the number of knots we should add in the x-direction.
        nplx = 1
        if(nx.eq.nminx) go to 180
        npl1 = nplusx*2
        rn = nplusx
        if(reducx.gt.acc) npl1 = rn*fpms/reducx
        nplx = min0(nplusx*2,max0(npl1,nplusx/2,1))
c  find nply, the number of knots we should add in the y-direction.
 180    nply = 1
        if(ny.eq.nminy) go to 190
        npl1 = nplusy*2
        rn = nplusy
        if(reducy.gt.acc) npl1 = rn*fpms/reducy
        nply = min0(nplusy*2,max0(npl1,nplusy/2,1))
 190    if(nplx-nply) 210,200,230
 200    if(lastdi.lt.0) go to 230
 210    if(nx.eq.nxe) go to 230
c  addition in the x-direction.
        lastdi = -1
        nplusx = nplx
        ifsx = 0
        do 220 l=1,nplusx
c  add a new knot in the x-direction
          call fpknot(x,mx,tx,nx,fpintx,nrdatx,nrintx,nxest,1)
c  test whether we cannot further increase the number of knots in the
c  x-direction.
          if(nx.eq.nxe) go to 250
 220    continue
        go to 250
 230    if(ny.eq.nye) go to 210
c  addition in the y-direction.
        lastdi = 1
        nplusy = nply
        ifsy = 0
        do 240 l=1,nplusy
c  add a new knot in the y-direction.
          call fpknot(y,my,ty,ny,fpinty,nrdaty,nrinty,nyest,1)
c  test whether we cannot further increase the number of knots in the
c  y-direction.
          if(ny.eq.nye) go to 250
 240    continue
c  restart the computations with the new set of knots.
 250  continue
c  test whether the least-squares polynomial is a solution of our
c  approximation problem.
 300  if(ier.eq.(-2)) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spline sp(x,y)                c
c *****************************************************                c
c  we have determined the number of knots and their position. we now   c
c  compute the b-spline coefficients of the smoothing spline sp(x,y).  c
c  this smoothing spline varies with the parameter p in such a way thatc
c    f(p) = sumi=1,mx(sumj=1,my((z(i,j)-sp(x(i),y(j)))**2)             c
c  is a continuous, strictly decreasing function of p. moreover the    c
c  least-squares polynomial corresponds to p=0 and the least-squares   c
c  spline to p=infinity. iteratively we then have to determine the     c
c  positive value of p such that f(p)=s. the process which is proposed c
c  here makes use of rational interpolation. f(p) is approximated by a c
c  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
c  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
c  are used to calculate the new value of p such that r(p)=s.          c
c  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = one
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p)=s.
      do 350 iter = 1,maxit
c  find the smoothing spline sp(x,y) and the corresponding sum of
c  squared residuals fp.
        call fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty,
     *  ny,p,c,nc,fp,fpintx,fpinty,mm,mynx,kx1,kx2,ky1,ky2,wrk(lsx),
     *  wrk(lsy),wrk(lri),wrk(lq),wrk(lax),wrk(lay),wrk(lbx),wrk(lby),
     *  nrx,nry)
c  test whether the approximation sp(x,y) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 320
        if((f2-f3).gt.acc) go to 310
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 350
 310    if(f2.lt.0.) ich3 = 1
 320    if(ich1.ne.0) go to 340
        if((f1-f2).gt.acc) go to 330
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 350
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 350
c  test whether the iteration process proceeds as theoretically
c  expected.
 330    if(f2.gt.0.) ich1 = 1
 340    if(f2.ge.f1 .or. f2.le.f3) go to 410
c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 350  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
      fp = 0.
 440  return
      end!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,!{{{
     * nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c given the set of values z(i,j) on the rectangular grid (x(i),y(j)),
c i=1,...,mx;j=1,...,my, subroutine regrid determines a smooth bivar-
c iate spline approximation s(x,y) of degrees kx and ky on the rect-
c angle xb <= x <= xe, yb <= y <= ye.
c if iopt = -1 regrid calculates the least-squares spline according
c to a given set of knots.
c if iopt >= 0 the total numbers nx and ny of these knots and their
c position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
c ally by the routine. the smoothness of s(x,y) is then achieved by
c minimalizing the discontinuity jumps in the derivatives of s(x,y)
c across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
c the amounth of smoothness is determined by the condition that f(p) =
c sum ((z(i,j)-s(x(i),y(j))))**2) be <= s, with s a given non-negative
c constant, called the smoothing factor.
c the fit is given in the b-spline representation (b-spline coefficients
c c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
c uated by means of subroutine bispev.
c
c calling sequence:
c     call regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
c    *  nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer flag. on entry iopt must specify whether a least-
c          squares spline (iopt=-1) or a smoothing spline (iopt=0 or 1)
c          must be determined.
c          if iopt=0 the routine will start with an initial set of knots
c          tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
c          1,...,ky+1. if iopt=1 the routine will continue with the set
c          of knots found at the last call of the routine.
c          attention: a call with iopt=1 must always be immediately pre-
c                     ceded by another call with iopt=1 or iopt=0 and
c                     s.ne.0.
c          unchanged on exit.
c  mx    : integer. on entry mx must specify the number of grid points
c          along the x-axis. mx > kx . unchanged on exit.
c  x     : double precision array of dimension at least (mx). before entry, x(i)
c          must be set to the x-co-ordinate of the i-th grid point
c          along the x-axis, for i=1,2,...,mx. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c  my    : integer. on entry my must specify the number of grid points
c          along the y-axis. my > ky . unchanged on exit.
c  y     : double precision array of dimension at least (my). before entry, y(j)
c          must be set to the y-co-ordinate of the j-th grid point
c          along the y-axis, for j=1,2,...,my. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c  z     : double precision array of dimension at least (mx*my).
c          before entry, z(my*(i-1)+j) must be set to the data value at
c          the grid point (x(i),y(j)) for i=1,...,mx and j=1,...,my.
c          unchanged on exit.
c  xb,xe : double precision values. on entry xb,xe,yb and ye must specify the bound-
c  yb,ye   aries of the rectangular approximation domain.
c          xb<=x(i)<=xe,i=1,...,mx; yb<=y(j)<=ye,j=1,...,my.
c          unchanged on exit.
c  kx,ky : integer values. on entry kx and ky must specify the degrees
c          of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
c          (kx=ky=3) splines. unchanged on exit.
c  s     : double precision. on entry (in case iopt>=0) s must specify the smoothing
c          factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  nxest : integer. unchanged on exit.
c  nyest : integer. unchanged on exit.
c          on entry, nxest and nyest must specify an upper bound for the
c          number of knots required in the x- and y-directions respect.
c          these numbers will also determine the storage space needed by
c          the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
c          in most practical situation nxest = mx/2, nyest=my/2, will
c          be sufficient. always large enough are nxest=mx+kx+1, nyest=
c          my+ky+1, the number of knots needed for interpolation (s=0).
c          see also further comments.
c  nx    : integer.
c          unless ier=10 (in case iopt >=0), nx will contain the total
c          number of knots with respect to the x-variable, of the spline
c          approximation returned. if the computation mode iopt=1 is
c          used, the value of nx should be left unchanged between sub-
c          sequent calls.
c          in case iopt=-1, the value of nx should be specified on entry
c  tx    : double precision array of dimension nmax.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the x-variable, i.e. the position of
c          the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
c          position of the additional knots tx(1)=...=tx(kx+1)=xb and
c          tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
c          if the computation mode iopt=1 is used, the values of tx(1),
c          ...,tx(nx) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tx(kx+2),
c          ...tx(nx-kx-1) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  ny    : integer.
c          unless ier=10 (in case iopt >=0), ny will contain the total
c          number of knots with respect to the y-variable, of the spline
c          approximation returned. if the computation mode iopt=1 is
c          used, the value of ny should be left unchanged between sub-
c          sequent calls.
c          in case iopt=-1, the value of ny should be specified on entry
c  ty    : double precision array of dimension nmax.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the y-variable, i.e. the position of
c          the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
c          position of the additional knots ty(1)=...=ty(ky+1)=yb and
c          ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
c          if the computation mode iopt=1 is used, the values of ty(1),
c          ...,ty(ny) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values ty(ky+2),
c          ...ty(ny-ky-1) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  c     : double precision array of dimension at least (nxest-kx-1)*(nyest-ky-1).
c          on succesful exit, c contains the coefficients of the spline
c          approximation s(x,y)
c  fp    : double precision. unless ier=10, fp contains the sum of squared
c          residuals of the spline approximation returned.
c  wrk   : double precision array of dimension (lwrk). used as workspace.
c          if the computation mode iopt=1 is used the values of wrk(1),
c          ...,wrk(4) should be left unchanged between subsequent calls.
c  lwrk  : integer. on entry lwrk must specify the actual dimension of
c          the array wrk as declared in the calling (sub)program.
c          lwrk must not be too small.
c           lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
c            my*(ky+1) +u
c           where u is the larger of my and nxest.
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c          if the computation mode iopt=1 is used the values of iwrk(1),
c          ...,iwrk(3) should be left unchanged between subsequent calls
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= 3+mx+my+nxest+nyest.
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the spline returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline returned is an interpolating
c            spline (fp=0).
c   ier=-2 : normal return. the spline returned is the least-squares
c            polynomial of degrees kx and ky. in this extreme case fp
c            gives the upper bound for the smoothing factor s.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters nxest and
c            nyest.
c            probably causes : nxest or nyest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the least-squares spline
c            according to the current set of knots. the parameter fp
c            gives the corresponding sum of squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration proces for finding a smoothing spline with
c            fp = s. probably causes : s too small.
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=3  : error. the maximal number of iterations maxit (set to 20
c            by the program) allowed for finding a smoothing spline
c            with fp=s has been reached. probably causes : s too small
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt<=1, 1<=kx,ky<=5, mx>kx, my>ky, nxest>=2*kx+2,
c            nyest>=2*ky+2, kwrk>=3+mx+my+nxest+nyest,
c            lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
c             my*(ky+1) +max(my,nxest),
c            xb<=x(i-1)<x(i)<=xe,i=2,..,mx,yb<=y(j-1)<y(j)<=ye,j=2,..,my
c            if iopt=-1: 2*kx+2<=nx<=min(nxest,mx+kx+1)
c                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
c                        2*ky+2<=ny<=min(nyest,my+ky+1)
c                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
c                    the schoenberg-whitney conditions, i.e. there must
c                    be subset of grid co-ordinates xx(p) and yy(q) such
c                    that   tx(p) < xx(p) < tx(p+kx+1) ,p=1,...,nx-kx-1
c                           ty(q) < yy(q) < ty(q+ky+1) ,q=1,...,ny-ky-1
c            if iopt>=0: s>=0
c                        if s=0 : nxest>=mx+kx+1, nyest>=my+ky+1
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c
c further comments:
c   regrid does not allow individual weighting of the data-values.
c   so, if these were determined to widely different accuracies, then
c   perhaps the general data set routine surfit should rather be used
c   in spite of efficiency.
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the least-squares polynomial (degrees kx,ky) if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the accuracy of the data values.
c   if the user has an idea of the statistical errors on the data, he
c   can also find a proper estimate for s. for, by assuming that, if he
c   specifies the right s, regrid will return a spline s(x,y) which
c   exactly reproduces the function underlying the data he can evaluate
c   the sum((z(i,j)-s(x(i),y(j)))**2) to find a good estimate for this s
c   for example, if he knows that the statistical errors on his z(i,j)-
c   values is not greater than 0.1, he may expect that a good s should
c   have a value not larger than mx*my*(0.1)**2.
c   if nothing is known about the statistical error in z(i,j), s must
c   be determined by trial and error, taking account of the comments
c   above. the best is then to start with a very large value of s (to
c   determine the least-squares polynomial and the corresponding upper
c   bound fp0 for s) and then to progressively decrease the value of s
c   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
c   and more carefully as the approximation shows more detail) to
c   obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if regrid is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   regrid once more with the selected value for s but now with iopt=0.
c   indeed, regrid may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds nxest and
c   nyest. indeed, if at a certain stage in regrid the number of knots
c   in one direction (say nx) has reached the value of its upper bound
c   (nxest), then from that moment on all subsequent knots are added
c   in the other (y) direction. this may indicate that the value of
c   nxest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c   for example, by setting nxest=2*kx+2 (the lowest allowable value for
c   nxest), the user can indicate that he wants an approximation which
c   is a simple polynomial of degree kx in the variable x.
c
c  other subroutines required:
c    fpback,fpbspl,fpregr,fpdisc,fpgivs,fpgrre,fprati,fprota,fpchec,
c    fpknot
c
c  references:
c   dierckx p. : a fast algorithm for smoothing data on a rectangular
c                grid while using spline functions, siam j.numer.anal.
c                19 (1982) 1286-1304.
c   dierckx p. : a fast algorithm for smoothing data on a rectangular
c                grid while using spline functions, report tw53, dept.
c                computer science,k.u.leuven, 1980.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1989
c
c  ..
c  ..scalar arguments..
      double precision xb,xe,yb,ye,s,fp
      integer iopt,mx,my,kx,ky,nxest,nyest,nx,ny,lwrk,kwrk,ier
c  ..array arguments..
      double precision x(mx),y(my),z(mx*my),tx(nxest),ty(nyest),
     * c((nxest-kx-1)*(nyest-ky-1)),wrk(lwrk)
      integer iwrk(kwrk)
c  ..local scalars..
      double precision tol
      integer i,j,jwrk,kndx,kndy,knrx,knry,kwest,kx1,kx2,ky1,ky2,
     * lfpx,lfpy,lwest,lww,maxit,nc,nminx,nminy,mz
c  ..function references..
      integer max0
c  ..subroutine references..
c    fpregr,fpchec
c  ..
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(kx.le.0 .or. kx.gt.5) go to 70
      kx1 = kx+1
      kx2 = kx1+1
      if(ky.le.0 .or. ky.gt.5) go to 70
      ky1 = ky+1
      ky2 = ky1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 70
      nminx = 2*kx1
      if(mx.lt.kx1 .or. nxest.lt.nminx) go to 70
      nminy = 2*ky1
      if(my.lt.ky1 .or. nyest.lt.nminy) go to 70
      mz = mx*my
      nc = (nxest-kx1)*(nyest-ky1)
      lwest = 4+nxest*(my+2*kx2+1)+nyest*(2*ky2+1)+mx*kx1+
     * my*ky1+max0(nxest,my)
      kwest = 3+mx+my+nxest+nyest
      if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 70
      if(xb.gt.x(1) .or. xe.lt.x(mx)) go to 70
      do 10 i=2,mx
        if(x(i-1).ge.x(i)) go to 70
  10  continue
      if(yb.gt.y(1) .or. ye.lt.y(my)) go to 70
      do 20 i=2,my
        if(y(i-1).ge.y(i)) go to 70
  20  continue
      if(iopt.ge.0) go to 50
      if(nx.lt.nminx .or. nx.gt.nxest) go to 70
      j = nx
      do 30 i=1,kx1
        tx(i) = xb
        tx(j) = xe
        j = j-1
  30  continue
      call fpchec(x,mx,tx,nx,kx,ier)
      if(ier.ne.0) go to 70
      if(ny.lt.nminy .or. ny.gt.nyest) go to 70
      j = ny
      do 40 i=1,ky1
        ty(i) = yb
        ty(j) = ye
        j = j-1
  40  continue
      call fpchec(y,my,ty,ny,ky,ier)
      if(ier) 70,60,70
  50  if(s.lt.0.) go to 70
      if(s.eq.0. .and. (nxest.lt.(mx+kx1) .or. nyest.lt.(my+ky1)) )
     * go to 70
      ier = 0
c  we partition the working space and determine the spline approximation
  60  lfpx = 5
      lfpy = lfpx+nxest
      lww = lfpy+nyest
      jwrk = lwrk-4-nxest-nyest
      knrx = 4
      knry = knrx+mx
      kndx = knry+my
      kndy = kndx+nxest
      call fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
     * tol,maxit,nc,nx,tx,ny,ty,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),
     * wrk(lfpx),wrk(lfpy),iwrk(1),iwrk(2),iwrk(3),iwrk(knrx),
     * iwrk(knry),iwrk(kndx),iwrk(kndy),wrk(lww),jwrk,ier)
  70  return
      end!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly)!{{{
c  ..scalar arguments..
      integer nx,ny,kx,ky,mx,my
c  ..array arguments..
      integer lx(mx),ly(my)
      double precision tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),
     * z(mx*my),wx(mx,kx+1),wy(my,ky+1)
c  ..local scalars..
      integer kx1,ky1,l,l1,l2,m,nkx1,nky1
      double precision arg,sp,tb,te
c  ..local arrays..
      double precision h(6)
c  ..subroutine references..
c    fpbspl
c  ..
      kx1 = kx+1
      nkx1 = nx-kx1
      tb = tx(kx1)
      te = tx(nkx1+1)
      l = kx1
      l1 = l+1
      do 40 i=1,mx
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  10    if(arg.lt.tx(l1) .or. l.eq.nkx1) go to 20
        l = l1
        l1 = l+1
        go to 10
  20    call fpbspl(tx,nx,kx,arg,l,h)
        lx(i) = l-kx1
        do 30 j=1,kx1
          wx(i,j) = h(j)
  30    continue
  40  continue
      ky1 = ky+1
      nky1 = ny-ky1
      tb = ty(ky1)
      te = ty(nky1+1)
      l = ky1
      l1 = l+1
      do 80 i=1,my
        arg = y(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  50    if(arg.lt.ty(l1) .or. l.eq.nky1) go to 60
        l = l1
        l1 = l+1
        go to 50
  60    call fpbspl(ty,ny,ky,arg,l,h)
        ly(i) = l-ky1
        do 70 j=1,ky1
          wy(i,j) = h(j)
  70    continue
  80  continue
      m = 0
      do 130 i=1,mx
        l = lx(i)*nky1
        do 90 i1=1,kx1
          h(i1) = wx(i,i1)
  90    continue
        do 120 j=1,my
          l1 = l+ly(j)
          sp = 0.
          do 110 i1=1,kx1
            l2 = l1
            do 100 j1=1,ky1
              l2 = l2+1
              sp = sp+c(l2)*h(i1)*wy(j,j1)
 100        continue
            l1 = l1+nky1
 110      continue
          m = m+1
          z(m) = sp
 120    continue
 130  continue
      return
      end!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,!{{{
     * iwrk,kwrk,ier)
c  subroutine bispev evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
c  ,my a bivariate spline s(x,y) of degrees kx and ky, given in the
c  b-spline representation.
c
c  calling sequence:
c     call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
c    * iwrk,kwrk,ier)
c
c  input parameters:
c   tx    : double precision array, length nx, which contains the position of the
c           knots in the x-direction.
c   nx    : integer, giving the total number of knots in the x-direction
c   ty    : double precision array, length ny, which contains the position of the
c           knots in the y-direction.
c   ny    : integer, giving the total number of knots in the y-direction
c   c     : double precision array, length (nx-kx-1)*(ny-ky-1), which contains the
c           b-spline coefficients.
c   kx,ky : integer values, giving the degrees of the spline.
c   x     : double precision array of dimension (mx).
c           before entry x(i) must be set to the x co-ordinate of the
c           i-th grid point along the x-axis.
c           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
c   mx    : on entry mx must specify the number of grid points along
c           the x-axis. mx >=1.
c   y     : double precision array of dimension (my).
c           before entry y(j) must be set to the y co-ordinate of the
c           j-th grid point along the y-axis.
c           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
c   my    : on entry my must specify the number of grid points along
c           the y-axis. my >=1.
c   wrk   : double precision array of dimension lwrk. used as workspace.
c   lwrk  : integer, specifying the dimension of wrk.
c           lwrk >= mx*(kx+1)+my*(ky+1)
c   iwrk  : integer array of dimension kwrk. used as workspace.
c   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
c
c  output parameters:
c   z     : double precision array of dimension (mx*my).
c           on succesful exit z(my*(i-1)+j) contains the value of s(x,y)
c           at the point (x(i),y(j)),i=1,...,mx;j=1,...,my.
c   ier   : integer error flag
c    ier=0 : normal return
c    ier=10: invalid input data (see restrictions)
c
c  restrictions:
c   mx >=1, my >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
c   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
c   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
c
c  other subroutines required:
c    fpbisp,fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer nx,ny,kx,ky,mx,my,lwrk,kwrk,ier
c  ..array arguments..
      integer iwrk(kwrk)
      double precision tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),
     * z(mx*my),wrk(lwrk)
c  ..local scalars..
      integer i,iw,lwest
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      lwest = (kx+1)*mx+(ky+1)*my
      if(lwrk.lt.lwest) go to 100
      if(kwrk.lt.(mx+my)) go to 100
      if(mx-1) 100,30,10
  10  do 20 i=2,mx
        if(x(i).lt.x(i-1)) go to 100
  20  continue
  30  if(my-1) 100,60,40
  40  do 50 i=2,my
        if(y(i).lt.y(i-1)) go to 100
  50  continue
  60  ier = 0
      iw = mx*(kx+1)+1
      call fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk(1),wrk(iw),
     * iwrk(1),iwrk(mx+1))
 100  return
      end!}}}
!---------------------------------------------------------------------

