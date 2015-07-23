C Kinder Interfaces to the FitPack spline functions, and some old 
C functions from Numerical Recipes used by the previous code
!******************************************************************
      SUBROUTINE SplFit(x,y,m,t,c,n)
!******************************************************************!{{{
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
!******************************************************************!}}}
      END SUBROUTINE ! SplFit
!******************************************************************

!******************************************************************
      SUBROUTINE TwoDimSplFit(x,y,z,mx,my,tx,ty,c,nx,ny)
!******************************************************************!{{{
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

      Parameter (iopt=0,s=0.0E0)

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
!******************************************************************!}}}
      END SUBROUTINE ! TwoDimSplFit
!******************************************************************

!******************************************************************
      DOUBLE PRECISION FUNCTION NIntegrate(x,y,m,a,b)
!******************************************************************!{{{
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
!******************************************************************!}}}
      END FUNCTION ! NIntegrate
!******************************************************************

!******************************************************************
      SUBROUTINE polint(xa,ya,n,x,y,dy)
!******************************************************************!{{{
      Implicit none!{{{
      Integer n,NMAX
      Double precision dy,x,y,xa(n),ya(n)
      Parameter (NMAX=10)
      Integer i,m,ns
      Double precision den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)!}}}
      ns=1
      dif=dabs(x-xa(1))
      do 11 i=1,n
        dift=dabs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
 11   continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.D0)then
            print *,ho,hp,xa,ya,n,x,y
            stop 'failure in polint'
          end if
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
 12     continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
 13   continue
      return
!******************************************************************!}}}
      END ! SUBROUTINE polint
!     (C) Copr. 1986-92 Numerical Recipes Software V,3.
!******************************************************************

!******************************************************************
      SUBROUTINE hunt(xx,n,x,jlo)
!******************************************************************!{{{
!     Takes in monotonic array xx of length n, and finds the index jlo
!      of the element Closest to x
      Integer jlo,n!{{{
      Double precision x,xx(n)
      Integer inc,jhi,jm
      LOGICAL ascnd!}}}
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
 1      jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
 2      jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
 3    if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
!******************************************************************!}}}
      END ! SUBROUTINE hunt
!     (C) Copr. 1986-92 Numerical Recipes Software V,3.
!******************************************************************

! To test this file     
!******************************************************************
      !Program main
!!******************************************************************!{{{

      !Implicit None
      !!Integer i
      !!Parameter (i=4002)
      !Integer mx, my
      !Parameter (mx=100,my=50)

      !!Double precision a(i), b(i)
      !!Double precision knot(i+4), coeff(i+4)
      !Double precision x(mx), y(my), z(mx*my)
      !Double precision tx(mx+4), ty(my+4), c(mx*my)
      !Double precision xeval, yeval, zeval
      !Double precision wrk(mx*4+my*4)
      !Integer lwrk, iwrk(mx+my), kwrk

      !Parameter (lwrk=mx*4+my*4, kwrk=mx+my)

      !!Integer j, n
      !Integer m, n, nx, ny, ier

      !!Do j=1, i
         !!a(j) = dble(j)
         !!b(j) = dble(j**2)
      !!End do
      !Do m=1, mx
         !x(m) = m
      !End do

      !Do n=1, my
         !y(n) = n
      !End do

      !Do m=1, mx
         !Do n=1, my
            !z((m-1)*my+n) = m*n
         !End do
      !End do            

      !!Call SplFit(a, b, i, knot, coeff, n)
      !Call TwoDimSplFit(x,y,z,mx,my,tx,ty,c,nx,ny)

      !!Do j=1, n
         !!Print *,knot(j), coeff(j)
      !!End do

      !xeval = 44.8D0
      !yeval = 34.65D0
      !Call bispev(tx,nx,ty,ny,c,3,3,xeval,1,yeval,1,zeval,wrk,lwrk,
      !* iwrk,kwrk,ier)

      !Print *,zeval

      !Stop
!!******************************************************************!}}}
      !End program !main
!******************************************************************
