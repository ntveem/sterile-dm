      program thermpop

C  corrected 3/21/02 for number density calculation

      integer ngauss,ntemp
      double precision pi,zeta3
      parameter(ngauss=100,ntemp=500)
      parameter(pi=3.141592654d0,zeta3=1.20206d0)

      double precision w(ngauss),absc(ngauss)
      double precision w2(ngauss),absc2(ngauss)

      double precision massmu,masstau,massc,massb  !mu,tau,charm,bottom masses
      double precision masss1,masss2 ! strange quark mass
      
      double precision lowtemp,hightemp,dtemp,temp,lambda,numdens,ke
      double precision edens
      double precision vtch
      double precision mass

      integer i,j

      massmu  = 105.658357d0   !MeV
      masstau = 1777.03d0
      massc   = 1.25d+3
      massb   = 4.20d+3
      masss1  = 70.d0
      masss2  = 175.d0


      call gaulag(absc,w,ngauss,0.d0)
      call gaulag(absc2,w2,ngauss,0.5d0)



C--- Muon Calculation

      mass = massmu

      lowtemp = 1.d0
      hightemp = 3.d0*mass
      dtemp = (hightemp-lowtemp)/(ntemp-1)

      
      open(unit=10,file='thermmu.dat',status='new')
      do i=1,ntemp
         numdens = 0.d0
         edens   = 0.d0
         ke      = 0.d0
         temp   = hightemp-dtemp*(i-1)
         lambda = temp/mass 
         do j=1,ngauss
C            numdens = numdens + (2.d0/(pi**2))*mass**2*temp*w(j)
C     $           *(lambda*absc(j)+1.d0)**2/(exp(1.d0/lambda)+exp(-absc(j
C     $           )))
            numdens = numdens + (2.d0/(pi**2))*lambda**2*mass**3*w(j)
     $           *sqrt(absc(j)**2+2.d0*absc(j)/lambda)*(absc(j)*lambda+1
     $           .d0)/(exp(1.d0/lambda)+exp(-absc(j)))
            edens  = edens + (2.d0/(pi**2))*mass**2*temp**2*w2(j)*(1.d0
     $           +absc2(j)*lambda)**2*dsqrt(absc2(j)+2.d0/lambda)/(exp(1
     $           .d0/lambda)+exp(-absc2(j)))

C********* vtch already has phase space/spin factors (2.d0/(pi**2))  *****
            ke = ke +  mass**2*temp**2*w2(j)*(1.d0
     $           +absc2(j)*lambda)**2*dsqrt(absc2(j)+2.d0/lambda)/(exp(1
     $           .d0/lambda)+exp(-absc2(j)))  
C********* vtch already has phase space/spin factors (2.d0/(pi**2))  *****
        end do
         vtch    = -1.48308d-21*ke
         numdens = numdens/(0.75d0*zeta3/pi**2*4.d0*temp**3)
         write(10,102)temp,numdens,vtch,edens
      end do


      close(unit=10)

 100  format(3(1pe15.7))
 102  format(4(1pe15.7))
 110  format(2(1pe15.7))
C --- Tauon Calculation

      open(unit=10,file='thermtau.dat',status='new')

      mass = masstau

      lowtemp = 10.d0
      hightemp = 3.d0*mass
      dtemp = (hightemp-lowtemp)/(ntemp-1)


      
      open(unit=10,file='thermtau.dat',status='new')
      do i=1,ntemp
         numdens = 0.d0
         ke      = 0.d0
         temp   = hightemp-dtemp*(i-1)
         lambda = temp/mass 
         do j=1,ngauss
            numdens = numdens + (2.d0/(pi**2))*lambda**2*mass**3*w(j)
     $           *sqrt(absc(j)**2+2.d0*absc(j)/lambda)*(absc(j)*lambda+1
     $           .d0)/(exp(1.d0/lambda)+exp(-absc(j)))

C            numdens = numdens + (2.d0/(pi**2))*mass**2*temp*w(j)
C     $           *(lambda*absc(j)+1.d0)**2/(exp(1.d0/lambda)+exp(-absc(j
C     $           )))
C********* vtch already has phase space/spin factors (2.d0/(pi**2))  *****
            ke = ke +  mass**2*temp**2*w2(j)*(1.d0
     $           +absc2(j)*lambda)**2*dsqrt(absc2(j)+2.d0/lambda)/(exp(1
     $           .d0/lambda)+exp(-absc2(j)))  
C********* vtch already has phase space/spin factors (2.d0/(pi**2))  *****
        end do
         vtch    = -1.48308d-21*ke
         numdens = numdens/(0.75d0*zeta3/pi**2*4.d0*temp**3)
         write(10,100)temp,numdens,vtch
      end do


      close(unit=10)

C--- Charm Quark Calculation : min at 50 MeV

      open(unit=10,file='thermc.dat',status='new')
      mass = massc

      lowtemp = 50.d0
      hightemp = 3.d0*mass
      dtemp = (hightemp-lowtemp)/(ntemp-1)


      
      open(unit=10,file='thermc.dat',status='new')
      do i=1,ntemp
         numdens = 0.d0
         temp   = hightemp-dtemp*(i-1)
         lambda = temp/mass 
         do j=1,ngauss
            numdens = numdens + (2.d0/(pi**2))*lambda**2*mass**3*w(j)
     $           *sqrt(absc(j)**2+2.d0*absc(j)/lambda)*(absc(j)*lambda+1
     $           .d0)/(exp(1.d0/lambda)+exp(-absc(j)))

C            numdens = numdens + (2.d0/(pi**2))*mass**2*temp*w(j)
C     $           *(lambda*absc(j)+1.d0)**2/(exp(1.d0/lambda)+exp(-absc(j
C     $           )))
         end do
         numdens = numdens/(0.75d0*zeta3/pi**2*4.d0*temp**3)
         write(10,110)temp,numdens
      end do


      close(unit=10)


C---Bottom Quark Calculation - cut off at 50 MeV

      open(unit=10,file='thermb.dat',status='new')

      mass = massb

      lowtemp = 50.d0
      hightemp = 3.d0*mass
      dtemp = (hightemp-lowtemp)/(ntemp-1)


      
      open(unit=10,file='thermb.dat',status='new')
      do i=1,ntemp
         numdens = 0.d0
         temp   = hightemp-dtemp*(i-1)
         lambda = temp/mass 
         do j=1,ngauss
            numdens = numdens + (2.d0/(pi**2))*lambda**2*mass**3*w(j)
     $           *sqrt(absc(j)**2+2.d0*absc(j)/lambda)*(absc(j)*lambda+1
     $           .d0)/(exp(1.d0/lambda)+exp(-absc(j)))
C            numdens = numdens + (2.d0/(pi**2))*mass**2*temp*w(j)
C     $           *(lambda*absc(j)+1.d0)**2/(exp(1.d0/lambda)+exp(-absc(j
C     $           )))
         end do
         numdens = numdens/(0.75d0*zeta3/pi**2*4.d0*temp**3)
         write(10,110)temp,numdens
      end do
      close(unit=10)

      end !program thermpop


      SUBROUTINE gaulag(x,w,n,alf)
      INTEGER n,MAXIT
      DOUBLE PRECISION alf,w(n),x(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.D-14,MAXIT=10)
CU    USES gammln
      INTEGER i,its,j
      DOUBLE PRECISION ai,gammln
      DOUBLE PRECISION p1,p2,p3,pp,z,z1
      do 13 i=1,n
        if(i.eq.1)then
          z=(1.d0+alf)*(3.d0+.92d0*alf)/(1.d0+2.4d0*n+1.8d0*alf)
        else if(i.eq.2)then
          z=z+(15.d0+6.25d0*alf)/(1.d0+.9d0*alf+2.5d0*n)
        else
          ai=i-2
          z=z+((1.d0+2.55d0*ai)/(1.9d0*ai)+1.26d0*ai*alf/
     $(1.d0+3.5d0*ai))*(z-x(i-2))/(1.d0+.3d0*alf)
        endif
        do 12 its=1,MAXIT
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
11        continue
          pp=(n*p1-(n+alf)*p2)/z
          z1=z
          z=z1-p1/pp
          if(dabs(z-z1).le.EPS)goto 1
12      continue
        pause 'too many iterations in gaulag'
1       x(i)=z
        w(i)=-dexp(gammln(alf+n)-gammln(dble(float(n))))/(pp*n*p2)
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.

      FUNCTION gammln(xx)
      DOUBLE PRECISION gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+dlog(stp*ser/x)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.
