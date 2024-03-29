c
c***************************************************************************
c


       subroutine maxwell_dist(u0, NUPAR, NUPER,
     .                       UminPara, UmaxPara,
     .                       UPERP, UPARA, DFDUPER, DFDUPAR)
       implicit none

       integer, intent(in) ::NUPAR,NUPER
       integer n,m

       real :: UminPara,UmaxPara
       real :: UPERP(NUPER)
       real :: UPARA(NUPAR)
       real, intent(inout) :: DFDUPER(NUPER,NUPAR), DFDUPAR(NUPER,NUPAR)
       real  pi, fnorm, f, u2, u0, u02

       parameter (PI = 3.141592653597932384)

       fnorm = u0**3 / pi**1.5
       u02 = u0**2

       do n = 1, NUPER
          do m = 1, NUPAR
             u2 = UPERP(n)**2 + UPARA(m)**2
             f = fnorm * exp(-u2 * u02)
             DFDUPER(n, m) = f *(-2.0 * UPERP(n) * u02)
             DFDUPAR(n, m) = f *(-2.0 * UPARA(m) * u02)
          end do
       end do


       end subroutine maxwell_dist

c
c***************************************************************************
c


      subroutine besjc (z,n,b,ier)
c
c subroutine besjc (z,n,b,ier)
c
c dimension of           b(n+1)
c arguments
c
c latest revision        october 1978
c
c purpose                to calculate j bessel functions for complex
c                        argument and integer order
c
c usage                  call besjc (z,n,b,ier)
c
c arguments
c
c on input               z
c                          complex argument for which j bessel functions
c                          are to be calculated.  abs(aimag(z)) must be
c                          less than the largest real argument that the
c                          fortran function exp can handle.
c
c                        n
c                          integer, the highest order to be calculated.
c                          n must be greater than or equal to zero.
c
c on output              b
c                          complex vector of length n+1 containing the
c                          bessel function values j-sub-0(z),j-sub-1(z),
c                          ...,j-sub-n(z) in b(1),b(2),...,b(n+1).
c                        ier
c                          an integer error flag.
c                          =0 if all desired orders have been calculated
c                             satisfactorily,
c                          =1 if abs(aimag(z)) is too large,
c                          =2 if n is less than zero,
c                          =2+k if only the first k results are correct.
c                             in the returned values b(m) for m greater
c                             than k, approximately the last
c                             alog10(abs(b(m)/b(k))) significant digits
c                             are in error.
c
c entry points           besjc
c
c special conditions     none
c
c common blocks          none
c
c i/o                    none, except for error messages produced by
c                        calling the error handling routine uliber.
c
c precision              single
c
c specialist             russ rew, ncar, boulder, colorado
c
c accuracy               in tests run on the 7600 with orders from 0
c                        through 10 using random values of the argument
c                        with absolute value less than 50 but
c                        concentrated around the origin, the maximum
c                        relative error (or absolute error when it was
c                        larger) observed was about 8.2e-14.
c
c timing                 on ncar"s control data 7600, besic takes about
c                        .32+.008*n milliseconds when z=(1.0,1.0).
c
c portability            ansi 1966 standard
c
c
c
c
c
c
      complex         z          ,b(1)
      data iorj/0/,xlarge/1000000./
c
      nb = n+1
      call b2slci (real(z),aimag(z),nb,iorj,b,ncalc)
      ier = 0
      if (ncalc .eq. nb) go to 103
      if (ncalc .ge. 0) go to 102
      if (n .ge. 0) go to 101
      ier = 2
      call uliber (ier,25h in besjc, n out of range,25)
      go to 103
  101 ier = 1
      call uliber (ier,25h in besjc, x out of range,25)
      go to 103
  102 ier = 2+ncalc
      call uliber (ier,40h in besjc, accuracy lost for some orders,40)
  103 return
      end


c
c********************************************************************
c

      subroutine deriv_x(f, id, jd, kd, i, j, k, imax, jmax, kmax,
     .                     dx, dfdx, d2fdx2)

      implicit none

      integer id, jd, kd, i, j, k, imax, jmax, kmax
      real f(id, jd, kd), dx, dfdx, d2fdx2

      if(i .ne. 1 .and. i .ne. imax)then
         dfdx = (f(i+1, j, k) - f(i-1, j, k)) / (2.0 * dx)
         d2fdx2 = (f(i+1,j,k) - 2. * f(i,j,k) + f(i-1,j,k)) / dx**2
      end if

      if(i .eq. 1)then
         dfdx = (f(i+1, j, k) - f(i, j, k)) / dx
         d2fdx2 = (f(i+2,j,k) - 2. * f(i+1,j,k) + f(i,j,k)) / dx**2
      end if

      if(i .eq. imax)then
         dfdx = (f(i, j, k) - f(i-1, j, k)) / dx
         d2fdx2 = (f(i,j,k) - 2. * f(i-1,j,k) + f(i-2,j,k)) / dx**2
      end if

      return
      end

c
c***************************************************************************
c
      subroutine deriv_y(f, id, jd, kd, i, j, k, imax, jmax, kmax,
     .                     dy, dfdy, d2fdy2)

      implicit none

      integer id, jd, kd, i, j, k, imax, jmax, kmax
      real f(id, jd, kd), dy, dfdy, d2fdy2

      if(j .ne. 1 .and. j .ne. jmax)then
         dfdy = (f(i, j+1, k) - f(i, j-1, k)) / (2.0 * dy)
         d2fdy2 = (f(i,j+1,k) - 2. * f(i,j,k) + f(i,j-1,k)) / dy**2
      end if

      if(j .eq. 1)then
         dfdy = (f(i, j+1, k) - f(i, j, k)) / dy
         d2fdy2 = (f(i,j+2,k) - 2. * f(i,j+1,k) + f(i,j,k)) / dy**2
      end if

      if(j .eq. jmax)then
         dfdy = (f(i, j, k) - f(i, j-1, k)) / dy
         d2fdy2 = (f(i,j,k) - 2. * f(i,j-1,k) + f(i,j-2,k)) / dy**2
      end if

      return
      end
c
c***************************************************************************
c
      subroutine deriv_xy(f, id, jd, kd, i, j, k, imax, jmax, kmax,
     .                     dx, dy, d2fdxy)

      implicit none

      integer id, jd, kd, i, j, k, imax, jmax, kmax
      real f(id, jd, kd), dx, dy, d2fdxy

      d2fdxy = 0.0

      if(i .ne. 1 .and. i .ne. imax .and.
     .   j .ne. 1 .and. j .ne. jmax) then

         d2fdxy = (f(i+1, j+1, k) - f(i-1, j+1, k)
     .           - f(i+1, j-1, k) + f(i-1, j-1, k) ) / (4. * dx * dy)
      end if


      return
      end

c
c***************************************************************************
c
      subroutine deriv_phi(f, id, jd, kd, i, j, k, imax, jmax, kmax,
     .                     dphi, dfdphi, d2fdphi2)

      implicit none

      integer id, jd, kd, i, j, k, imax, jmax, kmax, kp1, km1
      real f(id, jd, kd), dphi, dfdphi, d2fdphi2

      kp1 = k + 1
      km1 = k - 1

      if(k .eq. 1)    km1 = kmax
      if(k .eq. kmax) kp1 = 1

      dfdphi = (f(i, j, kp1) - f(i, j, km1)) / (2.0 * dphi)
      d2fdphi2 = (f(i,j,kp1) - 2.0 * f(i,j,k) + f(i,j,km1)) / dphi**2

      return
      end

c
c***************************************************************************
c

      subroutine deriv_xphi(f, id, jd, kd, i, j, k, imax, jmax, kmax,
     .                     dx, dphi, d2fdxphi)

      implicit none

      integer id, jd, kd, i, j, k, imax, jmax, kmax, kp1, km1
      real f(id, jd, kd), dx, dphi, d2fdxphi

      kp1 = k + 1
      km1 = k - 1

      if(k .eq. 1)    km1 = kmax
      if(k .eq. kmax) kp1 = 1

      if (i .ne. 1 .and. i .ne. imax)then

         d2fdxphi = (f(i+1, j, kp1) - f(i+1, j, km1)
     .             - f(i-1, j, kp1) + f(i-1, j, km1)) / (4. * dx * dphi)


      end if

      return
      end

c
c***************************************************************************
c
      subroutine deriv_yphi(f, id, jd, kd, i, j, k, imax, jmax, kmax,
     .                     dy, dphi, d2fdyphi)

      implicit none

      integer id, jd, kd, i, j, k, imax, jmax, kmax, kp1, km1
      real f(id, jd, kd), dy, dphi, d2fdyphi

      kp1 = k + 1
      km1 = k - 1

      if(k .eq. 1)    km1 = kmax
      if(k .eq. kmax) kp1 = 1

      if (j .ne. 1 .and. j .ne. jmax)then

         d2fdyphi = (f(i, j+1, kp1) - f(i, j+1, km1)
     .             - f(i, j-1, kp1) + f(i, j-1, km1)) / (4. * dy * dphi)


      end if

      return
      end
c
c***************************************************************************
c
      subroutine torint(x, y, phi, f,  ans,
     .   nx1, nx2, ny1, ny2, nphi1, nphi2,
     .   capr, nxmax, nymax, nphimax)

      implicit none

      integer i, j, k, nx1, nx2, ny1, ny2, nphi1, nphi2,
     .   nxmax, nymax, nphimax

      real f(nxmax, nxmax, nphimax),
     .     x(nxmax), y(nymax), phi(nphimax),
     .     capr(nxmax), ans, dx, dy, dphi

      ans = 0.0

      do k = nphi1, nphi2 - 1
         dphi = phi(k+1) - phi(k)

         do i = nx1, nx2 - 1
            dx = x(i+1) - x(i)

            do j = ny1, ny2 - 1
               dy = y(j+1) - y(j)

               ans = ans + dx * dy * dphi * (
     .                      capr(i) * (f(i,j,k)   + f(i,j+1,k)
     .                               + f(i,j,k+1) + f(i,j+1,k+1))

     .                  + capr(i+1) * (f(i+1,j,k)   + f(i+1,j+1,k)
     .                               + f(i+1,j,k+1) + f(i+1,j+1,k+1))
     .                                       )/8.0
            end do
         end do
      end do



      return
      end

c
c***************************************************************************
c


      SUBROUTINE ZFUN(Z,FU)
C
C     ROUTINE WHICH EVALUATES THE PLASMA DISPERSION FUNCTION.  USES
C     NUMERICAL INTEGRATION (ABSOLUTE VALUE OF THE COMPLEX ARGUMENT, Z,
C     LESS THAN 5) OR ASYMPTOTIC EXPANSION.
C
      COMPLEX Z,FU,TEMP1,TEMP2,Z2,TPIIOD
      DIMENSION C(21),D(21),E(21),F(21),W(21)
      DATA DELTA/5.0E-1/, YI/-1.0E+0/, N/21/, ITEST/0/
      DATA ZERO/0.0E+0/, OSQPI/5.6418958355E-1/, TPI/6.28318530718E+0/
C
      IF(ITEST .EQ. 1) GO TO 1
      ITEST=1
C***     DEFINE WEIGHTS AND STORE CONSTANTS USED FOR INTEGRATION.
C***     WEIGHTS ARE DERIVED FROM A 3 POINT INTEGRATION SCHEME.
      NM3=N-3
      CONOI=DELTA*OSQPI
      W(1)=CONOI*3.75E-1
      W(2)=CONOI*1.1666666667E+0
      W(3)=CONOI*9.5833333333E-1
      DO 20 I=1,3
        NMIP1=N-I+1
   20 W(NMIP1)=W(I)
      DO 30 I=4,NM3
   30 W(I)=CONOI
      NO2=N/2
      NO2P1=NO2+1
      X=ZERO
      Y=YI
      E(NO2P1)=X
      F(NO2P1)=Y
      TEMP1=CMPLX(X,Y)
      TEMP2=CEXP(-TEMP1*TEMP1)
      C(NO2P1)=REAL(TEMP2)*W(NO2P1)
      D(NO2P1)=AIMAG(TEMP2)*W(NO2P1)
      DO 200 I=1,NO2
      X=DELTA*FLOAT(I)
      NPI=NO2P1+I
      NMI=NO2P1-I
      TEMP1=CMPLX(X,Y)
      TEMP2=CEXP(-TEMP1*TEMP1)
      C(NPI)=REAL(TEMP2)*W(NPI)
      C(NMI)=REAL(TEMP2)*W(NMI)
      D(NPI)=AIMAG(TEMP2)*W(NPI)
      D(NMI)=AIMAG(TEMP2)*W(NMI)*(-1.0E+0)
      E(NMI)=-X
      E(NPI)=X
      F(NPI)=Y
      F(NMI)=Y
  200 CONTINUE
      TPIOD=TPI/DELTA
      TPIIOD=CMPLX(ZERO,TPIOD)
C***     BEGIN CALCULATIONS.
    1 G=REAL(Z)
      YY=AIMAG(Z)
      H=ABS(YY)
      ZA=Z*CONJG(Z)
      IF(ZA .GE. 2.5E+1) GO TO 5
      Z2=Z*Z
C***     NUMERICAL INTEGRATION.
C***     F=1/SQRT(PI)*SUM OF...W(I)*EXP(-X(I)**2)/(X(I)-Z)...I=1,N.
C***     INTEGRATION IS ALONG A LINE X(I) IN THE COMPLEX PLANE, WHERE
C***     THE IMAGINARY PART OF X(I)=YI AND THE DIFFERENCE BETWEEN
C***     SUCCESSIVE REAL PARTS OF X(I)=DELTA.  LIMITS OF INTEGRATION
C***     ARE FROM -DELTA*N/2 TO DELTA*N/2.
C***     COMPUTE THE INTEGRAL BY TAKING THE SUM FROM 1 TO N OF THE
C***     CONSTANTS DIVIDED BY X(I)-Z.  USES REAL ARITHMETIC.
      ZR=0.0E+0
      ZI=0.0E+0
      DO 7 I=1,N
      A=E(I)-G
      B=F(I)-H
      DEN=A*A+B*B
      ODEN=1.0E+0/DEN
      ZR=ZR+(A*C(I)+B*D(I))*ODEN
      ZI=ZI+(A*D(I)-B*C(I))*ODEN
    7 CONTINUE
C***     ADD THE CORRECTION TERM.
      FU=CMPLX(ZR,ZI)+(0.0E+0,-3.5449077018E+0)*CEXP(-Z2-TPIOD*
     *   (H-YI)+TPIIOD*G)
      IF(YY .GE. ZERO) GO TO 6
C***     IMAGINARY PART OF ARGUMENT IS NEGATIVE.
      FU=CONJG(FU)+(0.0E+0,3.5449077018E+0)*CEXP(-Z2)
      GO TO 6
C***     MAGNITUDE OF ARGUMENT IS GREATER THAN 5, USE
C***     ASYMPTOTIC EXPANSION.
    5 CALL AEXPAN(Z,G,H,YY,FU)
    6 RETURN
      END
      SUBROUTINE AEXPAN(Z,G,H,YY,FU)
C
C     ROUTINE WHICH COMPUTES THE PLASMA DISPERSION FUNCTION USING
C     ASYMPTOTIC EXPANSION.  IF THE IMAGINARY PART OF THE ARGUMENT, YY,
C     IS EQUAL TO ZERO REAL ARITHMETIC IS USED.
C
      COMPLEX Z,FU,A,Z2,OZ2
      DATA ZERO/0.0E+0/, N/8/
C
      IF(YY .EQ. ZERO) GO TO 1
C***     COMPLEX ARITHMETIC.
      Z2=Z*Z
      OZ2=1.0E+0/Z2
      FU=-1.0E+0/Z
      A=FU
      EN=5.0E-1
      DO 10 I=1,N
      A=EN*A*OZ2
      FU=FU+A
   10 EN=EN+1.0E+0
      IF(YY .GT. ZERO) GO TO 30
      IF(H .GT. SQRT(G*G+1.72E+2)) GO TO 20
      FU=FU+(0.0E+0,3.5449077018E+0)*CEXP(-Z2)
      GO TO 30
C***     ERROR STOP TO AVOID OVERFLOW.
   20 WRITE (51,100) Z
  100 FORMAT(//1X,'*** ERROR STOP IN FRDCNT ROUTINE, ARGUMENT IS',
     * ' TOO SMALL,  ARG =',1P2E14.6)
      STOP
C***     REAL ARITHMETIC.
    1 X2=G*G
      OX2=1.0E+0/X2
      F=-1.0E+0/G
      B=F
      EN=5.0E-1
      DO 11 I=1,N
      B=EN*B*OX2
      F=F+B
   11 EN=EN+1.0E+0
      C=1.7724538509E+0*EXP(-X2)
      FU=CMPLX(F,C)
   30 RETURN
      END
C
C
C TWO POL APPROXIMATION
C
      SUBROUTINE Z2P(Z,F)
      COMPLEX Z,F
      F=CMPLX(.5,.81)/(CMPLX(.51,-.81)-Z)-CMPLX(.5,-.81)
     1/(CMPLX(.51,.81)+Z)
      RETURN
      END

c
c THIS CODE IS SUPERIOR TO PDF AND PDFLL.  NO ERROR AT ZETA=(4.,0.1)
c  disp(z) is the plasma z function
c
      complex function disp(z)
      complex z
      real im
      x=real(z)
      y=aimag(z)
      call zzdisp(x,y,re,im)
      disp=cmplx(re,im)
      return
      end
c
c***************************************************************************
c
c
c  disp1(z) calculates first derivative of disp(z)
c
      complex function disp1(z)
      complex z,dp,disp
      complex zz
      zz=z
      dp=disp(z)
      disp1=-2.0-2.0*zz*dp
      return
      end
c
c***************************************************************************
c
c
c  calculates plasma z function
c
      subroutine zzdisp(x,y,zzr,zzi)
      x1=abs(x)
      y1=abs(y)
      call wzdisp(x1,y1,wzr1,wzi1)
      if(x.ge.0..and.y.ge.0.) go to 1
      if(x.le.0..and.y.ge.0.) go to 2
      a=2*x1*y1
      b=-(x1*x1-y1*y1)
      abr=2*exp(b)*cos(a)
      abi=-2*exp(b)*sin(a)
      wzr1=abr-wzr1
      wzi1=abi-wzi1
      if(x.le.0..and.y.le.0.) go to 1
    2 wzi1=-wzi1
    1 zzr=-1.7724538509055*wzi1
      zzi=1.7724538509055*wzr1
      return
      end
c
c***************************************************************************
c
cc   needed for z function
cc
      subroutine wzdisp(x,y,re,im)
      equivalence(epsh,epsl,epsy)
      real im,lambda
      integer capn
      logical b
      epsh=1.e-12
      if(y.lt.4.29 .and. x.lt.5.33) go to 10
c
c  (x,y) belongs to q1-r
c
      h=0.
      capn=0
      nu=8
      go to 20
c
c  (x,y) belongs to r
c
   10 s=(1.-y/4.29)*sqrt(1.-x*x/28.41)
      h=1.6*s
      h2=2.*h
      capn=6.+23*s+.5
      nu=9.+21*s+.5
      lambda=h2**capn
   20 b=(h.eq.0. .or. lambda.lt.epsl)
c
c  statement (lambda.lt.epsl) covers the underflow case
c  when h(.gt.0) is very small.
c
      rr=0.
      ri=0.
      sr=0.
      si=0.
      nup=nu+1
      do 100 i=1,nup
      n=nup-i
      np1=n+1
      tr=y+h+np1*rr
      ti=x-np1*ri
      c=.5/(tr*tr+ti*ti)
      rr=c*tr
      ri=c*ti
      if(.not.(h.gt.0..and.n.le.capn)) go to 100
      tr=lambda+sr
      sr=rr*tr-ri*si
      si=ri*tr+rr*si
      lambda=lambda/h2
  100 continue
      cc=1.12837916709551
      if(y.lt.epsy) go to 120
      if(b) go to 110
      re=sr*cc
      go to 130
  110 re=rr*cc
      go to 130
  120 re=exp(-x*x)
  130 if(b) go to 140
      im=si*cc
      go to 150
  140 im=ri*cc
  150 return
      end

c
c***************************************************************************
c
      real function second1(dummy)

      implicit none

      integer mtime, mclock
      real dummy

c*****FALCON:
      double precision mpi_wtime
      external mpi_wtime
      second1 = MPI_WTIME()


c*****EAGLE:
c      mtime = mclock()
c      second1  = 0.01 * mtime


      return
      end
