c
c***************************************************************************
c
       subroutine sigmah_stix_elect(i, j, n, m,
     .   bmod, gradprlb,
     .   xm, q, xn, xnuomg,
     .   xkt, omgc, omgp2,
     .   lmin, lmax, nzfun, ibessel,
     .   xkxsav, xkysav, xkzsav, capr,
     .   bx,    by,    bz,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   sigxx, sigxy, sigxz,
     .   sigyx, sigyy, sigyz,
     .   sigzx, sigzy, sigzz,
     .   iflag_gammab, delta0, xk0, damping, xkx_cutoff, xky_cutoff)

*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     and uses the old left handed U matrix in paper
*     ---------------------------------------------------------

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, iflag_gammab

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, sgn_kprl
      real dakbdkb, xnuomg, xkprl_cutoff
      real akprl, gammab, rho, alpha, eps0, omgrf, v0i
      real a, b, bmod, gradprlb
      real bx, by, bz
      real xkxsav, xkysav, xkzsav, capr
      real xkphi
      real xkalp, xkbet, xk0, rgamma

      real delta0, damping, xkx_cutoff, xky_cutoff, kr, step

      real uxx, uxy, uxz,
     .     uyx, uyy, uyz,
     .     uzx, uzy, uzz

      complex zi, zfunct, fzeta, omgrfc
      complex zfunct0, zeta0, sig3cold, z0, z1, z2
      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l


      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz

      parameter (lmaxdim = 99)

      complex xil(0: lmaxdim), xilp(0: lmaxdim)
      complex exil(0: lmaxdim), exilp(0: lmaxdim), exilpp(0: lmaxdim)

      complex zetal, zieps0, arg,
     .   al, bl, cl,
     .   gamma, zeta_eff

      complex dpldkb, dplpdkb, dbldkb
      complex dpldakb, dplpdakb, dbldakb


      common/sigcom/zi, eps0, v0i, omgrf



      xme = 9.11e-31
      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rho = alpha / omgc
      xkphi = xkzsav / capr
      omgrfc = omgrf * (1. + zi * xnuomg)

      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi

      xkperp = sqrt(xkalp**2 + xkbet**2)
      akprl = abs(xkprl)



      sgn_kprl = sign(1.0, xkprl)

      if (akprl .lt. 0.1) then

         xkprl = 0.1 * sgn_kprl
         akprl = abs(xkprl)

      end if


      gamma = 0.5 * xkperp**2 * rho**2
      rgamma = real(gamma)


      call besiexp(gamma, lmax, exil, exilp, exilpp, lmaxdim)

      sig0 = 0.0
      sig1 = 0.0
      sig2 = 0.0
      sig3 = 0.0
      sig4 = 0.0
      sig5 = 0.0


      do l = - 1, 1
         labs = abs(l)



         zetal = (omgrfc - l * omgc) / (xkprl * alpha)

*        --------------------------------
*        Set gammab to zero for electrons;
*        --------------------------------
         gammab = 0.0


         if(abs(gammab) .gt. 1000.0) gammab = 1000.0
         if(abs(gammab) .lt. .01)gammab = .01



         if (nzfun .eq. 0) call z_approx(sgn_kprl, zetal, 0.0,
     .                                                     z0, z1, z2)
         if (nzfun .eq. 1) call z_approx(sgn_kprl, zetal, gammab,
     .                                                     z0, z1, z2)
         if (nzfun .eq. 2) call z_smithe(sgn_kprl, zetal, gammab,
     .                                                     z0, z1, z2)
         if (nzfun .eq. 3) call z_table(sgn_kprl, zetal, gammab,
     .                                               .001, z0, z1, z2)


         al = 1.0 / (xkprl * alpha) * z0
         bl = 1.0 / (xkprl * alpha) * z1
         cl = 1.0 / (xkprl * alpha) * z2




         sig0l = - zieps0 * omgp2 * rho**2 * (exil(labs) - exilp(labs))
     .           * al
         sig1l = - zieps0 * omgp2 * l**2 * exil(labs)/ gamma * al
         sig2l = - eps0 * omgp2 * l * (exil(labs) - exilp(labs)) * al
         sig3l = - zieps0 * omgp2 * 2.0 * exil(labs) * cl
         sig4l = - zieps0 * omgp2 * rho * l * exil(labs)/ gamma  * bl
         sig5l = - eps0 * omgp2 * rho * (exil(labs) - exilp(labs)) * bl



         sig0 = sig0 + sig0l
         sig1 = sig1 + sig1l
         sig2 = sig2 + sig2l
         sig3 = sig3 + sig3l
         sig4 = sig4 + sig4l
         sig5 = sig5 + sig5l


      end do

*     ------------
*     Cold plasma:
*     ------------

c      sig0 = 0.0
c      sig1 = zieps0 * omgrfc * omgp2 / (omgrfc**2 - omgc**2)
c      sig2 = - eps0 * omgc   * omgp2 / (omgrfc**2 - omgc**2)
c      sig3 = zieps0 * omgp2 / omgrfc
c      sig4 = 0.0
c      sig5 = 0.0

      sig1 = sig1 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2
      sig3 = sig3 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2

      kr = sqrt((xkxsav / xkx_cutoff)**2 
     .     + (xkysav / xky_cutoff)**2)
      step = damping * kr**16 / (1. + kr**16)
      sig3 = sig3 * (1.0 + step)

*     -------------------
*     Swanson's rotation:
*     -------------------
      sigxx = sig1 + sig0 * xkbet**2
      sigxy = sig2 - sig0 * xkbet * xkalp
      sigxz = sig4 * xkalp + sig5 * xkbet

      sigyx = - sig2 - sig0 * xkbet * xkalp
      sigyy =   sig1 + sig0 * xkalp**2
      sigyz =   sig4 * xkbet - sig5 * xkalp

      sigzx = sig4 * xkalp - sig5 * xkbet
      sigzy = sig4 * xkbet + sig5 * xkalp
      sigzz = sig3


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
      end
c
c***************************************************************************
c

       subroutine sigmah_stix(i, j, n, m,
     .   bmod, gradprlb,
     .   xm, q, xn, xnuomg,
     .   xkt, omgc, omgp2,
     .   lmin, lmax, nzfun, ibessel,
     .   xkxsav, xkysav, xkzsav, capr,
     .   bx,    by,    bz,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   sigxx, sigxy, sigxz,
     .   sigyx, sigyy, sigyz,
     .   sigzx, sigzy, sigzz,
     .   iflag_gammab, delta0, xk0)

*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     and uses the old left handed U matrix in paper
*     ---------------------------------------------------------

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, iflag_gammab
      real delta0

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, sgn_kprl
      real dakbdkb, xnuomg
      real akprl, gammab, rho, alpha, eps0, omgrf, v0i
      real a, b, bmod, gradprlb
      real bx, by, bz
      real xkxsav, xkysav, xkzsav, capr
      real xkphi
      real xkalp, xkbet, xk0, rgamma

      real uxx, uxy, uxz,
     .     uyx, uyy, uyz,
     .     uzx, uzy, uzz

      complex zi, zfunct, fzeta, omgrfc
      complex zfunct0, zeta0, sig3cold, z0, z1, z2
      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l


      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz

      parameter (lmaxdim = 99)

      complex xil(0: lmaxdim), xilp(0: lmaxdim)
      complex exil(0: lmaxdim), exilp(0: lmaxdim), exilpp(0: lmaxdim)

      complex zetal, zieps0, arg,
     .   al, bl, cl,
     .   gamma, zeta_eff

      complex dpldkb, dplpdkb, dbldkb
      complex dpldakb, dplpdakb, dbldakb


      common/sigcom/zi, eps0, v0i, omgrf



      xme = 9.11e-31
      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rho = alpha / omgc
      xkphi = xkzsav / capr
      omgrfc = omgrf * (1. + zi * xnuomg)

      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi

      xkperp = sqrt(xkalp**2 + xkbet**2)
      akprl = abs(xkprl)



      sgn_kprl = sign(1.0, xkprl)

      if (akprl .lt. 0.1) then

         xkprl = 0.1 * sgn_kprl
         akprl = abs(xkprl)

      end if


      gamma = 0.5 * xkperp**2 * rho**2
      rgamma = real(gamma)


      call besiexp(gamma, lmax, exil, exilp, exilpp, lmaxdim)

      sig0 = 0.0
      sig1 = 0.0
      sig2 = 0.0
      sig3 = 0.0
      sig4 = 0.0
      sig5 = 0.0


      do l = lmin, lmax
         labs = abs(l)



            zetal = (omgrfc - l * omgc) / (xkprl * alpha)


            gammab = l * omgc / (2.0 * alpha * xkprl**2)
     .                                            * gradprlb / bmod



            if(abs(gammab) .gt. 1000.0) gammab = 1000.0
            if(abs(gammab) .lt. .01)gammab = .01



            if (nzfun .eq. 0) call z_approx(sgn_kprl, zetal, 0.0,
     .                                                     z0, z1, z2)
            if (nzfun .eq. 1) call z_approx(sgn_kprl, zetal, gammab,
     .                                                     z0, z1, z2)
            if (nzfun .eq. 2) call z_smithe(sgn_kprl, zetal, gammab,
     .                                                     z0, z1, z2)
            if (nzfun .eq. 3) call z_table(sgn_kprl, zetal, gammab,
     .                                               .001, z0, z1, z2)


            al = 1.0 / (xkprl * alpha) * z0
            bl = 1.0 / (xkprl * alpha) * z1
            cl = 1.0 / (xkprl * alpha) * z2




         sig0l = - zieps0 * omgp2 * rho**2 * (exil(labs) - exilp(labs))
     .           * al
         sig1l = - zieps0 * omgp2 * l**2 * exil(labs)/ gamma * al
         sig2l = - eps0 * omgp2 * l * (exil(labs) - exilp(labs)) * al
         sig3l = - zieps0 * omgp2 * 2.0 * exil(labs) * cl
         sig4l = - zieps0 * omgp2 * rho * l * exil(labs)/ gamma  * bl
         sig5l = - eps0 * omgp2 * rho * (exil(labs) - exilp(labs)) * bl



         sig0 = sig0 + sig0l
         sig1 = sig1 + sig1l
         sig2 = sig2 + sig2l
         sig3 = sig3 + sig3l
         sig4 = sig4 + sig4l
         sig5 = sig5 + sig5l


      end do

      sig1 = sig1 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2
      sig3 = sig3 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2

*     -------------------
*     Swanson's rotation:
*     -------------------
      sigxx = sig1 + sig0 * xkbet**2
      sigxy = sig2 - sig0 * xkbet * xkalp
      sigxz = sig4 * xkalp + sig5 * xkbet

      sigyx = - sig2 - sig0 * xkbet * xkalp
      sigyy =   sig1 + sig0 * xkalp**2
      sigyz =   sig4 * xkbet - sig5 * xkalp

      sigzx = sig4 * xkalp - sig5 * xkbet
      sigzy = sig4 * xkbet + sig5 * xkalp
      sigzz = sig3


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
      end
c
c***************************************************************************
c



      subroutine sigmac_stix(i, j, n, m,
     .   xm, q, xn, xnuomg,
     .   xkt, omgc, omgp2,
     .   lmin, lmax, nzfun, ibessel,
     .   xkxsav, xkysav, xkzsav, capr,
     .   bx,    by,    bz,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   sigxx, sigxy, sigxz,
     .   sigyx, sigyy, sigyz,
     .   sigzx, sigzy, sigzz)

*     -------------------------------------------------------
*     This routine uses the old left handed U matrix in paper
*     -------------------------------------------------------

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xkb, akb
      real dakbdkb
      real akprl, gammab, alpha, eps0, omgrf, v0i
      real a, b
      real bx, by, bz
      real xkxsav, xkysav, xkzsav, capr
      real xkphi
      real xkalp, xkbet, sqx, xnuomg

      real uxx, uxy, uxz,
     .     uyx, uyy, uyz,
     .     uzx, uzy, uzz

      complex zi, omgrfc
      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l


      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz


      parameter (lmaxdim = 99)

      complex xil(0: lmaxdim), xilp(0: lmaxdim)
      complex exil(0: lmaxdim), exilp(0: lmaxdim), exilpp(0: lmaxdim)

      complex tpl, tml, tplp, tplpp, zetalp, zetalm, zieps0,
     .   bl, gamma

      complex dpldkb, dplpdkb, dbldkb
      complex dpldakb, dplpdakb, dbldakb


      common/sigcom/zi, eps0, v0i, omgrf


      zieps0 = zi * eps0
      xkphi = xkzsav / capr


      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi


      omgrfc = omgrf * (1. + zi * xnuomg)


      sig0 = 0.0
      sig1 = zieps0 * omgrfc * omgp2 / (omgrfc**2 - omgc**2)
      sig2 = - eps0 * omgc   * omgp2 / (omgrfc**2 - omgc**2)
      sig3 = zieps0 * omgp2 / omgrfc
      sig4 = 0.0
      sig5 = 0.0


*     -------------------
*     Swanson's rotation:
*     -------------------
      sigxx = sig1 + sig0 * xkbet**2
      sigxy = sig2 - sig0 * xkbet * xkalp
      sigxz = sig4 * xkalp + sig5 * xkbet

      sigyx = - sig2 - sig0 * xkbet * xkalp
      sigyy =   sig1 + sig0 * xkalp**2
      sigyz =   sig4 * xkbet - sig5 * xkalp

      sigzx = sig4 * xkalp - sig5 * xkbet
      sigzy = sig4 * xkbet + sig5 * xkalp
      sigzz = sig3


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
      end
c
c*******************************************************************************
c

      subroutine besiexp(gamma, lmax, expbes, expbesp, expbespp,
     1    lmaxdim)
c*****   calculates exp(-gamma) times the modified bessel functions
c*****   (and derivatives) of order up to lmax

      implicit none

      integer lmax, nmax, ier, l, lmaxdim
      real gammod
      complex gamma, expbes(0:lmaxdim), expbesp(0:lmaxdim),
     1   expbespp(0:lmaxdim),
     1   xil(0:lmaxdim), xilp(0:lmaxdim), xilpp(0:lmaxdim), exgam

      complex b(100)

      exgam = exp(-gamma)
      gammod = cabs(gamma)

      if(gammod .le. 700.)then
         nmax = lmax + 1
         call besic(gamma, nmax, b, ier)
         if(ier .ne. 0)write(6,100) ier

         do l = 0, lmax
            xil(l) = b(l+1)
         end do

         do l = 0, lmax
           if(l .eq. 0) xilp(0) = xil(1)
           if(l .ne. 0) xilp(l) = xil(l-1) - l / gamma * xil(l)
           xilpp(l) = xil(l) * (1. + (l / gamma)**2) - xilp(l) / gamma
           expbes(l) = exgam * xil(l)
           expbesp(l) = exgam * xilp(l)
           expbespp(l) = exgam * xilpp(l)
         end do
      end if

      if(gammod .gt. 700.)then
         do l = 0, lmax
            call bes_asym(gamma, l, expbes(l), expbesp(l))
            expbespp(l) = expbes(l) * (1. + (l / gamma)**2)
     1          - expbesp(l) / gamma
         end do
      end if

  100 format('ier = ', i5, 'besic failed')
      return
      end


c
c***************************************************************************
c
      subroutine bes_asym(z, n, exil, exilp)

      implicit none

      integer mu, n
      real pi
      complex z, exil, exilp
      data pi/3.141592654/

      mu = 4 * n**2
      exil =  1.0 / csqrt(2.0 * pi * z)
     1   * (1.0
     1   - (mu - 1)/(8.0 * z)
     1   + (mu - 1) * (mu - 9) / (2.0 * (8.0 * z)**2)
     1   - (mu - 1) * (mu - 9) * (mu - 25) / (6.0 * (8.0 * z)**3)  )
      exilp = 1.0 / csqrt(2.0 * pi * z)
     1   * (1.0
     1   - (mu + 3)/(8.0 * z)
     1   + (mu - 1) * (mu + 15) / (2.0 * (8.0 * z)**2)
     1   - (mu - 1) * (mu - 9) * (mu + 35) / (6.0 * (8.0 * z)**3)  )

      return
      end

c
c***************************************************************************
c

      SUBROUTINE BESIC (Z,N,B,IER)
C
C PACKAGE BESC           (NOTE---DOCUMENTATION FOR INDIVIDUAL ROUTINES
C                        FOLLOWS THE GENERAL PACKAGE INFORMATION.)
C
C LATEST REVISION        NOVEMBER 1978
C
C PURPOSE                TO CALCULATE BESSEL FUNCTIONS I AND J FOR
C                        COMPLEX ARGUMENT AND INTEGER ORDERS.
C
C USAGE                  TO CALCULATE I BESSEL FUNCTIONS---
C
C                          CALL BESIC(Z,N,B,IER)
C
C                        TO CALCULATE J BESSEL FUNCTIONS---
C
C                          CALL BESJC(Z,N,B,IER)
C
C ENTRY POINTS           BESIC,BESJC,B2SLCI
C
C COMMON BLOCKS          NONE
C
C I/O                    NONE, EXCEPT FOR ERROR MESSAGES PRODUCED BY
C                        CALLING THE ERROR HANDLING ROUTINE ULIBER.
C
C REQUIRED LIBRARY       NONE
C ROUTINES
C
C SPECIALIST             RUSSELL K. REW
C
C LANGUAGE               FORTRAN
C
C HISTORY                DAVID J. SOOKNE WROTE THE CORE ROUTINE B2SLCI,
C                        WHICH WAS ORIGINALLY DOUBLE PRECISION.  THE
C                        SINGLE PRECISION VERSION, AND THE USER ENTRIES
C                        BESIC AND BESJC WERE WRITTEN BY RUSSELL K. REW
C                        IN SEPTEMBER, 1977.
C
C SPACE REQUIRED         1712 (OCTAL) = 970 (DECIMAL) ON THE 7600.
C
C PORTABILITY            THIS PACKAGE CONFORMS TO THE 1966 ANSI STANDARD
C                        AS CONFIRMED BY THE PFORT VERIFIER.  THERE ARE
C                        FOUR MACHINE DEPENDENT CONSTANTS IN B2SLCI THAT
C                        ARE DESCRIBED IN COMMENT CARDS.
C
C REQUIRED RESIDENT      ULIBER (AN ERROR ROUTINE),SQRT,EXP,COS,SIN
C ROUTINES
C
C REFERENCES             "BESSEL FUNCTIONS I AND J OF COMPLEX
C                        ARGUMENT AND INTEGER ORDER" AND
C                        "CERTIFICATION OF AN ALGORITHM FOR
C                        BESSEL FUNCTIONS OF COMPLEX ARGUMENT"
C                        BOTH BY DAVID J. SOOKNE, JOURNAL OF REASEARCH
C                        OF THE NATIONAL BUREAU OF STANDARDS-B. MATHE-
C                        MATICAL SCIENCES, VOL. 77A, NOS. 3 AND 4, JULY-
C                        DECEMBER, 1973.
C
C METHOD                 BACKWARD RECURSION WITH STRICT CONTROL OF
C                        ERROR.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE BESIC (Z,N,B,IER)
C
C DIMENSION OF           B(N+1)
C ARGUMENTS
C
C LATEST REVISION        OCTOBER 1978
C
C PURPOSE                TO CALCULATE I BESSEL FUNCTIONS FOR COMPLEX
C                        ARGUMENT AND INTEGER ORDER
C
C USAGE                  CALL BESIC (Z,N,B,IER)
C
C ARGUMENTS
C
C ON INPUT               Z
C                          COMPLEX ARGUMENT FOR WHICH I BESSEL FUNCTIONS
C                          ARE TO BE CALCULATED.  ABS(AIMAG(Z)) MUST BE
C                          LESS THAN THE LARGEST REAL ARGUMENT THAT THE
C                          FORTRAN FUNCTION EXP CAN HANDLE.
C
C                        N
C                          INTEGER, THE HIGHEST ORDER TO BE CALCULATED.
C                          N MUST BE GREATER THAN OR EQUAL TO ZERO.
C
C ON OUTPUT              B
C                          COMPLEX VECTOR OF LENGTH N+1 CONTAINING THE
C                          BESSEL FUNCTION VALUES I-SUB-0(Z),I-SUB-1(Z),
C                          ...,I-SUB-N(Z) IN B(1),B(2),...,B(N+1).
C                        IER
C                          AN INTEGER ERROR FLAG.
C                          =0 IF ALL DESIRED ORDERS HAVE BEEN CALCULATED
C                             SATISFACTORILY,
C                          =1 IF ABS(AIMAG(Z)) IS TOO LARGE,
C                          =2 IF N IS LESS THAN ZERO,
C                          =2+K IF ONLY THE FIRST K RESULTS ARE CORRECT.
C                             IN THE RETURNED VALUES B(M) FOR M GREATER
C                             THAN K, APPROXIMATELY THE LAST
C                             ALOG10(ABS(B(M)/B(K))) SIGNIFICANT DIGITS
C                             ARE IN ERROR.
C
C ENTRY POINTS           BESIC
C
C SPECIAL CONDITIONS     NONE
C
C COMMON BLOCKS          NONE
C
C I/O                    NONE, EXCEPT FOR ERROR MESSAGES PRODUCED BY
C                        CALLING THE ERROR HANDLING ROUTINE ULIBER.
C
C PRECISION              SINGLE
C
C SPECIALIST             RUSS REW, NCAR, BOULDER, COLORADO
C
C ACCURACY               IN TESTS RUN ON THE 7600 WITH ORDERS FROM 0
C                        THROUGH 10 USING RANDOM VALUES OF THE ARGUMENT
C                        WITH ABSOLUTE VALUE LESS THAN 50 BUT
C                        CONCENTRATED AROUND THE ORIGIN, THE MAXIMUM
C                        RELATIVE ERROR (OR ABSOLUTE ERROR WHEN IT WAS
C                        LARGER) OBSERVED WAS ABOUT 8.1E-14.
C
C TIMING                 ON NCAR"S CONTROL DATA 7600, BESIC TAKES ABOUT
C                        .32+.008*N MILLISECONDS WHEN Z=(1.0,1.0).
C
C PORTABILITY            ANSI 1966 STANDARD
C
C
C
C
C
C
      COMPLEX         Z          ,B(100)
      real b2(2, 100)


      DATA IORJ/1/,XLARGE/741.66/

      NB = N+1



      CALL B2SLCI (REAL(Z),AIMAG(Z),NB,IORJ,B2,NCALC)

      do n = 1, nb
         b(n) = cmplx(b2(1,n), b2(2,n))
      end do

      IER = 0
      IF (NCALC .EQ. NB) GO TO 103
      IF (NCALC .GE. 0) GO TO 102
      IF (N .GE. 0) GO TO 101
      IER = 2
      CALL ULIBER (IER,25H IN BESIC, N OUT OF RANGE,25)
      GO TO 103
  101 IER = 1
      CALL ULIBER (IER,25H IN BESIC, X OUT OF RANGE,25)
      GO TO 103
  102 IER = 2+NCALC
      CALL ULIBER (IER,40H IN BESIC, ACCURACY LOST FOR SOME ORDERS,40)
  103 RETURN
      END

c
c***************************************************************************
c

      function fzeta (arg)
      complex a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, arg, aux0
     1   , aux1, term, z, zz, fzeta
      data d1r/0.0/
      data d1i/1.77245385090551/
      data d2r/0.0/
      data d2i/3.54490770181103/
      data d3r/0.0/
      data d3i/7.08981540362206/
      data d4/0.33333333333333/
      data eps/1.0E-07/

c     data d4/0.33333333333333/
c     common/zetcom/eps
cray  code analysis
cray  optimize
c
      i = 0
      z = arg
      zz = z*z
      x = real(z)
      y = aimag(z)
      d1 = cmplx(d1r,d1i)
      d2 = cmplx(d2r,d2i)
      d3 = cmplx(d3r,d3i)
      ymag = abs(y)
      if (ymag - 1.0 .ge. 0.) then
c
c     continued fraction method: abs(y).ge.1.0
c
         y0 = y
         y = ymag
         aux1 = 1.5 - z*z
         aux2 = 0.0
         del = 1.5
         a1 = 0.0
         a2 = -1.0
         b1 = 1.0
         b2 = aux1
         c1 = a2/b2
c
  100    continue
         aux1 = aux1 + 2.0
         aux2 = aux2 - del
         del = del + 2.0
         a3 = aux1*a2 + aux2*a1
         b3 = aux1*b2 + aux2*b1
         c2 = a3/b3
         c3 = c2 - c1
         c3r = real(c3)
         c3i = aimag(c3)
         if (abs(c3r) + abs(c3i) .lt. eps) go to 110
         a1 = a2
         a2 = a3
         b1 = b2
         b2 = b3
         c1 = c2
         go to 100
  110    continue
         if (y0 .lt. 0.) then
            y = y0
            c2 = conjg(c2) - d3*z*exp(-zz)
         endif
         aux0 = -(0.5*c2 + 1.0)/z
      else
c
c     asymptotic series method: abs(x).ge.4.0 and abs(y).lt.1.0
c
         xmag = abs(x)
         if (xmag - 4.0 .lt. 0.) go to 130
         term = 1.0/z
         aux0 = -term
         aux1 = 0.5*term**2
         p = 1.0
         if (y .le. 0.) then
            if (y .ne. 0.) then
               aux0 = aux0 + d2*exp(-zz)
            else
               aux0 = aux0 + d1*exp(-zz)
            endif
         endif
  120    continue
         term = aux1*term*p
         aux0 = aux0 - term
         p = p + 2.0
         termr = real(term)
         termi = aimag(term)
c     if(abs(termr)+abs(termi).lt.eps)30,18
         if (abs(termr) + abs(termi) .lt. eps) go to 160
         go to 120
c
c     power series method: abs(x).lt.4.0 and abs(y).lt.1.0
c
  130    continue
         aux0 = 1.0
         aux1 = -(zz + zz)
         aux2 = eps/(eps + xmag + ymag)
         term = d4*aux1
         p = 3.0
  140    continue
         aux0 = aux0 + term
         termr = real(term)
         termi = aimag(term)
c     if(abs(termr)+abs(termi).lt.aux2)26,24
         if (abs(termr) + abs(termi) .lt. aux2) go to 150
         p = p + 2.0
         term = aux1*term/p
         go to 140
  150    continue
         aux0 = d1*exp(-zz) - 2.0*z*aux0
      endif
  160 continue
      fzeta = aux0
      if (i .le. 0) return
      fzeta = -2.0*(1.0 + arg*aux0)
      return
      end
c
c***************************************************************************
c


      subroutine b2slci (x,y,nb,ize,b,ncalc)
c
c this routine calculates bessel functions i and j of
c complex argument and integer order.
c
c
c      explanation of variables in the calling sequence
c
c x     real part of the complex argument
c       for which i*s or j*s are to be calculated.  if i*s
c       are to be calculated, abs(x) must not exceed exparg
c       (which see below).
c y     imaginary part of the argument.  if j*s are to be
c       calculated, abs(y) must not exceed exparg.
c nb    integer type.  1 + highest order to be calculated.
c       it must be positive.
c ize   integer type.  zero if j*s are to be calculated, 1
c       if i*s are to be calculated.
c b     real array dimensioned b(2,nb), need not be initialized by user.
c       if the routine terminates normally, (ncalc=nb), it returns
c       j(or i)-sub-zero through j(or i)-sub-nb-minus-one of z in b.
c       the real parts of the results are returned in b(1,n) and the
c       corresponding imaginary parts are returned in b(2,n) for n=1,2,
c       ...,nb.  (in the documentation below, br(n) and bi(n) refer to
c       b(1,n) and b(2,n), respectively.)
c ncalc integer type, need not be initialized by user.
c       before using the results, the user should check that
c       ncalc=nb, i.e. all orders have been calculated to
c       the desired accuracy.  see error returns below.
c
c
c       explanation of machine-dependent constants
c
c nsig  decimal significance desired.  should be set to
c       ifix(alog10(2)*nbit+1), where nbit is the number of
c       bits in the mantissa of a real variable.
c       setting nsig higher will increase cpu time without
c       increasing accuracy, while setting nsig lower will
c       decrease accuracy.  if only single-precision
c       accuracy is desired, replace nbit by the number of
c       bits in the mantissa of a single-precision variable.
c       the relative truncation error is limited to t=.5*10
c       **-nsig for order greater than abs(z), and for order
c       less than abs(z) (general test), the relative error
c       is limited to t for function values of magnitude at
c       least 1, and the absolute error is limited to t for
c       smaller values.
c nten  largest integer k such that 10**k is machine-
c       representable in real.
c largez upper limit on the magnitude of z.  bear in mind
c       that if abs(z)=n, then at least n iterations of the
c       backward recursion will be executed.
c exparg largest real argument that the library
c       exp routine can handle.
c
c
c                            error returns
c
c       let g denote either i or j.
c       in case of an error, ncalc.ne.nb, and not all g*s
c  are calculated to the desired accuracy.
c       if ncalc.lt.0, an argument is out of range.  nb.le.0
c  or ize is neither 0 nor 1 or ize=0 and abs(y).gt.exparg,
c  or ize=1 and abs(x).gt.exparg.  in this case, the vectors
c  br and bi are not calculated, and ncalc is set to
c  min0(nb,0)-1 so ncalc.ne.nb.
c       nb.gt.ncalc.gt.0 will occur if nb.gt.magz and abs(g-
c  sub-nb-of-z/g-sub-magx+np-of-z).lt.10.**(nten/2), i.e. nb
c  is much greater than magz.  in this case, br(n) and bi(n)
c  are calculated to the desired accuracy for n.le.ncalc,
c  but for ncalc.lt.n.le.nb, precision is lost.  if n.gt.
c  ncalc and abs(g(ncalc-1)/g(n-1)).eq.10**-k, then the last
c  k significant figures of g(n-1) (=br(n)+i*bi(n)) are
c  erroneous.  if the user wishes to calculate g(n-1) to
c  higher accuracy, he should use an asymptotic formula for
c  large order.
c
      real            x          ,y          ,b          ,pr         ,
     1                pi         ,plastr     ,plasti     ,poldr      ,
     2                poldi      ,psaver     ,psavei     ,exparg     ,
     3                test       ,tover      ,tempar     ,tempai     ,
     4                tempbr     ,tempbi     ,tempcr     ,tempci     ,
     5                sign       ,sumr       ,sumi       ,zinvr      ,
     6                zinvi
      dimension       b(2,nb)
c
c machine dependent data for 7600 and cray
c
c      data nsig,nten,largez,exparg
c     1    /15,293,100000,741.66e0/
c     data nsig/15/,nten/293/,largez/100000/,exparg/741.66/

C MACHINE DEPENDENT DATA FOR VAX
C
c     DATA NSIG,NTEN,LARGEZ,EXPARG
c    1    /7,38,10000.,88.02/

c
c machine dependent data for IBM 580 workstation
c
c     data nsig/15/,nten/308/,largez/100000/,exparg/709.78/
      data nsig/15/,nten/308/,largez/100000/,exparg/709.78/

      tempar = sqrt(x*x+y*y)
      magz = ifix((tempar))
      if (nb.gt.0 .and. magz.le.largez .and.
     1    ((ize.eq.0 .and. abs(y).le.exparg) .or.
     2                               (ize.eq.1 .and. abs(x).le.exparg)))
     3    go to 101
c
c error return -- z, nb, or ize is out of range
c
      ncalc = min0(nb,0)-1
      return
  101 sign = (float(1-2*ize))
      ncalc = nb
c
c use 2-term ascending series for small z
c
      if (tempar**4 .lt. .1e0**nsig) go to 136
c
c initialize the calculation of the p*s
c
      nbmz = nb-magz
      n = magz+1
      if (abs(x) .lt. abs(y)) go to 102
      zinvr = 1.e0/(x+y*y/x)
      zinvi = -y*zinvr/x
      go to 103
  102 zinvi = -1.e0/(y+x*x/y)
      zinvr = -x*zinvi/y
  103 plastr = 1.e0
      plasti = 0.e0
      pr = sign*(float(2*n))*zinvr
      pi = sign*(float(2*n))*zinvi
      test = 2.e0*1.e1**nsig
      m = 0
      if (nbmz .lt. 3) go to 105
c
c calculate p*s until n=nb-1.  check for possible overflow.
c
c the following devious computation replaces
c     tover=10.0**(nten-nsig)
c and is necessitated by the poor power algorithm on the ncar 7600
c in order to prevent overflow
c
      intemp = nten-nsig
      inhlf = intemp/2
      tover = 10.0**inhlf*10.0**(intemp-inhlf)
      nstart = magz+2
      nend = nb-1
      do 104 n=nstart,nend
         poldr = plastr
         poldi = plasti
         plastr = pr
         plasti = pi
         pr = sign*((float(2*n))*(plastr*zinvr-plasti*zinvi)-poldr)
         pi = sign*((float(2*n))*(plasti*zinvr+plastr*zinvi)-poldi)
         if ((pr/tover)**2+(pi/tover)**2-1.e0) 104,104,106
  104 continue
      n = nend
c
c calculate special significance test for nbmz.gt.2.
c
      tempbi = max(abs(pr),abs(pi))
      tempbi = tempbi*sqrt(2.e0*1.e1**nsig*
     1         sqrt(((pr/tempbi)**2+(pi/tempbi)**2)*
     2                         ((plastr/tempbi)**2+(plasti/tempbi)**2)))
      test = max(test,tempbi)
c
c calculate p*s until significance test is passed.
c
  105 n = n+1
      poldr = plastr
      poldi = plasti
      plastr = pr
      plasti = pi
      pr = sign*((float(2*n))*(plastr*zinvr-plasti*zinvi)-poldr)
      pi = sign*((float(2*n))*(plasti*zinvr+plastr*zinvi)-poldi)
      if ((pr/test)**2+(pi/test)**2 .lt. 1.e0) go to 105
      if (m .eq. 1) go to 110
c
c calculate strict variant of significance test, and
c calculate p*s until this test is passed.
c
      m = 1
      tempbi = max(abs(pr),abs(pi))
      tempbr = sqrt(((pr/tempbi)**2+(pi/tempbi)**2)/
     1                          ((plastr/tempbi)**2+(plasti/tempbi)**2))
      tempbi = (float(n+1))/tempar
      if (tempbr+1.e0/tempbr .gt. 2.e0*tempbi)
     1    tempbr = tempbi+sqrt(tempbi**2-1.e0)
      test = test/sqrt(tempbr-1.e0/tempbr)
      if ((pr/test)**2+(pi/test)**2-1.e0) 105,110,110
  106 nstart = n+1
c
c to avoid overflow, normalize p*s by dividing by tover.
c calculate p*s until unnormalized p would overflow.
c
      pr = pr/tover
      pi = pi/tover
      plastr = plastr/tover
      plasti = plasti/tover
      psaver = pr
      psavei = pi
      tempcr = plastr
      tempci = plasti
      test = 1.e1**(2*nsig)
  107 n = n+1
      poldr = plastr
      poldi = plasti
      plastr = pr
      plasti = pi
      pr = sign*((float(2*n))*(plastr*zinvr-plasti*zinvi)-poldr)
      pi = sign*((float(2*n))*(plasti*zinvr+plastr*zinvi)-poldi)
      if (pr**2+pi**2 .le. test) go to 107
c
c calculate backward test, and find ncalc, the highest n
c such that the test is passed.
c
      tempbr = sqrt((plastr**2+plasti**2)/(poldr**2+poldi**2))
      tempbi = (float(n))/tempar
      if (tempbr+1.e0/tempbr .gt. 2.e0*tempbi)
     1    tempbr = tempbi+sqrt(tempbi**2-1.e0)
      test = .5e0*(1.e0-1.e0/tempbr**2)/1.e1**nsig
      test = ((plastr**2+plasti**2)*test)*((poldr**2+poldi**2)*test)
      pr = plastr*tover
      pi = plasti*tover
      n = n-1
      nend = min0(nb,n)
      do 108 ncalc=nstart,nend
         poldr = tempcr
         poldi = tempci
         tempcr = psaver
         tempci = psavei
         psaver = sign*((float(2*n))*(tempcr*zinvr-tempci*zinvi)-poldr)
         psavei = sign*((float(2*n))*(tempci*zinvr+tempcr*zinvi)-poldi)
         if ((psaver**2+psavei**2)*(tempcr**2+tempci**2)-test)
     1       108,108,109
  108 continue
      ncalc = nend+1
  109 ncalc = ncalc-1
c
c the coefficient of b(n) in the normalization sum is
c m*sqrt(-1)**imag, where m=-2,0, or 2, and imag is 0 or 1.
c calculate recursion rules for m and imag, and initialize
c them.
c
  110 n = n+1
      tempbr = (float(ize))*x+(float(1-ize))*y
      ipos = 0
      if (tempbr) 111,112,111
  111 ipos = ifix((1.1e0*tempbr/abs(tempbr)))
  112 mrecur = 4*((2+ize+ipos)/2)-3-2*(ize+ipos)
      k = 2+ipos+2*ize*ipos**2-ize
      l = n-4*(n/4)
      mlast = 2+8*((k*l)/4)-4*((k*l)/2)
      if (ipos.eq.0 .and. (l.eq.1 .or. l.eq.3)) mlast = 0
      l = l+3-4*((l+3)/4)
      m = 2+8*((k*l)/4)-4*((k*l)/2)
      if (ipos.eq.0 .and. (l.eq.1 .or. l.eq.3)) m = 0
      imrecr = (1-ize)*ipos**2
      imag = imrecr*(l-2*(l/2))
c
c initialize the backward recursion and the normalization
c sum.
c
      tempbr = 0.e0
      tempbi = 0.e0
      if (abs(pi) .gt. abs(pr)) go to 113
      tempar = 1.e0/(pr+pi*(pi/pr))
      tempai = -(pi*tempar)/pr
      go to 114
  113 tempai = -1.e0/(pi+pr*(pr/pi))
      tempar = -(pr*tempai)/pi
  114 if (imag .ne. 0) go to 115
      sumr = (float(m))*tempar
      sumi = (float(m))*tempai
      go to 116
  115 sumr = -(float(m))*tempai
      sumi = (float(m))*tempar
  116 nend = n-nb
      if (nend) 123,120,117
c
c recur backward via difference equation calculating (but
c not storing) br(n) and bi(n) until n=nb.
c
  117 do 119 l=1,nend
         n = n-1
         tempcr = tempbr
         tempci = tempbi
         tempbr = tempar
         tempbi = tempai
         pr = (float(2*n))*zinvr
         pi = (float(2*n))*zinvi
         tempar = pr*tempbr-pi*tempbi-sign*tempcr
         tempai = pr*tempbi+pi*tempbr-sign*tempci
         imag = (1-imag)*imrecr
         k = mlast
         mlast = m
         m = k*mrecur
         if (imag .ne. 0) go to 118
         sumr = sumr+(float(m))*tempar
         sumi = sumi+(float(m))*tempai
         go to 119
  118    sumr = sumr-(float(m))*tempai
         sumi = sumi+(float(m))*tempar
  119 continue
c
c store br(nb), bi(nb)
c
  120 b(1,n) = tempar
      b(2,n) = tempai
      if (n .gt. 1) go to 121
c
c nb=1.  since 2*tempar and 2*tempai were added to sumr and
c sumi respectively, we must subtract tempar and tempai
c
      sumr = sumr-tempar
      sumi = sumi-tempai
      go to 130
c
c calculate and store br(nb-1),bi(nb-1)
c
  121 n = n-1
      pr = (float(2*n))*zinvr
      pi = (float(2*n))*zinvi
      b(1,n) = pr*tempar-pi*tempai-sign*tempbr
      b(2,n) = pr*tempai+pi*tempar-sign*tempbi
      if (n .eq. 1) go to 129
      imag = (1-imag)*imrecr
      k = mlast
      mlast = m
      m = k*mrecur
      if (imag .ne. 0) go to 122
      sumr = sumr+(float(m))*b(1,n)
      sumi = sumi+(float(m))*b(2,n)
      go to 125
  122 sumr = sumr-(float(m))*b(2,n)
      sumi = sumi+(float(m))*b(1,n)
      go to 125
c
c n.lt.nb, so store br(n), bi(n), and set higher orders zero
c
  123 b(1,n) = tempar
      b(2,n) = tempai
      nend = -nend
      do 124 l=1,nend
         b(1,n+l) = 0.e0
         b(2,n+l) = 0.e0
  124 continue
  125 nend = n-2
      if (nend .eq. 0) go to 128
c
c calculate via difference equation and store br(n),bi(n),
c until n=2
c
      do 127 l=1,nend
         n = n-1
         pr = (float(2*n))*zinvr
         pi = (float(2*n))*zinvi
         b(1,n) = pr*b(1,n+1)-pi*b(2,n+1)-sign*b(1,n+2)
         b(2,n) = pr*b(2,n+1)+pi*b(1,n+1)-sign*b(2,n+2)
         imag = (1-imag)*imrecr
         k = mlast
         mlast = m
         m = k*mrecur
         if (imag .ne. 0) go to 126
         sumr = sumr+(float(m))*b(1,n)
         sumi = sumi+(float(m))*b(2,n)
         go to 127
  126    sumr = sumr-(float(m))*b(2,n)
         sumi = sumi+(float(m))*b(1,n)
  127 continue
c
c calculate and store br(1), bi(1)
c
  128 b(1,1) = 2.e0*(b(1,2)*zinvr-b(2,2)*zinvi)-sign*b(1,3)
      b(2,1) = 2.e0*(b(1,2)*zinvi+b(2,2)*zinvr)-sign*b(2,3)
  129 sumr = sumr+b(1,1)
      sumi = sumi+b(2,1)
c
c calculate normalization factor, tempar +i*tempai
c
  130 if (ize .eq. 1) go to 131
      tempcr = (float(ipos))*y
      tempci = (float(-ipos))*x
      go to 132
  131 tempcr = (float(ipos))*x
      tempci = (float(ipos))*y
  132 tempcr = exp(tempcr)
      tempbr = cos(tempci)
      tempbi = sin(tempci)
      if (abs(sumr) .lt. abs(sumi)) go to 133
      tempci = sumi/sumr
      tempcr = (tempcr/sumr)/(1.e0+tempci*tempci)
      tempar = tempcr*(tempbr+tempbi*tempci)
      tempai = tempcr*(tempbi-tempbr*tempci)
      go to 134
  133 tempci = sumr/sumi
      tempcr = (tempcr/sumi)/(1.e0+tempci*tempci)
      tempar = tempcr*(tempbr*tempci+tempbi)
      tempai = tempcr*(tempbi*tempci-tempbr)
c
c normalize
c
  134 do 135 n=1,nb
         tempbr = b(1,n)*tempar-b(2,n)*tempai
         b(2,n) = b(1,n)*tempai+b(2,n)*tempar
         b(1,n) = tempbr
  135 continue
      return
c
c two-term ascending series for small z
c
  136 tempar = 1.e0
      tempai = 0.e0
      tempcr = .25e0*(x*x-y*y)
      tempci = .5e0*x*y
      b(1,1) = 1.e0-sign*tempcr
      b(2,1) = -sign*tempci
      if (nb .eq. 1) go to 138
      do 137 n=2,nb
         tempbr = (tempar*x-tempai*y)/(float(2*n-2))
         tempai = (tempar*y+tempai*x)/(float(2*n-2))
         tempar = tempbr
         tempbr = (float(n))
         b(1,n) = tempar*(1.e0-sign*tempcr/tempbr)+tempai*tempci/tempbr
         b(2,n) = tempai*(1.e0-sign*tempcr/tempbr)-tempar*tempci/tempbr
  137 continue
  138 return
c
c revision history---
c
c october 1978     first added to nssl
c november 1978    changed the value of the machine dependent constant
c                  nten from 322 to 293 to correct an undetected
c                  underflow and subsequent division of 0 by 0 when
c                  besir was called with large n (greater than 140).
c-----------------------------------------------------------------------
      end
c
c uliber     from portlib                                  12/27/80
      subroutine uliber (ierr,messg,nmessg)
c
c dimension of           messg(1+(nmessg-1)/ncpwd), where ncpwd is the
c arguments              number of characters that can be stored in one
c                        integer word.
c
c latest revision        april 1977.
c
c purpose                prints an error message.
c
c usage                  call uliber (ierr,messg,nmessg)
c
c arguments
c
c on input               ierr
c                          error number.  if ierr .lt. 32, the error is
c                          considered to be non-fatal and uliber returns
c                          after printing the error message.  if ierr
c                          .ge. 32 the error is considered fatal, and
c                          uliber stops after printing the message.
c                        messg
c                          the error message to be printed.
c                        nmessg
c                          the length of the message, in characters.
c
      dimension messg(1)
      integer erunit,prunit
c
c erunit should be set to the local unit number for error messages.
c
      data erunit/6/,prunit/6/
c
c      if (erunit .eq. 0) write(prunit,1000)
c
c 1000 format(50h1uliber, a subroutine to print error messages, has/
c     1       50h been called, but no local implementation of      /
c     2       50h uliber has been provided.  to properly implement /
c     3       50h this subroutine for the local machine and        /
c     4       50h environment, please see the comments in uliber.  /)
c
c replace the if statement and format statement above with the following
c code where $ncpwd should be replaced with the number of characters
c that can be stored in one integer word, and $n should be replaced by
c the value of 1+(79/$ncpwd).
c
       mm=1+(nmessg-1)/8
       write(erunit,1000) ierr,(messg(i),i=1,mm)
1000   format(' ****error #' ,i5,25h,error message follows.../
     1       (18 a8 ))
c
   10 if (ierr .ge. 32) stop
      return
      end
c
c**********************************************************************
c
      function bessi(n,x)

      implicit none

      integer n, iacc
      real bessi, x, bigno, bigni
      parameter (iacc = 40, bigno = 1.0e10, bigni = 1.0e-10)
c     Uses bessi0
c     Returns the modified Bessel function In(x) for any real x and n.ge.2
      integer j, m
      real bi, bim, bip, tox, bessi0
      if (n .lt. 2) pause 'bad argument n in bessi'
      if (x .eq. 0.) then
         bessi = 0.
      else
         tox = 2.0/abs(x)
         bip = 0.0
         bi = 1.0
         bessi = 0.
         m = 2 * ((n + int(sqrt(float(iacc * n)))))
         do j = m, 1, -1
            bim = bip + float(j) * tox * bi
            bip = bi
            bi = bim
            if (abs(bi) .gt. bigno) then
               bessi = bessi * bigni
               bi = bi * bigni
               bip = bip * bigni
            endif
            if (j .eq. n) bessi = bip
         enddo
         bessi = bessi * bessi0(x) / bi
         if(x .lt. 0. .and. mod(n,2) .eq. 1) bessi = -bessi
      endif
      return
      end
c
c**********************************************************************
c
      function bessi0(x)

      implicit none

      real bessi0, x
c      Returns the modified Bessel function I0(x) for any real x
      real ax
      double precision p1, p2, p3, p4, p5, p6, p7,
     1                 q1, q2, q3, q4, q5, q6, q7, q8, q9, y
      save  p1, p2, p3, p4, p5, p6, p7,
     1      q1, q2, q3, q4, q5, q6, q7, q8, q9

      data p1, p2, p3, p4, p5, p6, p7 / 1.0d0, 3.5156229d0, 3.0899424d0,
     1   1.2067492d0, 0.2659732d0, 0.360768d-1, 0.45813d-2 /

      data q1, q2, q3, q4, q5, q6, q7, q8, q9/0.39894228d0,0.1328592d-1,
     1   0.225319d-2,-0.157565d-2, 0.916281d-2, -0.2057706d-1,
     1   0.2635537d-1, -0.1647633d-1, 0.392377d-2 /

      if (abs(x) .lt. 3.75) then
         y = (x / 3.75)**2
         bessi0 = p1 + y * (p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      else
         ax = abs(x)
         y = 3.75 / ax
         bessi0 = (exp(ax) / sqrt(ax)) * (q1+y*(q2+y*(q3+y*(q4
     1      +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
      endif
      return
      end
c
c**********************************************************************
c
      function bessi1(x)

      implicit none

      real bessi1, x
c      Returns the modified Bessel function I1(x) for any real x
      real ax
      double precision p1, p2, p3, p4, p5, p6, p7,
     1                 q1, q2, q3, q4, q5, q6, q7, q8, q9, y
      save  p1, p2, p3, p4, p5, p6, p7,
     1      q1, q2, q3, q4, q5, q6, q7, q8, q9

      data p1, p2, p3, p4, p5, p6, p7 / 0.5d0,0.87890594d0,0.51498869d0,
     1   0.15084934d0, 0.2658733d-1, 0.301532d-2, 0.32411d-3 /

      data q1, q2, q3, q4, q5, q6, q7, q8, q9/0.39894228d0,-0.398802d-1,
     1   -0.362018d-2,0.163801d-2,-0.1031555d-1, 0.2282967d-1,
     1   -0.2895312d-1, 0.1787654d-1, -0.420059d-2 /

      if (abs(x) .lt. 3.75) then
         y = (x / 3.75)**2
         bessi1 = x * (p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
         ax = abs(x)
         y = 3.75 / ax
         bessi1 = (exp(ax) / sqrt(ax)) * (q1+y*(q2+y*(q3+y*(q4
     1      +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
      endif
      return
      end

c
c***************************************************************************
c
      function bessk0 (x)

      implicit none

      real bessk0, x
c--uses bessi0
c--returns the modified Bessel function Ko(x) for any positive real x
      real bessi0
      doubleprecision p1, p2, p3, p4, p5, p6, p7
      doubleprecision q1, q2, q3, q4, q5, q6, q7, y
c--accumulate polynomials in double precision
      save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      data p1,p2,p3,p4,p5,p6,p7/-0.57721566D0,0.42278420D0,
     1   0.23069756D0,0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      data q1,q2,q3,q4,q5,q6,q7/1.25331414D0,-0.7832358D-1,0.2189568D-1,
     1   -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      if (x .lt. 2.0) then
c  Polynomial fit.
         y = x*x/4.0
         bessk0 =(-log(x/2.0)*bessi0(x))+(p1+y*(p2+y*(p3+
     *               y*(p4+y*(p5+y*(p6+y*p7))))))
      else
         y = (2.0/x)
         bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+
     *          y*(q4+ y*(q5+y*(q6+y*q7))))))
      endif
      return
      end
c
c***************************************************************************
c

      function bessk1 (x)

      implicit none

      real bessk1, x
c--uses bessi1
c--returns the modified Bessel function K1(x) for any positive real x
      real bessi1
      doubleprecision p1, p2, p3, p4, p5, p6, p7
      doubleprecision q1, q2, q3, q4, q5, q6, q7, y
c--accumulate polynomials in double precision
      save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      data p1,p2,p3,p4,p5,p6,p7/1.0D0,0.15443144D0,-0.67278579D0,
     1  -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      data q1,q2,q3,q4,q5,q6,q7/1.25331414D0,0.23498619D0,-0.3655620D-1,
     1   0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
      if (x .lt. 2.0) then
c  Polynomial fit.
         y = x*x/4.0
         bessk1 =(log(x/2.0)*bessi1(x))+(1.0/x)*(p1+y*(p2+
     *         y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
         y = 2.0/x
         bessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+
     *          y*(q4+ y*(q5+y*(q6+y*q7))))))
      endif
      return
      end



c
c**********************************************************************
c


