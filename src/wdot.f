c
c***************************************************************************
c
      subroutine wl_fst2(xkperpl, xkperpr, xkprll, xkprlr,
     .   tn, tnp, tnpp, bn, zetan, lmaxdim,
     .   xkphi, xm, q, xkt, omgc, omgp2, omgrf,
     .   lmin, lmax, nzfun,
     .   dxdth, capr, xntau, btau,
     .   wxx, wxy, wxz, fyxx, fyxy, fyxz,
     .   wyx, wyy, wyz, fyyx, fyyy, fyyz,
     .   wzx, wzy, wzz, fyzx, fyzy, fyzz)

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs

      real xkperpl, xkperpr, xkprll, xkprlr,
     .     xkprl, gamleft, gamright,
     .     xm, q, xkt, omgc, omgp2, omgrf

      real akprl, gammab, rho, alpha, eps0

      real xkphi, dxdth, capr, xntau, btau

      complex gambar, gamtild
      complex wxx, wxy, wxz,
     1        wyx, wyy, wyz,
     2        wzx, wzy, wzz

      complex fyxx, fyxy, fyxz,
     1        fyyx, fyyy, fyyz,
     2        fyzx, fyzy, fyzz

      complex wxxl, wxyl, wxzl,
     1        wyxl, wyyl, wyzl,
     2        wzxl, wzyl, wzzl


      complex exil(0: lmaxdim), exilp(0: lmaxdim), exilpp(0: lmaxdim)

      complex zi, zieps0


      complex tn(-lmaxdim : lmaxdim),
     .        tnp(-lmaxdim : lmaxdim),
     .        tnpp(-lmaxdim : lmaxdim),
     .        zetan(-lmaxdim : lmaxdim),
     .        bn(-lmaxdim : lmaxdim)

      zi = cmplx(0.0, 1.0)
      eps0 = 8.85e-12

      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rho = alpha / omgc


      gamleft  = 0.5 * xkperpl**2 * rho**2
      gamright = 0.5 * xkperpr**2 * rho**2


c      gamtild = sqrt(gamleft * gamright)
c      The sqrt loses sign of k's so use instead
      gamtild = 0.5 * xkperpl * xkperpr * rho**2
      gambar = 0.5 * (gamleft + gamright)



c--   Calculate and store Bessel functions for l=0 to l = lmax (Il is ever in l):

      call besiexp2(gamtild, gambar, lmax, exil, exilp, exilpp, lmaxdim)

      wxx = 0.0
      wxy = 0.0
      wxz = 0.0
      wyx = 0.0
      wyy = 0.0
      wyz = 0.0
      wzx = 0.0
      wzy = 0.0
      wzz = 0.0

      fyxx = 0.0
      fyxy = 0.0
      fyxz = 0.0
      fyyx = 0.0
      fyyy = 0.0
      fyyz = 0.0
      fyzx = 0.0
      fyzy = 0.0
      fyzz = 0.0

      do l = lmin, lmax
         labs = abs(l)


         wxxl = - zieps0 * l**2 * exil(labs)/ gamtild * tn(l)
         wxyl = - eps0 * l * (gamright / gamtild * exil(labs)
     1              - exilp(labs)) * tn(l)
         wxzl = zieps0 * xkperpl * l * exil(labs) / gamleft * bn(l)



         wyxl =   eps0 * l * (gamleft / gamtild * exil(labs)
     1              - exilp(labs)) *  tn(l)
         wyyl = - zieps0 * (l**2 / gamtild * exil(labs)
     1                + 2. * gamtild * exil(labs)
     1                - 2. * gambar * exilp(labs) ) *  tn(l)
         wyzl = - eps0 * ( xkperpl * exil(labs)
     1                           - xkperpr * exilp(labs) ) * bn(l)



         wzxl = zieps0 * xkperpr * l * exil(labs)/ gamright * bn(l)
         wzyl =   eps0 * ( xkperpr * exil(labs)
     1                           - xkperpl * exilp(labs) ) * bn(l)
         wzzl = zieps0 * exil(labs) * zetan(l) * tnp(l)


         wxx = wxx + wxxl
         wxy = wxy + wxyl
         wxz = wxz + wxzl
         wyx = wyx + wyxl
         wyy = wyy + wyyl
         wyz = wyz + wyzl
         wzx = wzx + wzxl
         wzy = wzy + wzyl
         wzz = wzz + wzzl

         fyxx = fyxx + l / (2. * omgrf) * wxxl
         fyxy = fyxy + l / (2. * omgrf) * wxyl
         fyxz = fyxz + l / (2. * omgrf) * wxzl
         fyyx = fyyx + l / (2. * omgrf) * wyxl
         fyyy = fyyy + l / (2. * omgrf) * wyyl
         fyyz = fyyz + l / (2. * omgrf) * wyzl
         fyzx = fyzx + l / (2. * omgrf) * wzxl
         fyzy = fyzy + l / (2. * omgrf) * wzyl
         fyzz = fyzz + l / (2. * omgrf) * wzzl

      end do

      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
      end
c
c***************************************************************************
c
      subroutine wl_fst(xkperpl, xkperpr, xkprll, xkprlr,
     .   tn, tnp, tnpp, bn, zetan, lmaxdim,
     .   xkphi, xm, q, xkt, omgc, omgp2,
     .   lmin, lmax, nzfun,
     .   dxdth, capr, xntau, btau,
     .   wxx, wxy, wxz,
     .   wyx, wyy, wyz,
     .   wzx, wzy, wzz)

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs

      real xkperpl, xkperpr, xkprll, xkprlr,
     .     xkprl, gamleft, gamright,
     .     xm, q, xkt, omgc, omgp2

      real akprl, gammab, rho, alpha, eps0, omgrf, v0i

      real xkphi, dxdth, capr, xntau, btau

      complex gambar, gamtild
      complex wxx, wxy, wxz,
     1        wyx, wyy, wyz,
     2        wzx, wzy, wzz

      complex wxxl, wxyl, wxzl,
     1        wyxl, wyyl, wyzl,
     2        wzxl, wzyl, wzzl


      complex exil(0: lmaxdim), exilp(0: lmaxdim), exilpp(0: lmaxdim)

      complex zi, zieps0

      common/sigcom/zi, eps0, v0i, omgrf
      complex tn(-lmaxdim : lmaxdim),
     .        tnp(-lmaxdim : lmaxdim),
     .        tnpp(-lmaxdim : lmaxdim),
     .        zetan(-lmaxdim : lmaxdim),
     .        bn(-lmaxdim : lmaxdim)



      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rho = alpha / omgc


      gamleft  = 0.5 * xkperpl**2 * rho**2
      gamright = 0.5 * xkperpr**2 * rho**2


c      gamtild = sqrt(gamleft * gamright)
c      The sqrt loses sign of k's so use instead
      gamtild = 0.5 * xkperpl * xkperpr * rho**2
      gambar = 0.5 * (gamleft + gamright)



c--   Calculate and store Bessel functions for l=0 to l = lmax (Il is ever in l):

      call besiexp2(gamtild, gambar, lmax, exil, exilp, exilpp, lmaxdim)

      wxx = 0.0
      wxy = 0.0
      wxz = 0.0

      wyx = 0.0
      wyy = 0.0
      wyz = 0.0

      wzx = 0.0
      wzy = 0.0
      wzz = 0.0

      do l = lmin, lmax
         labs = abs(l)


         wxxl = - zieps0 * l**2 * exil(labs)/ gamtild * tn(l)
         wxyl = - eps0 * l * (gamright / gamtild * exil(labs)
     1              - exilp(labs)) * tn(l)
         wxzl = zieps0 * xkperpl * l * exil(labs) / gamleft * bn(l)



         wyxl =   eps0 * l * (gamleft / gamtild * exil(labs)
     1              - exilp(labs)) *  tn(l)
         wyyl = - zieps0 * (l**2 / gamtild * exil(labs)
     1                + 2. * gamtild * exil(labs)
     1                - 2. * gambar * exilp(labs) ) *  tn(l)
         wyzl = - eps0 * ( xkperpl * exil(labs)
     1                           - xkperpr * exilp(labs) ) * bn(l)



         wzxl = zieps0 * xkperpr * l * exil(labs)/ gamright * bn(l)
         wzyl =   eps0 * ( xkperpr * exil(labs)
     1                           - xkperpl * exilp(labs) ) * bn(l)
         wzzl = zieps0 * exil(labs) * zetan(l) * tnp(l)


         wxx = wxx + wxxl
         wxy = wxy + wxyl
         wxz = wxz + wxzl
         wyx = wyx + wyxl
         wyy = wyy + wyyl
         wyz = wyz + wyzl
         wzx = wzx + wzxl
         wzy = wzy + wzyl
         wzz = wzz + wzzl

      end do

      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
      end
c
c***************************************************************************
c
      subroutine estix(ex, ey, ez, ealpp, ebetp, ebp,
     .   bx, by, bz,
     .   xkx, xky, xkphi,
     .   xkprl, xkperp)

      implicit none

      real bx, by, bz, sqx, xkprl, xkalp, xkbet, xkperp,
     .   cosalp, sinalp, xkx, xky, xkphi

      real uxx, uxy, uxz,
     1     uyx, uyy, uyz,
     2     uzx, uzy, uzz

      complex ex, ey, ez, ealpp, ebetp, ebp, ealp, ebet, eb

      sqx = sqrt(1.0 - bx**2)

      xkprl = bx * xkx + by * xky + bz * xkphi
      xkalp = (xkx * (1.0 - bx**2) - xky * bx * by
     .       - xkphi * bx * bz) / sqx
      xkbet = (xky * bz - xkphi * by) / sqx

      xkperp = sqrt(xkalp**2 + xkbet**2)


      cosalp = xkalp / xkperp
      sinalp = xkbet / xkperp


      uxx =   sqx
      uxy = - bx * by / sqx
      uxz = - bx * bz / sqx
      uyx =   0.0
      uyy =   bz / sqx
      uyz = - by / sqx
      uzx =   bx
      uzy =   by
      uzz =   bz

      ealp = uxx * ex + uxy * ey + uxz * ez
      ebet = uyx * ex + uyy * ey + uyz * ez
      eb   = uzx * ex + uzy * ey + uzz * ez

      ealpp =   cosalp * ealp + sinalp * ebet
      ebetp = - sinalp * ealp + cosalp * ebet
      ebp   = eb

      return
      end
c
c***************************************************************************
c
      subroutine besiexp2(gamtild, gambar, lmax, expbes, expbesp,
     1   expbespp, lmaxdim)
c*****   calculates exp(-gambar) times Il(gamtild)
c*****   (and derivatives) of order up to lmax

      implicit none

      integer lmax, nmax, ier, l, lmaxdim
      real gammod
      complex gamtild, gambar, expbes(0:lmaxdim), expbesp(0:lmaxdim),
     1   expbespp(0:lmaxdim),
     1   b(100), xil(0:99), xilp(0:99), xilpp(0:99), exgam

      exgam = cexp(-gambar)
      gammod = cabs(gamtild)

      if(gammod .le. 700.)then
         nmax = lmax + 1

         call besic(gamtild, nmax, b, ier)
         if(ier .ne. 0)write(15,100) ier

         do l = 0, lmax
            xil(l) = b(l+1)
         end do

         do l = 0, lmax
           if(l .eq. 0) xilp(0) = xil(1)
           if(l .ne. 0) xilp(l) = xil(l-1) - l / gamtild * xil(l)
           xilpp(l) = - xilp(l) / gamtild
     1                 + (1. + (l / gamtild)**2)* xil(l)

           expbes(l) =  exgam * xil(l)
           expbesp(l) = exgam * xilp(l)
           expbespp(l) = exgam * xilpp(l)
         end do
      end if

      if(gammod .gt. 700.)then
         do l = 0, lmax
            call bes_asymp(gamtild, gambar, l, expbes(l), expbesp(l))
            expbespp(l) = - expbesp(l) / gamtild
     1             + (1. + (l / gamtild)**2)* expbes(l)
         end do
      end if

  100 format('ier = ', i5, 'besic failed')
  101 format(1p8e12.4)
  102 format(i10, 1p8e12.4)
      return
      end
c
c***************************************************************************
c
      subroutine bes_asymp(z, y, n, exil, exilp)

      implicit none

      integer mu, n
      real pi
      complex z, y, exil, exilp
      data pi/3.141592654/

      mu = 4 * n**2
      exil =  1.0 / csqrt(2.0 * pi * z) * cexp(z - y)
     1   * (1.0
     1   - (mu - 1)/(8.0 * z)
     1   + (mu - 1) * (mu - 9) / (2.0 * (8.0 * z)**2)
     1   - (mu - 1) * (mu - 9) * (mu - 25) / (6.0 * (8.0 * z)**3)  )

      exilp = 1.0 / csqrt(2.0 * pi * z) * cexp(z - y)
     1   * (1.0
     1   - (mu + 3)/(8.0 * z)
     1   + (mu - 1) * (mu + 15) / (2.0 * (8.0 * z)**2)
     1   - (mu - 1) * (mu - 9) * (mu + 35) / (6.0 * (8.0 * z)**3)  )


  102 format(i10, 1p8e12.4)
      return
      end
c
c***************************************************************************
c
      subroutine fluxavg(f, favg, rho, nxdim, nydim, nphidim, nrhodim,
     .   nnodex, nnodey, nnodephi, nnoderho, drho, dx, dy, dphi, capr,
     .   fvol, vol)

      implicit none

      integer nxdim, nydim, nphidim, nrhodim, nnodex, nnodey, nnodephi,
     .   nnoderho
      integer n, i, j, k

      real f(nxdim, nydim, nphidim), favg(nrhodim),
     .    rho(nxdim, nydim, nphidim),
     .   drho, dx, dy, dphi, fvol(nrhodim), vol(nrhodim), capr(nxdim)


      do n = 1, nnoderho
          fvol(n) = 0.0
          vol(n) = 0.0
      end do


      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nnodephi
               n = int(rho(i,j,k) / drho) + 1
               if(n .ge. 1 .and. n .le. nnoderho)then

                  fvol(n) = fvol(n) + dx * dy * capr(i) * dphi *f(i,j,k)
                  vol(n) =   vol(n) + dx * dy * capr(i) * dphi

               end if
            end do
         end do
      end do


      do n = 1, nnoderho
         if(vol(n) .ne. 0.0)
     .     favg(n) = fvol(n) / vol(n)
c         write(6, 100)n, fvol(n), vol(n), favg(n)
      end do



  100 format (1i10, 1p8e12.4)
  101 format (10i10)
      return
      end
c
c***************************************************************************
c
