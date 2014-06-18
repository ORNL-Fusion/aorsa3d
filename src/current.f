c
c***************************************************************************
c

      subroutine current_elect(xjpx, xjpy, xjpz, nxmx, nymx, nphimx,
     .   xjpx2d, xjpy2d, xjpz2d, gradprlb,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2,
     .   xm, q, xn, xkt, lmax,
     .   xkx, xky, xkz, nzfun, ibessel,
     .   x, y, phi,
     .   exk, eyk, ezk, capr,
     .   bxn, byn, bzn, bmod, dx, dy, dphi,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     .   xx, yy, zz, isigma, xnuomg, signb, psi, psilim,
     .   myid, nproc,
     .   iflag_gammab, delta0, xk0, damping, xkx_cutoff, xky_cutoff)

*-----------------------------------------------------------------------
*     This subroutine calculates the plasma current for a single species
*-----------------------------------------------------------------------

      implicit none

      logical ismine

      integer nproc, myid, ngrid, id, iflag_gammab
      real delta0, xk0, damping, xkx_cutoff, xky_cutoff
      integer  i, j, k, n, m, nphi, lmax, nzfun, ibessel, irnc, icnc
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,
     .   icontxt
      integer dlen_, desc_amat(dlen_)


      integer nxmx, nymx, nphimx,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2, isigma

      complex xx(nkdim1 : nkdim2, 1 : nxmx),
     .        yy(mkdim1 : mkdim2, 1 : nymx),
     .        zz(nphidim1 : nphidim2, 1 : nphimx)

      complex xjpx(nxmx, nymx, nphimx),
     .        xjpy(nxmx, nymx, nphimx),
     .        xjpz(nxmx, nymx, nphimx)

      real uxx(nxmx, nymx, nphimx),
     .     uxy(nxmx, nymx, nphimx),
     .     uxz(nxmx, nymx, nphimx),
     .     uyx(nxmx, nymx, nphimx),
     .     uyy(nxmx, nymx, nphimx),
     .     uyz(nxmx, nymx, nphimx),
     .     uzx(nxmx, nymx, nphimx),
     .     uzy(nxmx, nymx, nphimx),
     .     uzz(nxmx, nymx, nphimx)


      complex xjpx2d(nxmx, nymx),
     .        xjpy2d(nxmx, nymx),
     .        xjpz2d(nxmx, nymx)

      real bmod(nxmx, nymx, nphimx), signb, dum
      real xm, omgc, omgp2, eps0, xnuomg
      real capr(nxmx), dx, dy, q, xme, dphi


      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz


      complex cexpkxkykz
      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 :nphidim2),
     1        eyk(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 :nphidim2),
     1        ezk(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 :nphidim2)

      real xkx(nkdim1 : nkdim2),
     .     xky(mkdim1 : mkdim2),
     .     xkz(nphidim1 : nphidim2)
      real xn(nxmx, nymx, nphimx), xkt(nxmx, nymx, nphimx),
     .     x(nxmx), y(nymx), phi(nphimx)


      real bxn(nxmx, nymx, nphimx),
     .     byn(nxmx, nymx, nphimx),
     .     bzn(nxmx, nymx, nphimx)

      real gradprlb(nxmx, nymx, nphimx)
      real psi(nxmx, nymx, nphimx), psilim

      data eps0, xme / 8.85e-12, 9.11e-31 /

c      eps0 = 8.85e-12
c      zi = cmplx(0.0, 1.0)

      xjpx(:,:,:) = 0.0
      xjpy(:,:,:) = 0.0
      xjpz(:,:,:) = 0.0

*------------------------------------
*     Loop over mode numbers and mesh
*------------------------------------

      do k = 1, nnodephi
         do i = 1, nnodex
            do j = 1, nnodey

               ngrid = (k - 1) * nnodex * nnodey
     .               + (j - 1) * nnodex
     .               +  i

               id = mod(ngrid, nproc)
               if(id .eq. myid)then

                  omgc  = q * bmod(i, j, k) / xm * signb
                  omgp2 = xn(i,j,k) * q**2  / (eps0 * xm)


                  xjpx(i,j,k) = 0.0
                  xjpy(i,j,k) = 0.0
                  xjpz(i,j,k) = 0.0

                  if(psi(i,j,k) .le. psilim)then


                     do nphi = nphi1, nphi2
                        do n = nkx1, nkx2
                           do m = nky1, nky2


                              if (isigma .eq. 1)
     .                        call sigmah_stix_elect(i, j, n, m,
     .                        bmod(i,j,k), gradprlb(i,j,k),
     .                        xm, q, xn(i,j,k), xnuomg,
     .                        xkt(i,j,k), omgc, omgp2,
     .                        -lmax, lmax, nzfun, ibessel,
     .                        xkx(n), xky(m), xkz(nphi), capr(i),
     .                        bxn(i,j,k), byn(i,j,k), bzn(i,j,k),
     .                        uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .                        uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .                        uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .                        sigxx, sigxy, sigxz,
     .                        sigyx, sigyy, sigyz,
     .                        sigzx, sigzy, sigzz,
     .                        iflag_gammab, delta0, xk0,
     .                        damping, xkx_cutoff, xky_cutoff)


                              if (isigma .eq. 0)
     .                        call sigmac_stix(i, j, n, m,
     .                        xm, q, xn(i,j,k), xnuomg,
     .                        xkt(i,j,k), omgc, omgp2,
     .                        -lmax, lmax, nzfun, ibessel,
     .                        xkx(n), xky(m), xkz(nphi), capr(i),
     .                        bxn(i,j,k), byn(i,j,k), bzn(i,j,k),
     .                        uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .                        uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .                        uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .                        sigxx, sigxy, sigxz,
     .                        sigyx, sigyy, sigyz,
     .                        sigzx, sigzy, sigzz)

                              cexpkxkykz = xx(n,i) * yy(m,j) *zz(nphi,k)

                              xjpx(i,j,k) = xjpx(i,j,k) +
     .                              (sigxx * exk(n,m,nphi)
     1                             + sigxy * eyk(n,m,nphi)
     1                             + sigxz * ezk(n,m,nphi)) * cexpkxkykz
                              xjpy(i,j,k) = xjpy(i,j,k) +
     .                              (sigyx * exk(n,m,nphi)
     1                             + sigyy * eyk(n,m,nphi)
     1                             + sigyz * ezk(n,m,nphi)) * cexpkxkykz
                              xjpz(i,j,k) = xjpz(i,j,k) +
     .                              (sigzx * exk(n,m,nphi)
     1                             + sigzy * eyk(n,m,nphi)
     1                             + sigzz * ezk(n,m,nphi)) * cexpkxkykz



                           end do
                        end do
                     end do

                  end if
               end if

            end do
         end do
      end do


      call blacs_barrier(icontxt, 'All')


      do k = 1, nnodephi

         do i = 1, nnodex
            do j = 1, nnodey
               xjpx2d(i, j) = xjpx(i, j, k)
               xjpy2d(i, j) = xjpy(i, j, k)
               xjpz2d(i, j) = xjpz(i, j, k)
            end do
         end do

         call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx2d,
     .      nxmx, -1, -1)
         call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy2d,
     .      nxmx, -1, -1)
         call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz2d,
     .      nxmx, -1, -1)

         do i = 1, nnodex
            do j = 1, nnodey
               xjpx(i, j, k) = xjpx2d(i, j)
               xjpy(i, j, k) = xjpy2d(i, j)
               xjpz(i, j, k) = xjpz2d(i, j)
            end do
         end do

      end do


      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end

c
c***************************************************************************
c

      subroutine current(xjpx, xjpy, xjpz, nxmx, nymx, nphimx,
     .   xjpx2d, xjpy2d, xjpz2d, gradprlb,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2,
     .   xm, q, xn, xkt, lmax,
     .   xkx, xky, xkz, nzfun, ibessel,
     .   x, y, phi,
     .   exk, eyk, ezk, capr,
     .   bxn, byn, bzn, bmod, dx, dy, dphi,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     .   xx, yy, zz, isigma, xnuomg, signb, psi, psilim,
     .   myid, nproc, iflag_gammab, delta0, xk0)

*-----------------------------------------------------------------------
*     This subroutine calculates the plasma current for a single species
*-----------------------------------------------------------------------

      implicit none

      logical ismine

      integer nproc, myid, ngrid, id, iflag_gammab
      real delta0, xk0
      integer  i, j, k, n, m, nphi, lmax, nzfun, ibessel, irnc, icnc
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,
     .   icontxt
      integer dlen_, desc_amat(dlen_)


      integer nxmx, nymx, nphimx,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2, isigma

      complex xx(nkdim1 : nkdim2, 1 : nxmx),
     .        yy(mkdim1 : mkdim2, 1 : nymx),
     .        zz(nphidim1 : nphidim2, 1 : nphimx)

      complex xjpx(nxmx, nymx, nphimx),
     .        xjpy(nxmx, nymx, nphimx),
     .        xjpz(nxmx, nymx, nphimx)

      real uxx(nxmx, nymx, nphimx),
     .     uxy(nxmx, nymx, nphimx),
     .     uxz(nxmx, nymx, nphimx),
     .     uyx(nxmx, nymx, nphimx),
     .     uyy(nxmx, nymx, nphimx),
     .     uyz(nxmx, nymx, nphimx),
     .     uzx(nxmx, nymx, nphimx),
     .     uzy(nxmx, nymx, nphimx),
     .     uzz(nxmx, nymx, nphimx)


      complex xjpx2d(nxmx, nymx),
     .        xjpy2d(nxmx, nymx),
     .        xjpz2d(nxmx, nymx)

      real bmod(nxmx, nymx, nphimx), signb, dum
      real xm, omgc, omgp2, eps0, xnuomg
      real capr(nxmx), dx, dy, q, xme, dphi


      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz


      complex cexpkxkykz
      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 :nphidim2),
     1        eyk(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 :nphidim2),
     1        ezk(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 :nphidim2)

      real xkx(nkdim1 : nkdim2),
     .     xky(mkdim1 : mkdim2),
     .     xkz(nphidim1 : nphidim2)
      real xn(nxmx, nymx, nphimx), xkt(nxmx, nymx, nphimx),
     .     x(nxmx), y(nymx), phi(nphimx)


      real bxn(nxmx, nymx, nphimx),
     .     byn(nxmx, nymx, nphimx),
     .     bzn(nxmx, nymx, nphimx)

      real gradprlb(nxmx, nymx, nphimx)
      real psi(nxmx, nymx, nphimx), psilim

      data eps0, xme / 8.85e-12, 9.11e-31 /

c      eps0 = 8.85e-12
c      zi = cmplx(0.0, 1.0)

      xjpx(:,:,:) = 0.0
      xjpy(:,:,:) = 0.0
      xjpz(:,:,:) = 0.0

*------------------------------------
*     Loop over mode numbers and mesh
*------------------------------------

      do k = 1, nnodephi
         do i = 1, nnodex
            do j = 1, nnodey

               ngrid = (k - 1) * nnodex * nnodey
     .               + (j - 1) * nnodex
     .               +  i

               id = mod(ngrid, nproc)
               if(id .eq. myid)then

                  omgc  = q * bmod(i, j, k) / xm * signb
                  omgp2 = xn(i,j,k) * q**2  / (eps0 * xm)


                  xjpx(i,j,k) = 0.0
                  xjpy(i,j,k) = 0.0
                  xjpz(i,j,k) = 0.0

                  if(psi(i,j,k) .le. psilim)then


                     do nphi = nphi1, nphi2
                        do n = nkx1, nkx2
                           do m = nky1, nky2


                              if (isigma .eq. 1)
     .                        call sigmah_stix(i, j, n, m,
     .                        bmod(i,j,k), gradprlb(i,j,k),
     .                        xm, q, xn(i,j,k), xnuomg,
     .                        xkt(i,j,k), omgc, omgp2,
     .                        -lmax, lmax, nzfun, ibessel,
     .                        xkx(n), xky(m), xkz(nphi), capr(i),
     .                        bxn(i,j,k), byn(i,j,k), bzn(i,j,k),
     .                        uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .                        uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .                        uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .                        sigxx, sigxy, sigxz,
     .                        sigyx, sigyy, sigyz,
     .                        sigzx, sigzy, sigzz,
     .                        iflag_gammab, delta0, xk0)


                              if (isigma .eq. 0)
     .                        call sigmac_stix(i, j, n, m,
     .                        xm, q, xn(i,j,k), xnuomg,
     .                        xkt(i,j,k), omgc, omgp2,
     .                        -lmax, lmax, nzfun, ibessel,
     .                        xkx(n), xky(m), xkz(nphi), capr(i),
     .                        bxn(i,j,k), byn(i,j,k), bzn(i,j,k),
     .                        uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .                        uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .                        uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .                        sigxx, sigxy, sigxz,
     .                        sigyx, sigyy, sigyz,
     .                        sigzx, sigzy, sigzz)

                              cexpkxkykz = xx(n,i) * yy(m,j) *zz(nphi,k)

                              xjpx(i,j,k) = xjpx(i,j,k) +
     .                              (sigxx * exk(n,m,nphi)
     1                             + sigxy * eyk(n,m,nphi)
     1                             + sigxz * ezk(n,m,nphi)) * cexpkxkykz
                              xjpy(i,j,k) = xjpy(i,j,k) +
     .                              (sigyx * exk(n,m,nphi)
     1                             + sigyy * eyk(n,m,nphi)
     1                             + sigyz * ezk(n,m,nphi)) * cexpkxkykz
                              xjpz(i,j,k) = xjpz(i,j,k) +
     .                              (sigzx * exk(n,m,nphi)
     1                             + sigzy * eyk(n,m,nphi)
     1                             + sigzz * ezk(n,m,nphi)) * cexpkxkykz



                           end do
                        end do
                     end do

                  end if
               end if

            end do
         end do
      end do


      call blacs_barrier(icontxt, 'All')


      do k = 1, nnodephi

         do i = 1, nnodex
            do j = 1, nnodey
               xjpx2d(i, j) = xjpx(i, j, k)
               xjpy2d(i, j) = xjpy(i, j, k)
               xjpz2d(i, j) = xjpz(i, j, k)
            end do
         end do

         call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx2d,
     .      nxmx, -1, -1)
         call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy2d,
     .      nxmx, -1, -1)
         call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz2d,
     .      nxmx, -1, -1)

         do i = 1, nnodex
            do j = 1, nnodey
               xjpx(i, j, k) = xjpx2d(i, j)
               xjpy(i, j, k) = xjpy2d(i, j)
               xjpz(i, j, k) = xjpz2d(i, j)
            end do
         end do

      end do


      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end




