
c
c***************************************************************************
c
      subroutine checkx(myid, i, j, k, omgrf,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .   fdk, fek, ffk, u, v, w, xjx, xjy, xjz,
     .   nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .   nnodex, nnodey, nnodephi, nxmx, nymx, nphimx,
     .   xprime, yprime, phiprime)

      implicit none

      integer nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .    nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .    nnodex, nnodey, nnodephi, nxmx, nymx, nphimx,
     .    i, j, k, n, m, nphi, myid

      real xprime(nxmx), yprime(nymx), phiprime(nphimx)
      real eps0, omgrf

      complex xjx(nxmx, nymx, nphimx),
     .        xjy(nxmx, nymx, nphimx),
     .        xjz(nxmx, nymx, nphimx)

      complex fdk(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        fek(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        ffk(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2)

      complex term1, term2, term3


      complex xlhs, xrhs
      complex u(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        v(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        w(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2)
      complex zi


      zi = cmplx(0.0, 1.0)
      eps0 = 8.85e-12



      if (myid .eq. 0) then
         write(15, 910)
         write(15, 930)
         write(15, 920)
         write(15, 940)
         write(15, 920)

         write(6, 910)
         write(6, 930)
         write(6, 920)
         write(6, 940)
         write(6, 920)



         term1 = 0.0
         term2 = 0.0
         term3 = 0.0


         do nphi = nphi1, nphi2
            do n = nkx1, nkx2
               do m = nky1, nky2

                  term1 = term1 + fdk(n, m, nphi) * u(n, m, nphi)
                  term2 = term2 + fek(n, m, nphi) * v(n, m, nphi)
                  term3 = term3 + ffk(n, m, nphi) * w(n, m, nphi)

               end do
            end do
         end do

         xlhs = term1 + term2 + term3
         xrhs = -zi / omgrf / eps0 * xjx(i, j, k)
         write(15, 900) i, j, xprime(i), yprime(j), xlhs, xrhs
         write(6,  900) i, j, xprime(i), yprime(j), xlhs, xrhs
      endif


      return
  900 format(i5, i5, 1p9e12.3)
  940 format(1h , 4h   i, 5h    j, 4x, 6h x(i) , 6x, 6h y(j) ,
     1                    6x,  6hre lhs, 6x, 6him lhs, 6x,
     1                         6hre rhs, 6x, 6him rhs)
  910 format (1h1)
  920 format (1h0)
  930 format(3x, 14h u diagnostics)
      end

c
c***************************************************************************
c
      subroutine checky(myid, i, j, k, omgrf,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .   fgk, fak, fpk, u, v, w, xjx, xjy, xjz,
     .   nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .   nnodex, nnodey, nnodephi, nxmx, nymx, nphimx,
     .   xprime, yprime, phiprime)

      implicit none

      integer nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .    nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .    nnodex, nnodey, nnodephi, nxmx, nymx, nphimx,
     .    i, j, k, n, m, nphi, myid

      real xprime(nxmx), yprime(nymx), phiprime(nphimx)
      real eps0, omgrf

      complex xjx(nxmx, nymx, nphimx),
     .        xjy(nxmx, nymx, nphimx),
     .        xjz(nxmx, nymx, nphimx)

      complex fgk(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        fak(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        fpk(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2)

      complex term1, term2, term3


      complex xlhs, xrhs
      complex u(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        v(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        w(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2)
      complex zi


      zi = cmplx(0.0, 1.0)
      eps0 = 8.85e-12



      if (myid .eq. 0) then
         write(15, 910)
         write(15, 930)
         write(15, 920)
         write(15, 940)
         write(15, 920)

         write(6, 910)
         write(6, 930)
         write(6, 920)
         write(6, 940)
         write(6, 920)



         term1 = 0.0
         term2 = 0.0
         term3 = 0.0


         do nphi = nphi1, nphi2
            do n = nkx1, nkx2
               do m = nky1, nky2

                  term1 = term1 + fgk(n, m, nphi) * u(n, m, nphi)
                  term2 = term2 + fak(n, m, nphi) * v(n, m, nphi)
                  term3 = term3 + fpk(n, m, nphi) * w(n, m, nphi)

               end do
            end do
         end do

         xlhs = term1 + term2 + term3
         xrhs = -zi / omgrf / eps0 * xjy(i, j, k)
         write(15, 900) i, j, xprime(i), yprime(j), xlhs, xrhs
         write(6,  900) i, j, xprime(i), yprime(j), xlhs, xrhs
      endif


      return
  900 format(i5, i5, 1p9e12.3)
  940 format(1h , 4h   i, 5h    j, 4x, 6h x(i) , 6x, 6h y(j) ,
     1                    6x,  6hre lhs, 6x, 6him lhs, 6x,
     1                         6hre rhs, 6x, 6him rhs)
  910 format (1h1)
  920 format (1h0)
  930 format(3x, 14h v diagnostics)
      end

c
c***************************************************************************
c
      subroutine checkz(myid, i, j, k, omgrf,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .   frk, fqk, fsk, u, v, w, xjx, xjy, xjz,
     .   nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .   nnodex, nnodey, nnodephi, nxmx, nymx, nphimx,
     .   xprime, yprime, phiprime)

      implicit none

      integer nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .    nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .    nnodex, nnodey, nnodephi, nxmx, nymx, nphimx,
     .    i, j, k, n, m, nphi, myid

      real xprime(nxmx), yprime(nymx), phiprime(nphimx)
      real eps0, omgrf

      complex xjx(nxmx, nymx, nphimx),
     .        xjy(nxmx, nymx, nphimx),
     .        xjz(nxmx, nymx, nphimx)

      complex frk(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        fqk(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        fsk(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2)

      complex term1, term2, term3


      complex xlhs, xrhs
      complex u(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        v(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        w(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2)
      complex zi


      zi = cmplx(0.0, 1.0)
      eps0 = 8.85e-12



      if (myid .eq. 0) then
         write(15, 910)
         write(15, 930)
         write(15, 920)
         write(15, 940)
         write(15, 920)

         write(6, 910)
         write(6, 930)
         write(6, 920)
         write(6, 940)
         write(6, 920)



         term1 = 0.0
         term2 = 0.0
         term3 = 0.0


         do nphi = nphi1, nphi2
            do n = nkx1, nkx2
               do m = nky1, nky2

                  term1 = term1 + frk(n, m, nphi) * u(n, m, nphi)
                  term2 = term2 + fqk(n, m, nphi) * v(n, m, nphi)
                  term3 = term3 + fsk(n, m, nphi) * w(n, m, nphi)

               end do
            end do
         end do

         xlhs = term1 + term2 + term3
         xrhs = -zi / omgrf / eps0 * xjz(i, j, k)
         write(15, 900) i, j, xprime(i), yprime(j), xlhs, xrhs
         write(6,  900) i, j, xprime(i), yprime(j), xlhs, xrhs
      endif


      return
  900 format(i5, i5, 1p9e12.3)
  940 format(1h , 4h   i, 5h    j, 4x, 6h x(i) , 6x, 6h y(j) ,
     1                    6x,  6hre lhs, 6x, 6him lhs, 6x,
     1                         6hre rhs, 6x, 6him rhs)
  910 format (1h1)
  920 format (1h0)
  930 format(3x, 14h w diagnostics)
      end

c
c***************************************************************************
c

