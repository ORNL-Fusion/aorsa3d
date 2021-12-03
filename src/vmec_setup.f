c
c***************************************************************************
c

      subroutine vmec_eval_b(r, z, phi,
     .     mnmax, ns, m, n, rho0, rc, rs, zc, zs,
     .     rrc, rrs, zrc, zrs,
     .     mnmaxnyq, mnyq, nnyq, rhoh0, buc, bus, bvc, bvs,
     .     rho, theta,
     .     br, bz, bphi)

      implicit none

      real r, z, phi
      integer mnmax, ns
      integer m(mnmax), n(mnmax)
      real rho0(ns)
      real rc(mnmax, ns), rs(mnmax, ns), zc(mnmax, ns), zs(mnmax, ns),
     .     rrc(mnmax, ns), rrs(mnmax, ns),
     .     zrc(mnmax, ns), zrs(mnmax, ns)
      integer mnmaxnyq
      integer mnyq(mnmaxnyq), nnyq(mnmaxnyq)
      real rhoh0(ns - 1)
      real buc(mnmaxnyq, ns - 1), bus(mnmaxnyq, ns - 1),
     .     bvc(mnmaxnyq, ns - 1), bvs(mnmaxnyq, ns - 1)
      real rho, theta
      real br, bz, bphi

      real ph, c, s
      real rm, zm, rb, zb
      real, parameter :: rerr = 1.0e-4
      real, parameter :: zerr = 1.0e-4
      real, parameter :: rhomin = 1.0e-6
      real, parameter :: rhomax = 1.1
      real, parameter :: drhomax = 0.1
      real, parameter :: dthetamax = 0.1
      integer i, itr
      integer, parameter :: itrmax = 10
      integer ir
      real h, rci, rsi, zci, zsi, rrci, rrsi, zrci, zrsi
      real f1, f2, j11, j21, j12, j22, delta
      real dx
      real buci, busi, bvci, bvsi
      real ru, rv, zu, zv, bu, bv

!     initial guess
      rm = 0.0
      zm = 0.0
      do i = 1, mnmax
         ph = -n(i) * phi
         c = cos(ph)
         s = sin(ph)
         rm = rm + c * rc(i, 1) + s * rs(i, 1)
         zm = zm + c * zc(i, 1) + s * zs(i, 1)
      end do
      theta = atan2(z - zm, r - rm)
      rb = 0.0
      zb = 0.0
      do i = 1, mnmax
         ph = m(i) * theta - n(i) * phi
         c = cos(ph)
         s = sin(ph)
         rb = rb + c * rc(i, ns) + s * rs(i, ns)
         zb = zb + c * zc(i, ns) + s * zs(i, ns)
      end do
      rho = ((r - rm)**2 + (z - zm)**2)
     .     / ((rb - rm)**2 + (zb - zm)**2)
      rho = sqrt(rho)
      if (rho > 1.0) rho = 1.0

!     search solution
      itr = 0
      do while (.true.)
         f1 = 0.0               ! R
         f2 = 0.0               ! Z
         j11 = 0.0              ! dRdrho
         j21 = 0.0              ! dZdrho
         j12 = 0.0              ! dRdtheta (ru)
         j22 = 0.0              ! dZdtheta (zu)
         ir = int(rho**2 * (ns - 1)) + 1
         if (ir > (ns - 1)) ir = ns - 1
         h = (rho - rho0(ir)) / (rho0(ir + 1) - rho0(ir))
         do i = 1, mnmax
            rci = (1.0 - h) * rc(i, ir) + h * rc(i, ir + 1)
            rsi = (1.0 - h) * rs(i, ir) + h * rs(i, ir + 1)
            zci = (1.0 - h) * zc(i, ir) + h * zc(i, ir + 1)
            zsi = (1.0 - h) * zs(i, ir) + h * zs(i, ir + 1)
            rrci = (1.0 - h) * rrc(i, ir) + h * rrc(i, ir + 1)
            rrsi = (1.0 - h) * rrs(i, ir) + h * rrs(i, ir + 1)
            zrci = (1.0 - h) * zrc(i, ir) + h * zrc(i, ir + 1)
            zrsi = (1.0 - h) * zrs(i, ir) + h * zrs(i, ir + 1)

            ph = m(i) * theta - n(i) * phi
            c = cos(ph)
            s = sin(ph)

            f1 = f1 + c * rci + s * rsi
            f2 = f2 + c * zci + s * zsi
            j11 = j11 + c * rrci + s * rrsi
            j21 = j21 + c * zrci + s * zrsi
            j12 = j12 - m(i) * s * rci + m(i) * c * rsi
            j22 = j22 - m(i) * s * zci + m(i) * c * zsi
         end do
         f1 = f1 - r
         f2 = f2 - z

         s = (f1 / rerr)**2 + (f2 / zerr)**2
         if ((s < 1.0) .or. (rho > rhomax)) exit
         if (itr > itrmax) exit

         delta = j11 * j22 - j12 * j21
         dx = -(j22 * f1 - j12 * f2) / delta
         if (dx < -drhomax) then
            dx = -drhomax
         else if (dx > drhomax) then
            dx = drhomax
         end if
         rho = rho + dx
         if (rho < rhomin) rho = rhomin

         dx = -(-j21 * f1 + j11 * f2) / delta
         if (dx < -dthetamax) then
            dx = -dthetamax
         else if (dx > dthetamax) then
            dx = dthetamax
         end if
         theta = theta + dx

         itr = itr + 1
      end do

!     eval B
      f1 = 0.0
      f2 = 0.0
      ru = 0.0
      rv = 0.0
      zu = 0.0
      zv = 0.0
      do i = 1, mnmax
         rci = (1.0 - h) * rc(i, ir) + h * rc(i, ir + 1)
         rsi = (1.0 - h) * rs(i, ir) + h * rs(i, ir + 1)
         zci = (1.0 - h) * zc(i, ir) + h * zc(i, ir + 1)
         zsi = (1.0 - h) * zs(i, ir) + h * zs(i, ir + 1)

         ph = m(i) * theta - n(i) * phi
         c = cos(ph)
         s = sin(ph)

         f1 = f1 + c * rci + s * rsi
         f2 = f2 + c * zci + s * zsi
         ru = ru - m(i) * s * rci + m(i) * c * rsi
         rv = rv + n(i) * s * rci - n(i) * c * rsi
         zu = zu - m(i) * s * zci + m(i) * c * zsi
         zv = zv + n(i) * s * zci - n(i) * c * zsi
      end do

!     B is on half mesh
      ir = int(rho**2 * (ns - 1) + 0.5)
      if (ir < 1) then
         ir = 1
      else if (ir > (ns - 2)) then
         ir = ns - 2
      end if
      h = (rho - rhoh0(ir)) / (rhoh0(ir + 1) - rhoh0(ir))
      bu = 0.0
      bv = 0.0
      do i = 1, mnmaxnyq
         buci = (1.0 - h) * buc(i, ir) + h * buc(i, ir + 1)
         bvci = (1.0 - h) * bvc(i, ir) + h * bvc(i, ir + 1)
         busi = (1.0 - h) * bus(i, ir) + h * bus(i, ir + 1)
         bvsi = (1.0 - h) * bvs(i, ir) + h * bvs(i, ir + 1)

         ph = mnyq(i) * theta - nnyq(i) * phi
         c = cos(ph)
         s = sin(ph)

         bu = bu + c * buci + s * busi
         bv = bv + c * bvci + s * bvsi
      end do

      br = bu * ru + bv * rv
      bz = bu * zu + bv * zv
      bphi = r * bv

      return
      end

c
c***************************************************************************
c

      subroutine vmec_read_head_695(lun, rmaxsurf, rminsurf, zmaxsurf,
     .     nfp, ns, mnmax, lasymlogical, ierflag, nbsets)

      implicit none

      integer lun
      real wb, wp, gamma, pfac, rmaxsurf, rminsurf, zmaxsurf
      integer nfp, ns, mpol, ntor, mnmax, itfsq, niter,
     .     lasymlogical, lreconlogical, ierflag
      integer imse, itse, nbsets, nobd, nextcur, nstore_seq

 100  format(6(1x,i11))
 110  format(3(1x,e23.15))

      read(lun, 110) wb, wp, gamma, pfac,
     .     rmaxsurf, rminsurf, zmaxsurf
      read(lun, 100) nfp, ns, mpol, ntor, mnmax, itfsq, niter,
     .     lasymlogical, lreconlogical, ierflag
      read(lun, 100) imse, itse, nbsets, nobd, nextcur, nstore_seq

      return
      end

c
c***************************************************************************
c

      subroutine vmec_setup(icontxt, myid, nproc, iwout, wout,
     .     nxmx, nymx, nphimx, nnodex, nnodey, nnodephi,
     .     capr, y, phi, rwleft, rwright, ytop, ybottom, rt, b0,
     .     bxn, byn, bzn, bmod, psi, rho, thetap, nfp)

      use netcdf

      implicit none

      integer*4 istatus, ncid, mndimid, sdimid, nfpid,
     .     xmid, xnid, rmncid, zmnsid, bumncid, bvmncid
      integer*4 mnmax_r, ns_r

      integer myid, nproc, id
      integer icontxt
      integer iwout
      character*128 wout
      integer nxmx, nymx, nphimx
      integer nnodex, nnodey, nnodephi
      real capr(nxmx), y(nymx), phi(nphimx)
      real rwleft, rwright, ytop, ybottom, rt, b0
      real bxn(nxmx, nymx, nphimx),
     .     byn(nxmx, nymx, nphimx),
     .     bzn(nxmx, nymx, nphimx),
     .     bmod(nxmx, nymx, nphimx)
      real psi(nxmx, nymx, nphimx),
     .     thetap(nxmx, nymx, nphimx),
     .     rho(nxmx, nymx, nphimx)
      integer nfp

      integer :: lun = 40
      real wout_version
      real rmaxsurf, rminsurf, zmaxsurf
      integer ns, mnmax, lasymlogical, ierflag, nbsets
      real, allocatable, dimension(:) :: nbfld
      character*128 mgridfile
      integer, allocatable, dimension(:) :: xm, xn
      real, allocatable, dimension(:, :) :: rmnc, rmns, zmnc, zmns,
     .     lmnc, lmns, bmnc, gmnc, bsubumnc, bsubvmnc, bsubsmns,
     .     bsupumnc, bsupumns, bsupvmnc, bsupvmns,
     .     currvmnc,
     .     rrmnc, rrmns, zrmnc, zrmns
      real, allocatable, dimension(:) :: rho0, rhoh0

      integer i, j, k, ngrid
      real dx, dy, dphi, phimax

      real, parameter :: pi = 3.141592654

      if (iwout .eq. 0) then
         istatus = nf90_open(TRIM(wout), NF90_NOWRITE, ncid)

         istatus = nf90_inq_dimid(ncid, 'mn', mndimid)
         istatus = nf90_inq_dimid(ncid, 's', sdimid)
         istatus = nf90_inquire_dimension(ncid, mndimid, len=mnmax_r)
         istatus = nf90_inquire_dimension(ncid, sdimid, len=ns_r)
         mnmax = mnmax_r
         ns = ns_r

         istatus = nf90_inq_varid(ncid, 'nfp', nfpid)
         istatus = nf90_get_var(ncid, nfpid, nfp)
      else
         open(unit=lun, file=wout, status="old", form="formatted")

         read(lun, "(15x,f4.2)") wout_version
         if (myid .eq. 0) then
            write(6, "(' wout version = ',f4.2)") wout_version
         end if

         call vmec_read_head_695(lun, rmaxsurf, rminsurf, zmaxsurf,
     .        nfp, ns, mnmax, lasymlogical, ierflag, nbsets)
         if ((ierflag .ne. 0) .and. (ierflag .ne. 4)) then
           write(6, *) "ierflag = ", ierflag
         end if

         if (nbsets > 0) then
            allocate(nbfld(nbsets))
            read(lun, 110) (nbfld(i), i = 1, nbsets)
            deallocate(nbfld)
         end if

         read(lun, "(a)") mgridfile
      end if

      allocate(xm(mnmax))
      allocate(xn(mnmax))
      allocate(rmnc(mnmax, ns))
      allocate(rmns(mnmax, ns))
      allocate(zmnc(mnmax, ns))
      allocate(zmns(mnmax, ns))
      allocate(lmnc(mnmax, ns - 1))
      allocate(lmns(mnmax, ns - 1))
      allocate(bmnc(mnmax, ns - 1))
      allocate(gmnc(mnmax, ns - 1))
      allocate(bsubumnc(mnmax, ns - 1))
      allocate(bsubvmnc(mnmax, ns - 1))
      allocate(bsubsmns(mnmax, ns - 1))
      allocate(bsupumnc(mnmax, ns - 1))
      allocate(bsupumns(mnmax, ns - 1))
      allocate(bsupvmnc(mnmax, ns - 1))
      allocate(bsupvmns(mnmax, ns - 1))
      allocate(currvmnc(mnmax, ns - 1))

      if (iwout .eq. 0) then
         istatus = nf90_inq_varid(ncid, 'xm', xmid)
         istatus = nf90_get_var(ncid, xmid, xm)
         istatus = nf90_inq_varid(ncid, 'xn', xnid)
         istatus = nf90_get_var(ncid, xnid, xn)
         istatus = nf90_inq_varid(ncid, 'rmnc', rmncid)
         istatus = nf90_get_var(ncid, rmncid, rmnc)
         istatus = nf90_inq_varid(ncid, 'zmns', zmnsid)
         istatus = nf90_get_var(ncid, zmnsid, zmns)
         istatus = nf90_inq_varid(ncid, 'bsupumnc', bumncid)
         istatus = nf90_get_var(ncid, bumncid, bsupumnc)
         istatus = nf90_inq_varid(ncid, 'bsupvmnc', bvmncid)
         istatus = nf90_get_var(ncid, bvmncid, bsupvmnc)

         istatus = nf90_close(ncid)
      else
 100     format(6(1x,i11))
 110     format(3(1x,e23.15))

         bsupumns = 0.0
         bsupvmns = 0.0

         if (lasymlogical .gt. 0) then
            do i = 1, mnmax
               read(lun, 100) xm(i), xn(i)
               read(lun, 110) rmnc(i, 1), zmns(i, 1),
     .             lmns(i, 1), bmnc(i, 1), gmnc(i, 1),
     .             bsubumnc(i, 1), bsubvmnc(i, 1), bsubsmns(i, 1),
     .             bsupumnc(i, 1), bsupvmnc(i, 1), currvmnc(i, 1),
     .             rmns(i, 1), zmnc(i, 1), lmnc(i, 1)
            end do
            do j = 1, ns - 1
               do i = 1, mnmax
                  read(lun, 110) rmnc(i, j + 1), zmns(i, j + 1),
     .                lmns(i, j), bmnc(i, j), gmnc(i, j),
     .                bsubumnc(i, j), bsubvmnc(i, j), bsubsmns(i, j),
     .                bsupumnc(i, j), bsupvmnc(i, j), currvmnc(i, j),
     .                rmns(i, j + 1), zmnc(i, j + 1), lmnc(i, j)
               end do
            end do
         else
            rmns = 0.0
            zmnc = 0.0
            lmnc = 0.0

            do i = 1, mnmax
               read(lun, 100) xm(i), xn(i)
               read(lun, 110) rmnc(i, 1), zmns(i, 1),
     .             lmns(i, 1), bmnc(i, 1), gmnc(i, 1),
     .             bsubumnc(i, 1), bsubvmnc(i, 1), bsubsmns(i, 1),
     .             bsupumnc(i, 1), bsupvmnc(i, 1), currvmnc(i, 1)
            end do
            do j = 1, ns - 1
               do i = 1, mnmax
                  read(lun, 110) rmnc(i, j + 1), zmns(i, j + 1),
     .                lmns(i, j), bmnc(i, j), gmnc(i, j),
     .                bsubumnc(i, j), bsubvmnc(i, j), bsubsmns(i, j),
     .                bsupumnc(i, j), bsupvmnc(i, j), currvmnc(i, j)
               end do
            end do
         end if

         close(lun)
      end if

      allocate(rho0(ns))
      allocate(rhoh0(ns - 1))
      allocate(rrmnc(mnmax, ns))
      allocate(rrmns(mnmax, ns))
      allocate(zrmnc(mnmax, ns))
      allocate(zrmns(mnmax, ns))

!     mesh
      dx = (rmaxsurf - rminsurf) * 0.1
      if (rwleft .eq. 0.0) rwleft = rminsurf - dx
      if (rwright .eq. 0.0) rwright = rmaxsurf + dx
      dx = zmaxsurf * 0.2
      if (ybottom .eq. 0.0) ybottom = -zmaxsurf - dx
      if (ytop .eq. 0.0) ytop = zmaxsurf + dx

      dx = (rwright - rwleft) / nnodex
      do i = 1, nnodex
         capr(i) = (i-1) * dx  + dx / 2.0 + rwleft
      end do
      dy = (ytop - ybottom) / nnodey
      do j = 1, nnodey
         y(j) = (j-1) * dy  + dy / 2.0 + ybottom
      end do
      phimax = 2.0 * pi / nfp
      dphi = phimax / nnodephi
      do k = 1, nnodephi
         phi(k) = (k-1) * dphi - phimax / 2.0
      end do

      do i = 1, ns
         rho0(i) = sqrt(real(i - 1) / (ns - 1))
      end do
      do i = 1, ns - 1
         rhoh0(i) = sqrt((i - 0.5) / (ns - 1))
      end do

!     rho derivatives
      do i = 2, ns - 1
         rrmnc(:, i) = (rmnc(:, i + 1) - rmnc(:, i - 1))
     .        / (rho0(i + 1) - rho0(i - 1))
         rrmns(:, i) = (rmns(:, i + 1) - rmns(:, i - 1))
     .        / (rho0(i + 1) - rho0(i - 1))
         zrmnc(:, i) = (zmnc(:, i + 1) - zmnc(:, i - 1))
     .        / (rho0(i + 1) - rho0(i - 1))
         zrmns(:, i) = (zmns(:, i + 1) - zmns(:, i - 1))
     .        / (rho0(i + 1) - rho0(i - 1))
      end do
      rrmnc(:, 1) = (rmnc(:, 2) - rmnc(:, 1)) / (rho0(2) - rho0(1))
      rrmns(:, 1) = (rmns(:, 2) - rmns(:, 1)) / (rho0(2) - rho0(1))
      zrmnc(:, 1) = (zmnc(:, 2) - zmnc(:, 1)) / (rho0(2) - rho0(1))
      zrmns(:, 1) = (zmns(:, 2) - zmns(:, 1)) / (rho0(2) - rho0(1))
      rrmnc(:, ns) = (rmnc(:, ns) - rmnc(:, ns - 1))
     .     / (rho0(ns) - rho0(ns - 1))
      rrmns(:, ns) = (rmns(:, ns) - rmns(:, ns - 1))
     .     / (rho0(ns) - rho0(ns - 1))
      zrmnc(:, ns) = (zmnc(:, ns) - zmnc(:, ns - 1))
     .     / (rho0(ns) - rho0(ns - 1))
      zrmns(:, ns) = (zmns(:, ns) - zmns(:, ns - 1))
     .     / (rho0(ns) - rho0(ns - 1))

!     invert (rho, theta) & evaluate B
      rho = 0.0
      psi = 0.0                 ! =rho^2
      thetap = 0.0
      bmod = 0.0
      bxn = 0.0
      byn = 0.0
      bzn = 0.0
      do k = 1, nnodephi
         do i = 1, nnodex
            do j = 1, nnodey
               ngrid = (k - 1) * nnodex * nnodey
     .               + (j - 1) * nnodex
     .               +  i
               id = mod(ngrid, nproc)
               if (id .eq. myid) then
                  call vmec_eval_b(capr(i), y(j), phi(k),
     .                mnmax, ns, xm, xn, rho0, rmnc, rmns, zmnc, zmns,
     .                rrmnc, rrmns, zrmnc, zrmns,
     .                mnmax, xm, xn, rhoh0, bsupumnc, bsupumns,
     .                bsupvmnc, bsupvmns,
     .                rho(i, j, k), thetap(i, j, k),
     .                bxn(i, j, k), byn(i, j, k), bzn(i, j, k))

                  rho(i, j, k) = rho(i, j, k)
                  psi(i, j, k) = rho(i, j, k)**2
                  bmod(i, j, k) = bxn(i, j, k)**2 + byn(i, j, k)**2
     .                 + bzn(i, j, k)**2
                  bmod(i, j, k) = sqrt(bmod(i, j, k))
                  if (bmod(i, j, k) .gt. 0.0) then
                     bxn(i, j, k) = bxn(i, j, k) / bmod(i, j, k)
                     byn(i, j, k) = byn(i, j, k) / bmod(i, j, k)
                     bzn(i, j, k) = bzn(i, j, k) / bmod(i, j, k)
                  end if
               end if
            end do
         end do
      end do

      deallocate(xm)
      deallocate(xn)
      deallocate(rmnc)
      deallocate(rmns)
      deallocate(zmnc)
      deallocate(zmns)
      deallocate(lmnc)
      deallocate(lmns)
      deallocate(bmnc)
      deallocate(gmnc)
      deallocate(bsubumnc)
      deallocate(bsubvmnc)
      deallocate(bsubsmns)
      deallocate(bsupumnc)
      deallocate(bsupumns)
      deallocate(bsupvmnc)
      deallocate(bsupvmns)
      deallocate(currvmnc)

      deallocate(rho0)
      deallocate(rhoh0)
      deallocate(rrmnc)
      deallocate(rrmns)
      deallocate(zrmnc)
      deallocate(zrmns)

      call blacs_barrier(icontxt, 'All')

      call dgsum2d(icontxt, 'All', ' ', nxmx * nymx, nnodephi,
     .     rho, nxmx * nymx, -1, -1)
      call dgsum2d(icontxt, 'All', ' ', nxmx * nymx, nnodephi,
     .     psi, nxmx * nymx, -1, -1)
      call dgsum2d(icontxt, 'All', ' ', nxmx * nymx, nnodephi,
     .     thetap, nxmx * nymx, -1, -1)
      call dgsum2d(icontxt, 'All', ' ', nxmx * nymx, nnodephi,
     .     bmod, nxmx * nymx, -1, -1)
      call dgsum2d(icontxt, 'All', ' ', nxmx * nymx, nnodephi,
     .     bxn, nxmx * nymx, -1, -1)
      call dgsum2d(icontxt, 'All', ' ', nxmx * nymx, nnodephi,
     .     byn, nxmx * nymx, -1, -1)
      call dgsum2d(icontxt, 'All', ' ', nxmx * nymx, nnodephi,
     .     bzn, nxmx * nymx, -1, -1)

      rt = (rwleft + rwright) * 0.5
      b0 = bmod(nnodex / 2, nnodey / 2, nnodephi / 2)

      call blacs_barrier(icontxt, 'All')

!     write transform.out
      if (myid .eq. 0) then
 310     format(1p6e12.4)
 309     format(10i10)

         open(unit=40,file='transform.out',status='unknown',
     .      form='formatted')

         write(40, 309) nnodex, nnodey, nnodephi
         write(40, 310) rt, b0
         write(40, 310) (((bxn(i,j,k), i = 1, nnodex), j = 1, nnodey),
     .                    k = 1, nnodephi)
         write(40, 310) (((byn(i,j,k), i = 1, nnodex), j = 1, nnodey),
     .                    k = 1, nnodephi)
         write(40, 310) (((bzn(i,j,k), i = 1, nnodex), j = 1, nnodey),
     .                    k = 1, nnodephi)
         write(40, 310) (((bmod(i,j,k), i = 1, nnodex), j = 1, nnodey),
     .                    k = 1, nnodephi)
         write(40, 310) (((psi(i,j,k), i = 1, nnodex), j = 1, nnodey),
     .                    k = 1, nnodephi)
         write(40, 310) (((rho(i,j,k), i = 1, nnodex), j = 1, nnodey),
     .                    k = 1, nnodephi)
         write(40, 310) (((thetap(i,j,k), i = 1, nnodex), j=1, nnodey),
     .                    k = 1, nnodephi)
         write(40, 309) nfp

         close(40)
      end if

      return
      end
