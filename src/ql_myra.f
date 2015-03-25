      module ql_myra_mod
      use qlsum_myra_mod
      include 'mpif.h'
      contains
      
c
c***************************************************************************
c

      subroutine ql_myra_write(!bqlavg, cqlavg, eqlavg, fqlavg, 
     .   wdot_inout, fx0_inout, fy0_inout, fz0_inout,
     .   vol, nrhodim, nnoderho, drho,
     .   dx, dy, r0, nxdim, nydim, nphidim, nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, nkdim1, nkdim2,
     .   mkdim1, mkdim2, nphidim1, nphidim2,
     .   xm, q, xn, xkt, lmax,
     .   xkxsav, xkysav, xkzsav, nzfun, ibessel,
     .   exk, eyk, ezk, capr,
     .   bxn, byn, bzn,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     .   xx, yy, zz, isigma, xnuomg, signb, psi, psilim, nboundary,
     .   myid, nproc, delta0, gradprlb, bmod, ndist,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara_cql, UmaxPara_cql, UPERP_cql, UPARA_cql, UPERP, UPARA,
     .   vc_mks_cql, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj,
     .   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     .   nkperp, lmaxdim, nzeta_wdot, theta_,
     .   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg, 
     .   n_bin, upshift, i_write, xkx_cutoff, xky_cutoff, xkz_cutoff)

*----------------------------------------------------------------------
*     This subroutine calculates wdot and the quasi-linear operator for 
*     a single species
*----------------------------------------------------------------------

      implicit none

      integer splitrank, nloops, partition, start, finish
      integer pstart, pfinish, status(MPI_STATUS_SIZE), ierr
      integer recv_size, remainder
      real, allocatable :: bql_store(:,:,:)
      real, allocatable :: cql_store(:,:,:)
      real, allocatable :: eql_store(:,:,:)
      real, allocatable :: fql_store(:,:,:)
      logical iam_root
      integer left_neighbor, right_neighbor
      real token
      real upara_test
      integer mi_max, mi_min

      logical ismine
      
      real u2, fnorm, f_cql
      
      real, dimension(:,:), allocatable :: DFDUPER0, DFDUPAR0
      
      integer n_theta_max, n_psi_max, i_psi_eq, n_bin, upshift
      integer n_theta_(n_psi_max), n_theta, ntheta_giv, i_write

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp, lmaxdim
      integer  i, j, k, l, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,
     .   icontxt, nboundary, nzeta_wdot
      integer dlen_, desc_amat(dlen_)
      integer nxdim, nydim, nphidim, nnodex, nnodey, nnodephi,
     .   nrhodim, nnoderho,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, nkdim1, nkdim2,
     .   mkdim1, mkdim2, nphidim1, nphidim2, isigma
      integer ftrap, ni0, mi0
      
      real dfdth, dfdupar_check, dfduper_check, dfdth_check
      real uperp0_grid, upara0_grid, zeta, eta, ai, bi, ci, di
      real dfduper0_intplt, dfdupar0_intplt,
     .     xkx_cutoff, xky_cutoff, xkz_cutoff
      
      real theta_(n_theta_max, n_psi_max), thegiv, dthetag
      real deriv, dtheta, dtau_ratio_giv, tau_bounce_giv
      real uperp0, upara0, xlamda, derivb, dtheta0, factor, upara_mi
      real duperp, dupara

      real zeffcd, clight, rt, b0, damp, r1, xm1, ceta, epsa, xmut2,
     .   eta0, xnexp, xjtild, signkz, cfit, c1, afit, yt, bmaxa, vphase,
     .   omgrf, c_ehst, akprl, xkprl, xkphi,
     .   a, rmaxa, rmina,
     .   wphase, xlnlam, vth
      real x(nxdim), y(nydim), argd, drho!, bratio


      real xm, omgc, omgp2,
     .     psi(nxdim, nydim, nphidim), psilim, alpha
      real capr(nxdim), xnuomg, signb,
     .     delta0, dx, dy, r0
      real gradprlb(nxdim, nydim, nphidim), bmod(nxdim, nydim, nphidim),
     .     bmod_mid(nxdim, nydim, nphidim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim),
     .        zz(nphidim1 : nphidim2, 1 : nphidim)
     
      complex wdoti, fx0i, fy0i
      

!      real bqlavg(nuper, nupar, nnoderho)
!      real cqlavg(nuper, nupar, nnoderho)
!      real eqlavg(nuper, nupar, nnoderho)
!      real fqlavg(nuper, nupar, nnoderho)
      
      real wdot_inout(nxdim, nydim, nphidim)
      real  fx0_inout(nxdim, nydim, nphidim)
      real  fy0_inout(nxdim, nydim, nphidim)
      real  fz0_inout(nxdim, nydim, nphidim)
      
      real wdot(nnodex, nnodey, nnodephi)
      real fx0(nnodex, nnodey, nnodephi)
      real fy0(nnodex, nnodey, nnodephi)
      real fz0(nnodex, nnodey, nnodephi)

      real wdot2d(nnodex, nnodey)
      real fx02d(nnodex, nnodey)
      real fy02d(nnodex, nnodey)
      real fz02d(nnodex, nnodey)
      
!      real count(0 : 10000, 1), sum_count

      real vol(nrhodim), dldbavg(nrhodim)
	


      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz


      complex cexpkxky
      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 :nphidim2),
     1        eyk(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 :nphidim2),
     1        ezk(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 :nphidim2)

      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2),
     .     xkzsav(nphidim1 : nphidim2)
      real q, xn(nxdim, nydim, nphidim), xkt(nxdim, nydim, nphidim)
      real bxn(nxdim, nydim, nphidim), byn(nxdim, nydim, nphidim),
     .     bzn(nxdim, nydim, nphidim)

      real uxx(nxdim, nydim, nphidim),
     .     uxy(nxdim, nydim, nphidim),
     .     uxz(nxdim, nydim, nphidim),
     .     uyx(nxdim, nydim, nphidim),
     .     uyy(nxdim, nydim, nphidim),
     .     uyz(nxdim, nydim, nphidim),
     .     uzx(nxdim, nydim, nphidim),
     .     uzy(nxdim, nydim, nphidim),
     .     uzz(nxdim, nydim, nphidim)

      real rho(nxdim, nydim, nphidim)

      integer :: n_psi_dim, nuper, nupar, n_psi, mi, ni

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: UPERP_cql(NUPER), UPARA_cql(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)

      real bql, cql, eql, fql

      complex, dimension(:),  allocatable :: b_sum
      complex, dimension(:),  allocatable :: c_sum
      complex, dimension(:),  allocatable :: e_sum
      complex, dimension(:),  allocatable :: f_sum
      
      complex, dimension(:),  allocatable :: wdot_sum
      complex, dimension(:),  allocatable :: sum_fx0
      complex, dimension(:),  allocatable :: sum_fy0
      
      real, dimension(:,:,:), allocatable :: factvol
      real, dimension(:,:), allocatable :: factvol2d

!      real, dimension(:,:,:), allocatable :: bqlvol
!      real, dimension(:,:), allocatable :: bqlvol2d

!      real, dimension(:,:,:), allocatable :: cqlvol
!      real, dimension(:,:), allocatable :: cqlvol2d

!      real, dimension(:,:,:), allocatable :: eqlvol
!      real, dimension(:,:), allocatable :: eqlvol2d

!      real, dimension(:,:,:), allocatable :: fqlvol
!      real, dimension(:,:), allocatable :: fqlvol2d

      real :: W, ENORM, ZSPEC, ASPEC, BMAG
      real :: UminPara,UmaxPara
      real :: UminPara_cql,UmaxPara_cql
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, vc_mks_cql, rho_a(n_psi_dim)
      real :: eps0, pi, emax, u0, u_0, u_, costh, costh0, psic

      parameter (eps0 = 8.85e-12)
      parameter (PI = 3.141592653597932384)


      integer :: nwork, ip
      logical :: has_work
      integer, dimension(nnodex*nnodey*nnodephi) :: i_table, j_table,
     .     k_table
      
      allocate( dfduper0(nuper, nupar) )
      allocate( dfdupar0(nuper, nupar) )
      
      allocate(b_sum(NUPAR) )
      allocate(c_sum(NUPAR) )
      allocate(e_sum(NUPAR) )
      allocate(f_sum(NUPAR) )
      
      allocate(wdot_sum(NUPER) )
      allocate(sum_fx0(NUPER) )
      allocate(sum_fy0(NUPER) )      
      
      allocate(factvol(nuper, nupar, nnoderho) )
      allocate(factvol2d(nuper, nupar) )

!      allocate(bqlvol(nuper, nupar, nnoderho) )
!      allocate(bqlvol2d(nuper, nupar) )

!      allocate(cqlvol(nuper, nupar, nnoderho) )
!      allocate(cqlvol2d(nuper, nupar) )

!      allocate(eqlvol(nuper, nupar, nnoderho) )
!      allocate(eqlvol2d(nuper, nupar) )

!      allocate(fqlvol(nuper, nupar, nnoderho) )
!      allocate(fqlvol2d(nuper, nupar) )

!     -------------------------------------
!efd  initialize allocatable arrays to zero
!     -------------------------------------
      wdot = 0.0
      fx0 = 0.0
      fy0 = 0.0
      fz0 = 0.0

      b_sum = 0.0
      c_sum = 0.0
      e_sum = 0.0
      f_sum = 0.0

      wdot_sum = 0.0
      sum_fx0 = 0.0
      sum_fy0 = 0.0

!      count = 0.0
      
c      factvol = 0.0
c      factvol2d = 0.0      

!      bqlvol = 0.0
!      bqlvol2d = 0.0

!      cqlvol = 0.0
!      cqlvol2d = 0.0

!      eqlvol = 0.0
!      eqlvol2d = 0.0

!      fqlvol = 0.0
!      fqlvol2d = 0.0




      W = omgrf
      ZSPEC = q / 1.6e-19

      if(ndist .eq. 0)then   !--Maxwellian--!

         do ni = 1, NUPER
            UPERP(ni) = (real(ni-1)/real(NUPER-1))
         end do

         UminPara = -1.0
         UmaxPara =  1.0

         do mi = 1, NUPAR
            UPARA(mi) = (-1.0 + 2. * (real(mi-1) / real(NUPAR-1)))
         end do

      else   !--non-Maxwellian--!

            vc_mks = vc_mks_cql
            UminPara = UminPara_cql
            UmaxPara = UmaxPara_cql

            do ni = 1, nuper
               uperp(ni) = uperp_cql(ni)
            end do

            do mi = 1, nupar
               upara(mi) = upara_cql(mi)
            end do

      end if


      do n = 1, nnoderho

         do mi = 1, nupar
            do ni = 1, nuper
	    

!               bqlvol(ni, mi, n) = 0.0
!               cqlvol(ni, mi, n) = 0.0
!               eqlvol(ni, mi, n) = 0.0
!               fqlvol(ni, mi, n) = 0.0


!               bqlavg(ni, mi, n) = 0.0
!               cqlavg(ni, mi, n) = 0.0
!               eqlavg(ni, mi, n) = 0.0
!               fqlavg(ni, mi, n) = 0.0

            end do
         end do
      end do



      nwork = 0
      do k=1,nnodephi
      do j=1,nnodey
         do i=1,nnodex
           has_work = (psi(i,j,k) .le. psilim .and. nboundary .eq. 1
     .                                        .or. nboundary .eq. 0)
           if (has_work) then
              nwork = nwork + 1
              i_table(nwork) = i
              j_table(nwork) = j
              k_table(nwork) = k
           endif
         enddo
      enddo
      enddo
      

*     -----------------------
*     Loop over spatial mesh:
*     -----------------------

      nloops = nwork
      partition=int(nloops/nproc)
      remainder=mod(nloops,nproc)
      splitrank=nproc-remainder
      if(myid.lt.splitrank)then
         start=1+myid*partition
         finish=start+partition-1
      else
         start=1+myid*partition+myid-splitrank
         finish=start+partition
      endif
      
      if(ndist.eq.1) then
             allocate(bql_store(start:finish,1:nuper,1:nupar),
     .         cql_store(start:finish,1:nuper,1:nupar),
     .         eql_store(start:finish,1:nuper,1:nupar),
     .         fql_store(start:finish,1:nuper,1:nupar))
      else
      endif


      do ip = start,finish
      
            i = i_table(ip)
            j = j_table(ip)
            k = k_table(ip)
	
            alpha = sqrt(2.0 * xkt(i, j, k) / xm)
            n = int(rho(i,j,k) / drho) + 1
	
            if (ndist .eq. 0) then 
	       vc_mks =  3.0 * alpha    !--Maxwellian only--!
	       vc_mks_cql = vc_mks
	       UminPara_cql = -1.0 
	       UmaxPara_cql = 1.0
	    end if
	    
            u0 = vc_mks / alpha
	    	    
            Emax = 0.5 * xm * vc_mks**2
            Enorm = Emax / 1.6e-19
            ASPEC = xm / 1.67e-27
            omgc  = q * bmod(i, j, k) / xm * signb
            omgp2 = xn(i,j,k) * q**2  / (eps0 * xm)

            BMAG = omgc * (xm / q)
	       
!	    bratio = bmod_mid(i,j,k) / bmod(i,j,k)
!	    if (bratio .gt. 1.0) bratio = 1.0
	    
	    duperp = (uperp(nuper) - uperp(1)) / (nuper - 1)
            dupara = (upara(nupar) - upara(1)) / (nupar - 1)
	    
!	    psic = 1.0 / bratio

!           -----------------------------------------------------
!           get CQL3D distribution function on the midplane:
!           used for Wdot only - not the quasilinear coefficients
!           -----------------------------------------------------

            if(ndist .eq. 0)then   !--Maxwellian--!
	
	       call maxwell_dist(u0, NUPAR, NUPER,
     .              UminPara, UmaxPara,
     .              UPERP, UPARA, DFDUPER, DFDUPAR)

            else   !--non-Maxwellian--!
               call maxwell_dist(u0, NUPAR, NUPER,
     .              UminPara, UmaxPara,
     .              UPERP, UPARA, DFDUPER, DFDUPAR)
               if (myid.eq.0) then
                  write(*,*) 'non-Maxwellian not implemented'
               end if
            end if
	
!           ----------------------------------
!           Loop over perpendicular velocities
!           ----------------------------------

            do ni = 1, nuper
	    
	    if (ndist .eq. 0)
     .          call QLSUM_MAXWELLIAN(ni, b_sum, c_sum, e_sum, f_sum,
     .               wdot_sum(ni), sum_fx0(ni), sum_fy0(ni), W, ZSPEC, 
     .               ASPEC, BMAG, lmax, ENORM, UminPara, UmaxPara,
     .               NUPAR, NUPER, UPERP, UPARA, DFDUPER, DFDUPAR,
     .               exk, eyk, ezk, nkdim1, nkdim2, mkdim1, mkdim2,
     .               nphidim1, nphidim2,
     .               nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .               uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .               uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .               uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .               nxdim, nydim, nphidim, xkxsav, xkysav, xkzsav,
     .               capr(i), xx, yy, zz, i, j, k,
     .               lmaxdim, ndist, nzeta_wdot,
     .               gradprlb(i,j,k), bmod(i,j,k), omgc,
     .               alpha, xm,
     .               upshift, xkx_cutoff, xky_cutoff, xkz_cutoff, rt)
     
	    if (ndist .eq. 1)
     .         call QLSUM_NON_MAXWELLIAN(ni, b_sum, c_sum, e_sum, f_sum,
     .               wdot_sum(ni), sum_fx0(ni), sum_fy0(ni), W, ZSPEC, 
     .               ASPEC, BMAG, lmax, ENORM, UminPara, UmaxPara,
     .               NUPAR, NUPER, UPERP, UPARA, DFDUPER, DFDUPAR,
     .               exk, eyk, ezk, nkdim1, nkdim2, mkdim1, mkdim2,
     .               nphidim1, nphidim2,
     .               nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .               uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .               uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .               uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .               nxdim, nydim, nphidim, xkxsav, xkysav, xkzsav,
     .               capr(i), xx, yy, zz, i, j, k,
     .               lmaxdim, ndist, nzeta_wdot,
     .               gradprlb(i,j,k), bmod(i,j,k), omgc,
     .               alpha, xm,
     .               upshift, xkx_cutoff, xky_cutoff, xkz_cutoff, rt)
	       	    	            	           	       
	       
!              -----------------------------
!              Loop over parallel velocities
!              -----------------------------

c$$$               upara_test = sqrt(1.0 - (UPERP(ni))**2)
c$$$	       mi_max = ceiling(upara_test*(nupar-1)/2 + (nupar+1)/2)
c$$$	       mi_min = floor(-1.0*upara_test*(nupar-1)/2
c$$$     .                        + (nupar+1)/2)
c$$$
c$$$               do mi = 1, nupar
c$$$	       
c$$$	          bql = 0.0 
c$$$                  cql = 0.0
c$$$                  eql = 0.0
c$$$                  fql = 0.0
c$$$
c$$$	          if(mi.ge.mi_min .and. mi.le.mi_max)then
c$$$
c$$$                  bql = 1.0 / (8. * emax * dupara)
c$$$     .                  * eps0 * omgp2(i,j,k) / omgrf * real(b_sum(mi))
c$$$
c$$$                  cql = 1.0 / (8. * emax * dupara)
c$$$     .                  * eps0 * omgp2(i,j,k) / omgrf * real(c_sum(mi))
c$$$
c$$$                  eql = 1.0 / (8. * emax * dupara)
c$$$     .                  * eps0 * omgp2(i,j,k) / omgrf * real(e_sum(mi))
c$$$
c$$$                  fql = 1.0 / (8. * emax * dupara)
c$$$     .                  * eps0 * omgp2(i,j,k) / omgrf * real(f_sum(mi))
c$$$
c$$$
c$$$
c$$$!                  if(bql .ne. 0.0)count(myid, 1) = count(myid, 1) + 1.0
c$$$		     
c$$$		     
c$$$                  if(n .le. nnoderho)then
c$$$		     
c$$$*                    ----------------------
c$$$*                    calculate midplane u's
c$$$*                    ----------------------
c$$$                     argd = uperp(ni)**2 * (1.-bratio) + upara(mi)**2
c$$$                     if (argd .le. 0.0) argd = 1.0e-06
c$$$
c$$$                     uperp0 = uperp(ni) * sqrt(bratio)
c$$$                     upara0 = sign(1.0, upara(mi)) * sqrt(argd)
c$$$			 
c$$$                     u_  = sqrt(uperp(ni)**2 + upara(mi)**2)
c$$$                     if (u_  .eq. 0.0) u_  = 1.0e-08
c$$$                     u_0 = u_
c$$$			 
c$$$
c$$$                     ni0 = int((uperp0 - uperp(1)) / duperp) + 1
c$$$                     mi0 = int((upara0 - upara(1)) / dupara) + 1
c$$$
c$$$*                    --------------------------------------
c$$$*                    bounce average and map to midplane u's
c$$$*                    --------------------------------------
c$$$                     if(ni0 .ge. 1 .and. ni0 .le. nuper .and. 
c$$$     .                     mi0 .ge. 1 .and. mi0 .le. nupar)then
c$$$     
c$$$     
c$$$     
c$$$     			costh0 = upara0 / u_0
c$$$			costh = upara(mi) / u_
c$$$			
c$$$			if(costh0 .eq. 0.0)costh0 = 1.0e-08
c$$$			
c$$$			upara_mi = upara(mi)
c$$$			if(upara_mi .eq. 0.0)
c$$$     .                        upara_mi = (upara(mi) + upara(mi+1)) / 2.0 
c$$$        
c$$$
c$$$c                        factor = abs(psic * upara0 / upara_mi)
c$$$c		    	 factor = abs(sqrt(psic))
c$$$                        factor = 1.0
c$$$ 
c$$$			
c$$$c                        factvol(ni0, mi0, n) = factvol(ni0, mi0, n)
c$$$c     .                      + dx * dy * capr(i) / r0 * factor 			
c$$$
c$$$                        bqlvol(ni0, mi0, n) = bqlvol(ni0, mi0, n)
c$$$     .                      + dx * dy * capr(i) / r0 * bql * factor 
c$$$
c$$$                        cqlvol(ni0, mi0, n) = cqlvol(ni0, mi0, n)
c$$$     .                      + dx * dy * capr(i) / r0 * cql * factor
c$$$     .                      * costh / costh0 / sqrt(psic)  
c$$$
c$$$                        eqlvol(ni0, mi0, n) = eqlvol(ni0, mi0, n)
c$$$     .                      + dx * dy * capr(i) / r0 * eql * factor 
c$$$     .                      * costh / costh0 / psic
c$$$
c$$$                        fqlvol(ni0, mi0, n) = fqlvol(ni0, mi0, n)
c$$$     .                      + dx * dy * capr(i) / r0 * fql * factor 
c$$$     .                      * (costh / costh0)**2 / psic**1.5
c$$$
c$$$                     end if
c$$$			
c$$$
c$$$
c$$$		  end if
c$$$
c$$$                  else
c$$$                  endif
c$$$
c$$$                  if(ndist.eq.1)then
c$$$		  
c$$$                    bql_store(ip,ni,mi) = bql
c$$$	    	    cql_store(ip,ni,mi) = cql
c$$$		    eql_store(ip,ni,mi) = eql
c$$$		    fql_store(ip,ni,mi) = fql
c$$$		  endif
c$$$
c$$$   
c$$$               end do

            end do

	       
	    wdoti = 0.0
	    fx0i = 0.0
	    fy0i = 0.0
	       
*           --------------------------------------------
*           integrate wdot over perpendicular velocities
*           --------------------------------------------	       
	       
	    do ni = 2, nuper - 1
	       wdoti = wdoti + wdot_sum(ni) * duperp
	       fx0i = fx0i + sum_fx0(ni) * duperp
	       fy0i = fy0i + sum_fy0(ni) * duperp
	    end do
	    
	      
	    wdoti = wdoti + 0.5 * wdot_sum(1)     * duperp
     .                    + 0.5 * wdot_sum(nuper) * duperp
     
	    fx0i = fx0i + 0.5 * sum_fx0(1)     * duperp
     .                  + 0.5 * sum_fx0(nuper) * duperp
     
	    fy0i = fy0i + 0.5 * sum_fy0(1)     * duperp
     .                  + 0.5 * sum_fy0(nuper) * duperp 
              
	      
     
            wdot(i,j,k) = - pi / 2.0 * eps0 * omgp2
     .                                     / omgrf * real(wdoti)
               
            fx0(i,j,k)  = - pi / 2.0 * eps0 * omgp2
     .                      / omgrf * real(fx0i) / (2.0 * omgrf)
            fy0(i,j,k)  = - pi / 2.0 * eps0 * omgp2
     .                       / omgrf * real(fy0i)/ (2.0 * omgrf)
     	            
            xkphi = xkzsav(nphi1) / capr(i)
            if(xkphi .eq. 0.0)xkphi = 1.0e-05
            fz0(i,j,k)  = xkphi / omgrf * wdot(i,j,k)

	
!     --------------------------------
!     end loop over i,j,k spatial points
!     --------------------------------
      end do


      if(ndist.eq.1 .and. i_write .ne. 0)then
      
         iam_root = (myid .eq. 0)
         left_neighbor = mod( myid - 1 + nproc, nproc)
         right_neighbor = mod( myid + 1, nproc )
         token = 1.0

         if(iam_root)then
	 
           open(unit=43,file='out_orbitrf.coef',  status='replace',
     .                                             form='formatted')
	 
	   do ip=start,finish
	      i=i_table(ip)
	      j=j_table(ip)
              k=k_table(ip)
	      do ni=1,nuper
	        do mi=1,nupar
		  if(abs(bql_store(ip,ni,mi)) .gt. 10**(-10)) then

		     write(43,102) i, j, ni, mi,
     .                 bql_store(ip,ni,mi), cql_store(ip,ni,mi),
     .                 eql_store(ip,ni,mi), fql_store(ip,ni,mi)
	          endif
		enddo
	      enddo
           enddo
	   
	   close(43)
	   
	   call MPI_SEND(token,1, MPI_REAL,right_neighbor,
     .             2, MPI_COMM_WORLD, ierr)
           call MPI_RECV(token,1, MPI_REAL,left_neighbor,
     &             2, MPI_COMM_WORLD, status, ierr)
         else
	 
           call MPI_RECV(token,1, MPI_REAL,left_neighbor,
     &             2, MPI_COMM_WORLD, status, ierr)
	   
           open(unit=43,file='out_orbitrf.coef',  status='old',
     .                    form='formatted', position='append')
	 
	   do ip=start,finish
	      i=i_table(ip)
	      j=j_table(ip)
	      do ni=1,nuper
	        do mi=1,nupar
		  if(abs(bql_store(ip,ni,mi)) .gt. 10**(-10)) then

		     write(43,102) i, j, ni, mi,
     .                 bql_store(ip,ni,mi), cql_store(ip,ni,mi),
     .                 eql_store(ip,ni,mi), fql_store(ip,ni,mi)
	          endif
		enddo
	      enddo
           enddo
	   
	   close(43)
	   
	   call MPI_SEND(token,1, MPI_REAL,right_neighbor,
     .             2, MPI_COMM_WORLD, ierr)
       
         endif	 
      endif

      call blacs_barrier(icontxt, 'All')

!     -------------------
!     Sum over processors
!     -------------------

      call dgsum2d(icontxt, 'All', ' ', nnodex * nnodey, nnodephi,
     .     wdot, nnodex * nnodey, -1, -1)

      call dgsum2d(icontxt, 'All', ' ', nnodex * nnodey, nnodephi,
     .     fx0, nnodex * nnodey, -1, -1)

      call dgsum2d(icontxt, 'All', ' ', nnodex * nnodey, nnodephi,
     .     fy0, nnodex * nnodey, -1, -1)

      call dgsum2d(icontxt, 'All', ' ', nnodex * nnodey, nnodephi,
     .     fz0, nnodex * nnodey, -1, -1)
     
c$$$      do n = 1, nnoderho
c$$$
c$$$         do mi = 1, nupar
c$$$            do ni = 1, nuper
c$$$c	        factvol2d(ni, mi) = factvol(ni, mi, n)
c$$$               bqlvol2d(ni, mi) = bqlvol(ni, mi, n)
c$$$               cqlvol2d(ni, mi) = cqlvol(ni, mi, n)
c$$$               eqlvol2d(ni, mi) = eqlvol(ni, mi, n)
c$$$               fqlvol2d(ni, mi) = fqlvol(ni, mi, n)
c$$$            end do
c$$$         end do
c$$$
c$$$c         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, factvol2d,
c$$$c     .      nuper, -1, -1)
c$$$         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, bqlvol2d,
c$$$     .      nuper, -1, -1)
c$$$         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, cqlvol2d,
c$$$     .      nuper, -1, -1)
c$$$         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, eqlvol2d,
c$$$     .      nuper, -1, -1)
c$$$         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, fqlvol2d,
c$$$     .      nuper, -1, -1)
c$$$     
c$$$
c$$$
c$$$
c$$$         do mi = 1, nupar
c$$$            do ni = 1, nuper
c$$$c	       factvol(ni, mi, n) = factvol2d(ni, mi)
c$$$               bqlvol(ni, mi, n) = bqlvol2d(ni, mi)
c$$$               cqlvol(ni, mi, n) = cqlvol2d(ni, mi)
c$$$               eqlvol(ni, mi, n) = eqlvol2d(ni, mi)
c$$$               fqlvol(ni, mi, n) = fqlvol2d(ni, mi)
c$$$            end do
c$$$         end do
c$$$
c$$$
c$$$      end do
      
      
*     ------------------------------------
*     Divide by volume element (passed in)
*     ------------------------------------
c$$$      do n = 1, nnoderho
c$$$         do mi = 1, nupar
c$$$            do ni = 1, nuper
c$$$
c$$$*           --------------
c$$$*           bounce average
c$$$*           --------------
c$$$            if (vol(n) .ne. 0.0)then
c$$$            bqlavg(ni, mi, n) = bqlvol(ni, mi, n) / vol(n) * dldbavg(n) 
c$$$            cqlavg(ni, mi, n) = cqlvol(ni, mi, n) / vol(n) * dldbavg(n) 
c$$$            eqlavg(ni, mi, n) = eqlvol(ni, mi, n) / vol(n) * dldbavg(n) 
c$$$            fqlavg(ni, mi, n) = fqlvol(ni, mi, n) / vol(n) * dldbavg(n)
c$$$	    end if
c$$$	    
c$$$c	    if (bqlavg(ni, mi, n) .lt. 0.0)bqlavg(ni, mi, n) = 0.0
c$$$c	    if (cqlavg(ni, mi, n) .lt. 0.0)cqlavg(ni, mi, n) = 0.0
c$$$c	    if (eqlavg(ni, mi, n) .lt. 0.0)eqlavg(ni, mi, n) = 0.0
c$$$c	    if (fqlavg(ni, mi, n) .lt. 0.0)fqlavg(ni, mi, n) = 0.0	    
c$$$	    
c$$$            end do
c$$$         end do
c$$$      end do
      
      deallocate( dfduper0 )
      deallocate( dfdupar0 )


      deallocate(b_sum)
      deallocate(c_sum)
      deallocate(e_sum)
      deallocate(f_sum)
      
      deallocate(wdot_sum)
      deallocate(sum_fx0)
      deallocate(sum_fy0)      
      
      deallocate(factvol)
      deallocate(factvol2d)      
      
c$$$      deallocate(bqlvol)
c$$$      deallocate(bqlvol2d)
c$$$
c$$$      deallocate(cqlvol)
c$$$      deallocate(cqlvol2d)
c$$$
c$$$      deallocate(eqlvol)
c$$$      deallocate(eqlvol2d)
c$$$
c$$$      deallocate(fqlvol)
c$$$      deallocate(fqlvol2d)
      
      if(ndist.eq.1) then
         deallocate(bql_store, cql_store, eql_store,  fql_store)
      endif

      wdot_inout(1:nnodex,1:nnodey,1:nnodephi)
     .     = wdot(1:nnodex,1:nnodey,1:nnodephi)
      fx0_inout(1:nnodex,1:nnodey,1:nnodephi)
     .     = fx0(1:nnodex,1:nnodey,1:nnodephi)
      fy0_inout(1:nnodex,1:nnodey,1:nnodephi)
     .     = fy0(1:nnodex,1:nnodey,1:nnodephi)
      fz0_inout(1:nnodex,1:nnodey,1:nnodephi)
     .     = fz0(1:nnodex,1:nnodey,1:nnodephi)
            
!      if (myid .eq. 0)then
!         sum_count = 0.0
!         do n = 0, 5000
!            write(43, 101)n, count(n, 1)
!	    sum_count = sum_count + count(n, 1)
!         end do
!	 write(43, *)"sum_count = ", sum_count
!      end if

      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (1i6, 1p8e12.4)
  102 format (4i6, 1p8e12.4)

      end subroutine ql_myra_write

      
c
c***************************************************************************
c


c$$$      subroutine wdot_qlcheck(wdot_check, 
c$$$     .   nnoderho, nrhodim,
c$$$     .   bqlavg, cqlavg, xm, omgrf, xktavg, ndist,
c$$$     .   nupar, nuper, n_psi,
c$$$     .   n_psi_dim, dfduper, dfdupar,
c$$$     .   UminPara_cql, UmaxPara_cql, UPERP_cql, UPARA_cql, UPERP, UPARA,
c$$$     .   vc_mks_cql, df_cql_uprp, df_cql_uprl, rhon, rho_a, myid,
c$$$     .   dldbavg)
c$$$
c$$$*     ------------------------------------------------------------------
c$$$*     This subroutine calculates the flux averaged wdot for checking QL
c$$$*     ------------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$
c$$$      integer nnoderho, nrhodim, n, m, ndist, myid
c$$$      real bqlavg(nuper, nupar, nnoderho)
c$$$      real cqlavg(nuper, nupar, nnoderho)
c$$$      real xktavg(nrhodim), alpha, e
c$$$      real ans, dfdth, dfdu, u, xm, omgrf, jacobian, upara_mi
c$$$      real, dimension(:,:), allocatable :: wdot_int
c$$$      real wdot_check(nrhodim), rhon(nrhodim)
c$$$      real dldbavg(nrhodim)
c$$$
c$$$      integer  :: n_psi_dim, nuper, nupar, n_psi, mi0, ni0
c$$$
c$$$      real :: UPERP(NUPER), UPARA(NUPAR)
c$$$      real :: UPERP_cql(NUPER), UPARA_cql(NUPAR)
c$$$      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
c$$$      real :: DFDUPER0, DFDUPAR0, UPARA0
c$$$      real :: W, ENORM
c$$$      real :: UminPara,UmaxPara
c$$$      real :: UminPara_cql,UmaxPara_cql
c$$$
c$$$      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
c$$$      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
c$$$
c$$$
c$$$      real :: vc_mks, vc_mks_cql, rho_a(n_psi_dim)
c$$$      real :: eps0, pi, emax, u0, dfdu0, dfdth0, u_0
c$$$
c$$$      parameter (eps0 = 8.85e-12)
c$$$      parameter (PI = 3.141592653597932384)
c$$$
c$$$      allocate(wdot_int(nuper, nupar) )
c$$$
c$$$!     ------------------------------------
c$$$!efd  initialize allocatable array to zero
c$$$!     ------------------------------------
c$$$
c$$$      wdot_int = 0.0
c$$$
c$$$
c$$$
c$$$      e = 1.6e-19
c$$$      W = omgrf
c$$$
c$$$      if(ndist .eq. 0)then   !--Maxwellian--!
c$$$
c$$$         do n = 1, NUPER
c$$$            UPERP(n) = (real(n-1)/real(NUPER-1))
c$$$         end do
c$$$
c$$$         UminPara = -1.0
c$$$         UmaxPara =  1.0
c$$$
c$$$         do m = 1, NUPAR
c$$$            UPARA(m) = (-1.0 + 2. * (real(m-1) / real(NUPAR-1)))
c$$$         end do
c$$$
c$$$      else   !--non-Maxwellian--!
c$$$
c$$$         vc_mks = vc_mks_cql
c$$$         UminPara = UminPara_cql
c$$$         UmaxPara = UmaxPara_cql
c$$$
c$$$         do n = 1, nuper
c$$$            uperp(n) = uperp_cql(n)
c$$$         end do
c$$$
c$$$         do m = 1, nupar
c$$$            upara(m) = upara_cql(m)
c$$$         end do
c$$$
c$$$      end if
c$$$
c$$$
c$$$*     -------------------
c$$$*     Loop over rho mesh:
c$$$*     -------------------
c$$$
c$$$      do n = 1, nnoderho
c$$$
c$$$         alpha = sqrt(2.0 * xktavg(n) / xm)
c$$$         if (ndist .eq. 0) vc_mks =  3.0 * alpha    !--Maxwellian only--!
c$$$         u0 = vc_mks / alpha
c$$$
c$$$         Emax = 0.5 * xm * vc_mks**2
c$$$         Enorm = Emax / 1.6e-19
c$$$
c$$$!        ------------------------------------------------
c$$$!        get CQL3D distribution function on the midplane
c$$$!        ------------------------------------------------
c$$$
c$$$         if(ndist .eq. 0)then   !--Maxwellian--!
c$$$	
c$$$            call maxwell_dist(u0, NUPAR, NUPER,
c$$$     .                 UminPara, UmaxPara,
c$$$     .                 UPERP, UPARA, DFDUPER, DFDUPAR)
c$$$
c$$$         else   !--non-Maxwellian--!
c$$$	
c$$$            call cql3d_dist(nupar, nuper, n_psi,
c$$$     .                 n_psi_dim, rho_a, rhon(n),
c$$$     .                 UminPara,UmaxPara,
c$$$     .                 df_cql_uprp, df_cql_uprl,
c$$$     .                 UPERP, UPARA, DFDUPER, DFDUPAR)
c$$$
c$$$
c$$$
c$$$         end if
c$$$	 
c$$$!        -----------------------------
c$$$!        loop over MIDPLANE velocities
c$$$!        -----------------------------	 
c$$$
c$$$         do ni0 = 1, nuper
c$$$            do mi0 = 1, nupar
c$$$
c$$$               dfdupar0 = dfdupar(ni0, mi0)
c$$$               dfduper0 = dfduper(ni0, mi0)
c$$$
c$$$               u_0 = sqrt(uperp(ni0)**2 + upara(mi0)**2)
c$$$               if (u_0 .eq. 0.0) u_0 = 1.0e-08
c$$$	       	
c$$$               dfdu0 = (uperp(ni0) * dfduper0 + upara(mi0) * dfdupar0)
c$$$     .                / u_0
c$$$               dfdth0 = upara(mi0) * dfduper0 - uperp(ni0) * dfdupar0
c$$$	       
c$$$               wdot_int(ni0, mi0) = (bqlavg(ni0, mi0, n) * dfdu0
c$$$     .                             + cqlavg(ni0, mi0, n) * dfdth0) / u_0 
c$$$
c$$$            end do
c$$$         end do
c$$$
c$$$
c$$$	
c$$$!        ---------------------------------------------------
c$$$!        Do velocity space integral over midplane velocities
c$$$!        ---------------------------------------------------
c$$$
c$$$         wdot_check(n) = 0.0
c$$$
c$$$
c$$$         call ugrate(wdot_int, uperp, upara, nuper, nupar, ans, myid,xm)
c$$$
c$$$
c$$$         wdot_check(n) = - 4.0 * pi * e * enorm / dldbavg(n) * ans
c$$$
c$$$
c$$$      end do
c$$$
c$$$
c$$$      deallocate(wdot_int)
c$$$
c$$$
c$$$      return
c$$$
c$$$ 1311 format(1p9e12.4)
c$$$ 1312 format(i10, 1p9e12.4)
c$$$ 1313 format(2i10, 1p9e12.4)
c$$$  100 format (1p8e12.4)
c$$$  101 format (2i10, 1p8e12.4)
c$$$
c$$$      end subroutine 
c$$$      
c$$$
c$$$c
c$$$c***************************************************************************
c$$$c
c$$$
c$$$
c$$$
c$$$      subroutine ugrate(f, uperp, upara, nuper, nupar, fint, myid, xm)
c$$$
c$$$      implicit none
c$$$
c$$$      integer nuper, nupar, ni, mi, myid
c$$$
c$$$      real uperp(nuper), upara(nupar), f(nuper, nupar)
c$$$      real favg, fint, duperp, dupara, xm
c$$$
c$$$      duperp = (uperp(nuper) - uperp(1)) / (nuper - 1)
c$$$      dupara = (upara(nupar) - upara(1)) / (nupar - 1)
c$$$
c$$$      fint = 0.0
c$$$
c$$$      do ni = 1, nuper - 1
c$$$         do mi = 1, nupar - 1
c$$$            favg = (f(ni, mi)   + f(ni+1, mi)
c$$$     .            + f(ni, mi+1) + f(ni+1, mi+1)) / 4.0
c$$$
c$$$            fint = fint + favg * uperp(ni) * duperp * dupara
c$$$
c$$$c            if (myid.eq.0 .and. ni .eq. 32 .and. mi .eq. 95) then
c$$$c	         write(6 ,1313)ni, mi, xm, favg, fint
c$$$c	         write(15,1313)ni, mi, xm, favg, fint
c$$$c	    end if
c$$$
c$$$         end do
c$$$      end do
c$$$
c$$$ 1313 format(2i10, 1p9e12.4)
c$$$
c$$$      return
c$$$      end subroutine  ugrate

c
c***************************************************************************
c

       end module ql_myra_mod
