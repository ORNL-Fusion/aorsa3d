c
c***************************************************************************
c
       subroutine delta_(i, j, n, m,
     .   delta_x, delta_y, delta_z,
     .   bmod, gradprlb,
     .   xm, q, xn, xnuomg,
     .   xkt, omgc, omgp2,
     .   lmin, lmax, nzfun, ibessel,
     .   xkxsav, xkysav, xkzsav, capr,
     .   bx,    by,    bz,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz)

*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     and uses the old left handed U matrix in paper
*     ---------------------------------------------------------

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, sgn_kprl
      real dakbdkb, xnuomg, xkprl_cutoff
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

      complex delta_x, delta_y, delta_z
      complex delta_xl, delta_yl, delta_zl

      parameter (lmaxdim = 99)

      complex xil(0: lmaxdim), xilp(0: lmaxdim)
      complex exil(0: lmaxdim), exilp(0: lmaxdim),
     .     exilovergam(0: lmaxdim)

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

      if (xkprl  .eq. 0.0) xkprl  = 1.0e-08
      if (xkperp .eq. 0.0) xkperp = 1.0e-08

      sgn_kprl = sign(1.0, xkprl)
      akprl = abs(xkprl)


      gamma = 0.5 * xkperp**2 * rho**2
      rgamma = real(gamma)


      if(rgamma .ge. 1.0e-08)
     .     call besiexp(gamma, lmax, exil, exilp, lmaxdim, exilovergam)

      if(rgamma .lt. 1.0e-08)
     .     call bes_expand(gamma, lmax, exil, exilp, lmaxdim,
     .     exilovergam)

      delta_x = 0.0
      delta_y = 0.0
      delta_z = 0.0

      zeta0 = omgrfc / (xkprl * alpha)

      do l = lmin, lmax
         labs = abs(l)



         zetal = (omgrfc - l * omgc) / (xkprl * alpha)

         gammab = abs(l * omgc / (2.0 * alpha * xkprl**2)
     .        * gradprlb / bmod)


         gammab = 0.0
c         if(abs(gammab) .gt. 1000.0) gammab = 1000.0
         if(abs(gammab) .lt. .01)gammab = .01



         if (nzfun .eq. 0) call z_approx(sgn_kprl, zetal, 0.0,
     .                                                     z0, z1, z2)
         if (nzfun .eq. 1) call z_approx(sgn_kprl, zetal, gammab,
     .                                                     z0, z1, z2)
         if (nzfun .eq. 2) call z_smithe(sgn_kprl, zetal, gammab,
     .                                                     z0, z1, z2)
         if (nzfun .eq. 3) call z_table(sgn_kprl, zetal, gammab,
     .                                               .001, z0, z1, z2)


         delta_xl = - zieps0 / q * omgp2 / omgrf**2 * zeta0
     .        * xkperp * l * exilovergam(labs) * omgrf / omgc * z0
         
         delta_yl = - zieps0 / q * omgp2 / omgrf**2 * zeta0
     .        * zi * xkperp * omgrf /omgc * (exilp(labs) - exil(labs))
     .        * z0	
         
         delta_zl = - zieps0 / q * omgp2 / omgrf**2 * zeta0
     .        * xkprl * 2.0 * zeta0 * exil(labs) * z1	
         
         delta_x = delta_x + delta_xl 
         delta_y = delta_y + delta_yl 
         delta_z = delta_z + delta_zl 	    
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

      subroutine sigmad_cql3d(i, j, k, n, m, rho, rho_a,
     .   gradprlb, bmod, bmod0,
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
     .   delta0, ndist, nupar, nuper, n_psi,
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .   nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     .   i_sav, j_sav, k_sav, upshift, damping, xkx_cutoff, xky_cutoff,
     .   xkz_cutoff, rt, nkx2, nky2)


*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     No rotation is made.  Result is in the Stix frame.
*     ---------------------------------------------------------

      implicit none
      
      real, dimension(:,:), allocatable :: DFDUPER0, DFDUPAR0

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, k, n, m, nphi, ndist, myid, iflag, nproc, ni, mi,
     .    i_sav, j_sav, k_sav, ni0, mi0, upshift, nkx2, nky2
      integer n_upper, n_lower, m_upper, m_lower     
      
     
      real u0, u2, fnorm, f_cql, rt
      real fnmax, fmmax, nperp
      real fna, fma, fstepn, fstepm, fstep

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, y, sgn_kprl, reson, duperp, dupara
      real xkprl_eff0
      real dzetal(lmin:lmax), descrim
      real dakbdkb, xnuomg, gradprlb, bmod, bmod0, nu_coll
      real akprl,  rho, alpha, eps0, omgrf, v0i, emax
      real gammab(lmin:lmax), gamma_coll(lmin:lmax)
      real a, b, xnurf, pi, delta0, rhol
      real bx, by, bz, bratio, denom
      real dfdth, dfdupar_check, dfduper_check, dfdth_check
      real uperp0_grid, upara0_grid, zeta, eta, ai, bi, ci, di
      real dfduper0_intplt, dfdupar0_intplt

      real xkxsav, xkysav, xkzsav, capr, argd
      real xkphi
      real xkalp, xkbet, xk0, rgamma, xkx_cutoff, xky_cutoff,
     .     xkz_cutoff, damping, kr, step

      complex zi, zfunct, fzeta, omgrfc

      complex zfunct0, zeta0, sig3cold, z0, z1, z2, dz0, dz1, dz2

      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0_a, sig1_a, sig2_a, sig3_a, sig4_a, sig5_a
      complex sig0_h, sig1_h, sig2_h, sig3_h, sig4_h, sig5_h

      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
      logical :: l_interp   !new
      logical :: l_first  !new
      integer :: nkperp
      real :: kperp_max !new


      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz

      real uxx, uxy, uxz,
     1     uyx, uyy, uyz,
     2     uzx, uzy, uzz

      parameter (lmaxdim = 99)

      complex xil(0: lmaxdim), xilp(0: lmaxdim)
      complex exil(0: lmaxdim), exilp(0: lmaxdim),
     .                    exilovergam(0: lmaxdim)

      complex zetal(lmin:lmax)
      complex  zieps0, arg,
     .   al, bl, cl,
     .   gamma, zeta_eff

      integer  :: n_psi_dim
      integer :: nuper, nupar, n_psi

      integer NBESSJ

      real :: UPERP(NUPER),UPARA(NUPAR)
      real :: UPERP0, UPARA0

      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
      real :: W,K1(3),K2(3),KPER1,KPER2,ENORM,ZSPEC,ASPEC,BMAG,DensSPEC
      integer :: NSBESSJ,IFAIL
      COMPLEX WSPEC(3,3)
      complex :: factor
	real :: UminPara,UmaxPara
      real :: XI1(NUPER),JNXI1(NUPER, NBESSJ)
      real :: XI2(NUPER),JNXI2(NUPER, NBESSJ)


      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: rho_a(n_psi_dim), vc_mks


      integer :: i_uperp, i_upara, i_psi
      
      allocate( dfduper0(nuper, nupar) )
      allocate( dfdupar0(nuper, nupar) )

      nu_coll =  .01 * omgrf
      xme = 9.11e-31
      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rhol = alpha / omgc
      xkphi = xkzsav / capr
      omgrfc = omgrf * (1. + zi * xnuomg)


      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi
      xkperp = sqrt(xkalp**2 + xkbet**2)      
      
      
*     ------------------------------------
*     Optional: leave out upshift in xkprl
*     --------------------------------- --          
c      if (upshift .eq. 0)  xkprl = uzz * xkphi
      if (upshift .eq. 0)   xkprl = xkzsav / rt
      
      if (upshift .eq. -1) then
         kr = sqrt((xkxsav / xkx_cutoff)**2 
     .             + (xkysav / xky_cutoff)**2
     .             + (xkzsav / xkz_cutoff)**2)
         if (kr .gt. 1.0) xkprl = uzz * xkphi
      end if
      
      if (xkprl  .eq. 0.0) xkprl  = 1.0e-08
      if (xkperp .eq. 0.0) xkperp = 1.0e-08
                        
      sgn_kprl = sign(1.0, xkprl)
      akprl = abs(xkprl)       
      
      
      if(xkperp .gt. kperp_max)then
         write (6, *)"xkperp is gt kperp_max in sigmad"
         write (15, *)"xkperp is gt kperp_max in sigmad"
      end if
            
      
*     ---------------------------------
*     Calculate zetal(l) and gammab(l)
*     ---------------------------------      
      
      do l = lmin, lmax
         labs = abs(l)

         reson = (omgrf - l * real(omgc)) / omgrf
c         if (abs(reson) .lt. 0.02)then
c            zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
c            dzetal(l) = omgrf * xnuomg / (xkprl * alpha)
c         else
c            zetal(l) = (omgrf  - l * omgc) / (xkprl * alpha)
c            dzetal(l) = 0.0
c         end if
	 
	 zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
         dzetal(l) = omgrf * xnuomg / (xkprl * alpha)


         gammab(l) = abs(l * omgc / (2.0 * alpha * xkprl**2)
     .                                           * gradprlb / bmod)
         gamma_coll(l) = nu_coll / (akprl * alpha)


         if(xm .eq. xme)gammab(l) = 0.0
c         if(abs(gammab(l)) .gt. 1000.0) gammab(l) = 1000.0
         if(abs(gammab(l)) .lt. .01)gammab(l) = .01


      enddo

      
      
*     ------------------------------------------------
*     Calculate Brambilla's xkrpl_eff using l = 1 only
*     ------------------------------------------------
      y0 = 1.5
      y = y0
      

      if(sgn_kprl .ge. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            y = y0
            fgam = (sqrt(1. +  4. * gammab(1) * y) - 1.)
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if


      if(sgn_kprl .lt. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            descrim = 1. - 4. * gammab(1) * y0
            if (descrim .ge. 0.0) y =   y0
            if (descrim .lt. 0.0) y = - y0
            fgam = (1. - sqrt(1. -  4. * gammab(1) * y) )
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if
                    

*     -----------------------
*     Maxwellian distribution
*     -----------------------

      if(ndist .eq. 0)then

         gamma = 0.5 * xkperp**2 * rhol**2
         rgamma = real(gamma)


         if(rgamma .ge. 1.0e-08)
     .      call besiexp(gamma, lmax, exil, exilp, lmaxdim, exilovergam)

         if(rgamma .lt. 1.0e-08)
     .      call bes_expand(gamma, lmax, exil, exilp, lmaxdim,
     .                                                      exilovergam)


         sig0 = 0.0
         sig1 = 0.0
         sig2 = 0.0
         sig3 = 0.0
         sig4 = 0.0
         sig5 = 0.0



         do l = lmin, lmax
c	  do l = 0, 0
            labs = abs(l)

           if(nzfun .eq. 0) call z_approx(sgn_kprl, zetal(l), 0.0,
     .                                                     z0, z1, z2)
           if(nzfun .eq. 1) call z_approx(sgn_kprl,zetal(l),gammab(l),
     .                                                     z0, z1, z2)
           if(nzfun .eq. 2) call z_smithe(sgn_kprl,zetal(l),gammab(l),
     .                                                     z0, z1, z2)
           if(nzfun .eq. 3) call z_table(sgn_kprl,zetal(l),gammab(l),
     .                                      gamma_coll(l), z0, z1, z2)


            al = 1.0 / (xkprl * alpha) * z0
            bl = 1.0 / (xkprl * alpha) * z1
            cl = 1.0 / (xkprl * alpha) * z2


            sig0l = - zieps0 * omgp2 * rhol**2
     .                                 * (exil(labs) - exilp(labs)) * al
            sig1l = - zieps0 * omgp2 * l**2 * exilovergam(labs) * al
            sig2l = - eps0 * omgp2 * l * (exil(labs) - exilp(labs)) * al
            sig3l = - zieps0 * omgp2 * 2.0 * exil(labs) * cl
            sig4l = - zieps0 * omgp2 * rhol * l * exilovergam(labs) * bl
            sig5l = - eps0 * omgp2 * rhol
     .                                 * (exil(labs) - exilp(labs)) * bl


            sig0 = sig0 + sig0l
            sig1 = sig1 + sig1l
            sig2 = sig2 + sig2l
            sig3 = sig3 + sig3l
            sig4 = sig4 + sig4l
            sig5 = sig5 + sig5l

         end do

      end if

*     -----------------------------------
*     Funky METS calls for non Maxwellian:
*     -----------------------------------
      if (ndist .eq. 1) then
      
         if (upshift .ne. 0) xkprl = xkprl_eff
	 
      
         Emax = 0.5 * xm * vc_mks**2
         Enorm = Emax / 1.6e-19

         W = omgrf
         K1(1) = xkperp
         K1(2) = 0.0
         K1(3) = xkprl
         K2(1) = K1(1)
         K2(2) = K1(2)
         K2(3) = K1(3)
         KPER1 = xkperp
         KPER2 = xkperp


         ZSPEC = q / 1.6e-19
         ASPEC = xm / 1.67e-27
         BMAG = omgc * (xm / q)
         NSBESSJ = 2 * NBESSJ + 8
         DensSPEC = xn
         bratio = bmod0 / bmod
	 if(bratio .gt. 1.0) bratio = 1.0

	 duperp = uperp(nuper) / (nuper - 1)
	 dupara = 2.0 * upara(nupar) / (nupar - 1)

         if(i .ne. i_sav .or. j .ne. j_sav .or. k .ne. k_sav)then
	 
	    dfduper0 = 0.0
	    dfdupar0 = 0.0

!           ------------------------------------------------
!           get CQL3D distribution function on the midplane
!           ------------------------------------------------
            call cql3d_dist(nupar, nuper, n_psi,
     .                 n_psi_dim, rho_a, rho,
     .                 UminPara,UmaxPara,
     .                 df_cql_uprp, df_cql_uprl,
     .                 UPERP, UPARA, DFDUPER0, DFDUPAR0)
     

!           ------------------------------------------------
!           map CQL3D distribution function off the midplane
!           ------------------------------------------------

            if(bratio .ge. 0.0)then
	    
	       dfduper = 0.0
	       dfdupar = 0.0
	    
               do ni = 1, nuper
                  do mi = 1, nupar

                     argd = uperp(ni)**2 * (1. - bratio) 
     .                                               + upara(mi)**2
                     if (argd .le. 0.0) argd = 1.0e-06
		     
		     uperp0 = uperp(ni) * sqrt(bratio)
                     upara0 = sign(1.0, upara(mi)) * sqrt(argd)
		     
		     dfduper(ni, mi) = 0.0
                     dfdupar(ni, mi) = 0.0
		     
		     if(upara0 .ge. upara(1) .and. 
     .                                   upara0 .le. upara(nupar)) then
     
     		        ni0 = int((uperp0 - uperp(1)) / duperp) + 1
			mi0 = int((upara0 - upara(1)) / dupara) + 1
			
			dfduper0_intplt = dfduper0(ni0, mi0)
			dfdupar0_intplt = dfdupar0(ni0, mi0)
			
		        if (ni0 .lt. nuper .and. mi0 .lt. nupar) then
			
                        uperp0_grid = uperp(1) + (ni0 - 1) * duperp
			upara0_grid = upara(1) + (mi0 - 1) * dupara
						
			zeta = (uperp0 - uperp0_grid) / duperp
			eta  = (upara0 - upara0_grid) / dupara
			
                        ai = dfduper0(ni0, mi0)
                        bi = dfduper0(ni0+1 ,mi0) - dfduper0(ni0, mi0)
                        ci = dfduper0(ni0, mi0+1) - dfduper0(ni0, mi0)
                        di = dfduper0(ni0+1, mi0+1)+ dfduper0(ni0, mi0) 
     .                     - dfduper0(ni0+1, mi0) - dfduper0(ni0, mi0+1) 
			   
			dfduper0_intplt = ai + bi * zeta 
     .                                     + ci * eta + di * zeta * eta 			

                        ai = dfdupar0(ni0, mi0)
                        bi = dfdupar0(ni0+1 ,mi0) - dfdupar0(ni0, mi0)
                        ci = dfdupar0(ni0, mi0+1) - dfdupar0(ni0, mi0)
                        di = dfdupar0(ni0+1, mi0+1)+ dfdupar0(ni0, mi0) 
     .                     - dfdupar0(ni0+1, mi0) - dfdupar0(ni0, mi0+1) 
			   
			dfdupar0_intplt = ai + bi * zeta 
     .                                     + ci * eta + di * zeta * eta
     
                        end if 			

			
			if (upara0 .ne. 0.0)then
					     
			   dfdupar(ni, mi) = dfdupar0_intplt * 
     .                        upara(mi) / upara0
     
                           dfduper(ni, mi) = dfduper0_intplt * 
     .                        sqrt(bratio) + dfdupar0_intplt * 
     .                        uperp(ni) / upara0 * (1.0 - bratio)
     
c     			   dfdth = upara(mi) * dfduper(ni, mi)
c     .                           - uperp(ni) * dfdupar(ni, mi)


                        end if
			
		     end if
		     
		     
		     go to 5000
!                    ----------------------------
!                    optional analytic Maxwellian
!                    ----------------------------
		     pi = 3.141592654
		     
		     alpha = sqrt(2.0 * xkt / xm)
!                     vc_mks = 3.5 * alpha
                     u0 = vc_mks / alpha
		     
	             fnorm = u0**3 / pi**1.5 
			   
		     u2 = uperp(ni)**2 + upara(mi)**2
		     
                     f_cql = exp(-u2 * u0**2) * fnorm
	             dfduper(ni, mi) = -f_cql * 2. * uperp(ni) * u0**2
                     dfdupar(ni, mi) = -f_cql * 2. * upara(mi) * u0**2
 5000                continue		     


                  end do
               end do
	       
            end if
	    

!           --------------------------------------
!           Initialize the interpolation in k_perp:
!           --------------------------------------
            if (nkperp .ne. 0) then
 	       l_first  = .true.
               l_interp = .true.
	
               call GETNONMAXSIGMA_AORSA_NEWi(W,
     .                          ZSPEC,ASPEC,DensSPEC,BMAG,
     .                          K1,XI1,JNXI1,
     .                          K1,XI1,JNXI1,NBESSJ,
     .                          Enorm,UminPara,UmaxPara,
     .                          NUPAR,NUPER,UPERP,UPARA,
     .			        DFDUPER,DFDUPAR,
     .                          WSPEC,IFAIL,
     .                          l_first, l_interp, kperp_max, nkperp,
     .                          xkphi)
            end if



            i_sav = i
            j_sav = j
            k_sav = k

         end if






!        ------------------------------------
!        Complete integrals; no interpolation
!        ------------------------------------
         if (nkperp .eq. 0)then

            call WMATPRECALC_AORSA(ZSPEC,ASPEC,ENORM,BMAG,KPER1,UPERP,
     &                   NUPER,NBESSJ,NSBESSJ,XI1,JNXI1,IFAIL)

            call GETNONMAX_SIGMA_AORSA_NEW(W,
     .                          ZSPEC, ASPEC, DensSPEC, BMAG,
     .                          K1, XI1, JNXI1,
     .                          K1, XI1, JNXI1, NBESSJ,
     .                          Enorm, UminPara, UmaxPara,
     .                          NUPAR, NUPER, UPERP, UPARA,
     .			        DFDUPER, DFDUPAR,
     .                          WSPEC, IFAIL)

         else
!        -----------------
!        Use interpolation
!        -----------------

            l_first = .false.

	      call GETNONMAXSIGMA_AORSA_NEWi(W,
     .                       ZSPEC,ASPEC,DensSPEC,BMAG,
     .                       K1, XI1, JNXI1,
     .                       K1, XI1, JNXI1, NBESSJ,
     .                       Enorm, UminPara, UmaxPara,
     .                       NUPAR, NUPER, UPERP, UPARA,
     .			     DFDUPER, DFDUPAR,
     .                       WSPEC, IFAIL,
     .                       l_first, l_interp, kperp_max, nkperp,
     .                       xkphi)
         end if



         factor = cmplx(0.,-omgrf * eps0)

         sig1 = WSPEC(1,1) * factor
         sig2 = WSPEC(1,2) * factor
         sig3 = WSPEC(3,3) * factor
         sig4 = 0.0
         sig5 = 0.0
         sig0 = 0.0
         if (xkperp .gt. 0.01) then
            sig4 = WSPEC(3,1) / xkperp * factor
            sig5 = WSPEC(3,2) / xkperp * factor
            sig0 = (WSPEC(2,2) * factor - sig1) / xkperp**2
         endif

      end if

      sig1 = sig1 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2
      sig3 = sig3 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2
            
      if (xm .eq. xme) then
          kr = sqrt((xkxsav / xkx_cutoff)**2 
     .              + (xkysav / xky_cutoff)**2
     .              + (xkzsav / xkz_cutoff)**2)
          step = damping * kr**16 / (1. + kr**16)
	  
          sig3 = sig3 * (1.0 + step)
      end if

      go to 5500
*     --------------------------------------
*     Anti-aliasing filter (two-thirds rule):
*     -------------------------------------- 
      fnmax = 2. / 3. * real(nkx2) 
      fmmax = 2. / 3. * real(nky2) 

      nperp = sqrt((real(n) / fnmax)**2 + (real(m) / fmmax)**2)
c      if (nperp .gt. 1.0) then 
         step =  damping * nperp**16 / (1. + nperp**16) 
	                              	      	  
         sig0 = sig0 * (1.0 + step)      
         sig1 = sig1 * (1.0 + step) 
         sig2 = sig2 * (1.0 + step)   
         sig3 = sig3 * (1.0 + step)      
         sig4 = sig4 * (1.0 + step) 
         sig5 = sig5 * (1.0 + step)
	 	 
c      end if 

 5500 continue      
      
            
 	      	 	      	      
      
*     -----------------------------
*     Swanson's rotation (original):
*     -----------------------------
      sigxx = sig1 + sig0 * xkbet**2
      sigxy = sig2 - sig0 * xkbet * xkalp
      sigxz = sig4 * xkalp + sig5 * xkbet

      sigyx = - sig2 - sig0 * xkbet * xkalp
      sigyy =   sig1 + sig0 * xkalp**2
      sigyz =   sig4 * xkbet - sig5 * xkalp

      sigzx = sig4 * xkalp - sig5 * xkbet
      sigzy = sig4 * xkbet + sig5 * xkalp
      sigzz = sig3
      

      deallocate( dfduper0 )
      deallocate( dfdupar0 )


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
  102 format(2i10, 1p8e12.4)
  103 format(4i10, 1p8e12.4)
      end

c
c***************************************************************************
c




      subroutine cql3d_dist(nupar, nuper, n_psi,
     .                       n_psi_dim, rho_a, rho,
     .                       UminPara,UmaxPara,
     .                       df_cql_uprp, df_cql_uprl,
     .                       UPERP, UPARA, DFDUPER, DFDUPAR)

      implicit none
      integer   :: nupar, nuper, n_psi
      integer   :: n_psi_dim

      real   :: UminPara,UmaxPara
      real   :: UPERP(NUPER)
      real   :: UPARA(NUPAR)
      real   :: DFDUPER(NUPER,NUPAR),
     .                      DFDUPAR(NUPER,NUPAR)

      real   :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real   :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: f, u2, rho, rho_a(n_psi_dim)

      integer n, m, i_psi, i_psin

      i_psi = 1


      do i_psin = 1, n_psi - 1
         if(rho .ge. rho_a(i_psin)     .and.
     .      rho .lt. rho_a(i_psin + 1))
     .      i_psi = i_psin
      end do


      if(rho .ge. rho_a(n_psi)) i_psi = n_psi


      if(i_psi .lt. n_psi)then
         do n = 1, nuper
            do m = 1, nupar
               DFDUPER(n, m) = df_cql_uprp(n, m, i_psi)
     .            + (df_cql_uprp(n, m, i_psi + 1)
     .                                     - df_cql_uprp(n, m, i_psi))
     .                    * (rho              - rho_a(i_psi))
     .                    / (rho_a(i_psi + 1) - rho_a(i_psi))
               DFDUPAR(n, m) = df_cql_uprl(n, m, i_psi)
     .            + (df_cql_uprl(n, m, i_psi + 1)
     .                                     - df_cql_uprl(n, m, i_psi))
     .                    * (rho              - rho_a(i_psi))
     .                    / (rho_a(i_psi + 1) - rho_a(i_psi))
            end do
         end do
      end if



      if(i_psi .eq. n_psi)then
         do n = 1, nuper
            do m = 1, nupar
               DFDUPER(n, m) = df_cql_uprp(n, m, i_psi)
               DFDUPAR(n, m) = df_cql_uprl(n, m, i_psi)
            end do
         end do
      end if

  310 format(1p6e12.4)
  311 format(i10, 1p6e12.4)


      return
      end subroutine cql3d_dist

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
     .   delta0, xk0,
     .   damping, xkx_cutoff, xky_cutoff, xkz_cutoff)

*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     and uses the old left handed U matrix in paper
*     ---------------------------------------------------------

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, sgn_kprl, reson
      real dakbdkb, xnuomg, xkprl_cutoff
      real akprl, gammab, rho, alpha, eps0, omgrf, v0i
      real a, b, bmod, gradprlb
      real bx, by, bz
      real xkxsav, xkysav, xkzsav, capr
      real xkphi
      real xkalp, xkbet, xk0, rgamma

      real delta0, damping, xkx_cutoff, xky_cutoff, xkz_cutoff,
     .     kr, step

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
      complex exil(0: lmaxdim), exilp(0: lmaxdim),
     .     exilovergam(0: lmaxdim)

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

      if (xkprl  .eq. 0.0) xkprl  = 1.0e-08
      if (xkperp .eq. 0.0) xkperp = 1.0e-08

      sgn_kprl = sign(1.0, xkprl)
      akprl = abs(xkprl)


      gamma = 0.5 * xkperp**2 * rho**2
      rgamma = real(gamma)


      if(rgamma .ge. 1.0e-08)
     .     call besiexp(gamma, lmax, exil, exilp, lmaxdim, exilovergam)

      if(rgamma .lt. 1.0e-08)
     .     call bes_expand(gamma, lmax, exil, exilp, lmaxdim,
     .     exilovergam)

      sig0 = 0.0
      sig1 = 0.0
      sig2 = 0.0
      sig3 = 0.0
      sig4 = 0.0
      sig5 = 0.0


      do l = lmin, lmax
         labs = abs(l)

         reson = (omgrf - l * omgc) / omgrf
         if (abs(reson) .lt. 0.02)then
            zetal = (omgrfc - l * omgc) / (xkprl * alpha)
         else
            zetal = (omgrf - l * omgc) / (xkprl * alpha)
         end if

*        --------------------------------
*        Set gammab to zero for electrons;
*        --------------------------------
         gammab = 0.0


c         if(abs(gammab) .gt. 1000.0) gammab = 1000.0
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
         sig1l = - zieps0 * omgp2 * l**2 * exilovergam(labs) * al
         sig2l = - eps0 * omgp2 * l * (exil(labs) - exilp(labs)) * al
         sig3l = - zieps0 * omgp2 * 2.0 * exil(labs) * cl
         sig4l = - zieps0 * omgp2 * rho * l * exilovergam(labs) * bl
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
     .     + (xkysav / xky_cutoff)**2
     .     + (xkzsav / xkz_cutoff)**2)
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
     .   delta0, xk0)

*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     and uses the old left handed U matrix in paper
*     ---------------------------------------------------------

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m
      real delta0

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, sgn_kprl, reson
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
      complex exil(0: lmaxdim), exilp(0: lmaxdim),
     .     exilovergam(0: lmaxdim)

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

      if (xkprl  .eq. 0.0) xkprl  = 1.0e-08
      if (xkperp .eq. 0.0) xkperp = 1.0e-08

      sgn_kprl = sign(1.0, xkprl)
      akprl = abs(xkprl)


      gamma = 0.5 * xkperp**2 * rho**2
      rgamma = real(gamma)


      if(rgamma .ge. 1.0e-08)
     .     call besiexp(gamma, lmax, exil, exilp, lmaxdim, exilovergam)

      if(rgamma .lt. 1.0e-08)
     .     call bes_expand(gamma, lmax, exil, exilp, lmaxdim,
     .     exilovergam)

      sig0 = 0.0
      sig1 = 0.0
      sig2 = 0.0
      sig3 = 0.0
      sig4 = 0.0
      sig5 = 0.0


      do l = lmin, lmax
         labs = abs(l)

         reson = (omgrf - l * omgc) / omgrf
         if (abs(reson) .lt. 0.02)then
            zetal = (omgrfc - l * omgc) / (xkprl * alpha)
         else
            zetal = (omgrf - l * omgc) / (xkprl * alpha)
         end if


         gammab = abs(l * omgc / (2.0 * alpha * xkprl**2)
     .        * gradprlb / bmod)



c         if(abs(gammab) .gt. 1000.0) gammab = 1000.0
         if(abs(gammab) .lt. .01)gammab = .01



         if (nzfun .eq. 0) call z_approx(sgn_kprl, zetal, 0.0,
     .        z0, z1, z2)
         if (nzfun .eq. 1) call z_approx(sgn_kprl, zetal, gammab,
     .        z0, z1, z2)
         if (nzfun .eq. 2) call z_smithe(sgn_kprl, zetal, gammab,
     .        z0, z1, z2)
         if (nzfun .eq. 3) call z_table(sgn_kprl, zetal, gammab,
     .        .001, z0, z1, z2)


         al = 1.0 / (xkprl * alpha) * z0
         bl = 1.0 / (xkprl * alpha) * z1
         cl = 1.0 / (xkprl * alpha) * z2




         sig0l = - zieps0 * omgp2 * rho**2 * (exil(labs) - exilp(labs))
     .           * al
         sig1l = - zieps0 * omgp2 * l**2 * exilovergam(labs) * al
         sig2l = - eps0 * omgp2 * l * (exil(labs) - exilp(labs)) * al
         sig3l = - zieps0 * omgp2 * 2.0 * exil(labs) * cl
         sig4l = - zieps0 * omgp2 * rho * l * exilovergam(labs) * bl
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
      complex exil(0: lmaxdim), exilp(0: lmaxdim),
     .     exilovergam(0: lmaxdim)

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

      subroutine besiexp(gamma, lmax, expbes, expbesp, lmaxdim,
     .   expbesovergam)

*     ----------------------------------------------------------
*     Calculates exp(-gamma) times the modified bessel functions
*     (and derivatives) of order up to lmax
*     ----------------------------------------------------------

      implicit none

      integer lmax, nmax, ier, l, lmaxdim
      real gammod
      complex gamma, expbes(0: lmaxdim), expbesp(0: lmaxdim),
     .   expbesovergam(0: lmaxdim),
     .   xil(0: lmaxdim), xilp(0: lmaxdim), exgam

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
           expbes(l) = exgam * xil(l)
           expbesp(l) = exgam * xilp(l)
         end do
      end if

      if(gammod .gt. 700.)then
         do l = 0, lmax
            call bes_asym(gamma, l, expbes(l), expbesp(l))
         end do
      end if

      do l = 0, lmax
         expbesovergam(l) = expbes(l) / gamma
      end do

  100 format('ier = ', i5, 'besic failed')
      return
      end


c
c***************************************************************************
c


      subroutine bes_expand(gamma, lmax, expbes, expbesp, lmaxdim,
     .   expbesovergam)

*-------------------------------------------------------------------
*     Calculates exp(-gamma) times the modified bessel functions
*     (and derivatives) of order up to lmax using second order
*     expansion for small argument
*-------------------------------------------------------------------

      implicit none

      integer lmax, nmax, ier, l, lmaxdim
      real factrl, factl
      complex gamma, expbes(0:lmaxdim), expbesp(0:lmaxdim),
     .   expbesovergam(0:lmaxdim), xilovergam,
     .   xil(0:lmaxdim), xilp(0:lmaxdim), exgam

      complex b(100)

      exgam = 1.0 - gamma + gamma**2 / 2.0

      do l = 0, lmax
         factl = factrl(l)
         xil(l) = gamma**l / (2**l * factl) *
     .                                 ( 1. + gamma**2 / (4. * (l+1)))
         xilp(l) = gamma**(l-1) / (2**l * factl) *
     .                     (l + (l+2) * gamma**2 / (4. * (l+1)))

         xilovergam = gamma**(l-1) / (2**l * factl) *
     .                                 ( 1. + gamma**2 / (4. * (l+1)))

         expbes(l) = exgam * xil(l)
         expbesp(l) = exgam * xilp(l)
         expbesovergam(l) = exgam * xilovergam

c         write(6, 100)l, expbes(l), expbesp(l)
      end do

  100 format(i10, 1p8e12.4)

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


