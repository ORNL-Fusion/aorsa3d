!	module fourier_mod
!	contains

	subroutine convert3d_row_kron( row_in, row_out,
     &    nnodex,nnodey,nnodephi,
     &    nkx1,nkx2, nky1,nky2, nphi1,nphi2,
     &    inv_xx_t, inv_yy_t,inv_zz_t,
     &    xlen,ylen,philen, dx,dy,dphi )

	use kron3_mod
	implicit none
	

	integer nnodex,nnodey,nnodephi
	integer nkx1,nkx2,  nky1,nky2
	integer nphi1,nphi2

	complex inv_xx_t(1:nnodex,nkx1:nkx2)
	complex inv_yy_t(1:nnodey,nky1:nky2)
	complex inv_zz_t(1:nnodephi,nphi1:nphi2)

	real xlen,ylen,philen,   dx,dy,dphi
	
	complex row_in(3*nnodex*nnodey*nnodephi)
	complex row_out(3*nnodex*nnodey*nnodephi)

!	---------------
!	local variables
!	---------------
	complex*16, dimension(:,:,:), allocatable :: ealphak,ealpha
	complex*16, dimension(:,:,:), allocatable :: ebetak,ebeta
	complex*16, dimension(:,:,:), allocatable :: ebk,eb

	real*8 :: fact

	integer :: i,j,k,  m,n,nphi

        integer :: i_lo,i_hi,i_size
        integer :: j_lo,j_hi,j_size
        integer :: k_lo,k_hi,k_size
        integer :: m_lo,m_hi,m_size
        integer :: n_lo,n_hi,n_size
        integer :: nphi_lo,nphi_hi,nphi_size
	integer :: icnc , irnc

!       -----------------------------------------------------
!       Forward transformation  used in sft3d_kron from (i,j,k)
!	configuration space to (n,m,nphi) fourier space is
!
!	F( (n,m,nphi); (i,j,k) ) = 
!          kron( inv_zz(nphi,k), kron(inv_yy(m,j),inv_xx(n,i) )) 
!
!	ealphak(n,m,nphi) = F( (n,m,phi); (i,j,k) ) * ealpha(i,j,k)
!
!	Now we need the action of transpose(F) * vector
!
!	Let F_t( (i,j,k); (n,m,nphi) ) = transpose( F )
!
!	Let inv_xx_t(i,n) = transpose( inv_xx(n,i) )
!	Let inv_yy_t(j,m) = transpose( inv_yy(m,j) )
!	Let inv_zz_t(k,nphi) = transpose( inv_zz(nphi,k) )
!
!	Then F_t( (i,j,k); (n,m,nphi) ) = 
!	kron( inv_zz_t(k,nphi), kron(inv_yy_t(j,m), inv_xx_t(i,n)) )
!
!	ealpha(i,j,k) = F_t( (i,j,k); (n,m,nphi) ) * ealphak(n,m,nphi)
!
!       -----------------------------------------------------


	fact = (dx/xlen)*(dy/ylen)*(dphi/philen)




	nphi_lo = nphi1
	nphi_hi = nphi2

	n_lo = nkx1
	n_hi = nkx2

	m_lo = nky1
	m_hi = nky2
	
	k_lo = 1
	k_hi = nnodephi

	i_lo = 1
	i_hi = nnodex

	j_lo = 1
	j_hi = nnodey


	nphi_size = nphi_hi-nphi_lo+1
	n_size = n_hi-n_lo+1
	m_size = m_hi-m_lo+1

	k_size = k_hi-k_lo+1
	i_size = i_hi-i_lo+1
	j_size = j_hi-j_lo+1


!	---------------------------------------------------------
!	copy to ealphak(n,m,nphi), ebetak(n,m,nphi),ebk(n,m,nphi)
!	---------------------------------------------------------

	allocate( ealphak(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi) )
	allocate( ebetak(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi) )
	allocate( ebk(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi) )

	allocate( ealpha(i_lo:i_hi,j_lo:j_hi,k_lo:k_hi) )
	allocate( ebeta(i_lo:i_hi,j_lo:j_hi,k_lo:k_hi) )
	allocate( eb(i_lo:i_hi,j_lo:j_hi,k_lo:k_hi) )


	do nphi=nphi_lo,nphi_hi
	do m=m_lo,m_hi
	do n=n_lo,n_hi
	
                icnc = (nphi - nphi_lo) * 3 * m_size * n_size
     &             + (n - n_lo) * 3 * m_size
     &             + (m - m_lo) * 3 + 1

	  ealphak(n,m,nphi) = row_in(icnc)
	  ebetak(n,m,nphi)  = row_in(icnc+1)
	  ebk(n,m,nphi)     = row_in(icnc+2)

	enddo
	enddo
	enddo

        call kron3( 
     &              i_size,j_size,k_size,
     &              n_size, m_size, nphi_size,
     &              inv_zz_t,size(inv_zz_t,1),
     &              inv_yy_t,size(inv_yy_t,1),
     &              inv_xx_t,size(inv_xx_t,1),
     &              ealphak, size(ealphak,1), size(ealphak,2),
     &              ealpha, size(ealpha,1), size(ealpha,2) )

        call kron3( 
     &              i_size,j_size,k_size,
     &              n_size, m_size, nphi_size,
     &              inv_zz_t,size(inv_zz_t,1),
     &              inv_yy_t,size(inv_yy_t,1),
     &              inv_xx_t,size(inv_xx_t,1),
     &              ebetak, size(ebetak,1), size(ebetak,2),
     &              ebeta, size(ebeta,1), size(ebeta,2) )

        call kron3( 
     &              i_size,j_size,k_size,
     &              n_size, m_size, nphi_size,
     &              inv_zz_t,size(inv_zz_t,1),
     &              inv_yy_t,size(inv_yy_t,1),
     &              inv_xx_t,size(inv_xx_t,1),
     &              ebk, size(ebk,1), size(ebk,2),
     &              eb, size(eb,1), size(eb,2) )

!	---------------------------------------------------------
!	copy to ealpha(i,j,k), ebeta(i,j,k), eb(i,j,k) to row_out
!	---------------------------------------------------------

       do k = k_lo,k_hi
       do j = j_lo,j_hi
       do i = i_lo,i_hi

                irnc = (k-k_lo) * j_size * i_size * 3
     &             + (i-i_lo) * j_size * 3
     &             + (j-j_lo) * 3
     &             + 1

		row_out(irnc)   = fact * ealpha(i,j,k)
		row_out(irnc+1) = fact * ebeta(i,j,k)
		row_out(irnc+2) = fact * eb(i,j,k)
             end do
          end do
       end do

	deallocate( ealphak )
	deallocate( ebetak )
	deallocate( ebk )

	deallocate( ealpha )
	deallocate( ebeta )
	deallocate( eb )

	return
	end subroutine convert3d_row_kron




      subroutine sft3d_kron(f, fk,
     .   nxdim, nydim, nphidim,
     .   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz,
     .   xlen, ylen, philen, dx, dy, dphi)

      use kron3_mod
      implicit none

	integer, parameter :: idebug = 0

      integer nxdim, nydim, nphidim, n, m, nphi, i, j, k,
     .   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2

      real x(nxdim), y(nydim), phi(nphidim)
      real xlen, ylen, philen, dx, dy, dphi, fact

      real xkx(nkdim1 : nkdim2),
     .     xky(mkdim1 : mkdim2),
     .     xkz(nphidim1 : nphidim2)
      complex f(nxdim, nydim, nphidim),
     .   fk(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 : nphidim2)
      complex cexpkxkykz

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim),
     .        zz(nphidim1 : nphidim2, 1 : nphidim)

	integer :: i_lo,i_hi,i_size
	integer :: j_lo,j_hi,j_size
	integer :: k_lo,k_hi,k_size
	integer :: m_lo,m_hi,m_size
	integer :: n_lo,n_hi,n_size
	integer :: nphi_lo,nphi_hi,nphi_size
	complex*16, dimension(:,:), allocatable :: inv_xx
	complex*16, dimension(:,:), allocatable :: inv_yy
	complex*16, dimension(:,:), allocatable :: inv_zz
	complex*16, dimension(:,:,:), allocatable :: f_in
	complex*16, dimension(:,:,:), allocatable :: fk_kron
	real*8 :: err
	
	nphi_lo = nphi1
	nphi_hi = nphi2

	n_lo = nkx1
	n_hi = nkx2

	m_lo = nky1
	m_hi = nky2
	
	k_lo = 1
	k_hi = nnodephi

	i_lo = 1
	i_hi = nnodex

	j_lo = 1
	j_hi = nnodey


	nphi_size = nphi_hi-nphi_lo+1
	n_size = n_hi-n_lo+1
	m_size = m_hi-m_lo+1

	k_size = k_hi-k_lo+1
	i_size = i_hi-i_lo+1
	j_size = j_hi-j_lo+1

      fact = dx/xlen * dy/ylen * dphi/philen

      	allocate( fk_kron(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi) )
	allocate( f_in(i_lo:i_hi,j_lo:j_hi,k_lo:k_hi) )

	allocate( inv_xx(n_lo:n_hi,i_lo:i_hi) )
	allocate( inv_yy(m_lo:m_hi,j_lo:j_hi) )
	allocate( inv_zz(nphi_lo:nphi_hi,k_lo:k_hi) )

	do i=i_lo,i_hi
	do n=n_lo,n_hi
	  inv_xx(n,i) = 1.0d0/xx(n,i)
	enddo
	enddo

	do j=j_lo,j_hi
	do m=m_lo,m_hi
	  inv_yy(m,j) = 1.0d0/yy(m,j)
	enddo
	enddo

	do k=k_lo,k_hi
	do nphi=nphi_lo,nphi_hi
	  inv_zz(nphi,k) = 1.0d0/zz(nphi,k)
	enddo
	enddo


!	-----------------------------------------------------
!	fk(n,m,nphi) = 
!	kron( inv_zz(nphi,k), kron(inv_yy(m,j),inv_xx(n,i) )) * f(i,j,k)
!	-----------------------------------------------------
	f_in(i_lo:i_hi,j_lo:j_hi,k_lo:k_hi) = 
     &            f(i_lo:i_hi,j_lo:j_hi,k_lo:k_hi)

  	call kron3( n_size, m_size, nphi_size,  
     &              i_size,j_size,k_size,
     &		    inv_zz,size(inv_zz,1), 
     &		    inv_yy,size(inv_yy,1),
     &		    inv_xx,size(inv_xx,1),
     &		    f_in, size(f_in,1), size(f_in,2), 
     &              fk_kron, size(fk_kron,1), size(fk_kron,2) )

	fk(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi) = 
     &        fact * fk_kron(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi)


!debug-begin
	if (idebug.ge.2) then
!	------------
!	double check
!	------------
	fk_kron(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi) = 
     &          fk(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi)


      do nphi = nphi1, nphi2
         do n = nkx1, nkx2
            do m = nky1, nky2

               fk(n, m, nphi) = 0.0

               do k = 1, nnodephi
                  do i = 1, nnodex
                     do j = 1, nnodey

                        cexpkxkykz = xx(n, i) * yy(m, j) * zz(nphi, k)


                        fk(n, m, nphi) = fk(n, m, nphi)
     .                                 + f(i,j,k) / cexpkxkykz * fact

                     end do
                  end do
               end do

            end do
         end do
      end do

	err = maxval( abs( fk(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi) -
     &                 fk_kron(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi) ))
        write(*,*) 'sft3d_kron: err fk = ', err
        write(15,*) 'sft3d_kron: err fk = ', err

	endif
!debug-end

	deallocate( f_in )
	deallocate( fk_kron )
	deallocate( inv_xx )
	deallocate( inv_yy )
	deallocate( inv_zz )


      return
      end subroutine sft3d_kron



      subroutine fftinv3d_kron(x, y, phi, xkx, xky, xkz, f, a,
     .   nxdim, nydim, nphidim,
     1   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2,
     .   nnodex, nnodey, nnodephi,
     1   nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz)

	use kron3_mod

      implicit none

	integer, parameter :: idebug = 0

      integer nxdim, nydim, nphidim, n, m, nphi, i, j, k,
     1   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2,
     .   nnodex, nnodey, nnodephi,
     1   nkx1, nkx2, nky1, nky2, nphi1, nphi2

      real x(nxdim), y(nydim), phi(nphidim)
      real xkx(nkdim1 : nkdim2),
     .     xky(mkdim1 : mkdim2),
     .     xkz(nphidim1 : nphidim2)
      complex f(nxdim, nydim, nphidim),
     1   a(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 : nphidim2)
      complex cexpkxkykz

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim),
     .        zz(nphidim1 : nphidim2, 1 : nphidim)

	integer :: i_lo,i_hi, j_lo,j_hi, k_lo,k_hi
	integer :: n_lo,n_hi, m_lo,m_hi, nphi_lo,nphi_hi
	integer :: i_size,j_size,k_size
	integer :: m_size,n_size,nphi_size

	complex*16, dimension(:,:), allocatable :: xx_t,yy_t,zz_t
	complex*16, dimension(:,:,:), allocatable :: Ain,Fout
	real*8 :: err

	i_lo = 1
	i_hi = nnodex

	j_lo = 1
	j_hi = nnodey

	k_lo = 1
	k_hi = nnodephi

	nphi_lo = nphi1
	nphi_hi = nphi2
	
	n_lo = nkx1
	n_hi = nkx2

	m_lo = nky1
	m_hi = nky2

	i_size = i_hi-i_lo+1
	j_size = j_hi-j_lo+1
	k_size = k_hi-k_lo+1
	
	m_size = m_hi-m_lo+1
	n_size = n_hi-n_lo+1
	nphi_size = nphi_hi - nphi_lo + 1 

!	---------------------------
!	generate transpose matrices
!	---------------------------
	allocate(xx_t(i_lo:i_hi,n_lo:n_hi))
	allocate(yy_t(j_lo:j_hi,m_lo:m_hi))
	allocate(zz_t(k_lo:k_hi,nphi_lo:nphi_hi))

	do n=n_lo,n_hi
	do i=i_lo,i_hi
	   xx_t(i,n) = xx(n,i)
	enddo
	enddo

	do m=m_lo,m_hi
	do j=j_lo,j_hi
	  yy_t(j,m) = yy(m,j)
	enddo
	enddo

	do nphi=nphi_lo,nphi_hi
	do k=k_lo,k_hi
	   zz_t(k,nphi) = zz(nphi,k)
	enddo
	enddo

!	------------------------------------------------
!	The nested loop is performing  matrix-vector multiply
!	where the matrix is conceptually the kronecker product
!	of 3 matrices, zz_t, yy_t, xx_t
!
!	f(i,j,k) = 
!	  kron(zz_t(k,nphi), kron( yy_t(j,m), xx_t(i,n))*a(n,m,nphi)
!	------------------------------------------------
	allocate( Fout(i_lo:i_hi,j_lo:j_hi,k_lo:k_hi) )
	allocate( Ain(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi) )

	Ain(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi) = 
     &      a(n_lo:n_hi,m_lo:m_hi,nphi_lo:nphi_hi)

	call kron3(i_size,j_size,k_size, 
     &             n_size,m_size,nphi_size,
     &		zz_t, size(zz_t,1), 
     &          yy_t, size(yy_t,1), 
     &          xx_t, size(xx_t,1),
     &          Ain, size(Ain,1), size(Ain,2), 
     &          Fout, size(Fout,1), size(Fout,2) )

	f(i_lo:i_hi,j_lo:j_hi,k_lo:k_hi) = 
     &       Fout(i_lo:i_hi,j_lo:j_hi,k_lo:k_hi)



!debug-begin
      if (idebug.ge.2) then
      do k = 1, nnodephi
         do i = 1, nnodex
            do j = 1, nnodey
               f(i,j,k) = 0.0

               do nphi = nphi1, nphi2
c                 if(nphi .ne. 0) then
                     do n = nkx1, nkx2
                        do m = nky1, nky2

c                           cexpkxkykz = exp(zi * (xkx(n) * x(i)
c     .                                          + xky(m) * y(j)
c     .                                          + xkz(nphi) * phi(k)))

                           cexpkxkykz = xx(n, i) * yy(m, j) * zz(nphi,k)

                           f(i,j,k) = f(i,j,k) + a(n,m,nphi) *cexpkxkykz
                        end do
                     end do
c                 end if
               end do

            end do
         end do
      end do

	err = maxval( abs(f(i_lo:i_hi,j_lo:j_hi,k_lo:k_hi) -
     &                    Fout(i_lo:i_hi,j_lo:j_hi,k_lo:k_hi)))
	write(*,*) 'fftinv3d, err Fout ', err
	write(15,*) 'fftinv3d, err Fout ', err
      endif
!debug-end

	deallocate( xx_t )
	deallocate( yy_t )
	deallocate( zz_t )
	deallocate( Ain )
	deallocate( Fout )

      return
      end subroutine fftinv3d_kron
c*******************************************************************
c
      subroutine sft3d(f, fk,
     .   nxdim, nydim, nphidim,
     .   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz,
     .   xlen, ylen, philen, dx, dy, dphi)

      implicit none

      integer nxdim, nydim, nphidim, n, m, nphi, i, j, k,
     .   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2

      real x(nxdim), y(nydim), phi(nphidim)
      real xlen, ylen, philen, dx, dy, dphi, fact

      real xkx(nkdim1 : nkdim2),
     .     xky(mkdim1 : mkdim2),
     .     xkz(nphidim1 : nphidim2)
      complex f(nxdim, nydim, nphidim),
     .   fk(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 : nphidim2)
      complex cexpkxkykz

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim),
     .        zz(nphidim1 : nphidim2, 1 : nphidim)

      fact = dx/xlen * dy/ylen * dphi/philen

      do nphi = nphi1, nphi2
         do n = nkx1, nkx2
            do m = nky1, nky2

               fk(n, m, nphi) = 0.0

               do k = 1, nnodephi
                  do i = 1, nnodex
                     do j = 1, nnodey

                        cexpkxkykz = xx(n, i) * yy(m, j) * zz(nphi, k)


                        fk(n, m, nphi) = fk(n, m, nphi)
     .                                 + f(i,j,k) / cexpkxkykz * fact

                     end do
                  end do
               end do

            end do
         end do
      end do




      return
      end subroutine sft3d


c
c*******************************************************************
c
      subroutine convert3d_row(row_in, row_out,
     .   xlen, ylen, philen,
     .   nnodex, nnodey, nnodephi,
     .   nxdim, nydim, nphidim,
     .   nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .   xx, yy, zz, dx, dy, dphi,
     .   ndfmax, nmodesx, nmodesy, nmodesphi)

*     ---------------------------------------------------
*     performs a Fourier transform on one row of the 3-D
*     collocation matrix, p_amat
*        row_in(ncol)  = input  vector of dimension ndfmax
*        row_out(ncol) = output vector of dimension ndfmax
*     ----------------------------------------------------

      implicit none

      integer nnodex, nnodey, nnodephi, nxdim, nydim, nphidim,
     .   nkdim1, nkdim2,
     .   mkdim1, mkdim2, lkdim1, lkdim2, i, j, k, n, m, nphi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .   ndfmax, ncol, nmodesx, nmodesy, nmodesphi

      real xlen, ylen, philen, dx, dy, dphi, fact

      complex ealphak(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 :lkdim2),
     .        ebetak(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        ebk(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2)

      complex ealpha(nxdim, nydim, nphidim),
     .        ebeta(nxdim, nydim, nphidim),
     .        eb(nxdim, nydim, nphidim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim),
     .        zz(lkdim1 : lkdim2, 1 : nphidim)

      complex row_in(ndfmax), row_out(ndfmax)


      do nphi = nphi1, nphi2
         do n = nkx1, nkx2
            do m = nky1, nky2

               ncol = (nphi - nphi1) * 3 * nmodesx * nmodesy
     .             + (n - nkx1) * 3 * nmodesy
     .             + (m - nky1) * 3 + 1

               ealphak(n, m, nphi) = row_in(ncol)
               ebetak(n, m, nphi)  = row_in(ncol + 1)
               ebk(n, m, nphi)     = row_in(ncol + 2)

            end do
         end do
      end do


      call fftmat3d(ealphak, ealpha,
     .   nxdim, nydim, nphidim,
     .   nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz)

      call fftmat3d(ebetak, ebeta,
     .   nxdim, nydim, nphidim,
     .   nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz)

      call fftmat3d(ebk, eb,
     .   nxdim, nydim, nphidim,
     .   nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz)

      fact = dx / xlen * dy / ylen * dphi / philen

      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nnodephi

               ncol = (k - 1) * nnodex * nnodey * 3
     .              + (i - 1) * nnodey * 3
     .              + (j - 1) * 3
     .              + 1

               row_out(ncol)     = ealpha(i, j, k) * fact
               row_out(ncol + 1) = ebeta(i, j, k)  * fact
               row_out(ncol + 2) = eb(i, j, k)     * fact

            end do
         end do
      end do



      return
      end subroutine convert3d_row


c
c*******************************************************************
c
      subroutine fftmat3d(a, f,
     .   nxdim, nydim, nphidim,
     1   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2,
     .   nnodex, nnodey, nnodephi,
     1   nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz)

      implicit none

      integer nxdim, nydim, nphidim, n, m, nphi, i, j, k,
     1   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2,
     .   nnodex, nnodey, nnodephi,
     1   nkx1, nkx2, nky1, nky2, nphi1, nphi2

      real x(nxdim), y(nydim), phi(nphidim)
      real xkx(nkdim1 : nkdim2),
     .     xky(mkdim1 : mkdim2),
     .     xkz(nphidim1 : nphidim2)
      complex f(nxdim, nydim, nphidim),
     1   a(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 : nphidim2)
      complex cexpkxkykz

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim),
     .        zz(nphidim1 : nphidim2, 1 : nphidim)


      do k = 1, nnodephi
         do i = 1, nnodex
            do j = 1, nnodey

               f(i,j,k) = 0.0

               do nphi = nphi1, nphi2
                  do n = nkx1, nkx2
                     do m = nky1, nky2

                        cexpkxkykz = xx(n, i) * yy(m, j) * zz(nphi, k)
                        f(i,j,k) = f(i,j,k) + a(n,m,nphi) / cexpkxkykz

                     end do
                  end do
               end do

            end do
         end do
      end do


      return
      end subroutine fftmat3d
c
c*******************************************************************
c
      subroutine fftinv3d(x, y, phi, xkx, xky, xkz, f, a,
     .   nxdim, nydim, nphidim,
     1   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2,
     .   nnodex, nnodey, nnodephi,
     1   nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz)

      implicit none

      integer nxdim, nydim, nphidim, n, m, nphi, i, j, k,
     1   nkdim1, nkdim2, mkdim1, mkdim2, nphidim1, nphidim2,
     .   nnodex, nnodey, nnodephi,
     1   nkx1, nkx2, nky1, nky2, nphi1, nphi2

      real x(nxdim), y(nydim), phi(nphidim)
      real xkx(nkdim1 : nkdim2),
     .     xky(mkdim1 : mkdim2),
     .     xkz(nphidim1 : nphidim2)
      complex f(nxdim, nydim, nphidim),
     1   a(nkdim1 : nkdim2, mkdim1 : mkdim2, nphidim1 : nphidim2)
      complex cexpkxkykz

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim),
     .        zz(nphidim1 : nphidim2, 1 : nphidim)


      do k = 1, nnodephi
         do i = 1, nnodex
            do j = 1, nnodey
               f(i,j,k) = 0.0

               do nphi = nphi1, nphi2
c                 if(nphi .ne. 0) then
                     do n = nkx1, nkx2
                        do m = nky1, nky2

c                           cexpkxkykz = exp(zi * (xkx(n) * x(i)
c     .                                          + xky(m) * y(j)
c     .                                          + xkz(nphi) * phi(k)))

                           cexpkxkykz = xx(n, i) * yy(m, j) * zz(nphi,k)

                           f(i,j,k) = f(i,j,k) + a(n,m,nphi) *cexpkxkykz
                        end do
                     end do
c                 end if
               end do

            end do
         end do
      end do

      return
      end subroutine fftinv3d
c
c*******************************************************************
c



      subroutine sftinv2d(x, y, xkx, xky, f, a, nxdim, nydim,
     1   nkdim1, nkdim2, mkdim1, mkdim2, nnodex, nnodey,
     1   nkx1, nkx2, nky1, nky2)

      implicit none

      integer nxdim, nydim, n, m, i, j,
     1   nkdim1, nkdim2, mkdim1, mkdim2, nnodex, nnodey,
     1   nkx1, nkx2, nky1, nky2

      real x(nxdim), y(nydim)
      real xkx(nkdim1 : nkdim2), xky(mkdim1 : mkdim2)
      complex zi, f(nxdim, nydim),
     1   a(nkdim1 : nkdim2, mkdim1 : mkdim2)
      complex cexpkxky


      zi = cmplx(0.0, 1.0)

      do i = 1, nnodex
         do j = 1, nnodey
            f(i, j) = 0.0
            do n = nkx1, nkx2
               do m = nky1, nky2
                  cexpkxky = exp(zi * (xkx(n) * x(i) + xky(m) * y(j)))
                  f(i,j) = f(i,j) + a(n,m) * cexpkxky
               end do
            end do
         end do
      end do

      return
      end subroutine sftinv2d
c
c*******************************************************************
c

      subroutine sftinv2dp(x, y, xkx, xky, f, a, nxdim, nydim,
     1   nkdim1, nkdim2, mkdim1, mkdim2, nnodex, nnodey,
     1   nkx1, nkx2, nky1, nky2,
     1   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_)

      implicit none

      logical ismine

      integer nxdim, nydim, n, m, i, j,
     1   nkdim1, nkdim2, mkdim1, mkdim2, nnodex, nnodey,
     1   nkx1, nkx2, nky1, nky2, irnc, icnc

      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,
     .   icontxt
      integer dlen_, desc_amat(dlen_)

      real x(nxdim), y(nydim)
      real xkx(nkdim1 : nkdim2), xky(mkdim1 : mkdim2)
      complex zi, f(nxdim, nydim),
     1   a(nkdim1 : nkdim2, mkdim1 : mkdim2)
      complex cexpkxky


      zi = cmplx(0.0, 1.0)

c--   Loop over mode numbers and mesh:
      do i = 1, nnodex
         do j = 1, nnodey

            f(i, j) = 0.0

            irnc = (j-1) * 3 + (i-1) * nnodey * 3 + 1

            do n = nkx1, nkx2
               do m = nky1, nky2

                  icnc = (n - nkx1) * 3 * nnodey + (m - nky1) * 3 + 1

c                  call infog2l(irnc, icnc, desc_amat, nprow, npcol,
c     .                myrow, mycol, lrindx, lcindx, rsrc, csrc)

c                  ismine = (rsrc .eq. myrow) .and. (csrc .eq. mycol)


c                  if (ismine) then
                     cexpkxky = exp(zi * (xkx(n) * x(i) + xky(m) *y(j)))
                     f(i,j) = f(i,j) + a(n,m) * cexpkxky
c                  end if

               end do
            end do
         end do
      end do

c      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, f, nxdim,-1, -1)

      return 
      end subroutine sftinv2dp
c
c*******************************************************************
c


      subroutine sft(x, fx, xlen, nx, nxdim, nkdim1, nkdim2,
     1   nk1, nk2, an)

      implicit none

      integer nx, nxdim, nkdim1, nkdim2, n, nk, nk1, nk2
      real x(nxdim), xlen, xkn, pi
      complex fx(nxdim), an(nkdim1:nkdim2), zi, xint

      pi = 3.14159
      zi = cmplx(0.0, 1.0)

      do nk = nk1, nk2
         xkn = 2.0 * pi/xlen * nk
         an(nk) = 0.0
         do n = 1, nx-1
            xint = (fx(n)   * exp(-zi * xkn * x(n))
     1            + fx(n+1) * exp(-zi * xkn * x(n+1)) )/ 2.0
            an(nk) = an(nk) + xint * (x(n+1) - x(n)) / xlen
         end do
      end do

      return
      end subroutine sft
c
c*******************************************************************
c
      subroutine sftinv(x, fx, xlen, nx, nxdim, nkdim1, nkdim2,
     1    nk1, nk2, an)

      implicit none

      integer nx, nxdim, nkdim1, nkdim2, n, nk, nk1, nk2
      real x(nxdim), xlen, xkn, pi
      complex zi, fx(nxdim), an(nkdim1 : nkdim2)

      pi = 3.141592654
      zi = cmplx(0.0, 1.0)

      do n = 1, nx
         fx(n) = 0.0
         do nk = nk1, nk2
            xkn = 2.0 * pi * nk / xlen
            fx(n) = fx(n) + an(nk) * exp(zi * xkn * x(n))
         end do
      end do

      return
      end subroutine sftinv
c
c*******************************************************************
c
      subroutine fft(nx, work, nxdim, nkdim1, nkdim2, fx, an)

      implicit none

      integer nhalf, nx, nxdim, nkdim1, nkdim2, n, nk
      complex fx(nxdim), work(nxdim), an(nkdim1 : nkdim2)
      real x(nxdim)

      nhalf = nx / 2
      do n = 1, nx
         work(n) = fx(n)
      end do

      call four1(work, nx, -1)

      do n = 1, nx
         nk = n - 1
         if(nk .gt. nhalf) nk = nk - nx
         an(nk) = work(n) / nx
      end do
      return
      end subroutine fft
c
c***************************************************************************
c
      subroutine fftinv(nx, work, nxdim, nkdim1, nkdim2, fx, an)

      implicit none

      integer nhalf, nx, nxdim, nkdim1, nkdim2, n, nk
      complex fx(nxdim), work(nxdim), an(nkdim1 : nkdim2)
      real x(nxdim)

      nhalf = nx / 2

      do n = 1, nx
         nk = n - 1
         if(nk .gt. nhalf) nk = nk - nx
         work(n) = an(nk)
      end do

      call four1(work, nx, 1)

      do n = 1, nx
         fx(n) = work(n)
      end do

      return
      end subroutine fftinv
c
c***************************************************************************
c
      subroutine fftinvr(nx, work, nxdim, nkdim1, nkdim2, fx, an)

      implicit none

      integer nhalf, nx, nxdim, nkdim1, nkdim2, n, nk
      real fx(nxdim), an(nkdim1 : nkdim2)
      real x(nxdim)
      complex work(nxdim)

      nhalf = nx / 2

      do n = 1, nx
         nk = n - 1
         if(nk .gt. nhalf) nk = nk - nx
         work(n) = cmplx(an(nk), 0.0)
c         write(6, 100) n, work(n)
      end do

      call four1(work, nx, 1)

      do n = 1, nx
         fx(n) = real(work(n))
      end do
  100 format (i10, 1p8e12.4)
      return
      end subroutine fftinvr
c
c*******************************************************************
c
      subroutine four1(cdata, nn, isign)

      implicit none

      integer isign, nn
      complex cdata(nn)
      real data(2*nn)
c         Replaces data(1:2*nn) by its discrete Fourier transform, if isign
c         is input as 1; or replaces
c         data(1:2*nn) by nn times its inverse discrete Fourier transform,
c         if isign is input as -1.
c         data is a complex array of length nn or, equivatlently a real
c         array of length 2*nn.  nn MUST be an integer power of 2
c         (this is not checked for!).
      integer i, istep, j, m, mmax, n
      real tempi, tempr
c          Double precision for the trigonometic recurrences
      double precision theta, wi, wpi, wpr, wr, wtemp

      do j=1,nn
        data( 2*j-1 ) = real(cdata(j))
        data( 2*j   ) = aimag(cdata(j))
      enddo


      n = 2 * nn
      j = 1
c          This is bit reversal section of the routine
      do i = 1, n, 2
         if (j .gt. i) then
            tempr = data(j)
            tempi = data(j+1)
            data(j) = data(i)
            data(j+1) = data(i+1)
            data(i) = tempr
            data(i+1) = tempi
         endif
         m = n/2

    1    if (m.ge.2 .and. j.gt.m) then
            j = j - m
            m = m/2
            go to 1
         endif
         j = j + m
      end do
      mmax = 2
    2 if (n .gt. mmax) then
         istep = 2*mmax
         theta = 6.28318530717959D0/(isign*mmax)
         wpr=-2.d0*dsin(0.5d0*theta)**2

         wpi=dsin(theta)

         wr=1.d0

         wi=0.d0

         do m = 1, mmax, 2
            do i = m, n, istep
               j = i + mmax
               tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
               tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)

               data(j) = data(i) - tempr
               data(j+1) = data(i+1) - tempi
               data(i) = data(i) + tempr
               data(i+1) = data(i+1) + tempi
            end do
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
         end do
         mmax = istep
         go to 2
      endif

      do j=1,nn
        cdata(j) = cmplx( data(2*j-1), data(2*j) )
      enddo

      return
      end subroutine four1

!      end module fourier_mod

