	module kron3_mod
	implicit none

	private

	interface 
	subroutine zgemm(transA, transB, m,n,k,
     &		alpha, A, lda, B, ldb, beta, C, ldc )
	character transA, transB
	integer m,n,k,lda,ldb,ldc
	complex*16 alpha,beta, A(*),B(*),C(*)
	end subroutine zgemm
	end interface

	integer, parameter :: idebug = 0
	logical, parameter ::   use_gemm = .true.

	public :: kron2, kron3

	contains

	subroutine kron3( ny1, ny2, ny3, nx1,nx2,nx3,
     &		A, lda, B, ldb, C, ldc, 
     &          X, ldx1,ldx2, Y, ldy1,ldy2 )
	implicit none
	integer, intent(in) :: ny1,ny2,ny3, nx1,nx2,nx3
	integer, intent(in) :: lda,ldb,ldc, ldx1,ldx2, ldy1,ldy2

	complex*16, dimension(lda,*), intent(in) :: A
	complex*16, dimension(ldb,*), intent(in) :: B
	complex*16, dimension(ldc,*), intent(in) :: C
	complex*16, dimension(ldx1,ldx2,*), intent(in) :: X
	complex*16, dimension(ldy1,ldy2,*), intent(inout) :: Y

!	-----------------------------------------------------	
!	Y(ny1,ny2,ny3) = kron( A, kron(B,C)) * X(nx1,nx2,nx3)
!
!	A is (ny3,nx3), B is (ny2,nx2), C is (ny1,nx1)
!
!	kron(A, kron(B,C))*X = (kron(B,C)*X)*transpose(A)
!
!	Let be BCX(ny1,ny2,nx3) = 
!		kron(B(ny2,nx2), C(ny1,nx1))* X(nx1,nx2,nx3)
!
!	-----------------------------------------------------	
!	---------------
!	local variables
!	---------------
	complex*16, dimension(:,:), allocatable ::  Xin
	complex*16, dimension(:,:), allocatable ::  Yout, Youtb
	complex*16, dimension(:,:), allocatable ::  A_t
	complex*16, dimension(:,:,:), allocatable :: BCX, BCXb
	real*8 :: err
	integer :: ix1,ix2,ix3, iy1,iy2,iy3
	integer :: ldxin, ldyout 

	integer :: mm,nn,kk, ld1,ld2,ld3
	complex*16 :: alpha,beta
!


!	-------------------------------------------------------
!	Note computing kron(A,B,C)*X as
!
!	(kron(B,C)*X)*transpose(A)
!
!	is chosen for convenience in implementation and
!	may not be most efficient in terms of temporary storage
!	or amount of work.
!
!	-------------------------------------------------------

	allocate( BCX(ny1,ny2,nx3) )
	allocate( Xin(nx1,nx2) )
	allocate( Yout(ny1,ny2) )

	do ix3=1,nx3
	   do ix2=1,nx2
	   do ix1=1,nx1
		Xin(ix1,ix2) = X(ix1,ix2,ix3)
	   enddo
	   enddo
	
	   ldxin = size(Xin,1)
	   ldyout = size(Yout,1)

	   call kron2(  ny1,ny2, nx1,nx2, 
     &       B, ldb, C, ldc, Xin, ldxin, Yout,ldyout )

	   do iy2=1,ny2
	   do iy1=1,ny1
		BCX(iy1,iy2,ix3) = Yout(iy1,iy2)
	   enddo
	   enddo

	enddo
	deallocate( Xin )
	deallocate( Yout )

!debug-begin
	if (idebug.ge.2) then
	allocate( BCXb(ny1,ny2,nx3) )
	BCXb(:,:,:) = 0.0

	do ix3=1,nx3
	do iy2=1,ny2
	do iy1=1,ny1

	do ix2=1,nx2
	do ix1=1,nx1
	
	  BCXb(iy1,iy2,ix3) = BCXb(iy1,iy2,ix3) + 
     &      B(iy2,ix2) * C(iy1,ix1) * X(ix1,ix2,ix3)
	enddo
	enddo

	enddo
	enddo
	enddo

	err = maxval( abs(BCXb - BCX) )
	write(*,*) 'err BCXb ', err

        deallocate( BCXb )
	endif
!debug-end

!	----------------------------------------------------
!	Y(1:ny1,1:ny2,1:ny3) = 
!		BCX(1:ny1,1:ny2,1:nx3)*tranpose(A(ny3,nx3))
!	----------------------------------------------------
	allocate( A_t(1:nx3,1:ny3) )
	do ix3=1,nx3
	do iy3=1,ny3
	   A_t(ix3,iy3) = A(iy3,ix3)
	enddo
	enddo


	allocate( Xin(1:ny1,1:nx3) )
	allocate( Yout(1:ny1,1:ny3) )
	
	do iy2=1,ny2
!	   ---------------------------------------
!	   Xin(1:ny1,1:nx3) = BCX(1:ny1,iy2,1:nx3)
!	   ---------------------------------------
	   do ix3=1,nx3
	   do iy1=1,ny1
	     Xin(iy1,ix3) = BCX(iy1,iy2,ix3)
	   enddo
	   enddo

	   if (use_gemm) then
	     mm = ny1
	     nn = ny3
	     kk = nx3
	     ld1 = size(Xin,1)
	     ld2 = size(A_t,1)
	     ld3 = size(Yout,1)
	     alpha = 1.0d0
	     beta = 0.0d0

	     call zgemm( 'N','N', mm,nn,kk,
     &           alpha, Xin, ld1, A_t, ld2, beta, Yout, ld3 )

!debug-begin
	if (idebug.ge.2) then
	     allocate( Youtb(1:ny1,1:ny3) )

	     Youtb(1:ny1,1:ny3) = 
     &		matmul( Xin(1:ny1,1:nx3), A_t(1:nx3,1:ny3) )
	     err = maxval(abs(Youtb-Yout))
	     write(*,*) 'err Youtb ', err

	     deallocate( Youtb )
	endif	
!debug-end
	   else
	     Yout(1:ny1,1:ny3) = 
     &		matmul( Xin(1:ny1,1:nx3), A_t(1:nx3,1:ny3) )
	   endif
!	   ----------------------------------------
!	   Y(1:ny1, iy2, 1:ny3) = Yout(1:ny1,1:ny3)
!	   ----------------------------------------
	   do iy3=1,ny3
	   do iy1=1,ny1
	      Y(iy1,iy2,iy3) = Yout(iy1,iy3)
	   enddo
	   enddo

	enddo

	deallocate( Xin )
	deallocate( Yout )
	deallocate( A_t )

	deallocate( BCX )


	
	   

	return
	end subroutine kron3
	





	subroutine kron2( ny1, ny2, nx1,nx2,  
     &           A, lda, B, ldb, X, ldx, Y, ldy )
	implicit none
	integer, intent(in) :: ny1,ny2, nx1,nx2
	integer, intent(in) :: lda,ldb,ldx,ldy

	complex*16, dimension(lda,*), intent(in) :: A
	complex*16, dimension(ldb,*), intent(in) :: B
	complex*16, dimension(ldx,*), intent(in) :: X
	complex*16, dimension(ldy,*), intent(inout) :: Y
!	
!	compute Y = kron(A,B) * X
!
!	Y = B * X * transpose(A)
!
!	Y(ny1,ny2) = B(ny2,nx1) * X(nx1,nx2) * 
!                                transpose( A(ny2,nx2) )
!

!	---------------
!	local variables
!	---------------
	complex*16, dimension(:,:), allocatable :: BX, XA, A_t
	integer :: size_BX, size_XA
	integer :: flops_BX, flops_XA
	logical :: use_BX
	integer :: mm,nn,kk, ld1,ld2,ld3
	integer :: i,j
	complex*16 :: alpha,beta



	size_BX = ny1 * nx2
	size_XA = nx1 * ny2

	flops_BX = nx1*ny1*nx2 + nx2*ny1*ny2
	flops_XA = nx2*nx1*ny2 + nx1*ny1*nx2

	use_BX = (size_BX.lt.size_XA) .or.
     &      ((size_BX.eq.size_XA).and.(flops_BX.lt.flops_XA))

	if (use_BX) then
!	----------------------------------
!	compute as Y = (B*X) * tranpose(A)
!	BX(1:ny1,1:nx2) = B(1:ny1,1:nx1) * X(1:nx1, 1:nx2)
!	----------------------------------
	allocate( BX(ny1,nx2) )

	if (use_gemm) then
	   mm = ny1
	   nn = nx2
	   kk = nx1
	   ld1 = size(B,1)
	   ld2 = size(X,1)
	   ld3 = size(BX,1)
	   alpha = 1.0d0
	   beta = 0.0d0
	   call zgemm( 'N','N', mm,nn,kk,
     &		alpha, B,ld1, X, ld2, beta, BX, ld3 )

!	  ------------------------------------------
!	  Y(1:ny1,1:ny2) = matmul( BX(1:ny1, 1:nx2), 
!               transpose( A(1:ny2,1:nx2) ) )
!	  ------------------------------------------

	   mm = ny1
	   nn = ny2
	   kk = nx2
	   ld1 = size(BX,1)
	   ld2 = size(A,1)
	   ld3 = size(Y,1)
	   alpha = 1.0d0
	   beta = 0.0d0

	   call zgemm('N', 'T', mm,nn,kk,
     &          alpha, BX, ld1, A, ld2, beta, Y, ld3 )

	else
	  BX(1:ny1,1:nx2) = matmul( B(1:ny1,1:nx1), 
     &         X(1:nx1,1:nx2))

	  allocate( A_t(1:nx2,1:ny2) )
	  do j=1,ny2
	  do i=1,nx2
	    A_t(i,j) = A(j,i)
	  enddo
	  enddo
	  Y(1:ny1,1:ny2) = matmul( BX(1:ny1, 1:nx2), 
     &         A_t(1:nx2,1:ny2) )

	  deallocate( A_t )
	endif

	deallocate( BX )

	else
!	-----------------------------------
!	compute as Y = B * (X*transpose(A))
!	-----------------------------------
	allocate( XA(nx1,ny2) )

	if (use_gemm) then
	   mm = nx1
	   nn = ny2
	   kk = nx2
	   ld1 = size(X,1)
	   ld2 = size(A,1)
	   ld3 = size(XA,1)
	   alpha = 1.0d0
	   beta = 0.0d0

	   call zgemm( 'N', 'T', mm,nn,kk,
     &           alpha, X, ld1, A, ld2, beta, XA, ld3 )


	   mm = ny1
	   nn = ny2
	   kk = nx1
	   ld1 = size(B,1)
	   ld2 = size(XA,1)
	   ld3 = size(Y,1)
	   alpha = 1.0d0
	   beta = 0.0d0
	 
	   call zgemm( 'N','N', mm,nn,kk,
     &		alpha, B,ld1, XA, ld2, beta, Y, ld3 )
	   

	else

	  allocate( A_t(1:nx2,1:ny2) )
	  do j=1,ny2
	  do i=1,nx2
	    A_t(i,j) = A(j,i)
	  enddo
	  enddo

	  XA(1:nx1,1:ny2) = matmul( X(1:nx1,1:nx2),
     &        A_t(1:nx2,1:ny2)  )

	  deallocate( A_t )

	  Y(1:ny1,1:ny2) = matmul( B(1:ny1,1:nx1),
     &                    XA(1:nx1, 1:ny2) )
	endif

	deallocate( XA )
	endif

	return
	end subroutine kron2
	

	end module kron3_mod

