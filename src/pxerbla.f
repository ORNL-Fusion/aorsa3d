      SUBROUTINE PXERBLA( ICTXT, SRNAME, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 20, 1995
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, INFO
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      SRNAME
*     ..
*
*  Purpose
*  =======
*
*  PXERBLA is an error handler for the ScaLAPACK routines. It is called
*  by a ScaLAPACK routine if an input parameter has an invalid value.
*  A message is printed.  Installers may consider modifying this routine
*  in order to call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation. The context itself is global.
*
*  SRNAME  (global input) CHARACTER*(*)
*          The name of the routine which called PXERBLA.
*
*  INFO    (global input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      WRITE( *, FMT = 9999 ) MYROW, MYCOL, SRNAME, INFO
*
 9999 FORMAT( '{', I5, ',', I5, '}:  On entry to ', A,
     $        ' parameter number', I4, ' had an illegal value' )
*
*     End of PXERBLA
*
      END
