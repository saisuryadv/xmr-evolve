*> \brief \b DERRBSVR
*
*     SUBROUTINE DERRBSVR( PATH, NUNIT )
*
*     Exercise the input-validation branch of DBDSVR without depending
*     on LAPACK's TESTING infrastructure (CHKXER / XERBLA).  Each
*     invalid-argument case is dispatched and the returned INFO is
*     compared against the expected negative value.
*
      SUBROUTINE DERRBSVR( PATH, NUNIT )
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
      CHARACTER*(*)      PATH
      INTEGER            NUNIT
*     ..
*     .. Parameters ..
      INTEGER            NMAX, LW, LIW
      PARAMETER          ( NMAX = 4, LW = 4*NMAX*NMAX + 37*NMAX,
     $                     LIW = 20*NMAX )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, NFAIL, NS
      DOUBLE PRECISION   VL, VU
*     ..
*     .. Local Arrays ..
      INTEGER            IW( LIW )
      DOUBLE PRECISION   D( NMAX ), E( NMAX ), S( NMAX ), W( LW ),
     $                   Z( 2*NMAX, NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DBDSVR
*     ..
*     .. Executable Statements ..
      NFAIL = 0
      DO 5 I = 1, NMAX
         D( I ) = ONE
         E( I ) = ZERO
         S( I ) = ZERO
    5 CONTINUE
      VL = ZERO
      VU = ONE
*
      CALL DBDSVR( 'X', 'N', 'A', 2, D, E, VL, VU, 1, 2, NS, S,
     $             Z, 4, W, LW, IW, LIW, INFO )
      IF( INFO.NE.-1 ) NFAIL = NFAIL + 1
      CALL DBDSVR( 'U', 'X', 'A', 2, D, E, VL, VU, 1, 2, NS, S,
     $             Z, 4, W, LW, IW, LIW, INFO )
      IF( INFO.NE.-2 ) NFAIL = NFAIL + 1
      CALL DBDSVR( 'U', 'N', 'X', 2, D, E, VL, VU, 1, 2, NS, S,
     $             Z, 4, W, LW, IW, LIW, INFO )
      IF( INFO.NE.-3 ) NFAIL = NFAIL + 1
      CALL DBDSVR( 'U', 'N', 'A', -1, D, E, VL, VU, 1, 1, NS, S,
     $             Z, 4, W, LW, IW, LIW, INFO )
      IF( INFO.NE.-4 ) NFAIL = NFAIL + 1
      CALL DBDSVR( 'U', 'N', 'V', 2, D, E, -ONE, ONE, 1, 2, NS, S,
     $             Z, 4, W, LW, IW, LIW, INFO )
      IF( INFO.NE.-7 ) NFAIL = NFAIL + 1
      CALL DBDSVR( 'U', 'N', 'V', 2, D, E, ONE, ZERO, 1, 2, NS, S,
     $             Z, 4, W, LW, IW, LIW, INFO )
      IF( INFO.NE.-8 ) NFAIL = NFAIL + 1
      CALL DBDSVR( 'U', 'N', 'I', 2, D, E, VL, VU, 0, 2, NS, S,
     $             Z, 4, W, LW, IW, LIW, INFO )
      IF( INFO.NE.-9 ) NFAIL = NFAIL + 1
      CALL DBDSVR( 'U', 'N', 'I', 2, D, E, VL, VU, 2, 1, NS, S,
     $             Z, 4, W, LW, IW, LIW, INFO )
      IF( INFO.NE.-10 ) NFAIL = NFAIL + 1
      CALL DBDSVR( 'U', 'V', 'A', 2, D, E, VL, VU, 1, 2, NS, S,
     $             Z, 1, W, LW, IW, LIW, INFO )
      IF( INFO.NE.-14 ) NFAIL = NFAIL + 1
*
      IF( NFAIL.EQ.0 ) THEN
         WRITE( NUNIT, FMT = 9999 ) PATH
      ELSE
         WRITE( NUNIT, FMT = 9998 ) PATH, NFAIL
      END IF
*
 9999 FORMAT( 1X, A3, ' argument error checks passed.' )
 9998 FORMAT( 1X, A3, ' *** ', I3, ' ARGUMENT ERROR CHECKS FAILED ***' )
      RETURN
      END
*
*     Local XERBLA that lets DERRBSVR probe the illegal-argument branch
*     without aborting.  Overrides the LAPACK library XERBLA (which
*     STOPs) because this object is linked before the static archives.
*
      SUBROUTINE XERBLA( SRNAME, INFO )
      IMPLICIT NONE
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
*     Silently return; the caller reads INFO instead.
      RETURN
      END
