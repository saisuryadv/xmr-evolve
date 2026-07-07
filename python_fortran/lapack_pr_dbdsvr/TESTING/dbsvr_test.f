*> \brief \b DBSVRT  --  main program for the offline DBDSVR PR
*
*     Reads control parameters from stdin and runs the DBDSVR
*     canonical residual sweep.  Mimics the structure of LAPACK's
*     xeigtstd for the DBDSVDX section.
*
*     stdin format:
*         line 1: title (echoed to output)
*         line 2: NSIZES
*         line 3: NN(1..NSIZES)
*         line 4: THRESH
*         line 5: TSTERR   (T/F -- run DERRBSVR)
*         line 6: TSTREF   (T/F -- also run DBDSVDX on the same matrices)
*         line 7: ISEED(1..4)
*
      PROGRAM DBSVRT
      IMPLICIT NONE
*     ..
*     .. Parameters ..
      INTEGER            MAXN, MAXSIZ, NTYPES, LWORK, LIWORK
      PARAMETER          ( MAXN = 200, MAXSIZ = 20, NTYPES = 16,
     $                     LWORK  = 4*MAXN*MAXN + 40*MAXN + 4096,
     $                     LIWORK = 20*MAXN + 4096 )
      INTEGER            NIN, NOUT
      PARAMETER          ( NIN = 5, NOUT = 6 )
*     ..
*     .. Local Scalars ..
      CHARACTER*80       TITLE
      LOGICAL            TSTERR, TSTREF
      INTEGER            I, INFO, JSIZE, NSIZES
      DOUBLE PRECISION   THRESH
*     ..
*     .. Local Arrays ..
      LOGICAL            DOTYPE( NTYPES )
      INTEGER            ISEED( 4 ), NN( MAXSIZ )
*     ..
*     .. Allocatable workspaces ..
      DOUBLE PRECISION, ALLOCATABLE :: BD(:), BE(:), BD1(:), BE1(:),
     $                                  S1(:), S2(:), SA(:),
     $                                  U(:,:), VT(:,:),
     $                                  Z(:,:), A(:,:), DTMP(:), WORK(:)
      INTEGER, ALLOCATABLE :: IWORK(:)
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCHKBSVR, DERRBSVR
*     ..
*     .. Executable Statements ..
      READ( NIN, FMT = '(A)' ) TITLE
      READ( NIN, FMT = * ) NSIZES
      IF( NSIZES.GT.MAXSIZ ) THEN
         WRITE( NOUT, * ) 'NSIZES too large; increase MAXSIZ.'
         STOP 1
      END IF
      READ( NIN, FMT = * ) ( NN( I ), I = 1, NSIZES )
      READ( NIN, FMT = * ) THRESH
      READ( NIN, FMT = * ) TSTERR
      READ( NIN, FMT = * ) TSTREF
      READ( NIN, FMT = * ) ( ISEED( I ), I = 1, 4 )
*
      DO 5 I = 1, NTYPES
         DOTYPE( I ) = .TRUE.
    5 CONTINUE
*
      WRITE( NOUT, FMT = 9998 ) TITLE
      WRITE( NOUT, FMT = 9997 ) NSIZES, THRESH
      WRITE( NOUT, FMT = 9996 ) ( NN( I ), I = 1, NSIZES )
*
*     Verify max size fits our compile-time upper bound.
      DO 6 I = 1, NSIZES
         IF( NN( I ).GT.MAXN ) THEN
            WRITE( NOUT, * ) 'NN(', I, ') =', NN( I ),
     $           ' exceeds MAXN =', MAXN
            STOP 1
         END IF
    6 CONTINUE
*
*     Argument error checks.
      IF( TSTERR ) THEN
         CALL DERRBSVR( 'BSR', NOUT )
      END IF
*
*     Allocate workspaces sized to MAXN.
*
      ALLOCATE( BD( MAXN ), BE( MAXN ), BD1( MAXN ), BE1( MAXN ),
     $          S1( MAXN ), S2( MAXN ), SA( MAXN ),
     $          U( MAXN, MAXN ), VT( MAXN, MAXN ),
     $          Z( 2*MAXN, MAXN ),
     $          A( MAXN, MAXN ),
     $          DTMP( MAXN ),
     $          WORK( LWORK ),
     $          IWORK( LIWORK ) )
*
      CALL DCHKBSVR( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
     $               NOUT, TSTREF, BD, BE, BD1, BE1, S1, S2, SA,
     $               U, MAXN, VT, MAXN, Z, 2*MAXN,
     $               A, MAXN, DTMP, WORK, LWORK,
     $               IWORK, LIWORK, INFO )
*
      DEALLOCATE( BD, BE, BD1, BE1, S1, S2, SA, U, VT, Z, A, DTMP,
     $            WORK, IWORK )
*
      IF( INFO.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9995 )
         STOP 0
      ELSE
         WRITE( NOUT, FMT = 9994 ) INFO
         STOP 1
      END IF
*
 9998 FORMAT( /1X, 'DBDSVR tester: ', A )
 9997 FORMAT( 1X, 'Sizes: ', I3, '   THRESH = ', F8.2 )
 9996 FORMAT( 1X, 'N =', 20(1X, I5) )
 9995 FORMAT( /1X, '*** All tests passed ***' )
 9994 FORMAT( /1X, '*** ', I5, ' tests failed ***' )
      END
