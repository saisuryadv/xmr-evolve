*     Standalone probe: reproduce the type-4 N=100 U-ortho failure of
*     DBDSVDX in isolation, cycling ISEED forward through prior DBSVRMG
*     calls exactly as DCHKBSVR does.  Then call DBDSVDX with LDZ=2*N
*     (tight) and directly compute ‖I − UᵀU‖ without going through
*     DORT01.  Also test with LDZ=400 to check for LDZ sensitivity.
*
      PROGRAM PROBE
      IMPLICIT NONE
      INTEGER            N, MAXN, LWORK, LIWORK, JT, I, J, INFO, NS,
     $                   NNSIZE( 4 ), NNCUR
      PARAMETER          ( N = 100, MAXN = 200,
     $                     LWORK = 14*MAXN, LIWORK = 12*MAXN )
      CHARACTER          UPLO
      DOUBLE PRECISION   BD( MAXN ), BE( MAXN ), BD1( MAXN )
      DOUBLE PRECISION   BE1( MAXN ), S1( MAXN ), DTMP( MAXN )
      DOUBLE PRECISION   A( MAXN, MAXN ), U( MAXN, MAXN )
      DOUBLE PRECISION   UUT( MAXN, MAXN )
      DOUBLE PRECISION   ZTIGHT( 2*MAXN, MAXN )
      DOUBLE PRECISION   ZLOOSE( 2*MAXN, MAXN )
      DOUBLE PRECISION   WORK( LWORK ), RESID
      INTEGER            IWORK( LIWORK ), ISEED( 4 )
      DOUBLE PRECISION   ZERO, ONE, DLANGE
      EXTERNAL           DBSVRMG, DBDSVDX, DGEMM, DLANGE, DCOPY, DLASET
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*
      ISEED( 1 ) = 1
      ISEED( 2 ) = 2
      ISEED( 3 ) = 3
      ISEED( 4 ) = 4
      NNSIZE( 1 ) = 2
      NNSIZE( 2 ) = 10
      NNSIZE( 3 ) = 30
      NNSIZE( 4 ) = 60
*
*     Advance ISEED through the same sequence of DBSVRMG calls that
*     DCHKBSVR made before hitting (N=100, type 4): all types 1..12
*     at N=2, 10, 30, 60, then types 1..3 at N=100.
      DO 10 J = 1, 4
         NNCUR = NNSIZE( J )
         DO 20 JT = 1, 12
            CALL DBSVRMG( JT, NNCUR, ISEED, BD, BE, UPLO,
     $                    A, MAXN, DTMP, WORK, INFO )
   20    CONTINUE
   10 CONTINUE
      DO 30 JT = 1, 3
         CALL DBSVRMG( JT, 100, ISEED, BD, BE, UPLO,
     $                 A, MAXN, DTMP, WORK, INFO )
   30 CONTINUE
*
*     Now generate the target: N=100, type 4.
      CALL DBSVRMG( 4, N, ISEED, BD, BE, UPLO,
     $              A, MAXN, DTMP, WORK, INFO )
      WRITE(6, *) 'DBSVRMG type 4 INFO=', INFO, ' UPLO=', UPLO
      WRITE(6, *) 'BD(1)=', BD(1)
      WRITE(6, *) 'BD(50)=', BD(50)
      WRITE(6, *) 'BD(100)=', BD(100)
      WRITE(6, *) 'BE(1)=', BE(1)
      WRITE(6, *) 'BE(99)=', BE(99)
*
*     Run DBDSVDX with tight LDZ = 2*N.
      CALL DCOPY( N, BD, 1, BD1, 1 )
      CALL DCOPY( N-1, BE, 1, BE1, 1 )
      CALL DLASET( 'F', 2*N, N, ZERO, ZERO, ZTIGHT, 2*MAXN )
      CALL DBDSVDX( UPLO, 'V', 'A', N, BD1, BE1, ZERO, ZERO,
     $              0, 0, NS, S1, ZTIGHT, 2*N, WORK, IWORK, INFO )
      WRITE(6, *) 'DBDSVDX tight LDZ=2N=200 INFO=', INFO, ' NS=', NS
*
*     Unpack U columns (rows 1..N of ZTIGHT), compute UᵀU.
      DO 40 I = 1, NS
         CALL DCOPY( N, ZTIGHT(1, I), 1, U(1, I), 1 )
   40 CONTINUE
      CALL DGEMM( 'T', 'N', NS, NS, N, ONE, U, MAXN, U, MAXN, ZERO,
     $            UUT, MAXN )
      DO 50 I = 1, NS
         UUT( I, I ) = UUT( I, I ) - ONE
   50 CONTINUE
      RESID = DLANGE( 'F', NS, NS, UUT, MAXN, WORK ) / DBLE( N )
      WRITE(6, *) 'tight LDZ:  |I - UtU|_F / N =', RESID
*
*     Extra: mimic the tester's exact sequence -- RANGE='A' JOBZ='V',
*     then companion RANGE='A' JOBZ='N' (which might trash something),
*     then compute U-ortho.
      CALL DCOPY( N, BD, 1, BD1, 1 )
      CALL DCOPY( N-1, BE, 1, BE1, 1 )
      CALL DLASET( 'F', 2*MAXN, N, ZERO, ZERO, ZLOOSE, 2*MAXN )
      CALL DBDSVDX( UPLO, 'V', 'A', N, BD1, BE1, ZERO, ZERO,
     $              0, 0, NS, S1, ZLOOSE, 2*MAXN, WORK, IWORK, INFO )
      DO 45 I = 1, NS
         CALL DCOPY( N, ZLOOSE(1, I), 1, U(1, I), 1 )
   45 CONTINUE
      CALL DCOPY( N, BD, 1, BD1, 1 )
      CALL DCOPY( N-1, BE, 1, BE1, 1 )
      CALL DBDSVDX( UPLO, 'N', 'A', N, BD1, BE1, ZERO, ZERO,
     $              0, 0, NS, DTMP, ZLOOSE, 2*MAXN, WORK, IWORK, INFO )
      CALL DGEMM( 'T', 'N', NS, NS, N, ONE, U, MAXN, U, MAXN, ZERO,
     $            UUT, MAXN )
      DO 55 I = 1, NS
         UUT( I, I ) = UUT( I, I ) - ONE
   55 CONTINUE
      RESID = DLANGE( 'F', NS, NS, UUT, MAXN, WORK ) / DBLE( N )
      WRITE(6, *) 'tester seq: |I - UtU|_F / N =', RESID
*
*     Repeat with loose LDZ = 2*MAXN = 400 (what my tester uses).
      CALL DCOPY( N, BD, 1, BD1, 1 )
      CALL DCOPY( N-1, BE, 1, BE1, 1 )
      CALL DLASET( 'F', 2*MAXN, N, ZERO, ZERO, ZLOOSE, 2*MAXN )
      CALL DBDSVDX( UPLO, 'V', 'A', N, BD1, BE1, ZERO, ZERO,
     $              0, 0, NS, S1, ZLOOSE, 2*MAXN, WORK, IWORK, INFO )
      WRITE(6, *) 'DBDSVDX loose LDZ=400 INFO=', INFO, ' NS=', NS
      DO 60 I = 1, NS
         CALL DCOPY( N, ZLOOSE(1, I), 1, U(1, I), 1 )
   60 CONTINUE
      CALL DGEMM( 'T', 'N', NS, NS, N, ONE, U, MAXN, U, MAXN, ZERO,
     $            UUT, MAXN )
      DO 70 I = 1, NS
         UUT( I, I ) = UUT( I, I ) - ONE
   70 CONTINUE
      RESID = DLANGE( 'F', NS, NS, UUT, MAXN, WORK ) / DBLE( N )
      WRITE(6, *) 'loose LDZ:  |I - UtU|_F / N =', RESID
*
      STOP
      END
