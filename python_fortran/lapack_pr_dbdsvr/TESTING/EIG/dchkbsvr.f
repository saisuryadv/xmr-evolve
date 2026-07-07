*> \brief \b DCHKBSVR
*
*     SUBROUTINE DCHKBSVR( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
*                          NOUT, TSTREF, ...)
*
*     Bidiagonal-SVD tester -- lifts the DBDSVDX residual-and-metric
*     block from LAPACK's dchkbd.f *verbatim* and applies it to DBDSVR,
*     and optionally to DBDSVDX itself.
*
*     Metrics (per matrix, exactly as in dchkbd.f):
*       RANGE='A':  RESULT(20) = || SA - U^T B VT^T || / (n |B| ulp)  [DBDT04]
*                   RESULT(21) = || I - U^T U || / (n ulp)            [DORT01]
*                   RESULT(22) = || I - VT VT^T || / (n ulp)          [DORT01]
*                   RESULT(23) = 1/ulp if SA not non-increasing / >= 0
*                   RESULT(24) = max_j |SA(j)-S2(j)| / (dchkbd formula)
*       RANGE='I':  RESULT(25..29) same shape.
*       RANGE='V':  RESULT(30..34) same shape.
*
*     VL, VU for the RANGE='V' block are computed from the RANGE='A'
*     spectrum SA per dchkbd.f:1374-1388 (byte-identical formula).
*
*     A residual R(k) >= THRESH is a failure.  Reference default is
*     THRESH = 20.0 (mirrors dchkbd's bd.in).
*
      SUBROUTINE DCHKBSVR( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
     $                     NOUT, TSTREF, BD, BE, BD1, BE1, S1, S2, SA,
     $                     U, LDU, VT, LDVT, Z, LDZ, A, LDA, DTMP,
     $                     WORK, LWORK, IWORK, LIWORK, INFO )
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
      LOGICAL            TSTREF
      INTEGER            INFO, LDA, LDU, LDVT, LDZ, LIWORK, LWORK,
     $                   NOUT, NSIZES, NTYPES
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
      DOUBLE PRECISION   A( LDA, * ), BD( * ), BD1( * ), BE( * ),
     $                   BE1( * ), DTMP( * ), S1( * ), S2( * ),
     $                   SA( * ), U( LDU, * ), VT( LDVT, * ),
     $                   WORK( * ), Z( LDZ, * )
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, HALF, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0,
     $                     TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          UPLO
      INTEGER            I, IINFO, IL, IU, J, JSIZE, JTYPE, K, LWRK,
     $                   METH, N, NFAIL, NRUN, NSA, NS1, NS2
      DOUBLE PRECISION   ANORM, RTUNFL, TEMP1, TEMP2, ULP, ULPINV,
     $                   UNFL, VL, VU
*     ..
*     .. Local Arrays ..
      CHARACTER*8        METHNM( 2 )
      DOUBLE PRECISION   R( 20:34 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DBDSVDX, DBDSVR, DBDT04, DBSVRMG, DCOPY,
     $                   DORT01
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Data statements ..
      DATA               METHNM / 'DBDSVR  ', 'DBDSVDX ' /
*     ..
*     .. Executable Statements ..
      INFO   = 0
      NFAIL  = 0
      NRUN   = 0
      ULP    = DLAMCH( 'Precision' )
      UNFL   = DLAMCH( 'Safe minimum' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      LWRK   = LWORK
*
      DO 200 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         IF( N.LT.2 ) GO TO 200
*
         DO 190 JTYPE = 1, NTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 190
*
*           Generate B once and re-use for both solvers.
*
            CALL DBSVRMG( JTYPE, N, ISEED, BD, BE, UPLO, A, LDA, DTMP,
     $                    WORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9995 ) 'DBSVRMG', IINFO, JSIZE, JTYPE
               INFO = INFO + 1
               GO TO 190
            END IF
*
*           Compute ||B||.
            ANORM = ABS( BD( 1 ) )
            DO 10 I = 2, N
               ANORM = MAX( ANORM, ABS( BD( I ) ) + ABS( BE( I-1 ) ) )
   10       CONTINUE
            IF( ANORM.EQ.ZERO ) ANORM = ONE
*
*           Pick IL, IU (deterministic mid slice for reproducibility).
            IL = MAX( 1, N/4 + 1 )
            IU = MIN( N, ( 3*N ) / 4 )
            IF( IU.LT.IL ) THEN
               IL = 1
               IU = N
            END IF
*
*           For each solver (DBDSVR then optionally DBDSVDX): run the
*           full RANGE='A', RANGE='I', RANGE='V' residual battery.
*
            DO 180 METH = 1, 2
               IF( METH.EQ.2 .AND. .NOT.TSTREF ) GO TO 180
*
*              Zero out residuals.
               DO 15 I = 20, 34
                  R( I ) = ZERO
   15          CONTINUE
*
*              -------- RANGE='A': tests 20..24 --------
*
               CALL DCOPY( N, BD, 1, BD1, 1 )
               IF( N.GT.1 ) CALL DCOPY( N-1, BE, 1, BE1, 1 )
*
               IF( METH.EQ.1 ) THEN
                  CALL DBDSVR( UPLO, 'V', 'A', N, BD1, BE1, ZERO, ZERO,
     $                         0, 0, NSA, SA, Z, LDZ, WORK, LWRK,
     $                         IWORK, LIWORK, IINFO )
               ELSE
                  CALL DBDSVDX( UPLO, 'V', 'A', N, BD1, BE1, ZERO, ZERO,
     $                          0, 0, NSA, SA, Z, LDZ, WORK, IWORK,
     $                          IINFO )
               END IF
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9995 ) METHNM( METH )//'V/A',
     $                 IINFO, JSIZE, JTYPE
                  R( 20 ) = ULPINV
                  GO TO 150
               END IF
*
               DO 20 I = 1, NSA
                  CALL DCOPY( N, Z( 1  , I ), 1, U( 1, I ), 1 )
                  CALL DCOPY( N, Z( N+1, I ), 1, VT( I, 1 ), LDVT )
   20          CONTINUE
*
*              Companion N/A call for R(24) comparison.
               CALL DCOPY( N, BD, 1, BD1, 1 )
               IF( N.GT.1 ) CALL DCOPY( N-1, BE, 1, BE1, 1 )
               IF( METH.EQ.1 ) THEN
                  CALL DBDSVR( UPLO, 'N', 'A', N, BD1, BE1, ZERO, ZERO,
     $                         0, 0, NS2, S2, Z, LDZ, WORK, LWRK,
     $                         IWORK, LIWORK, IINFO )
               ELSE
                  CALL DBDSVDX( UPLO, 'N', 'A', N, BD1, BE1, ZERO, ZERO,
     $                          0, 0, NS2, S2, Z, LDZ, WORK, IWORK,
     $                          IINFO )
               END IF
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9995 ) METHNM( METH )//'N/A',
     $                 IINFO, JSIZE, JTYPE
                  R( 24 ) = ULPINV
                  GO TO 150
               END IF
*
               CALL DBDT04( UPLO, N, BD, BE, SA, NSA, U, LDU, VT, LDVT,
     $                      WORK, R( 20 ) )
               CALL DORT01( 'Columns', N, NSA, U, LDU, WORK, LWRK,
     $                      R( 21 ) )
               CALL DORT01( 'Rows',    NSA, N, VT, LDVT, WORK, LWRK,
     $                      R( 22 ) )
*
               DO 25 I = 1, NSA - 1
                  IF( SA( I ).LT.SA( I+1 ) ) R( 23 ) = ULPINV
                  IF( SA( I ).LT.ZERO )      R( 23 ) = ULPINV
   25          CONTINUE
               IF( NSA.GE.1 ) THEN
                  IF( SA( NSA ).LT.ZERO )    R( 23 ) = ULPINV
               END IF
*
               TEMP2 = ZERO
               IF( NSA.EQ.NS2 ) THEN
                  DO 27 J = 1, NSA
                     TEMP1 = ABS( SA( J ) - S2( J ) ) /
     $                       MAX( RTUNFL*MAX( SA( 1 ), ONE ),
     $                            ULP*MAX( ABS( SA( 1 ) ),
     $                                     ABS( S2( 1 ) ) ) )
                     TEMP2 = MAX( TEMP1, TEMP2 )
   27             CONTINUE
               ELSE
                  TEMP2 = ULPINV
               END IF
               R( 24 ) = TEMP2
*
*              -------- RANGE='I': tests 25..29 --------
*
               CALL DCOPY( N, BD, 1, BD1, 1 )
               IF( N.GT.1 ) CALL DCOPY( N-1, BE, 1, BE1, 1 )
*
               IF( METH.EQ.1 ) THEN
                  CALL DBDSVR( UPLO, 'V', 'I', N, BD1, BE1, ZERO, ZERO,
     $                         IL, IU, NS1, S1, Z, LDZ, WORK, LWRK,
     $                         IWORK, LIWORK, IINFO )
               ELSE
                  CALL DBDSVDX( UPLO, 'V', 'I', N, BD1, BE1, ZERO, ZERO,
     $                          IL, IU, NS1, S1, Z, LDZ, WORK, IWORK,
     $                          IINFO )
               END IF
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9995 ) METHNM( METH )//'V/I',
     $                 IINFO, JSIZE, JTYPE
                  R( 25 ) = ULPINV
                  GO TO 150
               END IF
*
*              Unpack Z -> U, VT (mirrors dchkbd:1296-1301).
               DO 30 I = 1, NS1
                  CALL DCOPY( N, Z( 1  , I ), 1, U( 1, I ), 1 )
                  CALL DCOPY( N, Z( N+1, I ), 1, VT( I, 1 ), LDVT )
   30          CONTINUE
*
*              Companion N/I call (values only) for S2.
               CALL DCOPY( N, BD, 1, BD1, 1 )
               IF( N.GT.1 ) CALL DCOPY( N-1, BE, 1, BE1, 1 )
               IF( METH.EQ.1 ) THEN
                  CALL DBDSVR( UPLO, 'N', 'I', N, BD1, BE1, ZERO, ZERO,
     $                         IL, IU, NS2, S2, Z, LDZ, WORK, LWRK,
     $                         IWORK, LIWORK, IINFO )
               ELSE
                  CALL DBDSVDX( UPLO, 'N', 'I', N, BD1, BE1, ZERO, ZERO,
     $                          IL, IU, NS2, S2, Z, LDZ, WORK, IWORK,
     $                          IINFO )
               END IF
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9995 ) METHNM( METH )//'N/I',
     $                 IINFO, JSIZE, JTYPE
                  R( 29 ) = ULPINV
                  GO TO 150
               END IF
*
*              R(25) = DBDT04 residual, R(26) = DORT01 U, R(27) = DORT01 VT.
               CALL DBDT04( UPLO, N, BD, BE, S1, NS1, U, LDU, VT, LDVT,
     $                      WORK, R( 25 ) )
               CALL DORT01( 'Columns', N, NS1, U, LDU, WORK, LWRK,
     $                      R( 26 ) )
               CALL DORT01( 'Rows',    NS1, N, VT, LDVT, WORK, LWRK,
     $                      R( 27 ) )
*
*              R(28) = sortedness / non-negativity (mirrors dchkbd:1347-57).
               DO 40 I = 1, NS1 - 1
                  IF( S1( I ).LT.S1( I+1 ) ) R( 28 ) = ULPINV
                  IF( S1( I ).LT.ZERO )      R( 28 ) = ULPINV
   40          CONTINUE
               IF( NS1.GE.1 ) THEN
                  IF( S1( NS1 ).LT.ZERO )    R( 28 ) = ULPINV
               END IF
*
*              R(29) = max |S1(j)-S2(j)| / (dchkbd:1360-66 formula).
               TEMP2 = ZERO
               IF( NS1.EQ.NS2 ) THEN
                  DO 50 J = 1, NS1
                     TEMP1 = ABS( S1( J ) - S2( J ) ) /
     $                       MAX( RTUNFL*MAX( S1( 1 ), ONE ),
     $                            ULP*MAX( ABS( S1( 1 ) ),
     $                                     ABS( S2( 1 ) ) ) )
                     TEMP2 = MAX( TEMP1, TEMP2 )
   50             CONTINUE
               ELSE
                  TEMP2 = ULPINV
               END IF
               R( 29 ) = TEMP2
*
*              -------- RANGE='V': tests 30..34 --------
*
*              VL, VU picked from SA (the RANGE='A' spectrum) exactly per
*              dchkbd.f:1374-1388.  SA is in descending order and here
*              SA(k) is the k-th largest singular value, so SA(IL) is
*              the IL-th largest.
*
               IF( IL.NE.1 ) THEN
                  TEMP1 = HALF*ABS( SA( IL ) - SA( IL-1 ) )
                  VU = SA( IL ) + MAX( TEMP1, ULP*ANORM, TWO*RTUNFL )
               ELSE
                  TEMP1 = HALF*ABS( SA( NSA ) - SA( 1 ) )
                  VU = SA( 1 ) + MAX( TEMP1, ULP*ANORM, TWO*RTUNFL )
               END IF
               IF( IU.NE.NSA ) THEN
                  TEMP1 = HALF*ABS( SA( IU+1 ) - SA( IU ) )
                  VL = SA( IU ) - MAX( ULP*ANORM, TWO*RTUNFL, TEMP1 )
               ELSE
                  TEMP1 = HALF*ABS( SA( NSA ) - SA( 1 ) )
                  VL = SA( NSA ) - MAX( ULP*ANORM, TWO*RTUNFL, TEMP1 )
               END IF
               VL = MAX( VL, ZERO )
               VU = MAX( VU, ZERO )
               IF( VL.GE.VU ) VU = MAX( VU*TWO, VU + VL + HALF )
*
               CALL DCOPY( N, BD, 1, BD1, 1 )
               IF( N.GT.1 ) CALL DCOPY( N-1, BE, 1, BE1, 1 )
               IF( METH.EQ.1 ) THEN
                  CALL DBDSVR( UPLO, 'V', 'V', N, BD1, BE1, VL, VU,
     $                         1, N, NS1, S1, Z, LDZ, WORK, LWRK,
     $                         IWORK, LIWORK, IINFO )
               ELSE
                  CALL DBDSVDX( UPLO, 'V', 'V', N, BD1, BE1, VL, VU,
     $                          0, 0, NS1, S1, Z, LDZ, WORK, IWORK,
     $                          IINFO )
               END IF
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9995 ) METHNM( METH )//'V/V',
     $                 IINFO, JSIZE, JTYPE
                  R( 30 ) = ULPINV
                  GO TO 150
               END IF
*
               DO 60 I = 1, NS1
                  CALL DCOPY( N, Z( 1  , I ), 1, U( 1, I ), 1 )
                  CALL DCOPY( N, Z( N+1, I ), 1, VT( I, 1 ), LDVT )
   60          CONTINUE
*
               CALL DCOPY( N, BD, 1, BD1, 1 )
               IF( N.GT.1 ) CALL DCOPY( N-1, BE, 1, BE1, 1 )
               IF( METH.EQ.1 ) THEN
                  CALL DBDSVR( UPLO, 'N', 'V', N, BD1, BE1, VL, VU,
     $                         1, N, NS2, S2, Z, LDZ, WORK, LWRK,
     $                         IWORK, LIWORK, IINFO )
               ELSE
                  CALL DBDSVDX( UPLO, 'N', 'V', N, BD1, BE1, VL, VU,
     $                          0, 0, NS2, S2, Z, LDZ, WORK, IWORK,
     $                          IINFO )
               END IF
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9995 ) METHNM( METH )//'N/V',
     $                 IINFO, JSIZE, JTYPE
                  R( 34 ) = ULPINV
                  GO TO 150
               END IF
*
               CALL DBDT04( UPLO, N, BD, BE, S1, NS1, U, LDU, VT, LDVT,
     $                      WORK, R( 30 ) )
               CALL DORT01( 'Columns', N, NS1, U, LDU, WORK, LWRK,
     $                      R( 31 ) )
               CALL DORT01( 'Rows',    NS1, N, VT, LDVT, WORK, LWRK,
     $                      R( 32 ) )
*
               DO 70 I = 1, NS1 - 1
                  IF( S1( I ).LT.S1( I+1 ) ) R( 33 ) = ULPINV
                  IF( S1( I ).LT.ZERO )      R( 33 ) = ULPINV
   70          CONTINUE
               IF( NS1.GE.1 ) THEN
                  IF( S1( NS1 ).LT.ZERO )    R( 33 ) = ULPINV
               END IF
*
               TEMP2 = ZERO
               IF( NS1.EQ.NS2 ) THEN
                  DO 80 J = 1, NS1
                     TEMP1 = ABS( S1( J ) - S2( J ) ) /
     $                       MAX( RTUNFL*MAX( S1( 1 ), ONE ),
     $                            ULP*MAX( ABS( S1( 1 ) ),
     $                                     ABS( S2( 1 ) ) ) )
                     TEMP2 = MAX( TEMP1, TEMP2 )
   80             CONTINUE
               ELSE
                  TEMP2 = ULPINV
               END IF
               R( 34 ) = TEMP2
*
  150          CONTINUE
*
*              Tally.
*
               DO 160 K = 20, 34
                  NRUN = NRUN + 1
                  IF( R( K ).GE.THRESH ) THEN
                     NFAIL = NFAIL + 1
                     WRITE( NOUT, FMT = 9999 ) METHNM( METH ),
     $                    N, JTYPE, K, R( K )
                  END IF
  160          CONTINUE
               WRITE( NOUT, FMT = 9998 ) METHNM( METH ), N, JTYPE,
     $              R( 20 ), R( 21 ), R( 22 ), R( 23 ), R( 24 ),
     $              R( 25 ), R( 26 ), R( 27 ), R( 28 ), R( 29 ),
     $              R( 30 ), R( 31 ), R( 32 ), R( 33 ), R( 34 )
*
  180       CONTINUE
*
  190    CONTINUE
  200 CONTINUE
*
      WRITE( NOUT, FMT = 9997 ) NRUN, NFAIL
      INFO = INFO + NFAIL
      RETURN
*
 9999 FORMAT( ' FAIL  ', A8, ' N=', I5, ' TYPE=', I2, ' TEST(',
     $        I2, ') =', ES10.2 )
 9998 FORMAT( ' ', A8, ' N=', I5, ' TYPE=', I2,
     $        ' A:', 5(1X,ES9.2), '  I:', 5(1X,ES9.2),
     $        '  V:', 5(1X,ES9.2) )
 9997 FORMAT( /1X, 'DCHKBSVR summary: ', I5, ' residuals tested, ',
     $        I5, ' failures at THRESH.' )
 9995 FORMAT( ' ERROR ', A, ' INFO=', I4, ' size#', I3, ' type', I3 )
*
      END
