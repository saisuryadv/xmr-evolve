*> \brief \b DBSVRMG -- best-recoverable 29-case LAW 166 generator
*
*  This generator reconstructs the 29-case robustness family behind the
*  bidiagonal-SVD work from the Groesser/Lang DMATGEN catalog: twelve
*  prescribed-spectrum cases (IDs 110--121) and seventeen entry-defined
*  cases (IDs 200--244).  LAW 166 does not publish the Figure 6.1 input
*  manifest, so this is deliberately described as a reconstruction rather
*  than as the unpublished original input deck.
*
*  Spectrum-defined and symmetric-tridiagonal cases are converted to an
*  upper bidiagonal B by forming T-shift*I = B**T B.  Direct bidiagonal
*  cases are returned unchanged.  The formulas and fixed LAPACK random
*  seed conventions follow the existing test-suite documentation.
*
      SUBROUTINE DBSVRMG( JTYPE, N, ISEED, BD, BE, UPLO,
     $                    A, LDA, D, WORK, INFO )
      IMPLICIT NONE
*
      CHARACTER          UPLO
      INTEGER            INFO, JTYPE, LDA, N
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   A( LDA, * ), BD( * ), BE( * ), D( * ),
     $                   WORK( * )
*
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF, BETA, SHIFT_FACTOR
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0,
     $                     TWO = 2.0D0, HALF = 0.5D0,
     $                     BETA = 0.5D0, SHIFT_FACTOR = 1.001D0 )
      INTEGER            I, HALF_N
      DOUBLE PRECISION   ALPHA, CENTER, EPS, EPS14, GL, RND,
     $                   SHIFT, SCALE, TEMP
      DOUBLE PRECISION   DLAMCH, DLARAN
      EXTERNAL           DLAMCH, DLARAN
      EXTERNAL           DLATMS, DPTTRF
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*
      INFO = 0
      UPLO = 'U'
      EPS = DLAMCH( 'Epsilon' )
      EPS14 = SQRT( SQRT( EPS ) )
      ALPHA = EPS
      IF( N.LT.2 ) THEN
         INFO = -2
         RETURN
      END IF
      IF( JTYPE.LT.1 .OR. JTYPE.GT.29 ) THEN
         INFO = -1
         RETURN
      END IF
      DO 5 I = 1, N
         BD( I ) = ZERO
         D( I ) = ZERO
    5 CONTINUE
      DO 6 I = 1, N-1
         BE( I ) = ZERO
    6 CONTINUE
*
*     Types 1--12 correspond to DMATGEN IDs 110--121.  D contains
*     eigenvalues of a symmetric tridiagonal T, including the signed
*     variants.  DLATMS reduces the diagonal spectrum to tridiagonal
*     form before the common Cholesky lift below.
      IF( JTYPE.LE.12 ) THEN
         IF( JTYPE.EQ.1 ) THEN
            DO 10 I = 1, N
               D( I ) = ONE
   10       CONTINUE
         ELSE IF( JTYPE.EQ.2 ) THEN
            DO 11 I = 1, N-1
               D( I ) = DBLE( I ) * EPS
   11       CONTINUE
            D( N ) = ONE
         ELSE IF( JTYPE.EQ.3 ) THEN
            D( 1 ) = EPS
            DO 12 I = 2, N-1
               D( I ) = ONE + DBLE( I ) * EPS14
   12       CONTINUE
            D( N ) = TWO
         ELSE IF( JTYPE.EQ.4 .OR. JTYPE.EQ.5 ) THEN
            DO 13 I = 1, N
               D( I ) = EPS + DBLE( I-1 ) * ( ONE-EPS ) /
     $                  DBLE( N-1 )
               IF( JTYPE.EQ.5 ) THEN
                  RND = DLARAN( ISEED )
                  IF( RND.LT.HALF ) D( I ) = -D( I )
               END IF
   13       CONTINUE
         ELSE IF( JTYPE.EQ.6 .OR. JTYPE.EQ.7 ) THEN
            DO 14 I = 1, N
               D( I ) = EPS**( DBLE( N-I ) / DBLE( N-1 ) )
               IF( JTYPE.EQ.7 ) THEN
                  RND = DLARAN( ISEED )
                  IF( RND.LT.HALF ) D( I ) = -D( I )
               END IF
   14       CONTINUE
         ELSE IF( JTYPE.EQ.8 ) THEN
            DO 15 I = 1, N
               D( I ) = DLARAN( ISEED )
   15       CONTINUE
         ELSE IF( JTYPE.EQ.9 ) THEN
            D( 1 ) = EPS*EPS
            SCALE = 10.0D0 * DBLE( N )
            DO 16 I = 2, N
               D( I ) = ( SCALE + DLARAN( ISEED ) ) / SCALE
   16       CONTINUE
         ELSE IF( JTYPE.EQ.10 ) THEN
            D( 1 ) = EPS
            DO 17 I = 2, N
               D( I ) = ONE
               IF( DLARAN( ISEED ).LT.HALF ) D( I ) = -ONE
   17       CONTINUE
         ELSE IF( JTYPE.EQ.11 ) THEN
            SCALE = 10.0D0 * DBLE( N )
            DO 18 I = 1, N-1
               D( I ) = EPS*EPS *
     $                  ( SCALE + DLARAN( ISEED ) ) / SCALE
   18       CONTINUE
            D( N ) = ONE
         ELSE
            DO 19 I = 1, N-1
               D( I ) = EPS
               IF( DLARAN( ISEED ).LT.HALF ) D( I ) = -EPS
   19       CONTINUE
            D( N ) = ONE
         END IF
         CALL DLATMS( N, N, 'S', ISEED, 'S', D, 0, ONE, ONE,
     $                1, 1, 'N', A, LDA, WORK, INFO )
         IF( INFO.NE.0 ) THEN
            INFO = 100 + INFO
            RETURN
         END IF
         DO 20 I = 1, N
            BD( I ) = A( I, I )
   20    CONTINUE
         DO 21 I = 1, N-1
            BE( I ) = A( I, I+1 )
   21    CONTINUE
         GO TO 500
      END IF
*
*     Types 13--16: DMATGEN IDs 200--203, ABCON tridiagonals.
      IF( JTYPE.GE.13 .AND. JTYPE.LE.16 ) THEN
         DO 30 I = 1, N
            BD( I ) = TWO
   30    CONTINUE
         DO 31 I = 1, N-1
            BE( I ) = ONE
   31    CONTINUE
         IF( JTYPE.EQ.14 .OR. JTYPE.EQ.15 ) BD( 1 ) = ONE
         IF( JTYPE.EQ.15 ) BD( N ) = 3.0D0
         IF( JTYPE.EQ.16 ) THEN
            BD( 1 ) = 3.0D0
            BD( N ) = 3.0D0
         END IF
         GO TO 500
      END IF
*
*     Type 17: ID 210, random bidiagonal entries.
      IF( JTYPE.EQ.17 ) THEN
         DO 40 I = 1, N
            BD( I ) = DLARAN( ISEED )
   40    CONTINUE
         DO 41 I = 1, N-1
            BE( I ) = DLARAN( ISEED )
   41    CONTINUE
         RETURN
      END IF
*
*     Types 18--19: IDs 220--221, graded direct bidiagonals.
      IF( JTYPE.EQ.18 ) THEN
         BD( N ) = ONE
         DO 50 I = N-1, 1, -1
            BD( I ) = BETA * BD( I+1 )
            BE( I ) = BD( I )
   50    CONTINUE
         RETURN
      ELSE IF( JTYPE.EQ.19 ) THEN
         BD( 1 ) = ONE
         DO 51 I = 1, N-1
            BD( I+1 ) = BETA * BD( I )
            BE( I ) = BD( I+1 )
   51    CONTINUE
         RETURN
      END IF
*
*     Types 20--21: IDs 222--223, Wilkinson +/- tridiagonals.
      IF( JTYPE.EQ.20 .OR. JTYPE.EQ.21 ) THEN
         CENTER = DBLE( N+1 ) / TWO
         DO 60 I = 1, N
            TEMP = CENTER - DBLE( I )
            IF( JTYPE.EQ.20 ) THEN
               BD( I ) = ABS( TEMP )
            ELSE
               BD( I ) = TEMP
            END IF
   60    CONTINUE
         DO 61 I = 1, N-1
            BE( I ) = ONE
   61    CONTINUE
         GO TO 500
      END IF
*
*     Types 22--23: IDs 224--225, W and double-W direct bidiagonals.
      IF( JTYPE.EQ.22 ) THEN
         DO 70 I = 1, N
            BD( I ) = ABS( ONE - DBLE( N ) / TWO + DBLE( I ) )
   70    CONTINUE
         DO 71 I = 1, N-1
            BE( I ) = ONE
   71    CONTINUE
         RETURN
      ELSE IF( JTYPE.EQ.23 ) THEN
         HALF_N = N / 2
         DO 72 I = 1, HALF_N
            BD( I ) = ABS( ONE - DBLE( N ) / 4.0D0 + DBLE( I ) )
   72    CONTINUE
         DO 73 I = HALF_N+1, N
            BD( I ) = BD( I-HALF_N )
   73    CONTINUE
         DO 74 I = 1, N-1
            BE( I ) = ONE
   74    CONTINUE
         RETURN
      END IF
*
*     Type 24: ID 230, Clement tridiagonal.
      IF( JTYPE.EQ.24 ) THEN
         DO 80 I = 1, N
            BD( I ) = ZERO
   80    CONTINUE
         DO 81 I = 1, N-1
            BE( I ) = SQRT( DBLE( I*(N-I) ) )
   81    CONTINUE
         GO TO 500
      END IF
*
*     Types 25--28: IDs 240--243, direct GRO bidiagonals.
      IF( JTYPE.GE.25 .AND. JTYPE.LE.28 ) THEN
         DO 90 I = 1, N
            BD( I ) = ALPHA
   90    CONTINUE
         DO 91 I = 1, N-1
            BE( I ) = ALPHA
   91    CONTINUE
         BD( 1 ) = ONE
         IF( JTYPE.GE.26 .AND. N.GE.2 ) BD( 2 ) = ONE
         IF( JTYPE.GE.27 ) THEN
            DO 92 I = 1, MIN( 4, N )
               BD( I ) = ONE
   92       CONTINUE
         END IF
         IF( JTYPE.EQ.28 ) THEN
            DO 93 I = 1, MIN( 4, N-1 )
               BE( I ) = ONE
   93       CONTINUE
         END IF
         RETURN
      END IF
*
*     Type 29: ID 244, shifted Wilkinson+ direct bidiagonal.
      CENTER = DBLE( N+1 ) / TWO
      DO 100 I = 1, N
         BD( I ) = ABS( CENTER - DBLE( I ) ) + ONE
  100 CONTINUE
      DO 101 I = 1, N-1
         BE( I ) = ONE
  101 CONTINUE
      RETURN
*
*     Common symmetric-tridiagonal to upper-bidiagonal conversion.
  500 CONTINUE
      GL = BD( 1 ) - ABS( BE( 1 ) )
      DO 510 I = 2, N-1
         GL = MIN( GL, BD( I )-ABS( BE( I-1 ) )-ABS( BE( I ) ) )
  510 CONTINUE
      GL = MIN( GL, BD( N )-ABS( BE( N-1 ) ) )
      SHIFT = ZERO
      IF( GL.LT.ZERO ) SHIFT = SHIFT_FACTOR * GL
      DO 511 I = 1, N
         D( I ) = BD( I )
         BD( I ) = BD( I ) - SHIFT
  511 CONTINUE
      DO 512 I = 1, N-1
         WORK( I ) = BE( I )
  512 CONTINUE
      CALL DPTTRF( N, BD, BE, INFO )
      IF( INFO.NE.0 ) THEN
*        Gershgorin can equal an eigenvalue for reducible/semidefinite
*        inputs.  Move one conservative unit farther left and retry.
         DO 513 I = 1, N
            BD( I ) = D( I ) - ( SHIFT-ONE )
  513    CONTINUE
         DO 514 I = 1, N-1
            BE( I ) = WORK( I )
  514    CONTINUE
         CALL DPTTRF( N, BD, BE, INFO )
         IF( INFO.NE.0 ) THEN
            INFO = 300 + INFO
            RETURN
         END IF
      END IF
      DO 520 I = 1, N
         BD( I ) = SQRT( ABS( BD( I ) ) )
  520 CONTINUE
      DO 521 I = 1, N-1
         BE( I ) = BE( I ) * BD( I )
  521 CONTINUE
      RETURN
      END
