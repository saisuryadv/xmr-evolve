*> \brief \b DBSVRMG — bidiagonal-SVD test matrix generator
*
*  =========================================================
*
*     SUBROUTINE DBSVRMG( JTYPE, N, ISEED, BD, BE, UPLO, A, LDA,
*                         D, WORK, INFO )
*
*     Generate a real N-by-N bidiagonal matrix B (returned as BD, BE)
*     for the DBDSVR canonical test suite.  Types 1-5 build a
*     bidiagonal directly; types 6-10 build a symmetric tridiagonal
*     first, then convert to bidiagonal via a shifted L*D*L^T
*     factorization (Cholesky) — B^T B = T_shift.  Types 11-12 apply
*     magnitude scaling to a base spectrum.
*
*     |Type| Source                             | Recipe                                     |
*     |----|------------------------------------|--------------------------------------------|
*     |  1 | dchkbd type 1                      | Zero bidiagonal                            |
*     |  2 | dchkbd type 2                      | Identity bidiagonal                        |
*     |  3 | dchkbd type 3  (DLATMS mode 4)     | Arithmetic-spectrum bidiag                 |
*     |  4 | dchkbd type 4  (DLATMS mode 3)     | Geometric-spectrum bidiag                  |
*     |  5 | dchkbd type 16                     | Log-distributed bidiag on [ulp, 1/ulp]     |
*     |  6 | dchkst type 8  (mode 4)            | Symm-tridiag arithmetic  -> Cholesky -> B  |
*     |  7 | dchkst type 9  (mode 3)            | Symm-tridiag geometric   -> Cholesky -> B  |
*     |  8 | dchkst type 10 (mode 1)            | Symm-tridiag clustered   -> Cholesky -> B  |
*     |  9 | dchkbd type 5                      | Arithmetic-spectrum bidiag x sqrt(overflow)|
*     | 10 | dchkbd type 6                      | Arithmetic-spectrum bidiag x sqrt(underflow)|
*     | 11 | dchkst type 3  (mode 4)            | Diagonal arithmetic                        |
*     | 12 | dchkst type 4  (mode 3)            | Diagonal geometric                         |
*     | 13 | dchkst type 5  (mode 1)            | Diagonal clustered (1, ulp, ..., ulp)      |
*     | 14 | dchkst type 6                      | Diagonal arithmetic x sqrt(overflow)       |
*     | 15 | dchkst type 7                      | Diagonal arithmetic x sqrt(underflow)      |
*     | 16 | dchkst type 21 (mode 3)            | SPD tridiag geometric, diag-dominant       |
*                                               |   -> Cholesky -> B                         |
*
*     Arguments
*     ---------
*     JTYPE (input)  INTEGER              matrix type 1..16
*     N     (input)  INTEGER              order (N >= 2 for interesting types)
*     ISEED (in/out) INTEGER(4)           random seed for DLATMS/DLARNV
*     BD    (output) DOUBLE(N)            bidiagonal main diagonal
*     BE    (output) DOUBLE(max(1,N-1))   bidiagonal off-diagonal
*     UPLO  (output) CHARACTER            'U' for upper (always here).
*     A     (workspace) DOUBLE(LDA,N)     scratch for DLATMS dense output
*     LDA   (input)  INTEGER              LDA >= N
*     D     (workspace) DOUBLE(N)         scratch for DLATMS prescribed spectrum
*     WORK  (workspace) DOUBLE(3*N)       scratch for DLATMS
*     INFO  (output) INTEGER              0 ok, +i => DPTTRF failed at pivot i,
*                                           -1 unknown JTYPE, -2 N too small
*
      SUBROUTINE DBSVRMG( JTYPE, N, ISEED, BD, BE, UPLO,
     $                    A, LDA, D, WORK, INFO )
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, JTYPE, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   A( LDA, * ), BD( * ), BE( * ), D( * ),
     $                   WORK( * )
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF, KTINY, SHIFT_FRAC
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                     HALF = 0.5D0, KTINY = 1.0D-6,
     $                     SHIFT_FRAC = 1.0D-3 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K, MODE
      DOUBLE PRECISION   ANORM, COND, DIAGMIN, DMAX, EPS, GL, LOGU,
     $                   OVFL, PVT, RNDVAL, SFMIN, SHIFT, TVAL, ULP
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RSEED( 4 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLATMS, DLASET, DLARNV, DPTTRF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, EXP, LOG, MAX, MIN, MOD, SQRT
*     ..
*     .. Executable Statements ..
      INFO = 0
      UPLO = 'U'
      EPS  = DLAMCH( 'Epsilon' )
      ULP  = DLAMCH( 'Precision' )
      SFMIN = DLAMCH( 'Safe Minimum' )
      OVFL = ONE / SFMIN
*
      IF( N.LT.2 ) THEN
         INFO = -2
         RETURN
      END IF
*
      DO 5 I = 1, N
         BD( I ) = ZERO
    5 CONTINUE
      DO 6 I = 1, MAX( 1, N-1 )
         BE( I ) = ZERO
    6 CONTINUE
*
*     ================= direct-bidiagonal types =========================
*
      IF( JTYPE.EQ.1 ) THEN
*        Type 1: zero bidiagonal.
         RETURN
*
      ELSE IF( JTYPE.EQ.2 ) THEN
*        Type 2: identity bidiagonal.
         DO 10 I = 1, N
            BD( I ) = ONE
   10    CONTINUE
         RETURN
*
      ELSE IF( JTYPE.EQ.3 .OR. JTYPE.EQ.4 .OR. JTYPE.EQ.5 .OR.
     $         JTYPE.EQ.9 .OR. JTYPE.EQ.10 ) THEN
*        Types 3, 4, 9, 10 : arithmetic/clustered spectrum via DLATMS,
*                            KL=0/KU=1 -> upper bidiagonal.
*        Type 5            : log-distributed bidiag on [ulp^2, ulp^-2].
         IF( JTYPE.EQ.5 ) THEN
*           Direct fill: BD(i) log-uniform on [ulp, 1/ulp].
            LOGU = LOG( ULP )
            DO 20 I = 1, N
               TVAL = TWO * DBLE( I - 1 ) / DBLE( N - 1 ) - ONE
               BD( I ) = EXP( TVAL * ( -LOGU ) )
   20       CONTINUE
            DO 21 I = 1, N - 1
               BE( I ) = KTINY * MIN( BD( I ), BD( I+1 ) )
   21       CONTINUE
            RETURN
         END IF
*
         COND = ONE / EPS
         DMAX = ONE
         IF( JTYPE.EQ.3 .OR. JTYPE.EQ.9 .OR. JTYPE.EQ.10 ) THEN
            MODE = 4
         ELSE
            MODE = 3
         END IF
         IF( JTYPE.EQ.9 ) DMAX = SQRT( OVFL ) * KTINY
         IF( JTYPE.EQ.10 ) DMAX = SQRT( SFMIN ) * ( ONE / KTINY )
*
         CALL DLATMS( N, N, 'S', ISEED, 'N', D, MODE, COND, DMAX,
     $                0, 1, 'N', A, LDA, WORK, INFO )
         IF( INFO.NE.0 ) THEN
            INFO = 100 + INFO
            RETURN
         END IF
         DO 30 I = 1, N
            BD( I ) = A( I, I )
   30    CONTINUE
         DO 31 I = 1, N - 1
            BE( I ) = A( I, I + 1 )
   31    CONTINUE
         RETURN
      END IF
*
*     ================= diagonal-spectrum types (dchkst 3-7) ============
*
      IF( JTYPE.EQ.11 .OR. JTYPE.EQ.12 .OR. JTYPE.EQ.13 .OR.
     $    JTYPE.EQ.14 .OR. JTYPE.EQ.15 ) THEN
*        Diagonal spectra via DLATMS(SYM='S', KL=KU=0):
*          11 = dchkst 3  (arith), 12 = dchkst 4  (geometric),
*          13 = dchkst 5  (clustered),
*          14 = dchkst 6  arith x sqrt(overflow),
*          15 = dchkst 7  arith x sqrt(underflow).
         IF( JTYPE.EQ.11 .OR. JTYPE.EQ.14 .OR. JTYPE.EQ.15 ) THEN
            MODE = 4
         ELSE IF( JTYPE.EQ.12 ) THEN
            MODE = 3
         ELSE
            MODE = 1
         END IF
         COND = ONE / EPS
         IF( JTYPE.EQ.14 ) THEN
            DMAX = SQRT( OVFL ) * ULP / DBLE( N )
         ELSE IF( JTYPE.EQ.15 ) THEN
            DMAX = SQRT( SFMIN ) * DBLE( N ) / ULP
         ELSE
            DMAX = ONE
         END IF
         CALL DLATMS( N, N, 'S', ISEED, 'S', D, MODE, COND, DMAX,
     $                0, 0, 'N', A, LDA, WORK, INFO )
         IF( INFO.NE.0 ) THEN
            INFO = 400 + INFO
            RETURN
         END IF
         DO 32 I = 1, N
            BD( I ) = A( I, I )
   32    CONTINUE
         RETURN
      END IF
*
*     ================= tridiagonal-plus-Cholesky types =================
*
      IF( JTYPE.EQ.6 .OR. JTYPE.EQ.7 .OR. JTYPE.EQ.8 .OR.
     $    JTYPE.EQ.16 ) THEN
*        Build symmetric tridiagonal via DLATMS(SYM='S' or 'P', KL=KU=1).
         IF( JTYPE.EQ.6 ) THEN
            MODE = 4
         ELSE IF( JTYPE.EQ.7 ) THEN
            MODE = 3
         ELSE IF( JTYPE.EQ.8 ) THEN
            MODE = 1
         ELSE
*           JTYPE = 16: dchkst 21 -- geometric SPD tridiag.
            MODE = 3
         END IF
         COND = ONE / EPS
         DMAX = ONE
         IF( JTYPE.EQ.16 ) THEN
            CALL DLATMS( N, N, 'S', ISEED, 'P', D, MODE, COND, DMAX,
     $                   1, 1, 'N', A, LDA, WORK, INFO )
         ELSE
            CALL DLATMS( N, N, 'S', ISEED, 'S', D, MODE, COND, DMAX,
     $                   1, 1, 'N', A, LDA, WORK, INFO )
         END IF
         IF( INFO.NE.0 ) THEN
            INFO = 200 + INFO
            RETURN
         END IF
*        For JTYPE=16 (dchkst 21) damp off-diagonals to keep the tridiag
*        diagonally dominant (mirror dchkst.f:879-887).
         IF( JTYPE.EQ.16 ) THEN
            DO 42 I = 2, N
               TVAL = ABS( A( I-1, I ) ) /
     $                SQRT( ABS( A( I-1, I-1 ) * A( I, I ) ) )
               IF( TVAL.GT.HALF ) THEN
                  A( I-1, I ) = HALF *
     $                  SQRT( ABS( A( I-1, I-1 ) * A( I, I ) ) )
                  A( I, I-1 ) = A( I-1, I )
               END IF
   42       CONTINUE
         END IF
         DO 40 I = 1, N
            BD( I ) = A( I, I )
   40    CONTINUE
         DO 41 I = 1, N - 1
            BE( I ) = A( I, I + 1 )
   41    CONTINUE
*
      ELSE
         INFO = -1
         RETURN
      END IF
*
*     -------- Convert symmetric tridiagonal T (in BD, BE) to bidiagonal
*     via LDL^T factorization (mirror test_dense_to_bidiag.py:104-138).
*
*     1) Gershgorin lower bound.
      GL = BD( 1 ) - ABS( BE( 1 ) )
      DO 70 I = 2, N - 1
         GL = MIN( GL, BD( I ) - ABS( BE( I-1 ) ) - ABS( BE( I ) ) )
   70 CONTINUE
      GL = MIN( GL, BD( N ) - ABS( BE( N-1 ) ) )
*
*     2) Shift so T - shift*I is safely positive definite.
      SHIFT = GL - SHIFT_FRAC * MAX( ABS( GL ), ONE )
      DO 71 I = 1, N
         BD( I ) = BD( I ) - SHIFT
   71 CONTINUE
*
*     3) Factor T_shift = L * D * L^T via DPTTRF (in-place: on exit BD
*        holds diag(D_chol) and BE holds subdiag(L)).
      CALL DPTTRF( N, BD, BE, INFO )
      IF( INFO.NE.0 ) THEN
*        A non-positive pivot appeared despite the shift.  Fall back to
*        a larger shift and retry once.
         DO 80 I = 1, N
            BD( I ) = BD( I ) - MAX( ABS( SHIFT ), ONE )
   80    CONTINUE
         CALL DPTTRF( N, BD, BE, INFO )
         IF( INFO.NE.0 ) THEN
            INFO = 300 + INFO
            RETURN
         END IF
      END IF
*
*     4) Extract bidiagonal factor B: BD(i) := sqrt(D_chol(i)),
*        BE(i) := E(i) * BD(i).
      DO 90 I = 1, N
         BD( I ) = SQRT( ABS( BD( I ) ) )
   90 CONTINUE
      DO 91 I = 1, N - 1
         BE( I ) = BE( I ) * BD( I )
   91 CONTINUE
*
      RETURN
*
*     End of DBSVRMG
*
      END
