      SUBROUTINE DLAR1V_TGK( REP, N, B1, BN, SIGMA, D, L, LD, LLD, EVAL,
     $                   GERSCH, Z, ZTZ, MINGMA, R, ISUPPZ, WORK )
*
*  -- New auxiliary routine for bidiagonal SVD via TGK-rooted MR^3 --
*     Based on LAPACK DLAR1V (Dhillon & Marques, Nov 11 2003).
*     stegr_ID/dlar1v.f is UNCHANGED; this is a separate copy that adds a
*     representation switch so the matrix M may be supplied either as an
*     L D L^T factorization (REP='L', identical to the original) or as a
*     symmetric tridiagonal -- in particular the Tridiagonal Golub-Kahan
*     (TGK) matrix -- with REP='T'.
*
*     .. Scalar Arguments ..
      CHARACTER          REP
      INTEGER            B1, BN, N, R
      DOUBLE PRECISION   EVAL, MINGMA, SIGMA, ZTZ
*     ..
*     .. Array Arguments ..
      INTEGER            ISUPPZ( * )
      DOUBLE PRECISION   D( * ), GERSCH( * ), L( * ), LD( * ), LLD( * ),
     $                   WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAR1V_TGK computes the (scaled) r-th column of the inverse of the
*  submatrix in rows B1 through BN of  M - sigma I,  where M is either
*  an L D L^T representation (REP='L') or a symmetric tridiagonal matrix
*  such as the TGK matrix (REP='T').  When sigma is close to an
*  eigenvalue, the computed vector is an accurate eigenvector.  The steps
*  are those of DLAR1V:
*  (a) stationary qd transform   M - sigma I = L(+) D(+) L(+)^T,
*  (b) progressive qd transform  M - sigma I = U(-) D(-) U(-)^T,
*  (c) the twist index r = argmin_k |gamma_k|, where gamma_k is the k-th
*      diagonal element of the twisted factorization and equals the
*      reciprocal of (M - sigma I)^{-1}_{kk},
*  (d) the scaled r-th column of the inverse from the twisted factor.
*
*  For REP='T' the qd transforms run on the tridiagonal entries directly
*  (diagonal in D, off-diagonals beta_i in L), so LD and LLD are not
*  referenced.  Using the identity  gamma_k = D(+)(k) + D(-)(k) - J_kk
*  with J_kk = D(k) - sigma (thesis eq. 3.1.12), the WORK slots are filled
*  so that the twist selection and back-substitution code below are
*  shared, unchanged, with the REP='L' path:
*       WORK(INDS + k-1) = D(+)(k)
*       WORK(INDP + k-1) = D(-)(k) - ( D(k) - sigma )
*       WORK(k)          = L(+)(k) = beta_k / D(+)(k)
*       WORK(INDUMN + k) = U(-)(k) = beta_k / D(-)(k+1)
*  so WORK(INDS+k-1) + WORK(INDP+k-1) = gamma_k.  For the TGK matrix
*  (D == 0) this is gamma_k = D(+)(k) + D(-)(k) + sigma.
*
*  Zero pivots are expected for the indefinite zero-diagonal TGK root.
*  In the direct recurrence a zero pivot self-heals through IEEE
*  arithmetic (D(+) -> -Inf, then L(+) -> 0), so no stationary/progressive
*  restart is needed; only the back-substitution uses a slow path (to
*  avoid Inf*0 = NaN), with beta-ratios per thesis Algorithm 3.3.1.
*
*  Arguments
*  =========
*
*  REP   (input) CHARACTER*1   'L' = L D L^T input; 'T' = tridiagonal/TGK.
*  N,B1,BN,SIGMA,EVAL,GERSCH,Z,ZTZ,MINGMA,R,ISUPPZ,WORK : as in DLAR1V.
*  D     (input) DOUBLE PRECISION (N)   REP='L': diag(D).  REP='T': diag(M).
*  L     (input) DOUBLE PRECISION (N-1) REP='L': subdiag(L). REP='T': beta_i.
*  LD,LLD(input) DOUBLE PRECISION (N-1) REP='L': L*D, L*L*D. REP='T': unused.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLKSIZ
      PARAMETER          ( BLKSIZ = 32 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            SAWNAN, TGK
      INTEGER            FROM, I, INDP, INDS, INDUMN, J, R1, R2, TO
      DOUBLE PRECISION   DMINUS, DPLUS, EPS, S, TMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      TGK = LSAME( REP, 'T' )
      EPS = DLAMCH( 'Precision' )
      IF( R.EQ.0 ) THEN
*
*        Eliminate the top and bottom indices from among the possible
*        values of R, using the Gerschgorin intervals.
*
         R1 = B1
         R2 = BN
         DO 10 I = B1, BN
            IF( EVAL.GE.GERSCH( 2*I-1 ) .AND. EVAL.LE.GERSCH( 2*I ) )
     $           THEN
               R1 = I
               GO TO 20
            END IF
   10    CONTINUE
         GO TO 40
   20    CONTINUE
         DO 30 I = BN, B1, -1
            IF( EVAL.GE.GERSCH( 2*I-1 ) .AND. EVAL.LE.GERSCH( 2*I ) )
     $           THEN
               R2 = I
               GO TO 40
            END IF
   30    CONTINUE
      ELSE
         R1 = R
         R2 = R
      END IF
*
   40 CONTINUE
      INDUMN = N
      INDS = 2*N + 1
      INDP = 3*N + 1
      SAWNAN = .FALSE.
*
*     ================================================================
*     Stationary transform (differential form) up to the index R2.
*     ================================================================
*
      IF( TGK ) THEN
*
*        Tridiagonal/TGK parent.  Direct recurrence; store D(+)(k).
*
         DPLUS = D( B1 ) - SIGMA
         WORK( INDS+B1-1 ) = DPLUS
         IF( DPLUS.EQ.ZERO ) SAWNAN = .TRUE.
         DO 45 I = B1, R2 - 1
            WORK( I ) = L( I ) / DPLUS
            DPLUS = D( I+1 ) - SIGMA - WORK( I )*L( I )
            WORK( INDS+I ) = DPLUS
            IF( DPLUS.EQ.ZERO ) SAWNAN = .TRUE.
   45    CONTINUE
      ELSE
*
*        L D L^T parent (original DLAR1V stationary loop).
*
         IF( B1.EQ.1 ) THEN
            WORK( INDS ) = ZERO
         ELSE
            WORK( INDS ) = LLD( B1-1 )
         END IF
         S = WORK( INDS ) - SIGMA
         DO 50 I = B1, R2 - 1
            DPLUS = D( I ) + S
            WORK( I ) = LD( I ) / DPLUS
            WORK( INDS+I ) = S*WORK( I )*L( I )
            S = WORK( INDS+I ) - SIGMA
   50    CONTINUE
*
         IF( .NOT.( S.GT.ZERO .OR. S.LT.ONE ) ) THEN
*           Slower version of the above loop if a NaN is detected.
            SAWNAN = .TRUE.
            J = B1 + 1
   60       CONTINUE
            IF( WORK( INDS+J ).GT.ZERO .OR. WORK( INDS+J ).LT.ONE ) THEN
               J = J + 1
               GO TO 60
            END IF
            WORK( INDS+J ) = LLD( J )
            S = WORK( INDS+J ) - SIGMA
            DO 70 I = J + 1, R2 - 1
               DPLUS = D( I ) + S
               WORK( I ) = LD( I ) / DPLUS
               IF( WORK( I ).EQ.ZERO ) THEN
                  WORK( INDS+I ) = LLD( I )
               ELSE
                  WORK( INDS+I ) = S*WORK( I )*L( I )
               END IF
               S = WORK( INDS+I ) - SIGMA
   70       CONTINUE
         END IF
      END IF
*
*     ================================================================
*     Progressive transform (differential form) down to the index R1.
*     ================================================================
*
      IF( TGK ) THEN
*
*        Tridiagonal/TGK parent.  Direct recurrence.  Store
*        WORK(INDP+i-1) = D(-)(i) - (D(i)-sigma)  and  U(-)(i).
*
         DMINUS = D( BN ) - SIGMA
         WORK( INDP+BN-1 ) = ZERO
         IF( DMINUS.EQ.ZERO ) SAWNAN = .TRUE.
         DO 75 I = BN - 1, R1, -1
            WORK( INDUMN+I ) = L( I ) / DMINUS
            DMINUS = D( I ) - SIGMA - WORK( INDUMN+I )*L( I )
            WORK( INDP+I-1 ) = DMINUS - ( D( I ) - SIGMA )
            IF( DMINUS.EQ.ZERO ) SAWNAN = .TRUE.
   75    CONTINUE
      ELSE
*
*        L D L^T parent (original DLAR1V progressive loop).
*
         WORK( INDP+BN-1 ) = D( BN ) - SIGMA
         DO 80 I = BN - 1, R1, -1
            DMINUS = LLD( I ) + WORK( INDP+I )
            TMP = D( I ) / DMINUS
            WORK( INDUMN+I ) = L( I )*TMP
            WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - SIGMA
   80    CONTINUE
         TMP = WORK( INDP+R1-1 )
         IF( .NOT.( TMP.GT.ZERO .OR. TMP.LT.ONE ) ) THEN
*           Slower version of the above loop if a NaN is detected.
            SAWNAN = .TRUE.
            J = BN - 3
   90       CONTINUE
            IF( WORK( INDP+J ).GT.ZERO .OR. WORK( INDP+J ).LT.ONE ) THEN
               J = J - 1
               GO TO 90
            END IF
            WORK( INDP+J ) = D( J+1 ) - SIGMA
            DO 100 I = J, R1, -1
               DMINUS = LLD( I ) + WORK( INDP+I )
               TMP = D( I ) / DMINUS
               WORK( INDUMN+I ) = L( I )*TMP
               IF( TMP.EQ.ZERO ) THEN
                  WORK( INDP+I-1 ) = D( I ) - SIGMA
               ELSE
                  WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - SIGMA
               END IF
  100       CONTINUE
         END IF
      END IF
*
*     ================================================================
*     Twist index: r = argmin_k |gamma_k|, gamma_k stored as
*     WORK(INDS+k-1) + WORK(INDP+k-1) for both representations.
*     ================================================================
*
      MINGMA = WORK( INDS+R1-1 ) + WORK( INDP+R1-1 )
      IF( MINGMA.EQ.ZERO )
     $   MINGMA = EPS*WORK( INDS+R1-1 )
      R = R1
      DO 110 I = R1, R2 - 1
         TMP = WORK( INDS+I ) + WORK( INDP+I )
         IF( TMP.EQ.ZERO )
     $      TMP = EPS*WORK( INDS+I )
         IF( ABS( TMP ).LT.ABS( MINGMA ) ) THEN
            MINGMA = TMP
            R = I + 1
         END IF
  110 CONTINUE
*
*     ================================================================
*     Compute the (scaled) r-th column of the inverse.  WORK(k)=L(+)(k)
*     and WORK(INDUMN+k)=U(-)(k) for both representations, so the fast
*     path is shared.  The slow (SAWNAN) path uses LD-ratios for REP='L'
*     and beta-ratios (L) for REP='T'.
*     ================================================================
*
      ISUPPZ( 1 ) = B1
      ISUPPZ( 2 ) = BN
      Z( R ) = ONE
      ZTZ = ONE
      IF( .NOT.SAWNAN ) THEN
*        No support truncation.  A TGK eigenvector can be BIMODAL -- the
*        left singular vector u occupies one end of z and the right vector
*        v the other, with a near-zero trough between -- so the usual
*        "two consecutive tiny entries => the eigenvector has died out,
*        stop" early-out would wrongly amputate the far hump.  This is
*        true for EVERY vector this routine computes, whether the
*        representation is the root TGK (REP='T') or a deeper L D L^T
*        child of it (REP='L'), since the eigenvectors are TGK
*        eigenvectors in both cases.  Hence we always form the full
*        vector and leave ISUPPZ = [B1,BN].
         DO 130 I = R - 1, B1, -1
            Z( I ) = -( WORK( I )*Z( I+1 ) )
            ZTZ = ZTZ + Z( I )*Z( I )
  130    CONTINUE
         DO 150 I = R + 1, BN
            Z( I ) = -( WORK( INDUMN+I-1 )*Z( I-1 ) )
            ZTZ = ZTZ + Z( I )*Z( I )
  150    CONTINUE
      ELSE
*
*        Slow path: some pivot vanished, so L(+)/U(-) may be +-Inf.
*        Use the recurrence of thesis Algorithm 3.3.1, switching the
*        coupling ratio between LD (REP='L') and beta=L (REP='T').
*
         DO 160 I = R - 1, B1, -1
            IF( Z( I+1 ).EQ.ZERO ) THEN
               IF( TGK ) THEN
                  Z( I ) = -( L( I+1 ) / L( I ) )*Z( I+2 )
               ELSE
                  Z( I ) = -( LD( I+1 ) / LD( I ) )*Z( I+2 )
               END IF
            ELSE
               Z( I ) = -( WORK( I )*Z( I+1 ) )
            END IF
            ZTZ = ZTZ + Z( I )*Z( I )
  160    CONTINUE
  170    CONTINUE
         DO 180 I = R, BN - 1
            IF( Z( I ).EQ.ZERO ) THEN
               IF( TGK ) THEN
                  Z( I+1 ) = -( L( I-1 ) / L( I ) )*Z( I-1 )
               ELSE
                  Z( I+1 ) = -( LD( I-1 ) / LD( I ) )*Z( I-1 )
               END IF
            ELSE
               Z( I+1 ) = -( WORK( INDUMN+I )*Z( I ) )
            END IF
            ZTZ = ZTZ + Z( I+1 )*Z( I+1 )
  180    CONTINUE
  190    CONTINUE
      END IF
*
      RETURN
*
*     End of DLAR1V_TGK
*
      END
