      PROGRAM TEST_STCOLL_ALLOC
*
*     Allocatable variant of dev/test_stcoll.f -- handles n up to ~10^4.
*     Reads an (upper) bidiagonal matrix from an STCollection-format .dat
*     file (line 1 = N, lines 2..N+1 = IDX  D(i)  E(i)), computes its SVD
*     via DBDSVDMR3, compares against DGESVD, and prints:
*       - max-norm metrics (existing behaviour):
*           rel.resid   = max|U S V^T - B|_max / sigma_max
*           orthU/orthV = max|U^T U - I|_max,  max|V^T V - I|_max
*       - Willems-Lang 2012 Table 5.1 paper-norm metrics:
*           paper.rel.resid = max_i max(||B v_i - u_i sigma_i||_2,
*                                       ||B^T u_i - v_i sigma_i||_2) / ||B||_2
*           paper.orthU/V  = ||U^T U - I||_2,  ||V^T V - I||_2   (spectral)
*     All values are also reported in units of N*eps.
*
      IMPLICIT NONE
      INTEGER N, M, INFO, INFO2, I, J, K, IDX, IOS
      DOUBLE PRECISION DD, EE, TMP, SDIFF, SMAX, SMIN, RELSV
      DOUBLE PRECISION RESID, ORTHU, ORTHV, NEPS
      DOUBLE PRECISION RESID_P, ORTHU_P, ORTHV_P
      DOUBLE PRECISION RNORM1, RNORM2
      CHARACTER*4096 FNAME
      INTEGER, ALLOCATABLE :: IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: D(:), E(:), S(:), SREF(:)
      DOUBLE PRECISION, ALLOCATABLE :: U(:,:), VT(:,:), WORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: AMAT(:,:), EVAL(:), RB(:)
      INTEGER LDU, LDVT, LWORK, LIWORK
      DOUBLE PRECISION DLAMCH, DNRM2
      EXTERNAL DBDSVDMR3, DGESVD, DSYEV, DLAMCH, DNRM2
      INTRINSIC ABS, MAX, MIN, DBLE, SQRT
*
      CALL GET_COMMAND_ARGUMENT( 1, FNAME )
      OPEN( UNIT=10, FILE=FNAME, STATUS='OLD', IOSTAT=IOS )
      IF( IOS.NE.0 ) THEN
         WRITE(*,*) 'cannot open ', TRIM(FNAME)
         STOP 1
      END IF
      READ(10,*) N
      IF( N.LT.1 ) THEN
         WRITE(*,*) 'bad N=', N
         STOP 1
      END IF
*
      LDU   = N
      LDVT  = N
      LWORK = 2*N*N + 100*N
      LIWORK= 30*N
      ALLOCATE( D(N), E(N), S(N), SREF(N) )
      ALLOCATE( U(LDU,N), VT(LDVT,N) )
      ALLOCATE( WORK(LWORK), IWORK(LIWORK) )
      ALLOCATE( AMAT(N,N), EVAL(N), RB(N) )
*
      DO I = 1, N
         READ(10,*) IDX, DD, EE
         D(I) = DD
         E(I) = EE
      END DO
      CLOSE(10)
      E(N) = 0.0D0
*
*     Build dense B in U as scratch for DGESVD reference (DGESVD destroys input).
      DO J = 1, N
         DO I = 1, N
            U(I,J) = 0.0D0
         END DO
      END DO
      DO I = 1, N
         U(I,I) = D(I)
      END DO
      DO I = 1, N-1
         U(I,I+1) = E(I)
      END DO
      CALL DGESVD( 'N','N', N, N, U, LDU, SREF, VT, LDVT, VT, LDVT,
     $             WORK, LWORK, INFO )
      IF( INFO.NE.0 ) WRITE(*,*) '   (DGESVD INFO=', INFO, ')'
*
      CALL DBDSVDMR3( 'V','U', N, D, E, S, U, LDU, VT, LDVT, M,
     $                WORK, LWORK, IWORK, LIWORK, INFO )
*
*     NaN/Inf audit
      K = 0
      DO 71 J = 1, M
         IF( .NOT.( S(J).EQ.S(J) .AND. ABS(S(J)).LE.1.0D300 ) ) K = K+1
         DO 70 I = 1, N
            IF( .NOT.( U(I,J).EQ.U(I,J) .AND.
     $                ABS(U(I,J)).LE.1.0D300 ) ) K = K + 1
            IF( .NOT.( VT(J,I).EQ.VT(J,I) .AND.
     $                ABS(VT(J,I)).LE.1.0D300 ) ) K = K + 1
   70    CONTINUE
   71 CONTINUE
      IF( K.GT.0 ) WRITE(*,*) '   *** NON-FINITE OUTPUT: ', K,
     $   ' bad entries in S/U/V ***'
*
      SMAX = 0.0D0
      SMIN = 1.0D30
      DO I = 1, M
         SMAX = MAX( SMAX, S(I) )
         IF( S(I).GT.0.0D0 ) SMIN = MIN( SMIN, S(I) )
      END DO
*
      SDIFF = 0.0D0
      RELSV = 0.0D0
      DO I = 1, M
         SDIFF = MAX( SDIFF, ABS(S(I)-SREF(M-I+1)) )
         IF( SREF(M-I+1).GT.0.0D0 )
     $      RELSV = MAX( RELSV, ABS(S(I)-SREF(M-I+1))/SREF(M-I+1) )
      END DO
*
*     ---- max-norm metrics (existing) ----
      RESID = 0.0D0
      DO I = 1, N
         DO J = 1, N
            TMP = 0.0D0
            DO K = 1, M
               TMP = TMP + U(I,K)*S(K)*VT(K,J)
            END DO
            DD = 0.0D0
            IF( I.EQ.J )      DD = D(I)
            IF( J.EQ.I+1 )    DD = E(I)
            RESID = MAX( RESID, ABS(TMP-DD) )
         END DO
      END DO
*
      ORTHU = 0.0D0
      ORTHV = 0.0D0
      DO I = 1, M
         DO J = 1, M
            TMP = 0.0D0
            DO K = 1, N
               TMP = TMP + U(K,I)*U(K,J)
            END DO
            IF( I.EQ.J ) TMP = TMP - 1.0D0
            ORTHU = MAX( ORTHU, ABS(TMP) )
            TMP = 0.0D0
            DO K = 1, N
               TMP = TMP + VT(I,K)*VT(J,K)
            END DO
            IF( I.EQ.J ) TMP = TMP - 1.0D0
            ORTHV = MAX( ORTHV, ABS(TMP) )
         END DO
      END DO
      IF( SMAX.GT.0.0D0 ) RESID = RESID / SMAX
*
*     ---- paper-norm metrics (Willems-Lang 2012 Table 5.1) ----
*
*     Per-triplet residual: max_i max(||B v_i - u_i s_i||_2,
*                                    ||B^T u_i - v_i s_i||_2) / ||B||_2
      RESID_P = 0.0D0
      DO 220 J = 1, M
*        r1 = B*v_j - s_j*u_j  ; v_j = VT(J,:),  u_j = U(:,J)
         DO 210 I = 1, N
            TMP = D(I)*VT(J,I)
            IF( I.LT.N ) TMP = TMP + E(I)*VT(J,I+1)
            RB(I) = TMP - S(J)*U(I,J)
  210    CONTINUE
         RNORM1 = DNRM2( N, RB, 1 )
*        r2 = B^T*u_j - s_j*v_j
         DO 215 I = 1, N
            TMP = D(I)*U(I,J)
            IF( I.GE.2 ) TMP = TMP + E(I-1)*U(I-1,J)
            RB(I) = TMP - S(J)*VT(J,I)
  215    CONTINUE
         RNORM2 = DNRM2( N, RB, 1 )
         RESID_P = MAX( RESID_P, MAX( RNORM1, RNORM2 ) )
  220 CONTINUE
      IF( SMAX.GT.0.0D0 ) RESID_P = RESID_P / SMAX
*
*     ortU_paper = ||U^T U - I||_2 = max |eigenvalue|  (DSYEV on M x M Gram)
      DO 235 J = 1, M
         DO 230 I = 1, M
            TMP = 0.0D0
            DO 225 K = 1, N
               TMP = TMP + U(K,I)*U(K,J)
  225       CONTINUE
            IF( I.EQ.J ) TMP = TMP - 1.0D0
            AMAT(I,J) = TMP
  230    CONTINUE
  235 CONTINUE
      CALL DSYEV( 'N','U', M, AMAT, N, EVAL, WORK, LWORK, INFO2 )
      ORTHU_P = 0.0D0
      IF( INFO2.EQ.0 ) THEN
         DO I = 1, M
            ORTHU_P = MAX( ORTHU_P, ABS(EVAL(I)) )
         END DO
      ELSE
         WRITE(*,*) '   (DSYEV U INFO=', INFO2, ')'
         ORTHU_P = ORTHU
      END IF
*
*     ortV_paper = ||V^T V - I||_2 ; VT stores V^T so V^T V = VT * VT^T
      DO 245 J = 1, M
         DO 240 I = 1, M
            TMP = 0.0D0
            DO 238 K = 1, N
               TMP = TMP + VT(I,K)*VT(J,K)
  238       CONTINUE
            IF( I.EQ.J ) TMP = TMP - 1.0D0
            AMAT(I,J) = TMP
  240    CONTINUE
  245 CONTINUE
      CALL DSYEV( 'N','U', M, AMAT, N, EVAL, WORK, LWORK, INFO2 )
      ORTHV_P = 0.0D0
      IF( INFO2.EQ.0 ) THEN
         DO I = 1, M
            ORTHV_P = MAX( ORTHV_P, ABS(EVAL(I)) )
         END DO
      ELSE
         WRITE(*,*) '   (DSYEV V INFO=', INFO2, ')'
         ORTHV_P = ORTHV
      END IF
*
*     ---- print ----
      WRITE(*,900) N, M, INFO, SMIN, SMAX, SDIFF, RELSV, RESID,
     $             ORTHU, ORTHV
  900 FORMAT(' N=',I6,' M=',I6,' INFO=',I3,' smin=',E10.3,
     $       ' smax=',E10.3,/,'   absdiff=',E10.3,' reldiff=',E10.3,
     $       ' rel.resid=',E10.3,' orthU=',E10.3,' orthV=',E10.3)
      NEPS = DBLE( N )*DLAMCH( 'Precision' )
      WRITE(*,910) RESID/NEPS, ORTHU/NEPS, ORTHV/NEPS
  910 FORMAT('   resid/(n.eps.||B||)=',F11.1,
     $       '   orthU/(n.eps)=',F11.1,'   orthV/(n.eps)=',F11.1)
      WRITE(*,920) RESID_P, ORTHU_P, ORTHV_P
  920 FORMAT('   paper.rel.resid=',E10.3,
     $       '   paper.orthU=',E10.3,'   paper.orthV=',E10.3)
      WRITE(*,930) RESID_P/NEPS, ORTHU_P/NEPS, ORTHV_P/NEPS
  930 FORMAT('   paper.resid/(n.eps.||B||)=',F11.1,
     $       '   paper.orthU/(n.eps)=',F11.1,
     $       '   paper.orthV/(n.eps)=',F11.1)
      DEALLOCATE( D, E, S, SREF, U, VT, WORK, IWORK )
      DEALLOCATE( AMAT, EVAL, RB )
      STOP
      END
