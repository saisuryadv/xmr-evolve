      PROGRAM TEST_DBDSQR_FULL
*
*     Full-SVD DBDSQR runner (singular values + U + V), used as the third
*     solver in eval_dbdsvdmr3.py (the n^3 baseline against which DBDSVDMR3
*     and self-contained MR3-GK are compared).
*
*     Reads an upper-bidiagonal matrix from an STCollection-format .dat
*     file (line 1 = N, lines 2..N+1 = IDX D(i) E(i)), computes B = U S V^T
*     via DBDSQR with NCVT = NRU = N (start from identity for both U^T and
*     VT, return U and VT), and prints the same paper-norm metric block as
*     test_stcoll_alloc.f so the orchestrator can parse one log format
*     across all three solvers.
*
*     Benchmark protocol (matches dbdsqr_ref.f and test_stcoll_alloc.f):
*       - Save D, E; DBDSQR overwrites both.
*       - One untimed warmup call.
*       - Adaptive bench loop: stop on (NREPS=1 && t>=0.5s) or
*         (NREPS>=3 && total>=0.2s) or NREPS>=30.
*       - Cache flush (32 MiB write+read) before every timed call.
*       - SYSTEM_CLOCK with INTEGER*8 ticks.  Report MIN-of-N.
*
      IMPLICIT NONE
      INTEGER, PARAMETER :: JUNK_SIZE = 4*1024*1024
      INTEGER N, INFO, INFO2, I, J, K, IDX, IOS, NREPS_EVAL, NREPS_REF
      INTEGER*8 TC1, TC2, TCR
      DOUBLE PRECISION DD, EE, TMP, SMAX, SMIN, NEPS, RELSV
      DOUBLE PRECISION RESID, ORTHU, ORTHV
      DOUBLE PRECISION RESID_P, ORTHU_P, ORTHV_P, RNORM1, RNORM2
      DOUBLE PRECISION T_EVAL, T_DBDSQR, T_THIS, T_TOTAL, JSUM
      CHARACTER*4096 FNAME
      DOUBLE PRECISION, ALLOCATABLE :: D(:), E(:), DSAVE(:), ESAVE(:)
      DOUBLE PRECISION, ALLOCATABLE :: ESCRATCH(:), S(:), SREF(:)
      DOUBLE PRECISION, ALLOCATABLE :: U(:,:), VT(:,:), WORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: AMAT(:,:), EVAL(:), RB(:)
      DOUBLE PRECISION, ALLOCATABLE :: JUNK(:)
      DOUBLE PRECISION DUMMY(1,1)
      INTEGER LDU, LDVT, LWORK
      DOUBLE PRECISION DLAMCH, DNRM2
      EXTERNAL DBDSQR, DSYEV, DLAMCH, DNRM2
      INTRINSIC ABS, MAX, MIN, DBLE
*
      CALL GET_COMMAND_ARGUMENT( 1, FNAME )
      OPEN( UNIT=10, FILE=FNAME, STATUS='OLD', IOSTAT=IOS )
      IF( IOS.NE.0 ) THEN
         WRITE(*,*) 'cannot open ', TRIM(FNAME)
         STOP 1
      END IF
      READ(10,*) N
      IF( N.LT.1 ) STOP 1
*
      LDU   = N
      LDVT  = N
      LWORK = MAX( 4*N, 100 )
      ALLOCATE( D(N), E(N), DSAVE(N), ESAVE(N), ESCRATCH(N) )
      ALLOCATE( S(N), SREF(N) )
      ALLOCATE( U(LDU,N), VT(LDVT,N), WORK(LWORK) )
      ALLOCATE( AMAT(N,N), EVAL(N), RB(N), JUNK(JUNK_SIZE) )
*
      DO I = 1, N
         READ(10,*) IDX, DD, EE
         DSAVE(I) = DD
         ESAVE(I) = EE
      END DO
      CLOSE(10)
      ESAVE(N) = 0.0D0
*
*     Initialize cache-flush scratch buffer.
      DO I = 1, JUNK_SIZE
         JUNK(I) = DBLE(I)
      END DO
*
*     -------- DBDSQR singular-values-only reference (same as test_stcoll_alloc) ----
      DO I = 1, N
         SREF(I) = DSAVE(I)
         ESCRATCH(I) = ESAVE(I)
      END DO
      CALL DBDSQR( 'U', N, 0, 0, 0, SREF, ESCRATCH, DUMMY, 1, DUMMY, 1,
     $             DUMMY, 1, WORK, INFO )
      NREPS_REF = 0
      T_DBDSQR  = 1.0D30
      T_TOTAL   = 0.0D0
      JSUM      = 0.0D0
   30 CONTINUE
         DO I = 1, N
            SREF(I) = DSAVE(I)
            ESCRATCH(I) = ESAVE(I)
         END DO
         DO J = 1, JUNK_SIZE
            JUNK(J) = JUNK(J) * 1.0000001D0 + 1.0D-30
            JSUM = JSUM + JUNK(J)
         END DO
         CALL SYSTEM_CLOCK( TC1, TCR )
         CALL DBDSQR( 'U', N, 0, 0, 0, SREF, ESCRATCH, DUMMY, 1,
     $                DUMMY, 1, DUMMY, 1, WORK, INFO )
         CALL SYSTEM_CLOCK( TC2 )
         T_THIS = DBLE( TC2 - TC1 ) / DBLE( TCR )
         IF( T_THIS.LT.T_DBDSQR ) T_DBDSQR = T_THIS
         T_TOTAL = T_TOTAL + T_THIS
         NREPS_REF = NREPS_REF + 1
         IF( NREPS_REF.EQ.1 .AND. T_THIS.GE.0.5D0 ) GOTO 31
         IF( NREPS_REF.GE.3 .AND. T_TOTAL.GE.0.2D0 ) GOTO 31
         IF( NREPS_REF.GE.30 ) GOTO 31
         GOTO 30
   31 CONTINUE
*
*     -------- Full DBDSQR with NCVT = NRU = N (the n^3 solver under eval) ----
*     Initialize VT = I (DBDSQR will overwrite with VT) and U = I (with U).
*     NRU here means "we pass an NRU x N matrix U into which DBDSQR will
*     left-multiply the right singular vectors, returning U^orig * Q_left".
*     Starting from U=I we recover Q_left, i.e. the left singular vectors.
*     Warmup.
      DO I = 1, N
         D(I) = DSAVE(I)
         E(I) = ESAVE(I)
         DO J = 1, N
            VT(I,J) = 0.0D0
            U(I,J) = 0.0D0
         END DO
         VT(I,I) = 1.0D0
         U(I,I) = 1.0D0
      END DO
      CALL DBDSQR( 'U', N, N, N, 0, D, E, VT, LDVT, U, LDU,
     $             DUMMY, 1, WORK, INFO )
*
      NREPS_EVAL = 0
      T_EVAL    = 1.0D30
      T_TOTAL   = 0.0D0
   40 CONTINUE
         DO I = 1, N
            D(I) = DSAVE(I)
            E(I) = ESAVE(I)
            DO J = 1, N
               VT(I,J) = 0.0D0
               U(I,J) = 0.0D0
            END DO
            VT(I,I) = 1.0D0
            U(I,I) = 1.0D0
         END DO
         DO J = 1, JUNK_SIZE
            JUNK(J) = JUNK(J) * 1.0000001D0 + 1.0D-30
            JSUM = JSUM + JUNK(J)
         END DO
         CALL SYSTEM_CLOCK( TC1, TCR )
         CALL DBDSQR( 'U', N, N, N, 0, D, E, VT, LDVT, U, LDU,
     $                DUMMY, 1, WORK, INFO )
         CALL SYSTEM_CLOCK( TC2 )
         T_THIS = DBLE( TC2 - TC1 ) / DBLE( TCR )
         IF( T_THIS.LT.T_EVAL ) T_EVAL = T_THIS
         T_TOTAL = T_TOTAL + T_THIS
         NREPS_EVAL = NREPS_EVAL + 1
         IF( NREPS_EVAL.EQ.1 .AND. T_THIS.GE.0.5D0 ) GOTO 41
         IF( NREPS_EVAL.GE.3 .AND. T_TOTAL.GE.0.2D0 ) GOTO 41
         IF( NREPS_EVAL.GE.30 ) GOTO 41
         GOTO 40
   41 CONTINUE
      IF( JSUM.EQ.1.2345D-300 ) WRITE(*,*) 'never', JUNK(1)
*
*     DBDSQR returns S = D in descending order, U has columns = u_i,
*     VT has rows = v_i^T (since we started from VT=I and DBDSQR applies
*     V^T from the LEFT to VT).  Store DBDSQR singular values into S
*     descending; for comparison, RELSV against SREF (which is identical
*     to D after this call) is trivially zero -- we keep it for completeness.
      DO I = 1, N
         S(I) = D(I)
      END DO
      RELSV = 0.0D0
      DO I = 1, N
         IF( SREF(I).GT.0.0D0 )
     $      RELSV = MAX( RELSV, ABS(S(I)-SREF(I))/SREF(I) )
      END DO
      SMAX = 0.0D0
      SMIN = 1.0D30
      DO I = 1, N
         SMAX = MAX( SMAX, S(I) )
         IF( S(I).GT.0.0D0 ) SMIN = MIN( SMIN, S(I) )
      END DO
*
*     ---- max-norm metrics ----
*     Note: DBDSQR's S is descending; we reconstruct U * diag(S) * VT.
      RESID = 0.0D0
      DO I = 1, N
         DO J = 1, N
            TMP = 0.0D0
            DO K = 1, N
               TMP = TMP + U(I,K)*S(K)*VT(K,J)
            END DO
            DD = 0.0D0
            IF( I.EQ.J )      DD = DSAVE(I)
            IF( J.EQ.I+1 )    DD = ESAVE(I)
            RESID = MAX( RESID, ABS(TMP-DD) )
         END DO
      END DO
      ORTHU = 0.0D0
      ORTHV = 0.0D0
      DO I = 1, N
         DO J = 1, N
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
*     ---- paper-norm metrics ----
      RESID_P = 0.0D0
      DO 220 J = 1, N
         DO 210 I = 1, N
            TMP = DSAVE(I)*VT(J,I)
            IF( I.LT.N ) TMP = TMP + ESAVE(I)*VT(J,I+1)
            RB(I) = TMP - S(J)*U(I,J)
  210    CONTINUE
         RNORM1 = DNRM2( N, RB, 1 )
         DO 215 I = 1, N
            TMP = DSAVE(I)*U(I,J)
            IF( I.GE.2 ) TMP = TMP + ESAVE(I-1)*U(I-1,J)
            RB(I) = TMP - S(J)*VT(J,I)
  215    CONTINUE
         RNORM2 = DNRM2( N, RB, 1 )
         RESID_P = MAX( RESID_P, MAX( RNORM1, RNORM2 ) )
  220 CONTINUE
      IF( SMAX.GT.0.0D0 ) RESID_P = RESID_P / SMAX
*     spectral ||U^T U - I||_2
      DO 235 J = 1, N
         DO 230 I = 1, N
            TMP = 0.0D0
            DO 225 K = 1, N
               TMP = TMP + U(K,I)*U(K,J)
  225       CONTINUE
            IF( I.EQ.J ) TMP = TMP - 1.0D0
            AMAT(I,J) = TMP
  230    CONTINUE
  235 CONTINUE
      CALL DSYEV( 'N','U', N, AMAT, N, EVAL, WORK, LWORK, INFO2 )
      ORTHU_P = 0.0D0
      IF( INFO2.EQ.0 ) THEN
         DO I = 1, N
            ORTHU_P = MAX( ORTHU_P, ABS(EVAL(I)) )
         END DO
      ELSE
         ORTHU_P = ORTHU
      END IF
      DO 245 J = 1, N
         DO 240 I = 1, N
            TMP = 0.0D0
            DO 238 K = 1, N
               TMP = TMP + VT(I,K)*VT(J,K)
  238       CONTINUE
            IF( I.EQ.J ) TMP = TMP - 1.0D0
            AMAT(I,J) = TMP
  240    CONTINUE
  245 CONTINUE
      CALL DSYEV( 'N','U', N, AMAT, N, EVAL, WORK, LWORK, INFO2 )
      ORTHV_P = 0.0D0
      IF( INFO2.EQ.0 ) THEN
         DO I = 1, N
            ORTHV_P = MAX( ORTHV_P, ABS(EVAL(I)) )
         END DO
      ELSE
         ORTHV_P = ORTHV
      END IF
*
      WRITE(*,900) N, N, INFO, SMIN, SMAX, 0.0D0, RELSV, RESID,
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
      WRITE(*,940) T_EVAL, T_DBDSQR, RELSV, NREPS_EVAL, NREPS_REF
  940 FORMAT('   t_eval=',ES12.5,'   t_dbdsqr=',ES12.5,
     $       '   sv_drift=',ES12.5,
     $       '   nreps_eval=',I0,'   nreps_ref=',I0)
      DEALLOCATE( D, E, DSAVE, ESAVE, ESCRATCH, S, SREF )
      DEALLOCATE( U, VT, WORK, AMAT, EVAL, RB, JUNK )
      STOP
      END
