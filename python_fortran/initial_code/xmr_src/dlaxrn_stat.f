      SUBROUTINE DLAXRN_STAT(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, NCOUNT, S
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  NCOUNT
      DOUBLE PRECISION, INTENT(OUT)  ::  S
*
*  Purpose
*  =======
*
*
*  ======================================================================
*
*     .. Parameters ..
*

*
*     PRBRK
*         Determines branch to choose for breaking a block ("R_brk").
*         Default: 8
*
      INTEGER, PARAMETER  ::  PRBRK = 8

*
*     .. Declarations ..
*
      DOUBLE PRECISION ZERO, ONE
      PARAMETER( ZERO=0.0D0, ONE=1.0D0 )
*
*     .. Locals ..
*
      DOUBLE PRECISION  SBOUND, NEGPIV
      DOUBLE PRECISION  SMT, GP, OLDGP, OLDSMT, X
      INTEGER           I, IPREV
*
*  ===== Executable Statements ==========================================
*
      SBOUND = MIN( PRBRK, N ) * ABS(TAU)
      NEGPIV = -PIVMIN
      IF( DIR.EQ.1 )THEN
         I    = 1
      ELSE
         I    = N
      ENDIF

      S = ZERO
      NCOUNT = 0
      DO
         DO
            IF( I.EQ.LBBEGK(IBB) )THEN
               EXIT
            ENDIF
            SMT = S - TAU
            GP = G(I) + SMT
            IF( ABS(GP).LT.PIVMIN )  GP = NEGPIV
            IF( GP.LT.ZERO )  NCOUNT = NCOUNT+1
            S = NGN(I) * (SMT / GP)
            I = I+DIR
         ENDDO
         IF( I.EQ.K )THEN
            EXIT
         ENDIF

         IBB = IBB+DIR

         SMT = S - TAU
         GP = G(I) + SMT
         IF( ABS(GP).LT.PIVMIN )  GP = NEGPIV
         IF( GP.LT.ZERO )  NCOUNT = NCOUNT+1
         S = - GNSQ(I) / GP
         OLDGP = GP
         OLDSMT = SMT
         IPREV = I
         I = I+DIR
         IF( I.EQ.K )THEN
            EXIT
         ENDIF

         SMT = S - TAU
         GP = G(I) + SMT
         IF( ABS(GP).LT.PIVMIN )  GP = NEGPIV
         IF( GP.LT.ZERO )  NCOUNT = NCOUNT+1
*
*        Compute next s after a block was broken
*        Notes:
*        (1) Recall that NGN(I-DIR) is undefined
*        (2) !!! SIGN(A,B) = Abs(A) if B>=0, -Abs(A) ow !!!
*
         IF( G(IPREV).EQ.ZERO )THEN
            S = - GNSQ(I) / GP
         ELSE
            IF( (ABS(S).LE.SBOUND) .OR.
     $          (SIGN(ONE,G(IPREV)).NE.SIGN(ONE,OLDGP)) )
     $      THEN
               X = G(IPREV)*SMT + GNSQ(IPREV)
            ELSE
               X = -S*OLDSMT - G(IPREV)*TAU
            ENDIF
            S = (GNSQ(I)*X) / (BDET(IPREV)*GP)
         ENDIF
         I = I+DIR
      ENDDO

      END SUBROUTINE DLAXRN_STAT
*
************************************************************************
