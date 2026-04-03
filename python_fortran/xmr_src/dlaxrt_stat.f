      SUBROUTINE DLAXRT_STAT(
     $             N, J, I0, I1, G, OMEGA, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, GPLUS, S
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, J, I0, I1
      INTEGER,          INTENT(IN)  ::  OMEGA(0:N+1)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      DOUBLE PRECISION, INTENT(OUT)  ::  GPLUS(I0:I1), S(N)
*
*  Purpose
*  =======
*
*     Computes stationary factorization from 1 to J or from N down to J.
*
*     I0=1 and I1=N-1
*       assumes that the twist of the source is >= J
*       produces
*         GPLUS(1),...,GPLUS(J-1) and S(1),...,S(J),
*       --> only S(1)=0 if J=1
*
*     I0=2 and I1=N
*       assumes that the twist of the source is <= J
*       produces
*         GPLUS(J+1),...,GPLUS(n) and S(J),...,S(N),
*       --> only S(N)=0 if J=N
*
*  ======================================================================
*
*     .. Declarations ..
*
      DOUBLE PRECISION ZERO, ONE
      PARAMETER( ZERO=0.0D0, ONE=1.0D0 )
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
*     .. Locals ..
*
      INTEGER I, IPREV, INEXT, DIR
      DOUBLE PRECISION SBOUND
      DOUBLE PRECISION SMT, OLDSMT, X

*
*     -- Executable Statements -----------------------------------------
*

      IF( I0.EQ.1 .AND. I1.EQ.N-1 )THEN
         DIR = +1
         I   = 1
      ELSEIF( I0.EQ.2 .AND. I1.EQ.N )THEN
         DIR = -1
         I   = N
      ELSE
         DIR = 0
         I   = -1
      ENDIF
      INEXT = I+DIR

      SBOUND = MIN( PRBRK, N ) * ABS(TAU)
      S(I) = ZERO
      DO
         DO
            IF( (I.EQ.J) .OR. (OMEGA(INEXT).NE.0) )THEN
               EXIT
            ENDIF
            SMT = S(I) - TAU
            GPLUS(I) = G(I) + SMT
            IF( ABS(GPLUS(I)).LT.PIVMIN )  GPLUS(I) = -PIVMIN
            S(INEXT) = NGN(I) * (SMT / GPLUS(I))
            I = INEXT
            INEXT = INEXT+DIR
         ENDDO
         IF( I.EQ.J )THEN
            EXIT
         ENDIF

         SMT = S(I) - TAU
         GPLUS(I) = G(I) + SMT
         IF( ABS(GPLUS(I)).LT.PIVMIN )  GPLUS(I) = -PIVMIN
         S(INEXT) = - GNSQ(I) / GPLUS(I)
         IF( INEXT.EQ.J )THEN
            EXIT
         ENDIF
         IPREV = I
         I     = INEXT
         INEXT = INEXT+DIR

         OLDSMT = SMT
         SMT = S(I) - TAU
         GPLUS(I)  = G(I) + SMT
         IF( ABS(GPLUS(I)).LT.PIVMIN )  GPLUS(I) = -PIVMIN
*
*        Compute next s after a block was broken
*        Note: Recall that NGN(IPREV) is undefined
*
         IF( G(IPREV).EQ.ZERO )THEN
            S(INEXT) = - GNSQ(I) / GPLUS(I)
         ELSE
            IF( (ABS(S(I)).LE.SBOUND) .OR.
     $          (SIGN(ONE,G(IPREV)).NE.SIGN(ONE,GPLUS(IPREV))) )
     $      THEN
               X = G(IPREV)*SMT + GNSQ(IPREV)
            ELSE
               X = -S(I)*OLDSMT - G(IPREV)*TAU
            ENDIF
            S(INEXT) = (GNSQ(I) * X) / (BDET(IPREV) * GPLUS(I))
         ENDIF
         I = INEXT
         INEXT = INEXT+DIR
      ENDDO
      END SUBROUTINE DLAXRT_STAT
*
************************************************************************
