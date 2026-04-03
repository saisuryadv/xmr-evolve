      SUBROUTINE DLAXRT_PROG(
     $             N, J1, J2, I0, I1, G, OMEGA, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, OFF, GPLUS, P
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, J1, J2, I0, I1
      INTEGER,          INTENT(IN)  ::  OMEGA(0:N+1)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU, OFF
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      DOUBLE PRECISION, INTENT(OUT)  ::  GPLUS(I0:I1), P(N)
*
*  Purpose
*  =======
*
*     Computes a progressive factorization from J1 to J2.
*
*     J1 < J2 and I0 = 1 and I1 = N-1
*        assumes J1 is twist in source
*        produces GPLUS(J1),...,GPLUS(J2-1) and P(J1+1),...,P(J2)
*
*     J1 > J2 and I0 = 2 and I1 = N
*        assumes J2 is twist in source
*        produces GPLUS(J1),...,GPLUS(J2+1) and P(J1-1),...,P(J2)
*
*  ======================================================================
*
*     .. Declarations ..
*
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO=0.0D0 )
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
*     .. Local ..
*
      INTEGER I, IPREV, INEXT, DIR, J
      DOUBLE PRECISION PBOUND, AUX
*
*     -- Executable Statements -----------------------------------------
*

      I = J1
      J = J2

      IF( J1.LT.J2 .AND. I0.EQ.1 .AND. I1.EQ.N-1 )THEN
         DIR = +1
      ELSEIF( J2.LT.J1 .AND. I0.EQ.2 .AND. I1.EQ.N )THEN
         DIR = -1
      ELSE
         I   = -1
         DIR = 0
      ENDIF

      PBOUND = MIN(PRBRK,N) * ABS(TAU)

      INEXT = I+DIR

      AUX = OFF + (G(I) - TAU)
      IF( OMEGA(I).EQ.1 )THEN
         GPLUS(I) = AUX
         IF( ABS(GPLUS(I)).LT.PIVMIN )  GPLUS(I) = -PIVMIN
         P(INEXT) = G(INEXT) - GNSQ(INEXT)/GPLUS(I)
         I = INEXT
         INEXT = INEXT+DIR
         AUX = P(I) - TAU
      ENDIF
      DO
         DO
            IF( (I.EQ.J) .OR. (OMEGA(INEXT).NE.0) )THEN
               EXIT
            ENDIF
            GPLUS(I) = NGN(INEXT) + AUX
            IF( ABS(GPLUS(I)).LT.PIVMIN )  GPLUS(I) = -PIVMIN
            P(INEXT) = (G(INEXT)/GPLUS(I)) * AUX
            I     = INEXT
            INEXT = INEXT+DIR
            AUX = P(I) - TAU
         ENDDO
         IF( I.EQ.J )THEN
            EXIT
         ENDIF

         GPLUS(I) = NGN(INEXT) + AUX
         IF( ABS(GPLUS(I)).LT.PIVMIN )  GPLUS(I) = -PIVMIN
         P(INEXT) = G(INEXT) - GNSQ(INEXT)/GPLUS(I)
         IF( INEXT.EQ.J )THEN
            EXIT
         ENDIF
         IPREV = I
         I     = INEXT
         INEXT = INEXT+DIR

         GPLUS(I) = P(I) - TAU
         IF( ABS(GPLUS(I)).LT.PIVMIN )  GPLUS(I) = -PIVMIN

         IF( (ABS(P(I)).LE.PBOUND) .OR. G(INEXT).EQ.ZERO ) THEN
            P(INEXT) = G(INEXT) - GNSQ(INEXT)/GPLUS(I)
         ELSE
            P(INEXT) = ( BDET(INEXT)*(AUX/GPLUS(IPREV))
     $                   - G(INEXT)*TAU )
     $                 / GPLUS(I)
         ENDIF
         I     = INEXT
         INEXT = INEXT+DIR
         AUX = P(I) - TAU
      ENDDO
      END SUBROUTINE DLAXRT_PROG
*
************************************************************************
