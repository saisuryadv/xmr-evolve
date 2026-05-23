      SUBROUTINE DLAXRM_STAT2(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, ANEGC, AAUX
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(2)
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  ANEGC(2)
      DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(2)
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
*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  ONE  = 1.0D0
*
*     .. Locals ..
*
      DOUBLE PRECISION  SGNGIP, NEGPIV
      INTEGER           I, IPREV

      DOUBLE PRECISION  AUX1, GPL1, SMT1, OLDSMT1, SBOUND1
      DOUBLE PRECISION  AUX2, GPL2, SMT2, OLDSMT2, SBOUND2
      INTEGER  NEGC1
      INTEGER  NEGC2
      LOGICAL  BRNCHI1
      LOGICAL  BRNCHI2
*
*   ====== Executable Statements ========================================
*
      NEGPIV = -PIVMIN
      IF( DIR.EQ.1 )THEN
         I    = 1
      ELSE
         I    = N
      ENDIF
      SBOUND1 = MIN( PRBRK, N ) * ABS(TAU(1))
      SBOUND2 = MIN( PRBRK, N ) * ABS(TAU(2))
      NEGC1 = 0
      NEGC2 = 0
      AUX1 = ZERO
      AUX2 = ZERO
      DO
         DO
            IF( I .EQ. LBBEGK(IBB) )THEN
               EXIT
            ENDIF
            AUX1 = AUX1 - TAU(1)
            AUX2 = AUX2 - TAU(2)
            GPL1 = G(I) + AUX1
            GPL2 = G(I) + AUX2
            IF( ABS(GPL1) .LT. PIVMIN )  GPL1 = NEGPIV
            IF( ABS(GPL2) .LT. PIVMIN )  GPL2 = NEGPIV
            IF( GPL1 .LT. ZERO )  NEGC1 = NEGC1 + 1
            IF( GPL2 .LT. ZERO )  NEGC2 = NEGC2 + 1
            AUX1 = NGN(I) * ( AUX1 / GPL1 )
            AUX2 = NGN(I) * ( AUX2 / GPL2 )
            I = I + DIR
         ENDDO
         IF( I .EQ. K )THEN
            EXIT
         ENDIF
*        Note: We did not use SMT yet

         IBB = IBB + DIR

         SMT1 = AUX1 - TAU(1)
         SMT2 = AUX2 - TAU(2)
         GPL1 = G(I) + SMT1
         GPL2 = G(I) + SMT2
         IF( ABS(GPL1) .LT. PIVMIN )  GPL1 = NEGPIV
         IF( ABS(GPL2) .LT. PIVMIN )  GPL2 = NEGPIV
         IF( GPL1 .LT. ZERO )  NEGC1 = NEGC1 + 1
         IF( GPL2 .LT. ZERO )  NEGC2 = NEGC2 + 1
         AUX1 = - GNSQ(I) / GPL1
         AUX2 = - GNSQ(I) / GPL2
         IPREV = I
         I = I+DIR
         IF( I .EQ. K )THEN
            EXIT
         ENDIF
         SGNGIP = SIGN(ONE,G(IPREV))
         OLDSMT1 = SMT1
         BRNCHI1 = ( ABS(AUX1).LE.SBOUND1 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL1) )
         OLDSMT2 = SMT2
         BRNCHI2 = ( ABS(AUX2).LE.SBOUND2 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL2) )
         SMT1 = AUX1 - TAU(1)
         SMT2 = AUX2 - TAU(2)
         GPL1 = G(I) + SMT1
         GPL2 = G(I) + SMT2
         IF( ABS(GPL1) .LT. PIVMIN )  GPL1 = NEGPIV
         IF( ABS(GPL2) .LT. PIVMIN )  GPL2 = NEGPIV
         IF( GPL1 .LT. ZERO )  NEGC1 = NEGC1 + 1
         IF( GPL2 .LT. ZERO )  NEGC2 = NEGC2 + 1
         IF( G(IPREV) .EQ. ZERO )THEN
            AUX1 = - GNSQ(I) / GPL1
            AUX2 = - GNSQ(I) / GPL2
         ELSE
            IF( BRNCHI1 )THEN
               AUX1 = G(IPREV)*SMT1 + GNSQ(IPREV)
            ELSE
               AUX1 = -AUX1*OLDSMT1 - G(IPREV)*TAU(1)
            ENDIF
            IF( BRNCHI2 )THEN
               AUX2 = G(IPREV)*SMT2 + GNSQ(IPREV)
            ELSE
               AUX2 = -AUX2*OLDSMT2 - G(IPREV)*TAU(2)
            ENDIF
            AUX1 = (GNSQ(I) * AUX1) / (BDET(IPREV) * GPL1)
            AUX2 = (GNSQ(I) * AUX2) / (BDET(IPREV) * GPL2)
         ENDIF
         I = I + DIR
      ENDDO

      ANEGC(1) = NEGC1
      ANEGC(2) = NEGC2
      AAUX(1) = AUX1
      AAUX(2) = AUX2
      END SUBROUTINE DLAXRM_STAT2
*
************************************************************************
