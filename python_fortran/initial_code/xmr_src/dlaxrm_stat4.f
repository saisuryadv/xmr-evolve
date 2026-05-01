      SUBROUTINE DLAXRM_STAT4(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, ANEGC, AAUX
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(4)
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  ANEGC(4)
      DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(4)
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
      DOUBLE PRECISION  AUX3, GPL3, SMT3, OLDSMT3, SBOUND3
      DOUBLE PRECISION  AUX4, GPL4, SMT4, OLDSMT4, SBOUND4
      INTEGER  NEGC1
      INTEGER  NEGC2
      INTEGER  NEGC3
      INTEGER  NEGC4
      LOGICAL  BRNCHI1
      LOGICAL  BRNCHI2
      LOGICAL  BRNCHI3
      LOGICAL  BRNCHI4
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
      SBOUND3 = MIN( PRBRK, N ) * ABS(TAU(3))
      SBOUND4 = MIN( PRBRK, N ) * ABS(TAU(4))
      NEGC1 = 0
      NEGC2 = 0
      NEGC3 = 0
      NEGC4 = 0
      AUX1 = ZERO
      AUX2 = ZERO
      AUX3 = ZERO
      AUX4 = ZERO
      DO
         DO
            IF( I .EQ. LBBEGK(IBB) )THEN
               EXIT
            ENDIF
            AUX1 = AUX1 - TAU(1)
            AUX2 = AUX2 - TAU(2)
            AUX3 = AUX3 - TAU(3)
            AUX4 = AUX4 - TAU(4)
            GPL1 = G(I) + AUX1
            GPL2 = G(I) + AUX2
            GPL3 = G(I) + AUX3
            GPL4 = G(I) + AUX4
            IF( ABS(GPL1) .LT. PIVMIN )  GPL1 = NEGPIV
            IF( ABS(GPL2) .LT. PIVMIN )  GPL2 = NEGPIV
            IF( ABS(GPL3) .LT. PIVMIN )  GPL3 = NEGPIV
            IF( ABS(GPL4) .LT. PIVMIN )  GPL4 = NEGPIV
            IF( GPL1 .LT. ZERO )  NEGC1 = NEGC1 + 1
            IF( GPL2 .LT. ZERO )  NEGC2 = NEGC2 + 1
            IF( GPL3 .LT. ZERO )  NEGC3 = NEGC3 + 1
            IF( GPL4 .LT. ZERO )  NEGC4 = NEGC4 + 1
            AUX1 = NGN(I) * ( AUX1 / GPL1 )
            AUX2 = NGN(I) * ( AUX2 / GPL2 )
            AUX3 = NGN(I) * ( AUX3 / GPL3 )
            AUX4 = NGN(I) * ( AUX4 / GPL4 )
            I = I + DIR
         ENDDO
         IF( I .EQ. K )THEN
            EXIT
         ENDIF
*        Note: We did not use SMT yet

         IBB = IBB + DIR

         SMT1 = AUX1 - TAU(1)
         SMT2 = AUX2 - TAU(2)
         SMT3 = AUX3 - TAU(3)
         SMT4 = AUX4 - TAU(4)
         GPL1 = G(I) + SMT1
         GPL2 = G(I) + SMT2
         GPL3 = G(I) + SMT3
         GPL4 = G(I) + SMT4
         IF( ABS(GPL1) .LT. PIVMIN )  GPL1 = NEGPIV
         IF( ABS(GPL2) .LT. PIVMIN )  GPL2 = NEGPIV
         IF( ABS(GPL3) .LT. PIVMIN )  GPL3 = NEGPIV
         IF( ABS(GPL4) .LT. PIVMIN )  GPL4 = NEGPIV
         IF( GPL1 .LT. ZERO )  NEGC1 = NEGC1 + 1
         IF( GPL2 .LT. ZERO )  NEGC2 = NEGC2 + 1
         IF( GPL3 .LT. ZERO )  NEGC3 = NEGC3 + 1
         IF( GPL4 .LT. ZERO )  NEGC4 = NEGC4 + 1
         AUX1 = - GNSQ(I) / GPL1
         AUX2 = - GNSQ(I) / GPL2
         AUX3 = - GNSQ(I) / GPL3
         AUX4 = - GNSQ(I) / GPL4
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
         OLDSMT3 = SMT3
         BRNCHI3 = ( ABS(AUX3).LE.SBOUND3 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL3) )
         OLDSMT4 = SMT4
         BRNCHI4 = ( ABS(AUX4).LE.SBOUND4 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL4) )
         SMT1 = AUX1 - TAU(1)
         SMT2 = AUX2 - TAU(2)
         SMT3 = AUX3 - TAU(3)
         SMT4 = AUX4 - TAU(4)
         GPL1 = G(I) + SMT1
         GPL2 = G(I) + SMT2
         GPL3 = G(I) + SMT3
         GPL4 = G(I) + SMT4
         IF( ABS(GPL1) .LT. PIVMIN )  GPL1 = NEGPIV
         IF( ABS(GPL2) .LT. PIVMIN )  GPL2 = NEGPIV
         IF( ABS(GPL3) .LT. PIVMIN )  GPL3 = NEGPIV
         IF( ABS(GPL4) .LT. PIVMIN )  GPL4 = NEGPIV
         IF( GPL1 .LT. ZERO )  NEGC1 = NEGC1 + 1
         IF( GPL2 .LT. ZERO )  NEGC2 = NEGC2 + 1
         IF( GPL3 .LT. ZERO )  NEGC3 = NEGC3 + 1
         IF( GPL4 .LT. ZERO )  NEGC4 = NEGC4 + 1
         IF( G(IPREV) .EQ. ZERO )THEN
            AUX1 = - GNSQ(I) / GPL1
            AUX2 = - GNSQ(I) / GPL2
            AUX3 = - GNSQ(I) / GPL3
            AUX4 = - GNSQ(I) / GPL4
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
            IF( BRNCHI3 )THEN
               AUX3 = G(IPREV)*SMT3 + GNSQ(IPREV)
            ELSE
               AUX3 = -AUX3*OLDSMT3 - G(IPREV)*TAU(3)
            ENDIF
            IF( BRNCHI4 )THEN
               AUX4 = G(IPREV)*SMT4 + GNSQ(IPREV)
            ELSE
               AUX4 = -AUX4*OLDSMT4 - G(IPREV)*TAU(4)
            ENDIF
            AUX1 = (GNSQ(I) * AUX1) / (BDET(IPREV) * GPL1)
            AUX2 = (GNSQ(I) * AUX2) / (BDET(IPREV) * GPL2)
            AUX3 = (GNSQ(I) * AUX3) / (BDET(IPREV) * GPL3)
            AUX4 = (GNSQ(I) * AUX4) / (BDET(IPREV) * GPL4)
         ENDIF
         I = I + DIR
      ENDDO

      ANEGC(1) = NEGC1
      ANEGC(2) = NEGC2
      ANEGC(3) = NEGC3
      ANEGC(4) = NEGC4
      AAUX(1) = AUX1
      AAUX(2) = AUX2
      AAUX(3) = AUX3
      AAUX(4) = AUX4
      END SUBROUTINE DLAXRM_STAT4
*
************************************************************************
