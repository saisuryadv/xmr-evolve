      SUBROUTINE DLAXRM_STAT8(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, ANEGC, AAUX
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(8)
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  ANEGC(8)
      DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(8)
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
      DOUBLE PRECISION  AUX5, GPL5, SMT5, OLDSMT5, SBOUND5
      DOUBLE PRECISION  AUX6, GPL6, SMT6, OLDSMT6, SBOUND6
      DOUBLE PRECISION  AUX7, GPL7, SMT7, OLDSMT7, SBOUND7
      DOUBLE PRECISION  AUX8, GPL8, SMT8, OLDSMT8, SBOUND8
      INTEGER  NEGC1
      INTEGER  NEGC2
      INTEGER  NEGC3
      INTEGER  NEGC4
      INTEGER  NEGC5
      INTEGER  NEGC6
      INTEGER  NEGC7
      INTEGER  NEGC8
      LOGICAL  BRNCHI1
      LOGICAL  BRNCHI2
      LOGICAL  BRNCHI3
      LOGICAL  BRNCHI4
      LOGICAL  BRNCHI5
      LOGICAL  BRNCHI6
      LOGICAL  BRNCHI7
      LOGICAL  BRNCHI8
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
      SBOUND5 = MIN( PRBRK, N ) * ABS(TAU(5))
      SBOUND6 = MIN( PRBRK, N ) * ABS(TAU(6))
      SBOUND7 = MIN( PRBRK, N ) * ABS(TAU(7))
      SBOUND8 = MIN( PRBRK, N ) * ABS(TAU(8))
      NEGC1 = 0
      NEGC2 = 0
      NEGC3 = 0
      NEGC4 = 0
      NEGC5 = 0
      NEGC6 = 0
      NEGC7 = 0
      NEGC8 = 0
      AUX1 = ZERO
      AUX2 = ZERO
      AUX3 = ZERO
      AUX4 = ZERO
      AUX5 = ZERO
      AUX6 = ZERO
      AUX7 = ZERO
      AUX8 = ZERO
      DO
         DO
            IF( I .EQ. LBBEGK(IBB) )THEN
               EXIT
            ENDIF
            AUX1 = AUX1 - TAU(1)
            AUX2 = AUX2 - TAU(2)
            AUX3 = AUX3 - TAU(3)
            AUX4 = AUX4 - TAU(4)
            AUX5 = AUX5 - TAU(5)
            AUX6 = AUX6 - TAU(6)
            AUX7 = AUX7 - TAU(7)
            AUX8 = AUX8 - TAU(8)
            GPL1 = G(I) + AUX1
            GPL2 = G(I) + AUX2
            GPL3 = G(I) + AUX3
            GPL4 = G(I) + AUX4
            GPL5 = G(I) + AUX5
            GPL6 = G(I) + AUX6
            GPL7 = G(I) + AUX7
            GPL8 = G(I) + AUX8
            IF( ABS(GPL1) .LT. PIVMIN )  GPL1 = NEGPIV
            IF( ABS(GPL2) .LT. PIVMIN )  GPL2 = NEGPIV
            IF( ABS(GPL3) .LT. PIVMIN )  GPL3 = NEGPIV
            IF( ABS(GPL4) .LT. PIVMIN )  GPL4 = NEGPIV
            IF( ABS(GPL5) .LT. PIVMIN )  GPL5 = NEGPIV
            IF( ABS(GPL6) .LT. PIVMIN )  GPL6 = NEGPIV
            IF( ABS(GPL7) .LT. PIVMIN )  GPL7 = NEGPIV
            IF( ABS(GPL8) .LT. PIVMIN )  GPL8 = NEGPIV
            IF( GPL1 .LT. ZERO )  NEGC1 = NEGC1 + 1
            IF( GPL2 .LT. ZERO )  NEGC2 = NEGC2 + 1
            IF( GPL3 .LT. ZERO )  NEGC3 = NEGC3 + 1
            IF( GPL4 .LT. ZERO )  NEGC4 = NEGC4 + 1
            IF( GPL5 .LT. ZERO )  NEGC5 = NEGC5 + 1
            IF( GPL6 .LT. ZERO )  NEGC6 = NEGC6 + 1
            IF( GPL7 .LT. ZERO )  NEGC7 = NEGC7 + 1
            IF( GPL8 .LT. ZERO )  NEGC8 = NEGC8 + 1
            AUX1 = NGN(I) * ( AUX1 / GPL1 )
            AUX2 = NGN(I) * ( AUX2 / GPL2 )
            AUX3 = NGN(I) * ( AUX3 / GPL3 )
            AUX4 = NGN(I) * ( AUX4 / GPL4 )
            AUX5 = NGN(I) * ( AUX5 / GPL5 )
            AUX6 = NGN(I) * ( AUX6 / GPL6 )
            AUX7 = NGN(I) * ( AUX7 / GPL7 )
            AUX8 = NGN(I) * ( AUX8 / GPL8 )
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
         SMT5 = AUX5 - TAU(5)
         SMT6 = AUX6 - TAU(6)
         SMT7 = AUX7 - TAU(7)
         SMT8 = AUX8 - TAU(8)
         GPL1 = G(I) + SMT1
         GPL2 = G(I) + SMT2
         GPL3 = G(I) + SMT3
         GPL4 = G(I) + SMT4
         GPL5 = G(I) + SMT5
         GPL6 = G(I) + SMT6
         GPL7 = G(I) + SMT7
         GPL8 = G(I) + SMT8
         IF( ABS(GPL1) .LT. PIVMIN )  GPL1 = NEGPIV
         IF( ABS(GPL2) .LT. PIVMIN )  GPL2 = NEGPIV
         IF( ABS(GPL3) .LT. PIVMIN )  GPL3 = NEGPIV
         IF( ABS(GPL4) .LT. PIVMIN )  GPL4 = NEGPIV
         IF( ABS(GPL5) .LT. PIVMIN )  GPL5 = NEGPIV
         IF( ABS(GPL6) .LT. PIVMIN )  GPL6 = NEGPIV
         IF( ABS(GPL7) .LT. PIVMIN )  GPL7 = NEGPIV
         IF( ABS(GPL8) .LT. PIVMIN )  GPL8 = NEGPIV
         IF( GPL1 .LT. ZERO )  NEGC1 = NEGC1 + 1
         IF( GPL2 .LT. ZERO )  NEGC2 = NEGC2 + 1
         IF( GPL3 .LT. ZERO )  NEGC3 = NEGC3 + 1
         IF( GPL4 .LT. ZERO )  NEGC4 = NEGC4 + 1
         IF( GPL5 .LT. ZERO )  NEGC5 = NEGC5 + 1
         IF( GPL6 .LT. ZERO )  NEGC6 = NEGC6 + 1
         IF( GPL7 .LT. ZERO )  NEGC7 = NEGC7 + 1
         IF( GPL8 .LT. ZERO )  NEGC8 = NEGC8 + 1
         AUX1 = - GNSQ(I) / GPL1
         AUX2 = - GNSQ(I) / GPL2
         AUX3 = - GNSQ(I) / GPL3
         AUX4 = - GNSQ(I) / GPL4
         AUX5 = - GNSQ(I) / GPL5
         AUX6 = - GNSQ(I) / GPL6
         AUX7 = - GNSQ(I) / GPL7
         AUX8 = - GNSQ(I) / GPL8
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
         OLDSMT5 = SMT5
         BRNCHI5 = ( ABS(AUX5).LE.SBOUND5 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL5) )
         OLDSMT6 = SMT6
         BRNCHI6 = ( ABS(AUX6).LE.SBOUND6 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL6) )
         OLDSMT7 = SMT7
         BRNCHI7 = ( ABS(AUX7).LE.SBOUND7 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL7) )
         OLDSMT8 = SMT8
         BRNCHI8 = ( ABS(AUX8).LE.SBOUND8 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL8) )
         SMT1 = AUX1 - TAU(1)
         SMT2 = AUX2 - TAU(2)
         SMT3 = AUX3 - TAU(3)
         SMT4 = AUX4 - TAU(4)
         SMT5 = AUX5 - TAU(5)
         SMT6 = AUX6 - TAU(6)
         SMT7 = AUX7 - TAU(7)
         SMT8 = AUX8 - TAU(8)
         GPL1 = G(I) + SMT1
         GPL2 = G(I) + SMT2
         GPL3 = G(I) + SMT3
         GPL4 = G(I) + SMT4
         GPL5 = G(I) + SMT5
         GPL6 = G(I) + SMT6
         GPL7 = G(I) + SMT7
         GPL8 = G(I) + SMT8
         IF( ABS(GPL1) .LT. PIVMIN )  GPL1 = NEGPIV
         IF( ABS(GPL2) .LT. PIVMIN )  GPL2 = NEGPIV
         IF( ABS(GPL3) .LT. PIVMIN )  GPL3 = NEGPIV
         IF( ABS(GPL4) .LT. PIVMIN )  GPL4 = NEGPIV
         IF( ABS(GPL5) .LT. PIVMIN )  GPL5 = NEGPIV
         IF( ABS(GPL6) .LT. PIVMIN )  GPL6 = NEGPIV
         IF( ABS(GPL7) .LT. PIVMIN )  GPL7 = NEGPIV
         IF( ABS(GPL8) .LT. PIVMIN )  GPL8 = NEGPIV
         IF( GPL1 .LT. ZERO )  NEGC1 = NEGC1 + 1
         IF( GPL2 .LT. ZERO )  NEGC2 = NEGC2 + 1
         IF( GPL3 .LT. ZERO )  NEGC3 = NEGC3 + 1
         IF( GPL4 .LT. ZERO )  NEGC4 = NEGC4 + 1
         IF( GPL5 .LT. ZERO )  NEGC5 = NEGC5 + 1
         IF( GPL6 .LT. ZERO )  NEGC6 = NEGC6 + 1
         IF( GPL7 .LT. ZERO )  NEGC7 = NEGC7 + 1
         IF( GPL8 .LT. ZERO )  NEGC8 = NEGC8 + 1
         IF( G(IPREV) .EQ. ZERO )THEN
            AUX1 = - GNSQ(I) / GPL1
            AUX2 = - GNSQ(I) / GPL2
            AUX3 = - GNSQ(I) / GPL3
            AUX4 = - GNSQ(I) / GPL4
            AUX5 = - GNSQ(I) / GPL5
            AUX6 = - GNSQ(I) / GPL6
            AUX7 = - GNSQ(I) / GPL7
            AUX8 = - GNSQ(I) / GPL8
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
            IF( BRNCHI5 )THEN
               AUX5 = G(IPREV)*SMT5 + GNSQ(IPREV)
            ELSE
               AUX5 = -AUX5*OLDSMT5 - G(IPREV)*TAU(5)
            ENDIF
            IF( BRNCHI6 )THEN
               AUX6 = G(IPREV)*SMT6 + GNSQ(IPREV)
            ELSE
               AUX6 = -AUX6*OLDSMT6 - G(IPREV)*TAU(6)
            ENDIF
            IF( BRNCHI7 )THEN
               AUX7 = G(IPREV)*SMT7 + GNSQ(IPREV)
            ELSE
               AUX7 = -AUX7*OLDSMT7 - G(IPREV)*TAU(7)
            ENDIF
            IF( BRNCHI8 )THEN
               AUX8 = G(IPREV)*SMT8 + GNSQ(IPREV)
            ELSE
               AUX8 = -AUX8*OLDSMT8 - G(IPREV)*TAU(8)
            ENDIF
            AUX1 = (GNSQ(I) * AUX1) / (BDET(IPREV) * GPL1)
            AUX2 = (GNSQ(I) * AUX2) / (BDET(IPREV) * GPL2)
            AUX3 = (GNSQ(I) * AUX3) / (BDET(IPREV) * GPL3)
            AUX4 = (GNSQ(I) * AUX4) / (BDET(IPREV) * GPL4)
            AUX5 = (GNSQ(I) * AUX5) / (BDET(IPREV) * GPL5)
            AUX6 = (GNSQ(I) * AUX6) / (BDET(IPREV) * GPL6)
            AUX7 = (GNSQ(I) * AUX7) / (BDET(IPREV) * GPL7)
            AUX8 = (GNSQ(I) * AUX8) / (BDET(IPREV) * GPL8)
         ENDIF
         I = I + DIR
      ENDDO

      ANEGC(1) = NEGC1
      ANEGC(2) = NEGC2
      ANEGC(3) = NEGC3
      ANEGC(4) = NEGC4
      ANEGC(5) = NEGC5
      ANEGC(6) = NEGC6
      ANEGC(7) = NEGC7
      ANEGC(8) = NEGC8
      AAUX(1) = AUX1
      AAUX(2) = AUX2
      AAUX(3) = AUX3
      AAUX(4) = AUX4
      AAUX(5) = AUX5
      AAUX(6) = AUX6
      AAUX(7) = AUX7
      AAUX(8) = AUX8
      END SUBROUTINE DLAXRM_STAT8
*
************************************************************************
