      SUBROUTINE DLAXRM_STAT16(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, ANEGC, AAUX
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(16)
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  ANEGC(16)
      DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(16)
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
      DOUBLE PRECISION  AUX9, GPL9, SMT9, OLDSMT9, SBOUND9
      DOUBLE PRECISION  AUX10, GPL10, SMT10, OLDSMT10, SBOUND10
      DOUBLE PRECISION  AUX11, GPL11, SMT11, OLDSMT11, SBOUND11
      DOUBLE PRECISION  AUX12, GPL12, SMT12, OLDSMT12, SBOUND12
      DOUBLE PRECISION  AUX13, GPL13, SMT13, OLDSMT13, SBOUND13
      DOUBLE PRECISION  AUX14, GPL14, SMT14, OLDSMT14, SBOUND14
      DOUBLE PRECISION  AUX15, GPL15, SMT15, OLDSMT15, SBOUND15
      DOUBLE PRECISION  AUX16, GPL16, SMT16, OLDSMT16, SBOUND16
      INTEGER  NEGC1
      INTEGER  NEGC2
      INTEGER  NEGC3
      INTEGER  NEGC4
      INTEGER  NEGC5
      INTEGER  NEGC6
      INTEGER  NEGC7
      INTEGER  NEGC8
      INTEGER  NEGC9
      INTEGER  NEGC10
      INTEGER  NEGC11
      INTEGER  NEGC12
      INTEGER  NEGC13
      INTEGER  NEGC14
      INTEGER  NEGC15
      INTEGER  NEGC16
      LOGICAL  BRNCHI1
      LOGICAL  BRNCHI2
      LOGICAL  BRNCHI3
      LOGICAL  BRNCHI4
      LOGICAL  BRNCHI5
      LOGICAL  BRNCHI6
      LOGICAL  BRNCHI7
      LOGICAL  BRNCHI8
      LOGICAL  BRNCHI9
      LOGICAL  BRNCHI10
      LOGICAL  BRNCHI11
      LOGICAL  BRNCHI12
      LOGICAL  BRNCHI13
      LOGICAL  BRNCHI14
      LOGICAL  BRNCHI15
      LOGICAL  BRNCHI16
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
      SBOUND9 = MIN( PRBRK, N ) * ABS(TAU(9))
      SBOUND10 = MIN( PRBRK, N ) * ABS(TAU(10))
      SBOUND11 = MIN( PRBRK, N ) * ABS(TAU(11))
      SBOUND12 = MIN( PRBRK, N ) * ABS(TAU(12))
      SBOUND13 = MIN( PRBRK, N ) * ABS(TAU(13))
      SBOUND14 = MIN( PRBRK, N ) * ABS(TAU(14))
      SBOUND15 = MIN( PRBRK, N ) * ABS(TAU(15))
      SBOUND16 = MIN( PRBRK, N ) * ABS(TAU(16))
      NEGC1 = 0
      NEGC2 = 0
      NEGC3 = 0
      NEGC4 = 0
      NEGC5 = 0
      NEGC6 = 0
      NEGC7 = 0
      NEGC8 = 0
      NEGC9 = 0
      NEGC10 = 0
      NEGC11 = 0
      NEGC12 = 0
      NEGC13 = 0
      NEGC14 = 0
      NEGC15 = 0
      NEGC16 = 0
      AUX1 = ZERO
      AUX2 = ZERO
      AUX3 = ZERO
      AUX4 = ZERO
      AUX5 = ZERO
      AUX6 = ZERO
      AUX7 = ZERO
      AUX8 = ZERO
      AUX9 = ZERO
      AUX10 = ZERO
      AUX11 = ZERO
      AUX12 = ZERO
      AUX13 = ZERO
      AUX14 = ZERO
      AUX15 = ZERO
      AUX16 = ZERO
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
            AUX9 = AUX9 - TAU(9)
            AUX10 = AUX10 - TAU(10)
            AUX11 = AUX11 - TAU(11)
            AUX12 = AUX12 - TAU(12)
            AUX13 = AUX13 - TAU(13)
            AUX14 = AUX14 - TAU(14)
            AUX15 = AUX15 - TAU(15)
            AUX16 = AUX16 - TAU(16)
            GPL1 = G(I) + AUX1
            GPL2 = G(I) + AUX2
            GPL3 = G(I) + AUX3
            GPL4 = G(I) + AUX4
            GPL5 = G(I) + AUX5
            GPL6 = G(I) + AUX6
            GPL7 = G(I) + AUX7
            GPL8 = G(I) + AUX8
            GPL9 = G(I) + AUX9
            GPL10 = G(I) + AUX10
            GPL11 = G(I) + AUX11
            GPL12 = G(I) + AUX12
            GPL13 = G(I) + AUX13
            GPL14 = G(I) + AUX14
            GPL15 = G(I) + AUX15
            GPL16 = G(I) + AUX16
            IF( ABS(GPL1) .LT. PIVMIN )  GPL1 = NEGPIV
            IF( ABS(GPL2) .LT. PIVMIN )  GPL2 = NEGPIV
            IF( ABS(GPL3) .LT. PIVMIN )  GPL3 = NEGPIV
            IF( ABS(GPL4) .LT. PIVMIN )  GPL4 = NEGPIV
            IF( ABS(GPL5) .LT. PIVMIN )  GPL5 = NEGPIV
            IF( ABS(GPL6) .LT. PIVMIN )  GPL6 = NEGPIV
            IF( ABS(GPL7) .LT. PIVMIN )  GPL7 = NEGPIV
            IF( ABS(GPL8) .LT. PIVMIN )  GPL8 = NEGPIV
            IF( ABS(GPL9) .LT. PIVMIN )  GPL9 = NEGPIV
            IF( ABS(GPL10) .LT. PIVMIN )  GPL10 = NEGPIV
            IF( ABS(GPL11) .LT. PIVMIN )  GPL11 = NEGPIV
            IF( ABS(GPL12) .LT. PIVMIN )  GPL12 = NEGPIV
            IF( ABS(GPL13) .LT. PIVMIN )  GPL13 = NEGPIV
            IF( ABS(GPL14) .LT. PIVMIN )  GPL14 = NEGPIV
            IF( ABS(GPL15) .LT. PIVMIN )  GPL15 = NEGPIV
            IF( ABS(GPL16) .LT. PIVMIN )  GPL16 = NEGPIV
            IF( GPL1 .LT. ZERO )  NEGC1 = NEGC1 + 1
            IF( GPL2 .LT. ZERO )  NEGC2 = NEGC2 + 1
            IF( GPL3 .LT. ZERO )  NEGC3 = NEGC3 + 1
            IF( GPL4 .LT. ZERO )  NEGC4 = NEGC4 + 1
            IF( GPL5 .LT. ZERO )  NEGC5 = NEGC5 + 1
            IF( GPL6 .LT. ZERO )  NEGC6 = NEGC6 + 1
            IF( GPL7 .LT. ZERO )  NEGC7 = NEGC7 + 1
            IF( GPL8 .LT. ZERO )  NEGC8 = NEGC8 + 1
            IF( GPL9 .LT. ZERO )  NEGC9 = NEGC9 + 1
            IF( GPL10 .LT. ZERO )  NEGC10 = NEGC10 + 1
            IF( GPL11 .LT. ZERO )  NEGC11 = NEGC11 + 1
            IF( GPL12 .LT. ZERO )  NEGC12 = NEGC12 + 1
            IF( GPL13 .LT. ZERO )  NEGC13 = NEGC13 + 1
            IF( GPL14 .LT. ZERO )  NEGC14 = NEGC14 + 1
            IF( GPL15 .LT. ZERO )  NEGC15 = NEGC15 + 1
            IF( GPL16 .LT. ZERO )  NEGC16 = NEGC16 + 1
            AUX1 = NGN(I) * ( AUX1 / GPL1 )
            AUX2 = NGN(I) * ( AUX2 / GPL2 )
            AUX3 = NGN(I) * ( AUX3 / GPL3 )
            AUX4 = NGN(I) * ( AUX4 / GPL4 )
            AUX5 = NGN(I) * ( AUX5 / GPL5 )
            AUX6 = NGN(I) * ( AUX6 / GPL6 )
            AUX7 = NGN(I) * ( AUX7 / GPL7 )
            AUX8 = NGN(I) * ( AUX8 / GPL8 )
            AUX9 = NGN(I) * ( AUX9 / GPL9 )
            AUX10 = NGN(I) * ( AUX10 / GPL10 )
            AUX11 = NGN(I) * ( AUX11 / GPL11 )
            AUX12 = NGN(I) * ( AUX12 / GPL12 )
            AUX13 = NGN(I) * ( AUX13 / GPL13 )
            AUX14 = NGN(I) * ( AUX14 / GPL14 )
            AUX15 = NGN(I) * ( AUX15 / GPL15 )
            AUX16 = NGN(I) * ( AUX16 / GPL16 )
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
         SMT9 = AUX9 - TAU(9)
         SMT10 = AUX10 - TAU(10)
         SMT11 = AUX11 - TAU(11)
         SMT12 = AUX12 - TAU(12)
         SMT13 = AUX13 - TAU(13)
         SMT14 = AUX14 - TAU(14)
         SMT15 = AUX15 - TAU(15)
         SMT16 = AUX16 - TAU(16)
         GPL1 = G(I) + SMT1
         GPL2 = G(I) + SMT2
         GPL3 = G(I) + SMT3
         GPL4 = G(I) + SMT4
         GPL5 = G(I) + SMT5
         GPL6 = G(I) + SMT6
         GPL7 = G(I) + SMT7
         GPL8 = G(I) + SMT8
         GPL9 = G(I) + SMT9
         GPL10 = G(I) + SMT10
         GPL11 = G(I) + SMT11
         GPL12 = G(I) + SMT12
         GPL13 = G(I) + SMT13
         GPL14 = G(I) + SMT14
         GPL15 = G(I) + SMT15
         GPL16 = G(I) + SMT16
         IF( ABS(GPL1) .LT. PIVMIN )  GPL1 = NEGPIV
         IF( ABS(GPL2) .LT. PIVMIN )  GPL2 = NEGPIV
         IF( ABS(GPL3) .LT. PIVMIN )  GPL3 = NEGPIV
         IF( ABS(GPL4) .LT. PIVMIN )  GPL4 = NEGPIV
         IF( ABS(GPL5) .LT. PIVMIN )  GPL5 = NEGPIV
         IF( ABS(GPL6) .LT. PIVMIN )  GPL6 = NEGPIV
         IF( ABS(GPL7) .LT. PIVMIN )  GPL7 = NEGPIV
         IF( ABS(GPL8) .LT. PIVMIN )  GPL8 = NEGPIV
         IF( ABS(GPL9) .LT. PIVMIN )  GPL9 = NEGPIV
         IF( ABS(GPL10) .LT. PIVMIN )  GPL10 = NEGPIV
         IF( ABS(GPL11) .LT. PIVMIN )  GPL11 = NEGPIV
         IF( ABS(GPL12) .LT. PIVMIN )  GPL12 = NEGPIV
         IF( ABS(GPL13) .LT. PIVMIN )  GPL13 = NEGPIV
         IF( ABS(GPL14) .LT. PIVMIN )  GPL14 = NEGPIV
         IF( ABS(GPL15) .LT. PIVMIN )  GPL15 = NEGPIV
         IF( ABS(GPL16) .LT. PIVMIN )  GPL16 = NEGPIV
         IF( GPL1 .LT. ZERO )  NEGC1 = NEGC1 + 1
         IF( GPL2 .LT. ZERO )  NEGC2 = NEGC2 + 1
         IF( GPL3 .LT. ZERO )  NEGC3 = NEGC3 + 1
         IF( GPL4 .LT. ZERO )  NEGC4 = NEGC4 + 1
         IF( GPL5 .LT. ZERO )  NEGC5 = NEGC5 + 1
         IF( GPL6 .LT. ZERO )  NEGC6 = NEGC6 + 1
         IF( GPL7 .LT. ZERO )  NEGC7 = NEGC7 + 1
         IF( GPL8 .LT. ZERO )  NEGC8 = NEGC8 + 1
         IF( GPL9 .LT. ZERO )  NEGC9 = NEGC9 + 1
         IF( GPL10 .LT. ZERO )  NEGC10 = NEGC10 + 1
         IF( GPL11 .LT. ZERO )  NEGC11 = NEGC11 + 1
         IF( GPL12 .LT. ZERO )  NEGC12 = NEGC12 + 1
         IF( GPL13 .LT. ZERO )  NEGC13 = NEGC13 + 1
         IF( GPL14 .LT. ZERO )  NEGC14 = NEGC14 + 1
         IF( GPL15 .LT. ZERO )  NEGC15 = NEGC15 + 1
         IF( GPL16 .LT. ZERO )  NEGC16 = NEGC16 + 1
         AUX1 = - GNSQ(I) / GPL1
         AUX2 = - GNSQ(I) / GPL2
         AUX3 = - GNSQ(I) / GPL3
         AUX4 = - GNSQ(I) / GPL4
         AUX5 = - GNSQ(I) / GPL5
         AUX6 = - GNSQ(I) / GPL6
         AUX7 = - GNSQ(I) / GPL7
         AUX8 = - GNSQ(I) / GPL8
         AUX9 = - GNSQ(I) / GPL9
         AUX10 = - GNSQ(I) / GPL10
         AUX11 = - GNSQ(I) / GPL11
         AUX12 = - GNSQ(I) / GPL12
         AUX13 = - GNSQ(I) / GPL13
         AUX14 = - GNSQ(I) / GPL14
         AUX15 = - GNSQ(I) / GPL15
         AUX16 = - GNSQ(I) / GPL16
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
         OLDSMT9 = SMT9
         BRNCHI9 = ( ABS(AUX9).LE.SBOUND9 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL9) )
         OLDSMT10 = SMT10
         BRNCHI10 = ( ABS(AUX10).LE.SBOUND10 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL10) )
         OLDSMT11 = SMT11
         BRNCHI11 = ( ABS(AUX11).LE.SBOUND11 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL11) )
         OLDSMT12 = SMT12
         BRNCHI12 = ( ABS(AUX12).LE.SBOUND12 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL12) )
         OLDSMT13 = SMT13
         BRNCHI13 = ( ABS(AUX13).LE.SBOUND13 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL13) )
         OLDSMT14 = SMT14
         BRNCHI14 = ( ABS(AUX14).LE.SBOUND14 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL14) )
         OLDSMT15 = SMT15
         BRNCHI15 = ( ABS(AUX15).LE.SBOUND15 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL15) )
         OLDSMT16 = SMT16
         BRNCHI16 = ( ABS(AUX16).LE.SBOUND16 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL16) )
         SMT1 = AUX1 - TAU(1)
         SMT2 = AUX2 - TAU(2)
         SMT3 = AUX3 - TAU(3)
         SMT4 = AUX4 - TAU(4)
         SMT5 = AUX5 - TAU(5)
         SMT6 = AUX6 - TAU(6)
         SMT7 = AUX7 - TAU(7)
         SMT8 = AUX8 - TAU(8)
         SMT9 = AUX9 - TAU(9)
         SMT10 = AUX10 - TAU(10)
         SMT11 = AUX11 - TAU(11)
         SMT12 = AUX12 - TAU(12)
         SMT13 = AUX13 - TAU(13)
         SMT14 = AUX14 - TAU(14)
         SMT15 = AUX15 - TAU(15)
         SMT16 = AUX16 - TAU(16)
         GPL1 = G(I) + SMT1
         GPL2 = G(I) + SMT2
         GPL3 = G(I) + SMT3
         GPL4 = G(I) + SMT4
         GPL5 = G(I) + SMT5
         GPL6 = G(I) + SMT6
         GPL7 = G(I) + SMT7
         GPL8 = G(I) + SMT8
         GPL9 = G(I) + SMT9
         GPL10 = G(I) + SMT10
         GPL11 = G(I) + SMT11
         GPL12 = G(I) + SMT12
         GPL13 = G(I) + SMT13
         GPL14 = G(I) + SMT14
         GPL15 = G(I) + SMT15
         GPL16 = G(I) + SMT16
         IF( ABS(GPL1) .LT. PIVMIN )  GPL1 = NEGPIV
         IF( ABS(GPL2) .LT. PIVMIN )  GPL2 = NEGPIV
         IF( ABS(GPL3) .LT. PIVMIN )  GPL3 = NEGPIV
         IF( ABS(GPL4) .LT. PIVMIN )  GPL4 = NEGPIV
         IF( ABS(GPL5) .LT. PIVMIN )  GPL5 = NEGPIV
         IF( ABS(GPL6) .LT. PIVMIN )  GPL6 = NEGPIV
         IF( ABS(GPL7) .LT. PIVMIN )  GPL7 = NEGPIV
         IF( ABS(GPL8) .LT. PIVMIN )  GPL8 = NEGPIV
         IF( ABS(GPL9) .LT. PIVMIN )  GPL9 = NEGPIV
         IF( ABS(GPL10) .LT. PIVMIN )  GPL10 = NEGPIV
         IF( ABS(GPL11) .LT. PIVMIN )  GPL11 = NEGPIV
         IF( ABS(GPL12) .LT. PIVMIN )  GPL12 = NEGPIV
         IF( ABS(GPL13) .LT. PIVMIN )  GPL13 = NEGPIV
         IF( ABS(GPL14) .LT. PIVMIN )  GPL14 = NEGPIV
         IF( ABS(GPL15) .LT. PIVMIN )  GPL15 = NEGPIV
         IF( ABS(GPL16) .LT. PIVMIN )  GPL16 = NEGPIV
         IF( GPL1 .LT. ZERO )  NEGC1 = NEGC1 + 1
         IF( GPL2 .LT. ZERO )  NEGC2 = NEGC2 + 1
         IF( GPL3 .LT. ZERO )  NEGC3 = NEGC3 + 1
         IF( GPL4 .LT. ZERO )  NEGC4 = NEGC4 + 1
         IF( GPL5 .LT. ZERO )  NEGC5 = NEGC5 + 1
         IF( GPL6 .LT. ZERO )  NEGC6 = NEGC6 + 1
         IF( GPL7 .LT. ZERO )  NEGC7 = NEGC7 + 1
         IF( GPL8 .LT. ZERO )  NEGC8 = NEGC8 + 1
         IF( GPL9 .LT. ZERO )  NEGC9 = NEGC9 + 1
         IF( GPL10 .LT. ZERO )  NEGC10 = NEGC10 + 1
         IF( GPL11 .LT. ZERO )  NEGC11 = NEGC11 + 1
         IF( GPL12 .LT. ZERO )  NEGC12 = NEGC12 + 1
         IF( GPL13 .LT. ZERO )  NEGC13 = NEGC13 + 1
         IF( GPL14 .LT. ZERO )  NEGC14 = NEGC14 + 1
         IF( GPL15 .LT. ZERO )  NEGC15 = NEGC15 + 1
         IF( GPL16 .LT. ZERO )  NEGC16 = NEGC16 + 1
         IF( G(IPREV) .EQ. ZERO )THEN
            AUX1 = - GNSQ(I) / GPL1
            AUX2 = - GNSQ(I) / GPL2
            AUX3 = - GNSQ(I) / GPL3
            AUX4 = - GNSQ(I) / GPL4
            AUX5 = - GNSQ(I) / GPL5
            AUX6 = - GNSQ(I) / GPL6
            AUX7 = - GNSQ(I) / GPL7
            AUX8 = - GNSQ(I) / GPL8
            AUX9 = - GNSQ(I) / GPL9
            AUX10 = - GNSQ(I) / GPL10
            AUX11 = - GNSQ(I) / GPL11
            AUX12 = - GNSQ(I) / GPL12
            AUX13 = - GNSQ(I) / GPL13
            AUX14 = - GNSQ(I) / GPL14
            AUX15 = - GNSQ(I) / GPL15
            AUX16 = - GNSQ(I) / GPL16
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
            IF( BRNCHI9 )THEN
               AUX9 = G(IPREV)*SMT9 + GNSQ(IPREV)
            ELSE
               AUX9 = -AUX9*OLDSMT9 - G(IPREV)*TAU(9)
            ENDIF
            IF( BRNCHI10 )THEN
               AUX10 = G(IPREV)*SMT10 + GNSQ(IPREV)
            ELSE
               AUX10 = -AUX10*OLDSMT10 - G(IPREV)*TAU(10)
            ENDIF
            IF( BRNCHI11 )THEN
               AUX11 = G(IPREV)*SMT11 + GNSQ(IPREV)
            ELSE
               AUX11 = -AUX11*OLDSMT11 - G(IPREV)*TAU(11)
            ENDIF
            IF( BRNCHI12 )THEN
               AUX12 = G(IPREV)*SMT12 + GNSQ(IPREV)
            ELSE
               AUX12 = -AUX12*OLDSMT12 - G(IPREV)*TAU(12)
            ENDIF
            IF( BRNCHI13 )THEN
               AUX13 = G(IPREV)*SMT13 + GNSQ(IPREV)
            ELSE
               AUX13 = -AUX13*OLDSMT13 - G(IPREV)*TAU(13)
            ENDIF
            IF( BRNCHI14 )THEN
               AUX14 = G(IPREV)*SMT14 + GNSQ(IPREV)
            ELSE
               AUX14 = -AUX14*OLDSMT14 - G(IPREV)*TAU(14)
            ENDIF
            IF( BRNCHI15 )THEN
               AUX15 = G(IPREV)*SMT15 + GNSQ(IPREV)
            ELSE
               AUX15 = -AUX15*OLDSMT15 - G(IPREV)*TAU(15)
            ENDIF
            IF( BRNCHI16 )THEN
               AUX16 = G(IPREV)*SMT16 + GNSQ(IPREV)
            ELSE
               AUX16 = -AUX16*OLDSMT16 - G(IPREV)*TAU(16)
            ENDIF
            AUX1 = (GNSQ(I) * AUX1) / (BDET(IPREV) * GPL1)
            AUX2 = (GNSQ(I) * AUX2) / (BDET(IPREV) * GPL2)
            AUX3 = (GNSQ(I) * AUX3) / (BDET(IPREV) * GPL3)
            AUX4 = (GNSQ(I) * AUX4) / (BDET(IPREV) * GPL4)
            AUX5 = (GNSQ(I) * AUX5) / (BDET(IPREV) * GPL5)
            AUX6 = (GNSQ(I) * AUX6) / (BDET(IPREV) * GPL6)
            AUX7 = (GNSQ(I) * AUX7) / (BDET(IPREV) * GPL7)
            AUX8 = (GNSQ(I) * AUX8) / (BDET(IPREV) * GPL8)
            AUX9 = (GNSQ(I) * AUX9) / (BDET(IPREV) * GPL9)
            AUX10 = (GNSQ(I) * AUX10) / (BDET(IPREV) * GPL10)
            AUX11 = (GNSQ(I) * AUX11) / (BDET(IPREV) * GPL11)
            AUX12 = (GNSQ(I) * AUX12) / (BDET(IPREV) * GPL12)
            AUX13 = (GNSQ(I) * AUX13) / (BDET(IPREV) * GPL13)
            AUX14 = (GNSQ(I) * AUX14) / (BDET(IPREV) * GPL14)
            AUX15 = (GNSQ(I) * AUX15) / (BDET(IPREV) * GPL15)
            AUX16 = (GNSQ(I) * AUX16) / (BDET(IPREV) * GPL16)
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
      ANEGC(9) = NEGC9
      ANEGC(10) = NEGC10
      ANEGC(11) = NEGC11
      ANEGC(12) = NEGC12
      ANEGC(13) = NEGC13
      ANEGC(14) = NEGC14
      ANEGC(15) = NEGC15
      ANEGC(16) = NEGC16
      AAUX(1) = AUX1
      AAUX(2) = AUX2
      AAUX(3) = AUX3
      AAUX(4) = AUX4
      AAUX(5) = AUX5
      AAUX(6) = AUX6
      AAUX(7) = AUX7
      AAUX(8) = AUX8
      AAUX(9) = AUX9
      AAUX(10) = AUX10
      AAUX(11) = AUX11
      AAUX(12) = AUX12
      AAUX(13) = AUX13
      AAUX(14) = AUX14
      AAUX(15) = AUX15
      AAUX(16) = AUX16
      END SUBROUTINE DLAXRM_STAT16
*
************************************************************************
