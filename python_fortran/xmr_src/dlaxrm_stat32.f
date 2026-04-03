      SUBROUTINE DLAXRM_STAT32(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, ANEGC, AAUX
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(32)
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  ANEGC(32)
      DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(32)
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
      DOUBLE PRECISION  AUX17, GPL17, SMT17, OLDSMT17, SBOUND17
      DOUBLE PRECISION  AUX18, GPL18, SMT18, OLDSMT18, SBOUND18
      DOUBLE PRECISION  AUX19, GPL19, SMT19, OLDSMT19, SBOUND19
      DOUBLE PRECISION  AUX20, GPL20, SMT20, OLDSMT20, SBOUND20
      DOUBLE PRECISION  AUX21, GPL21, SMT21, OLDSMT21, SBOUND21
      DOUBLE PRECISION  AUX22, GPL22, SMT22, OLDSMT22, SBOUND22
      DOUBLE PRECISION  AUX23, GPL23, SMT23, OLDSMT23, SBOUND23
      DOUBLE PRECISION  AUX24, GPL24, SMT24, OLDSMT24, SBOUND24
      DOUBLE PRECISION  AUX25, GPL25, SMT25, OLDSMT25, SBOUND25
      DOUBLE PRECISION  AUX26, GPL26, SMT26, OLDSMT26, SBOUND26
      DOUBLE PRECISION  AUX27, GPL27, SMT27, OLDSMT27, SBOUND27
      DOUBLE PRECISION  AUX28, GPL28, SMT28, OLDSMT28, SBOUND28
      DOUBLE PRECISION  AUX29, GPL29, SMT29, OLDSMT29, SBOUND29
      DOUBLE PRECISION  AUX30, GPL30, SMT30, OLDSMT30, SBOUND30
      DOUBLE PRECISION  AUX31, GPL31, SMT31, OLDSMT31, SBOUND31
      DOUBLE PRECISION  AUX32, GPL32, SMT32, OLDSMT32, SBOUND32
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
      INTEGER  NEGC17
      INTEGER  NEGC18
      INTEGER  NEGC19
      INTEGER  NEGC20
      INTEGER  NEGC21
      INTEGER  NEGC22
      INTEGER  NEGC23
      INTEGER  NEGC24
      INTEGER  NEGC25
      INTEGER  NEGC26
      INTEGER  NEGC27
      INTEGER  NEGC28
      INTEGER  NEGC29
      INTEGER  NEGC30
      INTEGER  NEGC31
      INTEGER  NEGC32
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
      LOGICAL  BRNCHI17
      LOGICAL  BRNCHI18
      LOGICAL  BRNCHI19
      LOGICAL  BRNCHI20
      LOGICAL  BRNCHI21
      LOGICAL  BRNCHI22
      LOGICAL  BRNCHI23
      LOGICAL  BRNCHI24
      LOGICAL  BRNCHI25
      LOGICAL  BRNCHI26
      LOGICAL  BRNCHI27
      LOGICAL  BRNCHI28
      LOGICAL  BRNCHI29
      LOGICAL  BRNCHI30
      LOGICAL  BRNCHI31
      LOGICAL  BRNCHI32
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
      SBOUND17 = MIN( PRBRK, N ) * ABS(TAU(17))
      SBOUND18 = MIN( PRBRK, N ) * ABS(TAU(18))
      SBOUND19 = MIN( PRBRK, N ) * ABS(TAU(19))
      SBOUND20 = MIN( PRBRK, N ) * ABS(TAU(20))
      SBOUND21 = MIN( PRBRK, N ) * ABS(TAU(21))
      SBOUND22 = MIN( PRBRK, N ) * ABS(TAU(22))
      SBOUND23 = MIN( PRBRK, N ) * ABS(TAU(23))
      SBOUND24 = MIN( PRBRK, N ) * ABS(TAU(24))
      SBOUND25 = MIN( PRBRK, N ) * ABS(TAU(25))
      SBOUND26 = MIN( PRBRK, N ) * ABS(TAU(26))
      SBOUND27 = MIN( PRBRK, N ) * ABS(TAU(27))
      SBOUND28 = MIN( PRBRK, N ) * ABS(TAU(28))
      SBOUND29 = MIN( PRBRK, N ) * ABS(TAU(29))
      SBOUND30 = MIN( PRBRK, N ) * ABS(TAU(30))
      SBOUND31 = MIN( PRBRK, N ) * ABS(TAU(31))
      SBOUND32 = MIN( PRBRK, N ) * ABS(TAU(32))
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
      NEGC17 = 0
      NEGC18 = 0
      NEGC19 = 0
      NEGC20 = 0
      NEGC21 = 0
      NEGC22 = 0
      NEGC23 = 0
      NEGC24 = 0
      NEGC25 = 0
      NEGC26 = 0
      NEGC27 = 0
      NEGC28 = 0
      NEGC29 = 0
      NEGC30 = 0
      NEGC31 = 0
      NEGC32 = 0
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
      AUX17 = ZERO
      AUX18 = ZERO
      AUX19 = ZERO
      AUX20 = ZERO
      AUX21 = ZERO
      AUX22 = ZERO
      AUX23 = ZERO
      AUX24 = ZERO
      AUX25 = ZERO
      AUX26 = ZERO
      AUX27 = ZERO
      AUX28 = ZERO
      AUX29 = ZERO
      AUX30 = ZERO
      AUX31 = ZERO
      AUX32 = ZERO
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
            AUX17 = AUX17 - TAU(17)
            AUX18 = AUX18 - TAU(18)
            AUX19 = AUX19 - TAU(19)
            AUX20 = AUX20 - TAU(20)
            AUX21 = AUX21 - TAU(21)
            AUX22 = AUX22 - TAU(22)
            AUX23 = AUX23 - TAU(23)
            AUX24 = AUX24 - TAU(24)
            AUX25 = AUX25 - TAU(25)
            AUX26 = AUX26 - TAU(26)
            AUX27 = AUX27 - TAU(27)
            AUX28 = AUX28 - TAU(28)
            AUX29 = AUX29 - TAU(29)
            AUX30 = AUX30 - TAU(30)
            AUX31 = AUX31 - TAU(31)
            AUX32 = AUX32 - TAU(32)
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
            GPL17 = G(I) + AUX17
            GPL18 = G(I) + AUX18
            GPL19 = G(I) + AUX19
            GPL20 = G(I) + AUX20
            GPL21 = G(I) + AUX21
            GPL22 = G(I) + AUX22
            GPL23 = G(I) + AUX23
            GPL24 = G(I) + AUX24
            GPL25 = G(I) + AUX25
            GPL26 = G(I) + AUX26
            GPL27 = G(I) + AUX27
            GPL28 = G(I) + AUX28
            GPL29 = G(I) + AUX29
            GPL30 = G(I) + AUX30
            GPL31 = G(I) + AUX31
            GPL32 = G(I) + AUX32
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
            IF( ABS(GPL17) .LT. PIVMIN )  GPL17 = NEGPIV
            IF( ABS(GPL18) .LT. PIVMIN )  GPL18 = NEGPIV
            IF( ABS(GPL19) .LT. PIVMIN )  GPL19 = NEGPIV
            IF( ABS(GPL20) .LT. PIVMIN )  GPL20 = NEGPIV
            IF( ABS(GPL21) .LT. PIVMIN )  GPL21 = NEGPIV
            IF( ABS(GPL22) .LT. PIVMIN )  GPL22 = NEGPIV
            IF( ABS(GPL23) .LT. PIVMIN )  GPL23 = NEGPIV
            IF( ABS(GPL24) .LT. PIVMIN )  GPL24 = NEGPIV
            IF( ABS(GPL25) .LT. PIVMIN )  GPL25 = NEGPIV
            IF( ABS(GPL26) .LT. PIVMIN )  GPL26 = NEGPIV
            IF( ABS(GPL27) .LT. PIVMIN )  GPL27 = NEGPIV
            IF( ABS(GPL28) .LT. PIVMIN )  GPL28 = NEGPIV
            IF( ABS(GPL29) .LT. PIVMIN )  GPL29 = NEGPIV
            IF( ABS(GPL30) .LT. PIVMIN )  GPL30 = NEGPIV
            IF( ABS(GPL31) .LT. PIVMIN )  GPL31 = NEGPIV
            IF( ABS(GPL32) .LT. PIVMIN )  GPL32 = NEGPIV
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
            IF( GPL17 .LT. ZERO )  NEGC17 = NEGC17 + 1
            IF( GPL18 .LT. ZERO )  NEGC18 = NEGC18 + 1
            IF( GPL19 .LT. ZERO )  NEGC19 = NEGC19 + 1
            IF( GPL20 .LT. ZERO )  NEGC20 = NEGC20 + 1
            IF( GPL21 .LT. ZERO )  NEGC21 = NEGC21 + 1
            IF( GPL22 .LT. ZERO )  NEGC22 = NEGC22 + 1
            IF( GPL23 .LT. ZERO )  NEGC23 = NEGC23 + 1
            IF( GPL24 .LT. ZERO )  NEGC24 = NEGC24 + 1
            IF( GPL25 .LT. ZERO )  NEGC25 = NEGC25 + 1
            IF( GPL26 .LT. ZERO )  NEGC26 = NEGC26 + 1
            IF( GPL27 .LT. ZERO )  NEGC27 = NEGC27 + 1
            IF( GPL28 .LT. ZERO )  NEGC28 = NEGC28 + 1
            IF( GPL29 .LT. ZERO )  NEGC29 = NEGC29 + 1
            IF( GPL30 .LT. ZERO )  NEGC30 = NEGC30 + 1
            IF( GPL31 .LT. ZERO )  NEGC31 = NEGC31 + 1
            IF( GPL32 .LT. ZERO )  NEGC32 = NEGC32 + 1
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
            AUX17 = NGN(I) * ( AUX17 / GPL17 )
            AUX18 = NGN(I) * ( AUX18 / GPL18 )
            AUX19 = NGN(I) * ( AUX19 / GPL19 )
            AUX20 = NGN(I) * ( AUX20 / GPL20 )
            AUX21 = NGN(I) * ( AUX21 / GPL21 )
            AUX22 = NGN(I) * ( AUX22 / GPL22 )
            AUX23 = NGN(I) * ( AUX23 / GPL23 )
            AUX24 = NGN(I) * ( AUX24 / GPL24 )
            AUX25 = NGN(I) * ( AUX25 / GPL25 )
            AUX26 = NGN(I) * ( AUX26 / GPL26 )
            AUX27 = NGN(I) * ( AUX27 / GPL27 )
            AUX28 = NGN(I) * ( AUX28 / GPL28 )
            AUX29 = NGN(I) * ( AUX29 / GPL29 )
            AUX30 = NGN(I) * ( AUX30 / GPL30 )
            AUX31 = NGN(I) * ( AUX31 / GPL31 )
            AUX32 = NGN(I) * ( AUX32 / GPL32 )
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
         SMT17 = AUX17 - TAU(17)
         SMT18 = AUX18 - TAU(18)
         SMT19 = AUX19 - TAU(19)
         SMT20 = AUX20 - TAU(20)
         SMT21 = AUX21 - TAU(21)
         SMT22 = AUX22 - TAU(22)
         SMT23 = AUX23 - TAU(23)
         SMT24 = AUX24 - TAU(24)
         SMT25 = AUX25 - TAU(25)
         SMT26 = AUX26 - TAU(26)
         SMT27 = AUX27 - TAU(27)
         SMT28 = AUX28 - TAU(28)
         SMT29 = AUX29 - TAU(29)
         SMT30 = AUX30 - TAU(30)
         SMT31 = AUX31 - TAU(31)
         SMT32 = AUX32 - TAU(32)
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
         GPL17 = G(I) + SMT17
         GPL18 = G(I) + SMT18
         GPL19 = G(I) + SMT19
         GPL20 = G(I) + SMT20
         GPL21 = G(I) + SMT21
         GPL22 = G(I) + SMT22
         GPL23 = G(I) + SMT23
         GPL24 = G(I) + SMT24
         GPL25 = G(I) + SMT25
         GPL26 = G(I) + SMT26
         GPL27 = G(I) + SMT27
         GPL28 = G(I) + SMT28
         GPL29 = G(I) + SMT29
         GPL30 = G(I) + SMT30
         GPL31 = G(I) + SMT31
         GPL32 = G(I) + SMT32
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
         IF( ABS(GPL17) .LT. PIVMIN )  GPL17 = NEGPIV
         IF( ABS(GPL18) .LT. PIVMIN )  GPL18 = NEGPIV
         IF( ABS(GPL19) .LT. PIVMIN )  GPL19 = NEGPIV
         IF( ABS(GPL20) .LT. PIVMIN )  GPL20 = NEGPIV
         IF( ABS(GPL21) .LT. PIVMIN )  GPL21 = NEGPIV
         IF( ABS(GPL22) .LT. PIVMIN )  GPL22 = NEGPIV
         IF( ABS(GPL23) .LT. PIVMIN )  GPL23 = NEGPIV
         IF( ABS(GPL24) .LT. PIVMIN )  GPL24 = NEGPIV
         IF( ABS(GPL25) .LT. PIVMIN )  GPL25 = NEGPIV
         IF( ABS(GPL26) .LT. PIVMIN )  GPL26 = NEGPIV
         IF( ABS(GPL27) .LT. PIVMIN )  GPL27 = NEGPIV
         IF( ABS(GPL28) .LT. PIVMIN )  GPL28 = NEGPIV
         IF( ABS(GPL29) .LT. PIVMIN )  GPL29 = NEGPIV
         IF( ABS(GPL30) .LT. PIVMIN )  GPL30 = NEGPIV
         IF( ABS(GPL31) .LT. PIVMIN )  GPL31 = NEGPIV
         IF( ABS(GPL32) .LT. PIVMIN )  GPL32 = NEGPIV
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
         IF( GPL17 .LT. ZERO )  NEGC17 = NEGC17 + 1
         IF( GPL18 .LT. ZERO )  NEGC18 = NEGC18 + 1
         IF( GPL19 .LT. ZERO )  NEGC19 = NEGC19 + 1
         IF( GPL20 .LT. ZERO )  NEGC20 = NEGC20 + 1
         IF( GPL21 .LT. ZERO )  NEGC21 = NEGC21 + 1
         IF( GPL22 .LT. ZERO )  NEGC22 = NEGC22 + 1
         IF( GPL23 .LT. ZERO )  NEGC23 = NEGC23 + 1
         IF( GPL24 .LT. ZERO )  NEGC24 = NEGC24 + 1
         IF( GPL25 .LT. ZERO )  NEGC25 = NEGC25 + 1
         IF( GPL26 .LT. ZERO )  NEGC26 = NEGC26 + 1
         IF( GPL27 .LT. ZERO )  NEGC27 = NEGC27 + 1
         IF( GPL28 .LT. ZERO )  NEGC28 = NEGC28 + 1
         IF( GPL29 .LT. ZERO )  NEGC29 = NEGC29 + 1
         IF( GPL30 .LT. ZERO )  NEGC30 = NEGC30 + 1
         IF( GPL31 .LT. ZERO )  NEGC31 = NEGC31 + 1
         IF( GPL32 .LT. ZERO )  NEGC32 = NEGC32 + 1
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
         AUX17 = - GNSQ(I) / GPL17
         AUX18 = - GNSQ(I) / GPL18
         AUX19 = - GNSQ(I) / GPL19
         AUX20 = - GNSQ(I) / GPL20
         AUX21 = - GNSQ(I) / GPL21
         AUX22 = - GNSQ(I) / GPL22
         AUX23 = - GNSQ(I) / GPL23
         AUX24 = - GNSQ(I) / GPL24
         AUX25 = - GNSQ(I) / GPL25
         AUX26 = - GNSQ(I) / GPL26
         AUX27 = - GNSQ(I) / GPL27
         AUX28 = - GNSQ(I) / GPL28
         AUX29 = - GNSQ(I) / GPL29
         AUX30 = - GNSQ(I) / GPL30
         AUX31 = - GNSQ(I) / GPL31
         AUX32 = - GNSQ(I) / GPL32
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
         OLDSMT17 = SMT17
         BRNCHI17 = ( ABS(AUX17).LE.SBOUND17 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL17) )
         OLDSMT18 = SMT18
         BRNCHI18 = ( ABS(AUX18).LE.SBOUND18 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL18) )
         OLDSMT19 = SMT19
         BRNCHI19 = ( ABS(AUX19).LE.SBOUND19 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL19) )
         OLDSMT20 = SMT20
         BRNCHI20 = ( ABS(AUX20).LE.SBOUND20 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL20) )
         OLDSMT21 = SMT21
         BRNCHI21 = ( ABS(AUX21).LE.SBOUND21 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL21) )
         OLDSMT22 = SMT22
         BRNCHI22 = ( ABS(AUX22).LE.SBOUND22 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL22) )
         OLDSMT23 = SMT23
         BRNCHI23 = ( ABS(AUX23).LE.SBOUND23 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL23) )
         OLDSMT24 = SMT24
         BRNCHI24 = ( ABS(AUX24).LE.SBOUND24 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL24) )
         OLDSMT25 = SMT25
         BRNCHI25 = ( ABS(AUX25).LE.SBOUND25 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL25) )
         OLDSMT26 = SMT26
         BRNCHI26 = ( ABS(AUX26).LE.SBOUND26 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL26) )
         OLDSMT27 = SMT27
         BRNCHI27 = ( ABS(AUX27).LE.SBOUND27 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL27) )
         OLDSMT28 = SMT28
         BRNCHI28 = ( ABS(AUX28).LE.SBOUND28 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL28) )
         OLDSMT29 = SMT29
         BRNCHI29 = ( ABS(AUX29).LE.SBOUND29 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL29) )
         OLDSMT30 = SMT30
         BRNCHI30 = ( ABS(AUX30).LE.SBOUND30 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL30) )
         OLDSMT31 = SMT31
         BRNCHI31 = ( ABS(AUX31).LE.SBOUND31 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL31) )
         OLDSMT32 = SMT32
         BRNCHI32 = ( ABS(AUX32).LE.SBOUND32 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL32) )
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
         SMT17 = AUX17 - TAU(17)
         SMT18 = AUX18 - TAU(18)
         SMT19 = AUX19 - TAU(19)
         SMT20 = AUX20 - TAU(20)
         SMT21 = AUX21 - TAU(21)
         SMT22 = AUX22 - TAU(22)
         SMT23 = AUX23 - TAU(23)
         SMT24 = AUX24 - TAU(24)
         SMT25 = AUX25 - TAU(25)
         SMT26 = AUX26 - TAU(26)
         SMT27 = AUX27 - TAU(27)
         SMT28 = AUX28 - TAU(28)
         SMT29 = AUX29 - TAU(29)
         SMT30 = AUX30 - TAU(30)
         SMT31 = AUX31 - TAU(31)
         SMT32 = AUX32 - TAU(32)
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
         GPL17 = G(I) + SMT17
         GPL18 = G(I) + SMT18
         GPL19 = G(I) + SMT19
         GPL20 = G(I) + SMT20
         GPL21 = G(I) + SMT21
         GPL22 = G(I) + SMT22
         GPL23 = G(I) + SMT23
         GPL24 = G(I) + SMT24
         GPL25 = G(I) + SMT25
         GPL26 = G(I) + SMT26
         GPL27 = G(I) + SMT27
         GPL28 = G(I) + SMT28
         GPL29 = G(I) + SMT29
         GPL30 = G(I) + SMT30
         GPL31 = G(I) + SMT31
         GPL32 = G(I) + SMT32
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
         IF( ABS(GPL17) .LT. PIVMIN )  GPL17 = NEGPIV
         IF( ABS(GPL18) .LT. PIVMIN )  GPL18 = NEGPIV
         IF( ABS(GPL19) .LT. PIVMIN )  GPL19 = NEGPIV
         IF( ABS(GPL20) .LT. PIVMIN )  GPL20 = NEGPIV
         IF( ABS(GPL21) .LT. PIVMIN )  GPL21 = NEGPIV
         IF( ABS(GPL22) .LT. PIVMIN )  GPL22 = NEGPIV
         IF( ABS(GPL23) .LT. PIVMIN )  GPL23 = NEGPIV
         IF( ABS(GPL24) .LT. PIVMIN )  GPL24 = NEGPIV
         IF( ABS(GPL25) .LT. PIVMIN )  GPL25 = NEGPIV
         IF( ABS(GPL26) .LT. PIVMIN )  GPL26 = NEGPIV
         IF( ABS(GPL27) .LT. PIVMIN )  GPL27 = NEGPIV
         IF( ABS(GPL28) .LT. PIVMIN )  GPL28 = NEGPIV
         IF( ABS(GPL29) .LT. PIVMIN )  GPL29 = NEGPIV
         IF( ABS(GPL30) .LT. PIVMIN )  GPL30 = NEGPIV
         IF( ABS(GPL31) .LT. PIVMIN )  GPL31 = NEGPIV
         IF( ABS(GPL32) .LT. PIVMIN )  GPL32 = NEGPIV
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
         IF( GPL17 .LT. ZERO )  NEGC17 = NEGC17 + 1
         IF( GPL18 .LT. ZERO )  NEGC18 = NEGC18 + 1
         IF( GPL19 .LT. ZERO )  NEGC19 = NEGC19 + 1
         IF( GPL20 .LT. ZERO )  NEGC20 = NEGC20 + 1
         IF( GPL21 .LT. ZERO )  NEGC21 = NEGC21 + 1
         IF( GPL22 .LT. ZERO )  NEGC22 = NEGC22 + 1
         IF( GPL23 .LT. ZERO )  NEGC23 = NEGC23 + 1
         IF( GPL24 .LT. ZERO )  NEGC24 = NEGC24 + 1
         IF( GPL25 .LT. ZERO )  NEGC25 = NEGC25 + 1
         IF( GPL26 .LT. ZERO )  NEGC26 = NEGC26 + 1
         IF( GPL27 .LT. ZERO )  NEGC27 = NEGC27 + 1
         IF( GPL28 .LT. ZERO )  NEGC28 = NEGC28 + 1
         IF( GPL29 .LT. ZERO )  NEGC29 = NEGC29 + 1
         IF( GPL30 .LT. ZERO )  NEGC30 = NEGC30 + 1
         IF( GPL31 .LT. ZERO )  NEGC31 = NEGC31 + 1
         IF( GPL32 .LT. ZERO )  NEGC32 = NEGC32 + 1
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
            AUX17 = - GNSQ(I) / GPL17
            AUX18 = - GNSQ(I) / GPL18
            AUX19 = - GNSQ(I) / GPL19
            AUX20 = - GNSQ(I) / GPL20
            AUX21 = - GNSQ(I) / GPL21
            AUX22 = - GNSQ(I) / GPL22
            AUX23 = - GNSQ(I) / GPL23
            AUX24 = - GNSQ(I) / GPL24
            AUX25 = - GNSQ(I) / GPL25
            AUX26 = - GNSQ(I) / GPL26
            AUX27 = - GNSQ(I) / GPL27
            AUX28 = - GNSQ(I) / GPL28
            AUX29 = - GNSQ(I) / GPL29
            AUX30 = - GNSQ(I) / GPL30
            AUX31 = - GNSQ(I) / GPL31
            AUX32 = - GNSQ(I) / GPL32
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
            IF( BRNCHI17 )THEN
               AUX17 = G(IPREV)*SMT17 + GNSQ(IPREV)
            ELSE
               AUX17 = -AUX17*OLDSMT17 - G(IPREV)*TAU(17)
            ENDIF
            IF( BRNCHI18 )THEN
               AUX18 = G(IPREV)*SMT18 + GNSQ(IPREV)
            ELSE
               AUX18 = -AUX18*OLDSMT18 - G(IPREV)*TAU(18)
            ENDIF
            IF( BRNCHI19 )THEN
               AUX19 = G(IPREV)*SMT19 + GNSQ(IPREV)
            ELSE
               AUX19 = -AUX19*OLDSMT19 - G(IPREV)*TAU(19)
            ENDIF
            IF( BRNCHI20 )THEN
               AUX20 = G(IPREV)*SMT20 + GNSQ(IPREV)
            ELSE
               AUX20 = -AUX20*OLDSMT20 - G(IPREV)*TAU(20)
            ENDIF
            IF( BRNCHI21 )THEN
               AUX21 = G(IPREV)*SMT21 + GNSQ(IPREV)
            ELSE
               AUX21 = -AUX21*OLDSMT21 - G(IPREV)*TAU(21)
            ENDIF
            IF( BRNCHI22 )THEN
               AUX22 = G(IPREV)*SMT22 + GNSQ(IPREV)
            ELSE
               AUX22 = -AUX22*OLDSMT22 - G(IPREV)*TAU(22)
            ENDIF
            IF( BRNCHI23 )THEN
               AUX23 = G(IPREV)*SMT23 + GNSQ(IPREV)
            ELSE
               AUX23 = -AUX23*OLDSMT23 - G(IPREV)*TAU(23)
            ENDIF
            IF( BRNCHI24 )THEN
               AUX24 = G(IPREV)*SMT24 + GNSQ(IPREV)
            ELSE
               AUX24 = -AUX24*OLDSMT24 - G(IPREV)*TAU(24)
            ENDIF
            IF( BRNCHI25 )THEN
               AUX25 = G(IPREV)*SMT25 + GNSQ(IPREV)
            ELSE
               AUX25 = -AUX25*OLDSMT25 - G(IPREV)*TAU(25)
            ENDIF
            IF( BRNCHI26 )THEN
               AUX26 = G(IPREV)*SMT26 + GNSQ(IPREV)
            ELSE
               AUX26 = -AUX26*OLDSMT26 - G(IPREV)*TAU(26)
            ENDIF
            IF( BRNCHI27 )THEN
               AUX27 = G(IPREV)*SMT27 + GNSQ(IPREV)
            ELSE
               AUX27 = -AUX27*OLDSMT27 - G(IPREV)*TAU(27)
            ENDIF
            IF( BRNCHI28 )THEN
               AUX28 = G(IPREV)*SMT28 + GNSQ(IPREV)
            ELSE
               AUX28 = -AUX28*OLDSMT28 - G(IPREV)*TAU(28)
            ENDIF
            IF( BRNCHI29 )THEN
               AUX29 = G(IPREV)*SMT29 + GNSQ(IPREV)
            ELSE
               AUX29 = -AUX29*OLDSMT29 - G(IPREV)*TAU(29)
            ENDIF
            IF( BRNCHI30 )THEN
               AUX30 = G(IPREV)*SMT30 + GNSQ(IPREV)
            ELSE
               AUX30 = -AUX30*OLDSMT30 - G(IPREV)*TAU(30)
            ENDIF
            IF( BRNCHI31 )THEN
               AUX31 = G(IPREV)*SMT31 + GNSQ(IPREV)
            ELSE
               AUX31 = -AUX31*OLDSMT31 - G(IPREV)*TAU(31)
            ENDIF
            IF( BRNCHI32 )THEN
               AUX32 = G(IPREV)*SMT32 + GNSQ(IPREV)
            ELSE
               AUX32 = -AUX32*OLDSMT32 - G(IPREV)*TAU(32)
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
            AUX17 = (GNSQ(I) * AUX17) / (BDET(IPREV) * GPL17)
            AUX18 = (GNSQ(I) * AUX18) / (BDET(IPREV) * GPL18)
            AUX19 = (GNSQ(I) * AUX19) / (BDET(IPREV) * GPL19)
            AUX20 = (GNSQ(I) * AUX20) / (BDET(IPREV) * GPL20)
            AUX21 = (GNSQ(I) * AUX21) / (BDET(IPREV) * GPL21)
            AUX22 = (GNSQ(I) * AUX22) / (BDET(IPREV) * GPL22)
            AUX23 = (GNSQ(I) * AUX23) / (BDET(IPREV) * GPL23)
            AUX24 = (GNSQ(I) * AUX24) / (BDET(IPREV) * GPL24)
            AUX25 = (GNSQ(I) * AUX25) / (BDET(IPREV) * GPL25)
            AUX26 = (GNSQ(I) * AUX26) / (BDET(IPREV) * GPL26)
            AUX27 = (GNSQ(I) * AUX27) / (BDET(IPREV) * GPL27)
            AUX28 = (GNSQ(I) * AUX28) / (BDET(IPREV) * GPL28)
            AUX29 = (GNSQ(I) * AUX29) / (BDET(IPREV) * GPL29)
            AUX30 = (GNSQ(I) * AUX30) / (BDET(IPREV) * GPL30)
            AUX31 = (GNSQ(I) * AUX31) / (BDET(IPREV) * GPL31)
            AUX32 = (GNSQ(I) * AUX32) / (BDET(IPREV) * GPL32)
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
      ANEGC(17) = NEGC17
      ANEGC(18) = NEGC18
      ANEGC(19) = NEGC19
      ANEGC(20) = NEGC20
      ANEGC(21) = NEGC21
      ANEGC(22) = NEGC22
      ANEGC(23) = NEGC23
      ANEGC(24) = NEGC24
      ANEGC(25) = NEGC25
      ANEGC(26) = NEGC26
      ANEGC(27) = NEGC27
      ANEGC(28) = NEGC28
      ANEGC(29) = NEGC29
      ANEGC(30) = NEGC30
      ANEGC(31) = NEGC31
      ANEGC(32) = NEGC32
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
      AAUX(17) = AUX17
      AAUX(18) = AUX18
      AAUX(19) = AUX19
      AAUX(20) = AUX20
      AAUX(21) = AUX21
      AAUX(22) = AUX22
      AAUX(23) = AUX23
      AAUX(24) = AUX24
      AAUX(25) = AUX25
      AAUX(26) = AUX26
      AAUX(27) = AUX27
      AAUX(28) = AUX28
      AAUX(29) = AUX29
      AAUX(30) = AUX30
      AAUX(31) = AUX31
      AAUX(32) = AUX32
      END SUBROUTINE DLAXRM_STAT32
*
************************************************************************
