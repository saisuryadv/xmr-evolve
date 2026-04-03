      SUBROUTINE DLAXRM_STAT64(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, ANEGC, AAUX
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(64)
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  ANEGC(64)
      DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(64)
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
      DOUBLE PRECISION  AUX33, GPL33, SMT33, OLDSMT33, SBOUND33
      DOUBLE PRECISION  AUX34, GPL34, SMT34, OLDSMT34, SBOUND34
      DOUBLE PRECISION  AUX35, GPL35, SMT35, OLDSMT35, SBOUND35
      DOUBLE PRECISION  AUX36, GPL36, SMT36, OLDSMT36, SBOUND36
      DOUBLE PRECISION  AUX37, GPL37, SMT37, OLDSMT37, SBOUND37
      DOUBLE PRECISION  AUX38, GPL38, SMT38, OLDSMT38, SBOUND38
      DOUBLE PRECISION  AUX39, GPL39, SMT39, OLDSMT39, SBOUND39
      DOUBLE PRECISION  AUX40, GPL40, SMT40, OLDSMT40, SBOUND40
      DOUBLE PRECISION  AUX41, GPL41, SMT41, OLDSMT41, SBOUND41
      DOUBLE PRECISION  AUX42, GPL42, SMT42, OLDSMT42, SBOUND42
      DOUBLE PRECISION  AUX43, GPL43, SMT43, OLDSMT43, SBOUND43
      DOUBLE PRECISION  AUX44, GPL44, SMT44, OLDSMT44, SBOUND44
      DOUBLE PRECISION  AUX45, GPL45, SMT45, OLDSMT45, SBOUND45
      DOUBLE PRECISION  AUX46, GPL46, SMT46, OLDSMT46, SBOUND46
      DOUBLE PRECISION  AUX47, GPL47, SMT47, OLDSMT47, SBOUND47
      DOUBLE PRECISION  AUX48, GPL48, SMT48, OLDSMT48, SBOUND48
      DOUBLE PRECISION  AUX49, GPL49, SMT49, OLDSMT49, SBOUND49
      DOUBLE PRECISION  AUX50, GPL50, SMT50, OLDSMT50, SBOUND50
      DOUBLE PRECISION  AUX51, GPL51, SMT51, OLDSMT51, SBOUND51
      DOUBLE PRECISION  AUX52, GPL52, SMT52, OLDSMT52, SBOUND52
      DOUBLE PRECISION  AUX53, GPL53, SMT53, OLDSMT53, SBOUND53
      DOUBLE PRECISION  AUX54, GPL54, SMT54, OLDSMT54, SBOUND54
      DOUBLE PRECISION  AUX55, GPL55, SMT55, OLDSMT55, SBOUND55
      DOUBLE PRECISION  AUX56, GPL56, SMT56, OLDSMT56, SBOUND56
      DOUBLE PRECISION  AUX57, GPL57, SMT57, OLDSMT57, SBOUND57
      DOUBLE PRECISION  AUX58, GPL58, SMT58, OLDSMT58, SBOUND58
      DOUBLE PRECISION  AUX59, GPL59, SMT59, OLDSMT59, SBOUND59
      DOUBLE PRECISION  AUX60, GPL60, SMT60, OLDSMT60, SBOUND60
      DOUBLE PRECISION  AUX61, GPL61, SMT61, OLDSMT61, SBOUND61
      DOUBLE PRECISION  AUX62, GPL62, SMT62, OLDSMT62, SBOUND62
      DOUBLE PRECISION  AUX63, GPL63, SMT63, OLDSMT63, SBOUND63
      DOUBLE PRECISION  AUX64, GPL64, SMT64, OLDSMT64, SBOUND64
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
      INTEGER  NEGC33
      INTEGER  NEGC34
      INTEGER  NEGC35
      INTEGER  NEGC36
      INTEGER  NEGC37
      INTEGER  NEGC38
      INTEGER  NEGC39
      INTEGER  NEGC40
      INTEGER  NEGC41
      INTEGER  NEGC42
      INTEGER  NEGC43
      INTEGER  NEGC44
      INTEGER  NEGC45
      INTEGER  NEGC46
      INTEGER  NEGC47
      INTEGER  NEGC48
      INTEGER  NEGC49
      INTEGER  NEGC50
      INTEGER  NEGC51
      INTEGER  NEGC52
      INTEGER  NEGC53
      INTEGER  NEGC54
      INTEGER  NEGC55
      INTEGER  NEGC56
      INTEGER  NEGC57
      INTEGER  NEGC58
      INTEGER  NEGC59
      INTEGER  NEGC60
      INTEGER  NEGC61
      INTEGER  NEGC62
      INTEGER  NEGC63
      INTEGER  NEGC64
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
      LOGICAL  BRNCHI33
      LOGICAL  BRNCHI34
      LOGICAL  BRNCHI35
      LOGICAL  BRNCHI36
      LOGICAL  BRNCHI37
      LOGICAL  BRNCHI38
      LOGICAL  BRNCHI39
      LOGICAL  BRNCHI40
      LOGICAL  BRNCHI41
      LOGICAL  BRNCHI42
      LOGICAL  BRNCHI43
      LOGICAL  BRNCHI44
      LOGICAL  BRNCHI45
      LOGICAL  BRNCHI46
      LOGICAL  BRNCHI47
      LOGICAL  BRNCHI48
      LOGICAL  BRNCHI49
      LOGICAL  BRNCHI50
      LOGICAL  BRNCHI51
      LOGICAL  BRNCHI52
      LOGICAL  BRNCHI53
      LOGICAL  BRNCHI54
      LOGICAL  BRNCHI55
      LOGICAL  BRNCHI56
      LOGICAL  BRNCHI57
      LOGICAL  BRNCHI58
      LOGICAL  BRNCHI59
      LOGICAL  BRNCHI60
      LOGICAL  BRNCHI61
      LOGICAL  BRNCHI62
      LOGICAL  BRNCHI63
      LOGICAL  BRNCHI64
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
      SBOUND33 = MIN( PRBRK, N ) * ABS(TAU(33))
      SBOUND34 = MIN( PRBRK, N ) * ABS(TAU(34))
      SBOUND35 = MIN( PRBRK, N ) * ABS(TAU(35))
      SBOUND36 = MIN( PRBRK, N ) * ABS(TAU(36))
      SBOUND37 = MIN( PRBRK, N ) * ABS(TAU(37))
      SBOUND38 = MIN( PRBRK, N ) * ABS(TAU(38))
      SBOUND39 = MIN( PRBRK, N ) * ABS(TAU(39))
      SBOUND40 = MIN( PRBRK, N ) * ABS(TAU(40))
      SBOUND41 = MIN( PRBRK, N ) * ABS(TAU(41))
      SBOUND42 = MIN( PRBRK, N ) * ABS(TAU(42))
      SBOUND43 = MIN( PRBRK, N ) * ABS(TAU(43))
      SBOUND44 = MIN( PRBRK, N ) * ABS(TAU(44))
      SBOUND45 = MIN( PRBRK, N ) * ABS(TAU(45))
      SBOUND46 = MIN( PRBRK, N ) * ABS(TAU(46))
      SBOUND47 = MIN( PRBRK, N ) * ABS(TAU(47))
      SBOUND48 = MIN( PRBRK, N ) * ABS(TAU(48))
      SBOUND49 = MIN( PRBRK, N ) * ABS(TAU(49))
      SBOUND50 = MIN( PRBRK, N ) * ABS(TAU(50))
      SBOUND51 = MIN( PRBRK, N ) * ABS(TAU(51))
      SBOUND52 = MIN( PRBRK, N ) * ABS(TAU(52))
      SBOUND53 = MIN( PRBRK, N ) * ABS(TAU(53))
      SBOUND54 = MIN( PRBRK, N ) * ABS(TAU(54))
      SBOUND55 = MIN( PRBRK, N ) * ABS(TAU(55))
      SBOUND56 = MIN( PRBRK, N ) * ABS(TAU(56))
      SBOUND57 = MIN( PRBRK, N ) * ABS(TAU(57))
      SBOUND58 = MIN( PRBRK, N ) * ABS(TAU(58))
      SBOUND59 = MIN( PRBRK, N ) * ABS(TAU(59))
      SBOUND60 = MIN( PRBRK, N ) * ABS(TAU(60))
      SBOUND61 = MIN( PRBRK, N ) * ABS(TAU(61))
      SBOUND62 = MIN( PRBRK, N ) * ABS(TAU(62))
      SBOUND63 = MIN( PRBRK, N ) * ABS(TAU(63))
      SBOUND64 = MIN( PRBRK, N ) * ABS(TAU(64))
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
      NEGC33 = 0
      NEGC34 = 0
      NEGC35 = 0
      NEGC36 = 0
      NEGC37 = 0
      NEGC38 = 0
      NEGC39 = 0
      NEGC40 = 0
      NEGC41 = 0
      NEGC42 = 0
      NEGC43 = 0
      NEGC44 = 0
      NEGC45 = 0
      NEGC46 = 0
      NEGC47 = 0
      NEGC48 = 0
      NEGC49 = 0
      NEGC50 = 0
      NEGC51 = 0
      NEGC52 = 0
      NEGC53 = 0
      NEGC54 = 0
      NEGC55 = 0
      NEGC56 = 0
      NEGC57 = 0
      NEGC58 = 0
      NEGC59 = 0
      NEGC60 = 0
      NEGC61 = 0
      NEGC62 = 0
      NEGC63 = 0
      NEGC64 = 0
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
      AUX33 = ZERO
      AUX34 = ZERO
      AUX35 = ZERO
      AUX36 = ZERO
      AUX37 = ZERO
      AUX38 = ZERO
      AUX39 = ZERO
      AUX40 = ZERO
      AUX41 = ZERO
      AUX42 = ZERO
      AUX43 = ZERO
      AUX44 = ZERO
      AUX45 = ZERO
      AUX46 = ZERO
      AUX47 = ZERO
      AUX48 = ZERO
      AUX49 = ZERO
      AUX50 = ZERO
      AUX51 = ZERO
      AUX52 = ZERO
      AUX53 = ZERO
      AUX54 = ZERO
      AUX55 = ZERO
      AUX56 = ZERO
      AUX57 = ZERO
      AUX58 = ZERO
      AUX59 = ZERO
      AUX60 = ZERO
      AUX61 = ZERO
      AUX62 = ZERO
      AUX63 = ZERO
      AUX64 = ZERO
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
            AUX33 = AUX33 - TAU(33)
            AUX34 = AUX34 - TAU(34)
            AUX35 = AUX35 - TAU(35)
            AUX36 = AUX36 - TAU(36)
            AUX37 = AUX37 - TAU(37)
            AUX38 = AUX38 - TAU(38)
            AUX39 = AUX39 - TAU(39)
            AUX40 = AUX40 - TAU(40)
            AUX41 = AUX41 - TAU(41)
            AUX42 = AUX42 - TAU(42)
            AUX43 = AUX43 - TAU(43)
            AUX44 = AUX44 - TAU(44)
            AUX45 = AUX45 - TAU(45)
            AUX46 = AUX46 - TAU(46)
            AUX47 = AUX47 - TAU(47)
            AUX48 = AUX48 - TAU(48)
            AUX49 = AUX49 - TAU(49)
            AUX50 = AUX50 - TAU(50)
            AUX51 = AUX51 - TAU(51)
            AUX52 = AUX52 - TAU(52)
            AUX53 = AUX53 - TAU(53)
            AUX54 = AUX54 - TAU(54)
            AUX55 = AUX55 - TAU(55)
            AUX56 = AUX56 - TAU(56)
            AUX57 = AUX57 - TAU(57)
            AUX58 = AUX58 - TAU(58)
            AUX59 = AUX59 - TAU(59)
            AUX60 = AUX60 - TAU(60)
            AUX61 = AUX61 - TAU(61)
            AUX62 = AUX62 - TAU(62)
            AUX63 = AUX63 - TAU(63)
            AUX64 = AUX64 - TAU(64)
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
            GPL33 = G(I) + AUX33
            GPL34 = G(I) + AUX34
            GPL35 = G(I) + AUX35
            GPL36 = G(I) + AUX36
            GPL37 = G(I) + AUX37
            GPL38 = G(I) + AUX38
            GPL39 = G(I) + AUX39
            GPL40 = G(I) + AUX40
            GPL41 = G(I) + AUX41
            GPL42 = G(I) + AUX42
            GPL43 = G(I) + AUX43
            GPL44 = G(I) + AUX44
            GPL45 = G(I) + AUX45
            GPL46 = G(I) + AUX46
            GPL47 = G(I) + AUX47
            GPL48 = G(I) + AUX48
            GPL49 = G(I) + AUX49
            GPL50 = G(I) + AUX50
            GPL51 = G(I) + AUX51
            GPL52 = G(I) + AUX52
            GPL53 = G(I) + AUX53
            GPL54 = G(I) + AUX54
            GPL55 = G(I) + AUX55
            GPL56 = G(I) + AUX56
            GPL57 = G(I) + AUX57
            GPL58 = G(I) + AUX58
            GPL59 = G(I) + AUX59
            GPL60 = G(I) + AUX60
            GPL61 = G(I) + AUX61
            GPL62 = G(I) + AUX62
            GPL63 = G(I) + AUX63
            GPL64 = G(I) + AUX64
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
            IF( ABS(GPL33) .LT. PIVMIN )  GPL33 = NEGPIV
            IF( ABS(GPL34) .LT. PIVMIN )  GPL34 = NEGPIV
            IF( ABS(GPL35) .LT. PIVMIN )  GPL35 = NEGPIV
            IF( ABS(GPL36) .LT. PIVMIN )  GPL36 = NEGPIV
            IF( ABS(GPL37) .LT. PIVMIN )  GPL37 = NEGPIV
            IF( ABS(GPL38) .LT. PIVMIN )  GPL38 = NEGPIV
            IF( ABS(GPL39) .LT. PIVMIN )  GPL39 = NEGPIV
            IF( ABS(GPL40) .LT. PIVMIN )  GPL40 = NEGPIV
            IF( ABS(GPL41) .LT. PIVMIN )  GPL41 = NEGPIV
            IF( ABS(GPL42) .LT. PIVMIN )  GPL42 = NEGPIV
            IF( ABS(GPL43) .LT. PIVMIN )  GPL43 = NEGPIV
            IF( ABS(GPL44) .LT. PIVMIN )  GPL44 = NEGPIV
            IF( ABS(GPL45) .LT. PIVMIN )  GPL45 = NEGPIV
            IF( ABS(GPL46) .LT. PIVMIN )  GPL46 = NEGPIV
            IF( ABS(GPL47) .LT. PIVMIN )  GPL47 = NEGPIV
            IF( ABS(GPL48) .LT. PIVMIN )  GPL48 = NEGPIV
            IF( ABS(GPL49) .LT. PIVMIN )  GPL49 = NEGPIV
            IF( ABS(GPL50) .LT. PIVMIN )  GPL50 = NEGPIV
            IF( ABS(GPL51) .LT. PIVMIN )  GPL51 = NEGPIV
            IF( ABS(GPL52) .LT. PIVMIN )  GPL52 = NEGPIV
            IF( ABS(GPL53) .LT. PIVMIN )  GPL53 = NEGPIV
            IF( ABS(GPL54) .LT. PIVMIN )  GPL54 = NEGPIV
            IF( ABS(GPL55) .LT. PIVMIN )  GPL55 = NEGPIV
            IF( ABS(GPL56) .LT. PIVMIN )  GPL56 = NEGPIV
            IF( ABS(GPL57) .LT. PIVMIN )  GPL57 = NEGPIV
            IF( ABS(GPL58) .LT. PIVMIN )  GPL58 = NEGPIV
            IF( ABS(GPL59) .LT. PIVMIN )  GPL59 = NEGPIV
            IF( ABS(GPL60) .LT. PIVMIN )  GPL60 = NEGPIV
            IF( ABS(GPL61) .LT. PIVMIN )  GPL61 = NEGPIV
            IF( ABS(GPL62) .LT. PIVMIN )  GPL62 = NEGPIV
            IF( ABS(GPL63) .LT. PIVMIN )  GPL63 = NEGPIV
            IF( ABS(GPL64) .LT. PIVMIN )  GPL64 = NEGPIV
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
            IF( GPL33 .LT. ZERO )  NEGC33 = NEGC33 + 1
            IF( GPL34 .LT. ZERO )  NEGC34 = NEGC34 + 1
            IF( GPL35 .LT. ZERO )  NEGC35 = NEGC35 + 1
            IF( GPL36 .LT. ZERO )  NEGC36 = NEGC36 + 1
            IF( GPL37 .LT. ZERO )  NEGC37 = NEGC37 + 1
            IF( GPL38 .LT. ZERO )  NEGC38 = NEGC38 + 1
            IF( GPL39 .LT. ZERO )  NEGC39 = NEGC39 + 1
            IF( GPL40 .LT. ZERO )  NEGC40 = NEGC40 + 1
            IF( GPL41 .LT. ZERO )  NEGC41 = NEGC41 + 1
            IF( GPL42 .LT. ZERO )  NEGC42 = NEGC42 + 1
            IF( GPL43 .LT. ZERO )  NEGC43 = NEGC43 + 1
            IF( GPL44 .LT. ZERO )  NEGC44 = NEGC44 + 1
            IF( GPL45 .LT. ZERO )  NEGC45 = NEGC45 + 1
            IF( GPL46 .LT. ZERO )  NEGC46 = NEGC46 + 1
            IF( GPL47 .LT. ZERO )  NEGC47 = NEGC47 + 1
            IF( GPL48 .LT. ZERO )  NEGC48 = NEGC48 + 1
            IF( GPL49 .LT. ZERO )  NEGC49 = NEGC49 + 1
            IF( GPL50 .LT. ZERO )  NEGC50 = NEGC50 + 1
            IF( GPL51 .LT. ZERO )  NEGC51 = NEGC51 + 1
            IF( GPL52 .LT. ZERO )  NEGC52 = NEGC52 + 1
            IF( GPL53 .LT. ZERO )  NEGC53 = NEGC53 + 1
            IF( GPL54 .LT. ZERO )  NEGC54 = NEGC54 + 1
            IF( GPL55 .LT. ZERO )  NEGC55 = NEGC55 + 1
            IF( GPL56 .LT. ZERO )  NEGC56 = NEGC56 + 1
            IF( GPL57 .LT. ZERO )  NEGC57 = NEGC57 + 1
            IF( GPL58 .LT. ZERO )  NEGC58 = NEGC58 + 1
            IF( GPL59 .LT. ZERO )  NEGC59 = NEGC59 + 1
            IF( GPL60 .LT. ZERO )  NEGC60 = NEGC60 + 1
            IF( GPL61 .LT. ZERO )  NEGC61 = NEGC61 + 1
            IF( GPL62 .LT. ZERO )  NEGC62 = NEGC62 + 1
            IF( GPL63 .LT. ZERO )  NEGC63 = NEGC63 + 1
            IF( GPL64 .LT. ZERO )  NEGC64 = NEGC64 + 1
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
            AUX33 = NGN(I) * ( AUX33 / GPL33 )
            AUX34 = NGN(I) * ( AUX34 / GPL34 )
            AUX35 = NGN(I) * ( AUX35 / GPL35 )
            AUX36 = NGN(I) * ( AUX36 / GPL36 )
            AUX37 = NGN(I) * ( AUX37 / GPL37 )
            AUX38 = NGN(I) * ( AUX38 / GPL38 )
            AUX39 = NGN(I) * ( AUX39 / GPL39 )
            AUX40 = NGN(I) * ( AUX40 / GPL40 )
            AUX41 = NGN(I) * ( AUX41 / GPL41 )
            AUX42 = NGN(I) * ( AUX42 / GPL42 )
            AUX43 = NGN(I) * ( AUX43 / GPL43 )
            AUX44 = NGN(I) * ( AUX44 / GPL44 )
            AUX45 = NGN(I) * ( AUX45 / GPL45 )
            AUX46 = NGN(I) * ( AUX46 / GPL46 )
            AUX47 = NGN(I) * ( AUX47 / GPL47 )
            AUX48 = NGN(I) * ( AUX48 / GPL48 )
            AUX49 = NGN(I) * ( AUX49 / GPL49 )
            AUX50 = NGN(I) * ( AUX50 / GPL50 )
            AUX51 = NGN(I) * ( AUX51 / GPL51 )
            AUX52 = NGN(I) * ( AUX52 / GPL52 )
            AUX53 = NGN(I) * ( AUX53 / GPL53 )
            AUX54 = NGN(I) * ( AUX54 / GPL54 )
            AUX55 = NGN(I) * ( AUX55 / GPL55 )
            AUX56 = NGN(I) * ( AUX56 / GPL56 )
            AUX57 = NGN(I) * ( AUX57 / GPL57 )
            AUX58 = NGN(I) * ( AUX58 / GPL58 )
            AUX59 = NGN(I) * ( AUX59 / GPL59 )
            AUX60 = NGN(I) * ( AUX60 / GPL60 )
            AUX61 = NGN(I) * ( AUX61 / GPL61 )
            AUX62 = NGN(I) * ( AUX62 / GPL62 )
            AUX63 = NGN(I) * ( AUX63 / GPL63 )
            AUX64 = NGN(I) * ( AUX64 / GPL64 )
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
         SMT33 = AUX33 - TAU(33)
         SMT34 = AUX34 - TAU(34)
         SMT35 = AUX35 - TAU(35)
         SMT36 = AUX36 - TAU(36)
         SMT37 = AUX37 - TAU(37)
         SMT38 = AUX38 - TAU(38)
         SMT39 = AUX39 - TAU(39)
         SMT40 = AUX40 - TAU(40)
         SMT41 = AUX41 - TAU(41)
         SMT42 = AUX42 - TAU(42)
         SMT43 = AUX43 - TAU(43)
         SMT44 = AUX44 - TAU(44)
         SMT45 = AUX45 - TAU(45)
         SMT46 = AUX46 - TAU(46)
         SMT47 = AUX47 - TAU(47)
         SMT48 = AUX48 - TAU(48)
         SMT49 = AUX49 - TAU(49)
         SMT50 = AUX50 - TAU(50)
         SMT51 = AUX51 - TAU(51)
         SMT52 = AUX52 - TAU(52)
         SMT53 = AUX53 - TAU(53)
         SMT54 = AUX54 - TAU(54)
         SMT55 = AUX55 - TAU(55)
         SMT56 = AUX56 - TAU(56)
         SMT57 = AUX57 - TAU(57)
         SMT58 = AUX58 - TAU(58)
         SMT59 = AUX59 - TAU(59)
         SMT60 = AUX60 - TAU(60)
         SMT61 = AUX61 - TAU(61)
         SMT62 = AUX62 - TAU(62)
         SMT63 = AUX63 - TAU(63)
         SMT64 = AUX64 - TAU(64)
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
         GPL33 = G(I) + SMT33
         GPL34 = G(I) + SMT34
         GPL35 = G(I) + SMT35
         GPL36 = G(I) + SMT36
         GPL37 = G(I) + SMT37
         GPL38 = G(I) + SMT38
         GPL39 = G(I) + SMT39
         GPL40 = G(I) + SMT40
         GPL41 = G(I) + SMT41
         GPL42 = G(I) + SMT42
         GPL43 = G(I) + SMT43
         GPL44 = G(I) + SMT44
         GPL45 = G(I) + SMT45
         GPL46 = G(I) + SMT46
         GPL47 = G(I) + SMT47
         GPL48 = G(I) + SMT48
         GPL49 = G(I) + SMT49
         GPL50 = G(I) + SMT50
         GPL51 = G(I) + SMT51
         GPL52 = G(I) + SMT52
         GPL53 = G(I) + SMT53
         GPL54 = G(I) + SMT54
         GPL55 = G(I) + SMT55
         GPL56 = G(I) + SMT56
         GPL57 = G(I) + SMT57
         GPL58 = G(I) + SMT58
         GPL59 = G(I) + SMT59
         GPL60 = G(I) + SMT60
         GPL61 = G(I) + SMT61
         GPL62 = G(I) + SMT62
         GPL63 = G(I) + SMT63
         GPL64 = G(I) + SMT64
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
         IF( ABS(GPL33) .LT. PIVMIN )  GPL33 = NEGPIV
         IF( ABS(GPL34) .LT. PIVMIN )  GPL34 = NEGPIV
         IF( ABS(GPL35) .LT. PIVMIN )  GPL35 = NEGPIV
         IF( ABS(GPL36) .LT. PIVMIN )  GPL36 = NEGPIV
         IF( ABS(GPL37) .LT. PIVMIN )  GPL37 = NEGPIV
         IF( ABS(GPL38) .LT. PIVMIN )  GPL38 = NEGPIV
         IF( ABS(GPL39) .LT. PIVMIN )  GPL39 = NEGPIV
         IF( ABS(GPL40) .LT. PIVMIN )  GPL40 = NEGPIV
         IF( ABS(GPL41) .LT. PIVMIN )  GPL41 = NEGPIV
         IF( ABS(GPL42) .LT. PIVMIN )  GPL42 = NEGPIV
         IF( ABS(GPL43) .LT. PIVMIN )  GPL43 = NEGPIV
         IF( ABS(GPL44) .LT. PIVMIN )  GPL44 = NEGPIV
         IF( ABS(GPL45) .LT. PIVMIN )  GPL45 = NEGPIV
         IF( ABS(GPL46) .LT. PIVMIN )  GPL46 = NEGPIV
         IF( ABS(GPL47) .LT. PIVMIN )  GPL47 = NEGPIV
         IF( ABS(GPL48) .LT. PIVMIN )  GPL48 = NEGPIV
         IF( ABS(GPL49) .LT. PIVMIN )  GPL49 = NEGPIV
         IF( ABS(GPL50) .LT. PIVMIN )  GPL50 = NEGPIV
         IF( ABS(GPL51) .LT. PIVMIN )  GPL51 = NEGPIV
         IF( ABS(GPL52) .LT. PIVMIN )  GPL52 = NEGPIV
         IF( ABS(GPL53) .LT. PIVMIN )  GPL53 = NEGPIV
         IF( ABS(GPL54) .LT. PIVMIN )  GPL54 = NEGPIV
         IF( ABS(GPL55) .LT. PIVMIN )  GPL55 = NEGPIV
         IF( ABS(GPL56) .LT. PIVMIN )  GPL56 = NEGPIV
         IF( ABS(GPL57) .LT. PIVMIN )  GPL57 = NEGPIV
         IF( ABS(GPL58) .LT. PIVMIN )  GPL58 = NEGPIV
         IF( ABS(GPL59) .LT. PIVMIN )  GPL59 = NEGPIV
         IF( ABS(GPL60) .LT. PIVMIN )  GPL60 = NEGPIV
         IF( ABS(GPL61) .LT. PIVMIN )  GPL61 = NEGPIV
         IF( ABS(GPL62) .LT. PIVMIN )  GPL62 = NEGPIV
         IF( ABS(GPL63) .LT. PIVMIN )  GPL63 = NEGPIV
         IF( ABS(GPL64) .LT. PIVMIN )  GPL64 = NEGPIV
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
         IF( GPL33 .LT. ZERO )  NEGC33 = NEGC33 + 1
         IF( GPL34 .LT. ZERO )  NEGC34 = NEGC34 + 1
         IF( GPL35 .LT. ZERO )  NEGC35 = NEGC35 + 1
         IF( GPL36 .LT. ZERO )  NEGC36 = NEGC36 + 1
         IF( GPL37 .LT. ZERO )  NEGC37 = NEGC37 + 1
         IF( GPL38 .LT. ZERO )  NEGC38 = NEGC38 + 1
         IF( GPL39 .LT. ZERO )  NEGC39 = NEGC39 + 1
         IF( GPL40 .LT. ZERO )  NEGC40 = NEGC40 + 1
         IF( GPL41 .LT. ZERO )  NEGC41 = NEGC41 + 1
         IF( GPL42 .LT. ZERO )  NEGC42 = NEGC42 + 1
         IF( GPL43 .LT. ZERO )  NEGC43 = NEGC43 + 1
         IF( GPL44 .LT. ZERO )  NEGC44 = NEGC44 + 1
         IF( GPL45 .LT. ZERO )  NEGC45 = NEGC45 + 1
         IF( GPL46 .LT. ZERO )  NEGC46 = NEGC46 + 1
         IF( GPL47 .LT. ZERO )  NEGC47 = NEGC47 + 1
         IF( GPL48 .LT. ZERO )  NEGC48 = NEGC48 + 1
         IF( GPL49 .LT. ZERO )  NEGC49 = NEGC49 + 1
         IF( GPL50 .LT. ZERO )  NEGC50 = NEGC50 + 1
         IF( GPL51 .LT. ZERO )  NEGC51 = NEGC51 + 1
         IF( GPL52 .LT. ZERO )  NEGC52 = NEGC52 + 1
         IF( GPL53 .LT. ZERO )  NEGC53 = NEGC53 + 1
         IF( GPL54 .LT. ZERO )  NEGC54 = NEGC54 + 1
         IF( GPL55 .LT. ZERO )  NEGC55 = NEGC55 + 1
         IF( GPL56 .LT. ZERO )  NEGC56 = NEGC56 + 1
         IF( GPL57 .LT. ZERO )  NEGC57 = NEGC57 + 1
         IF( GPL58 .LT. ZERO )  NEGC58 = NEGC58 + 1
         IF( GPL59 .LT. ZERO )  NEGC59 = NEGC59 + 1
         IF( GPL60 .LT. ZERO )  NEGC60 = NEGC60 + 1
         IF( GPL61 .LT. ZERO )  NEGC61 = NEGC61 + 1
         IF( GPL62 .LT. ZERO )  NEGC62 = NEGC62 + 1
         IF( GPL63 .LT. ZERO )  NEGC63 = NEGC63 + 1
         IF( GPL64 .LT. ZERO )  NEGC64 = NEGC64 + 1
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
         AUX33 = - GNSQ(I) / GPL33
         AUX34 = - GNSQ(I) / GPL34
         AUX35 = - GNSQ(I) / GPL35
         AUX36 = - GNSQ(I) / GPL36
         AUX37 = - GNSQ(I) / GPL37
         AUX38 = - GNSQ(I) / GPL38
         AUX39 = - GNSQ(I) / GPL39
         AUX40 = - GNSQ(I) / GPL40
         AUX41 = - GNSQ(I) / GPL41
         AUX42 = - GNSQ(I) / GPL42
         AUX43 = - GNSQ(I) / GPL43
         AUX44 = - GNSQ(I) / GPL44
         AUX45 = - GNSQ(I) / GPL45
         AUX46 = - GNSQ(I) / GPL46
         AUX47 = - GNSQ(I) / GPL47
         AUX48 = - GNSQ(I) / GPL48
         AUX49 = - GNSQ(I) / GPL49
         AUX50 = - GNSQ(I) / GPL50
         AUX51 = - GNSQ(I) / GPL51
         AUX52 = - GNSQ(I) / GPL52
         AUX53 = - GNSQ(I) / GPL53
         AUX54 = - GNSQ(I) / GPL54
         AUX55 = - GNSQ(I) / GPL55
         AUX56 = - GNSQ(I) / GPL56
         AUX57 = - GNSQ(I) / GPL57
         AUX58 = - GNSQ(I) / GPL58
         AUX59 = - GNSQ(I) / GPL59
         AUX60 = - GNSQ(I) / GPL60
         AUX61 = - GNSQ(I) / GPL61
         AUX62 = - GNSQ(I) / GPL62
         AUX63 = - GNSQ(I) / GPL63
         AUX64 = - GNSQ(I) / GPL64
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
         OLDSMT33 = SMT33
         BRNCHI33 = ( ABS(AUX33).LE.SBOUND33 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL33) )
         OLDSMT34 = SMT34
         BRNCHI34 = ( ABS(AUX34).LE.SBOUND34 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL34) )
         OLDSMT35 = SMT35
         BRNCHI35 = ( ABS(AUX35).LE.SBOUND35 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL35) )
         OLDSMT36 = SMT36
         BRNCHI36 = ( ABS(AUX36).LE.SBOUND36 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL36) )
         OLDSMT37 = SMT37
         BRNCHI37 = ( ABS(AUX37).LE.SBOUND37 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL37) )
         OLDSMT38 = SMT38
         BRNCHI38 = ( ABS(AUX38).LE.SBOUND38 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL38) )
         OLDSMT39 = SMT39
         BRNCHI39 = ( ABS(AUX39).LE.SBOUND39 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL39) )
         OLDSMT40 = SMT40
         BRNCHI40 = ( ABS(AUX40).LE.SBOUND40 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL40) )
         OLDSMT41 = SMT41
         BRNCHI41 = ( ABS(AUX41).LE.SBOUND41 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL41) )
         OLDSMT42 = SMT42
         BRNCHI42 = ( ABS(AUX42).LE.SBOUND42 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL42) )
         OLDSMT43 = SMT43
         BRNCHI43 = ( ABS(AUX43).LE.SBOUND43 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL43) )
         OLDSMT44 = SMT44
         BRNCHI44 = ( ABS(AUX44).LE.SBOUND44 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL44) )
         OLDSMT45 = SMT45
         BRNCHI45 = ( ABS(AUX45).LE.SBOUND45 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL45) )
         OLDSMT46 = SMT46
         BRNCHI46 = ( ABS(AUX46).LE.SBOUND46 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL46) )
         OLDSMT47 = SMT47
         BRNCHI47 = ( ABS(AUX47).LE.SBOUND47 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL47) )
         OLDSMT48 = SMT48
         BRNCHI48 = ( ABS(AUX48).LE.SBOUND48 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL48) )
         OLDSMT49 = SMT49
         BRNCHI49 = ( ABS(AUX49).LE.SBOUND49 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL49) )
         OLDSMT50 = SMT50
         BRNCHI50 = ( ABS(AUX50).LE.SBOUND50 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL50) )
         OLDSMT51 = SMT51
         BRNCHI51 = ( ABS(AUX51).LE.SBOUND51 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL51) )
         OLDSMT52 = SMT52
         BRNCHI52 = ( ABS(AUX52).LE.SBOUND52 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL52) )
         OLDSMT53 = SMT53
         BRNCHI53 = ( ABS(AUX53).LE.SBOUND53 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL53) )
         OLDSMT54 = SMT54
         BRNCHI54 = ( ABS(AUX54).LE.SBOUND54 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL54) )
         OLDSMT55 = SMT55
         BRNCHI55 = ( ABS(AUX55).LE.SBOUND55 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL55) )
         OLDSMT56 = SMT56
         BRNCHI56 = ( ABS(AUX56).LE.SBOUND56 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL56) )
         OLDSMT57 = SMT57
         BRNCHI57 = ( ABS(AUX57).LE.SBOUND57 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL57) )
         OLDSMT58 = SMT58
         BRNCHI58 = ( ABS(AUX58).LE.SBOUND58 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL58) )
         OLDSMT59 = SMT59
         BRNCHI59 = ( ABS(AUX59).LE.SBOUND59 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL59) )
         OLDSMT60 = SMT60
         BRNCHI60 = ( ABS(AUX60).LE.SBOUND60 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL60) )
         OLDSMT61 = SMT61
         BRNCHI61 = ( ABS(AUX61).LE.SBOUND61 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL61) )
         OLDSMT62 = SMT62
         BRNCHI62 = ( ABS(AUX62).LE.SBOUND62 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL62) )
         OLDSMT63 = SMT63
         BRNCHI63 = ( ABS(AUX63).LE.SBOUND63 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL63) )
         OLDSMT64 = SMT64
         BRNCHI64 = ( ABS(AUX64).LE.SBOUND64 .OR.
     $                  SGNGIP.NE.SIGN(ONE,GPL64) )
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
         SMT33 = AUX33 - TAU(33)
         SMT34 = AUX34 - TAU(34)
         SMT35 = AUX35 - TAU(35)
         SMT36 = AUX36 - TAU(36)
         SMT37 = AUX37 - TAU(37)
         SMT38 = AUX38 - TAU(38)
         SMT39 = AUX39 - TAU(39)
         SMT40 = AUX40 - TAU(40)
         SMT41 = AUX41 - TAU(41)
         SMT42 = AUX42 - TAU(42)
         SMT43 = AUX43 - TAU(43)
         SMT44 = AUX44 - TAU(44)
         SMT45 = AUX45 - TAU(45)
         SMT46 = AUX46 - TAU(46)
         SMT47 = AUX47 - TAU(47)
         SMT48 = AUX48 - TAU(48)
         SMT49 = AUX49 - TAU(49)
         SMT50 = AUX50 - TAU(50)
         SMT51 = AUX51 - TAU(51)
         SMT52 = AUX52 - TAU(52)
         SMT53 = AUX53 - TAU(53)
         SMT54 = AUX54 - TAU(54)
         SMT55 = AUX55 - TAU(55)
         SMT56 = AUX56 - TAU(56)
         SMT57 = AUX57 - TAU(57)
         SMT58 = AUX58 - TAU(58)
         SMT59 = AUX59 - TAU(59)
         SMT60 = AUX60 - TAU(60)
         SMT61 = AUX61 - TAU(61)
         SMT62 = AUX62 - TAU(62)
         SMT63 = AUX63 - TAU(63)
         SMT64 = AUX64 - TAU(64)
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
         GPL33 = G(I) + SMT33
         GPL34 = G(I) + SMT34
         GPL35 = G(I) + SMT35
         GPL36 = G(I) + SMT36
         GPL37 = G(I) + SMT37
         GPL38 = G(I) + SMT38
         GPL39 = G(I) + SMT39
         GPL40 = G(I) + SMT40
         GPL41 = G(I) + SMT41
         GPL42 = G(I) + SMT42
         GPL43 = G(I) + SMT43
         GPL44 = G(I) + SMT44
         GPL45 = G(I) + SMT45
         GPL46 = G(I) + SMT46
         GPL47 = G(I) + SMT47
         GPL48 = G(I) + SMT48
         GPL49 = G(I) + SMT49
         GPL50 = G(I) + SMT50
         GPL51 = G(I) + SMT51
         GPL52 = G(I) + SMT52
         GPL53 = G(I) + SMT53
         GPL54 = G(I) + SMT54
         GPL55 = G(I) + SMT55
         GPL56 = G(I) + SMT56
         GPL57 = G(I) + SMT57
         GPL58 = G(I) + SMT58
         GPL59 = G(I) + SMT59
         GPL60 = G(I) + SMT60
         GPL61 = G(I) + SMT61
         GPL62 = G(I) + SMT62
         GPL63 = G(I) + SMT63
         GPL64 = G(I) + SMT64
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
         IF( ABS(GPL33) .LT. PIVMIN )  GPL33 = NEGPIV
         IF( ABS(GPL34) .LT. PIVMIN )  GPL34 = NEGPIV
         IF( ABS(GPL35) .LT. PIVMIN )  GPL35 = NEGPIV
         IF( ABS(GPL36) .LT. PIVMIN )  GPL36 = NEGPIV
         IF( ABS(GPL37) .LT. PIVMIN )  GPL37 = NEGPIV
         IF( ABS(GPL38) .LT. PIVMIN )  GPL38 = NEGPIV
         IF( ABS(GPL39) .LT. PIVMIN )  GPL39 = NEGPIV
         IF( ABS(GPL40) .LT. PIVMIN )  GPL40 = NEGPIV
         IF( ABS(GPL41) .LT. PIVMIN )  GPL41 = NEGPIV
         IF( ABS(GPL42) .LT. PIVMIN )  GPL42 = NEGPIV
         IF( ABS(GPL43) .LT. PIVMIN )  GPL43 = NEGPIV
         IF( ABS(GPL44) .LT. PIVMIN )  GPL44 = NEGPIV
         IF( ABS(GPL45) .LT. PIVMIN )  GPL45 = NEGPIV
         IF( ABS(GPL46) .LT. PIVMIN )  GPL46 = NEGPIV
         IF( ABS(GPL47) .LT. PIVMIN )  GPL47 = NEGPIV
         IF( ABS(GPL48) .LT. PIVMIN )  GPL48 = NEGPIV
         IF( ABS(GPL49) .LT. PIVMIN )  GPL49 = NEGPIV
         IF( ABS(GPL50) .LT. PIVMIN )  GPL50 = NEGPIV
         IF( ABS(GPL51) .LT. PIVMIN )  GPL51 = NEGPIV
         IF( ABS(GPL52) .LT. PIVMIN )  GPL52 = NEGPIV
         IF( ABS(GPL53) .LT. PIVMIN )  GPL53 = NEGPIV
         IF( ABS(GPL54) .LT. PIVMIN )  GPL54 = NEGPIV
         IF( ABS(GPL55) .LT. PIVMIN )  GPL55 = NEGPIV
         IF( ABS(GPL56) .LT. PIVMIN )  GPL56 = NEGPIV
         IF( ABS(GPL57) .LT. PIVMIN )  GPL57 = NEGPIV
         IF( ABS(GPL58) .LT. PIVMIN )  GPL58 = NEGPIV
         IF( ABS(GPL59) .LT. PIVMIN )  GPL59 = NEGPIV
         IF( ABS(GPL60) .LT. PIVMIN )  GPL60 = NEGPIV
         IF( ABS(GPL61) .LT. PIVMIN )  GPL61 = NEGPIV
         IF( ABS(GPL62) .LT. PIVMIN )  GPL62 = NEGPIV
         IF( ABS(GPL63) .LT. PIVMIN )  GPL63 = NEGPIV
         IF( ABS(GPL64) .LT. PIVMIN )  GPL64 = NEGPIV
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
         IF( GPL33 .LT. ZERO )  NEGC33 = NEGC33 + 1
         IF( GPL34 .LT. ZERO )  NEGC34 = NEGC34 + 1
         IF( GPL35 .LT. ZERO )  NEGC35 = NEGC35 + 1
         IF( GPL36 .LT. ZERO )  NEGC36 = NEGC36 + 1
         IF( GPL37 .LT. ZERO )  NEGC37 = NEGC37 + 1
         IF( GPL38 .LT. ZERO )  NEGC38 = NEGC38 + 1
         IF( GPL39 .LT. ZERO )  NEGC39 = NEGC39 + 1
         IF( GPL40 .LT. ZERO )  NEGC40 = NEGC40 + 1
         IF( GPL41 .LT. ZERO )  NEGC41 = NEGC41 + 1
         IF( GPL42 .LT. ZERO )  NEGC42 = NEGC42 + 1
         IF( GPL43 .LT. ZERO )  NEGC43 = NEGC43 + 1
         IF( GPL44 .LT. ZERO )  NEGC44 = NEGC44 + 1
         IF( GPL45 .LT. ZERO )  NEGC45 = NEGC45 + 1
         IF( GPL46 .LT. ZERO )  NEGC46 = NEGC46 + 1
         IF( GPL47 .LT. ZERO )  NEGC47 = NEGC47 + 1
         IF( GPL48 .LT. ZERO )  NEGC48 = NEGC48 + 1
         IF( GPL49 .LT. ZERO )  NEGC49 = NEGC49 + 1
         IF( GPL50 .LT. ZERO )  NEGC50 = NEGC50 + 1
         IF( GPL51 .LT. ZERO )  NEGC51 = NEGC51 + 1
         IF( GPL52 .LT. ZERO )  NEGC52 = NEGC52 + 1
         IF( GPL53 .LT. ZERO )  NEGC53 = NEGC53 + 1
         IF( GPL54 .LT. ZERO )  NEGC54 = NEGC54 + 1
         IF( GPL55 .LT. ZERO )  NEGC55 = NEGC55 + 1
         IF( GPL56 .LT. ZERO )  NEGC56 = NEGC56 + 1
         IF( GPL57 .LT. ZERO )  NEGC57 = NEGC57 + 1
         IF( GPL58 .LT. ZERO )  NEGC58 = NEGC58 + 1
         IF( GPL59 .LT. ZERO )  NEGC59 = NEGC59 + 1
         IF( GPL60 .LT. ZERO )  NEGC60 = NEGC60 + 1
         IF( GPL61 .LT. ZERO )  NEGC61 = NEGC61 + 1
         IF( GPL62 .LT. ZERO )  NEGC62 = NEGC62 + 1
         IF( GPL63 .LT. ZERO )  NEGC63 = NEGC63 + 1
         IF( GPL64 .LT. ZERO )  NEGC64 = NEGC64 + 1
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
            AUX33 = - GNSQ(I) / GPL33
            AUX34 = - GNSQ(I) / GPL34
            AUX35 = - GNSQ(I) / GPL35
            AUX36 = - GNSQ(I) / GPL36
            AUX37 = - GNSQ(I) / GPL37
            AUX38 = - GNSQ(I) / GPL38
            AUX39 = - GNSQ(I) / GPL39
            AUX40 = - GNSQ(I) / GPL40
            AUX41 = - GNSQ(I) / GPL41
            AUX42 = - GNSQ(I) / GPL42
            AUX43 = - GNSQ(I) / GPL43
            AUX44 = - GNSQ(I) / GPL44
            AUX45 = - GNSQ(I) / GPL45
            AUX46 = - GNSQ(I) / GPL46
            AUX47 = - GNSQ(I) / GPL47
            AUX48 = - GNSQ(I) / GPL48
            AUX49 = - GNSQ(I) / GPL49
            AUX50 = - GNSQ(I) / GPL50
            AUX51 = - GNSQ(I) / GPL51
            AUX52 = - GNSQ(I) / GPL52
            AUX53 = - GNSQ(I) / GPL53
            AUX54 = - GNSQ(I) / GPL54
            AUX55 = - GNSQ(I) / GPL55
            AUX56 = - GNSQ(I) / GPL56
            AUX57 = - GNSQ(I) / GPL57
            AUX58 = - GNSQ(I) / GPL58
            AUX59 = - GNSQ(I) / GPL59
            AUX60 = - GNSQ(I) / GPL60
            AUX61 = - GNSQ(I) / GPL61
            AUX62 = - GNSQ(I) / GPL62
            AUX63 = - GNSQ(I) / GPL63
            AUX64 = - GNSQ(I) / GPL64
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
            IF( BRNCHI33 )THEN
               AUX33 = G(IPREV)*SMT33 + GNSQ(IPREV)
            ELSE
               AUX33 = -AUX33*OLDSMT33 - G(IPREV)*TAU(33)
            ENDIF
            IF( BRNCHI34 )THEN
               AUX34 = G(IPREV)*SMT34 + GNSQ(IPREV)
            ELSE
               AUX34 = -AUX34*OLDSMT34 - G(IPREV)*TAU(34)
            ENDIF
            IF( BRNCHI35 )THEN
               AUX35 = G(IPREV)*SMT35 + GNSQ(IPREV)
            ELSE
               AUX35 = -AUX35*OLDSMT35 - G(IPREV)*TAU(35)
            ENDIF
            IF( BRNCHI36 )THEN
               AUX36 = G(IPREV)*SMT36 + GNSQ(IPREV)
            ELSE
               AUX36 = -AUX36*OLDSMT36 - G(IPREV)*TAU(36)
            ENDIF
            IF( BRNCHI37 )THEN
               AUX37 = G(IPREV)*SMT37 + GNSQ(IPREV)
            ELSE
               AUX37 = -AUX37*OLDSMT37 - G(IPREV)*TAU(37)
            ENDIF
            IF( BRNCHI38 )THEN
               AUX38 = G(IPREV)*SMT38 + GNSQ(IPREV)
            ELSE
               AUX38 = -AUX38*OLDSMT38 - G(IPREV)*TAU(38)
            ENDIF
            IF( BRNCHI39 )THEN
               AUX39 = G(IPREV)*SMT39 + GNSQ(IPREV)
            ELSE
               AUX39 = -AUX39*OLDSMT39 - G(IPREV)*TAU(39)
            ENDIF
            IF( BRNCHI40 )THEN
               AUX40 = G(IPREV)*SMT40 + GNSQ(IPREV)
            ELSE
               AUX40 = -AUX40*OLDSMT40 - G(IPREV)*TAU(40)
            ENDIF
            IF( BRNCHI41 )THEN
               AUX41 = G(IPREV)*SMT41 + GNSQ(IPREV)
            ELSE
               AUX41 = -AUX41*OLDSMT41 - G(IPREV)*TAU(41)
            ENDIF
            IF( BRNCHI42 )THEN
               AUX42 = G(IPREV)*SMT42 + GNSQ(IPREV)
            ELSE
               AUX42 = -AUX42*OLDSMT42 - G(IPREV)*TAU(42)
            ENDIF
            IF( BRNCHI43 )THEN
               AUX43 = G(IPREV)*SMT43 + GNSQ(IPREV)
            ELSE
               AUX43 = -AUX43*OLDSMT43 - G(IPREV)*TAU(43)
            ENDIF
            IF( BRNCHI44 )THEN
               AUX44 = G(IPREV)*SMT44 + GNSQ(IPREV)
            ELSE
               AUX44 = -AUX44*OLDSMT44 - G(IPREV)*TAU(44)
            ENDIF
            IF( BRNCHI45 )THEN
               AUX45 = G(IPREV)*SMT45 + GNSQ(IPREV)
            ELSE
               AUX45 = -AUX45*OLDSMT45 - G(IPREV)*TAU(45)
            ENDIF
            IF( BRNCHI46 )THEN
               AUX46 = G(IPREV)*SMT46 + GNSQ(IPREV)
            ELSE
               AUX46 = -AUX46*OLDSMT46 - G(IPREV)*TAU(46)
            ENDIF
            IF( BRNCHI47 )THEN
               AUX47 = G(IPREV)*SMT47 + GNSQ(IPREV)
            ELSE
               AUX47 = -AUX47*OLDSMT47 - G(IPREV)*TAU(47)
            ENDIF
            IF( BRNCHI48 )THEN
               AUX48 = G(IPREV)*SMT48 + GNSQ(IPREV)
            ELSE
               AUX48 = -AUX48*OLDSMT48 - G(IPREV)*TAU(48)
            ENDIF
            IF( BRNCHI49 )THEN
               AUX49 = G(IPREV)*SMT49 + GNSQ(IPREV)
            ELSE
               AUX49 = -AUX49*OLDSMT49 - G(IPREV)*TAU(49)
            ENDIF
            IF( BRNCHI50 )THEN
               AUX50 = G(IPREV)*SMT50 + GNSQ(IPREV)
            ELSE
               AUX50 = -AUX50*OLDSMT50 - G(IPREV)*TAU(50)
            ENDIF
            IF( BRNCHI51 )THEN
               AUX51 = G(IPREV)*SMT51 + GNSQ(IPREV)
            ELSE
               AUX51 = -AUX51*OLDSMT51 - G(IPREV)*TAU(51)
            ENDIF
            IF( BRNCHI52 )THEN
               AUX52 = G(IPREV)*SMT52 + GNSQ(IPREV)
            ELSE
               AUX52 = -AUX52*OLDSMT52 - G(IPREV)*TAU(52)
            ENDIF
            IF( BRNCHI53 )THEN
               AUX53 = G(IPREV)*SMT53 + GNSQ(IPREV)
            ELSE
               AUX53 = -AUX53*OLDSMT53 - G(IPREV)*TAU(53)
            ENDIF
            IF( BRNCHI54 )THEN
               AUX54 = G(IPREV)*SMT54 + GNSQ(IPREV)
            ELSE
               AUX54 = -AUX54*OLDSMT54 - G(IPREV)*TAU(54)
            ENDIF
            IF( BRNCHI55 )THEN
               AUX55 = G(IPREV)*SMT55 + GNSQ(IPREV)
            ELSE
               AUX55 = -AUX55*OLDSMT55 - G(IPREV)*TAU(55)
            ENDIF
            IF( BRNCHI56 )THEN
               AUX56 = G(IPREV)*SMT56 + GNSQ(IPREV)
            ELSE
               AUX56 = -AUX56*OLDSMT56 - G(IPREV)*TAU(56)
            ENDIF
            IF( BRNCHI57 )THEN
               AUX57 = G(IPREV)*SMT57 + GNSQ(IPREV)
            ELSE
               AUX57 = -AUX57*OLDSMT57 - G(IPREV)*TAU(57)
            ENDIF
            IF( BRNCHI58 )THEN
               AUX58 = G(IPREV)*SMT58 + GNSQ(IPREV)
            ELSE
               AUX58 = -AUX58*OLDSMT58 - G(IPREV)*TAU(58)
            ENDIF
            IF( BRNCHI59 )THEN
               AUX59 = G(IPREV)*SMT59 + GNSQ(IPREV)
            ELSE
               AUX59 = -AUX59*OLDSMT59 - G(IPREV)*TAU(59)
            ENDIF
            IF( BRNCHI60 )THEN
               AUX60 = G(IPREV)*SMT60 + GNSQ(IPREV)
            ELSE
               AUX60 = -AUX60*OLDSMT60 - G(IPREV)*TAU(60)
            ENDIF
            IF( BRNCHI61 )THEN
               AUX61 = G(IPREV)*SMT61 + GNSQ(IPREV)
            ELSE
               AUX61 = -AUX61*OLDSMT61 - G(IPREV)*TAU(61)
            ENDIF
            IF( BRNCHI62 )THEN
               AUX62 = G(IPREV)*SMT62 + GNSQ(IPREV)
            ELSE
               AUX62 = -AUX62*OLDSMT62 - G(IPREV)*TAU(62)
            ENDIF
            IF( BRNCHI63 )THEN
               AUX63 = G(IPREV)*SMT63 + GNSQ(IPREV)
            ELSE
               AUX63 = -AUX63*OLDSMT63 - G(IPREV)*TAU(63)
            ENDIF
            IF( BRNCHI64 )THEN
               AUX64 = G(IPREV)*SMT64 + GNSQ(IPREV)
            ELSE
               AUX64 = -AUX64*OLDSMT64 - G(IPREV)*TAU(64)
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
            AUX33 = (GNSQ(I) * AUX33) / (BDET(IPREV) * GPL33)
            AUX34 = (GNSQ(I) * AUX34) / (BDET(IPREV) * GPL34)
            AUX35 = (GNSQ(I) * AUX35) / (BDET(IPREV) * GPL35)
            AUX36 = (GNSQ(I) * AUX36) / (BDET(IPREV) * GPL36)
            AUX37 = (GNSQ(I) * AUX37) / (BDET(IPREV) * GPL37)
            AUX38 = (GNSQ(I) * AUX38) / (BDET(IPREV) * GPL38)
            AUX39 = (GNSQ(I) * AUX39) / (BDET(IPREV) * GPL39)
            AUX40 = (GNSQ(I) * AUX40) / (BDET(IPREV) * GPL40)
            AUX41 = (GNSQ(I) * AUX41) / (BDET(IPREV) * GPL41)
            AUX42 = (GNSQ(I) * AUX42) / (BDET(IPREV) * GPL42)
            AUX43 = (GNSQ(I) * AUX43) / (BDET(IPREV) * GPL43)
            AUX44 = (GNSQ(I) * AUX44) / (BDET(IPREV) * GPL44)
            AUX45 = (GNSQ(I) * AUX45) / (BDET(IPREV) * GPL45)
            AUX46 = (GNSQ(I) * AUX46) / (BDET(IPREV) * GPL46)
            AUX47 = (GNSQ(I) * AUX47) / (BDET(IPREV) * GPL47)
            AUX48 = (GNSQ(I) * AUX48) / (BDET(IPREV) * GPL48)
            AUX49 = (GNSQ(I) * AUX49) / (BDET(IPREV) * GPL49)
            AUX50 = (GNSQ(I) * AUX50) / (BDET(IPREV) * GPL50)
            AUX51 = (GNSQ(I) * AUX51) / (BDET(IPREV) * GPL51)
            AUX52 = (GNSQ(I) * AUX52) / (BDET(IPREV) * GPL52)
            AUX53 = (GNSQ(I) * AUX53) / (BDET(IPREV) * GPL53)
            AUX54 = (GNSQ(I) * AUX54) / (BDET(IPREV) * GPL54)
            AUX55 = (GNSQ(I) * AUX55) / (BDET(IPREV) * GPL55)
            AUX56 = (GNSQ(I) * AUX56) / (BDET(IPREV) * GPL56)
            AUX57 = (GNSQ(I) * AUX57) / (BDET(IPREV) * GPL57)
            AUX58 = (GNSQ(I) * AUX58) / (BDET(IPREV) * GPL58)
            AUX59 = (GNSQ(I) * AUX59) / (BDET(IPREV) * GPL59)
            AUX60 = (GNSQ(I) * AUX60) / (BDET(IPREV) * GPL60)
            AUX61 = (GNSQ(I) * AUX61) / (BDET(IPREV) * GPL61)
            AUX62 = (GNSQ(I) * AUX62) / (BDET(IPREV) * GPL62)
            AUX63 = (GNSQ(I) * AUX63) / (BDET(IPREV) * GPL63)
            AUX64 = (GNSQ(I) * AUX64) / (BDET(IPREV) * GPL64)
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
      ANEGC(33) = NEGC33
      ANEGC(34) = NEGC34
      ANEGC(35) = NEGC35
      ANEGC(36) = NEGC36
      ANEGC(37) = NEGC37
      ANEGC(38) = NEGC38
      ANEGC(39) = NEGC39
      ANEGC(40) = NEGC40
      ANEGC(41) = NEGC41
      ANEGC(42) = NEGC42
      ANEGC(43) = NEGC43
      ANEGC(44) = NEGC44
      ANEGC(45) = NEGC45
      ANEGC(46) = NEGC46
      ANEGC(47) = NEGC47
      ANEGC(48) = NEGC48
      ANEGC(49) = NEGC49
      ANEGC(50) = NEGC50
      ANEGC(51) = NEGC51
      ANEGC(52) = NEGC52
      ANEGC(53) = NEGC53
      ANEGC(54) = NEGC54
      ANEGC(55) = NEGC55
      ANEGC(56) = NEGC56
      ANEGC(57) = NEGC57
      ANEGC(58) = NEGC58
      ANEGC(59) = NEGC59
      ANEGC(60) = NEGC60
      ANEGC(61) = NEGC61
      ANEGC(62) = NEGC62
      ANEGC(63) = NEGC63
      ANEGC(64) = NEGC64
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
      AAUX(33) = AUX33
      AAUX(34) = AUX34
      AAUX(35) = AUX35
      AAUX(36) = AUX36
      AAUX(37) = AUX37
      AAUX(38) = AUX38
      AAUX(39) = AUX39
      AAUX(40) = AUX40
      AAUX(41) = AUX41
      AAUX(42) = AUX42
      AAUX(43) = AUX43
      AAUX(44) = AUX44
      AAUX(45) = AUX45
      AAUX(46) = AUX46
      AAUX(47) = AUX47
      AAUX(48) = AUX48
      AAUX(49) = AUX49
      AAUX(50) = AUX50
      AAUX(51) = AUX51
      AAUX(52) = AUX52
      AAUX(53) = AUX53
      AAUX(54) = AUX54
      AAUX(55) = AUX55
      AAUX(56) = AUX56
      AAUX(57) = AUX57
      AAUX(58) = AUX58
      AAUX(59) = AUX59
      AAUX(60) = AUX60
      AAUX(61) = AUX61
      AAUX(62) = AUX62
      AAUX(63) = AUX63
      AAUX(64) = AUX64
      END SUBROUTINE DLAXRM_STAT64
*
************************************************************************
