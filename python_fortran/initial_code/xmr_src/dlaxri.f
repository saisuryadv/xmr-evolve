      SUBROUTINE DLAXRI(
     $             N, D, E, WIL, WIU, NBLCKS, ABINDS,
     $             ABGERS, ABWIND, ABVSEP,
     $             RWORK, IWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N, WIL, WIU, NBLCKS
      INTEGER,          INTENT(IN)     ::  ABINDS(2*NBLCKS)
      DOUBLE PRECISION, INTENT(IN)     ::  D(N), E(N)
*
      INTEGER,          INTENT(INOUT)  ::  IWORK( 5 * NBLCKS )
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( N )
*
      INTEGER,          INTENT(OUT)    ::  ABWIND(2*NBLCKS)
      DOUBLE PRECISION, INTENT(OUT)    ::  ABGERS(2*NBLCKS)
      DOUBLE PRECISION, INTENT(OUT)    ::  ABVSEP(2*NBLCKS)
*
*
*  Purpose
*  =======
*
*    For a symmetric tridiagonal matrix T of dimension N, given by its
*    entries in D and E, and split into irreducible blocks according to
*    NBLCKS and ABINDS, determine how the eigenpair indices WIL:WIU
*    are related to the individual blocks.
*
*    PARAMOUNT design guideline: For the partial case, we want consistent
*    eigensystems and a complexity of O(kn). The latter means that we
*    cannot use full information from the Gersgorin Discs, since that
*    would require sorting them. To guarantee consistency we need the
*    separators.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The dimension of the matrix
*
*  D       (input) DOUBLE PRECISION array, dimension ( N )
*          The diagonal entries of the matrix T.
*
*  E       (input) DOUBLE PRECISION array, dimension ( N )
*          The offdiagonal elements of the matrix, should be zero at
*          positions where a block ends and positive otherwise, just
*          as DLAXRA delivers them.
*
*  WIL     (input) INTEGER
*  WIU     (input) INTEGER
*          The index range of eigenpairs of the whole matrix T that
*          are to be mapped to local block indices.
*
*  NBLCKS  (input) INTEGER
*          How many irreducible blocks the matrix T was split into,
*          NBLCKS >= 1.
*
*  ABINDS  (input) INTEGER array, dimension ( 2*NBLCKS )
*          For 1 <= i <= NBLCKS, the i'th irreducible block is the
*          principal submatrix ABLCKS(2*i-1):ABLCKS(2*i) of T.
*
*  ABGERS  (output) DOUBLE PRECISION array, dimension ( 2*NBLCKS )
*          Upon exit, ( ABGERS(2*i-1), ABGERS(2*i) ) is the union
*          of Gersgorin Discs of the i'th irreducible block.
*
*  ABWIND  (output) INTEGER array, dimension ( 2*NBLCKS )
*          Upon exit, ABWIND(2*i-1):ABWIND(2*i) are the local indices
*          of eigenpairs in the i'th block that, taken with respect to
*          the full matrix, form a part of WIL:WIU.
*          These ranges may be empty for some blocks, in which case
*          they are set to 0:-1.
*
*  ABVSEP  (output) DOUBLE PRECISION array, dimension ( 2*NBLCKS )
*          Value separators for the blocks, mainly to be used by DSTEXR.
*          Upon exit, ( ABVSEP(2*i-1), ABVSEP(2*i) ) is set to an
*          interval to which computed eigenvalues should be capped
*          to guarantee monotonicity between different calls with
*          other index ranges. Note that these need not be contained
*          in the Gersgorin Bounds, in fact they may even be disjunct
*          to them.
*
*  ======================================================================
*
*     .. Declarations ..
*
      INTERFACE
      SUBROUTINE DLAXRK(
     $             N, D, ESQ, NBLCKS, ABINDS, PIVMIN, ABGERS,
     $             INDEX,
     $             LB, NNL, NZL, ABXIL,
     $             UB, NNU, NZU, ABXIU,
     $             IWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, NBLCKS, INDEX
      INTEGER,          INTENT(IN)  ::  ABINDS(2*NBLCKS)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN
      DOUBLE PRECISION, INTENT(IN)  ::  D(N), ESQ(N)
      DOUBLE PRECISION, INTENT(IN)  ::  ABGERS(2*NBLCKS)
*
      INTEGER,          INTENT(INOUT)  ::  NNL, NZL, NNU, NZU
      INTEGER,          INTENT(INOUT)  ::  ABXIL(NBLCKS), ABXIU(NBLCKS)
      INTEGER,          INTENT(INOUT)  ::  IWORK(NBLCKS)
      DOUBLE PRECISION, INTENT(INOUT)  ::  LB, UB
*
      END SUBROUTINE DLAXRK
      END INTERFACE

      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH

*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER  ::  ONE  = 1.0D0

*
*     .. Locals ..
*
      DOUBLE PRECISION  EPS, SFMIN
      DOUBLE PRECISION  GL, GU, BGL, BGU, RTMP
      DOUBLE PRECISION  ACCTOL, PIVMIN, ESQMAX
      DOUBLE PRECISION  XA, XB, XC, XD, SEPAB, SEPCD

      INTEGER           IBLCK, BBEG, BEND, BLEN
      INTEGER           NNA, NZA, NNB, NZB, NNC, NZC, NND, NZD
      INTEGER           IXESQ, IYBXIA, IYBXIB, IYBXIC, IYBXID, IYWORK
      INTEGER           NDROP, I, J, K, CNTAB, CNTCD

      LOGICAL           DODROP

*
*  ===== Executable Statements ==========================================
*


      EPS   = DLAMCH('Epsilon')
      SFMIN = DLAMCH('Safe Miniumum')

*     For each block
*      - determine union of Gershgorin Discs
*      - init wanted indices
      GL = D(1)
      GU = D(1)

      DO IBLCK = 1, NBLCKS
         BBEG = ABINDS(2*IBLCK-1)
         BEND = ABINDS(2*IBLCK)
         BLEN = BEND - BBEG + 1

         BGL = D(BBEG) - E(BBEG)
         BGU = D(BBEG) + E(BBEG)
         DO I = BBEG+1, BEND
            RTMP  = E(I-1) + E(I)
            BGL = MIN( BGL, D(I) - RTMP )
            BGU = MAX( BGU, D(I) + RTMP )
         ENDDO

         ABGERS(2*IBLCK-1) = BGL
         ABGERS(2*IBLCK)   = BGU
         GL = MIN( GL, BGL )
         GU = MAX( GU, BGU )

         ABVSEP(2*IBLCK-1) = BGL
         ABVSEP(2*IBLCK)   = BGU

         ABWIND(2*IBLCK-1) = 1
         ABWIND(2*IBLCK)   = BLEN
      ENDDO


*     If the full spectrum is desired, we are done.
      IF( WIL.EQ.1 .AND. WIU.EQ.N )THEN
         RETURN
      ENDIF

*     If the matrix has only one block, determining the index range
*     is trivial as well. However, we still need to set the separators
*     to guarantee monotonicity between different subsets.

*     -------------------------------------------------------------------

*     Prepare bisection
      ESQMAX = ZERO
      IXESQ   = 1
      DO I = 1, N-1
         RTMP = E(I)**2
         RWORK(IXESQ-1 + I) = RTMP
         ESQMAX = MAX( ESQMAX, RTMP )
      ENDDO
      RWORK(IXESQ-1 + N) = ZERO

*     We refine up to absolute accuracy
      ACCTOL = (N * EPS) * MAX( ABS(GL), ABS(GU), GU-GL )

*     Set minimum pivot in sturm sequence (for DLAXRK)
      PIVMIN = ( SFMIN / EPS ) * MAX( ESQMAX, ONE )

      IYBXIA = 1
      IYBXIB = 1 +   NBLCKS
      IYBXIC = 1 + 2*NBLCKS
      IYBXID = 1 + 3*NBLCKS
      IYWORK = 1 + 4*NBLCKS

      XA = GL
      XB = GU
*     init lb data
      DO IBLCK = 1, NBLCKS
         BBEG = ABINDS(2*IBLCK-1)
         BEND = ABINDS(2*IBLCK)
         BLEN = BEND - BBEG + 1
         IWORK(IYBXIA-1 + IBLCK) = 0
         IWORK(IYBXIB-1 + IBLCK) = 2*BLEN
      ENDDO
      NNA = 0
      NZA = 0
      NNB = N
      NZB = 0

*     Pass 1
*     Bisect this interval as long as it contains both WIL and WIU.
*     This is true initially, since we start with IL=1, IU=N.
*
      CNTAB = 0
      DO
*        mirror lb&ub data data to WIU
         XC = XA
         XD = XB
         DO IBLCK = 1, NBLCKS
            IWORK(IYBXIC-1 + IBLCK) = IWORK(IYBXIA-1 + IBLCK)
            IWORK(IYBXID-1 + IBLCK) = IWORK(IYBXIB-1 + IBLCK)
         ENDDO
         NNC = NNA
         NZC = NZA
         NND = NNB
         NZD = NZB

         IF( (XB-XA) .LE. ACCTOL .OR.
     $       (NNA.EQ.WIL-1 .AND. NNB.EQ.WIU) )
     $   THEN
            EXIT
         ENDIF

         CALL DLAXRK(
     $          N, D, RWORK(IXESQ), NBLCKS, ABINDS, PIVMIN, ABGERS, WIL,
     $          XA, NNA, NZA, IWORK(IYBXIA),
     $          XB, NNB, NZB, IWORK(IYBXIB),
     $          IWORK(IYWORK)
     $        )
         CNTAB = CNTAB + 1

*        check if we can split the interval
         IF( (NNB+NZB) .LT. WIU )THEN
*           WIU is in the half that was just split away
            XC = XB
            DO IBLCK = 1, NBLCKS
               IWORK(IYBXIC-1 + IBLCK) = IWORK(IYBXIB-1 + IBLCK)
            ENDDO
            NNC = NNB
            NZC = NZB
            EXIT
         ELSEIF( (NNB+1) .LE. WIU )THEN
*           The bound B is exact for WIU
            XC = XB
            XD = XB
            DO IBLCK = 1, NBLCKS
               IWORK(IYBXIC-1 + IBLCK) = IWORK(IYBXIB-1 + IBLCK)
               IWORK(IYBXID-1 + IBLCK) = IWORK(IYBXIB-1 + IBLCK)
            ENDDO
            NNC = NNB
            NZC = NZB
            NND = NNB
            NZD = NZB
            EXIT
         ENDIF
      ENDDO



*
*     Now the complete data for intervals [A,B] and [C,D] are set.
*     Cases: Either both intervals are identical and small, or
*     they share one bound (the previous midpoint).

*     Left Side
      CNTAB = 0
      DO
         IF( (XB-XA).LE.ACCTOL .OR. (NNA+NZA+1).GE.WIL )THEN
            EXIT
         ENDIF

         CALL DLAXRK(
     $          N, D, RWORK(IXESQ), NBLCKS, ABINDS, PIVMIN,
     $          ABGERS, WIL,
     $          XA, NNA, NZA, IWORK(IYBXIA),
     $          XB, NNB, NZB, IWORK(IYBXIB),
     $          IWORK(IYWORK)
     $        )
         CNTAB = CNTAB + 1
      ENDDO

*     Right Side
      CNTCD = 0
      DO
         IF( (XD-XC).LE.ACCTOL .OR. NND.LE.WIU )THEN
            EXIT
         ENDIF

         CALL DLAXRK(
     $          N, D, RWORK(IXESQ), NBLCKS, ABINDS, PIVMIN,
     $          ABGERS, WIU,
     $          XC, NNC, NZC, IWORK(IYBXIC),
     $          XD, NND, NZD, IWORK(IYBXID),
     $          IWORK(IYWORK)
     $        )
         CNTCD = CNTCD + 1
      ENDDO



*     Set index ranges & separators per block
*     If now the intervals contain ews left of WIL or right of WIU,
*     then the interval was refined to full precision and we have some
*     very close (or even multiple) eigenvalues. Then we have to drop
*     indices from the per-block ranges.
      DO IBLCK = 1, NBLCKS
         ABWIND(2*IBLCK-1)  = 0
         ABWIND(2*IBLCK)    = -1
         ABVSEP(2*IBLCK-1)  = ZERO
         ABVSEP(2*IBLCK)    = ZERO
      ENDDO

*     Left Side, [A,B] contains ews NNA+1:NNB+NZB of the full matrix
      NDROP  = MAX( 0, WIL - (NNA+1) )
      DODROP = (NDROP .GT. 0)
*     The Separator
      IF( .NOT. DODROP )THEN
         SEPAB = XA
      ELSE
         SEPAB = HALF*(XA + XB)
      ENDIF
*     (drop from first block to last)
      DO IBLCK = 1, NBLCKS
         I = IWORK(IYBXIA-1 + IBLCK) / 2 + 1
         J = (IWORK(IYBXIB-1 + IBLCK) + 1) / 2
         K = I
*        local block ews I:J lie within [A,B]
         IF( I .LE. J )THEN
            K = MIN( J+1, I + NDROP )
            NDROP = NDROP - (K - I)
         ENDIF
         ABWIND(2*IBLCK-1) = K
         ABVSEP(2*IBLCK-1) = SEPAB
      ENDDO


*     Right Side, [C,D]  contains ews NNC+1:NND+NZD of the full matrix
      NDROP  = MAX( 0, (NND+NZD) - WIU )
      DODROP = (NDROP .GT. 0)
*     The Separator
      IF( .NOT. DODROP )THEN
         SEPCD = XD
      ELSE
         SEPCD = HALF*(XC + XD)
      ENDIF
*     (drop from last block to first)
      DO IBLCK = NBLCKS, 1, -1
         I = IWORK(IYBXIC-1 + IBLCK) / 2 + 1
         J = (IWORK(IYBXID-1 + IBLCK) + 1) / 2
         K = J
*        local block ews I:J lie within [C,D]
         IF( I .LE. J )THEN
            K = MAX( I-1, J - NDROP )
            NDROP = NDROP - (J - K)
         ENDIF
         ABWIND(2*IBLCK) = K
         ABVSEP(2*IBLCK) = SEPCD
      ENDDO

*     establish well-defined state for blocks without wanted ews
      DO IBLCK = 1, NBLCKS
         IF( .NOT.(ABWIND(2*IBLCK-1).LE.ABWIND(2*IBLCK)) )THEN
            ABWIND(2*IBLCK-1) = 0
            ABWIND(2*IBLCK)   = -1
            ABVSEP(2*IBLCK-1) = ZERO
            ABVSEP(2*IBLCK)   = ZERO
         ENDIF
      ENDDO



      END SUBROUTINE DLAXRI
*
************************************************************************
