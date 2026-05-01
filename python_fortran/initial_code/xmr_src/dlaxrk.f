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
*
*  Purpose
*  =======
*
*    Do one bisection step for a symmetric tridiagonal matrix build
*    from irreducible subblocks.
*
*    For LB, NNL is the negcount and NZL the zero-count wrt the full
*    matrix, and ABXIL holds the blockwise inertias (these we can
*    define as 2*negc + issing since the blocks are irreducible,
*    cf. DLAXRN). Analogously for UB, NNU, NZU and ABXIU.
*    This should be true upon entry and will still be true upon exit.
*
*    Basically the routine computes the sturm count of T-mid(lb,ub)
*    and based on this makes the interval smaller while still containing
*    ew INDEX.
*
*    The Gershgorin Bounds are used to avoid unnessessary sturm count
*    computations, should for some block the shift fall outside them.
*
*    Preconditions:
*    (1)  NNL+NZL+1 <= INDEX <= NNU
*    (2)  LB < UB, there must be at least one fp-number between them.
*         Since the sturm counts are backward stable only in an absolute
*         sense, it is not recommended to call this routine with UB-LB
*         smaller than about N*EPS*||T||, because beyond that, you get
*         more or less random results.
*
*    Postconditions:
*    (1)  NNL+1 <= INDEX <= NNU+NZU
*    (2)  Either
*           LB = UB, NNL+1 <= INDEX <= NNL+NZL, NNL=NNU, NZL=NZU
*         or
*           LB < UB, NNL+NZL <= NNU
*
*  ======================================================================
*
*     .. Declarations ..
*
      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH

*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  HALF = 0.5D0
*
*     .. Locals ..
*
      DOUBLE PRECISION  EPS
      DOUBLE PRECISION  MID, AUX, DP

      INTEGER           IBLCK, BBEG, BEND, BLEN
      INTEGER           TNNEG, TNZER, BNNEG, BSING, BIL, BIU, I
*
*  ===== Executable Statements ==========================================
*
      EPS = DLAMCH('Epsilon')
      MID = HALF*(LB + UB)


*
*     =====================
*      Compute Sturm Count
*     =====================
*
      TNNEG = 0
      TNZER = 0
      DO IBLCK = 1, NBLCKS
         BBEG = ABINDS(2*IBLCK-1)
         BEND = ABINDS(2*IBLCK)
         BLEN = BEND - BBEG + 1

         BIL = ABXIL(IBLCK) / 2 + 1
         BIU = (ABXIU(IBLCK) + 1) / 2
*        Interval [LB,UB] contains local ews BIL:BIU of the block.

         BNNEG = 0
         BSING = 0
         IF( MID .GE. ABGERS(2*IBLCK) )THEN
            BNNEG = BLEN
         ELSEIF( MID .LE. ABGERS(2*IBLCK-1) )THEN
            BNNEG = 0
         ELSEIF( BIL .GT. BIU )THEN
            BNNEG = ABXIL(IBLCK) / 2
         ELSE
            AUX = ZERO
            I = BBEG
            DO
               DP = (D(I) - MID) - AUX
               IF( ABS(DP).LT.PIVMIN )  DP = -PIVMIN
               IF( DP .LT. ZERO )  BNNEG = BNNEG + 1
               IF( I .EQ. BEND )  EXIT
               AUX = ESQ(I) / DP
               I = I+1
            ENDDO
            IF( DP .EQ. ZERO )  BSING = 1
         ENDIF

         IWORK(IBLCK) = 2*BNNEG + BSING

         TNNEG = TNNEG + BNNEG
         TNZER = TNZER + BSING
      ENDDO
*
*     =================
*      Bisect Interval
*     =================
*
*      Based on the sturm count and the preconditions, we can partition
*      the interval (LB,UB) into three parts:
*       (LB,MID)  contains ews ...,TNNEG
*       [ MID ]   contains ews TNNEG+1:TNNEG+TNZER
*       (MID,UB)  contains ews TNNEG+TNZER+1,...
*
      IF( TNNEG .LT. INDEX )THEN
*        Restrict lower bound to mid
         LB  = MID
         NNL = TNNEG
         NZL = TNZER
         ABXIL(:) = IWORK(1:NBLCKS)
      ENDIF
      IF( INDEX .LE. TNNEG+TNZER )THEN
*        Restrict upper bound to mid
         UB  = MID
         NNU = TNNEG
         NZU = TNZER
         ABXIU(:) = IWORK(1:NBLCKS)
      ENDIF


      END SUBROUTINE DLAXRK
*
************************************************************************
