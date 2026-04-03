      SUBROUTINE DLAXRB_CLSSFY(
     $             N, REPR, REPI, DEPTH, SPDIAM,
     $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU,
     $             RGINFO,
     $             WIL, WIU, GIL, GIU, GAPTOL,
     $             RWORK, IWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N, DEPTH, IL, IU
      INTEGER,          INTENT(IN)     ::  WIL, WIU, GIL, GIU
      INTEGER,          INTENT(IN)     ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)     ::  SPDIAM, GAPTOL
      DOUBLE PRECISION, INTENT(IN)     ::  REPR(4*N+3)
*
      INTEGER,          INTENT(INOUT)  ::  EWL_AE( 2*IL-1 : 2*IU )
      INTEGER,          INTENT(INOUT)  ::  IWORK( 2*N + 4*(IU-IL+1) )
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU( 2*IL-1 : 2*IU )
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( IU - IL + 1 )
*
      INTEGER,          INTENT(OUT)    ::  RGINFO( IL : IU-1 )
*
*  Purpose
*  =======
*
*     Refine the eigenvalue approximation for the node [DEPTH,IL:IU]
*     in the ew-list and classify them into singletons and clusters.
*
*     Let JL:JU be the intersection of IL:IU and WIL:WIU. Then the
*     eigenvalues JL:JU are our main concern, and in principle we only
*     want to refine them. But to enable the construction of an ortho-
*     normal base , we have to stick to a conformal subtree.
*     In particular, the final classified bounds and gaps must be
*     independent of WIL:WIU.
*
*     It is helpful to think in terms of the full bisection tree be-
*     longing to each node in the representation tree: The tree of
*     intervals one gets by bisecting all initial bounds IL:IU to
*     full precision.
*
*     Let KL be maximal in IL:JL such that there is a relative gap
*     left of it or KL=IL, and analogously let KU be minimal in JU:IU
*     with a relative gap on its right or KU=IU.
*     Then we guarantee that the delivered bounds KL:KU, as well as
*     all gaps left, right and in-between, are the same as if the
*     routine were called with WIL=IL and WIU=IU. Note that the bounds
*     outside of KL:KU may be completely different.
*
*     The measure of relative gaps we use is |a-b| / absmax(a,b),
*     with one exception: relgap(0,0) = 1. This is because, if the
*     upper bound of one interval and the lower bound of the next
*     are both equal to zero, we can identify a gap regardless, since
*     only one can have zero after all (unreduced matrix).
*
*  Notes:
*  (1) Because the bounds may change after the gaps are fixed, the
*      quantity gap(gcbnds) / absmax( ub(i), lb(i+1) ) may become smaller
*      after that ( lb(i+1) can increas, for example ).
*      But the true relative gap is monotone in |ub-lb| (at least if
*      sign lb = sign ub), so that can only increase.
*
*  Preconditions:
*      IL < IU
*      GIL <= WIL <= WIU <= GIU
*      JL:JU is not empty
*
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The dimension of the matrix
*
*  REPR    (input) DOUBLE PRECISION array, dimension (4*N+3)
*  REPI    (input) INTEGER array, dimension (6+N+N/2)
*          Real and integer data used to represent the matrix
*
*  DEPTH   (input) INTEGER
*          Depth of the node in the representation tree
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          Index range belonging to the node
*
*  LGAP    (input/output) DOUBLE PRECISION
*  UGAP    (input/output) DOUBLE PRECISION
*  EWL_AE  (input/output) INTEGER array, dimension (1:2*N)
*  EWL_LU  (input/output) DOUBLE PRECISION array, dimension ( 1:2*N )
*          They constitute the initial list of eigenvalue bounds
*          (ew-list) and will hold the refinements, see DLAXRL
*          for details. Need only be initialized properly.
*
*  RGINFO  (output) INTEGER array, dimension ( IL, IU-1 )
*          RGINFO(i) holds a status flags to indicate the kind of
*          relative gap to the right of ew i, which can be
*             FULLGAP - if there is a full gap
*             NOFULL  - if there is definitely no full gap between them
*             UNKNOWN - if in the list ews i and i+1 still belong to
*                       the same interval
*          Corresponding constants are defined in gapinfo.inc.f.
*            It is guaranteed that not yet emerged gaps in the list
*          have info UNKNOWN, and that for all emerged gaps
*          the info is either NOFULL or FULLGAP.
*
*  WIL     (input) INTEGER
*  WIU     (input) INTEGER
*          The range of wanted eigenvalues.
*          Eigenvalues in the intersection of IL:IU and WIL:WIU will
*          be fully classified in any case, this intersection should
*          not be empty.
*
*  GIL     (input) INTEGER
*  GIU     (input) INTEGER
*          The global range of eigenvalues for which consistency is
*          desired, must be a superset of WIL:WIU.
*          Eigenvalues outiside of GIL:GIU will not be refined at all.
*
*  GAPTOL  (input) DOUBLE PRECISION
*          The gap tolerance, must be < 1
*
*  ======================================================================
c@extract -b parameters.inc.f procname=laxrb_clssfy
*
*     .. Parameters ..
*
*     Full gaps are also recognized if
*       absgap >= AVGAPFAC * SPDIAM / (N-1)
*     and
*       DEPTH <= MAXAVGAPDEPTH.
*     Set the latter to <0 to deactivate.
      DOUBLE PRECISION, PARAMETER  ::  AVGAPFAC = 0.3D0
      INTEGER, PARAMETER           ::  MAXAVGAPDEPTH = 10
*
*     .. Declarations ..
*
      INTERFACE
      SUBROUTINE DLAXRC(
     $            N, REPR, REPI,
     $            IL, IU, LGAP, UGAP, EWL_AE, EWL_LU,
     $            NUMEW, EWINDS,
     $            RELTOL, ABSTOL, RWORK, IWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N, IL, IU, NUMEW
      INTEGER,          INTENT(IN)     ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)     ::  RELTOL, ABSTOL
      DOUBLE PRECISION, INTENT(IN)     ::  REPR(4*N+3)
*
      INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU)
      INTEGER,          INTENT(INOUT)  ::  EWINDS(NUMEW)
      INTEGER,          INTENT(INOUT)  ::  IWORK( 2*NUMEW + IU-IL+1 )
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU)
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( NUMEW )
*
      END SUBROUTINE DLAXRC
      END INTERFACE

      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH

*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER  ::  ONE  = 1.0D0
      DOUBLE PRECISION, PARAMETER  ::  TWO  = 2.0D0

      INTEGER, PARAMETER  ::  GI_NOFULL  = -2
      INTEGER, PARAMETER  ::  GI_UNKNOWN = -1
      INTEGER, PARAMETER  ::  GI_NOGAP   =  0
      INTEGER, PARAMETER  ::  GI_INTGAP  =  1
      INTEGER, PARAMETER  ::  GI_FULLGAP =  2
*
*     .. Local Variables ..
*
      DOUBLE PRECISION  EPS, PREC, BISACC
      DOUBLE PRECISION  AVGTOL, AVGTHRESH, GAPTHRESH
      DOUBLE PRECISION  ABSTOL, RELTOL, REPELG
      DOUBLE PRECISION  LB, UB, ABSMAX, SPREAD

      INTEGER           JXWORK, IXWORK, JX2BIS, JXRBL, JXRBU
      INTEGER           ILEN, I, J, K, KL, KU, JL, JU, HL, HU
      INTEGER           KLX, KUX, SOFFL, SOFFU, IREFL, IREFU
      INTEGER           NTOBIS, JGAPL, JGAPU
      INTEGER           NRBL, NRBU, NNOCR, IOCR
      INTEGER           JXG, JXNGN, TWIST

      LOGICAL           SEEKL, SEEKU, DOAVG

      INTEGER           BISCNT
      COMMON /XMRSTATS/
     $       XTIME1, XTIME2, XTIME3,
     $       XDDDDD,
     $       XNBLCKS,
     $       XNNODES, XMAXDEPTH, XNUMFN, XNUMFT, XNUMGV, XNUMGV0,
     $       XNUMFS_2, XNUMFS_K,
     $       XNBIS_INIT,
     $       XNBIS_COB, XNBIS_IIB, XNBIS_CLASS, XNBIS_SNG, XNBIS_CLB,
     $       XNBIS_WASTED,
     $       XNRQI, XNRQIBIS, XMAXNRQI, XMAXNRQIBIS,
     $       XNENVGV, XNENVTF,
     $       XIIIII,
     $       XSTEALTHMODE
*
*           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*           !!!! SYNCHRONIZE ANY CHANGES HERE WITH xmr.h !!!!!
*           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*     Accumulative part: All these are accumulated over multiple
*     calls of stexr/laxrv.

*        Times for stages 1-3 (before/during/after block loop)
         DOUBLE PRECISION  XTIME1, XTIME2, XTIME3

*        Marker to verify data integrity
         DOUBLE PRECISION  XDDDDD

         INTEGER XNBLCKS
*
         INTEGER XNNODES, XMAXDEPTH

*        Number of calls to dlaxrn, dlaxrt, dlaxrg, and
*        dlaxrs (all twists/single twist)
         INTEGER XNUMFN, XNUMFT, XNUMGV, XNUMGV0, XNUMFS_2, XNUMFS_K

*        Number of bisection steps in various parts.
         INTEGER XNBIS_INIT
         INTEGER XNBIS_COB, XNBIS_IIB, XNBIS_CLASS, XNBIS_SNG, XNBIS_CLB

         INTEGER XNBIS_WASTED

*        Count/max rqi steps
         INTEGER XNRQI, XNRQIBIS
         INTEGER XMAXNRQI, XMAXNRQIBIS

*        For computing envelopes, number of full twisted factos
*        and gv-calls
         INTEGER XNENVGV, XNENVTF

*        Marker to verify data integrity
         INTEGER XIIIII

*    Temporaries and Configuration
*        If true, calls to internal routines dlaxr[n,t,s] are not
*        counted. Higher level routines like DLAXRB are not affected.
*        We need this to hide bisections steps done within aux
*        routines like xmr_estrc.
         LOGICAL XSTEALTHMODE

      SAVE /XMRSTATS/
*
*  ===== Executable Statements ==========================================
*
      IF( DEPTH .LE. MAXAVGAPDEPTH )THEN
         JXG = (0)
         JXNGN = (2*N+1)
         TWIST = REPI(2)

         REPELG = ZERO
         DO I = 1, N
            REPELG = MAX( REPELG, ABS(REPR(JXG+I)) )
         ENDDO
         REPELG = MAX( REPELG,
     $                 ABS(REPR(JXNGN+TWIST-1)),
     $                 ABS(REPR(JXNGN+TWIST+1)) )

*        Scale with DEPTH to account for larger residual effects
         AVGTOL = MAX(SPDIAM,REPELG) * ((DEPTH * AVGAPFAC) / (N-1))
         DOAVG  = .TRUE.
      ELSE
         AVGTOL = ZERO
         DOAVG  = .FALSE.
      ENDIF


      EPS    = DLAMCH('Epsilon')
      PREC   = DLAMCH('Precision')
      ILEN   = IU-IL+1

*     Overview of used index ranges:
*      JL:JU   are the indices of ews we want from this node
*              (intersection if IL:IU with WIL:WIU)
*      HL:HU   are the indices of ews from this node for which
*              consistency is desired, superset of JL:JU
*              (intersection of IL:IU with GIL:GIU)
*     NOTE: Neither of those need to be boundary indices of intervals.

      JL = MAX( IL, WIL )
      JU = MIN( IU, WIU )
      HL = MAX( IL, GIL )
      HU = MIN( IU, GIU )

      JXRBL  = 1
      JXRBU  = 1 + N
      JX2BIS = 1 + 2*N
      JXWORK = 1 + 2*N + ILEN
      IXWORK = 1

      IWORK(:) = 0
      RWORK(:) = ZERO


      BISCNT = XNUMFN
*
*   =====================================================================
*
*     Refine all eigenvalues enough to classify them

      ABSTOL = HALF * AVGTOL
      RELTOL = HALF * GAPTOL

      DO I = JL, JU
         IWORK(JX2BIS + I-JL) = I
      ENDDO

      CALL DLAXRC( N, REPR, REPI,
     $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU,
     $             JU-JL+1, IWORK(JX2BIS), RELTOL, ABSTOL,
     $             RWORK(IXWORK), IWORK(JXWORK)
     $           )

      XNBIS_CLASS = XNBIS_CLASS + (XNUMFN - BISCNT)
      BISCNT = XNUMFN

*   =====================================================================
*
*     Classify them

*     Note: This approach is an easy way to consistency, but means that
*     the gap tolerances are taken somewhat relaxed. With these settings,
*     we do guarantee that each classified gap has a relative width
*     exceeding GAPTOL or absolute width exceeding AVGTOL, but it may
*     be that gaps up to 2*GAPTOL or 2*ABSTOL, respectively, are not
*     recognized as such.

      AVGTHRESH = TWO * AVGTOL
      GAPTHRESH = TWO * GAPTOL
      IF( .NOT. DOAVG )  AVGTHRESH = SPDIAM

      RGINFO(:) = GI_UNKNOWN

      I = JL
      DO
         J = EWL_AE(2*I)
         IF( J .GE. JU )  EXIT
         I = J+1
*        Evaluate if there is a full gap in ]J,I[

         LB = EWL_LU(2*J-1)
         UB = EWL_LU(2*I)
*        Now LB is the lower bound of the interval containing ew J,
*        and UB is the upper bound of the interval with I.

         ABSMAX = MAX( ABS(LB), ABS(UB) )
         SPREAD = UB - LB

         IF( SPREAD .GE. MIN( AVGTHRESH, ABSMAX*GAPTHRESH ) )THEN
*           Yes, there must be a gap between J and I.
            RGINFO(J) = GI_FULLGAP
         ELSE
*           No, we don't recognize it.
*           However, due to the fact that the intervals were not refined
*           to full accuracy, there may actually be a gap here. But it
*           could not be a wide one, cf the note above.
            RGINFO(J) = GI_NOFULL
         ENDIF

         J = EWL_AE(2*I)
      ENDDO

*   =====================================================================
*
*     Handle the fringes

*     We have to find the extents of clusters that intersect JL:JU
*     but are not contained in it.


      KL   = EWL_AE(2*JL-1)
      KU   = EWL_AE(2*JU)

      NRBL = 0
      NRBU = 0

      SEEKL = .TRUE.
      SEEKU = .TRUE.
      KLX   = -1
      KUX   = -1

      SOFFL = 1
      SOFFU = 1

      DO

*        Determine actual range we are looking in
         IF( SEEKL )THEN
            DO
               IF( NRBL .EQ. 0 )THEN
                  KLX = HL-1
                  EXIT
               ENDIF
               KLX = EWL_AE( 2 * IWORK(JXRBL + NRBL-1) )
               IF( KLX .LT. KL )THEN
                  LB = EWL_LU(2*KLX-1)
                  UB = EWL_LU(2*KL)
                  ABSMAX = MAX( ABS(LB), ABS(UB) )
                  SPREAD = UB - LB
                  IF( SPREAD .GE. MIN( AVGTHRESH, ABSMAX*GAPTHRESH ) )
     $            THEN
*                    Yes, this still looks like it holds a gap
                     IF( KLX .EQ. KL-1 )THEN
                        RGINFO(KLX) = GI_FULLGAP
                        SEEKL = .FALSE.
                     ENDIF
                     EXIT
                  ENDIF
                  KL = MAX( HL, EWL_AE(2*KLX-1) )
               ENDIF
               NRBL = NRBL-1
            ENDDO
         ENDIF
         IF( SEEKU )THEN
            DO
               IF( NRBU .EQ. 0 )THEN
                  KUX = HU+1
                  EXIT
               ENDIF
               KUX = EWL_AE( 2 * IWORK(JXRBU + NRBU-1) - 1 )
               IF( KUX .GT. KU )THEN
                  LB = EWL_LU(2*KU-1)
                  UB = EWL_LU(2*KUX)
                  ABSMAX = MAX( ABS(LB), ABS(UB) )
                  SPREAD = UB - LB
                  IF( SPREAD .GE. MIN( AVGTHRESH, ABSMAX*GAPTHRESH ) )
     $            THEN
                     IF( KU .EQ. KUX-1 )THEN
                        RGINFO(KU) = GI_FULLGAP
                        SEEKU = .FALSE.
                     ENDIF
                     EXIT
                  ENDIF
                  KU = MIN( HU, EWL_AE(2*KUX) )
               ENDIF
               NRBU = NRBU-1
            ENDDO
         ENDIF

         SEEKL = SEEKL .AND. (HL .LT. KL)
         SEEKU = SEEKU .AND. (KU .LT. HU)

         IF( .NOT.SEEKL .AND. .NOT.SEEKU )  EXIT


*     We are looking for a gap between KLX AND KL, and between KU and KUX.
*     The intervals for KL and KU are refined already, and both are
*     boundary indices. Also, KLX is either HL-1 or is upper boundary
*     of an interval that was refined, and analogously for KUX.



*        Determine outside eigenvalues to refine
         NTOBIS = 0
*
         IREFL = -1
         IF( SEEKL )THEN
            IF( KLX .EQ. HL-1 )THEN
               IREFL = KL - SOFFL
               SOFFL = SOFFL * 2
            ELSE
               IREFL = (KLX + KL) / 2
            ENDIF
            IREFL = MAX(HL,IREFL)
            IWORK(JX2BIS + NTOBIS) = IREFL
            NTOBIS = NTOBIS + 1
            IWORK(JXRBL + NRBL) = IREFL
            NRBL = NRBL + 1
         ENDIF
*
         IREFU = -1
         IF( SEEKU )THEN
            IF( KUX .EQ. HU+1 )THEN
               IREFU = KU + SOFFU
               SOFFU = SOFFU * 2
            ELSE
               IREFU = (KU + KUX) / 2
            ENDIF
            IREFU = MIN(HU,IREFU)
            IWORK(JX2BIS + NTOBIS) = IREFU
            NTOBIS = NTOBIS + 1
            IWORK(JXRBU + NRBU) = IREFU
            NRBU = NRBU + 1
         ENDIF


         CALL DLAXRC( N, REPR, REPI,
     $                IL, IU, LGAP, UGAP, EWL_AE, EWL_LU,
     $                NTOBIS, IWORK(JX2BIS), RELTOL, ABSTOL,
     $                RWORK(IXWORK), IWORK(JXWORK)
     $              )


      ENDDO

*     Now KL and KU are correctly placed. Two possible scenarios, for
*     KL these are
*       KL = HL   meaning there is no gap between KL and JL, or
*       HL < KL <= JL  with RGINFO(KL-1) = GI_FULLGAP


      XNBIS_CLASS = XNBIS_CLASS + (XNUMFN - BISCNT)
      BISCNT = XNUMFN

*   =====================================================================

*     For consistency we now have to reset inner bounds for clusters
*     which concern us (intersect JL:JU) without being owned (not a
*     subset of JL:JU). The reason is that we searched for the (from our
*     perspective) outside end of the cluster starting from JL:JU, so
*     so for another call with different JL:JU, but still intersecting
*     the same cluster, the bounds inside the cluster will be different.
*     But note that the two intervals constituting the 'outer' gaps on
*     each side will always come out identical.


*     As a first step we have to determine the outermost full gaps inside
*     of JL:JU. Note that there may be none.

      K = EWL_AE(2*JL)
      DO
         IF( K .GE. JU )THEN
            K = KU
            EXIT
         ENDIF
         IF( RGINFO(K) .EQ. GI_FULLGAP )THEN
            EXIT
         ENDIF
         K = EWL_AE(2*(K+1))
      ENDDO
      JGAPL = K

      K = EWL_AE(2*JU-1)
      DO
         IF( K .LE. JGAPL )THEN
            K = KL
            EXIT
         ENDIF
         IF( RGINFO(K-1) .EQ. GI_FULLGAP )THEN
            EXIT
         ENDIF
         K = EWL_AE(2*(K-1)-1)
      ENDDO
      JGAPU = K

*     If JGAPL=KU and JGAPU=KL, then there is no gap  within JL:JU,
*     meaning the ews we want all belong to one cluster (which we may
*     or may not own).

      NNOCR = 0
      IF( JGAPL.EQ.KU .AND. JGAPU.EQ.KL )THEN
         IF( KL.LT.JL .OR. KU.GT.JU )THEN
            IWORK(1) = KL
            IWORK(2) = KU
            NNOCR = 1
         ENDIF
      ELSE
         IF( KL.LT.JL )THEN
            NNOCR = NNOCR + 1
            IWORK(1) = KL
            IWORK(2) = JGAPL
         ENDIF
         IF( KU.GT.JU )THEN
            NNOCR = NNOCR + 1
            IWORK(2*NNOCR-1) = JGAPU
            IWORK(2*NNOCR)   = KU
         ENDIF
      ENDIF

      DO IOCR = 1, NNOCR
         I  = EWL_AE( 2 * IWORK(2*IOCR-1) )
         J  = EWL_AE( 2 * IWORK(2*IOCR) - 1 )
         LB = EWL_LU( 2 * I )
         UB = EWL_LU( 2 * J - 1 )
         IF( I .LT. J )THEN
            RGINFO(I) = GI_NOFULL
            I = I+1
            J = J-1
            DO K = I, J
               EWL_AE(2*K-1) = I
               EWL_AE(2*K)   = J
               EWL_LU(2*K-1) = LB
               EWL_LU(2*K)   = UB
            ENDDO
            RGINFO(I:J-1) = GI_UNKNOWN
            RGINFO(J)     = GI_NOFULL
         ENDIF
      ENDDO
*
*   =====================================================================
*



      END SUBROUTINE DLAXRB_CLSSFY
*
*************************************************************************
