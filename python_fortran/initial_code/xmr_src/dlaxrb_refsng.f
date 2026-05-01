      SUBROUTINE DLAXRB_REFSNG(
     $             N, REPR, REPI, DEPTH,
     $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, RGINFO,
     $             RWORK, IWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N, IL, IU, DEPTH
      INTEGER,          INTENT(IN)     ::  RGINFO( IL : IU-1 )
      INTEGER,          INTENT(IN)     ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)     ::  REPR(4*N+3)
*
      INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU)
      INTEGER,          INTENT(INOUT)  ::  IWORK( 2*(IU-IL+1) )
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU)
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( 2*(IU-IL+1) )
*
*
*  Purpose
*  =======
*
*     Refine bounds for singletons some more after basic classification.
*
*     Note: This could equally well be done in DLAXRX before setting up
*           the RQI, but by doing them all together here we can exploit
*           the vectorized negcounts via DLAXRM.
*
*  Arguments
*  =========
*
*  RGINFO  (input) INTEGER, dimension ( IL : IU-1 )
*          RGINFO(i) holds a status flags to indicate the kind of
*          relative gap to the right of ew i, as they have been set
*          by DLAXRB_CLASS.
*          The possible values are listed in gapinfo.inc.f.
*
*  ======================================================================
*
*     .. Declarations ..
*
      INTERFACE
      SUBROUTINE DLAXRL_REFINE(
     $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, I, J, LAMBDA, XI
     $           )
      IMPLICIT NONE
*
      INTEGER         , INTENT(IN)  ::  IL, IU, I, J, XI
      DOUBLE PRECISION, INTENT(IN)  ::  LAMBDA
*
      INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU)
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU)
      END SUBROUTINE DLAXRL_REFINE
      END INTERFACE
      INTERFACE
      SUBROUTINE  DLAXRM( N, REPR, REPI, M, ATAU, AXI )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, M
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), ATAU(M)
*
      INTEGER,          INTENT(OUT)  ::  AXI(M)
      END SUBROUTINE DLAXRM
      END INTERFACE
*
*     .. Parameters ..
*
*     Singletons intervals are refined until
*        width <= Min( GAPACCSNG*MINGAP, RELACCSNG*mid ).
*
c      DOUBLE PRECISION, PARAMETER  ::  GAPACCSNG = 0.0001D0
c      DOUBLE PRECISION, PARAMETER  ::  RELACCSNG = 0.0001D0
      DOUBLE PRECISION GAPACCSNG, RELACCSNG

*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  HALF = 0.5D0

      INTEGER, PARAMETER  ::  GI_NOFULL  = -2
      INTEGER, PARAMETER  ::  GI_UNKNOWN = -1
      INTEGER, PARAMETER  ::  GI_NOGAP   =  0
      INTEGER, PARAMETER  ::  GI_INTGAP  =  1
      INTEGER, PARAMETER  ::  GI_FULLGAP =  2
*
*     .. Locals ..
*
      DOUBLE PRECISION  MID, LB, UB, SLGAP, SUGAP, WIDTH, TOLSNG
      INTEGER           ILEN, JXEIND, JXINER, IXTAUS, IXATOL
      INTEGER           I, J, K, M, NTOBIS, GILEFT, GIRGHT

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
*     -- Executable Statements -----------------------------------------
*
C     0.0001  FAC / N
*     Note: Far better would be to go with minimal initial refinement into
*     RQI and optionally do initial refines (based on dlaxrn) of the RQC is
*     too big. But needs some restructuring there, so only worthwhile if time
*     for initial refines of sng with this tols is significant.
      GAPACCSNG = 0.01 / N
      RELACCSNG = 0.01 / N

      BISCNT = XNUMFN

      ILEN = IU - IL + 1
*
*     -----------------------
*      Auxiliary data fields
*     -----------------------
*     IWORK:
*       EIND(1:M)   -  Indices of yet unconverged intervals
*       INER(1:M)   -  inertias of their midpoints
*     RWORK:
*       TAUS(1:M)   -  midpoints of unconverged intervals
*       ATOL(1:M)   -  absolute tolerance for the singleton
*
      JXEIND = 0
      JXINER = ILEN

      IXTAUS = 0
      IXATOL = ILEN

*     -----------------
*      Scan singletons
*     -----------------
      M = 0
      I = IL
      DO
         J = EWL_AE(2*I)

         GILEFT = GI_FULLGAP
         GIRGHT = GI_FULLGAP
         SLGAP  = LGAP
         SUGAP  = UGAP
         IF( I .NE. IL )THEN
            GILEFT = RGINFO(I-1)
            SLGAP  = EWL_LU(2*I-1) - EWL_LU(2*I-2)
         ENDIF
         IF( J .NE. IU )THEN
            GIRGHT = RGINFO(J)
            SUGAP  = EWL_LU(2*J+1) - EWL_LU(2*J)
         ENDIF

         IF( I.EQ.J .AND.
     $       GILEFT.EQ.GI_FULLGAP .AND. GIRGHT.EQ.GI_FULLGAP )
     $   THEN
*           This is a singleton
            M = M + 1
            IWORK(JXEIND + M) = I
            RWORK(IXATOL + M) = GAPACCSNG * MIN(SLGAP, SUGAP)
         ENDIF

         IF( J .EQ. IU )  EXIT

         I = J+1
      ENDDO

*
*     ===========
*      Main Loop
*     ===========
*
      DO
*
*        -----------------------
*         Check for convergence
*        -----------------------
*
         K = 1
         DO
            IF( K .GT. M )  EXIT

            I      = IWORK(JXEIND + K)
            LB     = EWL_LU(2*I-1)
            UB     = EWL_LU(2*I)
            MID    = HALF*(LB + UB)
            WIDTH  = UB - LB
            TOLSNG = MIN( RWORK(IXATOL + K), RELACCSNG*ABS(MID) )

            IF( MID.EQ.LB .OR. MID.EQ.UB .OR. WIDTH.LE.TOLSNG )THEN
*              this interval is converged
               IWORK(JXEIND + K) = IWORK(JXEIND + M)
               RWORK(IXATOL + K) = RWORK(IXATOL + M)
               M = M - 1
            ELSE
*              not converged, add midpoint for one bisection step
               NTOBIS = NTOBIS + 1
               RWORK( IXTAUS + K ) = MID
               K = K + 1
            ENDIF
         ENDDO
*
*        -----------
*         Loop Exit
*        -----------
*
         IF( M .EQ. 0 )  EXIT
*
*        ---------------------------------------
*         Bisect all unconverged intervals once
*        ---------------------------------------
*
         CALL DLAXRM( N, REPR, REPI, M,
     $                RWORK(IXTAUS + 1), IWORK(JXINER + 1) )

         DO K = 1, M
            I = IWORK(JXEIND + K)

            CALL DLAXRL_REFINE(
     $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU,
     $             I, I, RWORK(IXTAUS + K), IWORK(JXINER + K)
     $      )
         ENDDO

      ENDDO

      XNBIS_SNG = XNBIS_SNG + (XNUMFN - BISCNT)

      END SUBROUTINE DLAXRB_REFSNG
*
************************************************************************

