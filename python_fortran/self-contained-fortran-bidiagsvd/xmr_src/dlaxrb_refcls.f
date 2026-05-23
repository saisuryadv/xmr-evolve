      SUBROUTINE DLAXRB_REFCLS(
     $             N, REPR, REPI, DEPTH,
     $             ICBEG, ICEND, LGAP, UGAP, EWL_AE, EWL_LU,
     $             RGINFO,
     $             RWORK, IWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N, DEPTH, ICBEG, ICEND
      INTEGER,          INTENT(IN)     ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)     ::  REPR(4*N+3)
*
      INTEGER,          INTENT(INOUT)  ::  EWL_AE( 2*ICBEG-1 : 2*ICEND )
      INTEGER,          INTENT(INOUT)  ::  IWORK( 4*(ICEND-ICBEG+1) )
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU( 2*ICBEG-1 : 2*ICEND )
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( (ICEND-ICBEG+1) )
*
      INTEGER,          INTENT(INOUT)  ::  RGINFO( ICBEG : ICEND-1 )
*
*  Purpose
*  =======
*
*     Refine the eigenvalue approximations for a cluster, and reveal
*     internal gaps within the cluster.
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
*  ICBEG   (input) INTEGER
*  ICEND   (input) INTEGER
*          Index range constituting the (single) cluster,
*          1 <= ICBEG < ICEND <= N.
*
*  RGINFO  (intput/output) INTEGER array, dimension ( IL : IU-1 )
*          RGINFO(i) holds a status flags to indicate the kind of
*          relative gap to the right of ew i. Should have been set
*          by DLAXRB. Corresponding constants are defined in
*          gapinfo.inc.f.
*            Since we are supposed to be given a single cluster, the
*          entries should be UNKNOWN for gaps still hidden in an
*          interval, and not GI_FULLGAP otherwise.
*            Upon exit, some of the gaps that were previously set
*          as GI_NOFULL may be changed to GI_INTGAP or GI_NOGAP.
*
*  ======================================================================
c@extract -b parameters.inc.f procname=laxrb_refcls
*
*     .. Parameters ..
*
*     Gap Tolerance for subclusters is 2**LOGSCLTOLFAC * SQRT(PREC)
      INTEGER, PARAMETER  ::  LOGSCLTOLFAC = 0

*     Cluster boundaries will be refined to relative accuracy
*     RACLBFAC*N*EPS
      INTEGER, PARAMETER  ::  RACLBFAC = 4


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
      DOUBLE PRECISION  EPS, PREC, SCLTOL, RACLB
      DOUBLE PRECISION  ABSMAX, ABSGAP, LB, UB

      INTEGER           ICLEN, NTOBIS, I, J
      INTEGER           JXWORK, IXWORK, JX2BIS

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
      EPS    = DLAMCH('Epsilon')
      PREC   = DLAMCH('Precision')
      ICLEN  = ICEND - ICBEG + 1

      SCLTOL = 2**LOGSCLTOLFAC * SQRT(PREC)
      RACLB  = RACLBFAC * N * EPS

      JX2BIS = 1
      JXWORK = ICLEN + 1
      IXWORK = 1


      BISCNT = XNUMFN


*
*   =====================================================================
*
*     Refine outer bounds

      IWORK(JX2BIS)     = ICBEG
      IWORK(JX2BIS + 1) = ICEND

      CALL DLAXRC( N, REPR, REPI,
     $             ICBEG, ICEND, LGAP, UGAP, EWL_AE, EWL_LU,
     $             2, IWORK(JX2BIS), RACLB, ZERO,
     $             RWORK(IXWORK), IWORK(JXWORK)
     $           )

      XNBIS_CLB = XNBIS_CLB + (XNUMFN - BISCNT)
      BISCNT = XNUMFN
*
*   =====================================================================
*
*     Refine all eigenvalues enough to classify them

C Extra in-cluster refinement deactivated for now. Effectively disables
C the use of internal gaps for the subset case (efficiency concerns), but
C note that for the full spectrum on level one, where the ews are refined
C to full acc using dqds, internal gaps will be classified and used
C anyway.

C      RELTOL = HALF * SCLTOL
C
C      IF( .FALSE. )THEN
C
C         DO I = ICBEG, ICEND
C            IWORK(JX2BIS + I-ICBEG) = I
C         ENDDO
C
C         CALL DLAXRC( N, REPR, REPI,
C     $                ICBEG, ICEND, LGAP, UGAP, EWL_AE, EWL_LU,
C     $                ICLEN, IWORK(JX2BIS), RELTOL, ZERO,
C     $                RWORK(IXWORK), IWORK(JXWORK)
C     $              )
C
C@LOGGING on
C         WRITE(FIDBIS,*)
C     $      '  Did ', (XNUMFN-BISCNT)*INVLEN,
C     $      '*iclen steps for classification.'
C         WRITE(FIDBIS,*) ' '
C@LOGGING !
C@STATS on
C         XNBIS_CLASS = XNBIS_CLASS + (XNUMFN - BISCNT)
C         BISCNT = XNUMFN
C@STATS !
C      ENDIF

*   =====================================================================
*
*     Classify them

      I = ICBEG
      DO
         J = EWL_AE(2*I)
         IF( J .EQ. ICEND )  EXIT
         I = J+1
*        Evaluate ]J,I[

         IF( RGINFO(J).EQ.GI_NOFULL .OR. RGINFO(J).EQ.GI_UNKNOWN )THEN

            UB = EWL_LU(2*J)
            LB = EWL_LU(2*I-1)

            ABSMAX = MAX( ABS(LB), ABS(UB) )
            ABSGAP = LB - UB

            IF( ABSGAP .GE. ABSMAX*SCLTOL )THEN
               RGINFO(J) = GI_INTGAP
            ELSE
               RGINFO(J) = GI_NOGAP
            ENDIF
         ENDIF

         J = EWL_AE(2*I)
      ENDDO
*
*   =====================================================================
*
*     Refine boundaries of internal gaps
*
      NTOBIS = 0

      J = ICBEG - 1
      DO
         I = J+1
         DO
            J = EWL_AE(2*(J+1))
            IF( J .EQ. ICEND )  EXIT
            IF( RGINFO(J) .EQ. GI_INTGAP )  EXIT
         ENDDO

         IF( ICBEG .LT. I )THEN
            IWORK( JX2BIS + NTOBIS ) = I
            NTOBIS = NTOBIS + 1
         ENDIF
         IF( I .LT. J .AND. J .LT. ICEND )THEN
            IWORK( JX2BIS + NTOBIS ) = J
            NTOBIS = NTOBIS + 1
         ENDIF

         IF( J .EQ. ICEND )  EXIT
      ENDDO

      CALL DLAXRC( N, REPR, REPI,
     $             ICBEG, ICEND, LGAP, UGAP, EWL_AE, EWL_LU,
     $             NTOBIS, IWORK(JX2BIS), RACLB, ZERO,
     $             RWORK(IXWORK), IWORK(JXWORK)
     $           )

      XNBIS_CLB = XNBIS_CLB + (XNUMFN - BISCNT)

      END SUBROUTINE DLAXRB_REFCLS
*
*************************************************************************
