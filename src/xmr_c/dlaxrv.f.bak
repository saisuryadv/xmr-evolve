      SUBROUTINE WSREQ_XRV(N, REQR, REQI)
      IMPLICIT NONE
      INTEGER, INTENT(IN)   ::  N
      INTEGER, INTENT(OUT)  ::  REQR, REQI
*     Convenience routine to perform a workspace query for DLAXRV.
      EXTERNAL DLAXRV
      DOUBLE PRECISION  R, A(2)
      INTEGER           I, T(2)
      LOGICAL           L
      REQR = -1
      REQI = -1
      CALL DLAXRV( N, A, A, T, T, A,
     $             I, I, R, R,
     $             A, A, I, T,
     $             A, REQR, T, REQI, I )
      END SUBROUTINE WSREQ_XRV
*
************************************************************************
*
      SUBROUTINE DLAXRV(
     $             N, E, ROOTR, ROOTI, EWL_AE, EWL_LU,
     $             WIL, WIU, SPDIAM, GAPTOL,
     $             W, Z, LDZ, ISUPPZ,
     $             RWORK, LRWORK, IWORK, LIWORK, INFO
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, WIL, WIU, LDZ
      INTEGER,          INTENT(IN)  ::  ROOTI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  SPDIAM, GAPTOL
      DOUBLE PRECISION, INTENT(IN)  ::  E( N-1 ), ROOTR(4*N+3)
*
      INTEGER,          INTENT(INOUT)  ::  LRWORK, LIWORK
      INTEGER,          INTENT(INOUT)  ::  EWL_AE( 1:2*N )
      INTEGER,          INTENT(INOUT)  ::  IWORK( LIWORK )
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU( 1:2*N )
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( LRWORK )
*
      INTEGER,          INTENT(OUT)  ::  INFO
      INTEGER,          INTENT(OUT)  ::  ISUPPZ( 2*WIL-1 : 2*WIU )
      DOUBLE PRECISION, INTENT(OUT)  ::  Z( LDZ, WIL : WIU )
      DOUBLE PRECISION, INTENT(OUT)  ::  W( WIL : WIU )
*
*  Purpose
*  =======
*
*    Compute all eigenpairs WIL:WIU of the symmetric tridiagonal matrix
*    given by (N,E,ROOTR,ROOTI) using the MRRR algorithm.
*
*    This is an internal routine, designed to be called from DSTEXR.
*    In particular we expect that the matrix was properly scaled and
*    split.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The dimension of the matrix
*
*  E       (input) DOUBLE PRECISION array, dimension ( 1:N-1 )
*          The offdiagonal elements of the matrix, should be positive.
*
*  ROOTR   (input) DOUBLE PRECISION array, dimension (4*N+3)
*  ROOTI   (input) INTEGER array, dimension (6+N+N/2)
*          Real and integer data used to represent the matrix, cf
*          representation.txt
*
*  EWL_AE  (input/output) INTEGER array, dimension (1:2*N)
*  EWL_LU  (input/output) DOUBLE PRECISION array, dimension ( 1:2*N )
*          Together they constitute the initial list of eigenvalue
*          bounds (ew-list), see ewlist.txt and DLAXRL for details.
*             The list must be initialized, but other than that we do
*          not specify how good the approximations are - in principle
*          there may just be one interval in there (e.g. 0:SPDIAM for
*          a p.d. matrix).
*          Will be destroyed on exit, that is, the bounds will then not
*          be valid for the root representation anymore.
*
*  WIL     (input) INTEGER
*  WIU     (input) INTEGER
*          The range of wanted eigenpairs, 1 <= WIL <= WIU <= N.
*
*  SPDIAM  (input) DOUBLE PRECISION
*          The spectral diameter of the matrix.
*
*  GAPTOL  (input) DOUBLE PRECISION
*          Gap Tolerance used to classify eigenvalues into singletons
*          and clusters.
*
*  W       (output) DOUBLE PRECISION array, dimension ( WIL:WIU )
*          The computed eigenvalues WIL:WIU.
*            If the i'th eigenvector could not be computed (see ISUPPZ),
*          then W(i) is set regardless to something that makes sense at
*          least with regard to the other eigenvalues
*
*  Z       (output) DOUBLE PRECISION array, dimension ( LDZ, WIL:WIU )
*          The computed orthonormal eigenvectors WIL:WIU.
*            If the i'th eigenvector could not be computed (see ISUPPZ),
*          then the column Z(:,i) is set zo zero.
*
*  LDZ     (input) INTEGER
*          Leading dimension of the array Z.
*
*  ISUPPZ  (output) INTEGER array, dimension ( 2*WIL-1 : 2*WIU )
*          If the i'th vector was succsessfully computed (cf INFO), where
*          WIL <= i <= WIU, then the vector Z(i) has nonzero entries
*          only within ISUPPZ(2*i-1):ISUPPZ(2*i), which is a subset
*          of 1:N.
*          Failure to compute the vector is indicated by
*          ISUPPZ(2*i-1) = 0, ISUPPZ(2*i) < 0, the precise value of
*          the latter may give information about the kind of failure:
*          -1:  depth limit exceeded
*          -2:  no rep found in DLAXRF
*
*  RWORK   (workspace) DOUBLE PRECISION, dimension (LRWORK)
*
*  LRWORK  (input/output) INTEGER
*          Dimension of the array RWORK. Set to -1 for a workspace
*          query, then the routine only computes the required real
*          workspace, sets LRWORK to this value and returns.
*            For a workspace query only N needs to be set, no other
*          argument is referenced. You can do a real and integer ws
*          query at once.
*
*  IWORK   (workspace) INTEGER, dimension (LIWORK)
*
*  LIWORK  (input/output) INTEGER
*          Dimension of the array IWORK. Set to -1 for a workspace
*          query, then the routine only computes the required integer
*          workspace, sets LIWORK to this value and returns.
*            For a workspace query only N needs to be set, no other
*          argument is referenced. You can do a real and integer ws
*          query at once.
*
*  INFO    (output) INTEGER
*          -1 - Not enough workspace given
*           0 - Everything ok, all desired eigenpairs were successfully
*               computed
*           1 - Some eigenpairs were not computed, see ISUPPZ for details
*
*  ======================================================================

Cextract -b parameters.inc.f procname=laxrv
*
*     .. Parameters ..
*
*     By how much ULP the shift should be perturbed componentwise always.
*     Set to 0 to deactivate.
      INTEGER, PARAMETER  ::  NULPSHIFTPERT = 0

*     By how much ULP the primary data of the father rep should be
*     perturbed for a gap retry in DLAXRF.
*     Set to 0 to deactivate gap retries.
      INTEGER, PARAMETER  ::  NULPGRTPERT = 16

*     Maximal allowed depth of a node in the tree (root has depth 0).
*     Note: We need to set this to a small constant since we do a depth
*     first traversal, ie, MAXDEPTH influences how much workspace we need.
      INTEGER, PARAMETER  ::  MAXDEPTH = 10
*
*     .. Declarations ..
*
      INTERFACE
      SUBROUTINE WSREQ_XRF(N, ICBEG, ICEND, REQR, REQI)
      IMPLICIT NONE
      INTEGER, INTENT(IN)   ::  N, ICBEG, ICEND
      INTEGER, INTENT(OUT)  ::  REQR, REQI
      END SUBROUTINE WSREQ_XRF
*
************************************************************************
*
      SUBROUTINE DLAXRF(
     $            N, E, GAPTOL, SPDIAM, DEPTH, TAUBAR, ICBEG, ICEND,
     $            FREPR, FREPI, FLGAP, FUGAP, FEWL_AE, FEWL_LU, FRGINFO,
     $            DOGRT, FREPR_GRT, SHFMOD, ENV, MODE,
     $            SREPR, SREPI, SLGAP, SUGAP, SEWL_AE, SEWL_LU, TAU,
     $            RWORK, LRWORK, IWORK, LIWORK, INFO
     $           )
      IMPLICIT NONE
* note: we want to be able to call this routine both within a df or a
* bf tree traversion strategy. For bf, the EWI-List of the son would
* be copied over the father afterwards and cached quants discarded,
* for df however we can directly use them.
*
      INTEGER,          INTENT(IN)  ::  N, DEPTH, ICBEG, ICEND, MODE
      INTEGER,          INTENT(IN)  ::  FREPI(6+N+N/2)
      INTEGER,          INTENT(IN)  ::  FRGINFO(ICBEG-1:ICEND)
      DOUBLE PRECISION, INTENT(IN)  ::  GAPTOL, SPDIAM, TAUBAR
      DOUBLE PRECISION, INTENT(IN)  ::  E(N-1), FREPR(4*N+3)
      DOUBLE PRECISION, INTENT(IN)  ::  FREPR_GRT(4*N+3)
      DOUBLE PRECISION, INTENT(IN)  ::  SHFMOD(3*N), ENV(N)
      LOGICAL,          INTENT(IN)  ::  DOGRT
*
      INTEGER,          INTENT(INOUT)  ::  LRWORK, LIWORK
      INTEGER,          INTENT(INOUT)  ::  IWORK( LIWORK )
      INTEGER,          INTENT(INOUT)  ::  FEWL_AE(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(INOUT)  ::  FLGAP, FUGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  FEWL_LU(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( LRWORK )
*
      INTEGER,          INTENT(OUT)  ::  INFO
      INTEGER,          INTENT(OUT)  ::  SREPI(6+N+N/2)
      INTEGER,          INTENT(OUT)  ::  SEWL_AE(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(OUT)  ::  SLGAP, SUGAP, TAU
      DOUBLE PRECISION, INTENT(OUT)  ::  SREPR(4*N+3)
      DOUBLE PRECISION, INTENT(OUT)  ::  SEWL_LU(2*ICBEG-1:2*ICEND)
      END SUBROUTINE DLAXRF
      END INTERFACE
      INTERFACE
      SUBROUTINE WSREQ_XRX(N,REQR,REQI)
      IMPLICIT NONE
      INTEGER, INTENT(IN)     ::  N
      INTEGER, INTENT(INOUT)  ::  REQR, REQI
      END SUBROUTINE WSREQ_XRX
*
************************************************************************
*
      SUBROUTINE DLAXRX(
     $             N, E, REPR, REPI,
     $             INDEX, MINGAP, LAMBDA,
     $             Z, ISUPPZ, RWORK, LRWORK, IWORK, LIWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N, INDEX
      INTEGER,          INTENT(IN)     ::  LRWORK, LIWORK
      INTEGER,          INTENT(IN)     ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)     ::  MINGAP
      DOUBLE PRECISION, INTENT(IN)     ::  E(N-1), REPR(4*N+3)
*
      INTEGER,          INTENT(INOUT)  ::  IWORK( LIWORK )
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( LRWORK ), LAMBDA
*
      INTEGER,          INTENT(OUT)    ::  ISUPPZ(2)
      DOUBLE PRECISION, INTENT(OUT)    ::  Z(N)
      END SUBROUTINE DLAXRX
      END INTERFACE
      INTERFACE
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
      END SUBROUTINE DLAXRB_CLSSFY
      END INTERFACE
      INTERFACE
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
      END SUBROUTINE DLAXRB_REFSNG
      END INTERFACE
      INTERFACE
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
      END SUBROUTINE DLAXRB_REFCLS
      END INTERFACE

      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH

      EXTERNAL          DLARNV


*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER  ::  ONE  = 1.0D0




*
*     .. Local Variables ..
*
      DOUBLE PRECISION  ATAUBAR( 0 : MAXDEPTH )
      DOUBLE PRECISION  ANLGAP( 0 : MAXDEPTH ), ANUGAP( 0 : MAXDEPTH )
      DOUBLE PRECISION  TAUBAR, LGAP, UGAP, CLGAP, CUGAP, MINENV
      DOUBLE PRECISION  LAMBDA, TAU, PREC, FAC, PIV

      INTEGER           AIL( 0 : MAXDEPTH ), AIU( 0 : MAXDEPTH )
      INTEGER           ANC( 0 : MAXDEPTH )
      INTEGER           DEPTH, IL, IU, ILEN, I, J, K, WI, WJ
      INTEGER           IXREP, JXREP, IXROOT, JXROOT
      INTEGER           IXSHFP, IXENV, IXFGRT, IXPRTE, IXPREP
      INTEGER           JXRGI, IXGBND, IXCEWB, JXCEWB, IXWORK, JXWORK
      INTEGER           XRFSTA, STATUS, RSTRIDER, RSTRIDEI, NPAIRS
      INTEGER           WSI_XRB, WSR_XRB, WSI_XRF, WSR_XRF
      INTEGER           WSI_XRX, WSR_XRX
      INTEGER           WSREMR, WSREMI, WSREQI, WSREQR
      INTEGER           ISEED(4), XRFMODE
      INTEGER           TWIST, TYPE

      LOGICAL           DOGRT


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
      INTEGER IXF77A, IXF77B, IXF77C

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

*     Number of real/integer data items for a representation
      RSTRIDER = (4*N+3)
      RSTRIDEI = (6+N+N/2)

*     ----------------------
*      Workspace Management
*     ----------------------

      IXWORK = 1
      JXWORK = 1

*     Right gap infos, indexed 0:N
      JXRGI = JXWORK
      JXWORK = JXWORK + N+1

*     In CEWB we buffer the ewl for the child as given by DLAXRF,
*     before copying it again over EWL.
      IXCEWB = IXWORK
      JXCEWB = JXWORK
      IXWORK = IXWORK + 2*N
      JXWORK = JXWORK + 2*N

*     Perturbation factors for the shift, we use three different
*     sets, one for shifts to the left, one identical to 1 for np
*     perturbation, and one for shifts to the right, cf DLAXRF.
      IXSHFP = IXWORK
      IXWORK = IXWORK + 3*N

*     Gap boundaries from the times when the gaps were classified.
*     gcbnd(2*j) and gcbnd(2*j+1) constitute gap j:j+1
*     We define the index IXGBND such that
*        RWORK(IXGBND + k) = gcbnd(k),  2 <= k <= 2N-1,
*     that is, we start with gcbnd(2), which is to hold the bound right
*     of ew 1.
      IXGBND = IXWORK - 2
      IXWORK = IXWORK + 2*N-2

*     The representation stack.
      JXROOT = JXWORK
      IXROOT = IXWORK
      JXWORK = JXWORK + (MAXDEPTH+1)*RSTRIDEI
      IXWORK = IXWORK + (MAXDEPTH+1)*RSTRIDER

*     The envelope for the current representation
      IXENV = IXWORK
      IXWORK = IXWORK + N

*     Data for perturbed representations to use for gap retries.
*       FGRT: Perturbation factors for the primary data of the current
*             representation.
*       PRTE: Perturbed versions of the offdiagonal entries.
*       PREP: Buffer to hold the perturbed representation data
*             (real only).
      DOGRT  = ( NULPGRTPERT .GT. 0 )
      IF( DOGRT )THEN
         IXFGRT = IXWORK
         IXPRTE = IXWORK + N
         IXPREP = IXWORK + 2*N-1
         IXWORK = IXWORK + 2*N-1 + RSTRIDER
      ELSE
*        buffers won't be accessed
         IXFGRT = 1
         IXPRTE = 1
         IXPREP = 1
      ENDIF
*
*     --------------------------
*      Handle workspace queries
*     --------------------------
*
      WSREMI = LIWORK - JXWORK + 1
      WSREMR = LRWORK - IXWORK + 1

*     Combine ws for dlaxrb_[clssfy|refsng|refcls]
      WSR_XRB = N
      WSI_XRB = 7*N
      CALL WSREQ_XRX( N, WSR_XRX, WSI_XRX )
      CALL WSREQ_XRF( N, 1, N, WSR_XRF, WSI_XRF )

      WSREQI = JXWORK-1 + MAX( WSI_XRB, WSI_XRX, WSI_XRF )
      WSREQR = IXWORK-1 + MAX( WSR_XRB, WSR_XRX, WSR_XRF )
      IF( LIWORK.EQ.-1 .OR. LRWORK.EQ.-1 )THEN
         IF( LIWORK.EQ.-1 )  LIWORK = WSREQI
         IF( LRWORK.EQ.-1 )  LRWORK = WSREQR
         RETURN
      ENDIF
      IF( WSREQI.GT.LIWORK .OR. WSREQR.GT.LRWORK )THEN
         INFO = -1
         RETURN
      ENDIF


*
*     -------------------------------------------------------------------
*

      PREC = DLAMCH('Precision')


*
*     .. Prepare shift random perturbation ..
*
      ISEED(1) = 29
      ISEED(2) = 6
      ISEED(3) = 1978
      ISEED(4) = 31

      IF( NULPSHIFTPERT .EQ. 0 )THEN
         DO IXF77A = IXSHFP, IXSHFP+3*N-1
            RWORK(IXF77A) = ONE
         ENDDO
      ELSE
         FAC = NULPSHIFTPERT * PREC
         CALL DLARNV( 1, ISEED, N, RWORK(IXSHFP) )
         DO I = 0, N-1
            RWORK(IXSHFP + I) = ONE  -  RWORK(IXSHFP + I) * FAC
         ENDDO
         DO IXF77A = IXSHFP+N, IXSHFP+2*N-1
            RWORK(IXF77A) = ONE
         ENDDO
         CALL DLARNV( 1, ISEED, N, RWORK(IXSHFP+2*N) )
         DO I = 2*N, 3*N-1
            RWORK(IXSHFP + I) = ONE  +  RWORK(IXSHFP + I) * FAC
         ENDDO
      ENDIF

*
*     .. Prepare perturbation factors for the representation data ..
*
      IF( DOGRT )THEN
         FAC = NULPGRTPERT * PREC
         CALL DLARNV( 2, ISEED, N, RWORK(IXFGRT) )
         DO I = 0, N-1
            RWORK(IXFGRT + I) = ONE  +  RWORK(IXFGRT + I) * FAC
         ENDDO
         CALL DLARNV( 2, ISEED, N-1, RWORK(IXPRTE) )
         DO I = 0, N-2
            RWORK(IXPRTE + I) =
     $         E(I+1) * (ONE  +  RWORK(IXPRTE + I) * FAC)
         ENDDO
      ENDIF


      INFO = 0
      NPAIRS = 0
      DO I = 1, 2*(WIU-WIL+1)
         ISUPPZ(I) = 0
      ENDDO

*     ****************************************************************
*     *                        MAIN LOOP                             *
*     *                       -----------                            *
*     *  Perform a depth-first traversion of the representation tree *
*     ****************************************************************
*
      IWORK(JXRGI)     = 2
      IWORK(JXRGI + N) = 2

      DO IXF77A = 1, RSTRIDER
         RWORK(IXROOT + IXF77A - 1) = ROOTR(IXF77A)
      ENDDO
      DO IXF77A = 1, RSTRIDEI
         IWORK(JXROOT + IXF77A - 1) = ROOTI(IXF77A)
      ENDDO

      AIL(0) = 1
      AIU(0) = N
      ANC(0) = 1
      ATAUBAR(0) = ZERO
*     We do know not much about how accurate the initial bounds are.
*     The near-zero bound only needs to be set to sth > 0, as dlaxrb
*     will update it automatically to at lest |ew|.
      ANLGAP(0) = MAX( SPDIAM, ABS(EWL_LU(1)) )
      ANUGAP(0) = MAX( SPDIAM, ABS(EWL_LU(2*N)) )

      DEPTH = 0
      DO
         IL     = AIL(DEPTH)
         IU     = AIU(DEPTH)
         ILEN   = IU-IL+1
         TAUBAR = ATAUBAR(DEPTH)

         IXREP = IXROOT + DEPTH*RSTRIDER
         JXREP = JXROOT + DEPTH*RSTRIDEI


         IF( ANC(DEPTH) .EQ. IL )THEN
*
*           Classify eigenvalues
*
C            CALL DLAXRB(
C     $             N, RWORK(IXREP), IWORK(JXREP), DEPTH, SPDIAM,
C     $             IL, IU, ANLGAP(DEPTH), ANUGAP(DEPTH),
C     $             EWL_AE(2*IL-1), EWL_LU(2*IL-1),
C     $             IWORK(JXRGI + IL), RWORK(IXGBND + 2*IL),
C     $             WIL, WIU, GAPTOL, .TRUE.,
C     $             RWORK(IXWORK), WSR_XRB, IWORK(JXWORK), WSI_XRB
C     $           )

            CALL DLAXRB_CLSSFY(
     $             N, RWORK(IXREP), IWORK(JXREP), DEPTH, SPDIAM,
     $             IL, IU, ANLGAP(DEPTH), ANUGAP(DEPTH),
     $             EWL_AE(2*IL-1), EWL_LU(2*IL-1), IWORK(JXRGI + IL),
     $             WIL, WIU, 1, N, GAPTOL,
     $             RWORK(IXWORK), IWORK(JXWORK)
     $           )

*           Record bounds at time of classification
            DO IXF77A = IL, IU-1
               RWORK(IXGBND + 2*IXF77A) = EWL_LU(2*IXF77A)
            ENDDO

*           Initial refinement of singletons
            CALL DLAXRB_REFSNG(
     $             N, RWORK(IXREP), IWORK(JXREP), DEPTH,
     $             IL, IU, ANLGAP(DEPTH), ANUGAP(DEPTH),
     $             EWL_AE(2*IL-1), EWL_LU(2*IL-1), IWORK(JXRGI + IL),
     $             RWORK(IXWORK), IWORK(JXWORK)
     $           )

            XNNODES = XNNODES + 1
            XMAXDEPTH = MAX( XMAXDEPTH, DEPTH )
         ENDIF
*
*        Handle children, singletons right here, clusters onto the
*        stack.
*
         DO
            IF( ANC(DEPTH) .GT. IU )THEN
*              backtrack
               DEPTH = DEPTH-1
               EXIT
            ENDIF

            I = ANC(DEPTH)
            J = I
            DO
               IF( IWORK(JXRGI + J) .EQ. 2 )  EXIT
               J = J+1
            ENDDO
            ANC(DEPTH) = J+1

*
*           Prune children that do not contain wanted eigenpairs
*
            IF( J.LT.WIL .OR. I.GT.WIU )  CYCLE
            WI = MAX(WIL, I)
            WJ = MIN(WIU, J)

*           Determine left and right gaps for this child. Here, using
*           for the far bounds the ones from when the gap was first
*           identified leads to conformal subtrees.
            IF( I .EQ. IL )THEN
               LGAP = ANLGAP(DEPTH)
            ELSE
               LGAP = EWL_LU( 2*I-1 ) - RWORK( IXGBND + 2*(I-1) )
            ENDIF
            IF( J .EQ. IU )THEN
               UGAP = ANUGAP(DEPTH)
            ELSE
               UGAP = RWORK( IXGBND + 2*J+1 ) - EWL_LU( 2*J )
            ENDIF


            IF( I .EQ. J )THEN
*              -----------
*               Singleton
*              -----------
*
               LAMBDA = HALF * (EWL_LU(2*I-1) + EWL_LU(2*I))

               CALL DLAXRX( N, E,
     $                RWORK(IXREP), IWORK(JXREP),
     $                I, MIN(LGAP,UGAP), LAMBDA,
     $                Z(1,I), ISUPPZ(2*I-1),
     $                RWORK(IXWORK), WSR_XRX, IWORK(JXWORK), WSI_XRX
     $              )

               W(I) = LAMBDA + TAUBAR

               NPAIRS = NPAIRS + 1
            ELSE
*              ---------
*               Cluster
*              ---------
*
               STATUS = 0
               TAU = ZERO

               IF( DEPTH .EQ. MAXDEPTH )THEN
                  STATUS = -1
               ELSE
*                 1 means try outside and also exploit subgaps
                  XRFMODE = 1
C                  XRFMODE = 2

*                 Refine the cluster
                  CALL DLAXRB_REFCLS(
     $                   N, RWORK(IXREP), IWORK(JXREP), DEPTH, I, J,
     $                   LGAP, UGAP, EWL_AE(2*I-1), EWL_LU(2*I-1),
     $                   IWORK(JXRGI+I),
     $                   RWORK(IXWORK), IWORK(JXWORK)
     $                 )

*                 Compute Envelope
                  CALL DLAXRF_ENV(
     $                   N, E, RWORK(IXREP), IWORK(JXREP), DEPTH, I, J,
     $                   LGAP, UGAP, EWL_AE(2*I-1), EWL_LU(2*I-1),
     $                   IWORK(JXRGI+I-1),
     $                   XRFMODE, RWORK(IXENV), MINENV,
     $                   RWORK(IXWORK), IWORK(JXWORK)
     $                 )

*                 Setup perturbed rep for gap retries
*                 Note: Due to the depth first approach, we have to do
*                 this here, possibly multiple  times for the same rep.
                  IF( DOGRT )THEN
                     DO K = 0, N-1
                        RWORK(IXPREP + (0) + K)
     $                     =   RWORK(IXREP + (0) + K)
     $                       * RWORK(IXFGRT + K)
                     ENDDO
                     PIV   = RWORK(IXREP-1 + (4*N+3))
                     TWIST = IWORK(JXREP-1 + (2))
                     TYPE  = IWORK(JXREP-1 + (1))
                     CALL DLAXRR( N, TWIST, TYPE, RWORK(IXPRTE), PIV,
     $                            RWORK(IXPREP), IWORK(JXREP)
     $                    )
                  ENDIF
                  CALL DLAXRF(
     $                   N, E, GAPTOL, SPDIAM, DEPTH, TAUBAR, I, J,
     $                   RWORK(IXREP), IWORK(JXREP), LGAP, UGAP,
     $                   EWL_AE(2*I-1), EWL_LU(2*I-1), IWORK(JXRGI+I-1),
     $                   DOGRT, RWORK(IXPREP), RWORK(IXSHFP),
     $                   RWORK(IXENV), XRFMODE,
     $                   RWORK(IXREP+RSTRIDER), IWORK(JXREP+RSTRIDEI),
     $                   CLGAP, CUGAP, IWORK(JXCEWB), RWORK(IXCEWB),
     $                   TAU,
     $                   RWORK(IXWORK), WSR_XRF, IWORK(JXWORK), WSI_XRF,
     $                   XRFSTA
     $                 )

                  IF( XRFSTA .NE. 0 )THEN
*                    Could not find a child representation
                     STATUS = -2
                  ELSE
                     DO K = 0, 2*(J-I+1)-1
                        EWL_LU(2*I-1 + K) = RWORK(IXCEWB + K)
                        EWL_AE(2*I-1 + K) = IWORK(JXCEWB + K)
                     ENDDO
                  ENDIF
               ENDIF


               IF( STATUS .NE. 0 )THEN
                  DO IXF77A = 2*WI, 2*WJ, 2
                     ISUPPZ(IXF77A) = STATUS
                  ENDDO
                  DO IXF77A = 1, N
                  DO IXF77B = WI, WJ
                     Z(IXF77A, IXF77B) = ZERO
                  ENDDO
                  ENDDO
                  DO IXF77A = WI, WJ
                     W(IXF77A) = (EWL_LU(2*I-1) + EWL_LU(2*J))/2 + TAUBAR
                  ENDDO
               ELSE
                  DEPTH = DEPTH+1
                  ATAUBAR(DEPTH) = TAUBAR + TAU
                  ANLGAP(DEPTH)  = CLGAP
                  ANUGAP(DEPTH)  = CUGAP
                  AIL(DEPTH) = I
                  AIU(DEPTH) = J
                  ANC(DEPTH) = I
                  EXIT
               ENDIF
            ENDIF
         ENDDO

         IF( DEPTH .LT. 0 )THEN
            EXIT
         ENDIF

      ENDDO

      IF( NPAIRS .LT. (WIU-WIL+1) )THEN
         INFO = 1
      ENDIF


      END SUBROUTINE DLAXRV
*
*************************************************************************
