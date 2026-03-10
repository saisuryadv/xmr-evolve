/* dlaxrv_f77.f -- translated by f2c (version 20240504).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal xtime1, xtime2, xtime3, xddddd;
    integer xnblcks, xnnodes, xmaxdepth, xnumfn, xnumft, xnumgv, xnumgv0, 
	    xnumfs_2__, xnumfs_k__, xnbis_init__, xnbis_cob__, xnbis_iib__, 
	    xnbis_class__, xnbis_sng__, xnbis_clb__, xnbis_wasted__, xnrqi, 
	    xnrqibis, xmaxnrqi, xmaxnrqibis, xnenvgv, xnenvtf, xiiiii;
    logical xstealthmode;
} xmrstats_;

#define xmrstats_1 xmrstats_

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int wsreq_xrv__(integer *n, integer *reqr, integer *reqi)
{
    doublereal a[2];
    integer i__;
    doublereal r__;
    integer t[2];
    extern /* Subroutine */ int dlaxrv_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    integer *);

/*     IMPLICIT NONE */
/*     Convenience routine to perform a workspace query for DLAXRV. */
    *reqr = -1;
    *reqi = -1;
    dlaxrv_(n, a, a, t, t, a, &i__, &i__, &r__, &r__, a, a, &i__, t, a, reqr, 
	    t, reqi, &i__);
    return 0;
} /* wsreq_xrv__ */


/* *********************************************************************** */

/* Subroutine */ int dlaxrv_(integer *n, doublereal *e, doublereal *rootr, 
	integer *rooti, integer *ewl_ae__, doublereal *ewl_lu__, integer *wil,
	 integer *wiu, doublereal *spdiam, doublereal *gaptol, doublereal *w, 
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *rwork, 
	integer *lrwork, integer *iwork, integer *liwork, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    integer rstridei, rstrider, i__, j, k;
    extern /* Subroutine */ int wsreq_xrf__(integer *, integer *, integer *, 
	    integer *, integer *), wsreq_xrx__(integer *, integer *, integer *
	    );
    integer il, iu, wi, wj;
    extern /* Subroutine */ int dlaxrf_env__(integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *);
    doublereal fac;
    integer anc[11], ail[11], aiu[11];
    doublereal tau, piv, lgap;
    integer ilen;
    doublereal prec, ugap;
    integer type__;
    doublereal clgap;
    integer iseed[4];
    /* real tauba; -- removed, was f2c bug, should be taubar */
    doublereal cugap;
    integer depth, ixf77a, ixf77b;
    logical dogrt;
    integer jxrgi, ixrep, jxrep, ixenv;
    extern /* Subroutine */ int dlaxrb_refcls__(integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    integer twist;
    extern /* Subroutine */ int dlaxrb_refsng__(integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), dlaxrb_clssfy__(integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *);
    doublereal lambda;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal anlgap[11], anugap[11];
    integer ixgbnd;
    doublereal taubar;
    extern /* Subroutine */ int dlaxrf_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, logical *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *);
    integer ixcewb, jxcewb;
    extern /* Subroutine */ int dlarnv_(integer *, integer *, integer *, 
	    doublereal *);
    doublereal minenv;
    integer npairs;
    extern /* Subroutine */ int dlaxrr_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    integer ixshfp;
    extern /* Subroutine */ int dlaxrx_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *);
    integer ixfgrt, wsremi, ixprep, xrfsta, wsreqi, ixprte, wsremr, ixwork, 
	    ixroot, jxroot, jxwork, status, wsreqr;
    doublereal ataubar[11];
    integer xrfmode, wsi_xrb__, wsi_xrf__, wsr_xrb__, wsr_xrf__, wsi_xrx__, 
	    wsr_xrx__;

/*     IMPLICIT NONE */




/*  Purpose */
/*  ======= */

/*    Compute all eigenpairs WIL:WIU of the symmetric tridiagonal matrix */
/*    given by (N,E,ROOTR,ROOTI) using the MRRR algorithm. */

/*    This is an internal routine, designed to be called from DSTEXR. */
/*    In particular we expect that the matrix was properly scaled and */
/*    split. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The dimension of the matrix */

/*  E       (input) DOUBLE PRECISION array, dimension ( 1:N-1 ) */
/*          The offdiagonal elements of the matrix, should be positive. */

/*  ROOTR   (input) DOUBLE PRECISION array, dimension (4*N+3) */
/*  ROOTI   (input) INTEGER array, dimension (6+N+N/2) */
/*          Real and integer data used to represent the matrix, cf */
/*          representation.txt */

/*  EWL_AE  (input/output) INTEGER array, dimension (1:2*N) */
/*  EWL_LU  (input/output) DOUBLE PRECISION array, dimension ( 1:2*N ) */
/*          Together they constitute the initial list of eigenvalue */
/*          bounds (ew-list), see ewlist.txt and DLAXRL for details. */
/*             The list must be initialized, but other than that we do */
/*          not specify how good the approximations are - in principle */
/*          there may just be one interval in there (e.g. 0:SPDIAM for */
/*          a p.d. matrix). */
/*          Will be destroyed on exit, that is, the bounds will then not */
/*          be valid for the root representation anymore. */

/*  WIL     (input) INTEGER */
/*  WIU     (input) INTEGER */
/*          The range of wanted eigenpairs, 1 <= WIL <= WIU <= N. */

/*  SPDIAM  (input) DOUBLE PRECISION */
/*          The spectral diameter of the matrix. */

/*  GAPTOL  (input) DOUBLE PRECISION */
/*          Gap Tolerance used to classify eigenvalues into singletons */
/*          and clusters. */

/*  W       (output) DOUBLE PRECISION array, dimension ( WIL:WIU ) */
/*          The computed eigenvalues WIL:WIU. */
/*            If the i'th eigenvector could not be computed (see ISUPPZ), */
/*          then W(i) is set regardless to something that makes sense at */
/*          least with regard to the other eigenvalues */

/*  Z       (output) DOUBLE PRECISION array, dimension ( LDZ, WIL:WIU ) */
/*          The computed orthonormal eigenvectors WIL:WIU. */
/*            If the i'th eigenvector could not be computed (see ISUPPZ), */
/*          then the column Z(:,i) is set zo zero. */

/*  LDZ     (input) INTEGER */
/*          Leading dimension of the array Z. */

/*  ISUPPZ  (output) INTEGER array, dimension ( 2*WIL-1 : 2*WIU ) */
/*          If the i'th vector was succsessfully computed (cf INFO), where */
/*          WIL <= i <= WIU, then the vector Z(i) has nonzero entries */
/*          only within ISUPPZ(2*i-1):ISUPPZ(2*i), which is a subset */
/*          of 1:N. */
/*          Failure to compute the vector is indicated by */
/*          ISUPPZ(2*i-1) = 0, ISUPPZ(2*i) < 0, the precise value of */
/*          the latter may give information about the kind of failure: */
/*          -1:  depth limit exceeded */
/*          -2:  no rep found in DLAXRF */

/*  RWORK   (workspace) DOUBLE PRECISION, dimension (LRWORK) */

/*  LRWORK  (input/output) INTEGER */
/*          Dimension of the array RWORK. Set to -1 for a workspace */
/*          query, then the routine only computes the required real */
/*          workspace, sets LRWORK to this value and returns. */
/*            For a workspace query only N needs to be set, no other */
/*          argument is referenced. You can do a real and integer ws */
/*          query at once. */

/*  IWORK   (workspace) INTEGER, dimension (LIWORK) */

/*  LIWORK  (input/output) INTEGER */
/*          Dimension of the array IWORK. Set to -1 for a workspace */
/*          query, then the routine only computes the required integer */
/*          workspace, sets LIWORK to this value and returns. */
/*            For a workspace query only N needs to be set, no other */
/*          argument is referenced. You can do a real and integer ws */
/*          query at once. */

/*  INFO    (output) INTEGER */
/*          -1 - Not enough workspace given */
/*           0 - Everything ok, all desired eigenpairs were successfully */
/*               computed */
/*           1 - Some eigenpairs were not computed, see ISUPPZ for details */

/*  ====================================================================== */
/* extract -b parameters.inc.f procname=laxrv */

/*     .. Parameters .. */

/*     By how much ULP the shift should be perturbed componentwise always. */
/*     Set to 0 to deactivate. */
/*     By how much ULP the primary data of the father rep should be */
/*     perturbed for a gap retry in DLAXRF. */
/*     Set to 0 to deactivate gap retries. */
/*     Maximal allowed depth of a node in the tree (root has depth 0). */
/*     Note: We need to set this to a small constant since we do a depth */
/*     first traversal, ie, MAXDEPTH influences how much workspace we need. */

/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE WSREQ_XRF(N, ICBEG, ICEND, REQR, REQI) */
/*     IMPLICIT NONE */
/*     INTEGER, INTENT(IN)   ::  N, ICBEG, ICEND */
/*     INTEGER, INTENT(OUT)  ::  REQR, REQI */
/*     END SUBROUTINE WSREQ_XRF */

/* *********************************************************************** */

/*     SUBROUTINE DLAXRF( */
/*    $            N, E, GAPTOL, SPDIAM, DEPTH, TAUBAR, ICBEG, ICEND, */
/*    $            FREPR, FREPI, FLGAP, FUGAP, FEWL_AE, FEWL_LU, FRGINFO, */
/*    $            DOGRT, FREPR_GRT, SHFMOD, ENV, MODE, */
/*    $            SREPR, SREPI, SLGAP, SUGAP, SEWL_AE, SEWL_LU, TAU, */
/*    $            RWORK, LRWORK, IWORK, LIWORK, INFO */
/*    $           ) */
/*     IMPLICIT NONE */
/* note: we want to be able to call this routine both within a df or a */
/* bf tree traversion strategy. For bf, the EWI-List of the son would */
/* be copied over the father afterwards and cached quants discarded, */
/* for df however we can directly use them. */

/*     INTEGER,          INTENT(IN)  ::  N, DEPTH, ICBEG, ICEND, MODE */
/*     INTEGER,          INTENT(IN)  ::  FREPI(6+N+N/2) */
/*     INTEGER,          INTENT(IN)  ::  FRGINFO(ICBEG-1:ICEND) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  GAPTOL, SPDIAM, TAUBAR */
/*     DOUBLE PRECISION, INTENT(IN)  ::  E(N-1), FREPR(4*N+3) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  FREPR_GRT(4*N+3) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  SHFMOD(3*N), ENV(N) */
/*     LOGICAL,          INTENT(IN)  ::  DOGRT */

/*     INTEGER,          INTENT(INOUT)  ::  LRWORK, LIWORK */
/*     INTEGER,          INTENT(INOUT)  ::  IWORK( LIWORK ) */
/*     INTEGER,          INTENT(INOUT)  ::  FEWL_AE(2*ICBEG-1:2*ICEND) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  FLGAP, FUGAP */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  FEWL_LU(2*ICBEG-1:2*ICEND) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( LRWORK ) */

/*     INTEGER,          INTENT(OUT)  ::  INFO */
/*     INTEGER,          INTENT(OUT)  ::  SREPI(6+N+N/2) */
/*     INTEGER,          INTENT(OUT)  ::  SEWL_AE(2*ICBEG-1:2*ICEND) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  SLGAP, SUGAP, TAU */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  SREPR(4*N+3) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  SEWL_LU(2*ICBEG-1:2*ICEND) */
/*     END SUBROUTINE DLAXRF */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE WSREQ_XRX(N,REQR,REQI) */
/*     IMPLICIT NONE */
/*     INTEGER, INTENT(IN)     ::  N */
/*     INTEGER, INTENT(INOUT)  ::  REQR, REQI */
/*     END SUBROUTINE WSREQ_XRX */

/* *********************************************************************** */

/*     SUBROUTINE DLAXRX( */
/*    $             N, E, REPR, REPI, */
/*    $             INDEX, MINGAP, LAMBDA, */
/*    $             Z, ISUPPZ, RWORK, LRWORK, IWORK, LIWORK */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)     ::  N, INDEX */
/*     INTEGER,          INTENT(IN)     ::  LRWORK, LIWORK */
/*     INTEGER,          INTENT(IN)     ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)     ::  MINGAP */
/*     DOUBLE PRECISION, INTENT(IN)     ::  E(N-1), REPR(4*N+3) */

/*     INTEGER,          INTENT(INOUT)  ::  IWORK( LIWORK ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( LRWORK ), LAMBDA */

/*     INTEGER,          INTENT(OUT)    ::  ISUPPZ(2) */
/*     DOUBLE PRECISION, INTENT(OUT)    ::  Z(N) */
/*     END SUBROUTINE DLAXRX */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRB_CLSSFY( */
/*    $             N, REPR, REPI, DEPTH, SPDIAM, */
/*    $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, */
/*    $             RGINFO, */
/*    $             WIL, WIU, GIL, GIU, GAPTOL, */
/*    $             RWORK, IWORK */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)     ::  N, DEPTH, IL, IU */
/*     INTEGER,          INTENT(IN)     ::  WIL, WIU, GIL, GIU */
/*     INTEGER,          INTENT(IN)     ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)     ::  SPDIAM, GAPTOL */
/*     DOUBLE PRECISION, INTENT(IN)     ::  REPR(4*N+3) */

/*     INTEGER,          INTENT(INOUT)  ::  EWL_AE( 2*IL-1 : 2*IU ) */
/*     INTEGER,          INTENT(INOUT)  ::  IWORK( 2*N + 4*(IU-IL+1) ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU( 2*IL-1 : 2*IU ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( IU - IL + 1 ) */

/*     INTEGER,          INTENT(OUT)    ::  RGINFO( IL : IU-1 ) */
/*     END SUBROUTINE DLAXRB_CLSSFY */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRB_REFSNG( */
/*    $             N, REPR, REPI, DEPTH, */
/*    $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, RGINFO, */
/*    $             RWORK, IWORK */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)     ::  N, IL, IU, DEPTH */
/*     INTEGER,          INTENT(IN)     ::  RGINFO( IL : IU-1 ) */
/*     INTEGER,          INTENT(IN)     ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)     ::  REPR(4*N+3) */

/*     INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU) */
/*     INTEGER,          INTENT(INOUT)  ::  IWORK( 2*(IU-IL+1) ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( 2*(IU-IL+1) ) */

/*     END SUBROUTINE DLAXRB_REFSNG */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRB_REFCLS( */
/*    $             N, REPR, REPI, DEPTH, */
/*    $             ICBEG, ICEND, LGAP, UGAP, EWL_AE, EWL_LU, */
/*    $             RGINFO, */
/*    $             RWORK, IWORK */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)     ::  N, DEPTH, ICBEG, ICEND */
/*     INTEGER,          INTENT(IN)     ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)     ::  REPR(4*N+3) */

/*     INTEGER,          INTENT(INOUT)  ::  EWL_AE( 2*ICBEG-1 : 2*ICEND ) */
/*     INTEGER,          INTENT(INOUT)  ::  IWORK( 4*(ICEND-ICBEG+1) ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU( 2*ICBEG-1 : 2*ICEND ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( (ICEND-ICBEG+1) ) */

/*     INTEGER,          INTENT(INOUT)  ::  RGINFO( ICBEG : ICEND-1 ) */
/*     END SUBROUTINE DLAXRB_REFCLS */
/*     END INTERFACE */

/*     .. Constants .. */


/*     .. Local Variables .. */


/*           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*           !!!! SYNCHRONIZE ANY CHANGES HERE WITH xmr.h !!!!! */
/*           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*     Accumulative part: All these are accumulated over multiple */
/*     calls of stexr/laxrv. */
/*        Times for stages 1-3 (before/during/after block loop) */
/*        Marker to verify data integrity */

/*        Number of calls to dlaxrn, dlaxrt, dlaxrg, and */
/*        dlaxrs (all twists/single twist) */
/*        Number of bisection steps in various parts. */
/*        Count/max rqi steps */
/*        For computing envelopes, number of full twisted factos */
/*        and gv-calls */
/*        Marker to verify data integrity */
/*    Temporaries and Configuration */
/*        If true, calls to internal routines dlaxr[n,t,s] are not */
/*        counted. Higher level routines like DLAXRB are not affected. */
/*        We need this to hide bisections steps done within aux */
/*        routines like xmr_estrc. */

/*  ===== Executable Statements ========================================== */

/*     Number of real/integer data items for a representation */
    /* Parameter adjustments */
    --rooti;
    --rootr;
    --e;
    --ewl_ae__;
    --ewl_lu__;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --rwork;
    --iwork;

    /* Function Body */
    rstrider = (*n << 2) + 3;
    rstridei = *n + 6 + *n / 2;
/*     ---------------------- */
/*      Workspace Management */
/*     ---------------------- */
    ixwork = 1;
    jxwork = 1;
/*     Right gap infos, indexed 0:N */
    jxrgi = jxwork;
    jxwork = jxwork + *n + 1;
/*     In CEWB we buffer the ewl for the child as given by DLAXRF, */
/*     before copying it again over EWL. */
    ixcewb = ixwork;
    jxcewb = jxwork;
    ixwork += *n << 1;
    jxwork += *n << 1;
/*     Perturbation factors for the shift, we use three different */
/*     sets, one for shifts to the left, one identical to 1 for np */
/*     perturbation, and one for shifts to the right, cf DLAXRF. */
    ixshfp = ixwork;
    ixwork += *n * 3;
/*     Gap boundaries from the times when the gaps were classified. */
/*     gcbnd(2*j) and gcbnd(2*j+1) constitute gap j:j+1 */
/*     We define the index IXGBND such that */
/*        RWORK(IXGBND + k) = gcbnd(k),  2 <= k <= 2N-1, */
/*     that is, we start with gcbnd(2), which is to hold the bound right */
/*     of ew 1. */
    ixgbnd = ixwork - 2;
    ixwork = ixwork + (*n << 1) - 2;
/*     The representation stack. */
    jxroot = jxwork;
    ixroot = ixwork;
    jxwork += rstridei * 11;
    ixwork += rstrider * 11;
/*     The envelope for the current representation */
    ixenv = ixwork;
    ixwork += *n;
/*     Data for perturbed representations to use for gap retries. */
/*       FGRT: Perturbation factors for the primary data of the current */
/*             representation. */
/*       PRTE: Perturbed versions of the offdiagonal entries. */
/*       PREP: Buffer to hold the perturbed representation data */
/*             (real only). */
    dogrt = TRUE_;
    if (dogrt) {
	ixfgrt = ixwork;
	ixprte = ixwork + *n;
	ixprep = ixwork + (*n << 1) - 1;
	ixwork = ixwork + (*n << 1) - 1 + rstrider;
    } else {
/*        buffers won't be accessed */
	ixfgrt = 1;
	ixprte = 1;
	ixprep = 1;
    }

/*     -------------------------- */
/*      Handle workspace queries */
/*     -------------------------- */

    wsremi = *liwork - jxwork + 1;
    wsremr = *lrwork - ixwork + 1;
/*     Combine ws for dlaxrb_[clssfy|refsng|refcls] */
    wsr_xrb__ = *n;
    wsi_xrb__ = *n * 7;
    wsreq_xrx__(n, &wsr_xrx__, &wsi_xrx__);
    wsreq_xrf__(n, &c__1, n, &wsr_xrf__, &wsi_xrf__);
/* Computing MAX */
    i__1 = max(wsi_xrb__,wsi_xrx__);
    wsreqi = jxwork - 1 + max(i__1,wsi_xrf__);
/* Computing MAX */
    i__1 = max(wsr_xrb__,wsr_xrx__);
    wsreqr = ixwork - 1 + max(i__1,wsr_xrf__);
    if (*liwork == -1 || *lrwork == -1) {
	if (*liwork == -1) {
	    *liwork = wsreqi;
	}
	if (*lrwork == -1) {
	    *lrwork = wsreqr;
	}
	return 0;
    }
    if (wsreqi > *liwork || wsreqr > *lrwork) {
	*info = -1;
	return 0;
    }

/*     ------------------------------------------------------------------- */

    prec = dlamch_("Precision", (ftnlen)9);

/*     .. Prepare shift random perturbation .. */

    iseed[0] = 29;
    iseed[1] = 6;
    iseed[2] = 1978;
    iseed[3] = 31;
    if (TRUE_) {
	i__1 = ixshfp + *n * 3 - 1;
	for (ixf77a = ixshfp; ixf77a <= i__1; ++ixf77a) {
	    rwork[ixf77a] = 1.;
/* L90001: */
	}
/* L90002: */
    } else {
	fac = 0.;
	dlarnv_(&c__1, iseed, n, &rwork[ixshfp]);
	i__1 = *n - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    rwork[ixshfp + i__] = 1. - rwork[ixshfp + i__] * fac;
/* L90003: */
	}
/* L90004: */
	i__1 = ixshfp + (*n << 1) - 1;
	for (ixf77a = ixshfp + *n; ixf77a <= i__1; ++ixf77a) {
	    rwork[ixf77a] = 1.;
/* L90005: */
	}
/* L90006: */
	dlarnv_(&c__1, iseed, n, &rwork[ixshfp + (*n << 1)]);
	i__1 = *n * 3 - 1;
	for (i__ = *n << 1; i__ <= i__1; ++i__) {
	    rwork[ixshfp + i__] = rwork[ixshfp + i__] * fac + 1.;
/* L90007: */
	}
/* L90008: */
    }

/*     .. Prepare perturbation factors for the representation data .. */

    if (dogrt) {
	fac = prec * 16;
	dlarnv_(&c__2, iseed, n, &rwork[ixfgrt]);
	i__1 = *n - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    rwork[ixfgrt + i__] = rwork[ixfgrt + i__] * fac + 1.;
/* L90009: */
	}
/* L90010: */
	i__1 = *n - 1;
	dlarnv_(&c__2, iseed, &i__1, &rwork[ixprte]);
	i__1 = *n - 2;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    rwork[ixprte + i__] = e[i__ + 1] * (rwork[ixprte + i__] * fac + 
		    1.);
/* L90011: */
	}
/* L90012: */
    }
    *info = 0;
    npairs = 0;
    i__1 = *wiu - *wil + 1 << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	isuppz[i__] = 0;
/* L90013: */
    }
/* L90014: */
/*     **************************************************************** */
/*     *                        MAIN LOOP                             * */
/*     *                       -----------                            * */
/*     *  Perform a depth-first traversion of the representation tree * */
/*     **************************************************************** */

    iwork[jxrgi] = 2;
    iwork[jxrgi + *n] = 2;
    i__1 = rstrider;
    for (ixf77a = 1; ixf77a <= i__1; ++ixf77a) {
	rwork[ixroot + ixf77a - 1] = rootr[ixf77a];
/* L90015: */
    }
/* L90016: */
    i__1 = rstridei;
    for (ixf77a = 1; ixf77a <= i__1; ++ixf77a) {
	iwork[jxroot + ixf77a - 1] = rooti[ixf77a];
/* L90017: */
    }
/* L90018: */
    ail[0] = 1;
    aiu[0] = *n;
    anc[0] = 1;
    ataubar[0] = 0.;
/*     We do know not much about how accurate the initial bounds are. */
/*     The near-zero bound only needs to be set to sth > 0, as dlaxrb */
/*     will update it automatically to at lest |ew|. */
/* Computing MAX */
    d__1 = *spdiam, d__2 = abs(ewl_lu__[1]);
    anlgap[0] = max(d__1,d__2);
/* Computing MAX */
    d__2 = *spdiam, d__3 = (d__1 = ewl_lu__[*n * 2], abs(d__1));
    anugap[0] = max(d__2,d__3);
    depth = 0;
L90019:
    il = ail[depth];
    iu = aiu[depth];
    ilen = iu - il + 1;
    taubar = ataubar[depth];
    ixrep = ixroot + depth * rstrider;
    jxrep = jxroot + depth * rstridei;
    if (anc[depth] == il) {

/*           Classify eigenvalues */

/*            CALL DLAXRB( */
/*     $             N, RWORK(IXREP), IWORK(JXREP), DEPTH, SPDIAM, */
/*     $             IL, IU, ANLGAP(DEPTH), ANUGAP(DEPTH), */
/*     $             EWL_AE(2*IL-1), EWL_LU(2*IL-1), */
/*     $             IWORK(JXRGI + IL), RWORK(IXGBND + 2*IL), */
/*     $             WIL, WIU, GAPTOL, .TRUE., */
/*     $             RWORK(IXWORK), WSR_XRB, IWORK(JXWORK), WSI_XRB */
/*     $           ) */
	dlaxrb_clssfy__(n, &rwork[ixrep], &iwork[jxrep], &depth, spdiam, &il, 
		&iu, &anlgap[depth], &anugap[depth], &ewl_ae__[(il << 1) - 1],
		 &ewl_lu__[(il << 1) - 1], &iwork[jxrgi + il], wil, wiu, &
		c__1, n, gaptol, &rwork[ixwork], &iwork[jxwork]);
/*           Record bounds at time of classification */
	i__1 = iu - 1;
	for (ixf77a = il; ixf77a <= i__1; ++ixf77a) {
	    rwork[ixgbnd + (ixf77a << 1)] = ewl_lu__[ixf77a * 2];
/* L90021: */
	}
/* L90022: */
/*           Initial refinement of singletons */
	dlaxrb_refsng__(n, &rwork[ixrep], &iwork[jxrep], &depth, &il, &iu, &
		anlgap[depth], &anugap[depth], &ewl_ae__[(il << 1) - 1], &
		ewl_lu__[(il << 1) - 1], &iwork[jxrgi + il], &rwork[ixwork], &
		iwork[jxwork]);
	++xmrstats_1.xnnodes;
	xmrstats_1.xmaxdepth = max(xmrstats_1.xmaxdepth,depth);
    }

/*        Handle children, singletons right here, clusters onto the */
/*        stack. */

L90023:
    if (anc[depth] > iu) {
/*              backtrack */
	--depth;
	goto L90024;
    }
    i__ = anc[depth];
    j = i__;
L90025:
    if (iwork[jxrgi + j] == 2) {
	goto L90026;
    }
    ++j;
    goto L90025;
L90026:
    anc[depth] = j + 1;

/*           Prune children that do not contain wanted eigenpairs */

    if (j < *wil || i__ > *wiu) {
	goto L90023;
    }
    wi = max(*wil,i__);
    wj = min(*wiu,j);
/*           Determine left and right gaps for this child. Here, using */
/*           for the far bounds the ones from when the gap was first */
/*           identified leads to conformal subtrees. */
    if (i__ == il) {
	lgap = anlgap[depth];
    } else {
	lgap = ewl_lu__[(i__ << 1) - 1] - rwork[ixgbnd + (i__ - 1 << 1)];
    }
    if (j == iu) {
	ugap = anugap[depth];
    } else {
	ugap = rwork[ixgbnd + (j << 1) + 1] - ewl_lu__[j * 2];
    }
    if (i__ == j) {
/*              ----------- */
/*               Singleton */
/*              ----------- */

	lambda = (ewl_lu__[(i__ << 1) - 1] + ewl_lu__[i__ * 2]) * .5;
	d__1 = min(lgap,ugap);
	dlaxrx_(n, &e[1], &rwork[ixrep], &iwork[jxrep], &i__, &d__1, &lambda,
		&z__[i__ * z_dim1 + 1], &isuppz[(i__ << 1) - 1], &rwork[
		ixwork], &wsr_xrx__, &iwork[jxwork], &wsi_xrx__);
	w[i__] = lambda + taubar;
	++npairs;
    } else {
/*              --------- */
/*               Cluster */
/*              --------- */
	status = 0;
	tau = 0.;
	if (depth == 10) {
	    status = -1;
	} else {
/*                 1 means try outside and also exploit subgaps */
	    xrfmode = 1;
/*                  XRFMODE = 2 */
/*                 Refine the cluster */
	    dlaxrb_refcls__(n, &rwork[ixrep], &iwork[jxrep], &depth, &i__, &j,
		     &lgap, &ugap, &ewl_ae__[(i__ << 1) - 1], &ewl_lu__[(i__ 
		    << 1) - 1], &iwork[jxrgi + i__], &rwork[ixwork], &iwork[
		    jxwork]);
/*                 Compute Envelope */
	    dlaxrf_env__(n, &e[1], &rwork[ixrep], &iwork[jxrep], &depth, &i__,
		     &j, &lgap, &ugap, &ewl_ae__[(i__ << 1) - 1], &ewl_lu__[(
		    i__ << 1) - 1], &iwork[jxrgi + i__ - 1], &xrfmode, &rwork[
		    ixenv], &minenv, &rwork[ixwork], &iwork[jxwork]);
/*                 Setup perturbed rep for gap retries */
/*                 Note: Due to the depth first approach, we have to do */
/*                 this here, possibly multiple  times for the same rep. */
	    if (dogrt) {
		i__1 = *n - 1;
		for (k = 0; k <= i__1; ++k) {
		    rwork[ixprep + k] = rwork[ixrep + k] * rwork[ixfgrt + k];
/* L90027: */
		}
/* L90028: */
		piv = rwork[ixrep - 1 + ((*n << 2) + 3)];
		twist = iwork[jxrep + 1];
		type__ = iwork[jxrep];
		dlaxrr_(n, &twist, &type__, &rwork[ixprte], &piv, &rwork[
			ixprep], &iwork[jxrep]);
	    }
	    dlaxrf_(n, &e[1], gaptol, spdiam, &depth, &taubar, &i__, &j, &
		    rwork[ixrep], &iwork[jxrep], &lgap, &ugap, &ewl_ae__[(i__ 
		    << 1) - 1], &ewl_lu__[(i__ << 1) - 1], &iwork[jxrgi + i__ 
		    - 1], &dogrt, &rwork[ixprep], &rwork[ixshfp], &rwork[
		    ixenv], &xrfmode, &rwork[ixrep + rstrider], &iwork[jxrep 
		    + rstridei], &clgap, &cugap, &iwork[jxcewb], &rwork[
		    ixcewb], &tau, &rwork[ixwork], &wsr_xrf__, &iwork[jxwork],
		     &wsi_xrf__, &xrfsta);
	    if (xrfsta != 0) {
/*                    Could not find a child representation */
		status = -2;
	    } else {
		i__1 = (j - i__ + 1 << 1) - 1;
		for (k = 0; k <= i__1; ++k) {
		    ewl_lu__[(i__ << 1) - 1 + k] = rwork[ixcewb + k];
		    ewl_ae__[(i__ << 1) - 1 + k] = iwork[jxcewb + k];
/* L90029: */
		}
/* L90030: */
	    }
	}
	if (status != 0) {
	    i__1 = wj << 1;
	    for (ixf77a = wi << 1; ixf77a <= i__1; ixf77a += 2) {
		isuppz[ixf77a] = status;
/* L90031: */
	    }
/* L90032: */
	    i__1 = *n;
	    for (ixf77a = 1; ixf77a <= i__1; ++ixf77a) {
		i__2 = wj;
		for (ixf77b = wi; ixf77b <= i__2; ++ixf77b) {
		    z__[ixf77a + ixf77b * z_dim1] = 0.;
/* L90035: */
		}
/* L90036: */
/* L90033: */
	    }
/* L90034: */
	    i__1 = wj;
	    for (ixf77a = wi; ixf77a <= i__1; ++ixf77a) {
		w[ixf77a] = (ewl_lu__[(i__ << 1) - 1] + ewl_lu__[j * 2]) / 2
			+ taubar;
/* L90037: */
	    }
/* L90038: */
	} else {
	    ++depth;
	    ataubar[depth] = taubar + tau;
	    anlgap[depth] = clgap;
	    anugap[depth] = cugap;
	    ail[depth] = i__;
	    aiu[depth] = j;
	    anc[depth] = i__;
	    goto L90024;
	}
    }
    goto L90023;
L90024:
    if (depth < 0) {
	goto L90020;
    }
    goto L90019;
L90020:
    if (npairs < *wiu - *wil + 1) {
	*info = 1;
    }
    return 0;
} /* dlaxrv_ */

