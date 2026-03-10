/* dlaxrf_f77.f -- translated by f2c (version 20240504).
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
static integer c__4 = 4;

/* Subroutine */ int wsreq_xrf__(integer *n, integer *icbeg, integer *icend, 
	integer *reqr, integer *reqi)
{
    doublereal a[2];
    integer i__;
    logical l;
    doublereal r__;
    integer t[2];
    extern /* Subroutine */ int dlaxrf_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, logical *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *);

/*     IMPLICIT NONE */
    *reqr = -1;
    *reqi = -1;
    dlaxrf_(n, a, &r__, &r__, &i__, &r__, icbeg, icend, a, t, &r__, &r__, t, 
	    a, t, &l, a, a, a, &c__1, a, t, &r__, &r__, t, a, &r__, a, reqr, 
	    t, reqi, &i__);
    return 0;
} /* wsreq_xrf__ */


/* *********************************************************************** */

/* Subroutine */ int dlaxrf_(integer *n, doublereal *e, doublereal *gaptol, 
	doublereal *spdiam, integer *depth, doublereal *taubar, integer *
	icbeg, integer *icend, doublereal *frepr, integer *frepi, doublereal *
	flgap, doublereal *fugap, integer *fewl_ae__, doublereal *fewl_lu__, 
	integer *frginfo, logical *dogrt, doublereal *frepr_grt__, doublereal 
	*shfmod, doublereal *env, integer *mode, doublereal *srepr, integer *
	srepi, doublereal *slgap, doublereal *sugap, integer *sewl_ae__, 
	doublereal *sewl_lu__, doublereal *tau, doublereal *rwork, integer *
	lrwork, integer *iwork, integer *liwork, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;
    logical L__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer i_sign(integer *, integer *), pow_ii(integer *, integer *);

    /* Local variables */
    integer maxncand, nibweval;
    logical lasthope;
    integer i__, j, nibactive;
    doublereal bbesteval;
    integer bistrycnt, kl;
    doublereal xa, xb, xc, xd;
    integer ku, xi;
    extern /* Subroutine */ int dlaxrf_cob__(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *), 
	    dlaxrf_iib__(integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *);
    logical libnoneval, isgapretry;
    doublereal clb, cub, mid;
    integer loc, dir;
    doublereal eps;
    integer jxg, bbeg, bend;
    doublereal eval, prec;
    integer ixdp, ixrp;
    extern /* Subroutine */ int dlaxrf_seltw__(integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    integer abend[9], icand, ncand, iclen, maxnb, ixf77a, twist;
    extern /* Subroutine */ int dlaxrl_refine__(integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *);
    doublereal wntrw;
    extern /* Subroutine */ int dlaxrf_selshf__(integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *);
    integer indfbb;
    extern doublereal dlamch_(char *, ftnlen);
    integer ibatch, nbatch;
    logical cont_a__, cont_b__;
    integer iycfdg;
    logical cont_c__;
    integer nfudge, cobsta;
    logical dokill;
    integer iycloc;
    doublereal rclwid;
    integer ixgmap, jyomga, shpoff, nfbtrf, ixcevl, maxcpo, maxcpi, ixctau, 
	    shpsel;
    logical inerok;
    integer iyomdp, iycgrt, ixblup;
    extern /* Subroutine */ int dlaxrs_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *), dlaxrr_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    integer kminus, iyctwi;
    extern integer dlaxrn_(integer *, doublereal *, integer *, doublereal *
	    );
    integer wsreqi, iyomrp, iycxiz, xizero, wsreqr, ixwork, iywork, iytwok;
    logical hadeval;
    doublereal pivbase;
    logical success;

/*     IMPLICIT NONE */
/* note: we want to be able to call this routine both within a df or a */
/* bf tree traversion strategy. For bf, the EWI-List of the son would */
/* be copied over the father afterwards and cached quants discarded, */
/* for df however we can directly use them. */



/* # */
/*  Purpose */
/*  ======= */

/*     Find a shift such that the computed child representation is an RRR */
/*     (if possible). */

/*     Modes */
/*      1 - try both outside and inside locations */
/*      2 - place only shifts on the outside, or more precisely, that */
/*          depend on the outer bounds */
/*      3 - place only true inside shifts, that is, exactly those */
/*          not done in mode 2. */


/*     Return status */
/*      INFO=0  we found a shift TAU and the sons ewl-list and repdata are */
/*              set */

/*     For INFO != 0, we did not find an acceptable rep */
/*      INFO=1  no candidate found in select shift */
/*      INFO=2  none of the candidates we tried was acceptable */

/*  Arguments */
/*  ========= */

/*  SHFMOD   (input) DOUBLE PRECISION array, dimension(3*N) */
/*           Perturbations (1*alpha) used to modify the shift, with the */
/*           intent of making emergence of gaps in children more likely. */
/*           Three sets should be provided: */
/*             SHFMOD(1:N)     <= 1, for shifting left of a bound */
/*             SHFMOD(N+1:2N)   = 1, for no shift perturbation */
/*             SHFMOD(2N+1:3N) >= 1, for shifting right of a bound */
/*           These can all be one if perturbing the shift is deactivated. */

/*  DOGRT    (input) LOGICAL */
/*           Specifies if for children without a gap or at least wide */
/*           enough local spectrum, a gap retry is to be initiated, using */
/*           the same shift but a perturbed father representation. */

/*  FREPR_GRT  (input) DOUBLE PRECISION array, dimension( (4*N+3) ) */
/*           Perturbed real data of the father rep, to be used for gap */
/*           retrues. */
/*           Only accessed if DOGRT = .TRUE. */

/*  LRWORK   (input/output) INTEGER */
/*  LIWORK   (input/output) INTEGER */
/*           The lengths of the workspace arrays RWORK and IWORK, */
/*           respectively. They support a workspace query, if one or both */
/*           of them is -1, the routine just computes the required */
/*           workspace and sets LRWORK and LIWORK to it. */
/*           A ws query needs N, ICBEG, ICEND and MODE to be set, */
/*           no other argument is referenced. */
/*           Recommended: Don't use this directly, call WSREQ_XRF. */

/*  ====================================================================== */

/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRS( */
/*    $             N, IA, IE, REPR, REPI, */
/*    $             TAU, SHFPRT, */
/*    $             DPLUS, OMGADP, RPLUS, OMGARP, GAMMAP, TWISTOK, */
/*    $             RWORK */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, IA, IE */
/*     INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  TAU */
/*     DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), SHFPRT(N) */

/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK(2*N) */

/*     INTEGER,          INTENT(OUT)  ::  OMGADP(N), OMGARP(N) */
/*     INTEGER,          INTENT(OUT)  ::  TWISTOK(IA:IE) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  DPLUS(1:N-1), RPLUS(2:N) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  GAMMAP(IA:IE) */
/*     END SUBROUTINE DLAXRS */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRF_SELSHF( */
/*    $             N, REPR, REPI, */
/*    $             ICBEG, ICEND, LGAP, UGAP, EWL_AE, EWL_LU, RGINFO, */
/*    $             TAUBAR, GAPTOL, MAXCPO, MAXCPI, MAXNC, MAXNB, */
/*    $             NCAND, ACLOC, ACTAU, NBATCH, ABEND, */
/*    $             IWORK */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, ICBEG, ICEND */
/*     INTEGER,          INTENT(IN)  ::  MAXCPO, MAXCPI, MAXNC, MAXNB */
/*     INTEGER,          INTENT(IN)  ::  RGINFO(ICBEG-1:ICEND) */
/*     INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  TAUBAR, GAPTOL */

/*     INTEGER,          INTENT(INOUT) :: EWL_AE(2*ICBEG-1:2*ICEND) */
/*     INTEGER,          INTENT(INOUT) :: IWORK( 2*(ICEND-ICBEG+1) ) */
/*     DOUBLE PRECISION, INTENT(INOUT) :: LGAP, UGAP */
/*     DOUBLE PRECISION, INTENT(INOUT) :: EWL_LU(2*ICBEG-1:2*ICEND) */

/*     INTEGER,          INTENT(OUT) :: NCAND, NBATCH */
/*     INTEGER,          INTENT(OUT) :: ACLOC( MAXNC ) */
/*     INTEGER,          INTENT(OUT) :: ABEND( MAXNB ) */
/*     DOUBLE PRECISION, INTENT(OUT) :: ACTAU( MAXNC ) */
/*     END SUBROUTINE DLAXRF_SELSHF */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRF_SELTW( */
/*    $             N, KMINUS, D, OMEGAD, R, OMEGAR, E, GAMMA, TWISTOK, */
/*    $             ENV, SPDIAM, MINGAP, */
/*    $             DIR, LOC, TAU, WINOK, */
/*    $             K, EVAL, XIZERO, */
/*    $             RWORK, IWORK */
/*    $           ) */
/*     IMPLICIT NONE */

/*     LOGICAL,          INTENT(IN)     ::  WINOK */
/*     INTEGER,          INTENT(IN)     ::  N, KMINUS, DIR, LOC */
/*     INTEGER,          INTENT(IN)     ::  OMEGAD(N), OMEGAR(N) */
/*     INTEGER,          INTENT(IN)     ::  TWISTOK( N ) */
/*     DOUBLE PRECISION, INTENT(IN)     ::  SPDIAM, MINGAP, TAU */
/*     DOUBLE PRECISION, INTENT(IN)     ::  D(1:N-1), R(2:N), E(N-1) */
/*     DOUBLE PRECISION, INTENT(IN)     ::  ENV(N) */

/*     INTEGER,          INTENT(INOUT)  ::  IWORK( N+2 ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( 5*N ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  GAMMA(N) */

/*     INTEGER,          INTENT(OUT)    ::  K, XIZERO */
/*     DOUBLE PRECISION, INTENT(OUT)    ::  EVAL */
/*     END SUBROUTINE DLAXRf_SELTW */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRF_COB( */
/*    $             N, FLGAP, FLB, FUB, FUGAP, REPR, REPI, TAU, GAPTOL, */
/*    $             ICBEG, ICEND, XIZERO, SLGAP, SUGAP, EWL_AE, EWL_LU, */
/*    $             STATUS */
/*    $         ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, ICBEG, ICEND, XIZERO */
/*     INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  FLGAP, FLB, FUB, FUGAP */
/*     DOUBLE PRECISION, INTENT(IN)  ::  TAU, GAPTOL */
/*     DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3) */

/*     INTEGER,          INTENT(OUT)  ::  STATUS */
/*     INTEGER,          INTENT(OUT)  ::  EWL_AE(2*ICBEG-1:2*ICEND) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  SLGAP, SUGAP */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  EWL_LU(2*ICBEG-1:2*ICEND) */
/*     END SUBROUTINE DLAXRF_COB */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRF_IIB( */
/*    $             N, REPR, REPI, TAU, ICBEG, ICEND, */
/*    $             FEWL_AE, FEWL_LU, */
/*    $             SLGAP, SUGAP, SEWL_AE, SEWL_LU */
/*    $         ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, ICBEG, ICEND */
/*     INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2) */
/*     INTEGER,          INTENT(IN)  ::  FEWL_AE(2*ICBEG-1:2*ICEND) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  FEWL_LU(2*ICBEG-1:2*ICEND) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  TAU */
/*     DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3) */

/*     INTEGER,          INTENT(INOUT)  ::  SEWL_AE(2*ICBEG-1:2*ICEND) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  SLGAP, SUGAP */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  SEWL_LU(2*ICBEG-1:2*ICEND) */
/*     END SUBROUTINE DLAXRF_IIB */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRR( N, K, TYPE, E, PIVBASE, REPR, REPI ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, K, TYPE */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVBASE */
/*     DOUBLE PRECISION, INTENT(IN)  ::  E(1:N-1) */

/*     INTEGER,          INTENT(INOUT)  ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  REPR(4*N+3) */
/*     END SUBROUTINE DLAXRR */
/*     END INTERFACE */
/*     INTERFACE */
/*     FUNCTION DLAXRN(N, REPR, REPI, TAU) */
/*     IMPLICIT NONE */
/*     INTEGER  ::  DLAXRN */
/*     INTEGER,          INTENT(IN)  ::  N */
/*     INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), TAU */
/*     END FUNCTION DLAXRN */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRL_REFINE( */
/*    $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, I, J, LAMBDA, XI */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER         , INTENT(IN)  ::  IL, IU, I, J, XI */
/*     DOUBLE PRECISION, INTENT(IN)  ::  LAMBDA */

/*     INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU) */
/*     END SUBROUTINE DLAXRL_REFINE */
/*     END INTERFACE */

/*     .. Constants .. */


/*     .. Parameters .. */


/*     .. Parameters .. */

/*     The evaluation determines a real value EVAL for each candidate. */
/*     Based on this we categorize them broadly as follows */
/*        EVAL <= 1        passes a priori, will take it directly */
/*        EVAL <= EVALOK   is ok */
/*        EVAL >  EVALOK   only take as last resort */
/*     To transform a father bound between bisection and shifted zero */
/*     inertia semantics, we allow up to MAXNFBTRF tries, where the */
/*     relative backoff at the K'th try is increased by */
/*       FBTRFGRAN**K * N*EPS, */
/*     thus */
/*     For eigenvalue underflow with a consistent negcount, fudge the */
/*     shift by |tau|*FUDGEFAC*Eps away in the right direction and try */
/*     again. */
/*     Set either one to 0 to deactivate fudging. */
/*     If true, the eigenvalues of a candidate that has passed all other */
/*     tests are refined until we can see that the cluster is wide enough. */
/*     If that is not possible, then either a gap retry is initiated, if */
/*     that is activated, or the candidate is discarded. */
/*     Separate max number of candidates for inside locations, effective */
/*     will be minimum of this and MAXCPL */

/*     .. Constants .. */


/*     .. Locals .. */


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

    /* Parameter adjustments */
    --srepi;
    --srepr;
    --env;
    --shfmod;
    --frepr_grt__;
    --frepi;
    --frepr;
    --e;
    --fewl_ae__;
    --fewl_lu__;
    --frginfo;
    --sewl_ae__;
    --sewl_lu__;
    --rwork;
    --iwork;

    /* Function Body */
    iclen = *icend - *icbeg + 1;
    kminus = frepi[2];
    maxcpo = 3;
    maxcpi = 2;
    if (*mode == 2) {
	maxcpi = 0;
    } else if (*mode == 3) {
	maxcpo = 0;
    }
/* One idea may be to limit the backoffs for deeper levels, to prefer close */
/* ones to slightly better (but still bad) a priori criteria. */
/* Better solution:  Introduce relative offset from base bound into eval, maybe */
/*   like eval  <-  eval^2 * reloff */
/*      IF( DEPTH.GT.0 )THEN */
/*         MAXCPO = MIN( MAXCPO, 1 ) */
/*         MAXCPI = MIN( MAXCPI, 1 ) */
/*      ENDIF */
/* C~ */
    maxncand = (maxcpo << 1) + (*icend - *icbeg << 1) * maxcpi + 1;
    maxnb = max(maxcpo,maxcpi) * 3;

/*     ---------------------- */
/*      Workspace Allocation */
/*     ---------------------- */

    ixwork = 1;
    iywork = 1;
/*     Buffer arrays for DLAXRS */
/*       DP    real  dimension(N-1) */
/*       RP    real  dimension(N-1) */
/*       GMAP  real  dimension(N) */
/*       TWOK  int   dimension(N) */
    ixdp = ixwork;
    ixrp = ixdp + *n - 1;
    ixgmap = ixrp + *n - 1;
    ixwork = ixwork + *n * 3 - 2;
    iytwok = iywork;
    iyomdp = iytwok + *n;
    iyomrp = iyomdp + *n;
    iywork += *n * 3;
/*     Data per candidate (NCAND <= 2*ILEN), indexed from 0 */
/*     Set by SELSHF: */
/*       CLOC (int)   ew index we shift to, sign gives direction */
/*       COFF (real)  absolute offset from the ew bound */
/*     Set after evaluation: */
/*       CTWI (int)   optimal twist index */
/*       CXIZ (int)   zero inertia */
/*       CEVL (real)  evaluation */
/*     General: */
/*       CFDG (int)   counter for number of small fudges to avoid */
/*                    ew-underflow */
/*       CGRT (int)   set to one for a retry with the perturbed rep */
/*                    to get a wider cluster or even a gap */
    iycloc = iywork;
    iyctwi = iywork + maxncand;
    iycxiz = iywork + (maxncand << 1);
    iycfdg = iywork + maxncand * 3;
    iycgrt = iywork + (maxncand << 2);
    iywork += maxncand * 5;
    ixctau = ixwork;
    ixcevl = ixctau + maxncand;
    ixwork = ixcevl + maxncand;
/*     Blowup of father bounds to get consistent zero inertias. */
/*     Indexed (2*ICBEG-1:2*ICEND), conformingly with EWL. */
    ixblup = ixwork - ((*icbeg << 1) - 1);
    ixwork += iclen << 1;

/*     .. Handle Workspace Query .. */

/*      Additional ws for called routines: */
/*       2N real  for DLAXRS */
/*       5N real N+2 int  for _SELTW */
/*       2N int   for _SELSHF  [can use IWORK] */

    wsreqr = ixwork - 1 + *n * 5;
    wsreqi = iywork - 1 + *n + 2;
    if (*liwork == -1 || *lrwork == -1) {
	if (*liwork == -1) {
	    *liwork = wsreqi;
	}
	if (*lrwork == -1) {
	    *lrwork = wsreqr;
	}
	return 0;
    }
    if (*liwork < wsreqi || *lrwork < wsreqr) {
	*info = -1;
	return 0;
    }
/*     ------------------------------------------------------------- */
    eps = dlamch_("Epsilon", (ftnlen)7);
    prec = dlamch_("Precision", (ftnlen)9);
    jxg = 0;
    jyomga = 4;
    pivbase = frepr[(*n << 2) + 3];
    clb = fewl_lu__[(*icbeg << 1) - 1];
    cub = fewl_lu__[*icend * 2];
/* Computing MAX */
    d__1 = abs(clb), d__2 = abs(cub);
    rclwid = (cub - clb) / max(d__1,d__2);
/*     Setup tolerance for how wide the eigenvalues of the child have */
/*     to be so that we accept it (only used if FORCEGAP = .TRUE.). */
/*     [ Note: don't make this too small, otherwise the refinement */
/*       for revealing a gap becomes too expensive ] */
/*     We take the candidate immediately if its ews have relative width */
/*     exceeding WNTRW1, if it only exceeds WNTRW0 we may keep the */
/*     candidate as fallback. */
/* Computing MAX */
/* Computing 2nd power */
    d__3 = *gaptol;
/* Computing MIN */
    d__4 = *gaptol, d__5 = rclwid / *gaptol;
    d__1 = d__3 * d__3, d__2 = min(d__4,d__5);
    wntrw = max(d__1,d__2);

/*     .. Select which random perturbations for the shift to use .. */

/*      Set SHPOFF=0 to deactivate, SHPOFF=1 otherwise */
    shpoff = 0;
/* Computing MAX */
    d__1 = abs(clb), d__2 = abs(cub);
    if (cub - clb <= sqrt(prec) * max(d__1,d__2)) {
	shpoff = 1;
    }

/*     .. Prepare the candidates .. */

    dlaxrf_selshf__(n, &frepr[1], &frepi[1], icbeg, icend, flgap, fugap, &
	    fewl_ae__[1], &fewl_lu__[1], &frginfo[1], taubar, gaptol, &maxcpo,
	     &maxcpi, &maxncand, &maxnb, &ncand, &iwork[iycloc], &rwork[
	    ixctau], &nbatch, abend, &iwork[1]);
/*     Note: can use full IWORK, as not yet in use */
    if (ncand == 0) {
	*info = 1;
	return 0;
    }

/*     =================================================================== */
/*     =                          MAIN LOOP                              = */
/*     =================================================================== */

/*     Candidates go through several phases: */
/*       active if LOC != 0 */
/*       evaluated if TWI != 0 */
/*     Init all candidates as actuve but not evaluated and unfudged. */
    i__1 = iyctwi + maxncand - 1;
    for (ixf77a = iyctwi; ixf77a <= i__1; ++ixf77a) {
	iwork[ixf77a] = 0;
/* L90001: */
    }
/* L90002: */
    i__1 = iycfdg + maxncand - 1;
    for (ixf77a = iycfdg; ixf77a <= i__1; ++ixf77a) {
	iwork[ixf77a] = 0;
/* L90003: */
    }
/* L90004: */
    i__1 = iycgrt + maxncand - 1;
    for (ixf77a = iycgrt; ixf77a <= i__1; ++ixf77a) {
	iwork[ixf77a] = 0;
/* L90005: */
    }
/* L90006: */

/*     The blown up father bounds are zero as long as they are not set. */
    i__1 = ixblup + (*icend << 1);
    for (ixf77a = ixblup + (*icbeg << 1) - 1; ixf77a <= i__1; ++ixf77a) {
	rwork[ixf77a] = 0.;
/* L90007: */
    }
/* L90008: */

/*     BBEG:BEND hold the offsets (0-based) of the candidates in the */
/*     current batch. */
    bbeg = 0;
    bend = abend[0] - 1;
    ibatch = 1;
    lasthope = FALSE_;
    success = FALSE_;
L90009:
    if (ibatch > nbatch + 1) {
	goto L90010;
    }
/*        ---------------------------------------------------------------- */
/*                          Select Next Candidate */
/*        ---------------------------------------------------------------- */
/*        We evaluate the candidates in the order they are given. */
/*        If for the current batch all are evaluated we proceed with the */
/*        best one. */
    bbesteval = -1.;
    nibweval = 0;
    nibactive = 0;
    i__1 = bend;
    for (i__ = bbeg; i__ <= i__1; ++i__) {
	if (iwork[iycloc + i__] != 0) {
/*              candidate is active .. */
	    ++nibactive;
	    if (iwork[iyctwi + i__] != 0) {
/*                 .. and has been evaluated */
		++nibweval;
		if (nibweval == 1) {
		    bbesteval = rwork[ixcevl + i__];
		} else {
/* Computing MIN */
		    d__1 = bbesteval, d__2 = rwork[ixcevl + i__];
		    bbesteval = min(d__1,d__2);
		}
	    }
	}
/* L90011: */
    }
/* L90012: */
    icand = -1;
    i__1 = bend;
    for (i__ = bbeg; i__ <= i__1; ++i__) {
	if (iwork[iycloc + i__] != 0) {
/*              candidate is active .. */
	    if (iwork[iyctwi + i__] == 0) {
/*                 .. but has not yet been evaluated */
		icand = i__;
		goto L90014;
	    }
/*              .. and has already been evaluated .. */
	    if (rwork[ixcevl + i__] <= 5. || lasthope) {
/*                 .. and applies to current phase */
		if (icand == -1 || rwork[ixcevl + i__] < rwork[ixcevl + icand]
			) {
		    icand = i__;
		}
	    }
	}
/* L90013: */
    }
L90014:
    if (icand == -1) {
/*           no candidate found in current batch, proceed to the next */
	++ibatch;
	if (ibatch <= nbatch) {
	    bbeg = abend[ibatch - 2];
	    bend = abend[ibatch - 1] - 1;
	} else if (ibatch == nbatch + 1) {
	    bbeg = 0;
	    bend = ncand - 1;
/* Take fallback if we have one */
	    lasthope = TRUE_;
	} else {
	    bbeg = -1;
	    bend = -1;
	}
	goto L90009;
    }
    twist = 0;
    eval = -1.;
/*        ************************************************************** */
/*        *  Now a series of tests are performed to see if this is an  * */
/*        *  acceptable candidate. To keep the loop structure simple,  * */
/*        *  as soon as a test fails, we modify the candidate record   * */
/*        *  (eg to deactivate or signal as last resort) and cycle     * */
/*        *  back to the beginning of the loop.                        * */
/*        ************************************************************** */
    loc = (i__1 = iwork[iycloc + icand], abs(i__1));
    if (loc > *icend) {
/*           Special case: Shift not coupled to an eigenvalue */
	dir = 0;
	loc = 0;
    } else {
	dir = i_sign(&c__1, &iwork[iycloc + icand]);
    }
    *tau = rwork[ixctau + icand];
    nfudge = iwork[iycfdg + icand];
    shpsel = *n + 1 + dir * *n * shpoff;
    hadeval = iwork[iyctwi + icand] != 0;
    libnoneval = ! hadeval && nibweval == nibactive - 1;
    isgapretry = iwork[iycgrt + icand] > 0;
    dokill = FALSE_;
/*        Special cases: */
/*         LIBNONEVAL indicates that all other active candidates in the */
/*         batch have an evaluation already, only this candidate did not. */
/*         BBESTEVAL is the best evaluation of a candidate that we looked */
/*         at this round; this equals EVAL if HADEVAL is true. */
/*        If this is the first shift to this location, we may allow */
/*        backing off to get the zero inertia consistent. */
/*        Otherwise, if previous inflating of the bound came to close */
/*        to this shift, we kill it. */
    indfbb = 0;
    if (dir != 0) {
	dokill = FALSE_;
	if (dir == 1) {
	    indfbb = ixblup + (loc << 1);
	    if (rwork[indfbb] != 0. && ! hadeval && nfudge == 0 && ! 
		    isgapretry) {
		dokill = *tau <= rwork[indfbb];
	    }
	} else {
	    indfbb = ixblup + (loc << 1) - 1;
	    if (rwork[indfbb] != 0. && ! hadeval && nfudge == 0 && ! 
		    isgapretry) {
		dokill = *tau >= rwork[indfbb];
	    }
	}
	if (dokill) {
/*              too close, deactivate */
	    iwork[iycloc + icand] = 0;
	    goto L90009;
	}
    }

/*        ---------------------------------------------------------------- */
/*                            Perform the Shift */
/*        ---------------------------------------------------------------- */

/*        Retry until the inertia is consistent. */
    nfbtrf = 0;
L90015:
    if (hadeval) {
	kl = iwork[iyctwi + icand];
	ku = kl;
    } else {
	kl = 1;
	ku = *n;
    }
    if (! isgapretry) {
	dlaxrs_(n, &kl, &ku, &frepr[1], &frepi[1], tau, &shfmod[shpsel], &
		rwork[ixdp], &iwork[iyomdp], &rwork[ixrp], &iwork[iyomrp], &
		rwork[ixgmap - 1 + kl], &iwork[iytwok - 1 + kl], &rwork[
		ixwork]);
    } else {
/*              Same call as above, just with FREPR_GRT instead of FREPR */
	dlaxrs_(n, &kl, &ku, &frepr_grt__[1], &frepi[1], tau, &shfmod[shpsel],
		 &rwork[ixdp], &iwork[iyomdp], &rwork[ixrp], &iwork[iyomrp], &
		rwork[ixgmap - 1 + kl], &iwork[iytwok - 1 + kl], &rwork[
		ixwork]);
    }
/*           Select twist */
    if (hadeval) {
	twist = iwork[iyctwi + icand];
	xizero = iwork[iycxiz + icand];
	eval = rwork[ixcevl + icand];
    } else {
/*              Need full data here, KL=1, KU=N */
	d__1 = min(*flgap,*fugap);
	L__1 = dir != 0;
	dlaxrf_seltw__(n, &kminus, &rwork[ixdp], &iwork[iyomdp], &rwork[ixrp],
		 &iwork[iyomrp], &e[1], &rwork[ixgmap], &iwork[iytwok], &env[
		1], spdiam, &d__1, &dir, &loc, tau, &L__1, &twist, &eval, &
		xizero, &rwork[ixwork], &iwork[iywork]);
	if (twist == -1) {
	    iwork[iycloc + icand] = 0;
	    goto L90009;
	}
    }
/*           Determine if the inertia is consistent */
    if (twist == 0) {
	inerok = FALSE_;
    } else {
	inerok = dir == 0 || dir == -1 && xizero <= (loc << 1) - 1 || dir == 
		1 && xizero >= (loc << 1) - 1;
    }
    if (! inerok) {
/*              Initiate a retry if possible, or kill the candidate. */
/*              An fbtrf-retry is only considered if the transformed */
/*              father bound is not yet set. This implies that the */
/*              candidate has no evaluation yet. */
	if (rwork[indfbb] == 0. && nfbtrf <= 4) {
	    ++nfbtrf;
	    *tau += dir * abs(*tau) * pow_ii(&c__4, &nfbtrf) * *n * eps;
/*                  TAU = TAU + DIR*ABS(TAU)*((2**NFBTRF)*FBTRFGRAN*N*EPS) */
	    goto L90015;
	} else {
/*                 No further retry allowed. Kill the candidate. */
	    iwork[iycloc + icand] = 0;
	    goto L90009;
	}
    }

    goto L90016;
    goto L90015;
L90016:
/*        Now that the inertia is ok, update blown up father bounds */
    if (dir != 0) {
	rwork[indfbb] = *tau;
    }
    iwork[iyctwi + icand] = twist;
    iwork[iycxiz + icand] = xizero;
    rwork[ixcevl + icand] = eval;
    rwork[ixctau + icand] = *tau;
/*        Look at evaluation. */
/*        We only proceed with this candidate if */
/*          (a) it is a priori ok   (EVAL <= 1) */
/*        or */
/*          (b) it was the best in the batch (HADEVAL) and either */
/*              is acceptable (EVAL <= EVALOK) or we are desperate */
/*              (LASTHOPE) */
/*        or */
/*          (c) it was the only one without eval in the batch, but */
/*              after evaluating it we can see now that it is the */
/*              best */
/*              [special case to avoid unnessesary recomputation, if */
/*               we put it back now it is selected immediately again] */

    cont_a__ = eval <= 1.;
    cont_b__ = hadeval && (eval <= 5. || lasthope);
    cont_c__ = libnoneval && eval <= 5. && (bbesteval == -1. || eval <= min(
	    5.,bbesteval));
    if (! (cont_a__ || cont_b__ || cont_c__)) {
	goto L90009;
    }

/*        ---------------------------------------------------------------- */
/*                            Init cached quantities */
/*        ---------------------------------------------------------------- */

    i__1 = twist - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	srepr[jxg + i__] = rwork[ixdp - 1 + i__];
	srepi[jyomga + i__] = iwork[iyomdp - 1 + i__];
/* L90017: */
    }
/* L90018: */

    srepr[jxg + twist] = rwork[ixgmap - 1 + twist];
    if (twist == 1) {
	srepi[jyomga + twist] = iwork[iyomrp - 1 + twist];
    } else if (twist == *n) {
	srepi[jyomga + twist] = iwork[iyomdp - 1 + twist];
    } else {
	srepi[jyomga + twist] = 0;
    }

    i__1 = *n;
    for (i__ = twist + 1; i__ <= i__1; ++i__) {
	srepr[jxg + i__] = rwork[ixrp - 2 + i__];
	srepi[jyomga + i__] = iwork[iyomrp - 1 + i__];
/* L90019: */
    }
/* L90020: */
    dlaxrr_(n, &twist, &frepi[1], &e[1], &pivbase, &srepr[1], &srepi[1]);

/*        ---------------------------------------------------------------- */
/*                            Check Outer Bounds */
/*        ---------------------------------------------------------------- */

    bistrycnt = xmrstats_1.xnumfn;
    dlaxrf_cob__(n, flgap, &fewl_lu__[(*icbeg << 1) - 1], &fewl_lu__[*icend * 
	    2], fugap, &srepr[1], &srepi[1], tau, gaptol, icbeg, icend, &
	    xizero, slgap, sugap, &sewl_ae__[1], &sewl_lu__[1], &cobsta);
    if (abs(cobsta) == 1 || abs(cobsta) == 2) {
/*           Eigenvalue-"Undeflow". If possible, apply small fudge */
/*           to the shift and try again. */
	if (dir != 0 && nfudge < 1) {
/*              Reinit as new candidate, only increase its fudge */
/*              counter. */
	    *tau += abs(*tau) * ((dir << 1) * eps);
	    iwork[iyctwi + icand] = 0;
	    iwork[iycfdg + icand] = nfudge + 1;
/*              Normally we do expect that a small fudge won't change */
/*              the evaluation that much. Nevertheless, treating the */
/*              candidate again as new provides full flexibility and */
/*              avoids that the fudged shift becomes worse than */
/*              another candidate but is still preferred. */
	    xmrstats_1.xnbis_wasted__ += xmrstats_1.xnumfn - bistrycnt;
	    goto L90009;
	}
    }
    if (cobsta != 0) {
	xmrstats_1.xnbis_wasted__ += xmrstats_1.xnumfn - bistrycnt;
	iwork[iycloc + icand] = 0;
	goto L90009;
    }

/*        ---------------------------------------------------------------- */
/*                     Init Inner Bounds from Father */
/*        ---------------------------------------------------------------- */
    dlaxrf_iib__(n, &srepr[1], &srepi[1], tau, icbeg, icend, &fewl_ae__[1], &
	    fewl_lu__[1], slgap, sugap, &sewl_ae__[1], &sewl_lu__[1]);
/*         IF( IIBSTA .NE. 0 )THEN */
/*            Previously we discarded candidates where for some subclusters */
/*            of the father no corresponding bound for the son could be */
/*            found. This was switched off (its main use was for checking */
/*            the shift-relation in the coupled case, also it required */
/*            using the cluster internal gaps, and causes problems with */
/*            parallelization). */

/* *           Failed to obtain inner bounds, which means something is */
/* *           not right with the shift relation. */
/* @LOGGING on */
/*            WRITE(FIDRRF,*)' IIB failed, killing candidate',ICAND */
/* @LOGGING ! */
/* @STATS on */
/*            XNBIS_WASTED = XNBIS_WASTED + (XNUMFN - BISTRYCNT) */
/* @STATS ! */
/*            IWORK(IYCLOC + ICAND) = 0 */
/*            CYCLE candloop */
/*         ENDIF */
/*        ---------------------------------------------------------------- */
/*              [Optional]   See if the child has a gap */
/*        ---------------------------------------------------------------- */
    if (TRUE_) {
L90021:
	xa = sewl_lu__[(*icbeg << 1) - 1];
	xb = sewl_lu__[*icbeg * 2];
	xc = sewl_lu__[(*icend << 1) - 1];
	xd = sewl_lu__[*icend * 2];
/*              Note that we may have XA=XC, XB=XD */
/* Computing MAX */
	d__1 = abs(xb), d__2 = abs(xc);
	if (xc - xb >= wntrw * max(d__1,d__2)) {
/*                 Jep, wide enough, take it */
	    goto L90022;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    d__1 = abs(xa), d__2 = abs(xd);
	    if (xd - xa < wntrw * max(d__1,d__2)) {
/*                 No this child's ews are not wide enough and will */
/*                 never be. */
		if (*dogrt && ! isgapretry) {
/*                    Initiate a gap retry */
		    iwork[iycgrt + icand] = 1;
		    iwork[iyctwi + icand] = 0;
		    iwork[iycfdg + icand] = 0;
		} else {
/*                    Discard the candidate */
		    iwork[iycloc + icand] = 0;
		}
		xmrstats_1.xnbis_wasted__ += xmrstats_1.xnumfn - bistrycnt;
		goto L90009;
	    } else {
/*                 We cannot decide yet, bisect the larger interval. */
		if (xb - xa >= xd - xc) {
		    i__ = *icbeg;
		    j = sewl_ae__[i__ * 2];
		    mid = (xa + xb) * .5;
		} else {
		    j = *icend;
		    i__ = sewl_ae__[(j << 1) - 1];
		    mid = (xc + xd) * .5;
		}
		xi = dlaxrn_(n, &srepr[1], &srepi[1], &mid);
		++xmrstats_1.xnbis_class__;
		dlaxrl_refine__(icbeg, icend, slgap, sugap, &sewl_ae__[1], &
			sewl_lu__[1], &i__, &j, &mid, &xi);
	    }
	}
	goto L90021;
L90022:
	;
    }
/*        ---------------------------------------------------------------- */
/*        At this point the candidate passed all tests. */
/*        The son's data (SEWL,...) and TAU are set already. */
    success = TRUE_;
    goto L90010;
    goto L90009;
L90010:
    *info = 0;
    if (! success) {
	*info = 2;
    }
/* @LOGGING on */
/* 99         FORMAT(1X,I5,1X,@(FL)) */
/*            IF( N .LE. MAXLOGDIM )THEN */
/*               WRITE(FIDRRF,*) ' SUCCESS' */
/*               WRITE(FIDRRF,*) ' Bounds after COB:' */
/*               WRITE(FIDRRF,*) '  left gap  = ', SLGAP */
/*               I = ICBEG */
/*               DO */
/*                  J = SEWL_AE(2*I) */
/*                  WRITE(FIDRRF,FMT=99)  I, SEWL_LU(2*I-1) */
/*                  WRITE(FIDRRF,FMT=99)  J, SEWL_LU(2*J) */
/*                  IF( J .EQ. ICEND ) EXIT */
/*                  I = J+1 */
/*                  WRITE(FIDRRF,*) */
/*     $              '  [abs gap =',SEWL_LU(2*I-1)-SEWL_LU(2*J) */
/*               ENDDO */
/*               WRITE(FIDRRF,*) '  right gap = ', SUGAP */
/*            ENDIF */
/* @LOGGING ! */
    return 0;
} /* dlaxrf_ */

