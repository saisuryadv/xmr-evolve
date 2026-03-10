/* dlaxrb_refsng_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxrb_refsng__(integer *n, doublereal *repr, integer *
	repi, integer *depth, integer *il, integer *iu, doublereal *lgap, 
	doublereal *ugap, integer *ewl_ae__, doublereal *ewl_lu__, integer *
	rginfo, doublereal *rwork, integer *iwork)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__, j, k, m;
    doublereal gapaccsng, relaccsng, lb, ub, mid;
    integer ilen;
    doublereal slgap, sugap, width;
    extern /* Subroutine */ int dlaxrl_refine__(integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *);
    integer gileft, jxeind, biscnt, girght;
    extern /* Subroutine */ int dlaxrm_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *);
    integer ntobis, jxiner, ixatol;
    doublereal tolsng;
    integer ixtaus;

/*     IMPLICIT NONE */




/*  Purpose */
/*  ======= */

/*     Refine bounds for singletons some more after basic classification. */

/*     Note: This could equally well be done in DLAXRX before setting up */
/*           the RQI, but by doing them all together here we can exploit */
/*           the vectorized negcounts via DLAXRM. */

/*  Arguments */
/*  ========= */

/*  RGINFO  (input) INTEGER, dimension ( IL : IU-1 ) */
/*          RGINFO(i) holds a status flags to indicate the kind of */
/*          relative gap to the right of ew i, as they have been set */
/*          by DLAXRB_CLASS. */
/*          The possible values are listed in gapinfo.inc.f. */

/*  ====================================================================== */

/*     .. Declarations .. */

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
/*     INTERFACE */
/*     SUBROUTINE  DLAXRM( N, REPR, REPI, M, ATAU, AXI ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, M */
/*     INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), ATAU(M) */

/*     INTEGER,          INTENT(OUT)  ::  AXI(M) */
/*     END SUBROUTINE DLAXRM */
/*     END INTERFACE */

/*     .. Parameters .. */

/*     Singletons intervals are refined until */
/*        width <= Min( GAPACCSNG*MINGAP, RELACCSNG*mid ). */

/*      DOUBLE PRECISION, PARAMETER  ::  GAPACCSNG = 0.0001D0 */
/*      DOUBLE PRECISION, PARAMETER  ::  RELACCSNG = 0.0001D0 */

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

/*     -- Executable Statements ----------------------------------------- */

/*     0.0001  FAC / N */
/*     Note: Far better would be to go with minimal initial refinement into */
/*     RQI and optionally do initial refines (based on dlaxrn) of the RQC is */
/*     too big. But needs some restructuring there, so only worthwhile if time */
/*     for initial refines of sng with this tols is significant. */
    /* Parameter adjustments */
    --repi;
    --repr;
    --iwork;
    --rwork;
    --ewl_ae__;
    --ewl_lu__;
    --rginfo;

    /* Function Body */
    gapaccsng = .01f / *n;
    relaccsng = .01f / *n;
    biscnt = xmrstats_1.xnumfn;
    ilen = *iu - *il + 1;

/*     ----------------------- */
/*      Auxiliary data fields */
/*     ----------------------- */
/*     IWORK: */
/*       EIND(1:M)   -  Indices of yet unconverged intervals */
/*       INER(1:M)   -  inertias of their midpoints */
/*     RWORK: */
/*       TAUS(1:M)   -  midpoints of unconverged intervals */
/*       ATOL(1:M)   -  absolute tolerance for the singleton */

    jxeind = 0;
    jxiner = ilen;
    ixtaus = 0;
    ixatol = ilen;
/*     ----------------- */
/*      Scan singletons */
/*     ----------------- */
    m = 0;
    i__ = *il;
L90001:
    j = ewl_ae__[i__ * 2];
    gileft = 2;
    girght = 2;
    slgap = *lgap;
    sugap = *ugap;
    if (i__ != *il) {
	gileft = rginfo[i__ - 1];
	slgap = ewl_lu__[(i__ << 1) - 1] - ewl_lu__[(i__ << 1) - 2];
    }
    if (j != *iu) {
	girght = rginfo[j];
	sugap = ewl_lu__[(j << 1) + 1] - ewl_lu__[j * 2];
    }
    if (i__ == j && gileft == 2 && girght == 2) {
/*           This is a singleton */
	++m;
	iwork[jxeind + m] = i__;
	rwork[ixatol + m] = gapaccsng * min(slgap,sugap);
    }
    if (j == *iu) {
	goto L90002;
    }
    i__ = j + 1;
    goto L90001;
L90002:

/*     =========== */
/*      Main Loop */
/*     =========== */

L90003:

/*        ----------------------- */
/*         Check for convergence */
/*        ----------------------- */

    k = 1;
L90005:
    if (k > m) {
	goto L90006;
    }
    i__ = iwork[jxeind + k];
    lb = ewl_lu__[(i__ << 1) - 1];
    ub = ewl_lu__[i__ * 2];
    mid = (lb + ub) * .5;
    width = ub - lb;
/* Computing MIN */
    d__1 = rwork[ixatol + k], d__2 = relaccsng * abs(mid);
    tolsng = min(d__1,d__2);
    if (mid == lb || mid == ub || width <= tolsng) {
/*              this interval is converged */
	iwork[jxeind + k] = iwork[jxeind + m];
	rwork[ixatol + k] = rwork[ixatol + m];
	--m;
    } else {
/*              not converged, add midpoint for one bisection step */
	++ntobis;
	rwork[ixtaus + k] = mid;
	++k;
    }
    goto L90005;
L90006:

/*        ----------- */
/*         Loop Exit */
/*        ----------- */

    if (m == 0) {
	goto L90004;
    }

/*        --------------------------------------- */
/*         Bisect all unconverged intervals once */
/*        --------------------------------------- */

    dlaxrm_(n, &repr[1], &repi[1], &m, &rwork[ixtaus + 1], &iwork[jxiner + 1])
	    ;
    i__1 = m;
    for (k = 1; k <= i__1; ++k) {
	i__ = iwork[jxeind + k];
	dlaxrl_refine__(il, iu, lgap, ugap, &ewl_ae__[1], &ewl_lu__[1], &i__, 
		&i__, &rwork[ixtaus + k], &iwork[jxiner + k]);
/* L90007: */
    }
/* L90008: */
    goto L90003;
L90004:
    xmrstats_1.xnbis_sng__ += xmrstats_1.xnumfn - biscnt;
    return 0;
} /* dlaxrb_refsng__ */

