/* dlaxrn_f77.f -- translated by f2c (version 20240504).
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
static integer c_n1 = -1;

integer dlaxrn_(integer *n, doublereal *repr, integer *repi, doublereal *tau)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    integer i__, k, nb;
    doublereal sn, sp;
    integer ibb;
    doublereal gma;
    integer neg, ixg;
    extern /* Subroutine */ int dlaxrn_stat__(integer *, integer *, integer *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *)
	    ;
    integer negn, negp, ixngn, iylbbk, ixbdet, iyomga;
    doublereal pivmin;
    integer ixgnsq;
    doublereal pivbase;

/*     IMPLICIT NONE */

/*  Purpose */
/*  ======= */

/*    Perform a sturm count by doing a (suitably twisted) bidiagonal */
/*    factorization of T - TAU and then count the negative pivots. */
/*    However, for utmost information the function does not return the */
/*    negcount, but */
/*     xi := 2*negc + issing */
/*    where issing=1 if T-TAU is singular, and 0 otherwise. Thus it pro- */
/*    vides full information about the inertia. */

/*    Consequences / Usage */
/*    -------------------- */
/*         xi = 2k-1 (odd)  means we are on ew k */
/*         xi = 2k   (even) means we are in interval (k,k+1) */

/*    TAU is */
/*      upper bound for ew k         <=>  xi >= 2*k-1 */
/*      strict upper bound for ew k  <=>  xi >  2*k-1 */

/*      lower bound for ew k         <=>  xi <= 2*k-1 */
/*      strict lower bound for ew k  <=>  xi <  2*k-1 */

/*      in (ew k, ew k+1)  <=>  xi = 2k */
/*      in [ew k, ew k+1)  <=>  xi = 2k-1 or xi = 2k */
/*      in (ew k, ew k+1]  <=>  xi = 2k   or xi = 2k+1 */
/*      in [ew k, ew k+1]  <=>  2k-1 <= xi <= 2k+1 */

/*    For given inertia xi and sample lambda, */
/*      the nearest eigenvalue smaller than (not equal to) lambda is */
/*        xi / 2, */
/*      the nearest eigenvalue larger than (not equal to) lambda is */
/*        xi / 2 + MOD(xi,2) + 1  =  (xi+1) / 2 + 1 */

/*  ====================================================================== */


/*     .. Parameters .. */


/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRN_STAT( */
/*    $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET, */
/*    $             PIVMIN, TAU, NCOUNT, S */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, K, DIR */
/*     INTEGER,          INTENT(IN)  ::  LBBEGK(*) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU */
/*     DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N) */

/*     INTEGER,          INTENT(INOUT)  ::  IBB */

/*     INTEGER,          INTENT(OUT)  ::  NCOUNT */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  S */
/*     END SUBROUTINE DLAXRN_STAT */
/*     END INTERFACE */

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

/*  ----- Executable Statements ------------------------------------------ */

    /* Parameter adjustments */
    --repi;
    --repr;

    /* Function Body */
    if (! xmrstats_1.xstealthmode) {
	++xmrstats_1.xnumfn;
    }

/*     .. Decode Representation Data .. */

    ixg = 0;
    ixbdet = *n;
    ixngn = (*n << 1) + 1;
    ixgnsq = *n * 3 + 2;
    iylbbk = *n + 6;
    iyomga = 4;
    pivbase = repr[(*n << 2) + 3];
    k = repi[2];
    nb = repi[3];

/*     .. Special treatment for zero inertia .. */

    if (*tau == 0.) {
/*        NOTE: The source may have blocks, so here we need to */
/*              count them correctly (!) */
	neg = 0;

	ibb = 0;
	i__ = 1;
L90001:
L90003:
	if (i__ >= repi[iylbbk + ibb]) {
	    goto L90004;
	}
	if (repr[ixg + i__] < 0.) {
	    ++neg;
	}
	++i__;
	goto L90003;
L90004:
	if (i__ >= k) {
	    goto L90002;
	}
	++ibb;
	++neg;
	i__ += 2;
	goto L90001;
L90002:

	ibb = nb;
	i__ = *n;
L90005:
L90007:
	if (i__ <= repi[iylbbk + ibb]) {
	    goto L90008;
	}
	if (repr[ixg + i__] < 0.) {
	    ++neg;
	}
	--i__;
	goto L90007;
L90008:
	if (i__ <= k) {
	    goto L90006;
	}
	--ibb;
	++neg;
	i__ += -2;
	goto L90005;
L90006:

	ret_val = neg << 1;
	if (repi[iyomga + k] == 0) {
	    if (repr[ixg + k] < 0.) {
		ret_val += 2;
	    }
	    if (repr[ixg + k] == 0.) {
		++ret_val;
	    }
	}
	return ret_val;

    }

/*     General case, tau is nonzero */

    pivmin = pivbase;
    ibb = 1;
    negn = 0;
    sn = 0.;
    if (k > 1) {
	dlaxrn_stat__(n, &k, &c__1, &repr[ixg + 1], &ibb, &repi[iylbbk], &
		repr[ixgnsq + 1], &repr[ixngn + 1], &repr[ixbdet + 1], &
		pivmin, tau, &negn, &sn);
    }

    ibb = nb + 1;
    negp = 0;
    sp = 0.;
    if (k < *n) {
	dlaxrn_stat__(n, &k, &c_n1, &repr[ixg + 1], &ibb, &repi[iylbbk], &
		repr[ixgnsq + 1], &repr[ixngn + 1], &repr[ixbdet + 1], &
		pivmin, tau, &negp, &sp);
    }

    gma = repr[ixg + k] + (sn + sp - *tau);

    neg = negn + negp;
    if (gma < 0.) {
	++neg;
    }
    ret_val = neg << 1;
    if (gma == 0.) {
	++ret_val;
    }

    return ret_val;
} /* dlaxrn_ */

