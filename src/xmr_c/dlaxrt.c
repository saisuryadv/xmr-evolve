/* dlaxrt_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxrt_(integer *n, integer *j1, integer *j2, doublereal 
	*repr, integer *repi, doublereal *tau, doublereal *dplus, doublereal *
	rplus, doublereal *gammap, integer *twistok, doublereal *rwork)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    logical dogmastd;
    integer i__, k;
    doublereal tn, tp, gma, eps;
    integer ixg;
    extern /* Subroutine */ int dlaxrt_prog__(integer *, integer *, integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), dlaxrt_stat__(integer *
	    , integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    doublereal bigg, tmpa, tmpb, tmpf;
    integer ixtn, ixtp, twok;
    doublereal tmpg1, tmpg2;
    integer ixngn;
    doublereal pmtau, smtau;
    extern doublereal dlamch_(char *, ftnlen);
    integer ixbdet, iyomga;
    doublereal pivmin;
    integer ixgnsq;
    doublereal pivbase;

/*     IMPLICIT NONE */




/*  Purpose */
/*  ======= */

/*     Factorize twisted blocked to twisted non-blocked, for all twists */
/*     in J1:J2. */

/*  ====================================================================== */

/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRT_STAT( */
/*    $             N, J, I0, I1, G, OMEGA, GNSQ, NGN, BDET, */
/*    $             PIVMIN, TAU, GPLUS, S */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, J, I0, I1 */
/*     INTEGER,          INTENT(IN)  ::  OMEGA(0:N+1) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU */
/*     DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N) */

/*     DOUBLE PRECISION, INTENT(OUT)  ::  GPLUS(I0:I1), S(N) */
/*     END SUBROUTINE DLAXRT_STAT */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRT_PROG( */
/*    $             N, J1, J2, I0, I1, G, OMEGA, GNSQ, NGN, BDET, */
/*    $             PIVMIN, TAU, OFF, GPLUS, P */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, J1, J2, I0, I1 */
/*     INTEGER,          INTENT(IN)  ::  OMEGA(0:N+1) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU, OFF */
/*     DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N) */

/*     DOUBLE PRECISION, INTENT(OUT)  ::  GPLUS(I0:I1), P(N) */
/*     END SUBROUTINE DLAXRT_PROG */
/*     END INTERFACE */

/*     .. Constants .. */


/*     .. Parameters .. */


/*     PIVTAUFAC -- see laxrs  (might make sense to use different values) */


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

/*  ----- Executable Statements ----------------------------------------- */

    /* Parameter adjustments */
    --rwork;
    --repi;
    --repr;
    --dplus;
    /* RPLUS is declared as RPLUS(2:N) in Fortran, so adjust by 2 not 1 */
    rplus -= 2;
    --gammap;
    --twistok;

    /* Function Body */
    if (! xmrstats_1.xstealthmode) {
	++xmrstats_1.xnumft;
    }

/*     .. Decode Representation Data .. */

    ixg = 0;
    ixbdet = *n;
    ixngn = (*n << 1) + 1;
    ixgnsq = *n * 3 + 2;
    iyomga = 4;
    pivbase = repr[(*n << 2) + 3];
    k = repi[2];
    eps = dlamch_("Epsilon", (ftnlen)7);
    bigg = dlamch_("Overflow", (ftnlen)8);
/* Computing MAX */
    d__1 = pivbase, d__2 = eps * 0. * abs(*tau);
    pivmin = max(d__1,d__2);
    ixtn = 1;
    ixtp = *n + 1;
    i__1 = min(k,*j2);
    i__2 = *n - 1;
    dlaxrt_stat__(n, &i__1, &c__1, &i__2, &repr[ixg + 1], &repi[iyomga], &
	    repr[ixgnsq + 1], &repr[ixngn + 1], &repr[ixbdet + 1], &pivmin, 
	    tau, &dplus[1], &rwork[ixtn]);
    if (k < *j2) {
	tn = rwork[ixtn - 1 + k];
	i__1 = *n - 1;
	dlaxrt_prog__(n, &k, j2, &c__1, &i__1, &repr[ixg + 1], &repi[iyomga], 
		&repr[ixgnsq + 1], &repr[ixngn + 1], &repr[ixbdet + 1], &
		pivmin, tau, &tn, &dplus[1], &rwork[ixtn]);
    }
    i__1 = max(k,*j1);
    dlaxrt_stat__(n, &i__1, &c__2, n, &repr[ixg + 1], &repi[iyomga], &repr[
	    ixgnsq + 1], &repr[ixngn + 1], &repr[ixbdet + 1], &pivmin, tau, &
	    rplus[2], &rwork[ixtp]);
    if (*j1 < k) {
	tp = rwork[ixtp - 1 + k];
	dlaxrt_prog__(n, &k, j1, &c__2, n, &repr[ixg + 1], &repi[iyomga], &
		repr[ixgnsq + 1], &repr[ixngn + 1], &repr[ixbdet + 1], &
		pivmin, tau, &tp, &rplus[2], &rwork[ixtp]);
    }

/*     Compute gammas */

    if (*j1 <= k && k <= *j2) {
/*        k=1 or k=n also handled by this */
	tn = rwork[ixtn - 1 + k];
	tp = rwork[ixtp - 1 + k];
	gammap[k] = repr[ixg + k] + (tn + tp) - *tau;
	twistok[k] = 1;
    }
/*     Notes: */
/*     (1)  The following requires safe pivmin atm, otherwise */
/*          use safe mults */
/*     (2)  We access the next GAMMAP(i+1) as P before overwriting */
/*          it, so */
/*          DO NOT CHANGE THE LOOP ORDERS */

/*     Upper Part (analogous to below) */

/* Computing MIN */
    i__2 = *j2, i__3 = k - 1;
    i__1 = min(i__2,i__3);
    for (i__ = *j1; i__ <= i__1; ++i__) {
/*        we set gma just to soothe compilers */
	gma = bigg;
	twok = 0;
	dogmastd = repi[iyomga + i__] == 0 || i__ == 1 || i__ == *n;
	if (! dogmastd) {
	    dogmastd = repr[ixg + i__ - 1] == 0.;
	}
	if (dogmastd) {
/*           normal case (no inner block ending) or */
/*           block with zero d */
	    tn = rwork[ixtn - 1 + i__];
	    tp = rwork[ixtp - 1 + i__];
	    gma = tn + tp - *tau;
	    twok = 1;
	} else {
/*           twist where the source has a block-end */
	    if (i__ + 1 == k) {
		pmtau = rwork[ixtp - 1 + i__ + 1] + (repr[ixg + k] - *tau);
	    } else {
		pmtau = rwork[ixtp - 1 + i__ + 1] - *tau;
	    }
	    smtau = rwork[ixtn - 1 + i__ - 1] - *tau;
	    tmpa = repr[ixbdet + i__ - 1] * pmtau / (repr[ixg + i__ - 1] * 
		    rplus[i__ + 1]);
	    tmpb = repr[ixgnsq + i__ - 1] * smtau / (repr[ixg + i__ - 1] * 
		    dplus[i__ - 1]);
/*            TMPG1 = REPR(IXNGN + I) / RPLUS(I+1) */
/*            TMPG2 = REPR(IXG + I-1) / DPLUS(I-1) */
/* Computing MIN */
	    d__2 = abs(tmpa), d__3 = abs(tmpb);
	    if (abs(tmpa) + abs(tmpb) <= (d__1 = tmpa + tmpb, abs(d__1)) * 
		    10. || min(d__2,d__3) <= abs(*tau) * 5.) {
/*              method 1 */
		gma = tmpa + tmpb - *tau;
		twok = 1;
/* something is wrong with method 2, but turned out to be not absolutely */
/* required anyway */
/*            ELSEIF( ABS(TMPG1)+ABS(TMPG2) .LE. FIVE*(TMPG1-TMPG2) ) */
/*     $      THEN */
/* *              method 2 */
/*               TMPF  = REPR(IXG + I) * (PMTAU / RPLUS(I+1)) */
/*               TMPH  = REPR(IXGNSQ + I-1) / REPR(IXG + I-1) */
/*               GMA = ( TMPF + TMPH*(TMPG1-TMPG2) ) - TAU */
/*               TWOK = 1 */
	    }
	    if (twok == 0) {
		tmpg1 = rwork[ixtp - 1 + i__] - *tau;
/*              now TMPG1 = RP(i) as computed by progressive part */
		tmpg2 = repr[ixgnsq + i__ - 1] / dplus[i__ - 1];
		tmpf = tmpg1 - tmpg2;
		if (abs(tmpg1) + abs(tmpg2) <= abs(tmpf) * 10.) {
		    gma = tmpf;
		    twok = 1;
		}
	    }
	}
	twistok[i__] = twok;
	if (twok == 1) {
	    gammap[i__] = gma;
	}
/* L90001: */
    }
/* L90002: */

/*     Lower Part (analogous to above) */

/* Computing MAX */
    i__1 = *j1, i__2 = k + 1;
    i__3 = *j2;
    for (i__ = max(i__1,i__2); i__ <= i__3; ++i__) {
/*        set gma just to soothe compilers */
	gma = bigg;
	twok = 0;
	dogmastd = repi[iyomga + i__] == 0 || i__ == 1 || i__ == *n;
	if (! dogmastd) {
	    dogmastd = repr[ixg + i__ + 1] == 0.;
	}
	if (dogmastd) {
/*           normal case (no block ending) or */
/*           block with zero r */
	    tn = rwork[ixtn - 1 + i__];
	    tp = rwork[ixtp - 1 + i__];
	    gma = tn + tp - *tau;
	    twok = 1;
	} else {
/*           twist where the source has a block-end */
	    if (i__ - 1 == k) {
		pmtau = rwork[ixtp - 1 + i__ - 1] + (repr[ixg + k] - *tau);
	    } else {
		pmtau = rwork[ixtp - 1 + i__ - 1] - *tau;
	    }
	    smtau = rwork[ixtn - 1 + i__ + 1] - *tau;
	    tmpa = repr[ixbdet + i__ + 1] * pmtau / (repr[ixg + i__ + 1] * 
		    dplus[i__ - 1]);
	    tmpb = repr[ixgnsq + i__ + 1] * smtau / (repr[ixg + i__ + 1] * 
		    rplus[i__ + 1]);
/*            TMPG1 = REPR(IXNGN + I) / DPLUS(I-1) */
/*            TMPG2 = REPR(IXG + I+1) / RPLUS(I+1) */
/* Computing MIN */
	    d__2 = abs(tmpa), d__3 = abs(tmpb);
	    if (abs(tmpa) + abs(tmpb) <= (d__1 = tmpa + tmpb, abs(d__1)) * 
		    10. || min(d__2,d__3) <= abs(*tau) * 5.) {
/*              method 1 */
		gma = tmpa + tmpb - *tau;
		twok = 1;
/* see above */
/*            ELSEIF( ABS(TMPG1)+ABS(TMPG2) .LE. FIVE*(TMPG1-TMPG2) ) */
/*     $      THEN */
/* *              method 2 */
/*               TMPF  = REPR(IXG + I) * (PMTAU / DPLUS(I-1)) */
/*               TMPH  = REPR(IXGNSQ + I+1) / REPR(IXG + I+1) */
/*               GMA = ( TMPF + TMPH*(TMPG1-TMPG2) ) - TAU */
/*               TWOK = 1 */
	    }
	    if (twok == 0) {
		tmpg1 = rwork[ixtp - 1 + i__] - *tau;
/*              now TMPG1 = RP(i) as computed by progressive part */
		tmpg2 = repr[ixgnsq + i__ + 1] / rplus[i__ + 1];
		tmpf = tmpg1 - tmpg2;
		if (abs(tmpg1) + abs(tmpg2) <= abs(tmpf) * 10.) {
		    gma = tmpf;
		    twok = 1;
		}
	    }
	}
	twistok[i__] = twok;
	if (twok == 1) {
	    gammap[i__] = gma;
	}
/* L90003: */
    }
/* L90004: */

    return 0;
} /* dlaxrt_ */

