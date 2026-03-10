/* dlaxrg_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxrg0_(integer *n, doublereal *e, doublereal *repr, 
	integer *repi, doublereal *cuttol, doublereal *z__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, k;
    doublereal l, u, az;
    integer jxg;
    doublereal oldl, oldu;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    integer ixf77a, iactz, nzfix;
    doublereal prevz;
    extern doublereal dlamch_(char *, ftnlen);
    integer jxbdet, jyomga;
    doublereal invnrm, normsq, uthresh;

/*     IMPLICIT NONE */



/*  Purpose */
/*  ======= */

/*     Computes the FP-vector for a zero eigenvalue, normed to unity. */

/*     Concerning the block structure, there should be no block ending */
/*     at k, regardless of k=1, k=n or in-between. This could be fixed, */
/*     but atm we only compute vectors with blocks for a singular matrix, */
/*     and those should have a single pivot being zero. */

/*  ====================================================================== */

/*     .. Declarations .. */


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

    /* Parameter adjustments */
    --z__;
    --repi;
    --repr;
    --e;

    /* Function Body */
    if (! xmrstats_1.xstealthmode) {
	++xmrstats_1.xnumgv0;
    }
    jxg = 0;
    jxbdet = *n;
    jyomga = 4;
    k = repi[2];
    uthresh = dlamch_("Underflow", (ftnlen)9);
    z__[k] = 1.;
    normsq = 1.;
    nzfix = 0;

/*     Solve from twist on down to 1 */

    prevz = 1.;
    oldl = 0.;
    i__ = k - 1;
L90001:
    if (i__ == 0) {
	goto L90002;
    }
/* assume we have PREVZ such that Z(i+1) was computed as -OLDL*PREVZ */
/* if i+2 ends a block then PREVZ=Z(i+3), otherwise PREVZ=Z(i+2) */
/* CC */
    if (repi[jyomga + i__] != 0) {
/*           I ends a block, implies i < k here --> L(i) = E(i)/invD(i) */
	l = e[i__] * repr[jxg + i__ - 1] / repr[jxbdet + i__ - 1];
	iactz = i__ + 1;
    } else if (repi[jyomga + i__ + 1] != 0) {
/*           I starts a block, implies i < k-1 here */
	l = -(e[i__ + 1] * e[i__]) / repr[jxbdet + i__];
	iactz = i__ + 2;
    } else {
	l = e[i__] / repr[jxg + i__];
	iactz = i__ + 1;
    }

/*        Compute next vector entry, handle possible creeping underflow */

    if ((d__1 = z__[iactz], abs(d__1)) <= uthresh) {
/*           We should have I < K-1. */
/*           We also have |OLDL| < 1, |PREVZ| < 1, thus if there is */
/*           underflow anywhere in the following product, the exact */
/*           result would underflow as well. */
	z__[iactz] = 0.;
	z__[i__] = l * oldl * prevz;
	++nzfix;
    } else {
	z__[i__] = -l * z__[iactz];
    }

/*        support cut */

    az = (d__1 = z__[i__], abs(d__1)) + (d__2 = z__[i__ + 1], abs(d__2));
    if (az <= uthresh || az * e[i__] <= *cuttol) {
	i__1 = i__;
	for (ixf77a = 1; ixf77a <= i__1; ++ixf77a) {
	    z__[ixf77a] = 0.;
/* L90003: */
	}
/* L90004: */
	goto L90002;
    } else {
/* Computing 2nd power */
	d__1 = z__[i__];
	normsq += d__1 * d__1;
    }
    prevz = z__[iactz];
    oldl = l;
    --i__;
    goto L90001;
L90002:
    if (abs(z__[1]) <= uthresh) {
	z__[1] = 0.;
    }

/*     Solve from twist on down to n */
/*     (completely analogous to upper part) */

    prevz = 1.;
    oldu = 0.;
    i__ = k + 1;
L90005:
    if (i__ == *n + 1) {
	goto L90006;
    }
    if (repi[jyomga + i__] != 0) {
/*           I ends a block, implies i > k here --> U(i) = E(i-1)/invR(i) */
	u = e[i__ - 1] * repr[jxg + i__ + 1] / repr[jxbdet + i__ + 1];
	iactz = i__ - 1;
    } else if (repi[jyomga + i__ - 1] != 0) {
/*           I starts a block, implies i > k+1 here */
	u = -(e[i__ - 2] * e[i__ - 1]) / repr[jxbdet + i__];
	iactz = i__ - 2;
    } else {
	u = e[i__ - 1] / repr[jxg + i__];
	iactz = i__ - 1;
    }
    if ((d__1 = z__[iactz], abs(d__1)) <= uthresh) {
	z__[iactz] = 0.;
	z__[i__] = u * oldu * prevz;
	++nzfix;
    } else {
	z__[i__] = -u * z__[iactz];
    }
    az = (d__1 = z__[i__], abs(d__1)) + (d__2 = z__[i__ - 1], abs(d__2));
    if (az <= uthresh || az * e[i__ - 1] <= *cuttol) {
	i__1 = *n;
	for (ixf77a = i__; ixf77a <= i__1; ++ixf77a) {
	    z__[ixf77a] = 0.;
/* L90007: */
	}
/* L90008: */
	goto L90006;
    } else {
/* Computing 2nd power */
	d__1 = z__[i__];
	normsq += d__1 * d__1;
    }
    prevz = z__[iactz];
    oldu = u;
    ++i__;
    goto L90005;
L90006:
    if ((d__1 = z__[*n], abs(d__1)) <= uthresh) {
	z__[*n] = 0.;
    }

/*     Determine results */

    invnrm = 1. / sqrt(normsq);

/*     Normalize Z */

    dscal_(n, &invnrm, &z__[1], &c__1);
    return 0;
} /* dlaxrg0_ */


/* *********************************************************************** */
/* Subroutine */ int dlaxrg_(integer *n, integer *k, doublereal *d__, 
	doublereal *r__, doublereal *e, doublereal *gamma, doublereal *cuttol,
	 doublereal *z__, doublereal *normz, doublereal *resid, doublereal *
	rqcorr)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__;
    doublereal l, u, az;
    integer isa, ise;
    doublereal oldl, oldu;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    integer ixf77a, nzfix;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal invnrm, normsq, uthresh;

/*     IMPLICIT NONE */



/*  Purpose */
/*  ======= */

/*     Computes the FP-vector, normed to unity. */

/*     We require that the factorization exists, in particular the */
/*     elements L(i) and U(i) should be finite and nonzero. */

/*  ====================================================================== */

/*     .. Declarations .. */


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

    /* Parameter adjustments */
    --z__;
    --d__;
    /* R is declared as R(K+1:N) in Fortran, meaning R(K+1) should be the
       first element of the passed array. f2c's blind --r__ assumes 1-based
       indexing, but we need r__[K+1] = element[0], so adjust by K+1. */
    r__ -= (*k + 1);
    --e;

    /* Function Body */
    if (! xmrstats_1.xstealthmode) {
	++xmrstats_1.xnumgv;
    }
    uthresh = dlamch_("Underflow", (ftnlen)9);
    z__[*k] = 1.;
    normsq = 1.;
    nzfix = 0;
    isa = 1;
    ise = *n;

/*     Solve from twist on down to 1 */

    oldl = 0.;
    i__ = *k - 1;
L90009:
    if (i__ == 0) {
	goto L90010;
    }
    l = e[i__] / d__[i__];

/*        Compute next vector entry, handle possible creeping underflow */

    if ((d__1 = z__[i__ + 1], abs(d__1)) <= uthresh) {
/*           We should have I < K-1. */
/*           We also have |OLDL| < 1, |Z(i+2)| < 1, thus if there is */
/*           underflow anywhere in the following product, the exact */
/*           result would underflow as well. */
	z__[i__ + 1] = 0.;
	z__[i__] = l * oldl * z__[i__ + 2];
	++nzfix;
    } else {
	z__[i__] = -l * z__[i__ + 1];
    }

/*        support cut */

    az = (d__1 = z__[i__], abs(d__1)) + (d__2 = z__[i__ + 1], abs(d__2));
    if (az <= uthresh || az * e[i__] <= *cuttol) {
	i__1 = i__;
	for (ixf77a = 1; ixf77a <= i__1; ++ixf77a) {
	    z__[ixf77a] = 0.;
/* L90011: */
	}
/* L90012: */
	isa = i__ + 1;
	goto L90010;
    } else {
/* Computing 2nd power */
	d__1 = z__[i__];
	normsq += d__1 * d__1;
    }
    oldl = l;
    --i__;
    goto L90009;
L90010:

/*     Solve from twist on down to n */
/*     (completely analogous to upper part) */

    oldu = 0.;
    i__ = *k + 1;
L90013:
    if (i__ == *n + 1) {
	goto L90014;
    }
    u = e[i__ - 1] / r__[i__];
    if ((d__1 = z__[i__ - 1], abs(d__1)) <= uthresh) {
	z__[i__ - 1] = 0.;
	z__[i__] = u * oldu * z__[i__ - 2];
	++nzfix;
    } else {
	z__[i__] = -u * z__[i__ - 1];
    }
    az = (d__1 = z__[i__], abs(d__1)) + (d__2 = z__[i__ - 1], abs(d__2));
    if (az <= uthresh || az * e[i__ - 1] <= *cuttol) {
	i__1 = *n;
	for (ixf77a = i__; ixf77a <= i__1; ++ixf77a) {
	    z__[ixf77a] = 0.;
/* L90015: */
	}
/* L90016: */
	ise = i__ - 1;
	goto L90014;
    } else {
/* Computing 2nd power */
	d__1 = z__[i__];
	normsq += d__1 * d__1;
    }
    oldu = u;
    ++i__;
    goto L90013;
L90014:

/*     Determine results */

    *normz = sqrt(normsq);
    invnrm = 1. / *normz;
    *resid = abs(*gamma) * invnrm;
    *rqcorr = *gamma / normsq;

/*     Normalize Z */

    dscal_(n, &invnrm, &z__[1], &c__1);
    return 0;
} /* dlaxrg_ */

