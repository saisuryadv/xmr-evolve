/* dlar1v.f -- translated by f2c (version 20240504).
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
    doublereal timng[12];
} ttime_;

#define ttime_1 ttime_


/*       Reminder to Osni: */

/* 	Changes by Benedikt: */
/* 	* Argument list */
/* 	* call to DLATS1 */
/* 	* split up former SAWNAN into NANINS and NANINP */
/* 	* set FACTRL */

/*       Argument list changed: */
/* 	* dropped N */
/* 	* introduced FACTRL */
/* 	* replaced the former array WORK of length 4*N */
/* 	  by four arrays LPL, UMN, SPL, PMN of length N */

/* ****************************************************************************** */

/* Subroutine */ int dlar1v_hgb_(integer *b1, integer *bn, doublereal *sigma,
	doublereal *d__, doublereal *l, doublereal *ld, doublereal *lld, 
	doublereal *eval, doublereal *gersch, doublereal *z__, doublereal *
	ztz, doublereal *mingma, integer *r__, integer *isuppz, integer *
	factrl, doublereal *lpl, doublereal *umn, doublereal *spl, doublereal 
	*pmn)
{
    /* System generated locals */
    integer i__1;
    logical L__1;

    /* Local variables */
    integer i__, j;
    doublereal s;
    integer r1, r2;
    doublereal t1, t2, eps, tmp, dplus;
    extern /* Subroutine */ int dlats1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, logical *);
    extern doublereal dlamch_(char *, ftnlen), dsecnd_(void);
    logical naninp, nanins;
    doublereal dminus;


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1999 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */

/*  Purpose */
/*  ======= */

/*  DLAR1V computes the (scaled) r-th column of the inverse of */
/*  the sumbmatrix in rows B1 through BN of the tridiagonal matrix */
/*  L D L^T - sigma I. When sigma is close to an eigenvalue, the */
/*  computed vector is an accurate eigenvector. Usually, r corresponds */
/*  to the index where the eigenvector is largest in magnitude. */
/*  The following steps accomplish this computation : */
/*  (a) Stationary qd transform,  L D L^T - sigma I = L(+) D(+) L(+)^T, */
/*  (b) Progressive qd transform, L D L^T - sigma I = U(-) D(-) U(-)^T, */
/*  (c) Computation of the diagonal elements of the inverse of */
/*      L D L^T - sigma I by combining the above transforms, and choosing */
/*      r as the index where the diagonal of the inverse is (one of the) */
/*      largest in magnitude. */
/*  (d) Computation of the (scaled) r-th column of the inverse using the */
/*      twisted factorization obtained by combining the top part of the */
/*      the stationary and the bottom part of the progressive transform. */

/*  Arguments */
/*  ========= */

/*  B1       (input) INTEGER */
/*           First index of the submatrix of L D L^T. */

/*  BN       (input) INTEGER */
/*           Last index of the submatrix of L D L^T. */

/*  SIGMA    (input) DOUBLE PRECISION */
/*           The shift. In order to compute an accurate eigenvector, */
/*           SIGMA should be a good approximation to an eigenvalue */
/*           of L D L^T. */

/*  D        (input) DOUBLE PRECISION array, dimension (N) */
/*           The n diagonal elements of the diagonal matrix D. */

/*  L        (input) DOUBLE PRECISION array, dimension (N-1) */
/*           The (n-1) subdiagonal elements of the unit bidiagonal matrix */
/*           L, in elements 1 to N-1. */

/*  LD       (input) DOUBLE PRECISION array, dimension (N-1) */
/*           The n-1 elements L(i)*D(i). */

/*  LLD      (input) DOUBLE PRECISION array, dimension (N-1) */
/*           The n-1 elements L(i)*L(i)*D(i). */

/*  EVAL     (input) DOUBLE PRECISION */
/*           The eigenvalue corrsponding to SIGMA of a translate of */
/*           L D L^T. When the translate is zero, EVAL may be identical */
/*           to SIGMA. This variable facilitates calls to DLAR1V */
/*           from DLARRV. */

/*  GERSCH   (input) DOUBLE PRECISION array, dimension (2*N) */
/*           The N Gerschgorin intervals of a translate of L D L^T (see */
/*           comments for EVAL above). These are used to restrict */
/*           the initial search for R, when R is input as 0. */

/*  Z        (input/output) DOUBLE PRECISION array, dimension (N) */
/*           On input, all entries of Z must be set to 0. */
/*           On output, Z contains the (scaled) r-th column of the */
/*           inverse. The scaling is such that Z(R) equals 1. */

/*  ZTZ      (output) DOUBLE PRECISION */
/*           The square of the 2-norm of Z. */

/*  MINGMA   (output) DOUBLE PRECISION */
/*           The reciprocal of the largest (in magnitude) diagonal */
/*           element of the inverse of L D L^T - sigma I. */

/*  R        (input/output) INTEGER */
/*           The twist index for the twisted factorization used to */
/*           compute Z. */
/*           On input, 0 <= R <= N. If R is input as 0, R is set to */
/*           the index where (L D L^T - sigma I)^{-1} is largest */
/*           in magnitude. If 1 <= R <= N, R is unchanged. */
/*           On output, R contains the twist index used to compute Z. */

/*  ISUPPZ   (output) INTEGER array, dimension (2) */
/*           The support of the vector in Z, i.e., the vector Z is */
/*           nonzero only in elements ISUPPZ(1) through ISUPPZ( 2 ). */

/*  FACTRL   (output) INTEGER array, dimension (4) */
/*           FACTRL = FActorization ConTRol, is set with algorithmic */
/*           details from the factorizations. */
/* 	    FACTRL(1) = R1 */
/* 	    FACTRL(2) = R2 */
/* 	    FACTRL(3) = 0 or 1 depending on NANINS */
/* 	    FACTRL(4) = 0 or 1 depending on NANINP */

/*  LPL      (input) DOUBLE PRECISION array, dimension (N-1) */
/*           Off-diagonal elements of the lower unit bidiagonal */
/*           factor computed using the dstqds algorithm. */

/*  UMN      (input) DOUBLE PRECISION array, dimension (N-1) */
/*           Off-diagonal elements of the upper unit bidiagonal */
/*           factor computed using the dqds algorithm. */
/*           The arrays LPL and UMN are needed for a later call to */
/*           DLATS2. */

/*  SPL      (input) DOUBLE PRECISION array, dimension (N) */
/*           Auxiliary variables from the dstqds algorithm. */

/*  PMN      (input) DOUBLE PRECISION array, dimension (N) */
/*           Auxiliary variables from the dqds algorithm. */
/* 	    The arrays LPL, SPL, UMN, PMN are needed for a later call */
/* 	    to DLACSV if dealing with the bidiagonal SVD. */

/*  Further Details */
/*  =============== */

/* 	  LPL and SPL store the data generated by the */
/* 	  stationary qd transform for */
/*         L D L^T - sigma I = L(+) D(+) L(+)^T */
/* 	  LPL stores the lower subdiagonal of L(+). */
/* 	  SPL represents the auxiliary variables needed */
/* 	  in the differential stationary qd. */

/* 	  UMN and PMN store the data generated by the */
/* 	  progressive qd transform for */
/*         L D L^T - sigma I = U(-) D(-) U(-)^T */
/* 	  UMN stores the upper subdiagonal of U(-). */
/* 	  PMN represents the auxiliary variables needed */
/* 	  in the differential progressive qd. */

/* 	  The arrays LPL and UMN are needed for a later call to */
/* 	  DLATS2. */
/* 	  The arrays LPL, SPL, UMN, PMN are needed for a later call */
/* 	  to DLACSV if dealing with the bSVD. */

/* 	* In order to keep track of the computation of the qd transforms */
/* 	  there is an additional integer array of length 4 */
/* 	  FACTRL (FActorization ConTRoL) */

/* 	  FACTRL(1) = R1 */
/* 	  FACTRL(2) = R2 */
/* 	  FACTRL(3) = 0 or 1 depending on NANINS */
/* 	  FACTRL(4) = 0 or 1 depending on NANINP */

/* 	  R1 and R2 are used to safe operations: */
/* 	  Instead of computing all elements of L(+) and U(-), i.e., */
/* 	  L(+) = eye(n) + diag([LPL(1),...,LPL(N-1)],-1) */
/* 	  U(-) = eye(n) + diag([UMN(1),...,UMN(N-1)],1) */
/* 	  the computation is stopped earlier: */
/* 	  L(+) = eye(n) + diag([LPL(1),...,LPL(R2-1),0,...,0],-1) */
/* 	  U(-) = eye(n) + diag([0,...,0,UMN(R1),...,UMN(N-1)],1) */
/* 	  R1 and R2 are determined using the Gerschgorin intervals. */

/* 	  NANINS and NANINP indicate whether the factorization */
/* 	  using the stationary and progressive qd transform broke */
/* 	  down (a diagonal pivot element equals to zero). */
/* 	  In this case a different factorization procedure */
/* 	  is used, which in turn has effects when computing */
/* 	  the couplings in the bSVD algorithm. */



/*  Based on contributions by */
/*     Inderjit Dhillon, IBM Almaden, USA */
/*     Osni Marques, LBNL/NERSC, USA */
/*     Benedikt Grosser, BUGH Wuppertal, Germany */

/*     .. */
/*     .. Common block to return timings .. */
/*     .. */
/*     .. Array in Common .. */
/*     .. */
/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Subroutines .. */

/*     Reminder to Osni: */
/* 	new external subroutine: DLATS1 */
/* 	performs part (d) (see purpose) */

/*     .. */
/*     .. Local Scalars .. */

/*     Reminder to Osni: */
/* 	dropped: SAWNAN, FROM, TO, INDP, INDS, INDUMN */


/*     Reminder to Osni: */
/* 	introduced: NANINS, NANINP */
/* 	for couplings in the bSVD code it is necessary to */
/* 	differentiate if an NaN occured in the stationary qd */
/* 	(NANINS) or the progressive qd (NANINP). */

/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --pmn;
    --spl;
    --umn;
    --lpl;
    --factrl;
    --isuppz;
    --z__;
    --gersch;
    --lld;
    --ld;
    --l;
    --d__;

    /* Function Body */
    eps = dlamch_("Precision", (ftnlen)9);
    if (*r__ == 0) {

/*        Eliminate the top and bottom indices from among the possible */
/*        values of R. */

	r1 = *b1;
	r2 = *bn;
	i__1 = *bn;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    if (*eval >= gersch[(i__ << 1) - 1] && *eval <= gersch[i__ * 2]) {
		r1 = i__;
		goto L20;
	    }
/* L10: */
	}
	goto L40;
L20:
	i__1 = *b1;
	for (i__ = *bn; i__ >= i__1; --i__) {
	    if (*eval >= gersch[(i__ << 1) - 1] && *eval <= gersch[i__ * 2]) {
		r2 = i__;
		goto L40;
	    }
/* L30: */
	}
    } else {
	r1 = *r__;
	r2 = *r__;
    }

L40:

/*     Compute the stationary transform (using the differential form) */
/*     until the index R2 */

    t1 = dsecnd_();
    nanins = FALSE_;
    if (*b1 == 1) {
	spl[1] = 0.;
    } else {
	spl[1] = lld[*b1 - 1];
    }
    s = spl[1] - *sigma;
    i__1 = r2 - 1;
    for (i__ = *b1; i__ <= i__1; ++i__) {
	dplus = d__[i__] + s;
	lpl[i__] = ld[i__] / dplus;
	spl[i__ + 1] = s * lpl[i__] * l[i__];
	s = spl[i__ + 1] - *sigma;
/* L50: */
    }

    if (! (s > 0. || s < 1.)) {

/*        Run a slower version of the above loop if a NaN is detected */

	nanins = TRUE_;
	j = *b1 + 1;
L60:
	if (spl[j + 1] > 0. || spl[j + 1] < 1.) {
	    ++j;
	    goto L60;
	}
	spl[j + 1] = lld[j];
	s = spl[j + 1] - *sigma;
	i__1 = r2 - 1;
	for (i__ = j + 1; i__ <= i__1; ++i__) {
	    dplus = d__[i__] + s;
	    lpl[i__] = ld[i__] / dplus;
	    if (lpl[i__] == 0.) {
		spl[i__ + 1] = lld[i__];
	    } else {
		spl[i__ + 1] = s * lpl[i__] * l[i__];
	    }
	    s = spl[i__ + 1] - *sigma;
/* L70: */
	}
    }
    t2 = dsecnd_();
    ttime_1.timng[3] += t2 - t1;

/*     Compute the progressive transform (using the differential form) */
/*     until the index R1 */

    t1 = dsecnd_();
    naninp = FALSE_;
    pmn[*bn] = d__[*bn] - *sigma;
    i__1 = r1;
    for (i__ = *bn - 1; i__ >= i__1; --i__) {
	dminus = lld[i__] + pmn[i__ + 1];
	tmp = d__[i__] / dminus;
	umn[i__] = l[i__] * tmp;
	pmn[i__] = pmn[i__ + 1] * tmp - *sigma;
/* L80: */
    }
    tmp = pmn[r1];
    if (! (tmp > 0. || tmp < 1.)) {

/*        Run a slower version of the above loop if a NaN is detected */

	naninp = TRUE_;
	j = *bn - 3;
L90:
	if (pmn[j + 1] > 0. || pmn[j + 1] < 1.) {
	    --j;
	    goto L90;
	}
	pmn[j + 1] = d__[j + 1] - *sigma;
	i__1 = r1;
	for (i__ = j; i__ >= i__1; --i__) {
	    dminus = lld[i__] + pmn[i__ + 1];
	    tmp = d__[i__] / dminus;
	    umn[i__] = l[i__] * tmp;
	    if (tmp == 0.) {
		pmn[i__] = d__[i__] - *sigma;
	    } else {
		pmn[i__] = pmn[i__ + 1] * tmp - *sigma;
	    }
/* L100: */
	}
    }
    t2 = dsecnd_();
    ttime_1.timng[4] += t2 - t1;

/*     Find the index (from R1 to R2) of the largest (in magnitude) */
/*     diagonal element of the inverse */

    *mingma = spl[r1] + pmn[r1];
    if (*mingma == 0.) {
	*mingma = eps * spl[r1];
    }
    *r__ = r1;
    i__1 = r2 - 1;
    for (i__ = r1; i__ <= i__1; ++i__) {
	tmp = spl[i__ + 1] + pmn[i__ + 1];
	if (tmp == 0.) {
	    tmp = eps * spl[i__ + 1];
	}
	if (abs(tmp) < abs(*mingma)) {
	    *mingma = tmp;
	    *r__ = i__ + 1;
	}
/* L110: */
    }

/*     Compute the (scaled) r-th column of the inverse */


/*     Reminder to Osni: */
/*       modularized into routine DLATS1 */
/*       DLA + T(wisted matrix) S(olver) 1 */

    t1 = dsecnd_();
    L__1 = nanins || naninp;
    dlats1_(&z__[1], ztz, &lpl[1], &umn[1], &ld[1], b1, r__, bn, &isuppz[1], &
	    L__1);
    t2 = dsecnd_();
    ttime_1.timng[5] += t2 - t1;

/*     Set FACTRL */

    factrl[1] = r1;
    factrl[2] = r2;
    if (nanins) {
	factrl[3] = 1;
    } else {
	factrl[3] = 0;
    }
    if (naninp) {
	factrl[4] = 1;
    } else {
	factrl[4] = 0;
    }

    return 0;

/*     End of DLAR1V */

} /* dlar1v_hgb_ */

