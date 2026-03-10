/* dlarrc.f -- translated by f2c (version 20240504).
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

/* Table of constant values */

static integer c__1 = 1;


/*    Reminder to Osni: */
/*       This is an extension of DLARRF . */
/* 	Needs 4*N workspace compared to DLARRF . */

/* Subroutine */ int dlarrc_hgb_(integer *n, doublereal *d__, doublereal *l,
	doublereal *ld, integer *ifirst, integer *ilast, doublereal *w, 
	doublereal *sigma, doublereal *dplus, doublereal *lplus, integer *
	ndepth, doublereal *ab10, doublereal *dcoup, doublereal *lcoup, 
	doublereal *work, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    integer i__;
    doublereal s, sp1, eps, tmp;
    integer inds;
    doublereal dmax1, dmax2;
    integer inds1, inds2;
    doublereal delta;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1999 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Common block to return timings .. */
/*     .. */
/*     .. Array in Common .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  Given the initial representation L D L^T and its cluster of close */
/*  eigenvalues (in a relative measure), W( IFIRST ), W( IFIRST+1 ), ... */
/*  W( ILAST ), DLARRC finds a new relatively robust representation */
/*  L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the */
/*  eigenvalues of L(+) D(+) L(+)^T is relatively isolated. */

/*  For (NDEPTH.EQ.0) the corresponding coupled representation */
/*  LCOUP(+) DCOUP(+) LCOUP(+)^T is set up. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the matrix. */

/*  D       (input) DOUBLE PRECISION array, dimension (N) */
/*          The N diagonal elements of the diagonal matrix D. */

/*  L       (input) DOUBLE PRECISION array, dimension (N-1) */
/*          The (N-1) subdiagonal elements of the unit bidiagonal */
/*          matrix L. */

/*  LD      (input) DOUBLE PRECISION array, dimension (N-1) */
/*          The (N-1) elements L(i)*D(i). */

/*  LLD     (input) DOUBLE PRECISION array, dimension (N-1) */
/*          The (N-1) elements L(i)*L(i)*D(i). */

/*  IFIRST  (input) INTEGER */
/*          The index of the first eigenvalue in the cluster. */

/*  ILAST   (input) INTEGER */
/*          The index of the last eigenvalue in the cluster. */

/*  W       (input) DOUBLE PRECISION array, dimension (N) */
/*          The eigenvalues of L D L^T in ascending order. */
/*          W( IFIRST ) through W( ILAST ) form the cluster of relatively */
/*          close eigenalues. */

/*  SIGMA   (output) DOUBLE PRECISION */
/*          The shift used to form L(+) D(+) L(+)^T. */

/*  DPLUS   (output) DOUBLE PRECISION array, dimension (N) */
/*          The N diagonal elements of the diagonal matrix D(+). */

/*  LPLUS   (output) DOUBLE PRECISION array, dimension (N-1) */
/*          The first (N-1) elements of LPLUS contain the subdiagonal */
/*          elements of the unit bidiagonal matrix L(+). */

/*  NDEPTH  (input) INTEGER */
/*          Indicates if the initial matrix is positive definite */
/*          (NDEPTH.EQ.0) or not. */

/*  AB10    (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of B B^T . */

/*  DCOUP   (output) DOUBLE PRECISION array, dimension (N) */
/*          Diagonal elements of the diagonal matrix of the */
/* 	   coupled representation LCOUP DCOUP LCOUP^T . */

/*  LCOUP   (output) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of the lower unit bidiagonal matrix */
/* 	   of the coupled representation LCOUP DCOUP LCOUP^T . */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N) */
/*          Workspace. */

/*  INFO    (input) INTEGER */
/*          Error flag. */

/*  Further Details */
/*  =============== */

/*  Based on contributions by */
/*     Inderjit Dhillon, IBM Almaden, USA */
/*     Osni Marques, LBNL/NERSC, USA */
/*     Benedikt Grosser, BUGH Wuppertal, Germany */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --work;
    --lcoup;
    --dcoup;
    --ab10;
    --lplus;
    --dplus;
    --w;
    --ld;
    --l;
    --d__;

    /* Function Body */
    *info = 0;
    inds1 = (*n << 1) + 1;
    inds2 = *n * 3 + 1;
    eps = dlamch_("Precision", (ftnlen)9);
    *sigma = w[*ifirst];
    delta = eps * 2.;

/*     Compute the new relatively robust representation (RRR) */

L10:
    work[inds1] = -(*sigma);
    dplus[1] = d__[1] + work[inds1];
    dmax1 = abs(dplus[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lplus[i__] = ld[i__] / dplus[i__];
	work[inds1 + i__] = work[inds1 + i__ - 1] * lplus[i__] * l[i__] - *
		sigma;
	dplus[i__ + 1] = d__[i__ + 1] + work[inds1 + i__];
/* Computing MAX */
	d__2 = dmax1, d__3 = (d__1 = dplus[i__ + 1], abs(d__1));
	dmax1 = max(d__2,d__3);
/* L20: */
    }
    if (! (dmax1 > 0. || dmax1 < 1.)) {
	*sigma -= abs(*sigma) * delta;
	delta *= 2.;
	goto L10;
    }

    tmp = w[*ilast];
    delta = eps * 2.;
L30:
    work[inds2] = -tmp;
    work[1] = d__[1] + work[inds2];
    dmax2 = abs(work[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[*n + i__] = ld[i__] / work[i__];
	work[inds2 + i__] = work[inds2 + i__ - 1] * work[*n + i__] * l[i__] - 
		tmp;
	work[i__ + 1] = d__[i__ + 1] + work[inds2 + i__];
/* Computing MAX */
	d__2 = dmax2, d__3 = (d__1 = work[i__ + 1], abs(d__1));
	dmax2 = max(d__2,d__3);
/* L40: */
    }
    if (! (dmax2 > 0. || dmax2 < 1.)) {
	tmp += abs(tmp) * delta;
	delta *= 2.;
	goto L30;
    }

    inds = inds1;
    if (dmax2 < dmax1) {
	*sigma = tmp;
	dcopy_(n, &work[1], &c__1, &dplus[1], &c__1);
	i__1 = *n - 1;
	dcopy_(&i__1, &work[*n + 1], &c__1, &lplus[1], &c__1);
	inds = inds2;
    }

/*     Now compute the corresponding representation */
/*     LCOUP(+) DCOUP(+) LCOUP(+)^T */
/*     using coupling transformations. */

    if (*ndepth == 0) {
	sp1 = work[inds];
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s = sp1;
	    sp1 = work[inds + i__];
	    dcoup[i__] = sp1 * dplus[i__] / s;
	    lcoup[i__] = ab10[i__] / dcoup[i__];
/* L50: */
	}
	dcoup[*n] = -(*sigma) * dplus[*n] / sp1;
    } else {

/*       No backward stable couplings, flag an leave. */

	*info = 3;
    }

    return 0;

/*     End of DLARRC */

} /* dlarrc_hgb_ */

