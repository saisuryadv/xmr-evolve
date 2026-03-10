/* dlarri.f -- translated by f2c (version 20240504).
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
static doublereal c_b10 = 0.;
static integer c__0 = 0;

/* Subroutine */ int dlarri_(integer *n, logical *wantl, logical *wantr, 
	logical *wanta, doublereal *d__, doublereal *e, doublereal *tol, 
	integer *nsplit, integer *isplit, integer *m, doublereal *w, 
	doublereal *gersch, doublereal *d2, doublereal *e2, doublereal *de00, 
	doublereal *de10, doublereal *work, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    integer i__, j;
    doublereal t1, t2, gl;
    integer in;
    doublereal gu, tmp, offd;
    integer iend, jblk;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlasq2_(integer *, doublereal *, 
	    integer *);
    integer ibegin;
    extern doublereal dsecnd_(void);
    doublereal sgndef;


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

/*  Given the bidiagonal matrix B, DLARRI checks for "small" elements. */
/*  If no splittings occur either B^T B = L D L^T (WANTL or WANTA) */
/*  or B B^T = L D L^T (WANTR) is formed. Note that the decomposition */
/*  is always positive definite. The eigenvalues are found by the dqds */
/*  algorithm (subroutine DLASQ2). As an added benefit, DLARRI also outputs */
/*  the n Gerschgorin intervals for the matrices L D L^T. */

/*  Reminder to Osni: */
/*  Note 4 : Splittings are currently not supported, flagged by INFO = 33! */


/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the matrix. */

/*  WANTL   (input) LOGICAL */
/*          Only left singular vectors are desired. */

/*  WANTR   (input) LOGICAL */
/*          Only right singular vectors are desired. */

/*  WANTA   (input) LOGICAL */
/*          Both left and right singular vectors are desired. */

/*  D       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the N diagonal elements of the bidiagonal */
/*          matrix B. */
/*          On exit, the N diagonal elements of the diagonal */
/*          matrix D. */

/*  E       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the first (N-1) entries contain the off-diagonal */
/*          elements of the bidiagonal matrix B; E(N) need not be set. */
/*          On exit, E contains the subdiagonal elements of the unit */
/*          bidiagonal matrices L. */

/*  TOL     (input) DOUBLE PRECISION */
/*          The threshold for splitting. If on input |D(i)| < TOL or */
/*          |E(i)| < TOL, then the matrix B is split into smaller blocks. */

/*  NSPLIT  (output) INTEGER */
/*          The number of blocks B splits into. 1 <= NSPLIT <= N. */

/*  ISPLIT  (output) INTEGER array, dimension (2*N) */
/*          The splitting points, at which B breaks up into submatrices. */
/*          The first submatrix consists of rows/columns 1 to ISPLIT(1), */
/*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2), */
/*          etc., and the NSPLIT-th consists of rows/columns */
/*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N. */

/*  M       (output) INTEGER */
/*          The total number of eigenvalues of L D L^T found. */

/*  W       (output) DOUBLE PRECISION array, dimension (N) */
/*          The first M elements contain the eigenvalues */
/*          sorted in ascending order. */

/*  GERSCH  (output) DOUBLE PRECISION array, dimension (2*N) */
/*          The N Gerschgorin intervals. */

/*  D2      (output) DOUBLE PRECISION array, dimension (N) */
/*          D2(I) = D(I)*D(I) */

/*  E2      (output) DOUBLE PRECISION array, dimension (N-1) */
/*          E2(I) = E(I)*E(I) */

/*  DE00    (output) DOUBLE PRECISION array, dimension (N-1) */
/*          DE00(I) = D(I)*E(I) */

/*  DE01    (output) DOUBLE PRECISION array, dimension (N-1) */
/*          DE01(I) = D(I+1)*E(I) */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension MAX(1,4*N) */
/*          Workspace. */

/*  IWORK   (workspace) INTEGER array, dimension (2*N) */
/*          Workspace. */

/*  INFO    (output) INTEGER */
/*          Output error code from DLASQ2 */

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
    --de10;
    --de00;
    --e2;
    --d2;
    --gersch;
    --w;
    --isplit;
    --e;
    --d__;

    /* Function Body */
    *info = 0;

/*     Compute Splitting Points */

    *nsplit = 1;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = e[i__], abs(d__1)) <= *tol) {
	    isplit[*nsplit] = i__;
	    ++(*nsplit);
	}
/* L10: */
    }
    isplit[*nsplit] = *n;

/*     Check if matrix splits. */


/*     Reminder to Osni: */
/*     Cannot cope with splittings and negative entries. */

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (d__[i__] <= 0. || e[i__] <= 0.) {
	    *info = 33;
	    return 0;
	}
	if ((d__1 = d__[i__], abs(d__1)) <= *tol || (d__2 = e[i__], abs(d__2))
		 <= *tol) {
	    *info = 33;
	    return 0;
	}
/* L20: */
    }
    if (d__[*n] <= 0.) {
	*info = 33;
	return 0;
    }

    ibegin = 1;
    i__1 = *nsplit;
    for (jblk = 1; jblk <= i__1; ++jblk) {
	iend = isplit[jblk];
	if (ibegin == iend) {
	    w[ibegin] = d__[ibegin];
	    e[ibegin] = 0.;
	    ibegin = iend + 1;
	    goto L170;
	}
	in = iend - ibegin + 1;

/*        Reminder to Osni: */
/* 	 Assumption: IBEGIN = 1, IEND = N, IN = N */


/* 	 Compute some auxiliary variables. */

	i__2 = iend - 1;
	for (i__ = ibegin; i__ <= i__2; ++i__) {
	    d2[i__] = d__[i__] * d__[i__];
	    e2[i__] = e[i__] * e[i__];
	    de00[i__] = d__[i__] * e[i__];
	    de10[i__] = d__[i__ + 1] * e[i__];
/* L30: */
	}
	d2[iend] = d__[iend] * d__[iend];

/*        Form the IN Gerschgorin intervals */
/* 	 If (WANTR) then form T = B^T B explicitly */
/* 	   and find its Gerschgorin intervals. */
/* 	 If (WANTL) then form T = B B^T explicitly */
/* 	   and find its Gerschgorin intervals. */
/* 	 If (WANTA) then unify Gerschgorin intervals */
/* 	   of B^T B and B B^T. */

/* 	 Store input data: WORK(1:N) = D(1:N) */
/* 	                   WORK(N+1:2*N) = E(1:N) */

	dcopy_(&in, &d__[ibegin], &c__1, &work[1], &c__1);
	i__2 = in - 1;
	dcopy_(&i__2, &e[ibegin], &c__1, &work[*n + 1], &c__1);
	i__2 = in << 1;
	dcopy_(&i__2, &c_b10, &c__0, &gersch[(ibegin << 1) - 1], &c__1);

	if (*wantr || *wanta) {

/* 	   Compute T = B^T B and find its Gerschgorin intervals. */

	    d__[ibegin] = d2[ibegin];
	    i__2 = iend - 1;
	    for (i__ = ibegin; i__ <= i__2; ++i__) {
		d__[i__ + 1] = d2[i__ + 1] + e2[i__];
		e[i__] = de00[i__];
/* L40: */
	    }
	    gl = d__[ibegin] - (d__1 = e[ibegin], abs(d__1));
	    gu = d__[ibegin] + (d__1 = e[ibegin], abs(d__1));
	    gersch[(ibegin << 1) - 1] = gl;
	    gersch[ibegin * 2] = gu;
	    gersch[(iend << 1) - 1] = d__[iend] - (d__1 = e[iend - 1], abs(
		    d__1));
	    gersch[iend * 2] = d__[iend] + (d__1 = e[iend - 1], abs(d__1));
/* Computing MIN */
	    d__1 = gersch[(iend << 1) - 1];
	    gl = min(d__1,gl);
/* Computing MAX */
	    d__1 = gersch[iend * 2];
	    gu = max(d__1,gu);
	    i__2 = iend - 1;
	    for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
		offd = (d__1 = e[i__ - 1], abs(d__1)) + (d__2 = e[i__], abs(
			d__2));
		gersch[(i__ << 1) - 1] = d__[i__] - offd;
/* Computing MIN */
		d__1 = gersch[(i__ << 1) - 1];
		gl = min(d__1,gl);
		gersch[i__ * 2] = d__[i__] + offd;
/* Computing MAX */
		d__1 = gersch[i__ * 2];
		gu = max(d__1,gu);
/* L50: */
	    }
	}
	if (*wantl || *wanta) {

/* 	   Compute T = B B^T and find its Gerschgorin intervals. */

	    i__2 = iend - 1;
	    for (i__ = ibegin; i__ <= i__2; ++i__) {
		d__[i__] = d2[i__] + e2[i__];
		e[i__] = de10[i__];
/* L60: */
	    }
	    d__[iend] = d2[iend];
/* Computing MIN */
	    d__2 = gl, d__3 = d__[ibegin] - (d__1 = e[ibegin], abs(d__1));
	    gl = min(d__2,d__3);
/* Computing MAX */
	    d__2 = gu, d__3 = d__[ibegin] + (d__1 = e[ibegin], abs(d__1));
	    gu = max(d__2,d__3);
/* Computing MIN */
	    d__1 = gl, d__2 = gersch[(ibegin << 1) - 1];
	    gersch[(ibegin << 1) - 1] = min(d__1,d__2);
/* Computing MAX */
	    d__1 = gu, d__2 = gersch[ibegin * 2];
	    gersch[ibegin * 2] = max(d__1,d__2);
/* Computing MIN */
	    d__2 = gersch[(iend << 1) - 1], d__3 = d__[iend] - (d__1 = e[iend 
		    - 1], abs(d__1));
	    gersch[(iend << 1) - 1] = min(d__2,d__3);
/* Computing MAX */
	    d__2 = gersch[iend * 2], d__3 = d__[iend] + (d__1 = e[iend - 1], 
		    abs(d__1));
	    gersch[iend * 2] = max(d__2,d__3);
/* Computing MIN */
	    d__1 = gersch[(iend << 1) - 1];
	    gl = min(d__1,gl);
/* Computing MAX */
	    d__1 = gersch[iend * 2];
	    gu = max(d__1,gu);
	    i__2 = iend - 1;
	    for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
		offd = (d__1 = e[i__ - 1], abs(d__1)) + (d__2 = e[i__], abs(
			d__2));
/* Computing MIN */
		d__1 = gersch[(i__ << 1) - 1], d__2 = d__[i__] - offd;
		gersch[(i__ << 1) - 1] = min(d__1,d__2);
/* Computing MIN */
		d__1 = gersch[(i__ << 1) - 1];
		gl = min(d__1,gl);
/* Computing MAX */
		d__1 = gersch[i__ * 2], d__2 = d__[i__] + offd;
		gersch[i__ * 2] = max(d__1,d__2);
/* Computing MAX */
		d__1 = gersch[i__ * 2];
		gu = max(d__1,gu);
/* L70: */
	    }
	}

/* 	 Recover input data: */

	dcopy_(&in, &work[1], &c__1, &d__[ibegin], &c__1);
	i__2 = in - 1;
	dcopy_(&i__2, &work[*n + 1], &c__1, &e[ibegin], &c__1);

/* 	 Compute  B^T B = L D L^T or  B B^T = L D L^T */
/* 	 Since B^T B and B B^T are already definite, there is no */
/* 	 need to find an initial shift parameter as in DLARRE. */

	if (*wantl) {

/*          Compute B B^T = L D L^T . */
/* 	   Input: */
/* 		D(I), I=1,N - diagonal entries of B */
/* 	        E(1:N-1), I=1,N-1 - offdiagonal entries of B */
/* 	   Output: */
/* 		WORK(    I) = d_i , I=1,N - diagonal entries of D */
/* 		WORK(2*N+I) = 1/d_i , I=1,N-1, WORK(3*N) = ONE */
/* 		WORK(  N+I) = l_i , I=1,N-1 - offdiagonal entries of L */

	    tmp = d2[ibegin];
	    i__2 = iend - 1;
	    for (i__ = ibegin; i__ <= i__2; ++i__) {
		work[i__] = tmp + e2[i__];
		work[(*n << 1) + i__] = 1. / work[i__];
		work[*n + i__] = de10[i__] * work[(*n << 1) + i__];
		tmp = d2[i__ + 1] * work[(*n << 1) + i__] * tmp;
/* L80: */
	    }
	    work[iend] = tmp;
	    work[*n * 3] = 1.;

	} else {

/*          Compute B^T B = L D L^T . */
/* 	   Input: */
/* 		D(I), I=1,N - diagonal entries of B */
/* 	        E(1:N-1), I=1,N-1 - offdiagonal entries of B */
/* 	   Output: */
/* 		WORK(    I) = d_i , I=1,N - diagonal entries of D */
/* 		WORK(2*N+I) = 1/d_i , I=1,N-1, WORK(3*N) = ONE */
/* 		WORK(  N+I) = l_i , I=1,N-1 - offdiagonal entries of L */

	    i__2 = iend - 1;
	    for (i__ = ibegin; i__ <= i__2; ++i__) {
		work[i__] = d2[i__];
		work[(*n << 1) + i__] = 1. / work[i__];
		work[*n + i__] = e[i__] / d__[i__];
/* L90: */
	    }
	    work[in] = d2[in];
	    work[*n * 3] = 1.;
	}
	dcopy_(&in, &work[ibegin], &c__1, &d__[ibegin], &c__1);
	i__2 = in - 1;
	dcopy_(&i__2, &work[*n + ibegin], &c__1, &e[ibegin], &c__1);

/* 	 Initialize some relicts of DLARRE : */

	sgndef = 1.;
	e[iend] = 0.;

/*        Compute the eigenvalues of L D L^T. */

	j = ibegin;
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[(i__ << 1) - 1] = (d__1 = d__[j], abs(d__1));
	    work[i__ * 2] = e[j] * e[j] * work[(i__ << 1) - 1];
	    ++j;
/* L140: */
	}
	work[(in << 1) - 1] = (d__1 = d__[iend], abs(d__1));

	t1 = dsecnd_();
	dlasq2_(&in, &work[1], info);
	t2 = dsecnd_();
	ttime_1.timng[0] += t2 - t1;

	j = ibegin;
	if (sgndef > 0.) {
	    i__2 = in;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		w[j] = work[in - i__ + 1];
		++j;
/* L150: */
	    }
	} else {
	    i__2 = in;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		w[j] = -work[i__];
		++j;
/* L160: */
	    }
	}
	ibegin = iend + 1;
L170:
	;
    }
    *m = *n;

    return 0;

/*     End of DLARRI */

} /* dlarri_ */

