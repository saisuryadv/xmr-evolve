/* dlarrb.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlarrb_hgb_(integer *n, doublereal *d__, doublereal *l,
	doublereal *ld, doublereal *lld, integer *ifirst, integer *ilast, 
	doublereal *rtol1, doublereal *rtol2, doublereal *w, doublereal *wgap,
	 doublereal *werr, doublereal *work, integer *iwork, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__, j, k, p;
    doublereal s;
    integer i1, i2, kk;
    doublereal fac, gap, mid;
    integer cnt;
    doublereal eps, tmp, left;
    integer nint, prev, next, nleft;
    doublereal right, width, dplus, error;
    extern doublereal dlamch_(char *, ftnlen);
    integer nright, olnint;


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

/*  Given the relatively robust representation(RRR) L D L^T, DLARRB */
/*  does "limited" bisection to locate the eigenvalues of L D L^T, */
/*  W( IFIRST ) thru' W( ILAST ), to more accuracy. Intervals */
/*  [left, right] are maintained by storing their mid-points and */
/*  semi-widths in the arrays W and WERR respectively. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the matrix. */

/*  D       (input) DOUBLE PRECISION array, dimension (N) */
/*          The N diagonal elements of the diagonal matrix D. */

/*  L       (input) DOUBLE PRECISION array, dimension (N-1) */
/*          The (N-1) subdiagonal elements of the unit bidiagonal matrix L. */

/*  LD      (input) DOUBLE PRECISION array, dimension (N-1) */
/*          The (N-1) elements L(i)*D(i). */

/*  LLD     (input) DOUBLE PRECISION array, dimension (N-1) */
/*          The (N-1) elements L(i)*L(i)*D(i). */

/*  IFIRST  (input) INTEGER */
/*          The index of the first eigenvalue to be computed. */

/*  ILAST   (input) INTEGER */
/*          The index of the last eigenvalue to be computed. */

/*  RTOL1   (input) DOUBLE PRECISION */
/*          Tolerance for the convergence of the bisection intervals. */
/*          An interval [LEFT,RIGHT] has converged if */
/*          RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*ABS((LEFT+RIGHT)/2) ), */
/*          where GAP is the (estimated) distance to the neighboring */
/*          eigenvalue. */

/*  RTOL2   (input) DOUBLE PRECISION */
/*          Tolerance for the convergence of the bisection intervals. */
/*          An interval [LEFT,RIGHT] has converged if */
/*          RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*ABS((LEFT+RIGHT)/2) ), */
/*          where GAP is the (estimated) distance to the neighboring */
/*          eigenvalue. */

/*  W       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On input, W( IFIRST ) thru' W( ILAST ) are estimates of the */
/*          corresponding eigenvalues of L D L^T. */
/*          On output, these estimates are refined. */

/*  WGAP    (input/output) DOUBLE PRECISION array, dimension (N-1) */
/*          On input, the (estimated) gaps between consecutive */
/*          eigenvalues of L D L^T, i.e., WGAP(I) is the gap between */
/*          the (I+1)-st and I-th eigenvalues. Note that if */
/*          IFIRST.EQ.ILAST then WGAP(IFIRST) must be set to ZERO. */
/*          On output, these gaps are refined. */

/*  WERR    (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On input, WERR( IFIRST ) thru' WERR( ILAST ) are the errors */
/*          in the estimates W( IFIRST ) thru' W( ILAST ). */
/*          On output, these errors are refined. */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N) */
/*          Workspace. */

/*  IWORK   (workspace) INTEGER array, dimension (2*N) */
/*          Workspace. */

/* ****Reminder to Inder --- INFO is never set in this subroutine ****** */
/*  INFO    (output) INTEGER */
/*          Error flag. */

/*  Further Details */
/*  =============== */

/*  Based on contributions by */
/*     Inderjit Dhillon, IBM Almaden, USA */
/*     Osni Marques, LBNL/NERSC, USA */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */


/*     Check to see if any of the initial eigenvalue */
/*     estimates is acceptable. */

    /* Parameter adjustments */
    --iwork;
    --work;
    --werr;
    --wgap;
    --w;
    --lld;
    --ld;
    --l;
    --d__;

    /* Function Body */
    *info = 0;
    eps = dlamch_("Precision", (ftnlen)9);
    i1 = *ifirst;
    i2 = *ifirst;
    prev = 0;
    i__1 = *ilast;
    for (i__ = *ifirst; i__ <= i__1; ++i__) {
	if (i__ == *ifirst) {
	    gap = wgap[i__];
	} else if (i__ == *ilast) {
	    gap = wgap[i__ - 1];
	} else {
/* Computing MIN */
	    d__1 = wgap[i__ - 1], d__2 = wgap[i__];
	    gap = min(d__1,d__2);
	}
	error = werr[i__];
	k = i__ << 1;
	if (error < *rtol1 * gap) {
	    work[k - 1] = w[i__] - error;
	    work[k] = w[i__] + error;
	    iwork[k - 1] = -1;
	    if (i1 == i__) {
		++i1;
		prev = i__;
	    }
	} else {
	    iwork[k - 1] = 1;
	    i2 = i__;
	}
/* L10: */
    }

/*     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ]. */
/*     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while */
/*     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 ) */
/*     for an unconverged interval is set to the index of the next unconverged */
/*     interval, and is -1 or 0 for a converged interval. */

    i__ = i1;
    nint = 0;
L30:
    if (i__ <= i2) {
	if (iwork[(i__ << 1) - 1] == 1) {
	    fac = 1.;
	    left = w[i__] - werr[i__];

/*           Do while( CNT(LEFT).GT.I-1 ) */

L40:
	    if (i__ > i1 && left <= right) {
		left = right;
		cnt = i__ - 1;
	    } else {
		s = -left;
		cnt = 0;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    dplus = d__[j] + s;
		    s = s * lld[j] / dplus - left;
		    if (dplus < 0.) {
			++cnt;
		    }
/* L50: */
		}
		dplus = d__[*n] + s;
		if (dplus < 0.) {
		    ++cnt;
		}
		if (! (s > 0. || s < 1.)) {

/*                 Run a slower version of the above loop if a NaN is detected */

		    cnt = 0;
		    s = -left;
		    i__1 = *n - 1;
		    for (j = 1; j <= i__1; ++j) {
			dplus = d__[j] + s;
			if (dplus < 0.) {
			    ++cnt;
			}
			tmp = lld[j] / dplus;
			if (tmp == 0.) {
			    s = lld[j] - left;
			} else {
			    s = s * tmp - left;
			}
/* L55: */
		    }
		    dplus = d__[*n] + s;
		    if (dplus < 0.) {
			++cnt;
		    }
		}
		if (cnt > i__ - 1) {
		    left -= werr[i__] * fac;
		    fac *= 2.;
		    goto L40;
		}
	    }
	    nleft = cnt + 1;
	    i1 = min(i1,nleft);
	    fac = 1.;
	    right = w[i__] + werr[i__];

/*           Do while( CNT(RIGHT).LT.I ) */

L60:
	    s = -right;
	    cnt = 0;
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		dplus = d__[j] + s;
		s = s * lld[j] / dplus - right;
		if (dplus < 0.) {
		    ++cnt;
		}
/* L70: */
	    }
	    dplus = d__[*n] + s;
	    if (dplus < 0.) {
		++cnt;
	    }
	    if (! (s > 0. || s < 1.)) {

/*                 Run a slower version of the above loop if a NaN is detected */

		cnt = 0;
		s = -right;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    dplus = d__[j] + s;
		    if (dplus < 0.) {
			++cnt;
		    }
		    tmp = lld[j] / dplus;
		    if (tmp == 0.) {
			s = lld[j] - right;
		    } else {
			s = s * tmp - right;
		    }
/* L75: */
		}
		dplus = d__[*n] + s;
		if (dplus < 0.) {
		    ++cnt;
		}
	    }
	    if (cnt < i__) {
		right += werr[i__] * fac;
		fac *= 2.;
		goto L60;
	    }
	    ++nint;
	    k = nleft << 1;
	    work[k - 1] = left;
	    work[k] = right;
	    i__ = cnt + 1;
	    iwork[k - 1] = i__;
	    iwork[k] = cnt;
	    if (prev != nleft - 1) {
		work[k - 2] = left;
	    }
	    prev = nleft;
	} else {
	    right = work[i__ * 2];
	    ++iwork[k - 1];
	    prev = i__;
	    ++i__;
	}
	goto L30;
    }
    if (i__ <= *n && iwork[(i__ << 1) - 1] != -1) {
	work[(i__ << 1) - 1] = work[prev * 2];
    }

/*     Do while( NINT.GT.0 ) */

L80:
    prev = i1 - 1;
    olnint = nint;
    i__ = i1;
    i__1 = olnint;
    for (p = 1; p <= i__1; ++p) {
	k = i__ << 1;
	left = work[k - 1];
	right = work[k];
	next = iwork[k - 1];
	nright = iwork[k];
	mid = (left + right) * .5;
	width = right - mid;
/* Computing MAX */
	d__1 = abs(left), d__2 = abs(right);
	tmp = max(d__1,d__2);

/*        Check for convergence of intervals */

	gap = 0.;
	if (i__ == nright) {
	    if (prev > 0 && next <= *n) {
/* Computing MIN */
		d__1 = work[k + 1] - right, d__2 = left - work[k - 2];
		gap = min(d__1,d__2);
	    } else if (prev > 0) {
		gap = left - work[k - 2];
	    } else if (next <= *n) {
		gap = work[k + 1] - right;
	    }
	}
/* Computing MAX */
	d__1 = *rtol1 * gap, d__2 = *rtol2 * tmp;
	if (width < max(d__1,d__2)) {
	    --nint;
	    iwork[k - 1] = 0;
	    i__2 = nright;
	    for (j = i__ + 1; j <= i__2; ++j) {
		kk = j << 1;
		iwork[kk - 1] = 0;
		work[kk - 1] = left;
		work[kk] = right;
		wgap[j - 1] = 0.;
	    }
	    if (nright < *n) {
		wgap[nright] = work[(nright << 1) + 1] - right;
	    }
	    if (i1 == i__) {
		i1 = next;
	    } else if (next > 0) {
		iwork[(prev << 1) - 1] = next;
	    }
	    i__ = next;
	    goto L100;
	}
	prev = i__;

/*        Perform one bisection step */

	s = -mid;
	cnt = 0;
	i__2 = *n - 1;
	for (j = 1; j <= i__2; ++j) {
	    dplus = d__[j] + s;
	    s = s * lld[j] / dplus - mid;
	    if (dplus < 0.) {
		++cnt;
	    }
/* L90: */
	}
	dplus = d__[*n] + s;
	if (dplus < 0.) {
	    ++cnt;
	}
	if (! (s > 0. || s < 1.)) {

/*           Run a slower version of the above loop if a NaN is detected */

	    cnt = 0;
	    s = -mid;
	    i__2 = *n - 1;
	    for (j = 1; j <= i__2; ++j) {
		dplus = d__[j] + s;
		if (dplus < 0.) {
		    ++cnt;
		}
		tmp = lld[j] / dplus;
		if (tmp == 0.) {
		    s = lld[j] - mid;
		} else {
		    s = s * tmp - mid;
		}
/* L95: */
	    }
	    dplus = d__[*n] + s;
	    if (dplus < 0.) {
		++cnt;
	    }
	}
/* Computing MAX */
	i__2 = i__ - 1, i__3 = min(nright,cnt);
	cnt = max(i__2,i__3);
	if (cnt == i__ - 1) {
	    work[k - 1] = mid;
	} else if (cnt == nright) {
	    work[k] = mid;
	} else {
	    iwork[k] = cnt;
	    ++cnt;
	    iwork[k - 1] = cnt;
	    kk = cnt << 1;
	    iwork[kk - 1] = next;
	    iwork[kk] = nright;
	    work[k] = mid;
	    work[kk - 1] = mid;
	    work[kk] = right;
	    prev = cnt;
	    if (cnt - 1 > i__) {
		work[kk - 2] = mid;
	    }
	    if (cnt > *ifirst && cnt <= *ilast) {
		++nint;
	    } else if (cnt <= *ifirst) {
		i1 = cnt;
	    }
	}
	i__ = next;
L100:
	;
    }
    if (nint > 0) {
	goto L80;
    }
    i__1 = *ilast;
    for (i__ = *ifirst; i__ <= i__1; ++i__) {
	k = i__ << 1;
	if (iwork[k - 1] != -1) {
	    w[i__] = (work[k - 1] + work[k]) * .5;
	    werr[i__] = work[k] - w[i__];
	}
/* L110: */
    }

    return 0;

/*     End of DLARRB */

} /* dlarrb_hgb_ */

