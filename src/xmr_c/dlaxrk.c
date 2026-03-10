/* dlaxrk_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxrk_(integer *n, doublereal *d__, doublereal *esq, 
	integer *nblcks, integer *abinds, doublereal *pivmin, doublereal *
	abgers, integer *index, doublereal *lb, integer *nnl, integer *nzl, 
	integer *abxil, doublereal *ub, integer *nnu, integer *nzu, integer *
	abxiu, integer *iwork)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    doublereal dp;
    integer bil;
    doublereal mid;
    integer biu;
    doublereal eps, aux;
    integer bbeg, bend, blen, iblck, bnneg, bsing, tnneg, tnzer;
    extern doublereal dlamch_(char *, ftnlen);

/*     IMPLICIT NONE */




/*  Purpose */
/*  ======= */

/*    Do one bisection step for a symmetric tridiagonal matrix build */
/*    from irreducible subblocks. */

/*    For LB, NNL is the negcount and NZL the zero-count wrt the full */
/*    matrix, and ABXIL holds the blockwise inertias (these we can */
/*    define as 2*negc + issing since the blocks are irreducible, */
/*    cf. DLAXRN). Analogously for UB, NNU, NZU and ABXIU. */
/*    This should be true upon entry and will still be true upon exit. */

/*    Basically the routine computes the sturm count of T-mid(lb,ub) */
/*    and based on this makes the interval smaller while still containing */
/*    ew INDEX. */

/*    The Gershgorin Bounds are used to avoid unnessessary sturm count */
/*    computations, should for some block the shift fall outside them. */

/*    Preconditions: */
/*    (1)  NNL+NZL+1 <= INDEX <= NNU */
/*    (2)  LB < UB, there must be at least one fp-number between them. */
/*         Since the sturm counts are backward stable only in an absolute */
/*         sense, it is not recommended to call this routine with UB-LB */
/*         smaller than about N*EPS*||T||, because beyond that, you get */
/*         more or less random results. */

/*    Postconditions: */
/*    (1)  NNL+1 <= INDEX <= NNU+NZU */
/*    (2)  Either */
/*           LB = UB, NNL+1 <= INDEX <= NNL+NZL, NNL=NNU, NZL=NZU */
/*         or */
/*           LB < UB, NNL+NZL <= NNU */

/*  ====================================================================== */

/*     .. Declarations .. */


/*     .. Constants .. */


/*     .. Locals .. */


/*  ===== Executable Statements ========================================== */

    /* Parameter adjustments */
    --esq;
    --d__;
    --iwork;
    --abxiu;
    --abxil;
    --abgers;
    --abinds;

    /* Function Body */
    eps = dlamch_("Epsilon", (ftnlen)7);
    mid = (*lb + *ub) * .5;

/*     ===================== */
/*      Compute Sturm Count */
/*     ===================== */

    tnneg = 0;
    tnzer = 0;
    i__1 = *nblcks;
    for (iblck = 1; iblck <= i__1; ++iblck) {
	bbeg = abinds[(iblck << 1) - 1];
	bend = abinds[iblck * 2];
	blen = bend - bbeg + 1;
	bil = abxil[iblck] / 2 + 1;
	biu = (abxiu[iblck] + 1) / 2;
/*        Interval [LB,UB] contains local ews BIL:BIU of the block. */
	bnneg = 0;
	bsing = 0;
	if (mid >= abgers[iblck * 2]) {
	    bnneg = blen;
	} else if (mid <= abgers[(iblck << 1) - 1]) {
	    bnneg = 0;
	} else if (bil > biu) {
	    bnneg = abxil[iblck] / 2;
	} else {
	    aux = 0.;
	    i__ = bbeg;
L90003:
	    dp = d__[i__] - mid - aux;
	    if (abs(dp) < *pivmin) {
		dp = -(*pivmin);
	    }
	    if (dp < 0.) {
		++bnneg;
	    }
	    if (i__ == bend) {
		goto L90004;
	    }
	    aux = esq[i__] / dp;
	    ++i__;
	    goto L90003;
L90004:
	    if (dp == 0.) {
		bsing = 1;
	    }
	}
	iwork[iblck] = (bnneg << 1) + bsing;
	tnneg += bnneg;
	tnzer += bsing;
/* L90001: */
    }
/* L90002: */

/*     ================= */
/*      Bisect Interval */
/*     ================= */

/*      Based on the sturm count and the preconditions, we can partition */
/*      the interval (LB,UB) into three parts: */
/*       (LB,MID)  contains ews ...,TNNEG */
/*       [ MID ]   contains ews TNNEG+1:TNNEG+TNZER */
/*       (MID,UB)  contains ews TNNEG+TNZER+1,... */

    if (tnneg < *index) {
/*        Restrict lower bound to mid */
	*lb = mid;
	*nnl = tnneg;
	*nzl = tnzer;
	i__1 = *nblcks;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    abxil[i__] = iwork[i__];
/* L90005: */
	}
/* L90006: */
    }
    if (*index <= tnneg + tnzer) {
/*        Restrict upper bound to mid */
	*ub = mid;
	*nnu = tnneg;
	*nzu = tnzer;
	i__1 = *nblcks;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    abxiu[i__] = iwork[i__];
/* L90007: */
	}
/* L90008: */
    }
    return 0;
} /* dlaxrk_ */

