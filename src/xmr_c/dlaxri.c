/* dlaxri_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxri_(integer *n, doublereal *d__, doublereal *e, 
	integer *wil, integer *wiu, integer *nblcks, integer *abinds, 
	doublereal *abgers, integer *abwind, doublereal *abvsep, doublereal *
	rwork, integer *iwork)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__, j, k;
    doublereal gl, xa, xb, xc, gu, xd, bgl;
    integer nna;
    doublereal bgu;
    integer nnb, nnc, nnd;
    doublereal eps;
    integer nza, nzb, nzc, nzd, bbeg, bend, blen;
    doublereal rtmp;
    integer iblck, cntab;
    doublereal sepab;
    integer cntcd;
    doublereal sepcd, sfmin;
    integer ndrop, ixesq;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal acctol;
    integer iybxia, iybxib, iybxic, iybxid;
    logical dodrop;
    extern /* Subroutine */ int dlaxrk_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *);
    doublereal esqmax, pivmin;
    integer iywork;

/*     IMPLICIT NONE */





/*  Purpose */
/*  ======= */

/*    For a symmetric tridiagonal matrix T of dimension N, given by its */
/*    entries in D and E, and split into irreducible blocks according to */
/*    NBLCKS and ABINDS, determine how the eigenpair indices WIL:WIU */
/*    are related to the individual blocks. */

/*    PARAMOUNT design guideline: For the partial case, we want consistent */
/*    eigensystems and a complexity of O(kn). The latter means that we */
/*    cannot use full information from the Gersgorin Discs, since that */
/*    would require sorting them. To guarantee consistency we need the */
/*    separators. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The dimension of the matrix */

/*  D       (input) DOUBLE PRECISION array, dimension ( N ) */
/*          The diagonal entries of the matrix T. */

/*  E       (input) DOUBLE PRECISION array, dimension ( N ) */
/*          The offdiagonal elements of the matrix, should be zero at */
/*          positions where a block ends and positive otherwise, just */
/*          as DLAXRA delivers them. */

/*  WIL     (input) INTEGER */
/*  WIU     (input) INTEGER */
/*          The index range of eigenpairs of the whole matrix T that */
/*          are to be mapped to local block indices. */

/*  NBLCKS  (input) INTEGER */
/*          How many irreducible blocks the matrix T was split into, */
/*          NBLCKS >= 1. */

/*  ABINDS  (input) INTEGER array, dimension ( 2*NBLCKS ) */
/*          For 1 <= i <= NBLCKS, the i'th irreducible block is the */
/*          principal submatrix ABLCKS(2*i-1):ABLCKS(2*i) of T. */

/*  ABGERS  (output) DOUBLE PRECISION array, dimension ( 2*NBLCKS ) */
/*          Upon exit, ( ABGERS(2*i-1), ABGERS(2*i) ) is the union */
/*          of Gersgorin Discs of the i'th irreducible block. */

/*  ABWIND  (output) INTEGER array, dimension ( 2*NBLCKS ) */
/*          Upon exit, ABWIND(2*i-1):ABWIND(2*i) are the local indices */
/*          of eigenpairs in the i'th block that, taken with respect to */
/*          the full matrix, form a part of WIL:WIU. */
/*          These ranges may be empty for some blocks, in which case */
/*          they are set to 0:-1. */

/*  ABVSEP  (output) DOUBLE PRECISION array, dimension ( 2*NBLCKS ) */
/*          Value separators for the blocks, mainly to be used by DSTEXR. */
/*          Upon exit, ( ABVSEP(2*i-1), ABVSEP(2*i) ) is set to an */
/*          interval to which computed eigenvalues should be capped */
/*          to guarantee monotonicity between different calls with */
/*          other index ranges. Note that these need not be contained */
/*          in the Gersgorin Bounds, in fact they may even be disjunct */
/*          to them. */

/*  ====================================================================== */

/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRK( */
/*    $             N, D, ESQ, NBLCKS, ABINDS, PIVMIN, ABGERS, */
/*    $             INDEX, */
/*    $             LB, NNL, NZL, ABXIL, */
/*    $             UB, NNU, NZU, ABXIU, */
/*    $             IWORK */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, NBLCKS, INDEX */
/*     INTEGER,          INTENT(IN)  ::  ABINDS(2*NBLCKS) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN */
/*     DOUBLE PRECISION, INTENT(IN)  ::  D(N), ESQ(N) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  ABGERS(2*NBLCKS) */

/*     INTEGER,          INTENT(INOUT)  ::  NNL, NZL, NNU, NZU */
/*     INTEGER,          INTENT(INOUT)  ::  ABXIL(NBLCKS), ABXIU(NBLCKS) */
/*     INTEGER,          INTENT(INOUT)  ::  IWORK(NBLCKS) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  LB, UB */

/*     END SUBROUTINE DLAXRK */
/*     END INTERFACE */

/*     .. Constants .. */


/*     .. Locals .. */


/*  ===== Executable Statements ========================================== */

    /* Parameter adjustments */
    --rwork;
    --e;
    --d__;
    --iwork;
    --abvsep;
    --abwind;
    --abgers;
    --abinds;

    /* Function Body */
    eps = dlamch_("Epsilon", (ftnlen)7);
    sfmin = dlamch_("Safe Miniumum", (ftnlen)13);
/*     For each block */
/*      - determine union of Gershgorin Discs */
/*      - init wanted indices */
    gl = d__[1];
    gu = d__[1];
    i__1 = *nblcks;
    for (iblck = 1; iblck <= i__1; ++iblck) {
	bbeg = abinds[(iblck << 1) - 1];
	bend = abinds[iblck * 2];
	blen = bend - bbeg + 1;
	bgl = d__[bbeg] - e[bbeg];
	bgu = d__[bbeg] + e[bbeg];
	i__2 = bend;
	for (i__ = bbeg + 1; i__ <= i__2; ++i__) {
	    rtmp = e[i__ - 1] + e[i__];
/* Computing MIN */
	    d__1 = bgl, d__2 = d__[i__] - rtmp;
	    bgl = min(d__1,d__2);
/* Computing MAX */
	    d__1 = bgu, d__2 = d__[i__] + rtmp;
	    bgu = max(d__1,d__2);
/* L90003: */
	}
/* L90004: */
	abgers[(iblck << 1) - 1] = bgl;
	abgers[iblck * 2] = bgu;
	gl = min(gl,bgl);
	gu = max(gu,bgu);
	abvsep[(iblck << 1) - 1] = bgl;
	abvsep[iblck * 2] = bgu;
	abwind[(iblck << 1) - 1] = 1;
	abwind[iblck * 2] = blen;
/* L90001: */
    }
/* L90002: */
/*     If the full spectrum is desired, we are done. */
    if (*wil == 1 && *wiu == *n) {
	return 0;
    }
/*     If the matrix has only one block, determining the index range */
/*     is trivial as well. However, we still need to set the separators */
/*     to guarantee monotonicity between different subsets. */
/*     ------------------------------------------------------------------- */
/*     Prepare bisection */
    esqmax = 0.;
    ixesq = 1;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = e[i__];
	rtmp = d__1 * d__1;
	rwork[ixesq - 1 + i__] = rtmp;
	esqmax = max(esqmax,rtmp);
/* L90005: */
    }
/* L90006: */
    rwork[ixesq - 1 + *n] = 0.;
/*     We refine up to absolute accuracy */
/* Computing MAX */
    d__1 = abs(gl), d__2 = abs(gu), d__1 = max(d__1,d__2), d__2 = gu - gl;
    acctol = *n * eps * max(d__1,d__2);
/*     Set minimum pivot in sturm sequence (for DLAXRK) */
    pivmin = sfmin / eps * max(esqmax,1.);
    iybxia = 1;
    iybxib = *nblcks + 1;
    iybxic = (*nblcks << 1) + 1;
    iybxid = *nblcks * 3 + 1;
    iywork = (*nblcks << 2) + 1;
    xa = gl;
    xb = gu;
/*     init lb data */
    i__1 = *nblcks;
    for (iblck = 1; iblck <= i__1; ++iblck) {
	bbeg = abinds[(iblck << 1) - 1];
	bend = abinds[iblck * 2];
	blen = bend - bbeg + 1;
	iwork[iybxia - 1 + iblck] = 0;
	iwork[iybxib - 1 + iblck] = blen << 1;
/* L90007: */
    }
/* L90008: */
    nna = 0;
    nza = 0;
    nnb = *n;
    nzb = 0;
/*     Pass 1 */
/*     Bisect this interval as long as it contains both WIL and WIU. */
/*     This is true initially, since we start with IL=1, IU=N. */

    cntab = 0;
L90009:
/*        mirror lb&ub data data to WIU */
    xc = xa;
    xd = xb;
    i__1 = *nblcks;
    for (iblck = 1; iblck <= i__1; ++iblck) {
	iwork[iybxic - 1 + iblck] = iwork[iybxia - 1 + iblck];
	iwork[iybxid - 1 + iblck] = iwork[iybxib - 1 + iblck];
/* L90011: */
    }
/* L90012: */
    nnc = nna;
    nzc = nza;
    nnd = nnb;
    nzd = nzb;
    if (xb - xa <= acctol || nna == *wil - 1 && nnb == *wiu) {
	goto L90010;
    }
    dlaxrk_(n, &d__[1], &rwork[ixesq], nblcks, &abinds[1], &pivmin, &abgers[1]
	    , wil, &xa, &nna, &nza, &iwork[iybxia], &xb, &nnb, &nzb, &iwork[
	    iybxib], &iwork[iywork]);
    ++cntab;
/*        check if we can split the interval */
    if (nnb + nzb < *wiu) {
/*           WIU is in the half that was just split away */
	xc = xb;
	i__1 = *nblcks;
	for (iblck = 1; iblck <= i__1; ++iblck) {
	    iwork[iybxic - 1 + iblck] = iwork[iybxib - 1 + iblck];
/* L90013: */
	}
/* L90014: */
	nnc = nnb;
	nzc = nzb;
	goto L90010;
    } else if (nnb + 1 <= *wiu) {
/*           The bound B is exact for WIU */
	xc = xb;
	xd = xb;
	i__1 = *nblcks;
	for (iblck = 1; iblck <= i__1; ++iblck) {
	    iwork[iybxic - 1 + iblck] = iwork[iybxib - 1 + iblck];
	    iwork[iybxid - 1 + iblck] = iwork[iybxib - 1 + iblck];
/* L90015: */
	}
/* L90016: */
	nnc = nnb;
	nzc = nzb;
	nnd = nnb;
	nzd = nzb;
	goto L90010;
    }
    goto L90009;
L90010:

/*     Now the complete data for intervals [A,B] and [C,D] are set. */
/*     Cases: Either both intervals are identical and small, or */
/*     they share one bound (the previous midpoint). */
/*     Left Side */
    cntab = 0;
L90017:
    if (xb - xa <= acctol || nna + nza + 1 >= *wil) {
	goto L90018;
    }
    dlaxrk_(n, &d__[1], &rwork[ixesq], nblcks, &abinds[1], &pivmin, &abgers[1]
	    , wil, &xa, &nna, &nza, &iwork[iybxia], &xb, &nnb, &nzb, &iwork[
	    iybxib], &iwork[iywork]);
    ++cntab;
    goto L90017;
L90018:
/*     Right Side */
    cntcd = 0;
L90019:
    if (xd - xc <= acctol || nnd <= *wiu) {
	goto L90020;
    }
    dlaxrk_(n, &d__[1], &rwork[ixesq], nblcks, &abinds[1], &pivmin, &abgers[1]
	    , wiu, &xc, &nnc, &nzc, &iwork[iybxic], &xd, &nnd, &nzd, &iwork[
	    iybxid], &iwork[iywork]);
    ++cntcd;
    goto L90019;
L90020:
/*     Set index ranges & separators per block */
/*     If now the intervals contain ews left of WIL or right of WIU, */
/*     then the interval was refined to full precision and we have some */
/*     very close (or even multiple) eigenvalues. Then we have to drop */
/*     indices from the per-block ranges. */
    i__1 = *nblcks;
    for (iblck = 1; iblck <= i__1; ++iblck) {
	abwind[(iblck << 1) - 1] = 0;
	abwind[iblck * 2] = -1;
	abvsep[(iblck << 1) - 1] = 0.;
	abvsep[iblck * 2] = 0.;
/* L90021: */
    }
/* L90022: */
/*     Left Side, [A,B] contains ews NNA+1:NNB+NZB of the full matrix */
/* Computing MAX */
    i__1 = 0, i__2 = *wil - (nna + 1);
    ndrop = max(i__1,i__2);
    dodrop = ndrop > 0;
/*     The Separator */
    if (! dodrop) {
	sepab = xa;
    } else {
	sepab = (xa + xb) * .5;
    }
/*     (drop from first block to last) */
    i__1 = *nblcks;
    for (iblck = 1; iblck <= i__1; ++iblck) {
	i__ = iwork[iybxia - 1 + iblck] / 2 + 1;
	j = (iwork[iybxib - 1 + iblck] + 1) / 2;
	k = i__;
/*        local block ews I:J lie within [A,B] */
	if (i__ <= j) {
/* Computing MIN */
	    i__2 = j + 1, i__3 = i__ + ndrop;
	    k = min(i__2,i__3);
	    ndrop -= k - i__;
	}
	abwind[(iblck << 1) - 1] = k;
	abvsep[(iblck << 1) - 1] = sepab;
/* L90023: */
    }
/* L90024: */
/*     Right Side, [C,D]  contains ews NNC+1:NND+NZD of the full matrix */
/* Computing MAX */
    i__1 = 0, i__2 = nnd + nzd - *wiu;
    ndrop = max(i__1,i__2);
    dodrop = ndrop > 0;
/*     The Separator */
    if (! dodrop) {
	sepcd = xd;
    } else {
	sepcd = (xc + xd) * .5;
    }
/*     (drop from last block to first) */
    for (iblck = *nblcks; iblck >= 1; --iblck) {
	i__ = iwork[iybxic - 1 + iblck] / 2 + 1;
	j = (iwork[iybxid - 1 + iblck] + 1) / 2;
	k = j;
/*        local block ews I:J lie within [C,D] */
	if (i__ <= j) {
/* Computing MAX */
	    i__1 = i__ - 1, i__2 = j - ndrop;
	    k = max(i__1,i__2);
	    ndrop -= j - k;
	}
	abwind[iblck * 2] = k;
	abvsep[iblck * 2] = sepcd;
/* L90025: */
    }
/* L90026: */
/*     establish well-defined state for blocks without wanted ews */
    i__1 = *nblcks;
    for (iblck = 1; iblck <= i__1; ++iblck) {
	if (! (abwind[(iblck << 1) - 1] <= abwind[iblck * 2])) {
	    abwind[(iblck << 1) - 1] = 0;
	    abwind[iblck * 2] = -1;
	    abvsep[(iblck << 1) - 1] = 0.;
	    abvsep[iblck * 2] = 0.;
	}
/* L90027: */
    }
/* L90028: */
    return 0;
} /* dlaxri_ */

