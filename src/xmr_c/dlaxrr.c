/* dlaxrr_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxrr_(integer *n, integer *k, integer *type__, 
	doublereal *e, doublereal *pivbase, doublereal *repr, integer *repi)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    integer i__, nb, ibb, ixg;
    doublereal bdet, undef;
    integer ixngn;
    extern doublereal dlamch_(char *, ftnlen);
    integer iylbbk, ixbdet, iyomga, ixgnsq;

/*     IMPLICIT NONE */



/*  Purpose */
/*  ======= */

/*    Set up representation data structure. */
/*    Is called within dlarrf and for breadth-first traversal also */
/*    within dlaxrv. */
/*    Assumes that the following primary data elements have been set: */
/*      G(1:n), Omega(1:n). */

/*  ====================================================================== */

/*     .. Declarations .. */


/*     .. Constants .. */


/*     .. Locals .. */


/*  ----- Executable Statements ----------------------------------------- */

    /* Parameter adjustments */
    --repi;
    --repr;
    --e;

    /* Function Body */
    ixg = 0;
    ixbdet = *n;
    ixngn = (*n << 1) + 1;
    ixgnsq = *n * 3 + 2;
    iyomga = 4;
    iylbbk = *n + 6;
    undef = dlamch_("Overflow", (ftnlen)8);

/*     Analyse block structure */
/*     set squares of the offdiagonals */

    ibb = 0;
    i__1 = *k - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = e[i__];
	repr[ixgnsq + i__] = d__1 * d__1;
	if (repi[iyomga + i__] != 0) {
	    repi[iylbbk + ibb] = i__ - 1;
	    ++ibb;
	}
/* L90001: */
    }
/* L90002: */
    if (*k == *n) {
	if (repi[iyomga + *k] != 0) {
	    repi[iylbbk + ibb] = i__ - 1;
	    ++ibb;
	}
    }
    repi[iylbbk + ibb] = *k;
    repr[ixgnsq + *k] = undef;
    ++ibb;
    if (*k == 1) {
	if (repi[iyomga + *k] != 0) {
	    repi[iylbbk + ibb] = i__ + 1;
	    ++ibb;
	}
    }
    i__1 = *n;
    for (i__ = *k + 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = e[i__ - 1];
	repr[ixgnsq + i__] = d__1 * d__1;
	if (repi[iyomga + i__] != 0) {
	    repi[iylbbk + ibb] = i__ + 1;
	    ++ibb;
	}
/* L90003: */
    }
/* L90004: */
    nb = ibb - 1;
    i__1 = *n / 2;
    for (i__ = ibb; i__ <= i__1; ++i__) {
	repi[iylbbk + i__] = 0;
/* L90005: */
    }
/* L90006: */

/*     Upper part (analogous to below) */

    i__ = 1;
    ibb = 0;
L90007:
L90009:
    if (i__ == repi[iylbbk + ibb]) {
	goto L90010;
    }
    repr[ixngn + i__] = repr[ixgnsq + i__] / repr[ixg + i__];
    repr[ixbdet + i__] = 0.;
    ++i__;
    goto L90009;
L90010:
    if (i__ == *k) {
	goto L90008;
    }
    ++ibb;
    bdet = repr[ixg + i__] * repr[ixg + i__ + 1] - repr[ixgnsq + i__];
    repr[ixngn + i__] = 0.;
    repr[ixbdet + i__] = bdet;
    if (i__ == *n - 1) {
/*           I = N-1 implies K = N */
	goto L90008;
    }
    repr[ixngn + i__ + 1] = repr[ixgnsq + i__ + 1] * repr[ixg + i__] / bdet;
    repr[ixbdet + i__ + 1] = 0.;
    i__ += 2;
    goto L90007;
L90008:

/*     The Twist */

    repr[ixngn + *k] = undef;
    repr[ixbdet + *k] = 0.;

/*     Lower part (analogous to above) */

    i__ = *n;
    ibb = nb;
L90011:
L90013:
    if (i__ == repi[iylbbk + ibb]) {
	goto L90014;
    }
    repr[ixngn + i__] = repr[ixgnsq + i__] / repr[ixg + i__];
    repr[ixbdet + i__] = 0.;
    --i__;
    goto L90013;
L90014:
    if (i__ == *k) {
	goto L90012;
    }
    --ibb;
    bdet = repr[ixg + i__] * repr[ixg + i__ - 1] - repr[ixgnsq + i__];
    repr[ixngn + i__] = 0.;
    repr[ixbdet + i__] = bdet;
    if (i__ == 2) {
/*           I = 2 implies K = 1 */
	goto L90012;
    }
    repr[ixngn + i__ - 1] = repr[ixgnsq + i__ - 1] * repr[ixg + i__] / bdet;
    repr[ixbdet + i__ - 1] = 0.;
    i__ += -2;
    goto L90011;
L90012:
    repi[1] = *type__;
    repi[3] = nb;
    repi[2] = *k;
    repi[iyomga] = 0;
    repi[iyomga + *n + 1] = 0;
    repr[ixngn] = 0.;
    repr[ixngn + *n + 1] = 0.;
    repr[(*n << 2) + 3] = *pivbase;
    return 0;
} /* dlaxrr_ */

