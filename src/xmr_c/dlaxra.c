/* dlaxra_f77.f -- translated by f2c (version 20240504).
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

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int dlaxra_(integer *n, doublereal *d__, doublereal *e, 
	integer *s, doublereal *scale, doublereal *asptol, integer *nblcks, 
	integer *abinds)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__;
    doublereal gl, gu, off, eps, rtmp;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    doublereal oflow;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal tminnz, tmaxnz;

/*     IMPLICIT NONE */





/*  Purpose */
/*  ======= */

/*    For a symmetric tridiagonal matrix T of dimension N given by its */
/*    entries in D and E, apply some preparatory steps: */

/*    (1)  Find a diagonal sign matrix S = diag(+- 1) such that */
/*         the offdiagonal entries of T := S T S are nonnegative. */

/*    (2)  Scale the matrix entries by SCALE into proper numerical range. */
/*         NOTE: Here it would be tempting to scale each block */
/*               individually after step (3), but then DLAXRI cannot */
/*               determine index ranges anymore (would have to carry */
/*               the per-block scalings all the way through). */

/*    (3)  Split the matrix into irreducible blocks. */
/*         To this end, determine a suitable tolerance ASPTOL and set */
/*         all ofdiagonal entries below that threshold to zero. */
/*         The resulting block structure is returned in NLBKCS and ABINDS. */

/*    Upon exit, these transformations are reflected in the entries of */
/*    D and E. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The dimension of the matrix T. */

/*  D       (input/output) DOUBLE PRECISION array, dimension ( N ) */
/*          Upon entry, the diagonal entries of the matrix T. */
/*          Upon exit, the diagonal entries of the transformed matrix. */

/*  E       (input) DOUBLE PRECISION array, dimension ( N ) */
/*          Upon entry, the offdiagonal entries of the matrix T. */
/*          Upon exit, the offdiagonal entries of the transformed matrix. */
/*          These will be zero where a block ends, that is, at */
/*          ABINDS(2),...,ABINDS(2*NBLCKS), and > ASPTOL elsewhere. */

/*  S       (output)  INTEGER array, dimension ( N ) */
/*          The diagonal scaling matrix S, entries are +-1. */

/*  SCALE   (output)  DOUBLE PRECISION */
/*          The scaling factor used in step (2). */

/*  ASPTOL  (output)  DOUBLE PRECISION */
/*          The tolerance used for splitting in step (3). This will be */
/*          somewhere in the area of n*eps*norm, where norm is the */
/*          2-norm of the transformed matrix after steps (1) and (2). */

/*  NBLCKS  (output) INTEGER */
/*          How many irreducible blocks the matrix T was split into, */
/*          NBLCKS >= 1. */

/*  ABINDS  (output) INTEGER array, dimension ( 2*NBLCKS ) */
/*          Indices where the blocks start and end. */
/*          For 1 <= i <= NBLCKS, the i'th irreducible block is the */
/*          principal submatrix ABLCKS(2*i-1):ABLCKS(2*i). */
/*          Purely for increased ease-of-use there is some redundancy */
/*          here since always */
/*            ABINDS(1)=1, ABINDS(2*i)+1=ABINDS(2*i+1), */
/*          and */
/*            ABINDS(2*NBLCKS)=N. */

/*  ====================================================================== */

/*     .. Declarations .. */


/*     .. Parameters .. */


/*     .. Constants .. */


/*     .. Locals .. */


/*  ===== Executable Statements ========================================== */

    /* Parameter adjustments */
    --abinds;
    --s;
    --e;
    --d__;

    /* Function Body */
    eps = dlamch_("Epsilon", (ftnlen)7);
    oflow = dlamch_("Overflow", (ftnlen)8);
/*     ---------------------------- */
/*      Build S to make E positive */
/*     ---------------------------- */
    s[1] = 1;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s[i__] * e[i__] < 0.) {
	    s[i__ + 1] = -1;
	} else {
	    s[i__ + 1] = 1;
	}
	e[i__] = (d__1 = e[i__], abs(d__1));
/* L90001: */
    }
/* L90002: */
/*     ----------------------------------- */
/*      Compute union of Gershgorin Discs */
/*     ----------------------------------- */

    gl = d__[1] - e[1];
    gu = d__[1] + e[1];
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	off = e[i__ - 1] + e[i__];
/* Computing MIN */
	d__1 = gl, d__2 = d__[i__] - off;
	gl = min(d__1,d__2);
/* Computing MAX */
	d__1 = gu, d__2 = d__[i__] + off;
	gu = max(d__1,d__2);
/* L90003: */
    }
/* L90004: */
/* Computing MIN */
    d__1 = gl, d__2 = d__[*n] - e[*n - 1];
    gl = min(d__1,d__2);
/* Computing MAX */
    d__1 = gu, d__2 = d__[*n] + e[*n - 1];
    gu = max(d__1,d__2);

/*     ----------------------------- */
/*      Set tolerance for splitting */
/*     ----------------------------- */
/*      Note: Even if every offdiagonal element is <= neps and split away, */
/*            the norm of the additive error matrix is still bounded by */
/*            2neps. Hence, a tolerance of c*neps is fine. */
/* Computing MAX */
    d__1 = abs(gl), d__2 = abs(gu), d__1 = max(d__1,d__2), d__2 = gu - gl;
    *asptol = *n * 3 * eps * max(d__1,d__2);

/*     ------- */
/*      Do it */
/*     ------- */

    e[*n] = 0.;
    *nblcks = 1;
    abinds[1] = 1;
    i__ = 1;
L90005:
    if (e[i__] <= *asptol) {
	e[i__] = 0.;
	abinds[*nblcks * 2] = i__;
	if (i__ == *n) {
	    goto L90006;
	}
	++(*nblcks);
	abinds[(*nblcks << 1) - 1] = i__ + 1;
    }
    ++i__;
    goto L90005;
L90006:
/*     ------------------ */
/*           Scaling */
/*     ------------------ */
/*     Determine min/max nonzero entries */
    tminnz = e[1];
    tmaxnz = tminnz;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	d__1 = tminnz, d__2 = e[i__];
	tminnz = min(d__1,d__2);
/* Computing MAX */
	d__1 = tmaxnz, d__2 = e[i__];
	tmaxnz = max(d__1,d__2);
/* L90007: */
    }
/* L90008: */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (d__[i__] != 0.) {
	    rtmp = (d__1 = d__[i__], abs(d__1));
	    tminnz = min(tminnz,rtmp);
	    tmaxnz = max(tmaxnz,rtmp);
	}
/* L90009: */
    }
/* L90010: */
/*     Since the matrix is split already, we have the magnitudes of all */
/*     entries within [ EPS*TMAX, TMAX ], except possibly for smaller */
/*     entries on the diagonal. There is no reason not to move this */
/*     window such that it contains one. */
    *scale = 1.;
    if (tmaxnz > oflow / 4) {
	*scale = 1. / tmaxnz;
    } else if (tminnz > 1. || tmaxnz < 1.) {
	*scale = 2 / (tminnz + tmaxnz);
    }
    if (*scale != 1.) {
	dscal_(n, scale, &d__[1], &c__1);
	i__1 = *n - 1;
	dscal_(&i__1, scale, &e[1], &c__1);
    }
    return 0;
} /* dlaxra_ */

