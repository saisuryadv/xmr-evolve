/* dlats2.f -- translated by f2c (version 20240504).
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

/* ****************************************************************************** */

/* Subroutine */ int dlats2_(doublereal *z__, doublereal *ztz, doublereal *
	mingma, doublereal *l, doublereal *u, doublereal *ld, integer *b1, 
	integer *r__, integer *bn, logical *sawnan, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, im1;


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1999 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */

/*  Purpose */
/*  ======= */

/*  Solves N G N^T * x = z for x and copies x to z . */
/*  G is diagonal with G(R,R) = MINGMA, whereas N is a twisted matrix: */
/*  N = eye(n) + diag([L(1),...,L(R-1),0,...,0]) */
/*             + diag([0,...,0,U(R),...,U(n-1)]) */
/*  Thus only MINGMA, L and U is needed instead of the symbolic */
/*  matrices G and N. */

/*  DLATS2 is ashortcut for DLA + T(wisted matrix) S(olver) 2 . */

/*  Arguments */
/*  ========= */

/*  Z       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On input, the solution computed with DLATS1 . */
/*          On output, a refined solution, typically */
/*          an improved eigenvector approximation. */

/*  ZTZ     (output) DOUBLE PRECISION */
/*          Dot product ZTZ = Z^T Z . */

/*  MINGMA  (input) DOUBLE PRECISION */
/*          The twist element. */

/*  L       (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of the lower unit bidiagonal */
/*          factor computed using the dstqds algorithm. */

/*  U       (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of the upper unit bidiagonal */
/*          factor computed using the dqds algorithm. */

/*  LD      (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of L D L^T . */

/*  B1      (input) INTEGER */
/*          First index of splitted submatrix. */

/*  R       (input) INTEGER */
/*          The twist position. */

/*  BN      (input) INTEGER */
/*          Last index of splitted submatrix. */

/*  SAWNAN  (input) LOGICAL */
/*          Indicates if a breakdown occured. */

/*  X       (workspace) DOUBLE PRECISION array, dimension (N) */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */

/*     Do not solve anything if a NaN occured. */

    /* Parameter adjustments */
    --x;
    --ld;
    --u;
    --l;
    --z__;

    /* Function Body */
    if (*sawnan || *mingma == 0.) {
	return 0;
    }

/*     Solve N*x = z . */

    i__1 = *bn;
    for (i__ = *b1; i__ <= i__1; ++i__) {
	x[i__] = 0.;
/* L10: */
    }
    if (*r__ == *bn) {
	x[*b1] = z__[*b1];
	im1 = *b1;
	i__1 = *bn;
	for (i__ = *b1 + 1; i__ <= i__1; ++i__) {
	    x[i__] = z__[i__] - l[im1] * x[im1];
	    im1 = i__;
/* L20: */
	}
    }
    if (*r__ == *b1) {
	x[*bn] = z__[*bn];
	i__1 = *r__;
	for (i__ = *bn - 1; i__ >= i__1; --i__) {
	    x[i__] = z__[i__] - u[i__] * x[i__ + 1];
/* L30: */
	}
    }
    if (*b1 < *r__ && *r__ < *bn) {
	x[*b1] = z__[*b1];
	im1 = *b1;
	i__1 = *r__ - 1;
	for (i__ = *b1 + 1; i__ <= i__1; ++i__) {
	    x[i__] = z__[i__] - l[im1] * x[im1];
	    im1 = i__;
/* L40: */
	}
	x[*bn] = z__[*bn];
	i__1 = *r__ + 1;
	for (i__ = *bn - 1; i__ >= i__1; --i__) {
	    x[i__] = z__[i__] - u[i__] * x[i__ + 1];
/* L50: */
	}
	x[*r__] = z__[*r__] - l[*r__ - 1] * x[*r__ - 1] - u[*r__] * x[*r__ + 
		1];
    }

/*     Solve G*z = x . */

    i__1 = *r__ - 1;
    for (i__ = *b1; i__ <= i__1; ++i__) {
	z__[i__] = x[i__] * l[i__] / ld[i__];
/* L60: */
    }
    z__[*r__] = x[*r__] / *mingma;
    im1 = *r__;
    i__1 = *bn;
    for (i__ = *r__ + 1; i__ <= i__1; ++i__) {
	z__[i__] = x[i__] * u[im1] / ld[im1];
	im1 = i__;
/* L70: */
    }

/*     Solve NT*x = z . */

    x[*r__] = z__[*r__];
    i__1 = *b1;
    for (i__ = *r__ - 1; i__ >= i__1; --i__) {
	x[i__] = z__[i__] - l[i__] * x[i__ + 1];
/* L80: */
    }
    im1 = *r__;
    i__1 = *bn;
    for (i__ = *r__ + 1; i__ <= i__1; ++i__) {
	x[i__] = z__[i__] - u[im1] * x[im1];
	im1 = i__;
/* L90: */
    }

/*     Copy x to z and compute Z^T Z . */

    *ztz = 0.;
    i__1 = *bn;
    for (i__ = *b1; i__ <= i__1; ++i__) {
	z__[i__] = x[i__];
	*ztz += z__[i__] * z__[i__];
/* L100: */
    }

    return 0;

/*     END OF DLATS2 */

} /* dlats2_ */

