/* dlas2p.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlas2p_(integer *n, integer *b1, integer *r2, doublereal 
	*a2, doublereal *b2, doublereal *ab00, doublereal *ab10, doublereal *
	lfact, doublereal *sfact, doublereal *mu, logical *sawnan, doublereal 
	*lcoup, doublereal *pcoup)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    doublereal t, dfact, dcoup;


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1999 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */

/*  Purpose */
/*  ======= */

/*  Perform the coupling transformation */
/*  [LFACT, SFACT, MU] --> [LCOUP, PCOUP]. */
/*  LFACT and SFACT are the data form the explicit factorization */
/*  B^T B - MU * I = L D L^T - MU * I = LFACT DFACT LFACT^T */
/*  using the dstqds algorithm (SFACT contains auxiliary variables). */
/*  The entries of the corresponding representation */
/*  B B^T - MU * I = LCOUP DCOUP LCOUP^T */
/*  as well as the auxiliary variable PCOUP are computed using couplings. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the matrix. */

/*  B1      (input) INTEGER */
/*          First index of splitted submatrix. */

/*  R2      (input) INTEGER */
/*          Largest index for moderate twist elements. */

/*  A2      (input) DOUBLE PRECISION array, dimension (N) */
/*          Squared diagonal elements of the L D L^T representation. */

/*  B2      (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Squared off-diagonal elements of the L D L^T representation. */

/*  AB00    (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of B^T B . */

/*  AB10    (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of B B^T . */

/*  LFACT   (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of the lower unit bidiagonal */
/*          factor computed using the dstqds algorithm. */

/*  SFACT   (input) DOUBLE PRECISION array, dimension (N) */
/*          Auxiliary variables from the dstqds algorithm. */

/*  MU      (input) DOUBLE PRECISION */
/*          The shift parameter. */

/*  SAWNAN  (input) LOGICAL */
/*          Indicates if a breakdown occured. */

/*  LCOUP   (output) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of the corresponding lower unit */
/*          bidiagonal factor computed using coupling transformations. */

/*  PCOUP   (output) DOUBLE PRECISION array, dimension (N) */
/* 	   Corresponding auxiliary variables by coupling transformations. */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */
/*     .. */
    /* Parameter adjustments */
    --pcoup;
    --lcoup;
    --sfact;
    --lfact;
    --ab10;
    --ab00;
    --b2;
    --a2;

    /* Function Body */
    if (*sawnan) {

/*         Run slower, NaN-proof code */

	i__ = *b1;
L10:
	if (i__ <= *r2 - 1) {
	    dfact = ab00[i__] / lfact[i__];
	    if (dfact == 0.) {
		dcoup = b2[i__];
		pcoup[i__] = 0.;
		lcoup[i__] = ab10[i__] / dcoup;

/*             DFACT(R2-1).EQ.ZERO not possible */

		if (i__ < *r2 - 1) {
		    dcoup = sfact[i__ + 2] - *mu;
		    pcoup[i__ + 1] = -(*mu);
		    lcoup[i__ + 1] = ab10[i__ + 1] / dcoup;
		    i__ += 2;
		} else {
		    dcoup = -(*mu);
		    pcoup[*r2] = dcoup;
		    goto L20;
		}
	    } else {
		dfact = ab00[i__] / lfact[i__];
		t = dfact / (sfact[i__] - *mu);
		dcoup = (sfact[i__ + 1] - *mu) * t;
		pcoup[i__] = -(*mu) * t;
		lcoup[i__] = ab10[i__] / dcoup;
		++i__;
	    }
	    goto L10;
	}
L20:
	;
    } else {

/*         Run normal code */

	i__1 = *r2 - 1;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    dfact = ab00[i__] / lfact[i__];
	    t = dfact / (sfact[i__] - *mu);
	    dcoup = (sfact[i__ + 1] - *mu) * t;
	    pcoup[i__] = -(*mu) * t;
	    lcoup[i__] = ab10[i__] / dcoup;
/* L30: */
	}
	if (*r2 == *n) {
	    dfact = sfact[*n] - *mu + a2[*n];
	    pcoup[*n] = -(*mu) * dfact / (sfact[*n] - *mu);
	}
    }

    return 0;
} /* dlas2p_ */

