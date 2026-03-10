/* dlap2s.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlap2s_(integer *r1, integer *bn, doublereal *a2, 
	doublereal *ab00, doublereal *ab10, doublereal *ufact, doublereal *
	pfact, doublereal *mu, logical *sawnan, doublereal *ucoup, doublereal 
	*scoup)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    doublereal t, rfact, rcoup;


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
/*  [UFACT, PFACT, MU] --> [UCOUP, UCOUP]. */
/*  UFACT and PFACT are the data form the explicit factorization */
/*  B^T B - MU * I = L D L^T - MU * I = UFACT RFACT UFACT^T */
/*  using the dqds algorithm (UFACT contains auxiliary variables). */
/*  The entries of the corresponding representation */
/*  B B^T - MU * I = UCOUP RCOUP UCOUP^T */
/*  as well as the auxiliary variable SCOUP are computed using couplings. */

/*  Arguments */
/*  ========= */

/*  R1      (input) INTEGER */
/*          Smallest index for moderate twist elements. */

/*  BN      (input) INTEGER */
/*          Last index of splitted submatrix. */

/*  A2      (input) DOUBLE PRECISION array, dimension (N) */
/*          Squared diagonal elements of the L D L^T representation. */

/*  AB00    (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of B^T B . */

/*  AB10    (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of B B^T . */

/*  UFACT   (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of the upper unit bidiagonal */
/*          factor computed using the dqds algorithm. */

/*  PFACT   (input) DOUBLE PRECISION array, dimension (N) */
/*          Auxiliary variables from the dqds algorithm. */

/*  MU      (input) DOUBLE PRECISION */
/*          The shift parameter. */

/*  SAWNAN  (input) LOGICAL */
/*          Indicates if a breakdown occured. */

/*  UCOUP   (output) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of the corresponding upper unit */
/*          bidiagonal factor computed using coupling transformations. */

/*  SCOUP   (output) DOUBLE PRECISION array, dimension (N) */
/* 	   Corresponding auxiliary variables by coupling transformations. */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */
/*     .. */
    /* Parameter adjustments */
    --scoup;
    --ucoup;
    --pfact;
    --ufact;
    --ab10;
    --ab00;
    --a2;

    /* Function Body */
    if (*sawnan) {

/*         Run slower, NaN-proof code */

	rcoup = pfact[*bn];
	scoup[*bn] = -(*mu);
	t = 1. / pfact[*bn];
	i__ = *bn - 1;
L10:
	if (i__ >= *r1) {
	    rfact = ab00[i__] / ufact[i__];
	    if (rfact == 0.) {
		ucoup[i__] = ab10[i__] / rcoup;
		rcoup = a2[i__];
		scoup[i__] = 0.;

/*             RFACT(2).EQ.ZERO not possible */

		if (i__ == 1) {
		    goto L20;
		}
		--i__;
		ucoup[i__] = ab10[i__] / rcoup;
		rcoup = a2[i__] - *mu;
		t = 1. / pfact[i__];
		scoup[i__] = -(*mu);
		--i__;
	    } else {
		ucoup[i__] = ab10[i__] / rcoup;
		rfact = ab00[i__] / ufact[i__];
		rcoup = rfact * pfact[i__] * t;
		t = 1. / pfact[i__];
		scoup[i__] = -(*mu) * rcoup * t;
		--i__;
	    }
	    goto L10;
	}
L20:
	;
    } else {

/*         Run normal code */

	rcoup = pfact[*bn];
	scoup[*bn] = -(*mu);
	t = 1. / pfact[*bn];
	i__1 = *r1;
	for (i__ = *bn - 1; i__ >= i__1; --i__) {
	    ucoup[i__] = ab10[i__] / rcoup;
	    rfact = ab00[i__] / ufact[i__];
	    rcoup = rfact * pfact[i__] * t;
	    t = 1. / pfact[i__];
	    scoup[i__] = -(*mu) * rcoup * t;
/* L30: */
	}
    }

    return 0;
} /* dlap2s_ */

