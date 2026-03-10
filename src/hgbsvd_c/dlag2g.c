/* dlag2g.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlag2g_(integer *r1, integer *r2, integer *rf, 
	doublereal *sfact, doublereal *pfact, doublereal *scoup, doublereal *
	mu, doublereal *gcoumn, integer *rc)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    doublereal eps, gfact, gcoup;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal gfacmn;


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1999 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */

/*  Purpose */
/*  ======= */

/*  After performing the coupling transformation */
/*  [LFACT, SFACT, MU] --> [LCOUP, PCOUP] (in DLAS2P) and */
/*  [UFACT, PFACT, MU] --> [UCOUP, UCOUP] (in DLAP2S) */
/*  an appropriate coupled twisted representation of B^T B - MU * I */
/*  is generated yielding twist position RC and ist value GCOUMN . */

/*  Arguments */
/*  ========= */

/*  R1      (input) INTEGER */
/*          Smallest index for moderate twist elements. */

/*  R2      (input) INTEGER */
/*          Largest index for moderate twist elements. */

/*  RF      (input) INTEGER */
/*          Choosen index for the optimal twist element */
/*          of the explicit factorization. */

/*  SFACT   (input) DOUBLE PRECISION array, dimension (N) */
/*          Auxiliary variables from the dstqds algorithm. */

/*  PFACT   (input) DOUBLE PRECISION array, dimension (N) */
/*          Auxiliary variables from the dqds algorithm. */

/*  SCOUP   (input) DOUBLE PRECISION array, dimension (N) */
/* 	   Corresponding auxiliary variables generated */
/*          by coupling transformations in DLAP2S . */

/*  MU      (input) DOUBLE PRECISION */
/*          The shift parameter. */

/*  GCOUMN  (output) DOUBLE PRECISION */
/*          An acceptable twist element of the coupled twisted */
/*          representation. */

/*  RC      (output) INTEGER */
/*          An acceptable twist position of the coupled twisted */
/*          representation. */

/*     .. */
/*     .. External Functions .. */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */
/*     .. */
    /* Parameter adjustments */
    --scoup;
    --pfact;
    --sfact;

    /* Function Body */
    eps = dlamch_("Precision", (ftnlen)9);

/*     Find out if GCOUP(RF) denotes an acceptable twist position */
/*     by checking (ABS(GCOUP(RF)).LE.TWO*ABS(GFACT(RF)) . */

    gfacmn = sfact[*rf] + pfact[*rf];
    if (gfacmn == 0.) {
	*rc = *rf;
	*gcoumn = 0.;
	return 0;
    } else {
	*gcoumn = gfacmn * scoup[*rf] / (sfact[*rf] - *mu);
	if (*gcoumn == 0.) {
	    *gcoumn = eps * (scoup[*rf] + *mu);
	}
    }
    *rc = *rf;
    if (abs(*gcoumn) <= abs(gfacmn) * 2.) {
	return 0;
    }

/*     Determine an acceptable twist position from GFACT . */

    i__1 = *r2;
    for (i__ = *r1; i__ <= i__1; ++i__) {
	gfact = sfact[i__] + pfact[i__];
	gcoup = gfact * scoup[i__] / (sfact[i__] - *mu);
	if (gcoup == 0.) {
	    gcoup = eps * (scoup[i__] + *mu);
	}
	if (abs(gcoup) < abs(*gcoumn)) {
	    *rc = i__;
	    *gcoumn = gcoup;
	    if (abs(gcoup) <= abs(gfacmn) * 2.) {
		return 0;
	    }
	}
/* L10: */
    }

    return 0;
} /* dlag2g_ */

