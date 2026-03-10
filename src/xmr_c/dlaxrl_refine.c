/* dlaxrl_refine_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxrl_refine__(integer *il, integer *iu, doublereal *
	lgap, doublereal *ugap, integer *ewl_ae__, doublereal *ewl_lu__, 
	integer *i__, integer *j, doublereal *lambda, integer *xi)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer k, ja, ie;
    doublereal oldlb, oldub;

/*     IMPLICIT NONE */



/*  Purpose */
/*  ======= */

/*    Refine the interval [i,j] by bisecting it at lambda via the */
/*    inertia information xi. */

/*    Pre: [I,J] is an interval in the list, LAMBDA lies in the interior */
/*         of the current bounds. */

/*    Post: Eigenvalues [I,J] are covered by one or more intervals of */
/*         which the largest diameter is smaller than the diameter of */
/*         the original interval was. */

/*  ====================================================================== */

/*     .. Locals .. */


/*     -- Executable Statements ----------------------------------------- */


/*     Determine ew indices on the left and right of the sample */

    /* Parameter adjustments */
    --ewl_lu__;
    --ewl_ae__;

    /* Function Body */
    ie = *xi / 2;
    ja = (*xi + 1) / 2 + 1;

/*     Special handling of non-monotonicity. We could ignore the */
/*     inconsistencies, however, this might lead to an infinte */
/*     loop in bisection. */

    if (ja <= *i__ || ie < *i__ - 1) {
/*        In the second condition , the inertia indicates that an ew */
/*        from the left jumps in. */
/*        In any case this means LAMBDA can be seen as lower bound */
/*        for [I,J]. */
	ie = *i__ - 1;
	ja = *i__;
    } else if (ie >= *j || ja > *j + 1) {
/*        Analogously, we can regard LAMBDA as upper bound for [I,J]. */
	ja = *j + 1;
	ie = *j;
    }

/*     Now we can split the interval i:j with bounds lb:ub into three */
/*     intervals */
/*       [lb,lambda)     containing ews  i:ie, */
/*       [lambda,lambda] containing ew   ie+1:ja-1, and */
/*       (lambda,ub]     containing ews  ja:j, */
/*     where the respective index ranges may be empty, but */
/*     nevertheless we have */
/*       {i:ie} union {ie+1:ja-1} union {ja:j}  =  {i:j}. */

    oldlb = ewl_lu__[(*il << 1) - 1];
    oldub = ewl_lu__[*iu * 2];
    if (*i__ <= ie) {
	i__1 = ie;
	for (k = *i__; k <= i__1; ++k) {
	    ewl_ae__[k * 2] = ie;
	    ewl_lu__[k * 2] = *lambda;
/* L90001: */
	}
/* L90002: */
    }
    if (ie + 1 == ja - 1) {
	k = ie + 1;
	ewl_ae__[(k << 1) - 1] = k;
	ewl_ae__[k * 2] = k;
	ewl_lu__[(k << 1) - 1] = *lambda;
	ewl_lu__[k * 2] = *lambda;
    }
    if (ja <= *j) {
	i__1 = *j;
	for (k = ja; k <= i__1; ++k) {
	    ewl_ae__[(k << 1) - 1] = ja;
	    ewl_lu__[(k << 1) - 1] = *lambda;
/* L90003: */
	}
/* L90004: */
    }
    *lgap += ewl_lu__[(*il << 1) - 1] - oldlb;
    *ugap += oldub - ewl_lu__[*iu * 2];
    return 0;
} /* dlaxrl_refine__ */

