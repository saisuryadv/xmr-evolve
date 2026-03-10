/* dlaxrn0_f77.f -- translated by f2c (version 20240504).
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

integer dlaxrn0_(integer *n, integer *k, doublereal *d__, integer *omegad, 
	doublereal *r__, integer *omegar, doublereal *gamma)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    integer i__, xi, neg;

/*     IMPLICIT NONE */

/*  Purpose */
/*  ======= */

/*    Specialised routine to determine zero inertia from base data. */

/*  ====================================================================== */

/*     .. Constants .. */


/*     .. Declarations .. */


/*     .. Locals .. */


/*  ----- Executable Statements ------------------------------------------ */

    /* Parameter adjustments */
    --omegar;
    --r__;
    --omegad;
    --d__;

    /* Function Body */
    neg = 0;

    i__ = 1;
L90001:
    if (i__ >= *k) {
	goto L90002;
    }
    if (omegad[i__ + 1] != 0) {
	++neg;
	++i__;
    } else if (d__[i__] < 0.) {
	++neg;
    }
    ++i__;
    goto L90001;
L90002:

    i__ = *n;
L90003:
    if (i__ <= *k) {
	goto L90004;
    }
    if (omegar[i__ - 1] != 0) {
	++neg;
	--i__;
    } else if (r__[i__] < 0.) {
	++neg;
    }
    --i__;
    goto L90003;
L90004:

    xi = neg << 1;
    if (*k > 1 && omegad[*k] == 0 || *k < *n && omegar[*k] == 0) {
	if (*gamma < 0.) {
	    xi += 2;
	}
	if (*gamma == 0.) {
	    ++xi;
	}
    }

    ret_val = xi;
    return ret_val;
} /* dlaxrn0_ */

