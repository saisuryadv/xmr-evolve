/* dlaxrl_reset_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxrl_reset__(integer *il, integer *iu, doublereal *
	lgap, doublereal *ugap, integer *ewl_ae__, doublereal *ewl_lu__, 
	doublereal *l, integer *xil, doublereal *u, integer *xiu)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j;

/*     IMPLICIT NONE */



/*  Purpose */
/*  ======= */

/*    Reset the list. We need that L is a valid lower bound for ew IL, */
/*    and that U is a valid upper bound for ew IU. */

/*  ====================================================================== */

/*     .. Constants .. */


/*     .. Locals .. */


/*  ----- Executable Statements ----------------------------------------- */

    /* Parameter adjustments */
    --ewl_lu__;
    --ewl_ae__;

    /* Function Body */
    *lgap = 0.;
    *ugap = 0.;
    i__1 = *iu;
    for (i__ = *il; i__ <= i__1; ++i__) {
	ewl_lu__[(i__ << 1) - 1] = *l;
	ewl_lu__[i__ * 2] = *u;
	ewl_ae__[(i__ << 1) - 1] = *il;
	ewl_ae__[i__ * 2] = *iu;
/* L90001: */
    }
/* L90002: */
    if (*xil % 2 == 1) {
	j = ewl_ae__[*il * 2];
	ewl_ae__[*il * 2] = *il;
	ewl_lu__[*il * 2] = *l;
	i__1 = j;
	for (i__ = *il + 1; i__ <= i__1; ++i__) {
	    ewl_ae__[(i__ << 1) - 1] = *il + 1;
/* L90003: */
	}
/* L90004: */
    }
    if (*xiu % 2 == 1) {
	i__ = ewl_ae__[(*iu << 1) - 1];
	ewl_ae__[(*iu << 1) - 1] = *iu;
	ewl_lu__[(*iu << 1) - 1] = *u;
	i__1 = *iu - 1;
	for (j = i__; j <= i__1; ++j) {
	    ewl_ae__[j * 2] = *iu - 1;
/* L90005: */
	}
/* L90006: */
    }
    return 0;
} /* dlaxrl_reset__ */

