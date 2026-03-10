/* dlaxrl_update_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxrl_update__(integer *il, integer *iu, doublereal *
	lgap, doublereal *ugap, integer *ewl_ae__, doublereal *ewl_lu__, 
	doublereal *lambda, integer *xi)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__, j, ix, jx;
    extern /* Subroutine */ int dlaxrl_refine__(integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *);

/*     IMPLICIT NONE */



/*  Purpose */
/*  ======= */

/*    Incorporate the information given by sample LAMBDA with inertia */
/*    XI into the list. */
/*    Effectively this looks for an interval that contains LAMBDA and */
/*    refines it, but XI is employed to speed up the search. If LAMBDA */
/*    falls outside the intervals, the outer gaps will be updated if */
/*    possible. */

/*  ====================================================================== */

/*     INTERFACE */
/*     SUBROUTINE DLAXRL_REFINE( */
/*    $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, I, J, LAMBDA, XI */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER         , INTENT(IN)  ::  IL, IU, I, J, XI */
/*     DOUBLE PRECISION, INTENT(IN)  ::  LAMBDA */

/*     INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU) */
/*     END SUBROUTINE DLAXRL_REFINE */
/*     END INTERFACE */

/*     .. Constants .. */


/*     .. Locals .. */


/*     -- Executable Statements ----------------------------------------- */

/*     Look at all intervals that might benefit from this sample */
    /* Parameter adjustments */
    --ewl_lu__;
    --ewl_ae__;

    /* Function Body */
    if (*xi > (*il - 1 << 1) - 1 && *xi < (*iu + 1 << 1) - 1) {
/* Computing MAX */
	i__1 = *il, i__2 = *xi / 2;
	ix = ewl_ae__[(max(i__1,i__2) << 1) - 1];
/* Computing MIN */
	i__1 = *iu, i__2 = (*xi + 1) / 2 + 1;
	jx = ewl_ae__[min(i__1,i__2) * 2];
	i__ = ix;
L90001:
	j = ewl_ae__[i__ * 2];
	if (ewl_lu__[(i__ << 1) - 1] < *lambda && *lambda < ewl_lu__[j * 2]) {
	    dlaxrl_refine__(il, iu, lgap, ugap, &ewl_ae__[1], &ewl_lu__[1], &
		    i__, &j, lambda, xi);
	}
	if (j == jx || *lambda <= ewl_lu__[j * 2]) {
	    goto L90002;
	}
	i__ = j + 1;
	goto L90001;
L90002:
	;
    }
/*     Update the outer gaps */
/*     This must be done here, after refining the intervals, since by that */
/*     the outer bounds may have changed. */
    if (*xi == (*il - 1 << 1) - 1 || *xi == *il - 1 << 1) {
/*        is a valid left outer bound */
/* Computing MAX */
	d__1 = *lgap, d__2 = ewl_lu__[(*il << 1) - 1] - *lambda;
	*lgap = max(d__1,d__2);
    }
    if (*xi == *iu << 1 || *xi == (*iu << 1) + 1) {
/*        is a valid right outer bound */
/* Computing MAX */
	d__1 = *ugap, d__2 = *lambda - ewl_lu__[*iu * 2];
	*ugap = max(d__1,d__2);
    }
    return 0;
} /* dlaxrl_update__ */

