/* dlaxro_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxro_(integer *n, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, integer *isuppz, integer *revord)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, k, itmp;
    doublereal rtmp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);

/*     IMPLICIT NONE */





/*  Purpose */
/*  ======= */

/*    Sort the triplets ( W(i), Z(:,i), ISUPPZ(2*i-1:2*i) ) into */
/*    ascending order with respect to W. */

/*    The array REVORD is set to the reverse permutation: The eigenpairs */
/*    i in the output order was the REVORD(i)'th in the input order. */

/*  ====================================================================== */

/*     .. Delcarations .. */


/*     .. Locals .. */


/*  ===== Executable Statements ========================================== */

    /* Parameter adjustments */
    --revord;
    --isuppz;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;

    /* Function Body */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	revord[i__] = i__;
/* L90001: */
    }
/* L90002: */

    i__1 = *m - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rtmp = w[i__];
	j = i__;
	i__2 = *m;
	for (k = i__ + 1; k <= i__2; ++k) {
	    if (w[k] < rtmp) {
		rtmp = w[k];
		j = k;
	    }
/* L90005: */
	}
/* L90006: */

	if (j != i__) {
	    dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1], &
		    c__1);
	    rtmp = w[i__];
	    w[i__] = w[j];
	    w[j] = rtmp;
	    itmp = isuppz[(i__ << 1) - 1];
	    isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
	    isuppz[(j << 1) - 1] = itmp;
	    itmp = isuppz[i__ * 2];
	    isuppz[i__ * 2] = isuppz[j * 2];
	    isuppz[j * 2] = itmp;
	    itmp = revord[i__];
	    revord[i__] = revord[j];
	    revord[j] = itmp;
	}
/* L90003: */
    }
/* L90004: */
    return 0;
} /* dlaxro_ */

