/* reptools_f77.f -- translated by f2c (version 20240504).
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

/*  ====================================================================== */
/*   Convenience tools for working with our compressed representation */
/*   data structure, in particular to facilitate easy integration in */
/*   C++ Grailframe. */

/*   XMR_INITREP */
/*     to set up the representation given the primary data */
/*     explicitly (instead of assuming it is already set as */
/*     DLAXRR does it) */

/*   XMR_REPSIZE_REAL */
/*     given n, returns int giving absolute size of real part */

/*   XMR_REPSIZE_INT */
/*     given n, returns int giving absolute size of real part */

/*  ====================================================================== */
/* Subroutine */ int xmr_initrep__(integer *n, integer *k, integer *type__, 
	doublereal *g, integer *omega, doublereal *e, doublereal *pivbase, 
	doublereal *repr, integer *repi)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, ixg, iyomga;
    extern /* Subroutine */ int dlaxrr_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *);

/*     IMPLICIT NONE */
/*     .. */
/*     .. Scalar arguments .. */
/*     .. */
/*     .. */
/*     .. Array arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*    Set up representation data structure from the data given by */
/*    N, K, G, Omega and E. */

/*    This is a convenience routine, where no part of REPR or REPI */
/*    is assumed to ba initialized already. */

/*    Note that OMEGA is indexed 1:n, no padding here. */

/*  ====================================================================== */

/*     .. Declarations .. */


/*     .. Locals .. */


/*  ----- Executable Statements ----------------------------------------- */

    /* Parameter adjustments */
    --repi;
    --repr;
    --e;
    --omega;
    --g;

    /* Function Body */
    ixg = 0;
    iyomga = 4;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	repr[ixg + i__] = g[i__];
	repi[iyomga + i__] = omega[i__];
/* L90001: */
    }
/* L90002: */
    dlaxrr_(n, k, type__, &e[1], pivbase, &repr[1], &repi[1]);
    return 0;
} /* xmr_initrep__ */


/*     End of subroutine XMR_INITREP */

/* *********************************************************************** */

integer xmr_repsize_real__(integer *n)
{
    /* System generated locals */
    integer ret_val;

/*     IMPLICIT NONE */
    ret_val = (*n << 2) + 3;
    return ret_val;
} /* xmr_repsize_real__ */


/* *********************************************************************** */

integer xmr_repsize_int__(integer *n)
{
    /* System generated locals */
    integer ret_val;

/*     IMPLICIT NONE */
    ret_val = *n + 6 + *n / 2;
    return ret_val;
} /* xmr_repsize_int__ */

