/* dlats1.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlats1_(doublereal *z__, doublereal *ztz, doublereal *l, 
	doublereal *u, doublereal *ld, integer *b1, integer *r__, integer *bn,
	 integer *isuppz, logical *sawnan)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__, to;
    doublereal eps;
    integer from;
    extern doublereal dlamch_(char *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1999 */


/*  Reminder to Osni: */
/*  extracted from former dlar1v.f */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */

/*  Purpose */
/*  ======= */

/*  Solves N G N^T * Z = G(R,R) * E_R . */
/*  G is diagonal, whereas N is a twisted matrix: */
/*    N = eye(n) + diag([L(1),...,L(R-1),0,...,0]) */
/*               + diag([0,...,0,U(R),...,U(n-1)]) */
/*  The right hand side is given by G(R,R) times the R-th unit vector. */
/*  One proceeds as follows: */
/*    Z(R) = 1 */
/*    FOR I = R-1:-1:B1 */
/*      Z(I) = -( L( I )*Z( I+1 ) ) */
/*    END */
/*    FOR I = R+1:BN */
/*      Z( I ) = -( U( I-1 )*Z( I-1 ) ) */
/*    END */
/*  Thus only L and U are used instead of the symbolic matrices */
/*  G and N. */

/*  DLATS1 is ashortcut for DLA + T(wisted matrix) S(olver) 1 . */

/*  Arguments */
/*  ========= */

/*  Z       (output) DOUBLE PRECISION array, dimension (N) */
/*          The solution of the linear system, typically */
/*          a very good eigenvector approximation. */

/*  ZTZ     (output) DOUBLE PRECISION */
/*          Dot product ZTZ = Z^T Z . */

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

/*  ISUPPZ  (output) INTEGER array */
/*          Keeps the support information. The smallest and largest */
/*          indices of the non-zero elements of Z are stored in ISUPPZ(1) */
/*          and ISUPPZ(2) resp. */

/*  SAWNAN  (input) LOGICAL */
/*          Indicates if a breakdown occured. */

/*     .. Parameters .. */
/*     .. */
/*     .. External functions .. */
/*     .. */
/*     .. Local Scalars .. */

    /* Parameter adjustments */
    --isuppz;
    --ld;
    --u;
    --l;
    --z__;

    /* Function Body */
    eps = dlamch_("Precision", (ftnlen)9);
    isuppz[1] = *b1;
    isuppz[2] = *bn;
    z__[*r__] = 1.;
    *ztz = 1.;
    if (! (*sawnan)) {
	from = *r__ - 1;
/* Computing MAX */
	i__1 = *r__ - 32;
	to = max(i__1,*b1);
L120:
	if (from >= *b1) {
	    i__1 = to;
	    for (i__ = from; i__ >= i__1; --i__) {
		z__[i__] = -(l[i__] * z__[i__ + 1]);
		*ztz += z__[i__] * z__[i__];
/* L130: */
	    }
	    if ((d__1 = z__[to], abs(d__1)) <= eps && (d__2 = z__[to + 1], 
		    abs(d__2)) <= eps) {
		isuppz[1] = to + 2;
	    } else {
		from = to - 1;
/* Computing MAX */
		i__1 = to - 32;
		to = max(i__1,*b1);
		goto L120;
	    }
	}
	from = *r__ + 1;
/* Computing MIN */
	i__1 = *r__ + 32;
	to = min(i__1,*bn);
L140:
	if (from <= *bn) {
	    i__1 = to;
	    for (i__ = from; i__ <= i__1; ++i__) {
		z__[i__] = -(u[i__ - 1] * z__[i__ - 1]);
		*ztz += z__[i__] * z__[i__];
/* L150: */
	    }
	    if ((d__1 = z__[to], abs(d__1)) <= eps && (d__2 = z__[to - 1], 
		    abs(d__2)) <= eps) {
		isuppz[2] = to - 2;
	    } else {
		from = to + 1;
/* Computing MIN */
		i__1 = to + 32;
		to = min(i__1,*bn);
		goto L140;
	    }
	}
    } else {
	i__1 = *b1;
	for (i__ = *r__ - 1; i__ >= i__1; --i__) {
	    if (z__[i__ + 1] == 0.) {
		z__[i__] = -(ld[i__ + 1] / ld[i__]) * z__[i__ + 2];
	    } else {
		z__[i__] = -(l[i__] * z__[i__ + 1]);
	    }
	    if ((d__1 = z__[i__], abs(d__1)) <= eps && (d__2 = z__[i__ + 1], 
		    abs(d__2)) <= eps) {
		isuppz[1] = i__ + 2;
		goto L170;
	    }
	    *ztz += z__[i__] * z__[i__];
/* L160: */
	}
L170:
	i__1 = *bn - 1;
	for (i__ = *r__; i__ <= i__1; ++i__) {
	    if (z__[i__] == 0.) {
		z__[i__ + 1] = -(ld[i__ - 1] / ld[i__]) * z__[i__ - 1];
	    } else {
		z__[i__ + 1] = -(u[i__] * z__[i__]);
	    }
	    if ((d__1 = z__[i__], abs(d__1)) <= eps && (d__2 = z__[i__ + 1], 
		    abs(d__2)) <= eps) {
		isuppz[2] = i__ - 1;
		goto L190;
	    }
	    *ztz += z__[i__ + 1] * z__[i__ + 1];
/* L180: */
	}
L190:
	;
    }

    return 0;

/*     END OF DLATS1 */

} /* dlats1_ */

