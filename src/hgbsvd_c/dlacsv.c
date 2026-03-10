/* dlacsv.f -- translated by f2c (version 20240504).
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

static doublereal c_b3 = 1.;
static integer c__1 = 1;

/* ****************************************************************************** */

/* Subroutine */ int dlacsv_(integer *n, integer *b1, integer *r__, integer *
	bn, integer *ndepth, integer *factrl, doublereal *a, doublereal *b, 
	doublereal *a2, doublereal *b2, doublereal *ab00, doublereal *ab10, 
	doublereal *l, doublereal *umn, doublereal *s, doublereal *p, 
	doublereal *sigma, doublereal *eval, doublereal *gersch, doublereal *
	dcoup, doublereal *lcoup, doublereal *v, doublereal *u, integer *
	isuppu, doublereal *work)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    logical L__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    integer i__, r1, r2;
    doublereal tmp;
    integer uto;
    doublereal utu;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    integer indlc, indpc, rcind, indsc, induc, ufrom;
    extern /* Subroutine */ int dlag2g_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), dlap2s_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, logical *, doublereal *, doublereal *)
	    , dlas2p_(integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, logical *, doublereal *, doublereal *)
	    , dlats1_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, logical 
	    *), dlar1v_hgb_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dlats2_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, logical *, doublereal *);
    integer inddll;
    doublereal mingma;
    logical naninp, nanins;
    integer indwrk;


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1999 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */

/*  Purpose */
/*  ======= */

/*  DLACSV determines coupled singular vectors. */
/*  If (NDEPTH.EQ.0) the routines DLAS2P, DLAP2S and DLAG2G are called */
/*  to find an appropriate twisted representation for positive definite */
/*  initial matrices using certain coupling transformations. */
/*  For (NDEPTH.GT.0) the corresponding singular vector is determined */
/*  directly from the coupled representation LCOUP DCOUP LCOUP^T . */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the matrix. */

/*  B1      (input) INTEGER */
/*          First index of splitted submatrix. */

/*  R       (input) INTEGER */
/*          Choosen index for the optimal twist element */
/*          of the explicit factorization. */

/*  BN      (input) INTEGER */
/*          Last index of splitted submatrix. */

/*  NDEPTH  (input) INTEGER */
/*          Indicates if the initial matrix is positive definite */
/*          (NDEPTH.EQ.0) or not. */

/*  FACTRL  (input) INTEGER array, dimension (4) */
/*          Contains algorithmic details from the explicit factorization, */
/*          is set by DLAR1V . */

/*  A       (input) DOUBLE PRECISION array, dimension (N) */
/*          The diagonal entries of the bidiagonal matrix B. */

/*  B       (input) DOUBLE PRECISION array, dimension (N-1) */
/*          The off-diagonal entries of the bidiagonal matrix B. */

/*  A2      (input) DOUBLE PRECISION array, dimension (N) */
/*          Squared diagonal elements of the L D L^T representation. */

/*  B2      (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Squared off-diagonal elements of the L D L^T representation. */

/*  AB00    (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of B^T B . */

/*  AB10    (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of B B^T . */

/*  L       (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of the lower unit bidiagonal */
/*          factor computed using the dstqds algorithm. */

/*  UMN     (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of the upper unit bidiagonal */
/*          factor computed using the dqds algorithm. */

/*  S       (input) DOUBLE PRECISION array, dimension (N) */
/*          Auxiliary variables from the dstqds algorithm. */

/*  P       (input) DOUBLE PRECISION array, dimension (N) */
/*          Auxiliary variables from the dqds algorithm. */

/*  SIGMA   (input) DOUBLE PRECISION */
/*          The shift parameter. */

/*  EVAL    (input) DOUBLE PRECISION */
/*          The refined eigenvalue of the explicitly factorized */
/* 	   representation L D L^T . */

/*  GERSCH  (input) DOUBLE PRECISION array, dimension (N) */
/*          The (refined) Gerschgorin intervals. */

/*  DCOUP   (input) DOUBLE PRECISION array, dimension (N) */
/*          Diagonal elements of the diagonal matrix of the */
/* 	   coupled representation LCOUP DCOUP LCOUP^T . */

/*  LCOUP   (input) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of the lower unit bidiagonal matrix */
/* 	   of the coupled representation LCOUP DCOUP LCOUP^T . */

/*  V       (input) DOUBLE PRECISION array, dimension (N) */
/*          The right singular vector (needed for arranging signs). */

/*  U       (output) DOUBLE PRECISION array, dimension (N) */
/*          The coupled left singular vector. */

/*  ISUPPU  (output) INTEGER array */
/*          Keeps the support information. The smallest and largest */
/*          indices of the non-zero elements are stored in ISUPPU(1) */
/*          and ISUPPU(2) resp. */

/*  U       (output) DOUBLE PRECISION array, dimension (N) */
/*          The coupled left singular vector. */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (5*N) */
/*          Workspace. */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */
/*     .. */
    /* Parameter adjustments */
    --work;
    --isuppu;
    --u;
    --v;
    --lcoup;
    --dcoup;
    --gersch;
    --p;
    --s;
    --umn;
    --l;
    --ab10;
    --ab00;
    --b2;
    --a2;
    --b;
    --a;
    --factrl;

    /* Function Body */
    indlc = 1;
    indpc = *n + 1;
    induc = (*n << 1) + 1;
    indsc = *n * 3 + 1;
    indwrk = (*n << 2) + 1;
    inddll = indwrk;

/*     Extract information from FACTRL */

    r1 = factrl[1];
    r2 = factrl[2];
    nanins = factrl[3] == 1;
    naninp = factrl[4] == 1;

    if (*ndepth == 0) {

/*       Determine couplings for positive definite intitial matrices. */

	dlas2p_(n, b1, &r2, &a2[1], &b2[1], &ab00[1], &ab10[1], &l[1], &s[1], 
		sigma, &nanins, &work[indlc], &work[indpc]);
	dlap2s_(&r1, bn, &a2[1], &ab00[1], &ab10[1], &umn[1], &p[1], sigma, &
		naninp, &work[induc], &work[indsc]);
	dlag2g_(&r1, &r2, r__, &s[1], &p[1], &work[indsc], sigma, &mingma, &
		rcind);

/* 	Solve linear system with a matrix given as a twisted */
/* 	factorization with twist index RCIND. */

	L__1 = nanins || naninp;
	dlats1_(&u[1], &utu, &work[indlc], &work[induc], &ab10[1], b1, &rcind,
		 bn, &isuppu[1], &L__1);
	L__1 = nanins || naninp;
	dlats2_(&u[1], &utu, &mingma, &work[indlc], &work[induc], &ab10[1], &
		isuppu[1], &rcind, &isuppu[2], &L__1, &work[indwrk]);
    } else {

/*       Compute coupled singular vectors directly from the */
/* 	coupled representation LCOUP DCOUP LCOUP^T . */

	rcind = 0;
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[inddll + i__ - 1] = ab10[i__] * lcoup[i__];
/* L10: */
	}
	dlar1v_hgb_(b1, bn, sigma, &dcoup[1], &lcoup[1], &ab10[1], &work[inddll], 
		eval, &gersch[1], &u[1], &utu, &mingma, &rcind, &isuppu[1], &
		factrl[1], &work[indlc], &work[induc], &work[indsc], &work[
		indpc]);
	L__1 = factrl[3] == 1 || factrl[4] == 1;
	dlats2_(&u[1], &utu, &mingma, &work[indlc], &work[induc], &ab10[1], &
		isuppu[1], &rcind, &isuppu[2], &L__1, &work[indwrk]);
    }

/*     Arrange signs for U and V . */

/*     Compute temp = r-th component of bidiagonal * V . */

    if (*r__ < *n) {
	tmp = a[*r__] * v[*r__] + b[*r__] * v[*r__ + 1];
    } else {
	tmp = a[*n] * v[*n];
    }

/*     Compare signs of TMP and U(R) and scale appropriately. */

    ufrom = isuppu[1];
    uto = isuppu[2];
    d__1 = tmp * u[*r__];
    if (d_sign(&c_b3, &d__1) < 0.) {
	i__1 = uto - ufrom + 1;
	d__1 = -1. / sqrt(utu);
	dscal_(&i__1, &d__1, &u[ufrom], &c__1);
    } else {
	i__1 = uto - ufrom + 1;
	d__1 = 1. / sqrt(utu);
	dscal_(&i__1, &d__1, &u[ufrom], &c__1);
    }

    return 0;
} /* dlacsv_ */

