/* dbdsgr.f -- translated by f2c (version 20240504).
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

/* Common Block Declarations */

struct {
    doublereal timng[12];
} ttime_;

#define ttime_1 ttime_

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b15 = 1.;
static integer c__5 = 5;
static doublereal c_b30 = -1.;

/* Subroutine */ int dbdsgr_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info, 
	ftnlen jobz_len, ftnlen range_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j;
    doublereal t1, t2;
    integer jj;
    doublereal eps, tol, tmp;
    integer iend, itmp;
    doublereal tnrm;
    integer indd2, inde2, inddc;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    doublereal scale;
    integer indlc;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    logical wanta;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dcopy_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    logical wantl;
    integer lwmin;
    logical wantr, wantz;
    integer indde0, indde1;
    extern doublereal dlamch_(char *, ftnlen);
    logical alleig;
    integer ibegin;
    logical indeig;
    extern doublereal dsecnd_(void);
    integer iindbl;
    logical valeig;
    integer indaor, indbor;
    extern /* Subroutine */ int dlarri_(integer *, logical *, logical *, 
	    logical *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , xerbla_(char *, integer *, ftnlen);
    integer indsgl, iindwk, indgrs, indsgr;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int dlarrv_hgb_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, logical *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);
    integer ioffsu, ioffsv;
    doublereal thresh;
    integer iinspl, indwrk, liwmin, nsplit;
    logical lquery;

    /* Fortran I/O blocks */
    static cilist io___37 = { 0, 6, 0, 0, 0 };



/*  -- LAPACK computational routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1999 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Common block to return timings .. */
/*     .. */
/*     .. Array in Common .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DBDSGR computes selected singular values and, optionally, singular vectors */
/*  of a real symmetric bidiagonal matrix B. Singular values and */
/*  singular vectors can be selected by specifying either a range of values */
/*  or a range of indices for the desired singular values. The singular values */
/*  are computed by the dqds algorithm, while orthogonal eigenvectors are */
/*  computed from various ``good'' L D L^T representations (also known as */
/*  Relatively Robust Representations). Refer to DSTEGR for more information. */

/*  For more details, see "A new O(n^2) algorithm for the symmetric */
/*  tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon, */
/*  Computer Science Division Technical Report No. UCB/CSD-97-971, */
/*  UC Berkeley, May 1997. */

/*  The extension to the SVD is performed using a set of coupling */
/*  transformations relating the L D L^T representations of B^T B, B B^T */
/*  and the so-called Golub-Kahan matrix respectively. */
/*  The technique is described in "An O(n^2) algorithm for the bidiagonal */
/*  SVD", by Benedikt Grosser and Bruno Lang, to appear in LAA, 2001. */

/*  Note 1 : Currently DBDSGR is only set up to find ALL the n */
/*  singular values and corresponding singular vectors of B in O(n^2) time. */
/*  Note 2 : DSTEGR works only on machines which follow ieee-754 */
/*  floating-point standard in their handling of infinities and NaNs. */
/*  Normal execution of DSTEGR may create NaNs and infinities and hence */
/*  may abort due to a floating point exception in environments which */
/*  do not conform to the ieee standard. */
/*  Note 3 : Only couplings for so-called positive definite initial matrices */
/*  are implemented. Thus the routine may return no results for very */
/*  tight clusters of singular values. This is flagged by INFO=3, 4 or 5. */

/*  Reminder to Osni: */
/*  Note 4 : Splittings are currently not supported, flagged by INFO = 33! */


/*  Arguments */
/*  ========= */

/*  JOBZ    (input) CHARACTER*1 */
/*          = 'N':  Compute singular values only; */
/*          = 'L':  Compute singular values and left singular vectors. */
/*          = 'R':  Compute singular values and right singular vectors. */
/*          = 'A':  Compute singular values and both left and right */
/* 		   singular vectors. */

/*  RANGE   (input) CHARACTER*1 */
/*          = 'A': all singular values will be found. */
/*          = 'V': all singular values in the half-open interval (VL,VU] */
/*                 will be found. */
/*          = 'I': the IL-th through IU-th singular values will be found. */
/* ********* Only RANGE = 'A' is currently supported ********************* */

/*  N       (input) INTEGER */
/*          The order of the matrix.  N >= 0. */

/*  D       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the N diagonal elements of the bidiagonal matrix */
/*          B. On exit, D is overwritten. */

/*  E       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the (N-1) offdiagonal elements of the bidiagonal */
/*          matrix B in elements 1 to N-1 of E. E(N) need not be set on */
/*          input, but is used internally as workspace. */
/*          On exit, E is overwritten. */

/*  VL      (input) DOUBLE PRECISION */
/*  VU      (input) DOUBLE PRECISION */
/*          If RANGE='V', the lower and upper bounds of the interval to */
/*          be searched for singular values. VL < VU. */
/*          Not referenced if RANGE = 'A' or 'I'. */

/*  IL      (input) INTEGER */
/*  IU      (input) INTEGER */
/*          If RANGE='I', the indices (in decending order) of the */
/*          smallest and largest singular values to be returned. */
/*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/*          Not referenced if RANGE = 'A' or 'V'. */

/*  ABSTOL  (input) DOUBLE PRECISION */
/*          The absolute error tolerance for the */
/*          singular values/singular vectors. IF JOBZ = 'A', the singular */
/* 	   values and singular vectors output have residual norms bounded */
/* 	   ABSTOL, and the dot products between different singular vectors */
/*          are bounded by ABSTOL. For JOBZ = 'L' or JOBZ = 'R' the dot */
/* 	   products between different singular vectors are bounded by ABSTOL, */
/*          whereas residual norms cannot be computed. */
/* 	   Let T denote the Golub-Kahan matrix */
/* 	   T =   diag([D(1),E(1),...,E(N-1),D(N)],-1) */
/*              + diag([D(1),E(1),...,E(N-1),D(N)],1) */
/*                                                          [  0  B ] */
/*          resulting from a perfect shuffle permutation of [       ] . */
/*                                                          [ B^T 0 ] */
/* 	   If ABSTOL is less than N*EPS*|T|, then */
/*          N*EPS*|T| will be used in its place, where EPS is the */
/*          machine precision and |T| is the 1-norm of the tridiagonal */
/*          matrix. The singular values are computed to an accuracy of */
/*          EPS*|T| irrespective of ABSTOL. If high relative accuracy */
/*          is important, set ABSTOL to DLAMCH( 'Safe minimum' ). */
/*          See Barlow and Demmel "Computing Accurate Eigensystems of */
/*          Scaled Diagonally Dominant Matrices", LAPACK Working Note #7 */
/*          for a discussion of which matrices define their singular values */
/*          to high relative accuracy. */

/*  M       (output) INTEGER */
/*          The total number of singular values found.  0 <= M <= N. */
/*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. */

/*  W       (output) DOUBLE PRECISION array, dimension (N) */
/*          The first M elements contain the selected singular values in */
/*          descending order. */

/*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) ) */
/*          If JOBZ = 'L' or JOBZ = 'R', then if INFO = 0, the first M */
/* 	   columns of Z contain the orthonormal */
/*          left (JOBZ = 'L') or right (JOBZ = 'R') singular vectors */
/*          of the matrix  corresponding to the selected singular values, */
/*          with the i-th column of Z holding the singular vector */
/*          associated with W(i). */
/*          If JOBZ = 'A', then if INFO = 0, the first M */
/* 	   columns of Z contain the orthonormal */
/*          right singular vectors (internal matrix V) of the matrix B */
/*          in rows 1 to N, */
/* 	   whereas rows LDZ/2+1 to LDZ/2+N contain the orthonormal */
/* 	   left singular vectors (internal matrix U). */
/*          If JOBZ = 'N', then Z is not referenced. */
/*          Note: the user must ensure that at least max(1,M) columns are */
/*          supplied in the array Z; if RANGE = 'V', the exact value of M */
/*          is not known in advance and an upper bound must be used. */

/*  LDZ     (input) INTEGER */
/*          The leading dimension of the array Z.  LDZ >= 1, and if */
/*          JOBZ = 'L' or JOBZ = 'R', LDZ >= max(1,N). */
/* 	   If JOBZ = 'A', LDZ >= max(1,2*N) */
/*          Storage convention for JOBZ = 'A', LDZ = 16, M = 4, N = 6: */
/* 		row 1       VVVV:: */
/* 		            VVVV:: */
/* 		            VVVV:: */
/* 		            VVVV:: */
/* 		            VVVV:: */
/* 		row N       VVVV:: */
/* 		            :::::: */
/* 		row LDZ/2   :::::: */
/* 		row LDZ/2+1 UUUU:: */
/* 		            UUUU:: */
/* 		            UUUU:: */
/* 		            UUUU:: */
/* 		            UUUU:: */
/* 		row LDZ/2+N UUUU:: */
/* 		            :::::: */
/* 		row LDZ     :::::: */

/*  ISUPPZ  (output) INTEGER ARRAY, */
/* 	   If JOBZ = 'L' or 'R', dimension ( 2*max(1,LDZ) ) */
/* 	   If JOBZ = 'A', dimension ( 2*max(1,LDZ) ) */
/*          The support of the singular vectors in Z, i.e., the indices */
/*          indicating the nonzero elements in Z. If JOBZ = 'L' or 'R', */
/* 	   the i-th singular vector is nonzero only in elements */
/*          ISUPPZ( 2*i-1 ) through ISUPPZ( 2*i ). */
/*          If JOBZ = 'A', the i-th right singular vector is nonzero only */
/* 	   in elements ISUPPZ( 2*i-1 ) through ISUPPZ( 2*i ), and the */
/* 	   i-th left singular vector is nonzero onlyin elements */
/* 	   ISUPPZ( LDZ + 2*i-1 ) through ISUPPZ( LDZ + 2*i ). */

/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal */
/*          (and minimal) LWORK. */


/*  Reminder to Osni: */
/*  DBDSGR requires an workspace array of dimension 26*N . */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK.  LWORK >= max(1,26*N) */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK) */
/*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */

/*  LIWORK  (input) INTEGER */
/*          The dimension of the array IWORK.  LIWORK >= max(1,10*N) */

/*          If LIWORK = -1, then a workspace query is assumed; the */
/*          routine only calculates the optimal size of the IWORK array, */
/*          returns this value as the first entry of the IWORK array, and */
/*          no error message related to LIWORK is issued by XERBLA. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO = 1, internal error in DLARRI, */
/*                if INFO > 2, internal error in DLARRV. */

/* 	   Reminder to Osni: INFO = 33 - Splittings are not supported. */


/*  Further Details */
/*  =============== */

/*  Based on contributions by */
/*     Inderjit Dhillon, IBM Almaden, USA */
/*     Osni Marques, LBNL/NERSC, USA */
/*     Benedikt Grosser, BUGH Wuppertal, Germany */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     new LOGICAL variables */
/*     new INTEGER variables */
/*     changes for DOUBLE PRECISION variables */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --d__;
    --e;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    /* Function Body */
    wantl = lsame_(jobz, "L", (ftnlen)1, (ftnlen)1);
    wantr = lsame_(jobz, "R", (ftnlen)1, (ftnlen)1);
    wanta = lsame_(jobz, "A", (ftnlen)1, (ftnlen)1);
    wantz = wantl || wantr || wanta;
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

    lquery = *lwork == -1 || *liwork == -1;
    lwmin = *n * 26;
    liwmin = *n * 10;

    *info = 0;
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;

/*     The following two lines need to be removed once the */
/*     RANGE = 'V' and RANGE = 'I' options are provided. */

    } else if (valeig || indeig) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (valeig && *n > 0 && *vu <= *vl) {
	*info = -7;
    } else if (indeig && (*il < 1 || *il > *n)) {
	*info = -8;
    } else if (indeig && (*iu < *il || *iu > *n)) {
	*info = -9;
    } else if (*ldz < 1 || (wantz && *ldz < *n || wanta && *ldz < *n << 1)) {
	*info = -14;
    } else if (*lwork < lwmin && ! lquery) {
	*info = -17;
    } else if (*liwork < liwmin && ! lquery) {
	*info = -19;
    }
    if (*info == 0) {
	work[1] = (doublereal) lwmin;
	iwork[1] = liwmin;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DBDSGR", &i__1, (ftnlen)6);
	return 0;
    } else if (lquery) {
	return 0;
    }

    if (wanta) {
	ioffsv = 0;
	ioffsu = *ldz / 2;
    } else if (wantl) {
	ioffsu = 0;
	ioffsv = 0;
    } else {
	ioffsv = 0;
	ioffsu = 0;
    }

/*     Quick return if possible */

    *m = 0;
    if (*n == 0) {
	return 0;
    }

    if (*n == 1) {
	if (alleig || indeig) {
	    *m = 1;
	    w[1] = d__[1];
	} else {
	    if (*vl < d__[1] && *vu >= d__[1]) {
		*m = 1;
		w[1] = d__[1];
	    }
	}
	if (wantl || wanta) {
	    z__[ioffsu + 1 + z_dim1] = 1.;
	}
	if (wantr || wanta) {
	    z__[ioffsv + 1 + z_dim1] = 1.;
	}
	return 0;
    }

/*     Get machine constants. */

    eps = dlamch_("Precision", (ftnlen)9);

/*     Scale matrix. */

    tnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
    scale = 1. / tnrm;
    if (scale != 1.) {
	dscal_(n, &scale, &d__[1], &c__1);
	i__1 = *n - 1;
	dscal_(&i__1, &scale, &e[1], &c__1);
	tnrm *= scale;
    }

    t1 = dsecnd_();

/*     Set pointers. */

    indgrs = 1;
    indsgl = (*n << 1) + 1;
    indsgr = *n * 3 + 1;
    indaor = (*n << 2) + 1;
    indbor = *n * 5 + 1;
    indd2 = *n * 6 + 1;
    inde2 = *n * 7 + 1;
    indde0 = (*n << 3) + 1;
    indde1 = *n * 9 + 1;
    inddc = *n * 10 + 1;
    indlc = *n * 11 + 1;
    indwrk = *n * 12 + 1;

    iinspl = 1;
    iindbl = *n + 1;
    iindwk = (*n << 1) + 1;

/*     Compute signs of B's elements. */

    work[indsgr] = 1.;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[indsgl - 1 + i__] = d_sign(&c_b15, &d__[i__]);
	d__[i__] = work[indsgl - 1 + i__] * d__[i__];
	e[i__] = work[indsgl - 1 + i__] * e[i__];
	work[indsgr + i__] = d_sign(&c_b15, &e[i__]);
	e[i__] = work[indsgr + i__] * e[i__];
	d__[i__ + 1] = work[indsgr + i__] * d__[i__ + 1];
/* L5: */
    }
    work[indsgl - 1 + *n] = d_sign(&c_b15, &d__[*n]);
    d__[*n] = work[indsgl - 1 + i__] * d__[*n];

    dcopy_(n, &d__[1], &c__1, &work[indaor], &c__1);
    i__1 = *n - 1;
    dcopy_(&i__1, &e[1], &c__1, &work[indbor], &c__1);

/*     Compute the desired singular values of the bidiagonal B. */
/*     Form correponding Gerschgorin intervals  WORK( INDGRS ), */
/*     and some auxiliary variables WORK( INDD2 ), WORK( INDE2 ), */
/*     WORK( INDDE0 ), WORK( INDDE1 ) needed later. */


/*     Reminder to Osni: */
/*     New routine DLARRI instead of DLARRE. */
/*     Since splittings are not supported, this routine may return IINFO=33 */

    thresh = eps * tnrm;
    dlarri_(n, &wantl, &wantr, &wanta, &d__[1], &e[1], &thresh, &nsplit, &
	    iwork[iinspl], m, &w[1], &work[indgrs], &work[indd2], &work[inde2]
	    , &work[indde0], &work[indde1], &work[indwrk], &iinfo);
    t2 = dsecnd_();
    ttime_1.timng[6] += t2 - t1;
    s_wsle(&io___37);
    do_lio(&c__5, &c__1, (char *)&ttime_1.timng[6], (ftnlen)sizeof(doublereal)
	    );
    e_wsle();
    if (iinfo != 0) {
	*info = iinfo;
	return 0;
    }

    if (wantz) {

/*        Compute the desired singular vectors corresponding to the computed */
/*        singular values */

	ibegin = 1;
	i__1 = nsplit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iend = iwork[iinspl + i__ - 1];
	    i__2 = iend;
	    for (j = ibegin; j <= i__2; ++j) {
		iwork[iindbl + j - 1] = i__;
/* L10: */
	    }
	    ibegin = iend + 1;
/* L20: */
	}

/* Computing MAX */
	d__1 = *abstol, d__2 = (doublereal) (*n) * eps;
	tol = max(d__1,d__2);
	t1 = dsecnd_();

/*        Reminder to Osni: Interface changed. */
/*        Since DLARRV needs WORK( INDWRK:INDWRK+14*N-1 ) as workspace */
/* 	 and INDWRK = 12*N + 1, the overall workspace requirement */
/* 	 sums up to 26*N . */

	dlarrv_hgb_(n, &d__[1], &e[1], &iwork[iinspl], m, &w[1], &iwork[iindbl], &
		work[indgrs], &tol, &wanta, &z__[z_offset], ldz, &isuppz[1], &
		work[indaor], &work[indbor], &work[indd2], &work[inde2], &
		work[indde0], &work[indde1], &work[inddc], &work[indlc], &
		work[indwrk], &iwork[iindwk], &iinfo);
	t2 = dsecnd_();
	ttime_1.timng[7] += t2 - t1;
	if (iinfo != 0) {
	    *info = 2;
	    return 0;
	}

    }

    ibegin = 1;
    i__1 = nsplit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iend = iwork[iinspl + i__ - 1];
	i__2 = iend;
	for (j = ibegin; j <= i__2; ++j) {
	    w[j] = sqrt(w[j]);
/* L30: */
	}
	ibegin = iend + 1;
/* L40: */
    }

/*     If matrix was scaled, then rescale singular values appropriately. */

    if (scale != 1.) {
	d__1 = 1. / scale;
	dscal_(m, &d__1, &w[1], &c__1);
    }

/*     Arrange signs. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (work[indsgl - 1 + i__] != 1. && (wantl || wanta)) {
	    dscal_(n, &c_b30, &z__[ioffsu + i__ + z_dim1], ldz);
	}
	if (work[indsgr - 1 + i__] != 1. && (wantr || wanta)) {
	    dscal_(n, &c_b30, &z__[ioffsv + i__ + z_dim1], ldz);
	}
/* L45: */
    }

/*     Sort singular values to descending order, along with */
/*     singular vectors. */

    i__1 = *m - 1;
    for (j = 1; j <= i__1; ++j) {
	i__ = 0;
	tmp = w[j];
	i__2 = *m;
	for (jj = j + 1; jj <= i__2; ++jj) {
	    if (w[jj] > tmp) {
		i__ = jj;
		tmp = w[jj];
	    }
/* L50: */
	}
	if (i__ != 0) {
	    w[i__] = w[j];
	    w[j] = tmp;
	    if (wantz) {
		dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
		itmp = isuppz[(i__ << 1) - 1];
		isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
		isuppz[(j << 1) - 1] = itmp;
		itmp = isuppz[i__ * 2];
		isuppz[i__ * 2] = isuppz[j * 2];
		isuppz[j * 2] = itmp;
		if (wanta) {

/* 		    Assumes that left singular vectors */
/* 		    are stored in Z(IOFFSU+1:IOFFSU+N,1:N) */

		    dswap_(n, &z__[ioffsu + 1 + i__ * z_dim1], &c__1, &z__[
			    ioffsu + 1 + j * z_dim1], &c__1);
		    itmp = isuppz[ioffsu + (i__ << 1) - 1];
		    isuppz[ioffsu + (i__ << 1) - 1] = isuppz[ioffsu + (j << 1)
			     - 1];
		    isuppz[ioffsu + (j << 1) - 1] = itmp;
		    itmp = isuppz[ioffsu + (i__ << 1)];
		    isuppz[ioffsu + (i__ << 1)] = isuppz[ioffsu + (j << 1)];
		    isuppz[ioffsu + (j << 1)] = itmp;
		}
	    }
	}
/* L60: */
    }

    work[1] = (doublereal) lwmin;
    iwork[1] = liwmin;
    return 0;

/*     End of DBDSGR */

} /* dbdsgr_ */

