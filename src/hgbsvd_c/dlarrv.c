/* dlarrv.f -- translated by f2c (version 20240504).
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

static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int dlarrv_hgb_(integer *n, doublereal *d__, doublereal *l,
	integer *isplit, integer *m, doublereal *w, integer *iblock, 
	doublereal *gersch, doublereal *tol, logical *coup, doublereal *z__, 
	integer *ldz, integer *isuppz, doublereal *a, doublereal *b, 
	doublereal *a2, doublereal *b2, doublereal *ab00, doublereal *ab10, 
	doublereal *dcoup, doublereal *lcoup, doublereal *work, integer *
	iwork, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;
    logical L__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j, k, p, q;
    doublereal t1, t2;
    integer im, in;
    doublereal gap, eps, tmp;
    integer zto;
    doublereal ztz;
    integer iend, jblk;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    integer iter, temp[1], ktot;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    integer itmp1, itmp2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    integer indld;
    doublereal sigma;
    integer ndone, iinfo, iindr;
    doublereal resid;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    integer nclus;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    integer zfrom, iindc1, iindc2;
    extern /* Subroutine */ int dlar1v_hgb_(integer *, integer *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), dlats2_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, logical *, 
	    doublereal *);
    doublereal lambda;
    extern doublereal dlamch_(char *, ftnlen);
    integer ibegin;
    extern doublereal dsecnd_(void);
    integer indgap, indlld;
    extern /* Subroutine */ int dlarrb_hgb_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);
    doublereal mingma;
    extern /* Subroutine */ int dlarrc_hgb_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *);
    integer oldien, oldncl, factrl[4];
    doublereal relgap;
    extern /* Subroutine */ int dlarrf_hgb_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *);
    integer oldcls;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    integer ndepth, inderr, indlpl, iindwk, indpmn;
    extern /* Subroutine */ int dstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *);
    logical mgscls;
    integer lsbdpt, newcls, oldfst, indumn, indwrk, indspl, oldlst, ioffsu, 
	    ioffsv, newfrs, newftt;
    doublereal mgstol;
    integer nsplit;
    doublereal minrgp, nrminv;
    integer newlst;
    doublereal reltol;
    integer newsiz;
    doublereal rqcorr;
    extern /* Subroutine */ int dlacsv_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *)
	    ;


/*  -- LAPACK auxiliary routine (version 3.0) -- */
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

/*  Depending on COUP, DLARRV computes either eigenvectors of */
/*  a symmetric tridiagonal matrix T or sngular vector pairs */
/*  of a bidiagonal matrix B. If COUP is false, DLARRV computes */
/*  eigenvectors of the tridiagonal matrix */
/*  T = L D L^T given L, D and the eigenvalues of L D L^T. */
/*  The input eigenvalues should have high relative accuracy with */
/*  respect to the entries of L and D. The desired accuracy of the */
/*  eigenvectors can be specified by the input parameter TOL. */
/*  If COUP is true, DLARRV computes singular vector pairs of the */
/*  bidiagonal matrix B. The right singular vectors are computed */
/*  as eigenvectors of B^T B = L D L^T, wheras the left singular */
/*  vectors are determined using the coupled representation */
/*  B B^T = LCOUP DCOUP LCOUP^T . */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the matrix.  N >= 0. */

/*  D       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the N diagonal elements of the diagonal matrix D. */
/*          On exit, D may be overwritten. */

/*  L       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the (N-1) subdiagonal elements of the unit */
/*          bidiagonal matrix L are in elements 1 to N-1 of L. L(N) need */
/*          not be set on input, but is used internally as workspace. */
/*          On exit, L is overwritten. */

/*  ISPLIT  (input) INTEGER array, dimension (N) */
/*          The splitting points, at which T breaks up into submatrices. */
/*          The first submatrix consists of rows/columns 1 to */
/*          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1 */
/*          through ISPLIT( 2 ), etc. */

/*  M       (input) INTEGER */
/*          The total number of input eigenvalues.  0 <= M <= N. */

/*  W       (input) DOUBLE PRECISION array, dimension (N) */
/*          The first M elements of W contain the eigenvalues for */
/*          which eigenvectors are to be computed.  The eigenvalues */
/*          should be grouped by split-off block and ordered from */
/*          smallest to largest within the block ( The output array */
/*          W from DLARRE is expected here ). */
/*          Errors in W must be bounded by TOL (see above). */

/*  IBLOCK  (input) INTEGER array, dimension (N) */
/*          The submatrix indices associated with the corresponding */
/*          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to */
/*          the first submatrix from the top, =2 if W(i) belongs to */
/*          the second submatrix, etc. */

/*  TOL     (input) DOUBLE PRECISION */
/*          The absolute error tolerance for the eigenvectors. */
/*          The eigenvectors output have residual norms */
/*          bounded by TOL, and the dot products between different */
/*          eigenvectors are bounded by TOL. TOL must be at least */
/*          N*EPS*|T|, where EPS is the machine precision and |T| is */
/*          the 1-norm of the tridiagonal matrix. */

/*  COUP    (input) LOGICAL */
/*          If true, DLARRV computed singular vector pairs of */
/*          the bidiagonal matrix, otherwise the eigenvectors of */
/*          the symmetric tridiagonal matrix. */

/*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) ) */
/*          If INFO = 0 and COUP is false, the first M columns of Z */
/*          contain the orthonormal eigenvectors of the matrix T */
/*          corresponding to the input eigenvalues, with the i-th */
/*          column of Z holding the eigenvector associated with W(i). */
/*          If INFO = 0 and COUP is true, the first M columns of Z */
/*          contain the orthonormal singular vectors pairs of the */
/*          matrix B according to the following scheme: */
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
/*          Note: the user must ensure that at least max(1,M) columns are */
/*          supplied in the array Z. */

/*  LDZ     (input) INTEGER */
/*          The leading dimension of the array Z.  LDZ >= 1 */
/*          If (COUP=.false.), LDZ >= max(1,N), and */
/*          if (COUP=.true.), LDZ >= max(1,2*N). */


/*  Reminder to Osni: Dimension of ISUPPZ is changed. */

/*  ISUPPZ  (output) INTEGER ARRAY, dimension ( 2*max(1,LDZ) ) */
/*          The support of the eigenvectors in Z, i.e., the indices */
/*          indicating the nonzero elements in Z. The I-th eigenvector */
/*          is nonzero only in elements ISUPPZ( 2*I-1 ) through */
/*          ISUPPZ( 2*I ). */
/*          If COUP is true, then the I-th right singular vector */
/*          is nonzero only in elements ISUPPZ( 2*I-1 ) through */
/*          ISUPPZ( 2*I ), and the I-th left singular vector */
/*          is nonzero only in elements ISUPPZ( LDZ+2*I-1 ) through */
/*          ISUPPZ( LDZ+2*I ). */

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

/*  DCOUP   (output) DOUBLE PRECISION array, dimension (N) */
/*          Diagonal elements of the diagonal matrix of the */
/* 	   coupled representation LCOUP DCOUP LCOUP^T . */
/*          On exit, DCOUP is overwritten. */

/*  LCOUP   (output) DOUBLE PRECISION array, dimension (N-1) */
/*          Off-diagonal elements of the lower unit bidiagonal matrix */
/* 	   of the coupled representation LCOUP DCOUP LCOUP^T . */
/*          On exit, LCOUP is overwritten. */


/*  Reminder to Osni: */
/*  DLARRV now needs 14*N workspace compared to the former */
/*  version, which needed 13*N. */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (14*N) */

/*  IWORK   (workspace) INTEGER array, dimension (6*N) */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          > 0:  if INFO = 1, internal error in DLARRB */
/*                if INFO = 2, internal error in DSTEIN */
/*                if INFO = 3, NDEPTH to large */
/*                if INFO = 4, no couplings for DSTEIN branch */
/*                if INFO = 5, no couplings for Gram-Schmidt branch */

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
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Executable Statements .. */
/*     .. */
    /* Parameter adjustments */
    --d__;
    --l;
    --isplit;
    --w;
    --iblock;
    --gersch;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --a;
    --b;
    --a2;
    --b2;
    --ab00;
    --ab10;
    --dcoup;
    --lcoup;
    --work;
    --iwork;

    /* Function Body */
    inderr = *n;
    indld = *n << 1;
    indlld = *n * 3;
    indgap = *n << 2;
    indwrk = *n * 5 + 1;

    iindr = *n;
    iindc1 = *n << 1;
    iindc2 = *n * 3;
    iindwk = (*n << 2) + 1;

    eps = dlamch_("Precision", (ftnlen)9);

/*     Set some pointers for the singular vectors. */

    ioffsv = 0;
    ioffsu = *ldz / 2;
    i__1 = *n << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = 0;
/* L10: */
    }

/*     Initialize the storage for eigenvectors of L D L^T . */

    dlaset_("Full", n, n, &c_b5, &c_b5, &z__[ioffsv + 1 + z_dim1], ldz, (
	    ftnlen)4);

/*     vvvv Coupling vvvv */

    if (*coup) {

/*       Initialize the storage for eigenvectors of LCOUP DCOUP LCOUP^T . */

	dlaset_("Full", n, n, &c_b5, &c_b5, &z__[ioffsu + 1 + z_dim1], ldz, (
		ftnlen)4);
    }

/*     ^^^^ Coupling ^^^^ */


    mgstol = eps * 100.;

    nsplit = iblock[*m];
    ibegin = 1;
    i__1 = nsplit;
    for (jblk = 1; jblk <= i__1; ++jblk) {
	iend = isplit[jblk];

/*        Find the eigenvectors of the submatrix indexed IBEGIN */
/*        through IEND. */

	if (ibegin == iend) {
	    z__[ioffsv + ibegin + ibegin * z_dim1] = 1.;
	    isuppz[ioffsv + (ibegin << 1) - 1] = ibegin;
	    isuppz[ioffsv + (ibegin << 1)] = ibegin;

/* 	    vvvv Coupling vvvv */

	    if (*coup) {
		z__[ioffsu + ibegin + ibegin * z_dim1] = 1.;
		isuppz[(ioffsu << 1) + (ibegin << 1) - 1] = ibegin;
		isuppz[(ioffsu << 1) + (ibegin << 1)] = ibegin;
	    }

/* 	    ^^^^ Coupling ^^^^ */

	    ibegin = iend + 1;
	    goto L170;
	}
	oldien = ibegin - 1;
	in = iend - oldien;
/* Computing MIN */
	d__1 = .01, d__2 = 1. / (doublereal) in;
	reltol = min(d__1,d__2);
	im = in;
	dcopy_(&im, &w[ibegin], &c__1, &work[1], &c__1);
	i__2 = im - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[inderr + i__] = eps * (d__1 = work[i__], abs(d__1));
	    work[indgap + i__] = work[i__ + 1] - work[i__];
/* L30: */
	}
	work[inderr + im] = eps * (d__1 = work[im], abs(d__1));
/* Computing MAX */
	d__2 = (d__1 = work[im], abs(d__1));
	work[indgap + im] = max(d__2,eps);
	ndone = 0;

	ndepth = 0;
	lsbdpt = 1;
	nclus = 1;
	iwork[iindc1 + 1] = 1;
	iwork[iindc1 + 2] = im;

/*        While( NDONE.LT.IM ) do */

L40:
	if (ndone < im) {
	    oldncl = nclus;
	    nclus = 0;
	    lsbdpt = 1 - lsbdpt;
	    if (lsbdpt == 0) {
		oldcls = iindc1;
		newcls = iindc2;
	    } else {
		oldcls = iindc2;
		newcls = iindc1;
	    }
	    i__2 = oldncl;
	    for (i__ = 1; i__ <= i__2; ++i__) {

/*              If NDEPTH > 1, retrieve the relatively robust */
/*              representation (RRR) and perform limited bisection */
/*              (if necessary) to get approximate eigenvalues. */

		j = oldcls + (i__ << 1);
		oldfst = iwork[j - 1];
		oldlst = iwork[j];
		if (ndepth > 0) {
		    j = oldien + oldfst;
		    dcopy_(&in, &z__[ioffsv + ibegin + j * z_dim1], &c__1, &
			    d__[ibegin], &c__1);
		    i__3 = in - 1;
		    dcopy_(&i__3, &z__[ioffsv + ibegin + (j + 1) * z_dim1], &
			    c__1, &l[ibegin], &c__1);
		    dlaset_("Full", &in, &c__2, &c_b5, &c_b5, &z__[ioffsv + 
			    ibegin + j * z_dim1], ldz, (ftnlen)4);

/* 	          vvvv Coupling vvvv */

		    if (*coup) {
			dcopy_(n, &z__[ioffsu + ibegin + j * z_dim1], &c__1, &
				dcoup[ibegin], &c__1);
			i__3 = *n - 1;
			dcopy_(&i__3, &z__[ioffsu + ibegin + (j + 1) * z_dim1]
				, &c__1, &lcoup[ibegin], &c__1);
			dlaset_("Full", &in, &c__2, &c_b5, &c_b5, &z__[ioffsu 
				+ ibegin + j * z_dim1], ldz, (ftnlen)4);
		    }

/* 	          ^^^^ Coupling ^^^^ */

		}
		k = ibegin;
		i__3 = in - 1;
		for (j = 1; j <= i__3; ++j) {
		    tmp = d__[k] * l[k];
		    work[indld + j] = tmp;
		    work[indlld + j] = tmp * l[k];
		    ++k;
/* L50: */
		}
		if (ndepth > 0) {
		    t1 = dsecnd_();
		    d__1 = eps * 4.;
		    dlarrb_hgb_(&in, &d__[ibegin], &l[ibegin], &work[indld + 1], &
			    work[indlld + 1], &oldfst, &oldlst, &reltol, &
			    d__1, &work[1], &work[indgap + 1], &work[inderr + 
			    1], &work[indwrk], &iwork[iindwk], &iinfo);
		    if (iinfo != 0) {
			*info = 1;
			return 0;
		    }
		    t2 = dsecnd_();
		    ttime_1.timng[1] += t2 - t1;
		}

/*              Classify eigenvalues of the current representation (RRR) */
/*              as (i) isolated, (ii) loosely clustered or (iii) tightly */
/*              clustered */

		newfrs = oldfst;
		i__3 = oldlst;
		for (j = oldfst; j <= i__3; ++j) {
		    if (j == oldlst || work[indgap + j] >= reltol * (d__1 = 
			    work[j], abs(d__1))) {
			newlst = j;
		    } else {

/*                    continue (to the next loop) */

			relgap = work[indgap + j] / (d__1 = work[j], abs(d__1)
				);
			if (j == newfrs) {
			    minrgp = relgap;
			} else {
			    minrgp = min(minrgp,relgap);
			}
			goto L140;
		    }
		    newsiz = newlst - newfrs + 1;
		    newftt = oldien + newfrs;
		    mgscls = newsiz > 1 && newsiz <= 1 && minrgp >= mgstol;
		    if (newsiz > 1) {
			if (! mgscls) {

/*                       Find a new L D L^T representation if the cluster */
/*                       size is larger than MGSSIZ or the minimum */
/*                       relative gap within the cluster is too small. */

			    t1 = dsecnd_();

/* 	                vvvv Coupling vvvv */

			    if (*coup) {

/* 	Reminder to Osni: */
/* 	  DLARRC needs WORK( INDWRK:INDWRK+4*N-1 ) as workspace. */

				dlarrc_hgb_(&in, &d__[ibegin], &l[ibegin], &work[
					indld + 1], &newfrs, &newlst, &work[1]
					, &sigma, &z__[ioffsv + ibegin + 
					newftt * z_dim1], &z__[ioffsv + 
					ibegin + (newftt + 1) * z_dim1], &
					ndepth, &ab10[1], &z__[ioffsu + 
					ibegin + newftt * z_dim1], &z__[
					ioffsu + ibegin + (newftt + 1) * 
					z_dim1], &work[indwrk], info);
				if (*info != 0) {

/* 			    Cannot find a backward stable */
/* 			    coupling transformation. */

				    return 0;
				}
			    } else {
				dlarrf_hgb_(&in, &d__[ibegin], &l[ibegin], &work[
					indld + 1], &work[indlld + 1], &
					newfrs, &newlst, &work[1], &sigma, &
					z__[ioffsv + ibegin + newftt * z_dim1]
					, &z__[ioffsv + ibegin + (newftt + 1) 
					* z_dim1], &work[indwrk], &iwork[
					iindwk], info);
			    }

/* 	                ^^^^ Coupling ^^^^ */

			    t2 = dsecnd_();
			    ttime_1.timng[9] += t2 - t1;
			    if (*info == 0) {
				tmp = eps * abs(sigma);
				i__4 = newlst;
				for (k = newfrs; k <= i__4; ++k) {
				    work[k] -= sigma;
/* Computing MAX */
				    d__1 = work[indgap + k];
				    work[indgap + k] = max(d__1,tmp);
				    work[inderr + k] += tmp;
/* L51: */
				}
				++nclus;
				k = newcls + (nclus << 1);
				iwork[k - 1] = newfrs;
				iwork[k] = newlst;
			    } else {
				*info = 0;
				if (minrgp >= mgstol) {
				    mgscls = TRUE_;
				} else {

/* 	    		    vvvv Coupling vvvv */

				    if (*coup) {

/* 			      No couplings implemented for this */
/* 			      branch. */

					*info = 4;
					return 0;
				    }

/* 	    		    ^^^^ Coupling ^^^^ */


/*                             Call DSTEIN to process this tight cluster. */
/*                             This happens only if MINRGP <= MGSTOL */
/*                             and DLARRF returns INFO = 1. The latter */
/*                             means that a new RRR to "break" the */
/*                             cluster could not be found. */

				    work[indwrk] = d__[ibegin];
				    i__4 = in - 1;
				    for (k = 1; k <= i__4; ++k) {
					work[indwrk + k] = d__[ibegin + k] + 
						work[indlld + k];
/* L60: */
				    }
				    i__4 = newsiz;
				    for (k = 1; k <= i__4; ++k) {
					iwork[iindwk + k - 1] = 1;
/* L70: */
				    }
				    i__4 = newlst;
				    for (k = newfrs; k <= i__4; ++k) {
					isuppz[(ibegin + k << 1) - 3] = 1;
					isuppz[(ibegin + k << 1) - 2] = in;
/* L80: */
				    }
				    temp[0] = in;
				    t1 = dsecnd_();
				    dstein_(&in, &work[indwrk], &work[indld + 
					    1], &newsiz, &work[newfrs], &
					    iwork[iindwk], temp, &z__[ibegin 
					    + newftt * z_dim1], ldz, &work[
					    indwrk + in], &iwork[iindwk + in],
					     &iwork[iindwk + (in << 1)], &
					    iinfo);
				    t2 = dsecnd_();
				    ttime_1.timng[11] += t2 - t1;
				    if (iinfo != 0) {
					*info = 2;
					return 0;
				    }
				    ndone += newsiz;
				}
			    }
			}
		    }
		    if (newsiz == 1 || mgscls) {
			ktot = newftt;

/* 	Reminder to Osni: */
/* 	Interface of DLAR1V changed: Instead of */
/* 	using a single pointer to the workspace, we */
/* 	now differntiate WORK( INDLPL ), WORK( INDUMN ), */
/*       WORK( INDSPL ) and WORK( INDPMN ). */
/*       Pointers are choosen to work on */
/* 	WORK( INDWRK:INDWRK+4*N-1 ) */

			indlpl = indwrk;
			indumn = indwrk + *n;
			indspl = indwrk + (*n << 1);
			indpmn = indwrk + *n * 3;
			i__4 = newlst;
			for (k = newfrs; k <= i__4; ++k) {
			    iter = 0;
L90:
			    lambda = work[k];

/*                       Given LAMBDA, compute the eigenvector. */

			    t1 = dsecnd_();

/* 	Reminder to Osni: */
/* 	Interface of DLAR1V changed: Instead of */
/* 	using a single pointer to the workspace, we */
/* 	now differentiate WORK( INDLPL ), WORK( INDUMN ), */
/*       WORK( INDSPL ) and WORK( INDPMN ). */

			    dlar1v_hgb_(&c__1, &in, &lambda, &d__[ibegin], &l[
				    ibegin], &work[indld + 1], &work[indlld + 
				    1], &w[ibegin + k - 1], &gersch[(oldien <<
				     1) + 1], &z__[ioffsv + ibegin + ktot * 
				    z_dim1], &ztz, &mingma, &iwork[iindr + 
				    ktot], &isuppz[ioffsv + (ktot << 1) - 1], 
				    factrl, &work[indlpl], &work[indumn], &
				    work[indspl], &work[indpmn]);
			    t2 = dsecnd_();
			    ttime_1.timng[8] += t2 - t1;
			    tmp = 1. / ztz;
			    nrminv = sqrt(tmp);
			    resid = abs(mingma) * nrminv;
			    rqcorr = mingma * tmp;
			    if (k == in) {
				gap = work[indgap + k - 1];
			    } else if (k == 1) {
				gap = work[indgap + k];
			    } else {
/* Computing MIN */
				d__1 = work[indgap + k - 1], d__2 = work[
					indgap + k];
				gap = min(d__1,d__2);
			    }
			    ++iter;
			    if (resid > *tol * gap && abs(rqcorr) > eps * 4. *
				     abs(lambda)) {
				work[k] = lambda + rqcorr;
				if (iter < 8) {
				    goto L90;
				}
			    }

/* 	Reminder to Osni: */
/* 	A call to this routine heavily improves numerical orthogonality. */
/*       DLATS2 needs  WORK( INDWRK+4*N:INDWRK+5*N-1 ) as workspace. */

			    L__1 = factrl[2] == 1 || factrl[3] == 1;
			    dlats2_(&z__[ioffsv + ibegin + ktot * z_dim1], &
				    ztz, &mingma, &work[indlpl], &work[indumn]
				    , &work[indld + 1], &isuppz[ioffsv + (
				    ktot << 1) - 1], &iwork[iindr + ktot], &
				    isuppz[ioffsv + (ktot << 1)], &L__1, &
				    work[indwrk + (*n << 2)]);
			    tmp = 1. / ztz;
			    nrminv = sqrt(tmp);
			    iwork[ktot] = 1;
			    if (newsiz == 1) {
				++ndone;
			    }
			    zfrom = isuppz[(ktot << 1) - 1];
			    zto = isuppz[ktot * 2];
			    i__5 = zto - zfrom + 1;
			    dscal_(&i__5, &nrminv, &z__[ioffsv + ibegin + 
				    zfrom - 1 + ktot * z_dim1], &c__1);

/* 	    		vvvv Coupling vvvv */

			    if (*coup) {

/* 			  Determine the coupled left singular vectors. */


/* 	Reminder to Osni: */
/*       DLATS2 needs  WORK( INDWRK+4*IN:INDWRK+9*N-1 ) as workspace. */
/* 	Since INDWRK = 5*N+1, we end up with WORK( 1:14*N ) required */
/* 	workspace. */

				dlacsv_(&in, &c__1, &iwork[iindr + ktot], &in,
					 &ndepth, factrl, &a[ibegin], &b[
					ibegin], &a2[ibegin], &b2[ibegin], &
					ab00[ibegin], &ab10[ibegin], &work[
					indlpl], &work[indumn], &work[indspl],
					 &work[indpmn], &lambda, &w[ibegin + 
					k - 1], &gersch[(oldien << 1) + 1], &
					dcoup[ibegin], &lcoup[ibegin], &z__[
					ioffsv + ibegin + ktot * z_dim1], &
					z__[ioffsu + ibegin + ktot * z_dim1], 
					&isuppz[(ioffsu << 1) + (ktot << 1) - 
					1], &work[indwrk + (in << 2)]);
			    }

/* 	    		^^^^ Coupling ^^^^ */

			    ++ktot;
/* L100: */
			}
			if (newsiz > 1) {

/* 	    		vvvv Coupling vvvv */

			    if (*coup) {

/* 			  No couplings implemented for this */
/* 			  branch. */

				*info = 5;
				return 0;
			    }

/* 	    		^^^^ Coupling ^^^^ */

			    t1 = dsecnd_();
			    itmp1 = isuppz[(newftt << 1) - 1];
			    itmp2 = isuppz[newftt * 2];
			    ktot = oldien + newlst;
			    i__4 = ktot;
			    for (p = newftt + 1; p <= i__4; ++p) {
				i__5 = p - 1;
				for (q = newftt; q <= i__5; ++q) {
				    tmp = -ddot_(&in, &z__[ibegin + p * 
					    z_dim1], &c__1, &z__[ibegin + q * 
					    z_dim1], &c__1);
				    daxpy_(&in, &tmp, &z__[ibegin + q * 
					    z_dim1], &c__1, &z__[ibegin + p * 
					    z_dim1], &c__1);
/* L110: */
				}
				tmp = 1. / dnrm2_(&in, &z__[ibegin + p * 
					z_dim1], &c__1);
				dscal_(&in, &tmp, &z__[ibegin + p * z_dim1], &
					c__1);
/* Computing MIN */
				i__5 = itmp1, i__6 = isuppz[(p << 1) - 1];
				itmp1 = min(i__5,i__6);
/* Computing MAX */
				i__5 = itmp2, i__6 = isuppz[p * 2];
				itmp2 = max(i__5,i__6);
/* L120: */
			    }
			    i__4 = ktot;
			    for (p = newftt; p <= i__4; ++p) {
				isuppz[(p << 1) - 1] = itmp1;
				isuppz[p * 2] = itmp2;
/* L130: */
			    }
			    ndone += newsiz;
			    t2 = dsecnd_();
			    ttime_1.timng[2] += t2 - t1;
			}
		    }
		    newfrs = j + 1;
L140:
		    ;
		}
/* L150: */
	    }
	    ++ndepth;
	    goto L40;
	}
	j = ibegin << 1;
	i__2 = iend;
	for (i__ = ibegin; i__ <= i__2; ++i__) {
	    isuppz[ioffsv + j - 1] += oldien;
	    isuppz[ioffsv + j] += oldien;

/* 	    vvvv Coupling vvvv */

	    if (*coup) {
		isuppz[(ioffsu << 1) + j - 1] += oldien;
		isuppz[(ioffsu << 1) + j] += oldien;
	    }

/* 	    ^^^^ Coupling ^^^^ */

	    j += 2;
/* L160: */
	}
	ibegin = iend + 1;
L170:
	;
    }

    return 0;

/*     End of DLARRV */

} /* dlarrv_hgb_ */

