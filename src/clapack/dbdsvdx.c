/* dbdsvdx.f -- translated by f2c (version 20240504).
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

static doublereal c_b10 = 1.;
static doublereal c_b14 = -.125;
static integer c__1 = 1;
static doublereal c_b22 = 0.;
static integer c__2 = 2;

/* > \brief \b DBDSVDX */

/*  =========== DOCUMENTATION =========== */

/*  DBDSVDX computes the SVD of a bidiagonal matrix using TGK+DSTEVX. */
/*  Converted from LAPACK 3.6.0+ F90 source to F77 for f2c. */

/*  ===================================================================== */
/* Subroutine */ int dbdsvdx_(char *uplo, char *jobz, char *range, integer *n,
	 doublereal *d__, doublereal *e, doublereal *vl, doublereal *vu, 
	integer *il, integer *iu, integer *ns, doublereal *s, doublereal *z__,
	 integer *ldz, doublereal *work, integer *iwork, integer *info, 
	ftnlen uplo_len, ftnlen jobz_len, ftnlen range_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal), pow_dd(
	    doublereal *, doublereal *);

    /* Local variables */
    integer i__, j, k;
    doublereal mu, eps;
    integer nsl;
    doublereal tol, ulp;
    integer nru, nrv;
    doublereal emin;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    integer ntgk;
    doublereal smin, smax, nrmu, nrmv;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    logical sveq0;
    integer idbeg;
    doublereal sqrt2;
    integer idend;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    integer isbeg;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer idtgk, ietgk, iltgk, itemp;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    integer icolz;
    logical allsv;
    integer idptr;
    logical indsv;
    integer ieptr, iutgk;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    doublereal vltgk;
    logical lower;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    doublereal zjtji;
    logical split, valsv;
    integer isplt;
    doublereal ortol, vutgk;
    logical wantz;
    char rngvx[1];
    integer irowu, irowv, irowz;
    extern doublereal dlamch_(char *, ftnlen);
    integer iifail;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    doublereal abstol, thresh;
    integer iiwork;
    extern /* Subroutine */ int dstevx_(char *, char *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, integer *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK driver routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --d__;
    --e;
    --s;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;

    /* Function Body */
    allsv = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
    valsv = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
    indsv = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);

    *info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lower) {
	*info = -1;
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (allsv || valsv || indsv)) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*n > 0) {
	if (valsv) {
	    if (*vl < 0.) {
		*info = -7;
	    } else if (*vu <= *vl) {
		*info = -8;
	    }
	} else if (indsv) {
	    if (*il < 1 || *il > max(1,*n)) {
		*info = -9;
	    } else if (*iu < min(*n,*il) || *iu > *n) {
		*info = -10;
	    }
	}
    }
    if (*info == 0) {
	if (*ldz < 1 || wantz && *ldz < *n << 1) {
	    *info = -14;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DBDSVDX", &i__1, (ftnlen)7);
	return 0;
    }

/*     Quick return if possible (N.LE.1) */

    *ns = 0;
    if (*n == 0) {
	return 0;
    }

    if (*n == 1) {
	if (allsv || indsv) {
	    *ns = 1;
	    s[1] = abs(d__[1]);
	} else {
	    if (*vl < abs(d__[1]) && *vu >= abs(d__[1])) {
		*ns = 1;
		s[1] = abs(d__[1]);
	    }
	}
	if (wantz) {
	    z__[z_dim1 + 1] = d_sign(&c_b10, &d__[1]);
	    z__[z_dim1 + 2] = 1.;
	}
	return 0;
    }

    abstol = dlamch_("Safe Minimum", (ftnlen)12) * 2;
    ulp = dlamch_("Precision", (ftnlen)9);
    eps = dlamch_("Epsilon", (ftnlen)7);
    sqrt2 = sqrt(2.);
    ortol = sqrt(ulp);

/*     Criterion for splitting is taken from DBDSQR when singular */
/*     values are computed to relative accuracy TOL. */

/* Computing MAX */
/* Computing MIN */
    d__3 = 100., d__4 = pow_dd(&eps, &c_b14);
    d__1 = 10., d__2 = min(d__3,d__4);
    tol = max(d__1,d__2) * eps;

/*     Compute approximate maximum, minimum singular values. */

    i__ = idamax_(n, &d__[1], &c__1);
    smax = (d__1 = d__[i__], abs(d__1));
    i__1 = *n - 1;
    i__ = idamax_(&i__1, &e[1], &c__1);
/* Computing MAX */
    d__2 = smax, d__3 = (d__1 = e[i__], abs(d__1));
    smax = max(d__2,d__3);

/*     Compute threshold for neglecting D's and E's. */

    smin = abs(d__[1]);
    if (smin != 0.) {
	mu = smin;
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    mu = (d__2 = d__[i__], abs(d__2)) * (mu / (mu + (d__1 = e[i__ - 1]
		    , abs(d__1))));
	    smin = min(smin,mu);
	    if (smin == 0.) {
		goto L15;
	    }
/* L10: */
	}
L15:
	;
    }
    smin /= sqrt((doublereal) (*n));
    thresh = tol * smin;

/*     Check for zeros in D and E (splits), i.e. submatrices. */

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = d__[i__], abs(d__1)) <= thresh) {
	    d__[i__] = 0.;
	}
	if ((d__1 = e[i__], abs(d__1)) <= thresh) {
	    e[i__] = 0.;
	}
/* L20: */
    }
    if ((d__1 = d__[*n], abs(d__1)) <= thresh) {
	d__[*n] = 0.;
    }

/*     Pointers for arrays used by DSTEVX. */

    idtgk = 1;
    ietgk = idtgk + (*n << 1);
    itemp = ietgk + (*n << 1);
    iifail = 1;
    iiwork = iifail + (*n << 1);

/*     Set RNGVX, which corresponds to RANGE for DSTEVX in TGK mode. */
/*     VL,VU or IL,IU are redefined to conform to implementation a) */
/*     described in the leading comments. */

    iltgk = 0;
    iutgk = 0;
    vltgk = 0.;
    vutgk = 0.;

    if (allsv) {

/*        All singular values will be found. We aim at -s (see */
/*        leading comments) with RNGVX = 'I'. IL and IU are set */
/*        later (as ILTGK and IUTGK) according to the dimension */
/*        of the active submatrix. */

	*(unsigned char *)rngvx = 'I';
	if (wantz) {
	    i__1 = *n << 1;
	    i__2 = *n + 1;
	    dlaset_("F", &i__1, &i__2, &c_b22, &c_b22, &z__[z_offset], ldz, (
		    ftnlen)1);
	}
    } else if (valsv) {

/*        Find singular values in a half-open interval. We aim */
/*        at -s (see leading comments) and we swap VL and VU */
/*        (as VUTGK and VLTGK), changing their signs. */

	*(unsigned char *)rngvx = 'V';
	vltgk = -(*vu);
	vutgk = -(*vl);
	i__1 = (*n << 1) - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    work[idtgk + i__] = 0.;
/* L30: */
	}
	dcopy_(n, &d__[1], &c__1, &work[ietgk], &c__2);
	i__1 = *n - 1;
	dcopy_(&i__1, &e[1], &c__1, &work[ietgk + 1], &c__2);
	i__1 = *n << 1;
	dstevx_("N", "V", &i__1, &work[idtgk], &work[ietgk], &vltgk, &vutgk, &
		iltgk, &iltgk, &abstol, ns, &s[1], &z__[z_offset], ldz, &work[
		itemp], &iwork[iiwork], &iwork[iifail], info, (ftnlen)1, (
		ftnlen)1);
	if (*ns == 0) {
	    return 0;
	} else {
	    if (wantz) {
		i__1 = *n << 1;
		dlaset_("F", &i__1, ns, &c_b22, &c_b22, &z__[z_offset], ldz, (
			ftnlen)1);
	    }
	}
    } else if (indsv) {

/*        Find the IL-th through the IU-th singular values. We aim */
/*        at -s (see leading comments) and indices are mapped into */
/*        values, therefore mimicking DSTEBZ. */

	iltgk = *il;
	iutgk = *iu;
	*(unsigned char *)rngvx = 'V';
	i__1 = (*n << 1) - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    work[idtgk + i__] = 0.;
/* L40: */
	}
	dcopy_(n, &d__[1], &c__1, &work[ietgk], &c__2);
	i__1 = *n - 1;
	dcopy_(&i__1, &e[1], &c__1, &work[ietgk + 1], &c__2);
	i__1 = *n << 1;
	dstevx_("N", "I", &i__1, &work[idtgk], &work[ietgk], &vltgk, &vltgk, &
		iltgk, &iltgk, &abstol, ns, &s[1], &z__[z_offset], ldz, &work[
		itemp], &iwork[iiwork], &iwork[iifail], info, (ftnlen)1, (
		ftnlen)1);
	vltgk = s[1] - smax * 2. * ulp * *n;
	i__1 = (*n << 1) - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    work[idtgk + i__] = 0.;
/* L50: */
	}
	dcopy_(n, &d__[1], &c__1, &work[ietgk], &c__2);
	i__1 = *n - 1;
	dcopy_(&i__1, &e[1], &c__1, &work[ietgk + 1], &c__2);
	i__1 = *n << 1;
	dstevx_("N", "I", &i__1, &work[idtgk], &work[ietgk], &vutgk, &vutgk, &
		iutgk, &iutgk, &abstol, ns, &s[1], &z__[z_offset], ldz, &work[
		itemp], &iwork[iiwork], &iwork[iifail], info, (ftnlen)1, (
		ftnlen)1);
	vutgk = s[1] + smax * 2. * ulp * *n;
	vutgk = min(vutgk,0.);

/*        If VLTGK=VUTGK, DSTEVX returns an error message, */
/*        so if needed we change VUTGK slightly. */

	if (vltgk == vutgk) {
	    vltgk -= tol;
	}

	if (wantz) {
	    i__1 = *n << 1;
	    i__2 = *iu - *il + 1;
	    dlaset_("F", &i__1, &i__2, &c_b22, &c_b22, &z__[z_offset], ldz, (
		    ftnlen)1);
	}
    }

/*     Initialize variables and pointers for S, Z, and WORK. */

/*     NRU, NRV: number of rows in U and V for the active submatrix */
/*     IDBEG, ISBEG: offsets for the entries of D and S */
/*     IROWZ, ICOLZ: offsets for the rows and columns of Z */
/*     IROWU, IROWV: offsets for the rows of U and V */

    *ns = 0;
    nru = 0;
    nrv = 0;
    idbeg = 1;
    isbeg = 1;
    irowz = 1;
    icolz = 1;
    irowu = 2;
    irowv = 1;
    split = FALSE_;
    sveq0 = FALSE_;

/*     Form the tridiagonal TGK matrix. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s[i__] = 0.;
/* L60: */
    }
    work[ietgk + (*n << 1) - 1] = 0.;
    i__1 = (*n << 1) - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	work[idtgk + i__] = 0.;
/* L70: */
    }
    dcopy_(n, &d__[1], &c__1, &work[ietgk], &c__2);
    i__1 = *n - 1;
    dcopy_(&i__1, &e[1], &c__1, &work[ietgk + 1], &c__2);


/*     Check for splits in two levels, outer level */
/*     in E and inner level in D. */

    i__1 = *n << 1;
    for (ieptr = 2; ieptr <= i__1; ieptr += 2) {
	if (work[ietgk + ieptr - 1] == 0.) {

/*           Split in E (this piece of B is square) or bottom */
/*           of the (input bidiagonal) matrix. */

	    isplt = idbeg;
	    idend = ieptr - 1;
	    i__2 = idend;
	    for (idptr = idbeg; idptr <= i__2; idptr += 2) {
		if (work[ietgk + idptr - 1] == 0.) {

/*                 Split in D (rectangular submatrix). Set the number */
/*                 of rows in U and V (NRU and NRV) accordingly. */

		    if (idptr == idbeg) {

/*                    D=0 at the top. */

			sveq0 = TRUE_;
			if (idbeg == idend) {
			    nru = 1;
			    nrv = 1;
			}
		    } else if (idptr == idend) {

/*                    D=0 at the bottom. */

			sveq0 = TRUE_;
			nru = (idend - isplt) / 2 + 1;
			nrv = nru;
			if (isplt != idbeg) {
			    ++nru;
			}
		    } else {
			if (isplt == idbeg) {

/*                       Split: top rectangular submatrix. */

			    nru = (idptr - idbeg) / 2;
			    nrv = nru + 1;
			} else {

/*                       Split: middle square submatrix. */

			    nru = (idptr - isplt) / 2 + 1;
			    nrv = nru;
			}
		    }
		} else if (idptr == idend) {

/*                 Last entry of D in the active submatrix. */

		    if (isplt == idbeg) {

/*                    No split (trivial case). */

			nru = (idend - idbeg) / 2 + 1;
			nrv = nru;
		    } else {

/*                    Split: bottom rectangular submatrix. */

			nrv = (idend - isplt) / 2 + 1;
			nru = nrv + 1;
		    }
		}

		ntgk = nru + nrv;

		if (ntgk > 0) {

/*                 Compute eigenvalues/vectors of the active */
/*                 submatrix according to RANGE: */
/*                 if RANGE='A' (ALLSV) then RNGVX = 'I' */
/*                 if RANGE='V' (VALSV) then RNGVX = 'V' */
/*                 if RANGE='I' (INDSV) then RNGVX = 'V' */

		    iltgk = 1;
		    iutgk = ntgk / 2;
		    if (allsv || vutgk == 0.) {
			if (sveq0 || smin < eps || ntgk % 2 > 0) {
/*                        Special case: eigenvalue equal to zero or very */
/*                        small, additional eigenvector is needed. */
			    ++iutgk;
			}
		    }

/*                 Workspace needed by DSTEVX: */
/*                 WORK( ITEMP: ): 2*5*NTGK */
/*                 IWORK( 1: ): 2*6*NTGK */

		    dstevx_(jobz, rngvx, &ntgk, &work[idtgk + isplt - 1], &
			    work[ietgk + isplt - 1], &vltgk, &vutgk, &iltgk, &
			    iutgk, &abstol, &nsl, &s[isbeg], &z__[irowz + 
			    icolz * z_dim1], ldz, &work[itemp], &iwork[iiwork]
			    , &iwork[iifail], info, (ftnlen)1, (ftnlen)1);
		    if (*info != 0) {
/*                    Exit with the error code from DSTEVX. */
			return 0;
		    }

/*                 Compute EMIN = abs(max(S(ISBEG:ISBEG+NSL-1))) */
/*                 (replaces F90 MAXVAL intrinsic) */

		    emin = 0.;
		    i__3 = nsl - 1;
		    for (i__ = 0; i__ <= i__3; ++i__) {
			if ((d__1 = s[isbeg + i__], abs(d__1)) > emin) {
			    emin = (d__2 = s[isbeg + i__], abs(d__2));
			}
/* L80: */
		    }

		    if (nsl > 0 && wantz) {

/*                    Normalize u=Z([2,4,...],:) and v=Z([1,3,...],:), */
/*                    changing the sign of v as discussed in the leading */
/*                    comments. The norms of u and v may be (slightly) */
/*                    different from 1/sqrt(2) if the corresponding */
/*                    eigenvalues are very small or too close. We check */
/*                    those norms and, if needed, reorthogonalize the */
/*                    vectors. */

			if (nsl > 1 && vutgk == 0. && ntgk % 2 == 0 && emin ==
				 0. && ! split) {

/*                       D=0 at the top or bottom of the active submatrix: */
/*                       one eigenvalue is equal to zero; concatenate the */
/*                       eigenvectors corresponding to the two smallest */
/*                       eigenvalues. */

			    i__3 = ntgk - 1;
			    for (i__ = 0; i__ <= i__3; ++i__) {
				z__[irowz + i__ + (icolz + nsl - 2) * z_dim1] 
					+= z__[irowz + i__ + (icolz + nsl - 1)
					 * z_dim1];
				z__[irowz + i__ + (icolz + nsl - 1) * z_dim1] 
					= 0.;
/* L90: */
			    }
			}

/* Computing MIN */
			i__4 = nsl - 1, i__5 = nru - 1;
			i__3 = min(i__4,i__5);
			for (i__ = 0; i__ <= i__3; ++i__) {
			    nrmu = dnrm2_(&nru, &z__[irowu + (icolz + i__) * 
				    z_dim1], &c__2);
			    if (nrmu == 0.) {
				*info = (*n << 1) + 1;
				return 0;
			    }
			    d__1 = 1. / nrmu;
			    dscal_(&nru, &d__1, &z__[irowu + (icolz + i__) * 
				    z_dim1], &c__2);
			    if (nrmu != 1. && (d__1 = nrmu - ortol, abs(d__1))
				     * sqrt2 > 1.) {
				i__4 = i__ - 1;
				for (j = 0; j <= i__4; ++j) {
				    zjtji = -ddot_(&nru, &z__[irowu + (icolz 
					    + j) * z_dim1], &c__2, &z__[irowu 
					    + (icolz + i__) * z_dim1], &c__2);
				    daxpy_(&nru, &zjtji, &z__[irowu + (icolz 
					    + j) * z_dim1], &c__2, &z__[irowu 
					    + (icolz + i__) * z_dim1], &c__2);
/* L100: */
				}
				nrmu = dnrm2_(&nru, &z__[irowu + (icolz + i__)
					 * z_dim1], &c__2);
				d__1 = 1. / nrmu;
				dscal_(&nru, &d__1, &z__[irowu + (icolz + i__)
					 * z_dim1], &c__2);
			    }
/* L120: */
			}
/* Computing MIN */
			i__4 = nsl - 1, i__5 = nrv - 1;
			i__3 = min(i__4,i__5);
			for (i__ = 0; i__ <= i__3; ++i__) {
			    nrmv = dnrm2_(&nrv, &z__[irowv + (icolz + i__) * 
				    z_dim1], &c__2);
			    if (nrmv == 0.) {
				*info = (*n << 1) + 1;
				return 0;
			    }
			    d__1 = -1. / nrmv;
			    dscal_(&nrv, &d__1, &z__[irowv + (icolz + i__) * 
				    z_dim1], &c__2);
			    if (nrmv != 1. && (d__1 = nrmv - ortol, abs(d__1))
				     * sqrt2 > 1.) {
				i__4 = i__ - 1;
				for (j = 0; j <= i__4; ++j) {
				    zjtji = -ddot_(&nrv, &z__[irowv + (icolz 
					    + j) * z_dim1], &c__2, &z__[irowv 
					    + (icolz + i__) * z_dim1], &c__2);
				    daxpy_(&nru, &zjtji, &z__[irowv + (icolz 
					    + j) * z_dim1], &c__2, &z__[irowv 
					    + (icolz + i__) * z_dim1], &c__2);
/* L140: */
				}
				nrmv = dnrm2_(&nrv, &z__[irowv + (icolz + i__)
					 * z_dim1], &c__2);
				d__1 = 1. / nrmv;
				dscal_(&nrv, &d__1, &z__[irowv + (icolz + i__)
					 * z_dim1], &c__2);
			    }
/* L160: */
			}
			if (vutgk == 0. && idptr < idend && ntgk % 2 > 0) {

/*                       D=0 in the middle of the active submatrix (one */
/*                       eigenvalue is equal to zero): save the */
/*                       corresponding eigenvector for later use (when */
/*                       bottom of the active submatrix is reached). */

			    split = TRUE_;
			    i__3 = ntgk - 1;
			    for (i__ = 0; i__ <= i__3; ++i__) {
				z__[irowz + i__ + (*n + 1) * z_dim1] = z__[
					irowz + i__ + (*ns + nsl) * z_dim1];
				z__[irowz + i__ + (*ns + nsl) * z_dim1] = 0.;
/* L170: */
			    }
			}
		    }

		    nsl = min(nsl,nru);
		    sveq0 = FALSE_;

/*                 Absolute values of the eigenvalues of TGK. */

		    i__3 = nsl - 1;
		    for (i__ = 0; i__ <= i__3; ++i__) {
			s[isbeg + i__] = (d__1 = s[isbeg + i__], abs(d__1));
/* L180: */
		    }

/*                 Update pointers for TGK, S and Z. */

		    isbeg += nsl;
		    irowz += ntgk;
		    icolz += nsl;
		    irowu = irowz;
		    irowv = irowz + 1;
		    isplt = idptr + 1;
		    *ns += nsl;
		    nru = 0;
		    nrv = 0;
		}
		if (irowz < *n << 1 && wantz) {
		    i__3 = irowz - 1;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			z__[i__ + icolz * z_dim1] = 0.;
/* L190: */
		    }
		}
/* L400: */
	    }
	    if (split && wantz) {

/*              Bring back eigenvector corresponding */
/*              to eigenvalue equal to zero. */

		i__2 = idend - ntgk + 1;
		for (i__ = idbeg; i__ <= i__2; ++i__) {
		    z__[i__ + (isbeg - 1) * z_dim1] += z__[i__ + (*n + 1) * 
			    z_dim1];
		    z__[i__ + (*n + 1) * z_dim1] = 0.;
/* L410: */
		}
	    }
	    --irowv;
	    ++irowu;
	    idbeg = ieptr + 1;
	    sveq0 = FALSE_;
	    split = FALSE_;
	}
/* L500: */
    }

/*     Sort the singular values into decreasing order (insertion sort on */
/*     singular values, but only one transposition per singular vector) */

    i__1 = *ns - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = 1;
	smin = s[1];
	i__2 = *ns + 1 - i__;
	for (j = 2; j <= i__2; ++j) {
	    if (s[j] <= smin) {
		k = j;
		smin = s[j];
	    }
/* L510: */
	}
	if (k != *ns + 1 - i__) {
	    s[k] = s[*ns + 1 - i__];
	    s[*ns + 1 - i__] = smin;
	    if (wantz) {
		i__2 = *n << 1;
		dswap_(&i__2, &z__[k * z_dim1 + 1], &c__1, &z__[(*ns + 1 - 
			i__) * z_dim1 + 1], &c__1);
	    }
	}
/* L520: */
    }

/*     If RANGE=I, check for singular values/vectors to be discarded. */

    if (indsv) {
	k = *iu - *il + 1;
	if (k < *ns) {
	    i__1 = *ns;
	    for (i__ = k + 1; i__ <= i__1; ++i__) {
		s[i__] = 0.;
/* L530: */
	    }
	    if (wantz) {
		i__1 = *ns;
		for (j = k + 1; j <= i__1; ++j) {
		    i__2 = *n << 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			z__[i__ + j * z_dim1] = 0.;
/* L540: */
		    }
/* L550: */
		}
	    }
	    *ns = k;
	}
    }

/*     Reorder Z: U = Z( 1:N,1:NS ), V = Z( N+1:N*2,1:NS ). */
/*     If B is a lower diagonal, swap U and V. */

    if (wantz) {
	i__1 = *ns;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *n << 1;
	    dcopy_(&i__2, &z__[i__ * z_dim1 + 1], &c__1, &work[1], &c__1);
	    if (lower) {
		dcopy_(n, &work[2], &c__2, &z__[*n + 1 + i__ * z_dim1], &c__1)
			;
		dcopy_(n, &work[1], &c__2, &z__[i__ * z_dim1 + 1], &c__1);
	    } else {
		dcopy_(n, &work[2], &c__2, &z__[i__ * z_dim1 + 1], &c__1);
		dcopy_(n, &work[1], &c__2, &z__[*n + 1 + i__ * z_dim1], &c__1)
			;
	    }
/* L600: */
	}
    }

    return 0;

/*     End of DBDSVDX */

} /* dbdsvdx_ */

