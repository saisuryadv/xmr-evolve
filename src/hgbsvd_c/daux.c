/* daux.f -- translated by f2c (version 20240504).
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

static doublereal c_b4 = 0.;
static doublereal c_b5 = 1.;
static doublereal c_b9 = -1.;
static integer c__1 = 1;
static integer c__9 = 9;
static integer c__3 = 3;
static integer c__5 = 5;

/* #define VERBOSE */
/* #define INFO_VERBOSE */
/* #define DT */
/* #define VS */
/* #define PR */
/* #define HL */
/* #define EVAL */
/* #define HLVAL */
/* #define HLVEC */
/* #define BLAS_BASED */
/* #define NaN_Check */
/* #define Assign_Check */

/*       ftncheck:                               okay */
/*       1111:                                   okay */
/*       invalid input token:                    okay */
/*       Praeprozessor:                          BLAS_BASED,NaN_Check,Assign_Check */
/*       Kommentare:                             okay */
/*       Semantik:                               okay */
/*       Liste subroutinen, ueberfluesige:       okay */
/*       Auskommentierte Anweisungen:            okay */
/*       NaN-Proof:                              okay */
/*       workarounds:                            okay */
/*       optimiert:                              okay */
/*       info:                                   okay */
/*       init:                                   okay */
/*       Dokumentation purpose / parameter */

/* ****************************************************************************** */

/*       daux.f : auxiliary-routines for dtest.f */

/*       This file contains a set of useful tools for the driver dtest : */

/*       Date:   Sep 10, 1999 */
/*       Author: Benedikt Grosser */

/*       dcheck_by_max: ( show ) */
/*       dcheck_by_norm:         ( nicht benutzt ) */
/*       dminiorthcheck:         ( nicht benutzt ) */
/*       dmview:                 ( nicht benutzt ) */
/*       dvview:                 ( nicht benutzt ) */
/*       ivview:                 ( nicht benutzt ) */
/*       dlv:                    ( nicht benutzt ) */
/*       dlv1:                   ( nicht benutzt ) */
/*       dlv2:                   ( nicht benutzt ) */
/*       dlv2rd:                 ( nicht benutzt ) */

/* ****************************************************************************** */

/* Subroutine */ int dcheck_by_max__(integer *m, integer *n, doublereal *a, 
	doublereal *sigma, integer *ddim, doublereal *b, integer *edim, 
	doublereal *u, integer *mmax, doublereal *v, integer *nmax, logical *
	updu, logical *updv, logical *transu, logical *transv, doublereal *
	res, doublereal *orthu, doublereal *orthv, doublereal *maxu, 
	doublereal *maxv, doublereal *maxr, doublereal *sigmx, doublereal *
	maxuc, doublereal *maxvc, doublereal *maxrc, doublereal *work, 
	integer *lwork, doublereal *work1, integer *lwork1, integer *info)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    integer i__, j, k;
    doublereal dn, tmp, tmp1, fact;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    integer wind;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    doublereal dummy;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    doublereal macheps;



/*       purpose */

/*       Check accuracy: */
/*       Deviation from orthogonality and residual. */

/*       orthu = norm(UT*U-I,1)  maxu = max(max(UT*U-I)) maxuc = maxu/(m*eps) */
/*       orthv = norm(VT*V-I,1)  maxu = max(max(VT*V-I)) maxvc = maxv/(n*eps) */
/*       res = norm(U*Sigma*VT-B,1) or norm(UT*B*V-Sigma,1) */
/*       maxr = max(max(U*Sigma*VT-B)) or max(max(UT*B*V-Sigma)) */
/*       maxrc = maxc/(n*eps*sigmx) */

/*       .. */
/*       .. parameters .. */
/*       .. */

/*       .. */
/*       .. external functions .. */
/*       .. */
/*       .. */
/*       .. local variables .. */
/*       .. */
/*       .. */
/*       .. executable statements .. */
/*       .. */
    /* Parameter adjustments */
    --sigma;
    --a;
    --b;
    u_dim1 = *mmax;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *nmax;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --work;
    --work1;

    /* Function Body */
    *info = 0;
    dummy = 0.;
    macheps = dlamch_("P", (ftnlen)1);
    dn = (doublereal) (*n);
    fact = dn * macheps;

/*       dummy statments to make ftnchek happy */

    tmp1 = 0.;
    if (*transu) {
	tmp1 = 1.;
    }
    if (*transv) {
	tmp1 = 1.;
    }
    dummy = tmp1;

/*       determine orthu */

    if (! (*updu)) {
	*orthu = -1.;
    } else {
	dlaset_("F", m, m, &c_b4, &c_b5, &work[1], m, (ftnlen)1);
	dgemm_("T", "N", m, m, m, &c_b5, &u[u_offset], mmax, &u[u_offset], 
		mmax, &c_b9, &work[1], m, (ftnlen)1, (ftnlen)1);
	*orthu = dlange_("1", m, m, &work[1], m, &dummy, (ftnlen)1);
	*maxu = 0.;
	*maxuc = 0.;
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *m;
	    for (j = 1; j <= i__2; ++j) {
		wind = (i__ - 1) * *m + j;
		tmp = (d__1 = work[wind], abs(d__1));
		if (tmp > *maxu) {
		    *maxu = tmp;
		}
		if (tmp > fact) {
		    work[wind] = (d__1 = work[wind], abs(d__1)) / fact;
/* Computing MAX */
		    d__2 = *maxu, d__3 = (d__1 = work[wind], abs(d__1));
		    *maxuc = max(d__2,d__3);
		} else {
		    work[wind] = 0.;
		}
	    }
	}
    }

/*       determine orthv */

    if (! (*updv)) {
	*orthv = -1.;
    } else {
	dlaset_("F", n, n, &c_b4, &c_b5, &work[1], n, (ftnlen)1);
	dgemm_("T", "N", n, n, n, &c_b5, &v[v_offset], nmax, &v[v_offset], 
		nmax, &c_b9, &work[1], n, (ftnlen)1, (ftnlen)1);
	*orthv = dlange_("1", n, n, &work[1], n, &dummy, (ftnlen)1);
	*maxv = 0.;
	*maxvc = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		wind = (i__ - 1) * *n + j;
		tmp = (d__1 = work[wind], abs(d__1));
		if (tmp > *maxv) {
		    *maxv = tmp;
		}
		if (tmp > fact) {
		    work[wind] = tmp / fact;
/* Computing MAX */
		    d__2 = *maxv, d__3 = (d__1 = work[wind], abs(d__1));
		    *maxvc = max(d__2,d__3);
		} else {
		    work[wind] = 0.;
		}
	    }
	}
    }

/*       determine res */

    if (! (*updu && *updv)) {
	*res = -1.;
    } else {
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {

/*           compute B*V(:,i) */

		wind = (i__ - 1) * *n;
		i__3 = *n - 1;
		for (j = 1; j <= i__3; ++j) {
		    work[j + wind] = a[j] * v[j + i__ * v_dim1] + b[j] * v[j 
			    + 1 + i__ * v_dim1];
		}
		work[*n + wind] = a[*n] * v[*n + i__ * v_dim1];
		work1[k + wind] = ddot_(n, &u[k * u_dim1 + 1], &c__1, &work[
			wind + 1], &c__1);
	    }
	}
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    work1[k + (k - 1) * *n] -= sigma[k];
	}
	*res = dlange_("1", m, n, &work1[1], m, &dummy, (ftnlen)1);
	*maxr = 0.;
	*maxrc = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *m;
	    for (j = 1; j <= i__2; ++j) {
		wind = (i__ - 1) * *m + j;
		if ((d__1 = work1[wind], abs(d__1)) > *maxr) {
		    *maxr = (d__1 = work1[wind], abs(d__1));
		}
		if ((d__1 = work1[wind], abs(d__1)) > fact * *sigmx) {
		    work1[wind] = (d__1 = work1[wind], abs(d__1)) / (fact * *
			    sigmx);
/* Computing MAX */
		    d__2 = *maxr, d__3 = (d__1 = work1[wind], abs(d__1));
		    *maxrc = max(d__2,d__3);
		    work1[wind] = (d__1 = work1[wind], abs(d__1)) * (fact * *
			    sigmx);
		} else {
		    work1[wind] = 0.;
		}
	    }
	}
    }

    return 0;
} /* dcheck_by_max__ */


/* ****************************************************************************** */

/* Subroutine */ int dcheck_results__(integer *m, integer *n, doublereal *a, 
	doublereal *sigma, integer *ddim, doublereal *b, integer *edim, 
	doublereal *u, integer *mmax, doublereal *v, integer *nmax, logical *
	updu, logical *updv, logical *transu, logical *transv, doublereal *
	res, doublereal *orthu, doublereal *orthv, doublereal *maxu, 
	doublereal *maxv, doublereal *maxr, doublereal *sigmx, doublereal *
	maxuc, doublereal *maxvc, doublereal *maxrc, doublereal *work, 
	integer *lwork, doublereal *work1, integer *lwork1, integer *info)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    integer i__, j, k;
    doublereal dn, tmp, tmp1, fact;
    integer wind;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    doublereal dummy;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    doublereal macheps;

    /* Fortran I/O blocks */
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };




/*       purpose */

/*       Check accuracy: */
/*       Deviation from orthogonality and residual. */

/*       orthu = norm(UT*U-I,1)  maxu = max(max(UT*U-I)) maxuc = maxu/(m*eps) */
/*       orthv = norm(VT*V-I,1)  maxu = max(max(VT*V-I)) maxvc = maxv/(n*eps) */
/*       res = norm(U*Sigma*VT-B,1) or norm(UT*B*V-Sigma,1) */
/*       maxr = max(max(U*Sigma*VT-B)) or max(max(UT*B*V-Sigma)) */
/*       maxrc = maxc/(n*eps*sigmx) */

/*       .. */
/*       .. parameters .. */
/*       .. */

/*       .. */
/*       .. external functions .. */
/*       .. */
/*       .. */
/*       .. local variables .. */
/*       .. */
/*       .. */
/*       .. executable statements .. */
/*       .. */
    /* Parameter adjustments */
    --sigma;
    --a;
    --b;
    u_dim1 = *mmax;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *nmax;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --work;
    --work1;

    /* Function Body */
    *info = 0;
    dummy = 0.;
    macheps = dlamch_("P", (ftnlen)1);
    dn = (doublereal) (*n);
    fact = dn * macheps;

/*       dummy statments to make ftnchek happy */

    tmp1 = 0.;
    if (*transu) {
	tmp1 = 1.;
    }
    if (*transv) {
	tmp1 = 1.;
    }
    dummy = tmp1;

/*       determine orthu */

    if (! (*updu)) {
	*orthu = -1.;
    } else {
	dlaset_("F", m, m, &c_b4, &c_b5, &work[1], m, (ftnlen)1);
	dgemm_("T", "N", m, m, m, &c_b5, &u[u_offset], mmax, &u[u_offset], 
		mmax, &c_b9, &work[1], m, (ftnlen)1, (ftnlen)1);
	*orthu = dlange_("1", m, m, &work[1], m, &dummy, (ftnlen)1);
	*maxu = 0.;
	*maxuc = 0.;
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *m;
	    for (j = 1; j <= i__2; ++j) {
		wind = (i__ - 1) * *m + j;
		tmp = (d__1 = work[wind], abs(d__1));
		if (tmp > *maxu) {
		    *maxu = tmp;
		}
		if (tmp > fact) {
		    s_wsle(&io___20);
		    do_lio(&c__9, &c__1, "orthu ", (ftnlen)6);
		    do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer)
			    );
		    do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(integer));
		    do_lio(&c__5, &c__1, (char *)&tmp, (ftnlen)sizeof(
			    doublereal));
		    e_wsle();
		    work[wind] = (d__1 = work[wind], abs(d__1)) / fact;
/* Computing MAX */
		    d__2 = *maxu, d__3 = (d__1 = work[wind], abs(d__1));
		    *maxuc = max(d__2,d__3);
		} else {
		    work[wind] = 0.;
		}
	    }
	}
    }

/*       determine orthv */

    if (! (*updv)) {
	*orthv = -1.;
    } else {
	dlaset_("F", n, n, &c_b4, &c_b5, &work[1], n, (ftnlen)1);
	dgemm_("T", "N", n, n, n, &c_b5, &v[v_offset], nmax, &v[v_offset], 
		nmax, &c_b9, &work[1], n, (ftnlen)1, (ftnlen)1);
	*orthv = dlange_("1", n, n, &work[1], n, &dummy, (ftnlen)1);
	*maxv = 0.;
	*maxvc = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		wind = (i__ - 1) * *n + j;
		tmp = (d__1 = work[wind], abs(d__1));
		if (tmp > *maxv) {
		    *maxv = tmp;
		}
		if (tmp > fact) {
/* 	write(*,*) 'orthv ',i,j,tmp */
		    work[wind] = tmp / fact;
/* Computing MAX */
		    d__2 = *maxv, d__3 = (d__1 = work[wind], abs(d__1));
		    *maxvc = max(d__2,d__3);
		} else {
		    work[wind] = 0.;
		}
	    }
	}
    }

/*       determine res */

    if (! (*updu && *updv)) {
	*res = -1.;
    } else {
	*res = 0.;
	*maxr = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           compute B*V(:,i) */

	    i__2 = *n - 1;
	    for (j = 1; j <= i__2; ++j) {
		work[j] = a[j] * v[j + i__ * v_dim1] + b[j] * v[j + 1 + i__ * 
			v_dim1];
	    }
	    work[*n] = a[*n] * v[*n + i__ * v_dim1];
/* 	do k = 1,n */
/* 	  write(*,*) k,work(k),u(k,i) */
/* 	enddo */
/* 	read(*,*) */
	    if (FALSE_) {
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
		    *res = 0.;
		    i__3 = *n;
		    for (k = 1; k <= i__3; ++k) {
			*res += work[k] * u[k + j * u_dim1];
		    }
		    if (i__ == j) {
			*res -= sigma[i__];
		    }
		    if (abs(*res) > fact * sigma[*n]) {
			s_wsle(&io___22);
			do_lio(&c__9, &c__1, "res ", (ftnlen)4);
			do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(
				integer));
			do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(
				integer));
			do_lio(&c__5, &c__1, (char *)&(*res), (ftnlen)sizeof(
				doublereal));
			do_lio(&c__5, &c__1, (char *)&sigma[i__], (ftnlen)
				sizeof(doublereal));
			e_wsle();
		    }
		}
	    } else {
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
		    work[j] -= sigma[i__] * u[j + i__ * u_dim1];
/* Computing MAX */
		    d__2 = *maxr, d__3 = (d__1 = work[j], abs(d__1));
		    *maxr = max(d__2,d__3);
/* 	if (ABS(WORK(J)).gt.fact) then */
/* 	  write(*,*) 'res ',i,j,ABS(WORK(J)) */
/* 	endif */
		}
	    }
	    work1[k] = dnrm2_(n, &work[1], &c__1);
/* Computing MAX */
	    d__2 = *res, d__3 = (d__1 = work1[k], abs(d__1));
	    *res = max(d__2,d__3);
/*       write(*,*) i,res,abs(work1(k)) */
	}
    }

    return 0;
} /* dcheck_results__ */


/* ****************************************************************************** */

/* Subroutine */ int dmview_(doublereal *a, integer *m, integer *n, integer *
	lld, integer *mb, integer *nb, char *name__, ftnlen name_len)
{
    /* Format strings */
    static char fmt_101[] = "(d13.4)";
    static char fmt_102[] = "(d13.4,d13.4)";
    static char fmt_103[] = "(d13.4,d13.4,d13.4)";
    static char fmt_104[] = "(d13.4,d13.4,d13.4,d13.4)";
    static char fmt_105[] = "(d13.4,d13.4,d13.4,d13.4,d13.4)";
    static char fmt_106[] = "(d13.4,d13.4,d13.4,d13.4,d13.4,d13.4)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsle(cilist *), e_wsle(void), do_lio(integer *, integer *, char 
	    *, ftnlen), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), 
	    e_wsfe(void);

    /* Local variables */
    integer i__, j, numcol;

    /* Fortran I/O blocks */
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };
    static cilist io___26 = { 0, 6, 0, 0, 0 };
    static cilist io___29 = { 0, 6, 0, 0, 0 };
    static cilist io___30 = { 0, 6, 0, 0, 0 };
    static cilist io___31 = { 0, 6, 0, 0, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_101, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_102, 0 };
    static cilist io___35 = { 0, 6, 0, fmt_103, 0 };
    static cilist io___36 = { 0, 6, 0, fmt_104, 0 };
    static cilist io___37 = { 0, 6, 0, fmt_105, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_106, 0 };
    static cilist io___39 = { 0, 6, 0, 0, 0 };
    static cilist io___40 = { 0, 6, 0, 0, 0 };
    static cilist io___41 = { 0, 6, 0, 0, 0 };
    static cilist io___42 = { 0, 6, 0, 0, 0 };






/* ****** */

    /* Parameter adjustments */
    --a;

    /* Function Body */
    s_wsle(&io___23);
    e_wsle();
    s_wsle(&io___24);
    do_lio(&c__9, &c__1, "+++++ Matrix ", (ftnlen)13);
    do_lio(&c__9, &c__1, name__, (ftnlen)5);
    do_lio(&c__9, &c__1, " ++++++++++++++++++++++++", (ftnlen)25);
    e_wsle();
    s_wsle(&io___25);
    e_wsle();
    s_wsle(&io___26);
    do_lio(&c__9, &c__1, "==================================", (ftnlen)34);
    e_wsle();
    i__1 = *n;
    i__2 = *nb;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	numcol = *nb;
	if (i__ + *nb > *n) {
	    numcol = *n - i__ + 1;
	}
	s_wsle(&io___29);
	e_wsle();
	s_wsle(&io___30);
	do_lio(&c__9, &c__1, " Spalte ", (ftnlen)8);
	do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " bis ", (ftnlen)5);
	i__3 = i__ + numcol - 1;
	do_lio(&c__3, &c__1, (char *)&i__3, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___31);
	do_lio(&c__9, &c__1, "==================================", (ftnlen)34)
		;
	e_wsle();
	i__3 = *m;
	for (j = 1; j <= i__3; ++j) {
	    if (numcol == 1) {
		s_wsfe(&io___33);
		do_fio(&c__1, (char *)&a[*lld * (i__ - 1) + j], (ftnlen)
			sizeof(doublereal));
		e_wsfe();
	    }
	    if (numcol == 2) {
		s_wsfe(&io___34);
		do_fio(&c__1, (char *)&a[*lld * (i__ - 1) + j], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&a[*lld * i__ + j], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	    if (numcol == 3) {
		s_wsfe(&io___35);
		do_fio(&c__1, (char *)&a[*lld * (i__ - 1) + j], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&a[*lld * i__ + j], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&a[*lld * (i__ + 1) + j], (ftnlen)
			sizeof(doublereal));
		e_wsfe();
	    }
	    if (numcol == 4) {
		s_wsfe(&io___36);
		do_fio(&c__1, (char *)&a[*lld * (i__ - 1) + j], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&a[*lld * i__ + j], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&a[*lld * (i__ + 1) + j], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&a[*lld * (i__ + 2) + j], (ftnlen)
			sizeof(doublereal));
		e_wsfe();
	    }
	    if (numcol == 5) {
		s_wsfe(&io___37);
		do_fio(&c__1, (char *)&a[*lld * (i__ - 1) + j], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&a[*lld * i__ + j], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&a[*lld * (i__ + 1) + j], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&a[*lld * (i__ + 2) + j], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&a[*lld * (i__ + 3) + j], (ftnlen)
			sizeof(doublereal));
		e_wsfe();
	    }
	    if (numcol == 6) {
		s_wsfe(&io___38);
		do_fio(&c__1, (char *)&a[*lld * (i__ - 1) + j], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&a[*lld * i__ + j], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&a[*lld * (i__ + 1) + j], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&a[*lld * (i__ + 2) + j], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&a[*lld * (i__ + 3) + j], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&a[*lld * (i__ + 4) + j], (ftnlen)
			sizeof(doublereal));
		e_wsfe();
	    }
	    if (j % *mb == 0) {
		s_wsle(&io___39);
		e_wsle();
	    }
	}
    }
    s_wsle(&io___40);
    e_wsle();
    s_wsle(&io___41);
    do_lio(&c__9, &c__1, "+++++ Ende ", (ftnlen)11);
    do_lio(&c__9, &c__1, name__, (ftnlen)5);
    do_lio(&c__9, &c__1, " +++++++++++++++++++++++++", (ftnlen)26);
    e_wsle();
    s_wsle(&io___42);
    e_wsle();
    return 0;
} /* dmview_ */

