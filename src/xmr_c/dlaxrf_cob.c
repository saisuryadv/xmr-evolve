/* dlaxrf_cob_f77.f -- translated by f2c (version 20240504).
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
    doublereal xtime1, xtime2, xtime3, xddddd;
    integer xnblcks, xnnodes, xmaxdepth, xnumfn, xnumft, xnumgv, xnumgv0, 
	    xnumfs_2__, xnumfs_k__, xnbis_init__, xnbis_cob__, xnbis_iib__, 
	    xnbis_class__, xnbis_sng__, xnbis_clb__, xnbis_wasted__, xnrqi, 
	    xnrqibis, xmaxnrqi, xmaxnrqibis, xnenvgv, xnenvtf, xiiiii;
    logical xstealthmode;
} xmrstats_;

#define xmrstats_1 xmrstats_

/* Subroutine */ int dlaxrf_cob__(integer *n, doublereal *flgap, doublereal *
	flb, doublereal *fub, doublereal *fugap, doublereal *repr, integer *
	repi, doublereal *tau, doublereal *gaptol, integer *icbeg, integer *
	icend, integer *xizero, doublereal *slgap, doublereal *sugap, integer 
	*ewl_ae__, doublereal *ewl_lu__, integer *status)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    integer i__;
    doublereal l, freltight;
    logical ok;
    integer xi;
    doublereal lx, fbx, eps;
    integer isl, isu;
    doublereal sbx1, sbx2;
    integer ineg, ipos, type__;
    doublereal utol;
    extern /* Subroutine */ int dlaxrl_reset__(integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *);
    doublereal abnds[15];
    integer nbnds, ainrs[15];
    extern /* Subroutine */ int dlaxrl_update__(integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *);
    extern doublereal dlamch_(char *, ftnlen);
    doublereal fbound;
    logical gotbnd;
    doublereal frelax;
    extern integer dlaxrn_(integer *, doublereal *, integer *, doublereal *
	    );
    doublereal sbound, pivbase, wrelgap;

/*     IMPLICIT NONE */



/*  Purpose */
/*  ======= */

/*    Check Outer Bounds: Determine initial outer bounds, verify that */
/*    there is no eigenvalue underflow, and initialize a sons ewi-list. */

/*    Failure of some checks is indicated by STATUS != 0, the value then */
/*    indicates the kind of problem: */
/*        -1|+1  ew underflow on negative|positive side */
/*        -2|+2  two very small (about eps*tau) nonzero ews */
/*        -3|+3  could not verify left|right near bound */
/*        -5|+5  could not verify left|right far bound (the gaps) */

/*    If STATUS=0, all test are passed. Then the sons ew-list will be */
/*    initialized properly. */

/*    Independent of STATUS, XIZERO is always set to the sons zero */
/*    inertia. */

/*  ====================================================================== */

/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRL_RESET( */
/*    $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, L, XIL, U, XIU */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  IL, IU, XIL, XIU */
/*     DOUBLE PRECISION, INTENT(IN)  ::  L, U */

/*     INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU) */
/*     END SUBROUTINE DLAXRL_RESET */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRL_UPDATE( */
/*    $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, LAMBDA, XI */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  IL, IU, XI */
/*     DOUBLE PRECISION, INTENT(IN)  ::  LAMBDA */

/*     INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU) */
/*     END SUBROUTINE DLAXRL_UPDATE */
/*     END INTERFACE */
/*     INTERFACE */
/*     FUNCTION DLAXRN(N, REPR, REPI, TAU) */
/*     IMPLICIT NONE */
/*     INTEGER  ::  DLAXRN */
/*     INTEGER,          INTENT(IN)  ::  N */
/*     INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), TAU */
/*     END FUNCTION DLAXRN */
/*     END INTERFACE */

/*     .. Constants .. */

/*     The maximal number of bounds we generate to init the sons ew-list. */

/*     .. Parameters .. */

/*     Father and son bounds are relaxed before checking them by */
/*         2**LOGFBNDRELAXFAC * (N*EPS) */


/*     .. Locals .. */


/*           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*           !!!! SYNCHRONIZE ANY CHANGES HERE WITH xmr.h !!!!! */
/*           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*     Accumulative part: All these are accumulated over multiple */
/*     calls of stexr/laxrv. */
/*        Times for stages 1-3 (before/during/after block loop) */
/*        Marker to verify data integrity */

/*        Number of calls to dlaxrn, dlaxrt, dlaxrg, and */
/*        dlaxrs (all twists/single twist) */
/*        Number of bisection steps in various parts. */
/*        Count/max rqi steps */
/*        For computing envelopes, number of full twisted factos */
/*        and gv-calls */
/*        Marker to verify data integrity */
/*    Temporaries and Configuration */
/*        If true, calls to internal routines dlaxr[n,t,s] are not */
/*        counted. Higher level routines like DLAXRB are not affected. */
/*        We need this to hide bisections steps done within aux */
/*        routines like xmr_estrc. */

/*  ===== Executable Statements ========================================== */

    /* Parameter adjustments */
    --repi;
    --repr;
    --ewl_ae__;
    --ewl_lu__;

    /* Function Body */
    *status = 0;
    nbnds = 0;
    eps = dlamch_("Epsilon", (ftnlen)7);
    pivbase = repr[(*n << 2) + 3];
    utol = pivbase * (*n / eps);
    l = eps * abs(*tau);
    lx = l * 2;
    frelax = *n * eps * 128;
    wrelgap = *gaptol;
    freltight = *n * eps * 4;

/*     Determine the type of this shift based on the zero inertia. */

    ++nbnds;
    abnds[nbnds - 1] = 0.;
    ainrs[nbnds - 1] = *xizero;
/*     set IPOS as index of first positive ew */
/*     set INEG as index of first negative ew */
    ipos = (*xizero + 1) / 2 + 1;
    ineg = *xizero / 2;
    if (*xizero <= (*icbeg << 1) - 1) {
	type__ = 1;
    } else if (*xizero >= (*icend << 1) - 1) {
	type__ = -1;
    } else {
	type__ = 0;
    }

/*     ------------------------------- */
/*      Check positive very small ews */
/*     ------------------------------- */

    if (type__ != -1 && *status == 0) {
	xi = dlaxrn_(n, &repr[1], &repi[1], &l);
	++nbnds;
	abnds[nbnds - 1] = l;
	ainrs[nbnds - 1] = xi;
	if (! (xi <= (ipos << 1) - 1)) {

/*           L is not a lower bound for the first positive ew */

	    xi = dlaxrn_(n, &repr[1], &repi[1], &utol);
	    ++nbnds;
	    abnds[nbnds - 1] = utol;
	    ainrs[nbnds - 1] = xi;
	    if (! (xi == (ipos << 1) - 2 || xi == (ipos << 1) - 1)) {
/*              underflow on positive side */
		*status = 1;
	    } else if (ipos < *icend) {
		xi = dlaxrn_(n, &repr[1], &repi[1], &lx);
		++nbnds;
		abnds[nbnds - 1] = lx;
		ainrs[nbnds - 1] = xi;
		if (! (xi <= (ipos + 1 << 1) - 1)) {
/*                 we more than one very small ew between L and LX */
		    *status = 2;
		}
	    }
	}
    }

/*     ------------------------------- */
/*      Check negative very small ews */
/*     ------------------------------- */

    if (type__ != 1 && *status == 0) {
	d__1 = -l;
	xi = dlaxrn_(n, &repr[1], &repi[1], &d__1);
	++nbnds;
	abnds[nbnds - 1] = -l;
	ainrs[nbnds - 1] = xi;
	if (! (xi >= (ineg << 1) - 1)) {

/*          -L is not an upper bound for first negative ew */

	    d__1 = -utol;
	    xi = dlaxrn_(n, &repr[1], &repi[1], &d__1);
	    ++nbnds;
	    abnds[nbnds - 1] = -utol;
	    ainrs[nbnds - 1] = xi;
	    if (! (xi == (ineg << 1) - 1 || xi == ineg << 1)) {
/*              underflow on negative side */
		*status = -1;
	    } else if (ineg > *icbeg) {
		d__1 = -lx;
		xi = dlaxrn_(n, &repr[1], &repi[1], &d__1);
		++nbnds;
		abnds[nbnds - 1] = -lx;
		ainrs[nbnds - 1] = xi;
		if (! (xi >= (ineg - 1 << 1) - 1)) {
/*                 more than one very small ew between -LX and -L */
		    *status = -2;
		}
	    }
	}
    }

/*     ------------------ */
/*      Check left bound */
/*     ------------------ */

    if (*status == 0) {
	gotbnd = FALSE_;
	fbound = *flb - abs(*flb) * freltight;
	if (type__ == 1 && *tau >= fbound) {
	    gotbnd = TRUE_;
	} else {
	    sbound = fbound - *tau - (d__1 = fbound - *tau, abs(d__1)) * 
		    freltight;
	    xi = dlaxrn_(n, &repr[1], &repi[1], &sbound);
	    ++nbnds;
	    abnds[nbnds - 1] = sbound;
	    ainrs[nbnds - 1] = xi;
	    if (xi == (*icbeg << 1) - 2 || xi == (*icbeg << 1) - 1) {
		gotbnd = TRUE_;
	    }
	}
	fbound = *flb - abs(*flb) * frelax;
	if (! gotbnd && ! (type__ == 1 && *tau >= fbound)) {
	    sbound = fbound - *tau - (d__1 = fbound - *tau, abs(d__1)) * 
		    frelax;
	    xi = dlaxrn_(n, &repr[1], &repi[1], &sbound);
	    ++nbnds;
	    abnds[nbnds - 1] = sbound;
	    ainrs[nbnds - 1] = xi;
	    if (! (xi == (*icbeg << 1) - 2 || xi == (*icbeg << 1) - 1)) {
/*              Cannot verify left inner far bound. If this shift is */
/*              to the left we can live with that, since 0 is a bound. */
		if (type__ != 1) {
		    *status = -3;
		}
	    }
	}
/*         FBOUND = FLB - ABS(FLB)*FRELAX */

/*         IF( .NOT.(TYPE.EQ.+1 .AND. TAU.GE.FBOUND) )THEN */

/*            SBOUND = (FBOUND-TAU) - ABS(FBOUND-TAU)*FRELAX */
/*            XI = DLAXRN( N, REPR, REPI, SBOUND ) */
/*            NBNDS = NBNDS + 1 */
/*            ABNDS(NBNDS) = SBOUND */
/*            AINRS(NBNDS) = XI */

/*            IF( .NOT.(XI.EQ.2*ICBEG-2 .OR. XI.EQ.2*ICBEG-1) )THEN */
/* *              cannot verify left inner far bound */
/*               STATUS = -3 */
/*            ENDIF */
/*         ENDIF */
    }

/*     ------------------- */
/*      Check right bound */
/*     ------------------- */

    if (*status == 0) {
	gotbnd = FALSE_;
	fbound = *fub + abs(*fub) * freltight;
	if (type__ == -1 && *tau <= fbound) {
	    gotbnd = TRUE_;
	} else {
	    sbound = fbound - *tau + (d__1 = fbound - *tau, abs(d__1)) * 
		    freltight;
	    xi = dlaxrn_(n, &repr[1], &repi[1], &sbound);
	    ++nbnds;
	    abnds[nbnds - 1] = sbound;
	    ainrs[nbnds - 1] = xi;
	    if (xi == (*icend << 1) - 1 || xi == *icend << 1) {
		gotbnd = TRUE_;
	    }
	}
	fbound = *fub + abs(*fub) * frelax;
	if (! gotbnd && ! (type__ == -1 && *tau <= fbound)) {
	    sbound = fbound - *tau + (d__1 = fbound - *tau, abs(d__1)) * 
		    frelax;
	    xi = dlaxrn_(n, &repr[1], &repi[1], &sbound);
	    ++nbnds;
	    abnds[nbnds - 1] = sbound;
	    ainrs[nbnds - 1] = xi;
	    if (! (xi == (*icend << 1) - 1 || xi == *icend << 1)) {
/*              Cannot verifY right inner bound. If this shift is */
/*              to the right we can live with that, since 0 is a bound. */
		if (type__ != 1) {
		    *status = 3;
		}
	    }
	}
/*         FBOUND = FUB + ABS(FUB)*FRELAX */
/*         IF( .NOT.(TYPE.EQ.-1 .AND. TAU.LE.FBOUND) )THEN */

/*            SBOUND = (FBOUND-TAU) + ABS(FBOUND-TAU)*FRELAX */
/*            XI = DLAXRN( N, REPR, REPI, SBOUND ) */
/*            NBNDS = NBNDS + 1 */
/*            ABNDS(NBNDS) = SBOUND */
/*            AINRS(NBNDS) = XI */

/*            IF( .NOT.(XI.EQ.2*ICEND-1 .OR. XI.EQ.2*ICEND) )THEN */
/* *              cannot verifY right inner far bound */
/*               STATUS = +3 */
/*            ENDIF */
/*         ENDIF */
    }

/*     ----------------------- */
/*      Check left outer gap */
/*     ----------------------- */

    if (*status == 0 && type__ != 1) {
	fbound = *flb - abs(*flb) * frelax;
	fbx = *flb - *flgap + (d__1 = *flb - *flgap, abs(d__1)) * frelax;
	sbound = fbound - *tau - (d__1 = fbound - *tau, abs(d__1)) * frelax;
/* Computing MAX */
	d__2 = sbound * 2, d__3 = fbx - *tau + (d__1 = fbx - *tau, abs(d__1)) 
		* frelax;
	sbx1 = max(d__2,d__3);
	sbx2 = sbound - wrelgap * abs(sbound);
	ok = FALSE_;
/*        1st try: go wide */
	if (sbx1 < sbx2) {
	    xi = dlaxrn_(n, &repr[1], &repi[1], &sbx1);
	    ++nbnds;
	    abnds[nbnds - 1] = sbx1;
	    ainrs[nbnds - 1] = xi;
	    ok = xi == (*icbeg - 1 << 1) - 1 || xi == *icbeg - 1 << 1;
	}
/*        2nd try: go tight */
	if (! ok) {
	    xi = dlaxrn_(n, &repr[1], &repi[1], &sbx2);
	    ++nbnds;
	    abnds[nbnds - 1] = sbx2;
	    ainrs[nbnds - 1] = xi;
	    ok = xi == (*icbeg - 1 << 1) - 1 || xi == *icbeg - 1 << 1;
	}
	if (! ok) {
/*           Could not establish sons left outer gap */
	    *status = -5;
	}
    }

/*     ----------------------- */
/*      Check right outer gap */
/*     ----------------------- */

    if (*status == 0 && type__ != -1) {
	fbound = *fub + abs(*fub) * frelax;
	fbx = *fub + *fugap - (d__1 = *fub + *fugap, abs(d__1)) * frelax;
	sbound = fbound - *tau + (d__1 = fbound - *tau, abs(d__1)) * frelax;
/* Computing MIN */
	d__2 = sbound * 2, d__3 = fbx - *tau - (d__1 = fbx - *tau, abs(d__1)) 
		* frelax;
	sbx1 = min(d__2,d__3);
	sbx2 = sbound + wrelgap * abs(sbound);
	ok = FALSE_;
/*        1st try: go wide */
	if (sbx1 > sbx2) {
	    xi = dlaxrn_(n, &repr[1], &repi[1], &sbx1);
	    ++nbnds;
	    abnds[nbnds - 1] = sbx1;
	    ainrs[nbnds - 1] = xi;
	    ok = xi == *icend << 1 || xi == (*icend << 1) + 1;
	}
/*        2nd try: go tight */
	if (! ok) {
	    xi = dlaxrn_(n, &repr[1], &repi[1], &sbx2);
	    ++nbnds;
	    abnds[nbnds - 1] = sbx2;
	    ainrs[nbnds - 1] = xi;
	    ok = xi == *icend << 1 || xi == (*icend << 1) + 1;
	}
	if (! ok) {
/*           Could not establish sons right outer gap */
	    *status = 5;
	}
    }
/*     Don't count xizero */
    xmrstats_1.xnbis_cob__ = xmrstats_1.xnbis_cob__ + nbnds - 1;

/*     ------------------------ */
/*      Init the sons ewi-list */
/*     ------------------------ */

    if (*status == 0) {
	isl = 0;
	isu = 0;
	i__1 = nbnds;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xi = ainrs[i__ - 1];
	    if (xi <= (*icbeg << 1) - 1) {
		isl = i__;
	    }
	    if (xi >= (*icend << 1) - 1) {
		isu = i__;
	    }
/* L90001: */
	}
/* L90002: */
	dlaxrl_reset__(icbeg, icend, slgap, sugap, &ewl_ae__[1], &ewl_lu__[1],
		 &abnds[isl - 1], &ainrs[isl - 1], &abnds[isu - 1], &ainrs[
		isu - 1]);
	i__1 = nbnds;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (i__ != isl && i__ != isu) {
		dlaxrl_update__(icbeg, icend, slgap, sugap, &ewl_ae__[1], &
			ewl_lu__[1], &abnds[i__ - 1], &ainrs[i__ - 1]);
	    }
/* L90003: */
	}
/* L90004: */
/*        If we have a zero bound, the outer gaps were not checked since */
/*        we know there is a relative gap, the outer gaps for the son */
/*        will then be zero. To avoid zero gaps, we just copy the one */
/*        from the father. The gap are not used anyway, since dlaxrx has */
/*        special handling for zero eigenvalues. */
	if (ewl_lu__[(*icbeg << 1) - 1] == 0.) {
	    *slgap = *flgap;
	}
	if (ewl_lu__[*icend * 2] == 0.) {
	    *sugap = *fugap;
	}
    }
    return 0;
} /* dlaxrf_cob__ */

