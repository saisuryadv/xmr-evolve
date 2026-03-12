/* dlaxre_f77.f -- translated by f2c (version 20240504).
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

/* Global flag to control GK-form: 0 = auto (use GK when applicable), 1 = force disable */
int xmr_disable_gkform_ = 0;

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

/* Table of constant values */

static doublereal c_b16 = 1.;
static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b57 = 0.;

/* Subroutine */ int dlaxre_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *gl, doublereal *gu, doublereal *abserr, doublereal *
	gaptol, doublereal *repr, integer *repi, doublereal *tau, integer *
	ewl_ae__, doublereal *ewl_lu__, char *emode, doublereal *rwork, 
	integer *info, ftnlen emode_len)
{
    /* System generated locals */
    integer i__1;
    integer r__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    integer i__, j;
    doublereal dl, dp, du, vl, vu, fac, gma, mgl;
    extern /* Subroutine */ int dlaxre_initewldqds__(integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, doublereal *);
    doublereal eps, mgu;
    integer jxg;
    doublereal aux;
    integer try__;
    doublereal diag, eoff;
    integer negl;
    doublereal dmin__, prec, emax, dmax__;
    integer negu, ixdp;
    doublereal auxl;
    integer itmp, type__;
    doublereal auxu, rtmp;
    integer ixrp, iseed[4];
    doublereal signd;
    integer iinfo, ixf77a;
    doublereal sfmin;
    integer sqdim, jxngn, ixsqd, ixsqe, twist, ixsqw;
    extern /* Subroutine */ int dlasq1_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *), dlaxrl_update__(integer *, integer *,
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *);
    extern doublereal dlamch_(char *, ftnlen);
    logical disdef, risdef;
    extern /* Subroutine */ int dlarnv_(integer *, integer *, integer *, 
	    doublereal *);
    doublereal offset;
    integer jyomga;
    extern integer dlaxrn_(integer *, doublereal *, integer *, doublereal *
	    );
    extern /* Subroutine */ int dlaxrr_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    logical gotrep;
    doublereal vstart, pivbase;

/*     IMPLICIT NONE */





/*  Purpose */
/*  ======= */

/*    Build root representation and init ew-list. */

/*      Given a symmetric tridiagonal matrix T by its entries in (D,E), */
/*    we determine a representation REP of a matrix M with */
/*        (T+err)  -  TAU  =  M, */
/*    where TAU is a shift and err is a matrix with norm bounded by */
/*    ABSERR (used for matrices with a nearly constant diagonal, to get */
/*    them into GK form and then use blocks). */

/*    The specific nature of M is kept hidden, it may be positive or */
/*    negative definite, or some other kind of good representation (RRR). */

/*    It is assumed that the input tridiagonal T is (numerically) */
/*    irreducible and properly scaled, that is, one block produced */
/*    by DLAXRA. */

/*    Upon exit, REPR, REPI and E combined define the representation. */
/*    (E may be modified, since the representation data does not include */
/*    the offdiagonal entries itself) */

/*    Once the representation has been found, the ew-list is initialized. */
/*    How this is done is determined by EMODE: */
/*     EMODE = 'd' or 'D'    Init ews based on ews that were computed */
/*                           to full accuracy by dqds. */
/*                           NOTE: Complexity O(n^2) */

/*     EMODE = 'o' or 'O'    Use Gershgorin Discs for initial setup, in */
/*                           particular to obtain outer bounds. Maybe */
/*                           get some more info, but not exceeding linear */
/*                           complexity. */
/*                           NOTE: Complexity O(n) */

/*     EMODE = 'g' or 'G'    Use full information from the Gershgorin */
/*                           Discs. */
/*                           NOTE: Complexity O(n logn) */

/*     EMODE = 'b' or 'B'    Start with Gershgorin discs (g), then refine */
/*                           all to full precision using bisection. */
/*                           NOTE: Complexity O(n^2) */

/* !! At the moment only modes 'd' and 'o' are implemented. */

/*  Arguments */
/*  ========= */

/*  ====================================================================== */

/*     .. Declarations .. */

/*     INTERFACE */
/*     FUNCTION DLAXRN(N, REPR, REPI, TAU) */
/*     IMPLICIT NONE */
/*     INTEGER  ::  DLAXRN */
/*     INTEGER,          INTENT(IN)  ::  N */
/*     INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), TAU */
/*     END FUNCTION DLAXRN */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRR( N, K, TYPE, E, PIVBASE, REPR, REPI ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, K, TYPE */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVBASE */
/*     DOUBLE PRECISION, INTENT(IN)  ::  E(1:N-1) */

/*     INTEGER,          INTENT(INOUT)  ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  REPR(4*N+3) */
/*     END SUBROUTINE DLAXRR */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRE_INITEWLDQDS( */
/*    $             N, REPR, REPI, IL, IU, */
/*    $             QDVALS, WIL, WIU, EWL_AE, EWL_LU */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)   ::  N, IL, IU, WIL, WIU */
/*     INTEGER,          INTENT(IN)   ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)   ::  QDVALS(IL:IU) */
/*     DOUBLE PRECISION, INTENT(IN)   ::  REPR(4*N+3) */

/*     INTEGER,          INTENT(OUT)  ::  EWL_AE(2*IL-1:2*IU) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  EWL_LU(2*IL-1:2*IU) */

/*     END SUBROUTINE DLAXRE_INITEWLDQDS */
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

/*     .. Parameters .. */

/*     How many tries are allowed to find a root representation */
/*     This is just a stopgap for situations where something is very wrong, */
/*     normally 2 or at most 3 tries should always suffice. */
/*     By how much ULP to perturb the primary data of the root rep. */
/*     Set to 0 for no perturbation. */

/*     .. Constants .. */


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
    --rwork;
    --repi;
    --repr;
    --e;
    --d__;
    --ewl_ae__;
    --ewl_lu__;

    /* Function Body */
    eps = dlamch_("Epsilon", (ftnlen)7);
    prec = dlamch_("Precision", (ftnlen)9);
    sfmin = dlamch_("Safe Minimum", (ftnlen)12);
    jxg = 0;
    jyomga = 4;
    jxngn = (*n << 1) + 1;
/*     ------------- */
/*      Set PIVBASE */
/*     ------------- */
    emax = e[1];
    dmin__ = d__[1];
    dmax__ = d__[1];
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MIN */
	d__1 = dmin__, d__2 = d__[i__];
	dmin__ = min(d__1,d__2);
/* Computing MAX */
	d__1 = dmax__, d__2 = d__[i__];
	dmax__ = max(d__1,d__2);
/* Computing MAX */
	d__1 = emax, d__2 = e[i__];
	emax = max(d__1,d__2);
/* L90001: */
    }
/* L90002: */
    pivbase = sqrt(sfmin) * 4 * max(emax,1.);
    gotrep = FALSE_;
    if (xmr_disable_gkform_ == 0 && dmax__ - dmin__ <= *abserr * 2) {
/*        GK-type matrix support for bidiagonal SVD. */
/*        For TGK matrices: DMAX=DMIN=0, so condition is trivially */
/*        satisfied. TAU=0, TYPE=1 with alternating OMEGA. */
/*        ----------------------------------- */
/*         Matrices with a constant diagonal */
/*        ----------------------------------- */
/*        Using an absolute perturbation bounded by ABSERR gives a */
/*        constant diagonal. We can then shift it away to obtain a */
/*        nicely blockable matrix in GK-form, possibly with an odd */
/*        dimension. */
/*        A bidiagonal matrix can have tiny singular values even */
/*        if the entries itself are not small, even below PIVBASE. */
/*        This is a problem since the shifting algorithms (DLAXR[N,T,S]) */
/*        require a shift > PIVBASE/EPS, otherwise they cannot guarantee */
/*        relative accuracy. Right now we discard the option to interpret */
/*        the root as GK-matrix if that happens. */
/* Computing MAX */
	d__1 = abs(*gl), d__2 = abs(*gu), d__1 = max(d__1,d__2), d__2 = *gu - 
		*gl;
/* Computing 2nd power */
	d__3 = eps;
	vl = max(d__1,d__2) * (d__3 * d__3);
	negl = 0;
	auxl = 0.;
	i__ = 1;
L90003:
	dl = d__[i__] - vl - auxl;
	if (abs(dl) < pivbase) {
	    dl = -pivbase;
	}
	if (dl < 0.) {
	    ++negl;
	}
	if (i__ == *n) {
	    goto L90004;
	}
/* Computing 2nd power */
	d__1 = e[i__];
	rtmp = d__1 * d__1;
	auxl = rtmp / dl;
	++i__;
	goto L90003;
L90004:
/*        Only proceed if there are no ews between 0 and VL */
	if (negl == (*n + 1) / 2) {
	    *tau = (dmin__ + dmax__) * .5;
	    if (e[1] <= e[*n - 1]) {
		itmp = 0;
		twist = *n;
	    } else {
		itmp = 1 - *n % 2;
		twist = 1;
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		repr[jxg + i__] = 0.;
		repi[jyomga + i__] = itmp;
		itmp = 1 - itmp;
/* L90005: */
	    }
/* L90006: */
	    type__ = 1;
	    gotrep = TRUE_;
	}
/*        Note */
/*        ==== */
/*          The above support fot this special case is rudimentary at */
/*          best. One should then really better handle the problem just */
/*          as a bSVD routine would to it: */
/*          (1) For tiny singular values, use a couple of sweeps of the */
/*              implicit zero-shift QR to get them away. */
/*          (2) Compute at most one half of the eigenpairs for GK, then */
/*              negate half of the vector entries to get the rest. */

    }
    if (! gotrep) {
/*        -------------------------------------- */
/*         Treat as standard Tridiagonal Matrix */
/*        -------------------------------------- */
/*        Determine preferred side for shifting */
	vl = *gl + (*gu - *gl) / 4;
	vu = *gu - (*gu - *gl) / 4;
	negl = 0;
	negu = 0;
	auxl = 0.;
	auxu = 0.;
	i__ = 1;
L90007:
	dl = d__[i__] - vl - auxl;
	du = d__[i__] - vu - auxu;
	if (abs(dl) < pivbase) {
	    dl = -pivbase;
	}
	if (abs(du) < pivbase) {
	    du = -pivbase;
	}
	if (dl < 0.) {
	    ++negl;
	}
	if (du < 0.) {
	    ++negu;
	}
	if (i__ == *n) {
	    goto L90008;
	}
/* Computing 2nd power */
	d__1 = e[i__];
	rtmp = d__1 * d__1;
	auxl = rtmp / dl;
	auxu = rtmp / du;
	++i__;
	goto L90007;
L90008:
	if (negl >= *n - negu) {
	    vstart = *gl;
	    signd = 1.;
	} else {
	    vstart = *gu;
	    signd = -1.;
	}
/* Computing MAX */
	d__1 = abs(vstart), d__2 = *gu - *gl;
	offset = (*n << 1) * eps * max(d__1,d__2);
/*        Setup workspace for factorization data */
/*         RWORK(1:N)    = DP(1:N) */
/*         RWORK(N+1:2N) = RP(1:N) */
	ixdp = 0;
	ixrp = *n;
	gotrep = FALSE_;
	try__ = 0;
L90009:
	*tau = vstart - try__ * signd * offset;
/*           Top to Bottom */
	disdef = TRUE_;
	aux = 0.;
	i__ = 1;
L90011:
	dp = d__[i__] - *tau - aux;
	if (abs(dp) < pivbase) {
	    dp = signd * pivbase;
	}
	rwork[ixdp + i__] = dp;
	if (i__ == *n || ! disdef) {
	    goto L90012;
	}
	disdef = disdef && d_sign(&c_b16, &dp) == signd;
/* Computing 2nd power */
	d__1 = e[i__];
	aux = d__1 * d__1 / dp;
	++i__;
	goto L90011;
L90012:
	if (i__ == *n) {
/* Computing MAX */
	    d__2 = *n * eps * abs(*tau);
	    if ((d__1 = rwork[ixdp + i__], abs(d__1)) <= max(d__2,*abserr)) {
		rwork[ixdp + i__] = 0.;
	    }
	}
/*           Bottom to Top */
	risdef = TRUE_;
	aux = 0.;
	i__ = *n;
L90013:
	dp = d__[i__] - *tau - aux;
	if (abs(dp) < pivbase) {
	    dp = signd * pivbase;
	}
	rwork[ixrp + i__] = dp;
	if (i__ == 1 || ! risdef) {
	    goto L90014;
	}
	risdef = risdef && d_sign(&c_b16, &dp) == signd;
	--i__;
/* Computing 2nd power */
	d__1 = e[i__];
	aux = d__1 * d__1 / dp;
	goto L90013;
L90014:
	if (i__ == 1) {
/* Computing MAX */
	    d__2 = *n * eps * abs(*tau);
	    if ((d__1 = rwork[ixrp + 1], abs(d__1)) <= max(d__2,*abserr)) {
		rwork[ixrp + 1] = 0.;
	    }
	}
/*           We want that both directions are semi-definite, otherwise */
/*           the results may differ when called with the flipped */
/*           matrix. */
	if (disdef && risdef) {
/*              Take min gamma */
	    if ((d__1 = rwork[ixdp + 1], abs(d__1)) < (d__2 = rwork[ixrp + *n]
		    , abs(d__2))) {
		gma = rwork[ixrp + 1];
		twist = 1;
	    } else {
		gma = rwork[ixdp + *n];
		twist = *n;
	    }
/*              Set repdata if successful */
	    i__1 = twist - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		repr[jxg + i__] = rwork[ixdp + i__];
/* L90015: */
	    }
/* L90016: */
	    repr[jxg + twist] = gma;
	    i__1 = *n;
	    for (i__ = twist + 1; i__ <= i__1; ++i__) {
		repr[jxg + i__] = rwork[ixrp + i__];
/* L90017: */
	    }
/* L90018: */
	    i__1 = jyomga + *n;
	    for (ixf77a = jyomga + 1; ixf77a <= i__1; ++ixf77a) {
		repi[ixf77a] = 0;
/* L90019: */
	    }
/* L90020: */
	    type__ = 0;
	    gotrep = TRUE_;
	}
	++try__;
	if (try__ == 10 || gotrep) {
	    goto L90010;
	}
	goto L90009;
L90010:
	;
    }
    if (! gotrep) {
	*info = 1;
	return 0;
    }
/*     -------------- */
/*      Perturb Root */
/*     -------------- */
/*     RWORK is not in use right now, so we can employ it freely. */
    if (TRUE_) {
	iseed[0] = 377;
	iseed[1] = 610;
	iseed[2] = 987;
	iseed[3] = 1597;
	fac = prec * 8;
	i__1 = (*n << 1) - 1;
	dlarnv_(&c__2, iseed, &i__1, &rwork[1]);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    e[i__] *= rwork[i__] * fac + 1.;
/* L90021: */
	}
/* L90022: */
	i__ = 1;
	j = *n;
L90023:
	if (i__ > *n) {
	    goto L90024;
	}
	repr[jxg + i__] *= rwork[j] * fac + 1.;
	++i__;
	++j;
	goto L90023;
L90024:
	;
    }
    dlaxrr_(n, &twist, &type__, &e[1], &pivbase, &repr[1], &repi[1]);
/*     =================================================================== */
/*     =                        Init EW-List                             = */
/*     =================================================================== */
/*     RWORK is not in use */
    twist = repi[2];
    type__ = repi[1];
    if (*(unsigned char *)emode == 'd' || *(unsigned char *)emode == 'D') {
/*        Use DQDS */
	ixsqd = 1;
	ixsqe = *n + 1;
	ixsqw = ixsqe + *n;
	if (type__ == 1) {
/*           We can ignore if the matrix was flipped, just take the */
/*           bidiagonal for the GK-form from top to bottom. */
	    sqdim = (*n + 1) / 2;
	    i__ = 1;
	    j = 1;
L90025:
	    rwork[ixsqd - 1 + i__] = e[j];
	    ++j;
	    if (i__ == sqdim) {
		goto L90026;
	    }
	    rwork[ixsqe - 1 + i__] = e[j];
	    ++j;
	    ++i__;
	    goto L90025;
L90026:
	    dlasq1_(&sqdim, &rwork[ixsqd], &rwork[ixsqe], &rwork[ixsqw], &
		    iinfo);
	    if (iinfo != 0) {
		*info = 2;
		return 0;
	    }
/*           The singular values of the bidiagonal factor are in SQD, */
/*           but in descending order. We have to mirror them to obtain */
/*           eigenvalues of M. */
	    i__ = 1;
	    j = *n;
L90027:
	    if (i__ >= j) {
		goto L90028;
	    }
	    rwork[ixsqd - 1 + j] = rwork[ixsqd - 1 + i__];
	    rwork[ixsqd - 1 + i__] = -rwork[ixsqd - 1 + i__];
	    --j;
	    ++i__;
	    goto L90027;
L90028:
	    ;
	} else {
	    if (twist == 1) {
		signd = d_sign(&c_b16, &repr[jxg + *n]);
		i__ = *n;
L90029:
		rwork[ixsqd - 1 + i__] = sqrt(signd * repr[jxg + i__]);
		if (i__ == 1) {
		    goto L90030;
		}
		rwork[ixsqe - 1 + i__ - 1] = sqrt(signd * repr[jxngn + i__]);
		--i__;
		goto L90029;
L90030:
		;
	    } else {
		signd = d_sign(&c_b16, &repr[jxg + 1]);
		i__ = 1;
L90031:
		rwork[ixsqd - 1 + i__] = sqrt(signd * repr[jxg + i__]);
		if (i__ == *n) {
		    goto L90032;
		}
		rwork[ixsqe - 1 + i__] = sqrt(signd * repr[jxngn + i__]);
		++i__;
		goto L90031;
L90032:
		;
	    }
	    dlasq1_(n, &rwork[ixsqd], &rwork[ixsqe], &rwork[ixsqw], &iinfo);
	    if (iinfo != 0) {
		*info = 2;
		return 0;
	    }
/*           Now SQD holds the singular values of the bidiagonal factor */
/*           in descending order. We have to square them, adjust to a */
/*           negative definite M (SIGND=-1) and invert their order. */
	    if (signd == 1.) {
/*              swap order */
		i__ = 1;
		j = *n;
L90033:
		if (i__ >= j) {
		    goto L90034;
		}
		rtmp = rwork[ixsqd - 1 + i__];
		rwork[ixsqd - 1 + i__] = rwork[ixsqd - 1 + j];
		rwork[ixsqd - 1 + j] = rtmp;
		++i__;
		--j;
		goto L90033;
L90034:
		;
	    }
/*           square and sign */
	    i__1 = *n - 1;
	    for (i__ = 0; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		d__1 = rwork[ixsqd + i__];
		rwork[ixsqd + i__] = signd * (d__1 * d__1);
/* L90035: */
	    }
/* L90036: */
	}


	dlaxre_initewldqds__(n, &repr[1], &repi[1], &c__1, n, &rwork[ixsqd], &
		c__1, n, &ewl_ae__[1], &ewl_lu__[1]);
    } else if (*(unsigned char *)emode == 'o' || *(unsigned char *)emode ==
	    'O') {
/*        Use Gersgorin Discs for outer bounds only. */
/*        Recompute them for the RRR, instead of shifting GL,GU from */
/*        the original T. */
	eoff = e[1];
	diag = repr[jxg + 1];
	if (1 == twist) {
	    diag += repr[jxngn + 2];
	}
	mgl = diag - eoff;
	mgu = diag + eoff;
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    eoff = e[i__ - 1] + e[i__];
	    diag = repr[jxg + i__];
	    if (i__ <= twist) {
		diag += repr[jxngn + i__ - 1];
	    }
	    if (i__ >= twist) {
		diag += repr[jxngn + i__ + 1];
	    }
/* Computing MIN */
	    d__1 = mgl, d__2 = diag - eoff;
	    mgl = min(d__1,d__2);
/* Computing MAX */
	    d__1 = mgu, d__2 = diag + eoff;
	    mgu = max(d__1,d__2);
/* L90037: */
	}
/* L90038: */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ewl_ae__[(i__ << 1) - 1] = 1;
	    ewl_ae__[i__ * 2] = *n;
	    ewl_lu__[(i__ << 1) - 1] = mgl;
	    ewl_lu__[i__ * 2] = mgu;
/* L90039: */
	}
/* L90040: */
    } else {
	*info = -12;
	return 0;
    }
/*     Always include information from the zero inertia in the ewlist */
    vl = 1.;
    vu = 1.;
    r__1 = dlaxrn_(n, &repr[1], &repi[1], &c_b57);
    dlaxrl_update__(&c__1, n, &vl, &vu, &ewl_ae__[1], &ewl_lu__[1], &c_b57, &
	    r__1);
    ++xmrstats_1.xnbis_init__;
    *info = 0;
    return 0;
} /* dlaxre_ */

