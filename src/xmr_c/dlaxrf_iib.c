/* dlaxrf_iib_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxrf_iib__(integer *n, doublereal *repr, integer *repi,
	 doublereal *tau, integer *icbeg, integer *icend, integer *fewl_ae__, 
	doublereal *fewl_lu__, doublereal *slgap, doublereal *sugap, integer *
	sewl_ae__, doublereal *sewl_lu__)
{
    integer i__, j, xi;
    doublereal eps;
    integer xnn0, ison, json;
    extern /* Subroutine */ int dlaxrl_refine__(integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    doublereal arxfac[3];
    integer irxfac;
    doublereal fbound;
    logical gotbnd;
    extern integer dlaxrn_(integer *, doublereal *, integer *, doublereal *
	    );
    doublereal sbound;

/*     IMPLICIT NONE */



/*  Purpose */
/*  ======= */

/*    Use the father's inner bounds in FEWL to initialize bounds for */
/*    the son. */
/*    For each father bound of an interval, a couple of relative */
/*    relaxations are tried and the shift is applied to get a bound */
/*    for the son. */

/*  ====================================================================== */

/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRL_REFINE( */
/*    $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, I, J, LAMBDA, XI */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER         , INTENT(IN)  ::  IL, IU, I, J, XI */
/*     DOUBLE PRECISION, INTENT(IN)  ::  LAMBDA */

/*     INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU) */
/*     END SUBROUTINE DLAXRL_REFINE */
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


/*     .. Parameters .. */

/*         Intervals will be refined via bisection up to a relative */
/*         accuracy of BISACCFAC*N*EPS */
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
    --fewl_ae__;
    --fewl_lu__;
    --sewl_ae__;
    --sewl_lu__;

    /* Function Body */
    xnn0 = xmrstats_1.xnumfn;
    eps = dlamch_("Epsilon", (ftnlen)7);

    arxfac[0] = eps * 2;
    arxfac[1] = *n * eps * 4;
    arxfac[2] = *n * eps * 128;
/*     Loop over inner gaps */
    j = *icbeg - 1;
L90001:
    j = fewl_ae__[(j + 1) * 2];
    if (j == *icend) {
	goto L90002;
    }
    i__ = j + 1;

/*        Upper bound for J */

    gotbnd = FALSE_;
    irxfac = 1;
L90003:
    if (gotbnd || irxfac > 3) {
	goto L90004;
    }

    fbound = fewl_lu__[j * 2];
    fbound += arxfac[irxfac - 1] * abs(fbound);
    sbound = fbound - *tau;
    sbound += arxfac[irxfac - 1] * abs(sbound);
    ison = sewl_ae__[(j << 1) - 1];
    json = sewl_ae__[j * 2];

    if (sewl_lu__[(j << 1) - 1] < sbound && sbound < sewl_lu__[j * 2]) {
	xi = dlaxrn_(n, &repr[1], &repi[1], &sbound);
	dlaxrl_refine__(icbeg, icend, slgap, sugap, &sewl_ae__[1], &sewl_lu__[
		1], &ison, &json, &sbound, &xi);
    }

    gotbnd = sewl_lu__[j * 2] <= sbound;
    ++irxfac;
    goto L90003;
L90004:

/*        Lower bound for I */

    gotbnd = FALSE_;
    irxfac = 1;
L90005:
    if (gotbnd || irxfac > 3) {
	goto L90006;
    }

    fbound = fewl_lu__[(i__ << 1) - 1];
    fbound -= arxfac[irxfac - 1] * abs(fbound);
    sbound = fbound - *tau;
    sbound -= arxfac[irxfac - 1] * abs(sbound);
    ison = sewl_ae__[(i__ << 1) - 1];
    json = sewl_ae__[i__ * 2];

    if (sewl_lu__[(i__ << 1) - 1] < sbound && sbound < sewl_lu__[i__ * 2]) {
	xi = dlaxrn_(n, &repr[1], &repi[1], &sbound);
	dlaxrl_refine__(icbeg, icend, slgap, sugap, &sewl_ae__[1], &sewl_lu__[
		1], &ison, &json, &sbound, &xi);
    }

    gotbnd = sbound <= sewl_lu__[(i__ << 1) - 1];
    ++irxfac;
    goto L90005;
L90006:
    goto L90001;
L90002:
    xmrstats_1.xnbis_iib__ += xmrstats_1.xnumfn - xnn0;

    return 0;
} /* dlaxrf_iib__ */

