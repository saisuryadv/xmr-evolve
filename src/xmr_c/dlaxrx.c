/* dlaxrx_f77.f -- translated by f2c (version 20240504).
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

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b9 = 1.;

/* Subroutine */ int wsreq_xrx__(integer *n, integer *reqr, integer *reqi)
{
/*     IMPLICIT NONE */
    *reqr = *n * 7;
    *reqi = *n + 1;
    return 0;
} /* wsreq_xrx__ */


/* *********************************************************************** */

/* Subroutine */ int dlaxrx_(integer *n, doublereal *e, doublereal *repr, 
	integer *repi, integer *index, doublereal *mingap, doublereal *lambda,
	 doublereal *z__, integer *isuppz, doublereal *rwork, integer *lrwork,
	 integer *iwork, integer *liwork)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    integer biscount, i__;
    doublereal lb, ub, dir, eps, tmp;
    integer nbis;
    doublereal zeta;
    integer ixdp, iter, nrqi, ixrp, ixgma;
    doublereal resid;
    integer ixf77a, ixbuf[2];
    doublereal rtmpa, rtmpb, normz;
    integer ixwrk, twist;
    doublereal bcklda;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal absgma;
    logical havelb;
    doublereal acctol, mingma;
    logical haveub;
    integer ntwbad, bufsel;
    extern /* Subroutine */ int dlaxrg_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    doublereal optlda;
    extern /* Subroutine */ int dlaxrt_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *);
    doublereal rqcorr, restol, cuttol, optres;
    integer status, iytwok;
    extern /* Subroutine */ int dlaxrg0_(integer *, doublereal *, doublereal *
	    , integer *, doublereal *, doublereal *);
    logical haveacc;
    doublereal xlambda;
    logical haveopt;

/*     IMPLICIT NONE */




/*  Purpose */
/*  ======= */

/*     Compute an eigenvector for a singleton via RQI, backed up by */
/*     (improved) bisection. */

/*     Note:  We do not use initial bounds for the eigenvalue, as */
/*            these would be based on normal bisection. Trusting them */
/*            can lose to much accuracy. */
/*            This means we implicitly assume LAMBDA to have some */
/*            initial accuracy, and in particular that it is placed far */
/*            enough from other ews so that the RQI can converge. */
/*            If that is not the case, the routine will detect this */
/*            and produce a runtime abort (in debug mode only, this */
/*            should become an error code). */

/*  ====================================================================== */


/*     .. Parameters .. */


/*    ACCTOLFAC */
/*        Accept eigenpair as soon as change in eigenvalue is less than */
/*        ACCTOLFAC*Eps*|lambda| */


/*    RESTOLFAC */
/*        Accept eigenpair as soon as resiudal falls below */
/*        RESTOL=RESTOLFAC*n*Eps*mingap. */


/*     CUTTOLFAC */
/*        Allowed change to the residual through cutting the support is */
/*        CUTTOLFAC*RESTOL. */


/*     PBBLEN ('bisection burst length') */
/*        Number of steps in a bisection burst, >= 1 */

/*     RESDECFAC ('residual (norm) decrease factor') */
/*        Factor in (0,1) by how much we would like the residual norm to */
/*        decrease per iteration at least. */


/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRG0( */
/*    $             N, E, REPR, REPI, CUTTOL, Z */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)   ::  N */
/*     INTEGER,          INTENT(IN)   ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)   ::  REPR(4*N+3) */
/*     DOUBLE PRECISION, INTENT(IN)   ::  E(N-1) */
/*     DOUBLE PRECISION, INTENT(IN)   ::  CUTTOL */

/*     DOUBLE PRECISION, INTENT(OUT)  ::  Z(N) */
/*     END SUBROUTINE DLAXRG0 */
/*     SUBROUTINE DLAXRG( */
/*    $             N, K, D, R, E, GAMMA, CUTTOL, */
/*    $             Z, NORMZ, RESID, RQCORR */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, K */
/*     DOUBLE PRECISION, INTENT(IN)  ::  GAMMA, CUTTOL */
/*     DOUBLE PRECISION, INTENT(IN)  ::  D(1:K-1), R(K+1:N), E(1:N-1) */

/*     DOUBLE PRECISION, INTENT(OUT)  ::  NORMZ, RESID, RQCORR */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  Z(N) */
/*     END SUBROUTINE DLAXRG */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRT( */
/*    $             N, J1, J2, REPR, REPI, */
/*    $             TAU, DPLUS, RPLUS, GAMMAP, TWISTOK, RWORK */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, J1, J2 */
/*     INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  TAU */
/*     DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3) */

/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK(2*N) */

/*     INTEGER,          INTENT(OUT)  ::  TWISTOK(J1:J2) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  DPLUS(1:N-1), RPLUS(2:N) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  GAMMAP(J1:J2) */
/*     END SUBROUTINE DLAXRT */
/*     END INTERFACE */

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

/*  ----- Executable Statements ----------------------------------------- */

    /* Parameter adjustments */
    --z__;
    --repi;
    --repr;
    --e;
    --isuppz;
    --rwork;
    --iwork;

    /* Function Body */
    eps = dlamch_("Epsilon", (ftnlen)7);
    ixgma = 1;
/*  GMA(1:N) */
    ixdp = ixgma + *n;
/*  DPLUS(1:N-1) */
    ixrp = ixdp + *n - 1;
/*  RPLUS(2:N) */
    ixbuf[0] = ixrp + *n - 1;
    ixbuf[1] = ixbuf[0] + *n;
    ixwrk = ixbuf[1] + *n;
/*     Further workspace:  2N for DLAXRT */
/*     We hold two buffers to temporarily store the computed vectors, */
/*     indicated by IXBUFA and IXBUFB. The index IXACTBUF points to the */
/*     one we are currently using. Once a vector was found (HAVEOPT= */
/*     .TRUE.) the other buffer will hold the so-far best one. */
    bufsel = 1;
    iytwok = 1;
    acctol = eps * 2.;
    restol = *n * eps * 2. * *mingap;
    cuttol = restol * .5;
    havelb = FALSE_;
    haveub = FALSE_;
    haveopt = FALSE_;
    optres = -1.;
/*     Just to soothe compilers, these settings will never be used. */
    lb = 0.;
    ub = 0.;
    optlda = 0.;
    biscount = 0;
    haveacc = FALSE_;
    iter = 0;
    nbis = 0;
    nrqi = 0;
/*     ======================== */
/*      Special case: Lambda=0 */
/*     ======================== */
    if (*lambda == 0.) {
	dlaxrg0_(n, &e[1], &repr[1], &repi[1], &cuttol, &z__[1]);
	optlda = *lambda;
	optres = 0.;
	goto L100;
    }
/*     =============== */
/*        Main Loop */
/*     =============== */
L90001:

/*        LAMBDA is set, either LB and UB are still undefined, or LAMBDA */
/*        is consistent with them. */

    dlaxrt_(n, &c__1, n, &repr[1], &repi[1], lambda, &rwork[ixdp], &rwork[
	    ixrp], &rwork[ixgma], &iwork[iytwok], &rwork[ixwrk]);
    ++iter;

/*        Find minimal gamma */

    ntwbad = 0;
    twist = 0;
    absgma = -1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iwork[iytwok - 1 + i__] != 0) {
	    tmp = (d__1 = rwork[ixgma - 1 + i__], abs(d__1));
	    if (twist == 0 || tmp < absgma) {
		absgma = tmp;
		twist = i__;
	    }
	} else {
	    ++ntwbad;
	}
/* L90003: */
    }
/* L90004: */
    mingma = rwork[ixgma - 1 + twist];

/*        Update guidance interval */

    if (mingma == 0.) {
	lb = *lambda;
	havelb = TRUE_;
	ub = *lambda;
	haveub = TRUE_;
    } else if (mingma < 0.) {
	ub = *lambda;
	haveub = TRUE_;
    } else {
	lb = *lambda;
	havelb = TRUE_;
    }
    if (havelb && haveub) {
	haveacc = haveacc || ub - lb <= acctol * abs(*lambda);
    }

/*        Cancel a bisection request if full accuracy is reached or if */
/*        we can see that this is a good one. */

    if (biscount > 0 && (haveacc || absgma <= restol * 5 || haveopt && absgma 
	    <= optres)) {
	biscount = 0;
    }
    if (biscount != 0) {

/*           ----------- */
/*            Bisection */
/*           ----------- */

	++nbis;
	*lambda = (lb + ub) * .5;
	--biscount;
    } else {

/*           --------------- */
/*            Full RQI step */
/*           --------------- */

	dlaxrg_(n, &twist, &rwork[ixdp], &rwork[ixrp - 2 + twist + 1], &e[1], 
		&mingma, &cuttol, &rwork[ixbuf[bufsel - 1]], &normz, &resid, &
		rqcorr);
	++nrqi;
/*           Note: The current optimum is needed to evaluate the quality */
/*           of this step. This is the reason why we record the optimum */
/*           only at the end of the loop body below, and the loop exit is */
/*           there as well. */
	bcklda = *lambda;
	xlambda = *lambda + rqcorr;
	status = 0;
	if (resid <= restol) {

/*              converged: residual tolerance is reached */

	    status = 1;
/* @CHECKS on */
/*               IF( INDEX .EQ. 150 )THEN */
/*                  CALL FABORT('','') */
/*               ENDIF */
/* @CHECKS ! */
	} else if (xlambda == *lambda || abs(rqcorr) <= acctol * abs(*lambda))
		 {

/*              converged: rqc tolerance reached */

	    status = 2;
	} else if (haveacc) {

/*              converged: guidance interval has full acc */
/*                (although this rqc might jump out of it) */

	    status = 3;
	}
/*           If STATUS=0 still, we are not yet converged. The remainder of */
/*           the loop body centers around setting the next iterate LAMBDA. */
/*           Once we have one, STATUS is set to sth < 0. */
	if (status == 0 && havelb && haveub) {

/*              Check if we stay in the guidance interval, if not */
/*              switch to bisection */

	    if (xlambda <= lb || xlambda >= ub) {
		biscount = 4;
		status = -1;
	    }
	}
	if (status == 0) {

/*              Check if the residual norm is sufficiently better than */
/*              the current optimum. If so, then continue normally. */

	    if (! haveopt || resid <= optres * .7) {
		*lambda = xlambda;
		status = -2;
	    }
	}
	if (status == 0) {

/*              The residual is not good enough, try alternatives to avoid */
/*              slow convergence. */

	    dir = d_sign(&c_b9, &rqcorr);
/* Computing MAX */
	    d__1 = abs(rqcorr) * 3;
	    zeta = *lambda + dir * max(d__1,resid);
	    if (! havelb || ! haveub) {
		*lambda = zeta;
	    } else {
		if (dir == -1.) {
		    rtmpa = (lb + ub) * .5;
		    rtmpb = ub;
		} else {
		    rtmpa = lb;
		    rtmpb = (lb + ub) * .5;
		}
		if (rtmpa <= zeta && zeta <= rtmpb) {
		    *lambda = zeta;
		} else {
		    rtmpa = lb * .6875 + ub * .3125;
		    rtmpb = lb * .3125 + ub * .6875;
		    if (rtmpa <= xlambda && xlambda <= rtmpb) {
			*lambda = xlambda;
		    } else {
			biscount = 4;
		    }
		}
	    }
	}

/*           Record optimal shift */

	if (! haveopt || resid < optres) {
	    optlda = bcklda;
	    optres = resid;
	    haveopt = TRUE_;
/*              Keep this one as optimum, switch buffers */
	    bufsel = 3 - bufsel;
	}

/*           =========== */
/*            LOOP EXIT */
/*           =========== */

	if (status > 0) {
	    goto L90002;
	}

/*           Setup a bisection iterate */

	if (biscount != 0) {
	    *lambda = (lb + ub) * .5;
	}

    }

    goto L90001;
L90002:

/*     The loop was only exited after a vector was computed into the active */
/*     buffer. However, the last computed one need not be the optimum if */
/*     the iteration stopped due to maximum accuracy being reached. The */
/*     vector with optimal residual is kept in the backup buffer. */

    *lambda = optlda;
    bufsel = 3 - bufsel;
    i__1 = *n;
    for (ixf77a = 1; ixf77a <= i__1; ++ixf77a) {
	z__[ixf77a] = rwork[ixbuf[bufsel - 1] - 1 + ixf77a];
/* L90005: */
    }
/* L90006: */

/*     Determine the support */

L100:
    i__ = 1;
L90007:
    if (z__[i__] != 0.) {
	goto L90008;
    }
    ++i__;
    goto L90007;
L90008:
    isuppz[1] = i__;
    i__ = *n;
L90009:
    if (z__[i__] != 0.) {
	goto L90010;
    }
    --i__;
    goto L90009;
L90010:
    isuppz[2] = i__;
    xmrstats_1.xnrqi += nrqi;
    xmrstats_1.xnrqibis += nbis;
    xmrstats_1.xmaxnrqi = max(xmrstats_1.xmaxnrqi,nrqi);
    xmrstats_1.xmaxnrqibis = max(xmrstats_1.xmaxnrqibis,nbis);
    return 0;
} /* dlaxrx_ */

