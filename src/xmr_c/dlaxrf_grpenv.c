/* dlaxrf_grpenv_f77.f -- translated by f2c (version 20240504).
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
static doublereal c_b8 = 0.;

/* Subroutine */ int dlaxrf_grpenv__(integer *n, doublereal *e, doublereal *
	repr, integer *repi, integer *i__, integer *j, doublereal *lgap, 
	doublereal *lb, doublereal *ub, doublereal *ugap, doublereal *env, 
	doublereal *rwork, integer *iwork)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal envforce;
    integer k;
    doublereal sl, mu, su, invminrgap, eps;
    integer ixdp;
    doublereal xenv, rtmp;
    integer ixrp;
    doublereal delta, resid, normz;
    integer twist;
    doublereal bisacc;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal absgma, bckoff, mingap, sinbnd;
    integer ixgmal;
    extern /* Subroutine */ int dlaxrg_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    integer ixgmau;
    extern /* Subroutine */ int dlaxrt_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *);
    doublereal invrho, rqcorr;
    integer jxtokl, ixwork, jxtoku, jxwork;

/*     IMPLICIT NONE */





/*  Purpose */
/*  ======= */

/*     Compute envelope for (sub-)cluster I:J. */
/*     It is *not* checked beforehand if it makes sense to try for an */
/*     envelope here, that should have been done already. */
/*     We use the Parlett/Voemel envelope strategy, or a plain */
/*     twisted facto for a singleton subcluster. */

/*  ====================================================================== */

/*     .. Declarations .. */

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

/*     .. Constants .. */


/*     .. Parameters .. */

/*         Intervals will be refined via bisection up to a relative */
/*         accuracy of BISACCFAC*N*EPS */

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
    --iwork;
    --rwork;
    --env;
    --repi;
    --repr;
    --e;

    /* Function Body */
    eps = dlamch_("Precision", (ftnlen)9);
    bisacc = (*n << 2) * eps;
    mingap = min(*lgap,*ugap);

/*     Integer Workspace */

    jxtokl = 1;
    jxtoku = jxtokl + *n;
    jxwork = jxtoku + *n;
/*       --------- */
/*              2N */

/*     Real Workspace */

    ixdp = 1;
    ixrp = ixdp + *n;
    ixgmal = ixrp + *n;
    ixgmau = ixgmal + *n;
    ixwork = ixgmau + *n;
/*       --------- */
/*              4N */
/*            + 2N for DLAXRT */
/*       ========= */
/*              6N */
/*     -------------------- */
/*      Build the envelope */
/*     -------------------- */
    if (*i__ == *j) {
/*        ---------------------- */
/*         Singleton Subcluster */
/*        ---------------------- */
	d__1 = (*lb + *ub) * .5;
	dlaxrt_(n, &c__1, n, &repr[1], &repi[1], &d__1, &rwork[ixdp], &rwork[
		ixrp], &rwork[ixgmal], &iwork[1], &rwork[ixwork]);
/*        get min gamma */
	iwork[*n + 1] = -1;
	twist = 1;
L90001:
	if (iwork[twist] != 0) {
	    goto L90002;
	}
	++twist;
	goto L90001;
L90002:
	absgma = (d__1 = rwork[ixgmal - 1 + twist], abs(d__1));
	i__1 = *n;
	for (k = twist + 1; k <= i__1; ++k) {
	    if (iwork[k] != 0) {
		rtmp = (d__1 = rwork[ixgmal - 1 + k], abs(d__1));
		if (rtmp < absgma) {
		    absgma = rtmp;
		    twist = k;
		}
	    }
/* L90003: */
	}
/* L90004: */
/*        Compute eigenvector approximation */
	dlaxrg_(n, &twist, &rwork[ixdp], &rwork[ixrp - 2 + twist + 1], &e[1], 
		&rwork[ixgmal - 1 + twist], &c_b8, &rwork[ixwork], &normz, &
		resid, &rqcorr);
	sinbnd = resid / mingap;
/*        Extract envelope */
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	    d__2 = 1., d__3 = (d__1 = rwork[ixwork - 1 + k], abs(d__1)) + 
		    sinbnd;
	    env[k] = min(d__2,d__3);
/* L90005: */
	}
/* L90006: */
/* @LOGGING on */
/*         WRITE(FIDRRF,*) '   -- Singleton --' */
/*         WRITE(FIDRRF,*) '     absgma =',ABSGMA */
/*         WRITE(FIDRRF,*) '     twist  =',TWIST */
/*         WRITE(FIDRRF,*) '     sinbnd =',SINBND */
/* @LOGGING ! */
    } else {
/*        ---------- */
/*         Subgroup */
/*        ---------- */
/*        We use one parameter BCKOFF in place of TAU*DELTA in the */
/*        Parlett/Voemel paper. To avoid zero gammas we backoff at */
/*        least a bit. */
	mu = (*lb + *ub) * .5;
	delta = (*ub - *lb) * .5;
/* Computing MAX */
/* Computing MAX */
	d__3 = abs(*lb), d__4 = abs(*ub);
	d__1 = delta * 2, d__2 = delta + bisacc * max(d__3,d__4);
	bckoff = max(d__1,d__2);
/* Computing 2nd power */
	d__1 = (mingap + delta) / bckoff;
	invrho = 1. / (d__1 * d__1 - 1.);
	sl = mu - bckoff;
	su = mu + bckoff;
/*        Compute twisted factorizations at both ends. */
/*        Only the gammas are needed. */
	dlaxrt_(n, &c__1, n, &repr[1], &repi[1], &sl, &rwork[ixdp], &rwork[
		ixrp], &rwork[ixgmal], &iwork[jxtokl], &rwork[ixwork]);
	dlaxrt_(n, &c__1, n, &repr[1], &repi[1], &su, &rwork[ixdp], &rwork[
		ixrp], &rwork[ixgmau], &iwork[jxtoku], &rwork[ixwork]);
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    if (iwork[jxtokl - 1 + k] == 0 || iwork[jxtoku - 1 + k] == 0 || 
		    rwork[ixgmal - 1 + k] == 0. || rwork[ixgmau - 1 + k] == 
		    0.) {
		xenv = 1.;
	    } else {
/*              The smallest weight inside the cluster is at least */
/*              2 / BCKOFF. */
		rtmp = bckoff * .5;
		rtmp = (d__1 = rtmp / rwork[ixgmal - 1 + k] + rtmp / rwork[
			ixgmau - 1 + k], abs(d__1)) + (*n - (*j - *i__ + 1)) *
			 invrho;
/* Computing MIN */
		d__1 = 1., d__2 = sqrt(rtmp);
		xenv = min(d__1,d__2);
	    }
	    env[k] = xenv;
/* L90007: */
	}
/* L90008: */
    }

/*     Safeguard by forcing a minimal envelope entry. The rationale */
/*     here is that for a small change in the representation, a */
/*     componentwise change of */
/*       n eps relcond / relgap */
/*     has to be anticipated in each vector, so the envelope entries */
/*     really cannot fall below that. */

/* Computing MAX */
    d__1 = abs(*lb) / *lgap, d__2 = abs(*ub) / *ugap;
    invminrgap = max(d__1,d__2);
/* Computing MIN */
    d__1 = 1., d__2 = *n * 10 * eps * invminrgap;
    envforce = min(d__1,d__2);
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
	d__1 = envforce, d__2 = env[k];
	env[k] = max(d__1,d__2);
/* L90009: */
    }
/* L90010: */
/* @LOGGING on */
/*      ENVMIN = ONE */
/*      ENVMAX = ZERO */
/*      ENVAVG = ZERO */
/*      DO K = 1, N */
/*         ENVMIN = MIN(ENVMIN, ENV(K)) */
/*         ENVMAX = MAX(ENVMAX, ENV(K)) */
/*         ENVAVG = ENVAVG + ENV(K) */
/*      ENDDO */
/*      ENVAVG = ENVAVG / N */
/*      WRITE(FIDRRF,*) '   -- Summary --' */
/*      WRITE(FIDRRF,*) '     envforce = ',ENVFORCE */
/*      WRITE(FIDRRF,*) '     min = ',ENVMIN */
/*      WRITE(FIDRRF,*) '     avg = ',ENVAVG */
/*      WRITE(FIDRRF,*) '     max = ',ENVMAX */
/*      WRITE(FIDRRF,*) ' ' */
/* @LOGGING ! */

    if (*i__ == *j) {
	++xmrstats_1.xnenvgv;
	++xmrstats_1.xnenvtf;
    } else {
	xmrstats_1.xnenvtf += 2;
    }
    return 0;
} /* dlaxrf_grpenv__ */

