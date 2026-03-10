/* dlaxrf_env_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxrf_env__(integer *n, doublereal *e, doublereal *repr,
	 integer *repi, integer *depth, integer *icbeg, integer *icend, 
	doublereal *lgap, doublereal *ugap, integer *ewl_ae__, doublereal *
	ewl_lu__, integer *rginfo, integer *mode, doublereal *env, doublereal 
	*minenv, doublereal *rwork, integer *iwork)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j, k;
    doublereal lb, ub, gapl, gapu;
    integer ixf77a;
    doublereal mienv;
    integer ngrps;
    extern /* Subroutine */ int dlaxrf_grpenv__(integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    integer ixgenv, ixienv;
    logical grpsok;
    integer ixwork;
    logical toowide, toomany;

/*     IMPLICIT NONE */





/*  Purpose */
/*  ======= */

/*     Determine an envelope for the given cluster ICBEG:ICEND. */

/*     Modes: */
/*      1 - do full env (outer and inner, if possible) */
/*      2 - do only outer */
/*      3 - improve outer, meaning ENV is set to outer envelope already, */
/*          but now inner envs shall improve it */

/*  Notes */
/*  ===== */

/*     We use ideas from Parlett&voemel (submatrix-method); check */
/*     their paper about envelope localization. */

/*     For even trying to do an envelope we require that a (sub-)cluster */
/*     I:J fulfills */
/*     (1)  |J-I+1| < N/3 */
/*     (2)  width < mingap / ENVGAPFAC */

/*  ====================================================================== */

/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRF_GRPENV( */
/*    $            N, E, REPR, REPI, I, J, */
/*    $            LGAP, LB, UB, UGAP, */
/*    $            ENV, */
/*    $            RWORK, IWORK */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, I, J */
/*     INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  LGAP, LB, UB, UGAP */
/*     DOUBLE PRECISION, INTENT(IN)  ::  E(N-1), REPR(4*N+3) */

/*     INTEGER,          INTENT(INOUT)  ::  IWORK( 2*N ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( 6*N ) */

/*     DOUBLE PRECISION, INTENT(OUT)  ::  ENV(N) */

/*     END SUBROUTINE DLAXRF_GRPENV */
/*     END INTERFACE */

/*     .. Constants .. */


/*     .. Parameters .. */

/*         Cluster boundaries will be refined to relative accuracy */
/*         RACLBFAC*N*EPS */

/*     .. Locals .. */


/*  ===== Executable Statements ========================================== */

    /* Parameter adjustments */
    --iwork;
    --rwork;
    --env;
    --repi;
    --repr;
    --e;
    --ewl_ae__;
    --ewl_lu__;
    --rginfo;

    /* Function Body */
    if (FALSE_) {
	i__1 = *n;
	for (ixf77a = 1; ixf77a <= i__1; ++ixf77a) {
	    env[ixf77a] = 1.;
/* L90001: */
	}
/* L90002: */
	*minenv = 1.;
	return 0;
    }

/*     Real Workspace */

    ixienv = 1;
    ixgenv = *n + 1;
    ixwork = (*n << 1) + 1;
/*       -------- */
/*             2N */
/*           + 6N for GRPENV */
/*     Integer Workspace */
/*             2N for GRPENV */
/*     ------------------------ */
/*      Compute outer envelope */
/*     ------------------------ */
    if (*mode == 1 || *mode == 2) {
	lb = ewl_lu__[(*icbeg << 1) - 1];
	ub = ewl_lu__[*icend * 2];
	toomany = (*icend - *icbeg + 1) * 3 > *n;
	toowide = (ub - lb) * 100 > min(*lgap,*ugap);
	if (toomany || toowide) {
	    i__1 = *n;
	    for (ixf77a = 1; ixf77a <= i__1; ++ixf77a) {
		env[ixf77a] = 1.;
/* L90003: */
	    }
/* L90004: */
	    *minenv = 1.;
	} else {
	    dlaxrf_grpenv__(n, &e[1], &repr[1], &repi[1], icbeg, icend, lgap, 
		    &lb, &ub, ugap, &env[1], &rwork[ixwork], &iwork[1]);
	    *minenv = 1.;
	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
		d__1 = *minenv, d__2 = env[k];
		*minenv = min(d__1,d__2);
/* L90005: */
	    }
/* L90006: */
	}
    }
    if (*mode == 2) {
	return 0;
    }
/*     ----------------------------------- */
/*      Compute envelope from subclusters */
/*     ----------------------------------- */
/*     First pass */
/*     Determine number of subgroups and check that all satisfy the */
/*     necessary criteria to make an envelope useful */
    ngrps = 0;
    grpsok = TRUE_;
    i__ = *icbeg;
L90007:
    j = ewl_ae__[i__ * 2];
L90009:
    if (rginfo[j] > 0) {
	goto L90010;
    }
    j = ewl_ae__[(j + 1) * 2];
    goto L90009;
L90010:
    ++ngrps;

    lb = ewl_lu__[(i__ << 1) - 1];
    ub = ewl_lu__[j * 2];
    if (i__ == *icbeg) {
	gapl = *lgap;
    } else {
	gapl = lb - ewl_lu__[(i__ << 1) - 2];
    }
    if (j == *icend) {
	gapu = *ugap;
    } else {
	gapu = ewl_lu__[(j << 1) + 1] - ub;
    }
    toomany = (j - i__ + 1) * 3 > *n;
    toowide = (ub - lb) * 100 > min(gapl,gapu);
    grpsok = grpsok && ! (toomany || toowide);

    if (j == *icend || ! grpsok) {
	goto L90008;
    }
    i__ = j + 1;
    goto L90007;
L90008:
    if (! grpsok) {
	return 0;
    }
    if (ngrps == 1) {
	return 0;
    }
/*     Second pass */
/*     Accumulate envelopes for the subgroups */
    i__1 = ixienv + *n - 1;
    for (ixf77a = ixienv; ixf77a <= i__1; ++ixf77a) {
	rwork[ixf77a] = 0.;
/* L90011: */
    }
/* L90012: */
    mienv = 0.;

    i__ = *icbeg;
L90013:
    j = ewl_ae__[i__ * 2];
L90015:
    if (rginfo[j] > 0) {
	goto L90016;
    }
    j = ewl_ae__[(j + 1) * 2];
    goto L90015;
L90016:
    lb = ewl_lu__[(i__ << 1) - 1];
    ub = ewl_lu__[j * 2];
    if (i__ == *icbeg) {
	gapl = *lgap;
    } else {
	gapl = lb - ewl_lu__[(i__ << 1) - 2];
    }
    if (j == *icend) {
	gapu = *ugap;
    } else {
	gapu = ewl_lu__[(j << 1) + 1] - ub;
    }
    dlaxrf_grpenv__(n, &e[1], &repr[1], &repi[1], &i__, &j, &gapl, &lb, &ub, &
	    gapu, &rwork[ixgenv], &rwork[ixwork], &iwork[1]);

/*        accumulate */

    mienv = 0.;
    i__1 = *n - 1;
    for (k = 0; k <= i__1; ++k) {
/* Computing 2nd power */
	d__1 = rwork[ixgenv + k];
	rwork[ixienv + k] += d__1 * d__1;
/* Computing MIN */
	d__1 = mienv, d__2 = rwork[ixienv + k];
	mienv = min(d__1,d__2);
/* L90017: */
    }
/* L90018: */
    if (mienv >= 1. || j == *icend) {
	goto L90014;
    }
    i__ = j + 1;
    goto L90013;
L90014:
    if (mienv >= 1.) {
	return 0;
    }
/*     at this point, the inner envelope holds the squares of the */
/*     entries of a valid envelope approximation */
    *minenv = 1.;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	d__1 = env[k], d__2 = sqrt(rwork[ixienv - 1 + k]);
	env[k] = min(d__1,d__2);
/* Computing MIN */
	d__1 = *minenv, d__2 = env[k];
	*minenv = min(d__1,d__2);
/* L90019: */
    }
/* L90020: */
    return 0;
} /* dlaxrf_env__ */

