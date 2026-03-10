/* dlaxrf_seltw_f77.f -- translated by f2c (version 20240504).
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

static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int dlaxrf_seltw__(integer *n, integer *kminus, doublereal *
	d__, integer *omegad, doublereal *r__, integer *omegar, doublereal *e,
	 doublereal *gamma, integer *twistok, doublereal *env, doublereal *
	spdiam, doublereal *mingap, integer *dir, integer *loc, doublereal *
	tau, logical *winok, integer *k, doublereal *eval, integer *xizero, 
	doublereal *rwork, integer *iwork)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal invgrbnd, invmaxrc, failxrrr;
    integer i__;
    logical actinerok;
    doublereal gmathresh;
    integer xi;
    extern /* Subroutine */ int dlaxrf_seltw_part__(integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    doublereal belg_noenv__;
    integer fok, lok;
    doublereal eps, growthbound, belg, delg, feats[6];
    integer ixf77a, ninok;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal brcond, drcond;
    integer ixbeta;
    logical inerok;
    integer ixgmax, iyiner, ixrcsq, twinok;
    doublereal acteval;
    integer ixmxbrc, ixgsmax;

/*     IMPLICIT NONE */




/*  Purpose */
/*  ======= */

/*     Selects a twist with optimal element growth and relative condition */
/*     number. For this twist, the characteristics of the resulting rep */
/*     are determined and returned in ELG, RCOND, NCD and combined as */
/*     EVAL. */

/*     If all entries in TWISTOK are 0, the routine sets K=-1 and returns. */

/*  Arguments */
/*  ========= */

/*  KMINUS  (input)  INTEGER */
/*          Twist in the source. */

/*  TWISTOK (input)  INTEGER array, dimension (N) */
/*          Flags where a twist is allowed, as delivered by DLAXRS. */

/*  TAU     (input)  DOUBLE PRECISION */
/*  LOC     (input)  INTEGER, in {1,N} */
/*  DIR     (input)  INTEGER, in {-1,+1} */
/*          Specify the eigenvalue we are close to, on which side, and */
/*          the shift itself. */

/*  WINOK   (input)  LOGICAL */
/*          Set to true if the resulting shift shall have a consistent */
/*          inertia wrt LOC and DIR. The routine will try to get a */
/*          twist with consistent inertia anyway, but the impact of */
/*          WINOK is that it causes only twists with a consistent */
/*          inertia to be considered in the first place. */

/*  GAMMA   (input)  DOUBLE PRECISION array, dimension (N) */
/*          In concert with looking for consistent inertias, the routine */
/*          may set very small gammas (wrt |tau|) to zero. */

/*  K       (output)  INTEGER */
/*          The selected twist. */
/*          Set to <= 0 if no suitable twist was found: */
/*           0 if no twist had a consistent inertia when WINOK=.TRUE. */
/*          -1 if no twist was possible */

/*  EVAL    (output)  DOUBLE PRECISION, >= 0 */
/*          Evaluation, combined from the normalised single features. */
/*          Is <= 1 iff the rep passes the a priori criteria. */
/*          Only set if K!=-1. */

/*          The detailed features of the representation that constitute */
/*          eval are returned in the first entries of RWORK. */
/*          At the moment these are: */
/*            RWORK(1) = element growth without envelope / spdiam */
/*            RWORK(2) = element growth with envelope / spdiam */
/*            RWORK(3) = -1 [reserved for relcond without envelope] */
/*            RWORK(4) = relcond with envelope */
/*            RWORK(5) = -1 [reserved for ncd] */
/*            RWORK(6) = -1 [reserved for ncd] */

/*  ====================================================================== */

/*     .. Constants .. */


/*     .. Parameters .. */

/*     Restrict possible choice of twist for node reps: */
/*      0 - no restriction */
/*      1 - always twist at 1 */
/*      2 - always twist at n */

/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRF_SELTW_PART( */
/*    $             N, IA, IE, DIR, G, GN, OMEGA, TWISTOK, S, */
/*    $             ARCSQ, ABETA, AGMAX, AGSMAX, ABRCFAC */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER, INTENT(IN)  ::  N, IA, IE, DIR */
/*     INTEGER, INTENT(IN)  :: */
/*    $   TWISTOK( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ), */
/*    $   OMEGA(   MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ) */
/*     DOUBLE PRECISION, INTENT(IN) :: */
/*    $   G(  MIN(IA,IE) : MAX(IA,IE) ), */
/*    $   GN( MIN(IA,IE) : MAX(IA,IE) ), */
/*    $   S( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ) */

/*     DOUBLE PRECISION, INTENT(INOUT)  :: */
/*    $   ARCSQ ( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ), */
/*    $   ABETA ( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ), */
/*    $   AGMAX ( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ), */
/*    $   AGSMAX( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ), */
/*    $   ABRCFAC( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ) */
/*     END SUBROUTINE DLAXRF_SELTW_PART */
/*     END INTERFACE */

/*     .. Local Variables .. */


/*  ===== Executable Statements ========================================== */

/*     Determine range of possible twists */

    /* Parameter adjustments */
    --iwork;
    --rwork;
    --env;
    --twistok;
    --gamma;
    --e;
    --omegar;
    --omegad;
    --d__;
    /* R is declared as R(2:N) in Fortran, adjust by 2 not 1 */
    r__ -= 2;

    /* Function Body */
    fok = 1;
L90001:
    if (twistok[fok] != 0) {
	goto L90002;
    }
    ++fok;
    if (fok == *n + 1) {
	goto L90002;
    }
    goto L90001;
L90002:
    lok = *n;
L90003:
    if (twistok[lok] != 0) {
	goto L90004;
    }
    --lok;
    if (lok == 0) {
	goto L90004;
    }
    goto L90003;
L90004:

/*     Apply twist restriction */

    if (FALSE_) {
	lok = 1;
    }
    if (FALSE_) {
	fok = *n;
    }

/*     Return if no twist possible */

    if (fok > lok) {
	*k = -1;
	return 0;
    }
/*     ---------------------------- */
/*      Compute per-twist Inertias */
/*     ---------------------------- */
    iyiner = 1;
    i__1 = iyiner + *n + 1;
    for (ixf77a = iyiner; ixf77a <= i__1; ++ixf77a) {
	iwork[ixf77a] = 0;
/* L90005: */
    }
/* L90006: */
    eps = dlamch_("Epsilon", (ftnlen)7);
    gmathresh = eps * abs(*tau);
/*     We compute in two sweeps the inertias for the parts above */
/*     and below the twist candidate in IWORK(IYINER + [FOK:LOK]). */
/*     Entries where a block ends are not set, except if the block */
/*     is at the end. */
    xi = 0;
    i__ = 1;
L90007:
    iwork[iyiner + i__] = xi;
    if (i__ >= lok) {
	goto L90008;
    }
    if (omegad[i__ + 1] != 0) {
	xi += 2;
	++i__;
    } else if (d__[i__] < 0.) {
	xi += 2;
    }
    ++i__;
    goto L90007;
L90008:
    if (lok == *n && omegad[*n] != 0) {
	iwork[iyiner + *n] = xi;
    }

    xi = 0;
    i__ = *n;
L90009:
    iwork[iyiner + i__] += xi;
    if (i__ <= fok) {
	goto L90010;
    }
    if (omegar[i__ - 1] != 0) {
	xi += 2;
	--i__;
    } else if (r__[i__] < 0.) {
	xi += 2;
    }
    --i__;
    goto L90009;
L90010:
    if (fok == 1 && omegar[1] != 0) {
	iwork[iyiner + 1] = xi;
    }

/*     Note: The two special case handlers after the loops above already */
/*     set correct inertias for blocks at end, so we can ignore indices */
/*     where blocks end now. */
    twinok = -1;
    ninok = 0;
    i__1 = lok;
    for (i__ = fok; i__ <= i__1; ++i__) {
	if (omegad[i__] == 0 && omegar[i__] == 0) {
	    if ((d__1 = gamma[i__], abs(d__1)) <= gmathresh) {
/*              Set tiny gammas to zero (dlaxrs may or may not have done */
/*              that already). */
/*              Note that the inertia does not yet incorporate gamma. */
		if (*dir == 0 || iwork[iyiner + i__] == (*loc << 1) - 2) {
/*                 We will land exactly on ew LOC */
		    gamma[i__] = 0.;
		} else {
/*                 May already be to far inside, try to rescue what */
/*                 we can */
		    gamma[i__] = -(*dir) * gmathresh;
		}
	    }
	    if (gamma[i__] < 0.) {
		iwork[iyiner + i__] += 2;
	    } else if (gamma[i__] == 0.) {
		++iwork[iyiner + i__];
	    }
	}
	if (twistok[i__] != 0) {
	    xi = iwork[iyiner + i__];
	    if (*dir == -1 && xi <= (*loc << 1) - 1 || *dir == 1 && xi >= (*
		    loc << 1) - 1) {
		++ninok;
		if (twinok == -1 || (i__2 = i__ - *kminus, abs(i__2)) < (i__3 
			= twinok - *kminus, abs(i__3))) {
		    twinok = i__;
		}
	    }
	}
/* L90011: */
    }
/* L90012: */
/*     Now, if TWINOK != -1 then it is set to a twist with consistent */
/*     inertia that is as close as possible to the source twist. */
    if (*winok && ninok == 0) {
/*        None of the twists give a consistent inertia */
	*k = 0;
	return 0;
    }
/*     ------------------------------------ */
/*      Special case: Take guaranteed RRRs */
/*     ------------------------------------ */
    if (twinok != -1 && (*loc == 1 && *dir == -1 || *loc == *n && *dir == 1)) 
	    {
/*        Rationale for eval: If definite then entries are all bounded */
/*        by the norm. */
	*k = twinok;
	*eval = .125;
	*xizero = iwork[iyiner + twinok];
	feats[0] = 1. / *spdiam;
	feats[1] = feats[0];
	feats[3] = 1.;
	return 0;
    }

/*     ---------- */
/*      Evaluate */
/*     ---------- */

    ixrcsq = 0;
    ixbeta = *n;
    ixgmax = *n << 1;
    ixgsmax = *n * 3;
    ixmxbrc = *n << 2;

    i__1 = *n * 5;
    for (ixf77a = 1; ixf77a <= i__1; ++ixf77a) {
	rwork[ixf77a] = 0.;
/* L90013: */
    }
/* L90014: */

/*     --------------- */
/*      Top-to-Bottom */
/*     --------------- */

    if (lok > 1) {
	i__1 = lok - 1;
	dlaxrf_seltw_part__(n, &c__1, &i__1, &c__1, &d__[1], &e[1], &omegad[1]
		, &twistok[1], &env[1], &rwork[ixrcsq + 1], &rwork[ixbeta + 1]
		, &rwork[ixgmax + 1], &rwork[ixgsmax + 1], &rwork[ixmxbrc + 1]
		);
    }

/*     --------------- */
/*      Bottom-to-Top */
/*     --------------- */

    if (fok < *n) {
	i__1 = fok + 1;
	dlaxrf_seltw_part__(n, n, &i__1, &c_n1, &r__[fok + 1], &e[fok], &
		omegar[fok], &twistok[fok], &env[fok], &rwork[ixrcsq + fok], &
		rwork[ixbeta + fok], &rwork[ixgmax + fok], &rwork[ixgsmax + 
		fok], &rwork[ixmxbrc + fok]);
    }

/*     --------- */
/*      Combine */
/*     --------- */

    growthbound = *spdiam * 8;
    invgrbnd = 1. / growthbound;
    invmaxrc = .10000000000000001;
    failxrrr = (*n - 1) * *mingap / (*spdiam * sqrt(eps));
    for (ixf77a = 1; ixf77a <= 6; ++ixf77a) {
	feats[ixf77a - 1] = -1.;
/* L90015: */
    }
/* L90016: */
    *k = -1;
    *eval = -1.;
    *xizero = -1;
    inerok = FALSE_;
    i__1 = lok;
    for (i__ = fok; i__ <= i__1; ++i__) {
	if (twistok[i__] != 0) {
/*           Compute rc and elg that would result from choosing this */
/*           twist. */
	    brcond = rwork[ixrcsq + i__] + (d__1 = env[i__] * (rwork[ixbeta + 
		    i__] + env[i__]), abs(d__1));
/* Computing MAX */
	    d__1 = brcond * rwork[ixmxbrc + i__];
	    brcond = max(d__1,1.);
	    drcond = brcond * invmaxrc;
/* Computing MAX */
	    d__2 = rwork[ixgmax + i__], d__3 = (d__1 = gamma[i__], abs(d__1));
	    belg_noenv__ = max(d__2,d__3);
/* Computing MAX */
	    d__2 = rwork[ixgsmax + i__], d__3 = (d__1 = env[i__] * gamma[i__],
		     abs(d__1));
	    belg = max(d__2,d__3);
/*            IF( GROWTHBOUND.LT.BELG .AND. BELG.LT.FAILXRRR ) */
/*     $      THEN */
	    delg = belg;
/*            ELSE */
/*               DELG = BELG_NOENV */
/*            ENDIF */
	    delg *= invgrbnd;
	    acteval = drcond + delg;
	    xi = iwork[iyiner + i__];
	    actinerok = *dir == -1 && xi <= (*loc << 1) - 1 || *dir == 1 && 
		    xi >= (*loc << 1) - 1;
	    if (*k == -1 || acteval < *eval && (! inerok || actinerok || 
		    acteval * 2 < *eval)) {
		*k = i__;
		*eval = acteval;
		*xizero = xi;
		inerok = actinerok;
		feats[0] = belg_noenv__ / *spdiam;
		feats[1] = belg / *spdiam;
		feats[3] = brcond;
	    }
	}
/* L90017: */
    }
/* L90018: */
    if (*k != -1) {
	for (ixf77a = 1; ixf77a <= 6; ++ixf77a) {
	    rwork[ixf77a] = feats[ixf77a - 1];
/* L90019: */
	}
/* L90020: */
    }
    return 0;
} /* dlaxrf_seltw__ */

