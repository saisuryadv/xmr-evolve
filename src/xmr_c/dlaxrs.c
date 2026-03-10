/* dlaxrs_f77.f -- translated by f2c (version 20240504).
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
static integer c__2 = 2;

/* Subroutine */ int dlaxrs_(integer *n, integer *ia, integer *ie, doublereal 
	*repr, integer *repi, doublereal *tau, doublereal *shfprt, doublereal 
	*dplus, integer *omgadp, doublereal *rplus, integer *omgarp, 
	doublereal *gammap, integer *twistok, doublereal *rwork)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__, k, j1, j2, nb;
    doublereal tn, tp;
    integer ibb;
    doublereal gma, eps;
    integer ixg;
    extern /* Subroutine */ int dlaxrs_prog__(integer *, integer *, integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *)
	    , dlaxrs_stat__(integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , doublereal *);
    doublereal bigg;
    integer ixtn, ixtp, ixf77a, ixngn;
    extern doublereal dlamch_(char *, ftnlen);
    integer iylbbk;
    doublereal abstau;
    integer ixbdet, iyomga;
    doublereal pivmin;
    integer ixgnsq;
    doublereal pivbase;

/*     IMPLICIT NONE */




/*  Purpose */
/*  ======= */

/*     Factorize N G N' - TTAU = Np Gp Np', where both NGN and NpGpNp' are */
/*     twisted block factorizations, and TTAU is a diagonal matrix. */
/*     The union of data for all twists in IA:IE in the target is computed. */
/*   ! At the moment progressive par only implemented for non-blocked ! */

/*     There are two reasons why some twist may not be ok: a block would */
/*     end there (we don't allow that for now), or we cannot compute the */
/*     twist element stably. 1 and N are always ok. */

/*     Post conditions */
/*     --------------- */
/*       if IA = 1 then TWISTOK(1) = 1 */
/*       if IE = N then TWISTOK(N) = 1 */
/*       OMGADP(1) = 0 */
/*       OMGARP(N) = 0 */

/*  ====================================================================== */

/*     .. Parameters .. */


/*     PIVTAUFAC */
/*         Adjust PIVMIN to PIVTAUFAC*EPS*TAU, set to zero to deactivate */
/*         a shift-dependent pivmin. */


/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRS_STAT( */
/*    $             N, J, I0, I1, G, OMEGA, GNSQ, NGN, BDET, */
/*    $             PIVMIN, TAU, SHFPRT, GPLUS, OPLUS, S */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, J, I0, I1 */
/*     INTEGER,          INTENT(IN)  ::  OMEGA(0:N+1) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU */
/*     DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  SHFPRT(N) */

/*     INTEGER,          INTENT(OUT)  ::  OPLUS(N) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  GPLUS(I0:I1), S(N) */
/*     END SUBROUTINE DLAXRS_STAT */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRS_PROG( */
/*    $             N, J1, J2, I0, I1, G, OMEGA, GNSQ, NGN, BDET, */
/*    $             PIVMIN, TAU, SHFPRT, OFF, GPLUS, OPLUS, P */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, J1, J2, I0, I1 */
/*     INTEGER,          INTENT(IN)  ::  OMEGA(0:N+1) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU, OFF */
/*     DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  SHFPRT(N) */

/*     INTEGER,          INTENT(OUT)  ::  OPLUS(N) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  GPLUS(I0:I1), P(N) */
/*     END SUBROUTINE DLAXRS_PROG */
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
    --rwork;
    --omgarp;
    --omgadp;
    --shfprt;
    --repi;
    --repr;
    --dplus;
    --rplus;
    /* GAMMAP is declared as GAMMAP(IA:IE) in Fortran, adjust by IA */
    gammap -= *ia;
    /* TWISTOK is declared as TWISTOK(IA:IE) in Fortran, adjust by IA */
    twistok -= *ia;

    /* Function Body */
    if (! xmrstats_1.xstealthmode) {
	if (*ia != *ie) {
	    ++xmrstats_1.xnumfs_2__;
	} else {
	    ++xmrstats_1.xnumfs_k__;
	}
    }

/*     .. Decode Representation Data .. */

    ixg = 0;
    ixbdet = *n;
    ixngn = (*n << 1) + 1;
    ixgnsq = *n * 3 + 2;
    iyomga = 4;
    iylbbk = *n + 6;
    pivbase = repr[(*n << 2) + 3];
    k = repi[2];
    nb = repi[3];
    eps = dlamch_("Epsilon", (ftnlen)7);
    bigg = dlamch_("Overflow", (ftnlen)8);
    abstau = abs(*tau);
/*     Restrict range for twists if we have blocks .. */
/*     (progressive only implemented for no blocks in source yet) */
    j1 = *ia;
    j2 = *ie;
    if (nb > 0) {
/*        We cannot do a progressive part that would range over a */
/*        block in the source. */
	i__1 = nb;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    ibb = repi[iylbbk + i__];
	    if (ibb < k) {
/* Computing MAX */
/* Computing MIN */
		i__4 = k, i__5 = ibb + 2;
		i__2 = j1, i__3 = min(i__4,i__5);
		j1 = max(i__2,i__3);
	    }
	    if (ibb > k) {
/* Computing MIN */
/* Computing MAX */
		i__4 = k, i__5 = ibb - 2;
		i__2 = j2, i__3 = max(i__4,i__5);
		j2 = min(i__2,i__3);
	    }
/* L90001: */
	}
/* L90002: */
    }
/* Computing MAX */
    d__1 = pivbase, d__2 = eps * 0. * abstau;
    pivmin = max(d__1,d__2);
    ixtn = 1;
    ixtp = *n + 1;
    i__1 = min(k,j2);
    i__2 = *n - 1;
    dlaxrs_stat__(n, &i__1, &c__1, &i__2, &repr[ixg + 1], &repi[iyomga], &
	    repr[ixgnsq + 1], &repr[ixngn + 1], &repr[ixbdet + 1], &pivmin, 
	    tau, &shfprt[1], &dplus[1], &omgadp[1], &rwork[ixtn]);
    if (k < j2) {
	tn = rwork[ixtn - 1 + k];
	i__1 = *n - 1;
	dlaxrs_prog__(n, &k, &j2, &c__1, &i__1, &repr[ixg + 1], &repi[iyomga],
		 &repr[ixgnsq + 1], &repr[ixngn + 1], &repr[ixbdet + 1], &
		pivmin, tau, &shfprt[1], &tn, &dplus[1], &omgadp[1], &rwork[
		ixtn]);
    }
    i__1 = max(k,j1);
    dlaxrs_stat__(n, &i__1, &c__2, n, &repr[ixg + 1], &repi[iyomga], &repr[
	    ixgnsq + 1], &repr[ixngn + 1], &repr[ixbdet + 1], &pivmin, tau, &
	    shfprt[1], &rplus[1], &omgarp[1], &rwork[ixtp]);
    if (j1 < k) {
	tp = rwork[ixtp - 1 + k];
	dlaxrs_prog__(n, &k, &j1, &c__2, n, &repr[ixg + 1], &repi[iyomga], &
		repr[ixgnsq + 1], &repr[ixngn + 1], &repr[ixbdet + 1], &
		pivmin, tau, &shfprt[1], &tp, &rplus[1], &omgarp[1], &rwork[
		ixtp]);
    }

/*     Compute gammas */

/*     Notes: */
/*     (1)  If one of the parts places a block end at i, 1 < i < n, then */
/*          we do not want to twist there, as that would again expose */
/*          the elg that was hidden by the block. */
/*      ==> it may well be that no twist is ok */
/*     (2)  For twisting at i where a block ends in the source, the */
/*          same special handling as in dlaxrt does apply. */

    i__1 = *ie;
    for (ixf77a = *ia; ixf77a <= i__1; ++ixf77a) {
	twistok[ixf77a] = 0;
/* L90003: */
    }
/* L90004: */
    i__1 = j2;
    for (i__ = j1; i__ <= i__1; ++i__) {
	gma = bigg;
	if (i__ != *n && omgadp[i__] != 0 || i__ != 1 && omgarp[i__] != 0) {
	    twistok[i__] = 0;
	} else if (repi[iyomga + i__] == 0 || i__ == 1 || i__ == *n) {
/*           neither source nor target have a block here */
	    tn = rwork[ixtn - 1 + i__];
	    tp = rwork[ixtp - 1 + i__];
	    if (i__ != k) {
		gma = tn + tp - *tau * shfprt[i__];
	    } else {
		gma = repr[ixg + k] + (tn + tp - *tau * shfprt[i__]);
	    }
	    twistok[i__] = 1;
	} else {
/*           block ends in the source */
/*            TODO, see dlaxrt */
	    twistok[i__] = 0;
	}
	if (twistok[i__] == 1) {
	    gammap[i__] = gma;
	}
/* L90005: */
    }
/* L90006: */
/*      DO I = IA,IE */
/*         IF( (I.LT.J1) .OR. (I.GT.J2) .OR. */
/*     $       (I.NE.N .AND. OMGADP(I).NE.0) .OR. */
/*     $       (I.NE.1 .AND. OMGARP(I).NE.0) ) */
/*     $   THEN */
/*            TWISTOK(I) = 0 */
/*         ELSE */
/* *           normal case (no block ending) */
/*            TN = RWORK(IXTN-1 + I) */
/*            TP = RWORK(IXTP-1 + I) */
/*            GMA = (TN + TP) - TAU*SHFPRT(I) */
/*            IF( I.EQ.K )THEN */
/*               GMA = GMA + REPR(IXG + K) */
/*            ENDIF */
/*            GAMMAP(I)  = GMA */
/*            TWISTOK(I) = 1 */
/*         ENDIF */
/*      ENDDO */
    return 0;
} /* dlaxrs_ */

