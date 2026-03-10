/* dlaxrm_f77.f -- translated by f2c (version 20240504).
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
static integer c_n1 = -1;

/* Subroutine */ int dlaxrm_(integer *n, doublereal *repr, integer *repi, 
	integer *m, doublereal *atau, integer *axi)
{
    integer i__, j, k, nb, xi, ibb;
    doublereal gma;
    integer ixg, negn[64], negp[64];
    extern /* Subroutine */ int dlaxrm_stat2__(integer *, integer *, integer *
	    , doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *)
	    , dlaxrm_stat4__(integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *), 
	    dlaxrm_stat8__(integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *);
    doublereal auxn[64], auxp[64];
    extern /* Subroutine */ int dlaxrm_stat32__(integer *, integer *, integer 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, doublereal *
	    ), dlaxrm_stat16__(integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *), 
	    dlaxrm_stat64__(integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *);
    integer ixngn;
    doublereal gmarep;
    integer jylbbk, ixbdet;
    extern integer dlaxrn_(integer *, doublereal *, integer *, doublereal *
	    );
    doublereal pivmin;
    integer ixgnsq;
    doublereal pivbase;

/*     IMPLICIT NONE */



/*  Purpose */
/*  ======= */

/*    Compute multiple sturm counts with respect to the same */
/*    representation. The routine computes for each shift tau */
/*    in ATAU the inertia of T - tau into AXI---see dlaxrn for */
/*    information on how the inertias are to be interpreted. */

/*    Pre: No zero shift. */

/*  ====================================================================== */

/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRM_STAT64( */
/*    $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET, */
/*    $             PIVMIN, TAU, ANEGC, AAUX */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, K, DIR */
/*     INTEGER,          INTENT(IN)  ::  LBBEGK(*) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(64) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N) */

/*     INTEGER,          INTENT(INOUT)  ::  IBB */

/*     INTEGER,          INTENT(OUT)  ::  ANEGC(64) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(64) */
/*     END SUBROUTINE DLAXRM_STAT64 */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRM_STAT32( */
/*    $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET, */
/*    $             PIVMIN, TAU, ANEGC, AAUX */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, K, DIR */
/*     INTEGER,          INTENT(IN)  ::  LBBEGK(*) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(32) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N) */

/*     INTEGER,          INTENT(INOUT)  ::  IBB */

/*     INTEGER,          INTENT(OUT)  ::  ANEGC(32) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(32) */
/*     END SUBROUTINE DLAXRM_STAT32 */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRM_STAT16( */
/*    $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET, */
/*    $             PIVMIN, TAU, ANEGC, AAUX */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, K, DIR */
/*     INTEGER,          INTENT(IN)  ::  LBBEGK(*) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(16) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N) */

/*     INTEGER,          INTENT(INOUT)  ::  IBB */

/*     INTEGER,          INTENT(OUT)  ::  ANEGC(16) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(16) */
/*     END SUBROUTINE DLAXRM_STAT16 */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRM_STAT8( */
/*    $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET, */
/*    $             PIVMIN, TAU, ANEGC, AAUX */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, K, DIR */
/*     INTEGER,          INTENT(IN)  ::  LBBEGK(*) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(8) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N) */

/*     INTEGER,          INTENT(INOUT)  ::  IBB */

/*     INTEGER,          INTENT(OUT)  ::  ANEGC(8) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(8) */
/*     END SUBROUTINE DLAXRM_STAT8 */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRM_STAT4( */
/*    $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET, */
/*    $             PIVMIN, TAU, ANEGC, AAUX */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, K, DIR */
/*     INTEGER,          INTENT(IN)  ::  LBBEGK(*) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(4) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N) */

/*     INTEGER,          INTENT(INOUT)  ::  IBB */

/*     INTEGER,          INTENT(OUT)  ::  ANEGC(4) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(4) */
/*     END SUBROUTINE DLAXRM_STAT4 */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRM_STAT2( */
/*    $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET, */
/*    $             PIVMIN, TAU, ANEGC, AAUX */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, K, DIR */
/*     INTEGER,          INTENT(IN)  ::  LBBEGK(*) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(2) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N) */

/*     INTEGER,          INTENT(INOUT)  ::  IBB */

/*     INTEGER,          INTENT(OUT)  ::  ANEGC(2) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(2) */
/*     END SUBROUTINE DLAXRM_STAT2 */
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
/* extract -b parameters.inc.f procname=laxrm */

/*     The maximal number of negcounts that are to be performed at */
/*     once. Controls which of the routines dlaxrm_statXX are allowed. */
/*     As such, sensible values are those for which a corresponding */
/*     routine exists in the first place (1,2,4,8,16,32,64 atm). */
/*     Must be >= 1. */


/*     .. Constants .. */


/*     .. Local Variables .. */


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

/*     -- Executable Statements ----------------------------------------- */

/*     Note: Here we use absolute positions to the beginning of data. */
    /* Parameter adjustments */
    --repi;
    --repr;
    --axi;
    --atau;

    /* Function Body */
    ixg = 1;
    ixbdet = *n + 1;
    ixngn = (*n << 1) + 2;
    ixgnsq = *n * 3 + 3;
    jylbbk = *n + 6;
    pivbase = repr[(*n << 2) + 3];
    k = repi[2];
    nb = repi[3];
    gmarep = repr[ixg - 1 + k];
    pivmin = pivbase;
    j = 1;
    if (FALSE_) {
L90001:
	if (j + 63 > *m) {
	    goto L90002;
	}
	ibb = 1;
	dlaxrm_stat64__(n, &k, &c__1, &repr[ixg], &ibb, &repi[jylbbk], &repr[
		ixgnsq], &repr[ixngn], &repr[ixbdet], &pivmin, &atau[j], negn,
		 auxn);
	ibb = nb + 1;
	dlaxrm_stat64__(n, &k, &c_n1, &repr[ixg], &ibb, &repi[jylbbk], &repr[
		ixgnsq], &repr[ixngn], &repr[ixbdet], &pivmin, &atau[j], negp,
		 auxp);
	for (i__ = 1; i__ <= 64; ++i__) {
	    gma = gmarep + (auxn[i__ - 1] + auxp[i__ - 1] - atau[j]);
	    xi = negn[i__ - 1] + negp[i__ - 1] << 1;
	    if (gma < 0.) {
		xi += 2;
	    } else if (gma == 0.) {
		++xi;
	    }
	    axi[j] = xi;
	    ++j;
/* L90003: */
	}
/* L90004: */
	goto L90001;
L90002:
	;
    }
    if (FALSE_) {
L90005:
	if (j + 31 > *m) {
	    goto L90006;
	}
	ibb = 1;
	dlaxrm_stat32__(n, &k, &c__1, &repr[ixg], &ibb, &repi[jylbbk], &repr[
		ixgnsq], &repr[ixngn], &repr[ixbdet], &pivmin, &atau[j], negn,
		 auxn);
	ibb = nb + 1;
	dlaxrm_stat32__(n, &k, &c_n1, &repr[ixg], &ibb, &repi[jylbbk], &repr[
		ixgnsq], &repr[ixngn], &repr[ixbdet], &pivmin, &atau[j], negp,
		 auxp);
	for (i__ = 1; i__ <= 32; ++i__) {
	    gma = gmarep + (auxn[i__ - 1] + auxp[i__ - 1] - atau[j]);
	    xi = negn[i__ - 1] + negp[i__ - 1] << 1;
	    if (gma < 0.) {
		xi += 2;
	    } else if (gma == 0.) {
		++xi;
	    }
	    axi[j] = xi;
	    ++j;
/* L90007: */
	}
/* L90008: */
	goto L90005;
L90006:
	;
    }
    if (FALSE_) {
L90009:
	if (j + 15 > *m) {
	    goto L90010;
	}
	ibb = 1;
	dlaxrm_stat16__(n, &k, &c__1, &repr[ixg], &ibb, &repi[jylbbk], &repr[
		ixgnsq], &repr[ixngn], &repr[ixbdet], &pivmin, &atau[j], negn,
		 auxn);
	ibb = nb + 1;
	dlaxrm_stat16__(n, &k, &c_n1, &repr[ixg], &ibb, &repi[jylbbk], &repr[
		ixgnsq], &repr[ixngn], &repr[ixbdet], &pivmin, &atau[j], negp,
		 auxp);
	for (i__ = 1; i__ <= 16; ++i__) {
	    gma = gmarep + (auxn[i__ - 1] + auxp[i__ - 1] - atau[j]);
	    xi = negn[i__ - 1] + negp[i__ - 1] << 1;
	    if (gma < 0.) {
		xi += 2;
	    } else if (gma == 0.) {
		++xi;
	    }
	    axi[j] = xi;
	    ++j;
/* L90011: */
	}
/* L90012: */
	goto L90009;
L90010:
	;
    }
    if (FALSE_) {
L90013:
	if (j + 7 > *m) {
	    goto L90014;
	}
	ibb = 1;
	dlaxrm_stat8__(n, &k, &c__1, &repr[ixg], &ibb, &repi[jylbbk], &repr[
		ixgnsq], &repr[ixngn], &repr[ixbdet], &pivmin, &atau[j], negn,
		 auxn);
	ibb = nb + 1;
	dlaxrm_stat8__(n, &k, &c_n1, &repr[ixg], &ibb, &repi[jylbbk], &repr[
		ixgnsq], &repr[ixngn], &repr[ixbdet], &pivmin, &atau[j], negp,
		 auxp);
	for (i__ = 1; i__ <= 8; ++i__) {
	    gma = gmarep + (auxn[i__ - 1] + auxp[i__ - 1] - atau[j]);
	    xi = negn[i__ - 1] + negp[i__ - 1] << 1;
	    if (gma < 0.) {
		xi += 2;
	    } else if (gma == 0.) {
		++xi;
	    }
	    axi[j] = xi;
	    ++j;
/* L90015: */
	}
/* L90016: */
	goto L90013;
L90014:
	;
    }
    if (FALSE_) {
L90017:
	if (j + 3 > *m) {
	    goto L90018;
	}
	ibb = 1;
	dlaxrm_stat4__(n, &k, &c__1, &repr[ixg], &ibb, &repi[jylbbk], &repr[
		ixgnsq], &repr[ixngn], &repr[ixbdet], &pivmin, &atau[j], negn,
		 auxn);
	ibb = nb + 1;
	dlaxrm_stat4__(n, &k, &c_n1, &repr[ixg], &ibb, &repi[jylbbk], &repr[
		ixgnsq], &repr[ixngn], &repr[ixbdet], &pivmin, &atau[j], negp,
		 auxp);
	for (i__ = 1; i__ <= 4; ++i__) {
	    gma = gmarep + (auxn[i__ - 1] + auxp[i__ - 1] - atau[j]);
	    xi = negn[i__ - 1] + negp[i__ - 1] << 1;
	    if (gma < 0.) {
		xi += 2;
	    } else if (gma == 0.) {
		++xi;
	    }
	    axi[j] = xi;
	    ++j;
/* L90019: */
	}
/* L90020: */
	goto L90017;
L90018:
	;
    }
    if (TRUE_) {
L90021:
	if (j + 1 > *m) {
	    goto L90022;
	}
	ibb = 1;
	dlaxrm_stat2__(n, &k, &c__1, &repr[ixg], &ibb, &repi[jylbbk], &repr[
		ixgnsq], &repr[ixngn], &repr[ixbdet], &pivmin, &atau[j], negn,
		 auxn);
	ibb = nb + 1;
	dlaxrm_stat2__(n, &k, &c_n1, &repr[ixg], &ibb, &repi[jylbbk], &repr[
		ixgnsq], &repr[ixngn], &repr[ixbdet], &pivmin, &atau[j], negp,
		 auxp);
	for (i__ = 1; i__ <= 2; ++i__) {
	    gma = gmarep + (auxn[i__ - 1] + auxp[i__ - 1] - atau[j]);
	    xi = negn[i__ - 1] + negp[i__ - 1] << 1;
	    if (gma < 0.) {
		xi += 2;
	    } else if (gma == 0.) {
		++xi;
	    }
	    axi[j] = xi;
	    ++j;
/* L90023: */
	}
/* L90024: */
	goto L90021;
L90022:
	;
    }
    if (! xmrstats_1.xstealthmode) {
	xmrstats_1.xnumfn = xmrstats_1.xnumfn + j - 1;
    }
/*     Do all remaining with the standard single negcount */
L90025:
    if (j > *m) {
	goto L90026;
    }
    axi[j] = dlaxrn_(n, &repr[1], &repi[1], &atau[j]);
    ++j;
    goto L90025;
L90026:
    return 0;
} /* dlaxrm_ */

