/* dstexr_f77.f -- translated by f2c (version 20240504).
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

static doublereal c_b16 = .001;

/* Subroutine */ int dstexr_(integer *n, doublereal *d__, doublereal *e, 
	integer *wil, integer *wiu, doublereal *w, doublereal *z__, integer *
	ldz, integer *isuppz, doublereal *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, j, m, ib;
    extern /* Subroutine */ int wsreq_xrv__(integer *, integer *, integer *);
    integer za, ze, ixd, ixe;
    doublereal eps;
    integer iys, bbeg, bend, blen, bwil, bwiu;
    doublereal scale;
    char emode[1];
    doublereal bgerl, bgeru, bsepl;
    integer iinfo, ixf77a;
    doublereal bsepu, invsc;
    integer ixrep, iyrep;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal bsdiam;
    extern doublereal dsecnd_(void);
    extern /* Subroutine */ int dlaxra_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, integer *);
    integer nblcks, iyblck;
    doublereal bshift;
    extern /* Subroutine */ int dlaxre_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    char *, doublereal *, integer *, ftnlen);
    logical gotall;
    integer iyewae;
    extern /* Subroutine */ int dlaxri_(integer *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), dlaxro_(
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *);
    integer liwrem;
    extern /* Subroutine */ int dlaxrv_(integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    integer *);
    integer ixgers;
    doublereal asptol;
    integer iywind, lrwrem, wsreqi, ixewlu, ixvsep, wsreqr, iywork, ixwork;
    doublereal tstamp0, tstamp1;
    integer wsi_xrv__, wsr_xrv__;

/*     IMPLICIT NONE */





/*  Purpose */
/*  ======= */

/*    DSTEXR computes selected eigenpairs with indices wil:wiu for a symmetric */
/*    tridiagonal matrix T. The indexing of eigenpairs is with respect to an */
/*    ascending order of the eigenvalues, and the delivered eigenvalues will be */
/*    ordered ascendingly. */

/*    The computed eigenpairs will be consistent in the following sense: If the */
/*    routine is called twice with index sets wil1:wiu1 and wil2:wiu2, and those */
/*    index sets are non-overlapping, say wiu1 < wil2, then the results will obey */
/*    (1)  The computed eigenvalues from the first call are all <= the computed */
/*         eigenvalues from the second call. */
/*    (2)  The computed vectors from both calls will be numerically orthogonal */
/*          to each other. */

/*    Esp. (2) is the key to an easy 'naive' parallelization approach, where */
/*    each processor computes one about equally sized chunk of the desired */
/*    eigenpairs separately from all other processors. */
/*    The benefit of (1) is that the computed eigenpairs won't need to be */
/*    sorted across processors afterwards, that is, the whole computation can */
/*    be performed without *any* communication. */

/*  Arguments */
/*  ========= */

/*  D       (input) DOUBLE PRECISION ARRAY, dimension (N). */
/*          The diagonal entries of T. */

/*  E       (input) DOUBLE PRECISION ARRAY, dimension (N-1). */
/*          The offdiagonal entries of T. */

/*  WIL     (input) INTEGER */
/*  WIU     (input) INTEGER */
/*          The range of wanted eigenpairs, 1 <= WIL <= WIU <= N. */

/*  W       (output) DOUBLE PRECISION array, dimension ( WIL:WIU ) */
/*          The computed eigenvalues WIL:WIU. */
/*          The i'th entry is only set if the eigenpair could be computed, */
/*          as indicated by ISUPPZ. */

/*  Z       (output) DOUBLE PRECISION array, dimension ( LDZ, WIL:WIU ) */
/*          The computed orthonormal eigenvectors WIL:WIU. */
/*          The i'th column is only set if the eigenpair could be */
/*          computed, as indicated by ISUPPZ. */

/*  LDZ     (input) INTEGER */
/*          Leading dimension of array Z. */

/*  ISUPPZ  (output) INTEGER array, dimension (2*WIL-1:2*WIU) */
/*          If the i'th vector was succsessfully computed (cf INFO), where */
/*          WIL <= i <= WIU, then the vector Z(i) has nonzero entries */
/*          only within ISUPPZ(2*i-1):ISUPPZ(2*i), which is a subset */
/*          of 1:N. */
/*          Failure to compute the vector is indicated by */
/*          ISUPPZ(2*i-1) = 0, ISUPPZ(2*i) < 0, the precise value of */
/*          the latter may give information about the kind of failure: */
/*          -1:  depth limit exceeded */
/*          -2:  no rep found in DLAXRF */
/*          -3:  could not find a root rep */

/*  RWORK   (workspace) DOUBLE PRECISION, dimension (LRWORK) */

/*  LRWORK  (input/output) INTEGER */
/*          Dimension of the array RWORK. Set to -1 for a workspace */
/*          query, then the routine only computes the required real */
/*          workspace, sets LRWORK to this value and returns. */
/*            For a workspace query only N needs to be set, no other */
/*          argument is referenced. You can do a real and integer ws */
/*          query at once. */

/*  IWORK   (workspace) INTEGER, dimension (LIWORK) */

/*  LIWORK  (input/output) INTEGER */
/*          Dimension of the array IWORK. Set to -1 for a workspace */
/*          query, then the routine only computes the required integer */
/*          workspace, sets LIWORK to this value and returns. */
/*            For a workspace query only N needs to be set, no other */
/*          argument is referenced. You can do a real and integer ws */
/*          query at once. */

/*  INFO    (output) INTEGER */
/*          Indicate succes of the computation. For INFO=0 or INFO=1, the */
/*          result data structures are set, for other values this need not */
/*          be the case. */
/*          < 0 : Some argument had an illegal value */
/*            0 : Everything ok, all desired eigenpairs were computed. */
/*            1 : At least some eigenpairs could not be computed, refer */
/*                to ISUPPZ to see which. */
/*          > 1 : Some other error. */

/*  ====================================================================== */
/* extract -b parameters.inc.f procname=stexr */

/*     .. Parameters .. */

/*     For blocks of size <= QRDIM, the QR-Algorithm is invoked. */
/*     Activate the use of dqds to compute root-eigenvalues. Regardless */
/*     of this setting DQDS is only empleyed if all ews belonging to */
/*     a block are wanted. */

/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRA( */
/*    $             N, D, E, S, SCALE, ASPTOL, NBLCKS, ABINDS */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)     ::  N */

/*     DOUBLE PRECISION, INTENT(INOUT)  ::  D(N), E(N) */

/*     INTEGER,          INTENT(OUT)    ::  NBLCKS */
/*     INTEGER,          INTENT(OUT)    ::  S(N), ABINDS(2*N) */
/*     DOUBLE PRECISION, INTENT(OUT)    ::  ASPTOL, SCALE */

/*     END SUBROUTINE DLAXRA */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRE( */
/*    $             N, D, E, GL, GU, ABSERR, GAPTOL, */
/*    $             REPR, REPI, TAU, */
/*    $             EWL_AE, EWL_LU, EMODE, */
/*    $             RWORK, */
/*    $             INFO */
/*    $           ) */
/*     IMPLICIT NONE */

/*     CHARACTER,        INTENT(IN)     ::  EMODE */
/*     INTEGER,          INTENT(IN)     ::  N */
/*     DOUBLE PRECISION, INTENT(IN)     ::  GL, GU, ABSERR, GAPTOL */
/*     DOUBLE PRECISION, INTENT(IN)     ::  D(N) */

/*     DOUBLE PRECISION, INTENT(INOUT)  ::  E( N ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( 6*N ) */

/*     INTEGER,          INTENT(OUT)    ::  INFO */
/*     INTEGER,          INTENT(OUT)    ::  REPI( (6+N+N/2) ) */
/*     INTEGER,          INTENT(OUT)    ::  EWL_AE(1:2*N) */
/*     DOUBLE PRECISION, INTENT(OUT)    ::  TAU */
/*     DOUBLE PRECISION, INTENT(OUT)    ::  REPR( (4*N+3) ) */
/*     DOUBLE PRECISION, INTENT(OUT)    ::  EWL_LU(1:2*N) */

/*     END SUBROUTINE DLAXRE */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRI( */
/*    $             N, D, E, WIL, WIU, NBLCKS, ABINDS, */
/*    $             ABGERS, ABWIND, ABVSEP, */
/*    $             RWORK, IWORK */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)     ::  N, WIL, WIU, NBLCKS */
/*     INTEGER,          INTENT(IN)     ::  ABINDS(2*NBLCKS) */
/*     DOUBLE PRECISION, INTENT(IN)     ::  D(N), E(N) */

/*     INTEGER,          INTENT(INOUT)  ::  IWORK( 5 * NBLCKS ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( N ) */

/*     INTEGER,          INTENT(OUT)    ::  ABWIND(2*NBLCKS) */
/*     DOUBLE PRECISION, INTENT(OUT)    ::  ABGERS(2*NBLCKS) */
/*     DOUBLE PRECISION, INTENT(OUT)    ::  ABVSEP(2*NBLCKS) */

/*     END SUBROUTINE DLAXRI */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE DLAXRO( N, M, W, Z, LDZ, ISUPPZ, REVORD ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)     ::  N, M, LDZ */

/*     INTEGER,          INTENT(INOUT)  ::  ISUPPZ(2*M) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  W(M), Z(LDZ,M) */

/*     INTEGER,          INTENT(OUT)    ::  REVORD(M) */

/*     END SUBROUTINE DLAXRO */
/*     END INTERFACE */
/*     INTERFACE */
/*     SUBROUTINE WSREQ_XRV(N, REQR, REQI) */
/*     IMPLICIT NONE */
/*     INTEGER, INTENT(IN)   ::  N */
/*     INTEGER, INTENT(OUT)  ::  REQR, REQI */
/*     END SUBROUTINE WSREQ_XRV */

/* *********************************************************************** */

/*     SUBROUTINE DLAXRV( */
/*    $             N, E, ROOTR, ROOTI, EWL_AE, EWL_LU, */
/*    $             WIL, WIU, SPDIAM, GAPTOL, */
/*    $             W, Z, LDZ, ISUPPZ, */
/*    $             RWORK, LRWORK, IWORK, LIWORK, INFO */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, WIL, WIU, LDZ */
/*     INTEGER,          INTENT(IN)  ::  ROOTI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  SPDIAM, GAPTOL */
/*     DOUBLE PRECISION, INTENT(IN)  ::  E( N-1 ), ROOTR(4*N+3) */

/*     INTEGER,          INTENT(INOUT)  ::  LRWORK, LIWORK */
/*     INTEGER,          INTENT(INOUT)  ::  EWL_AE( 1:2*N ) */
/*     INTEGER,          INTENT(INOUT)  ::  IWORK( LIWORK ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU( 1:2*N ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( LRWORK ) */

/*     INTEGER,          INTENT(OUT)  ::  INFO */
/*     INTEGER,          INTENT(OUT)  ::  ISUPPZ( 2*WIL-1 : 2*WIU ) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  Z( LDZ, WIL : WIU ) */
/*     DOUBLE PRECISION, INTENT(OUT)  ::  W( WIL : WIU ) */
/*     END SUBROUTINE DLAXRV */
/*     END INTERFACE */

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

/*     =================================================================== */
/*     =                     Check input parameters                      = */
/*     =================================================================== */
/*     ---------------------- */
/*      Workspace Management */
/*     ---------------------- */
    /* Parameter adjustments */
    --e;
    --d__;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --rwork;
    --iwork;

    /* Function Body */
    ixwork = 1;
    iywork = 1;
/*     Copies of the matrix data to modify */
    ixd = ixwork;
    ixe = ixwork + *n;
    ixwork += *n << 1;
/*     Outer scaling by DLAXRA to get offdiagonals positive */
    iys = iywork;
    iywork += *n;
/*     Begin and end of blocks in the splitted matrix */
/*     ( we assume here N as upper bound for number of blocks, however */
/*       then we would not call the kernel, so this could be optimized ) */
    iyblck = iywork;
    iywork += *n << 1;
/*     Index ranges of desired eigenvalues for each block */
    iywind = iywork;
    iywork += *n << 1;
/*     Union of Gershgorin Discs and Value Separators for each block */
    ixgers = ixwork;
    ixvsep = ixwork + (*n << 1);
    ixwork += *n << 2;
/*     EW-List for the current block */
    iyewae = iywork;
    ixewlu = ixwork;
    iywork += *n << 1;
    ixwork += *n << 1;
/*     Root Representation for the current block */
    iyrep = iywork;
    ixrep = ixwork;
    iywork += *n + 6 + *n / 2;
    ixwork += (*n << 2) + 3;
    lrwrem = *lrwork - ixwork + 1;
    liwrem = *liwork - iywork + 1;

/*     -------------------------- */
/*      Handle workspace queries */
/*     -------------------------- */

    wsreq_xrv__(n, &wsr_xrv__, &wsi_xrv__);
    wsreqr = ixwork - 1 + wsr_xrv__;
    wsreqi = iywork - 1 + wsi_xrv__;
    if (*lrwork == -1 || *liwork == -1) {
	if (*lrwork == -1) {
	    *lrwork = wsreqr;
	}
	if (*liwork == -1) {
	    *liwork = wsreqi;
	}
	return 0;
    }
    tstamp0 = dsecnd_();
    eps = dlamch_("Epsilon", (ftnlen)7);
/*     Copy D and E so that we can modify them */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rwork[ixd + i__ - 1] = d__[i__];
/* L90001: */
    }
/* L90002: */
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rwork[ixe + i__ - 1] = e[i__];
/* L90003: */
    }
/* L90004: */
/*     =================================================================== */
/*     =                          Splitting                              = */
/*     =================================================================== */

/*     Split the matrix into irreducible blocks, record where they begin */
/*     and end, and scale them individually into proper numerical range. */

    dlaxra_(n, &rwork[ixd], &rwork[ixe], &iwork[iys], &scale, &asptol, &
	    nblcks, &iwork[iyblck]);
    invsc = 1. / scale;
/*     ======================== */
/*      Determine Index ranges */
/*     ======================== */
    dlaxri_(n, &rwork[ixd], &rwork[ixe], wil, wiu, &nblcks, &iwork[iyblck], &
	    rwork[ixgers], &iwork[iywind], &rwork[ixvsep], &rwork[ixwork], &
	    iwork[iywork]);
/*     =================================================================== */
/*     =                       Build the Vectors                         = */
/*     =================================================================== */
    tstamp1 = dsecnd_();
    xmrstats_1.xtime1 += tstamp1 - tstamp0;
    tstamp0 = tstamp1;
    xmrstats_1.xnblcks += nblcks;
    i__1 = *wiu - *wil + 1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *ldz;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__ + j * z_dim1] = 0.;
/* L90007: */
	}
/* L90008: */
/* L90005: */
    }
/* L90006: */
    i__1 = *wiu - *wil + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = 0.;
/* L90009: */
    }
/* L90010: */
    m = *wiu - *wil + 1;
    gotall = TRUE_;
    ze = *wil - 1;
    i__1 = nblcks;
    for (ib = 1; ib <= i__1; ++ib) {
	bbeg = iwork[iyblck - 1 + (ib << 1) - 1];
	bend = iwork[iyblck - 1 + (ib << 1)];
	blen = bend - bbeg + 1;
	bwil = iwork[iywind - 1 + (ib << 1) - 1];
	bwiu = iwork[iywind - 1 + (ib << 1)];
	bgerl = rwork[ixgers - 1 + (ib << 1) - 1];
	bgeru = rwork[ixgers - 1 + (ib << 1)];
	bsdiam = bgeru - bgerl;
	bsepl = rwork[ixvsep - 1 + (ib << 1) - 1];
	bsepu = rwork[ixvsep - 1 + (ib << 1)];
	if (bwil > bwiu) {
	    goto L90011;
	}
	za = ze + 1;
	ze = za + (bwiu - bwil);
	if (blen == 1) {
	    z__[bbeg + za * z_dim1] = 1.;
	    w[za] = rwork[ixd - 1 + bbeg];
	    isuppz[(za << 1) - 1] = 1;
	    isuppz[za * 2] = 1;
/* These special handlers are not implemented yet, since esp for debugging */
/* and testing we want to call the kernel on small matrices. */
/*         ELSEIF( BLEN .EQ. 2 )THEN */

/*            call special 2x2 handler */

/*         ELSEIF( BLEN .LE. QRDIM )THEN */

/*            call QR */

	} else {
/*           Do MRRR */
/*           Build the representation and init eigenvalues */
/*           The decision to compute the eigenvalues with dqds or not */
/*           may not depend on the local index range of wanted eigen- */
/*           pairs, otherwise we lose constistency. The sole exception */
/*           is if all ews in the block are wanted (even if for the full */
/*           matrix only a subset is desired). */
	    if (bwil == 1 && bwiu == blen) {
		*(unsigned char *)emode = 'd';
	    } else {
		*(unsigned char *)emode = 'o';
	    }
	    dlaxre_(&blen, &rwork[ixd - 1 + bbeg], &rwork[ixe - 1 + bbeg], &
		    bgerl, &bgeru, &asptol, &c_b16, &rwork[ixrep], &iwork[
		    iyrep], &bshift, &iwork[iyewae], &rwork[ixewlu], emode, &
		    rwork[ixwork], &iinfo, (ftnlen)1);
	    if (iinfo != 0) {
		i__2 = (ze << 1) - 1;
		for (ixf77a = (za << 1) - 1; ixf77a <= i__2; ixf77a += 2) {
		    isuppz[ixf77a] = 0;
/* L90013: */
		}
/* L90014: */
		i__2 = ze << 1;
		for (ixf77a = za << 1; ixf77a <= i__2; ixf77a += 2) {
		    isuppz[ixf77a] = -3;
/* L90015: */
		}
/* L90016: */
		goto L90011;
	    }
/*           Compute the vectors */
	    dlaxrv_(&blen, &rwork[ixe - 1 + bbeg], &rwork[ixrep], &iwork[
		    iyrep], &iwork[iyewae], &rwork[ixewlu], &bwil, &bwiu, &
		    bsdiam, &c_b16, &w[za], &z__[bbeg + za * z_dim1], ldz, &
		    isuppz[(za << 1) - 1], &rwork[ixwork], &lrwrem, &iwork[
		    iywork], &liwrem, &iinfo);
	    gotall = iinfo == 0;
	    if (iinfo != 0 && iinfo != 1) {
		*info = 10;
		return 0;
	    }
/*           Undo the shift from DLAXRE */
	    i__2 = ze;
	    for (j = za; j <= i__2; ++j) {
		w[j] += bshift;
/* L90017: */
	    }
/* L90018: */
	}

/*        ----------------------------------- */
/*         Adjust results to original matrix */
/*        ----------------------------------- */

/*        Eigenvalues: */
/*        - Cap eigenvalues by the separators */
/*        - Undo the scaling from DLAXRA */

/*        Eigenvectors: */
/*        - Undo the outer scaling by signature matrix in DLAXRA */
/*        - Adjust support to full matrix */
	i__2 = ze;
	for (j = za; j <= i__2; ++j) {
	    if (w[j] < bsepl) {
		w[j] = bsepl;
	    } else if (w[j] > bsepu) {
		w[j] = bsepu;
	    }
	    w[j] *= invsc;

	    if (isuppz[(j << 1) - 1] > 0) {
		isuppz[(j << 1) - 1] = isuppz[(j << 1) - 1] + bbeg - 1;
		isuppz[j * 2] = isuppz[j * 2] + bbeg - 1;
	    }
/* L90019: */
	}
/* L90020: */
	i__2 = bend;
	for (i__ = bbeg; i__ <= i__2; ++i__) {
	    if (iwork[iys - 1 + i__] != 1) {
/*              negate row of Z */
		i__3 = ze;
		for (j = za; j <= i__3; ++j) {
		    z__[i__ + j * z_dim1] = -z__[i__ + j * z_dim1];
/* L90023: */
		}
/* L90024: */
	    }
/* L90021: */
	}
/* L90022: */
L90011:
	;
    }
/* L90012: */
    tstamp1 = dsecnd_();
    xmrstats_1.xtime2 += tstamp1 - tstamp0;
    tstamp0 = tstamp1;
/*     =================================================================== */
/*                Sort the eigenpairs into ascending order */
/*     =================================================================== */
    dlaxro_(n, &m, &w[1], &z__[z_offset], ldz, &isuppz[1], &iwork[1]);
    *info = 0;
    if (! gotall) {
	*info = 1;
    }

    tstamp1 = dsecnd_();
    xmrstats_1.xtime3 += tstamp1 - tstamp0;
    tstamp0 = tstamp1;
    return 0;
} /* dstexr_ */

