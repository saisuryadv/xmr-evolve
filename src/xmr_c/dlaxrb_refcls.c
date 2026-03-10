/* dlaxrb_refcls_f77.f -- translated by f2c (version 20240504).
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

static integer c__2 = 2;
static doublereal c_b5 = 0.;

/* Subroutine */ int dlaxrb_refcls__(integer *n, doublereal *repr, integer *
	repi, integer *depth, integer *icbeg, integer *icend, doublereal *
	lgap, doublereal *ugap, integer *ewl_ae__, doublereal *ewl_lu__, 
	integer *rginfo, doublereal *rwork, integer *iwork)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j;
    doublereal lb, ub, eps, prec, raclb;
    integer iclen, jx2bis;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal absgap, absmax;
    extern /* Subroutine */ int dlaxrc_(integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    integer biscnt, ntobis;
    doublereal scltol;
    integer ixwork, jxwork;

/*     IMPLICIT NONE */




/*  Purpose */
/*  ======= */

/*     Refine the eigenvalue approximations for a cluster, and reveal */
/*     internal gaps within the cluster. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The dimension of the matrix */

/*  REPR    (input) DOUBLE PRECISION array, dimension (4*N+3) */
/*  REPI    (input) INTEGER array, dimension (6+N+N/2) */
/*          Real and integer data used to represent the matrix */

/*  ICBEG   (input) INTEGER */
/*  ICEND   (input) INTEGER */
/*          Index range constituting the (single) cluster, */
/*          1 <= ICBEG < ICEND <= N. */

/*  RGINFO  (intput/output) INTEGER array, dimension ( IL : IU-1 ) */
/*          RGINFO(i) holds a status flags to indicate the kind of */
/*          relative gap to the right of ew i. Should have been set */
/*          by DLAXRB. Corresponding constants are defined in */
/*          gapinfo.inc.f. */
/*            Since we are supposed to be given a single cluster, the */
/*          entries should be UNKNOWN for gaps still hidden in an */
/*          interval, and not GI_FULLGAP otherwise. */
/*            Upon exit, some of the gaps that were previously set */
/*          as GI_NOFULL may be changed to GI_INTGAP or GI_NOGAP. */

/*  ====================================================================== */
/* @extract -b parameters.inc.f procname=laxrb_refcls */

/*     .. Parameters .. */

/*     Gap Tolerance for subclusters is 2**LOGSCLTOLFAC * SQRT(PREC) */
/*     Cluster boundaries will be refined to relative accuracy */
/*     RACLBFAC*N*EPS */

/*     .. Declarations .. */

/*     INTERFACE */
/*     SUBROUTINE DLAXRC( */
/*    $            N, REPR, REPI, */
/*    $            IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, */
/*    $            NUMEW, EWINDS, */
/*    $            RELTOL, ABSTOL, RWORK, IWORK */
/*    $           ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)     ::  N, IL, IU, NUMEW */
/*     INTEGER,          INTENT(IN)     ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)     ::  RELTOL, ABSTOL */
/*     DOUBLE PRECISION, INTENT(IN)     ::  REPR(4*N+3) */

/*     INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU) */
/*     INTEGER,          INTENT(INOUT)  ::  EWINDS(NUMEW) */
/*     INTEGER,          INTENT(INOUT)  ::  IWORK( 2*NUMEW + IU-IL+1 ) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU) */
/*     DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( NUMEW ) */

/*     END SUBROUTINE DLAXRC */
/*     END INTERFACE */

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

/*  ===== Executable Statements ========================================== */

    /* Parameter adjustments */
    --repi;
    --repr;
    --iwork;
    --rwork;
    --ewl_ae__;
    --ewl_lu__;
    --rginfo;

    /* Function Body */
    eps = dlamch_("Epsilon", (ftnlen)7);
    prec = dlamch_("Precision", (ftnlen)9);
    iclen = *icend - *icbeg + 1;
    scltol = sqrt(prec);
    raclb = (*n << 2) * eps;
    jx2bis = 1;
    jxwork = iclen + 1;
    ixwork = 1;
    biscnt = xmrstats_1.xnumfn;

/*   ===================================================================== */

/*     Refine outer bounds */
    iwork[jx2bis] = *icbeg;
    iwork[jx2bis + 1] = *icend;
    dlaxrc_(n, &repr[1], &repi[1], icbeg, icend, lgap, ugap, &ewl_ae__[1], &
	    ewl_lu__[1], &c__2, &iwork[jx2bis], &raclb, &c_b5, &rwork[ixwork],
	     &iwork[jxwork]);
    xmrstats_1.xnbis_clb__ += xmrstats_1.xnumfn - biscnt;
    biscnt = xmrstats_1.xnumfn;

/*   ===================================================================== */

/*     Refine all eigenvalues enough to classify them */
/* Extra in-cluster refinement deactivated for now. Effectively disables */
/* the use of internal gaps for the subset case (efficiency concerns), but */
/* note that for the full spectrum on level one, where the ews are refined */
/* to full acc using dqds, internal gaps will be classified and used */
/* anyway. */
/*      RELTOL = HALF * SCLTOL */

/*      IF( .FALSE. )THEN */

/*         DO I = ICBEG, ICEND */
/*            IWORK(JX2BIS + I-ICBEG) = I */
/*         ENDDO */

/*         CALL DLAXRC( N, REPR, REPI, */
/*     $                ICBEG, ICEND, LGAP, UGAP, EWL_AE, EWL_LU, */
/*     $                ICLEN, IWORK(JX2BIS), RELTOL, ZERO, */
/*     $                RWORK(IXWORK), IWORK(JXWORK) */
/*     $              ) */

/* @LOGGING on */
/*         WRITE(FIDBIS,*) */
/*     $      '  Did ', (XNUMFN-BISCNT)*INVLEN, */
/*     $      '*iclen steps for classification.' */
/*         WRITE(FIDBIS,*) ' ' */
/* @LOGGING ! */
/* @STATS on */
/*         XNBIS_CLASS = XNBIS_CLASS + (XNUMFN - BISCNT) */
/*         BISCNT = XNUMFN */
/* @STATS ! */
/*      ENDIF */
/*   ===================================================================== */

/*     Classify them */
    i__ = *icbeg;
L90001:
    j = ewl_ae__[i__ * 2];
    if (j == *icend) {
	goto L90002;
    }
    i__ = j + 1;
/*        Evaluate ]J,I[ */
    if (rginfo[j] == -2 || rginfo[j] == -1) {
	ub = ewl_lu__[j * 2];
	lb = ewl_lu__[(i__ << 1) - 1];
/* Computing MAX */
	d__1 = abs(lb), d__2 = abs(ub);
	absmax = max(d__1,d__2);
	absgap = lb - ub;
	if (absgap >= absmax * scltol) {
	    rginfo[j] = 1;
	} else {
	    rginfo[j] = 0;
	}
    }
    j = ewl_ae__[i__ * 2];
    goto L90001;
L90002:

/*   ===================================================================== */

/*     Refine boundaries of internal gaps */

    ntobis = 0;
    j = *icbeg - 1;
L90003:
    i__ = j + 1;
L90005:
    j = ewl_ae__[(j + 1) * 2];
    if (j == *icend) {
	goto L90006;
    }
    if (rginfo[j] == 1) {
	goto L90006;
    }
    goto L90005;
L90006:
    if (*icbeg < i__) {
	iwork[jx2bis + ntobis] = i__;
	++ntobis;
    }
    if (i__ < j && j < *icend) {
	iwork[jx2bis + ntobis] = j;
	++ntobis;
    }
    if (j == *icend) {
	goto L90004;
    }
    goto L90003;
L90004:
    dlaxrc_(n, &repr[1], &repi[1], icbeg, icend, lgap, ugap, &ewl_ae__[1], &
	    ewl_lu__[1], &ntobis, &iwork[jx2bis], &raclb, &c_b5, &rwork[
	    ixwork], &iwork[jxwork]);
    xmrstats_1.xnbis_clb__ += xmrstats_1.xnumfn - biscnt;
    return 0;
} /* dlaxrb_refcls__ */

