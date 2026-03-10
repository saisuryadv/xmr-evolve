/* dlaxrb_clssfy_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxrb_clssfy__(integer *n, doublereal *repr, integer *
	repi, integer *depth, doublereal *spdiam, integer *il, integer *iu, 
	doublereal *lgap, doublereal *ugap, integer *ewl_ae__, doublereal *
	ewl_lu__, integer *rginfo, integer *wil, integer *wiu, integer *gil, 
	integer *giu, doublereal *gaptol, doublereal *rwork, integer *iwork)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    integer i__, j, k;
    doublereal gapthresh, avgthresh, lb;
    integer hl, jl;
    doublereal ub;
    integer kl, hu, ju, ku;
    doublereal eps;
    integer jxg, klx, kux, ilen;
    doublereal prec;
    integer nrbl, iocr, nrbu, jgapl;
    logical doavg;
    integer irefl;
    logical seekl;
    integer jgapu, soffl, irefu;
    logical seeku;
    integer nnocr, jxrbl, soffu, jxngn, itmp77, jxrbu, twist, jx2bis;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal absmax, repelg, spread;
    extern /* Subroutine */ int dlaxrc_(integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    integer biscnt;
    doublereal abstol, avgtol;
    integer ntobis;
    doublereal reltol;
    integer ixwork, jxwork;

/*     IMPLICIT NONE */




/*  Purpose */
/*  ======= */

/*     Refine the eigenvalue approximation for the node [DEPTH,IL:IU] */
/*     in the ew-list and classify them into singletons and clusters. */

/*     Let JL:JU be the intersection of IL:IU and WIL:WIU. Then the */
/*     eigenvalues JL:JU are our main concern, and in principle we only */
/*     want to refine them. But to enable the construction of an ortho- */
/*     normal base , we have to stick to a conformal subtree. */
/*     In particular, the final classified bounds and gaps must be */
/*     independent of WIL:WIU. */

/*     It is helpful to think in terms of the full bisection tree be- */
/*     longing to each node in the representation tree: The tree of */
/*     intervals one gets by bisecting all initial bounds IL:IU to */
/*     full precision. */

/*     Let KL be maximal in IL:JL such that there is a relative gap */
/*     left of it or KL=IL, and analogously let KU be minimal in JU:IU */
/*     with a relative gap on its right or KU=IU. */
/*     Then we guarantee that the delivered bounds KL:KU, as well as */
/*     all gaps left, right and in-between, are the same as if the */
/*     routine were called with WIL=IL and WIU=IU. Note that the bounds */
/*     outside of KL:KU may be completely different. */

/*     The measure of relative gaps we use is |a-b| / absmax(a,b), */
/*     with one exception: relgap(0,0) = 1. This is because, if the */
/*     upper bound of one interval and the lower bound of the next */
/*     are both equal to zero, we can identify a gap regardless, since */
/*     only one can have zero after all (unreduced matrix). */

/*  Notes: */
/*  (1) Because the bounds may change after the gaps are fixed, the */
/*      quantity gap(gcbnds) / absmax( ub(i), lb(i+1) ) may become smaller */
/*      after that ( lb(i+1) can increas, for example ). */
/*      But the true relative gap is monotone in |ub-lb| (at least if */
/*      sign lb = sign ub), so that can only increase. */

/*  Preconditions: */
/*      IL < IU */
/*      GIL <= WIL <= WIU <= GIU */
/*      JL:JU is not empty */


/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The dimension of the matrix */

/*  REPR    (input) DOUBLE PRECISION array, dimension (4*N+3) */
/*  REPI    (input) INTEGER array, dimension (6+N+N/2) */
/*          Real and integer data used to represent the matrix */

/*  DEPTH   (input) INTEGER */
/*          Depth of the node in the representation tree */

/*  IL      (input) INTEGER */
/*  IU      (input) INTEGER */
/*          Index range belonging to the node */

/*  LGAP    (input/output) DOUBLE PRECISION */
/*  UGAP    (input/output) DOUBLE PRECISION */
/*  EWL_AE  (input/output) INTEGER array, dimension (1:2*N) */
/*  EWL_LU  (input/output) DOUBLE PRECISION array, dimension ( 1:2*N ) */
/*          They constitute the initial list of eigenvalue bounds */
/*          (ew-list) and will hold the refinements, see DLAXRL */
/*          for details. Need only be initialized properly. */

/*  RGINFO  (output) INTEGER array, dimension ( IL, IU-1 ) */
/*          RGINFO(i) holds a status flags to indicate the kind of */
/*          relative gap to the right of ew i, which can be */
/*             FULLGAP - if there is a full gap */
/*             NOFULL  - if there is definitely no full gap between them */
/*             UNKNOWN - if in the list ews i and i+1 still belong to */
/*                       the same interval */
/*          Corresponding constants are defined in gapinfo.inc.f. */
/*            It is guaranteed that not yet emerged gaps in the list */
/*          have info UNKNOWN, and that for all emerged gaps */
/*          the info is either NOFULL or FULLGAP. */

/*  WIL     (input) INTEGER */
/*  WIU     (input) INTEGER */
/*          The range of wanted eigenvalues. */
/*          Eigenvalues in the intersection of IL:IU and WIL:WIU will */
/*          be fully classified in any case, this intersection should */
/*          not be empty. */

/*  GIL     (input) INTEGER */
/*  GIU     (input) INTEGER */
/*          The global range of eigenvalues for which consistency is */
/*          desired, must be a superset of WIL:WIU. */
/*          Eigenvalues outiside of GIL:GIU will not be refined at all. */

/*  GAPTOL  (input) DOUBLE PRECISION */
/*          The gap tolerance, must be < 1 */

/*  ====================================================================== */
/* @extract -b parameters.inc.f procname=laxrb_clssfy */

/*     .. Parameters .. */

/*     Full gaps are also recognized if */
/*       absgap >= AVGAPFAC * SPDIAM / (N-1) */
/*     and */
/*       DEPTH <= MAXAVGAPDEPTH. */
/*     Set the latter to <0 to deactivate. */

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
    if (*depth <= 10) {
	jxg = 0;
	jxngn = (*n << 1) + 1;
	twist = repi[2];
	repelg = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__2 = repelg, d__3 = (d__1 = repr[jxg + i__], abs(d__1));
	    repelg = max(d__2,d__3);
/* L90001: */
	}
/* L90002: */
/* Computing MAX */
	d__3 = repelg, d__4 = (d__1 = repr[jxngn + twist - 1], abs(d__1)), 
		d__3 = max(d__3,d__4), d__4 = (d__2 = repr[jxngn + twist + 1],
		 abs(d__2));
	repelg = max(d__3,d__4);
/*        Scale with DEPTH to account for larger residual effects */
	avgtol = max(*spdiam,repelg) * (*depth * .3 / (*n - 1));
	doavg = TRUE_;
    } else {
	avgtol = 0.;
	doavg = FALSE_;
    }
    eps = dlamch_("Epsilon", (ftnlen)7);
    prec = dlamch_("Precision", (ftnlen)9);
    ilen = *iu - *il + 1;
/*     Overview of used index ranges: */
/*      JL:JU   are the indices of ews we want from this node */
/*              (intersection if IL:IU with WIL:WIU) */
/*      HL:HU   are the indices of ews from this node for which */
/*              consistency is desired, superset of JL:JU */
/*              (intersection of IL:IU with GIL:GIU) */
/*     NOTE: Neither of those need to be boundary indices of intervals. */
    jl = max(*il,*wil);
    ju = min(*iu,*wiu);
    hl = max(*il,*gil);
    hu = min(*iu,*giu);
    jxrbl = 1;
    jxrbu = *n + 1;
    jx2bis = (*n << 1) + 1;
    jxwork = (*n << 1) + 1 + ilen;
    ixwork = 1;
    i__1 = (*n << 1) + (*iu - *il + 1 << 2);
    for (itmp77 = 1; itmp77 <= i__1; ++itmp77) {
	iwork[itmp77] = 0;
/* L90003: */
    }
/* L90004: */
    i__1 = *iu - *il + 1;
    for (itmp77 = 1; itmp77 <= i__1; ++itmp77) {
	rwork[itmp77] = 0.;
/* L90005: */
    }
/* L90006: */
    biscnt = xmrstats_1.xnumfn;

/*   ===================================================================== */

/*     Refine all eigenvalues enough to classify them */
    abstol = avgtol * .5;
    reltol = *gaptol * .5;
    i__1 = ju;
    for (i__ = jl; i__ <= i__1; ++i__) {
	iwork[jx2bis + i__ - jl] = i__;
/* L90007: */
    }
/* L90008: */
    i__1 = ju - jl + 1;
    dlaxrc_(n, &repr[1], &repi[1], il, iu, lgap, ugap, &ewl_ae__[1], &
	    ewl_lu__[1], &i__1, &iwork[jx2bis], &reltol, &abstol, &rwork[
	    ixwork], &iwork[jxwork]);
    xmrstats_1.xnbis_class__ += xmrstats_1.xnumfn - biscnt;
    biscnt = xmrstats_1.xnumfn;
/*   ===================================================================== */

/*     Classify them */
/*     Note: This approach is an easy way to consistency, but means that */
/*     the gap tolerances are taken somewhat relaxed. With these settings, */
/*     we do guarantee that each classified gap has a relative width */
/*     exceeding GAPTOL or absolute width exceeding AVGTOL, but it may */
/*     be that gaps up to 2*GAPTOL or 2*ABSTOL, respectively, are not */
/*     recognized as such. */
    avgthresh = avgtol * 2.;
    gapthresh = *gaptol * 2.;
    if (! doavg) {
	avgthresh = *spdiam;
    }
    i__1 = *iu - *il;
    for (itmp77 = 1; itmp77 <= i__1; ++itmp77) {
	rginfo[itmp77] = -1;
/* L90009: */
    }
/* L90010: */
    i__ = jl;
L90011:
    j = ewl_ae__[i__ * 2];
    if (j >= ju) {
	goto L90012;
    }
    i__ = j + 1;
/*        Evaluate if there is a full gap in ]J,I[ */
    lb = ewl_lu__[(j << 1) - 1];
    ub = ewl_lu__[i__ * 2];
/*        Now LB is the lower bound of the interval containing ew J, */
/*        and UB is the upper bound of the interval with I. */
/* Computing MAX */
    d__1 = abs(lb), d__2 = abs(ub);
    absmax = max(d__1,d__2);
    spread = ub - lb;
/* Computing MIN */
    d__1 = avgthresh, d__2 = absmax * gapthresh;
    if (spread >= min(d__1,d__2)) {
/*           Yes, there must be a gap between J and I. */
	rginfo[j] = 2;
    } else {
/*           No, we don't recognize it. */
/*           However, due to the fact that the intervals were not refined */
/*           to full accuracy, there may actually be a gap here. But it */
/*           could not be a wide one, cf the note above. */
	rginfo[j] = -2;
    }
    j = ewl_ae__[i__ * 2];
    goto L90011;
L90012:
/*   ===================================================================== */

/*     Handle the fringes */
/*     We have to find the extents of clusters that intersect JL:JU */
/*     but are not contained in it. */
    kl = ewl_ae__[(jl << 1) - 1];
    ku = ewl_ae__[ju * 2];
    nrbl = 0;
    nrbu = 0;
    seekl = TRUE_;
    seeku = TRUE_;
    klx = -1;
    kux = -1;
    soffl = 1;
    soffu = 1;
L90013:
/*        Determine actual range we are looking in */
    if (seekl) {
L90015:
	if (nrbl == 0) {
	    klx = hl - 1;
	    goto L90016;
	}
	klx = ewl_ae__[iwork[jxrbl + nrbl - 1] * 2];
	if (klx < kl) {
	    lb = ewl_lu__[(klx << 1) - 1];
	    ub = ewl_lu__[kl * 2];
/* Computing MAX */
	    d__1 = abs(lb), d__2 = abs(ub);
	    absmax = max(d__1,d__2);
	    spread = ub - lb;
/* Computing MIN */
	    d__1 = avgthresh, d__2 = absmax * gapthresh;
	    if (spread >= min(d__1,d__2)) {
/*                    Yes, this still looks like it holds a gap */
		if (klx == kl - 1) {
		    rginfo[klx] = 2;
		    seekl = FALSE_;
		}
		goto L90016;
	    }
/* Computing MAX */
	    i__1 = hl, i__2 = ewl_ae__[(klx << 1) - 1];
	    kl = max(i__1,i__2);
	}
	--nrbl;
	goto L90015;
L90016:
	;
    }
    if (seeku) {
L90017:
	if (nrbu == 0) {
	    kux = hu + 1;
	    goto L90018;
	}
	kux = ewl_ae__[(iwork[jxrbu + nrbu - 1] << 1) - 1];
	if (kux > ku) {
	    lb = ewl_lu__[(ku << 1) - 1];
	    ub = ewl_lu__[kux * 2];
/* Computing MAX */
	    d__1 = abs(lb), d__2 = abs(ub);
	    absmax = max(d__1,d__2);
	    spread = ub - lb;
/* Computing MIN */
	    d__1 = avgthresh, d__2 = absmax * gapthresh;
	    if (spread >= min(d__1,d__2)) {
		if (ku == kux - 1) {
		    rginfo[ku] = 2;
		    seeku = FALSE_;
		}
		goto L90018;
	    }
/* Computing MIN */
	    i__1 = hu, i__2 = ewl_ae__[kux * 2];
	    ku = min(i__1,i__2);
	}
	--nrbu;
	goto L90017;
L90018:
	;
    }
    seekl = seekl && hl < kl;
    seeku = seeku && ku < hu;
    if (! seekl && ! seeku) {
	goto L90014;
    }
/*     We are looking for a gap between KLX AND KL, and between KU and KUX. */
/*     The intervals for KL and KU are refined already, and both are */
/*     boundary indices. Also, KLX is either HL-1 or is upper boundary */
/*     of an interval that was refined, and analogously for KUX. */
/*        Determine outside eigenvalues to refine */
    ntobis = 0;

    irefl = -1;
    if (seekl) {
	if (klx == hl - 1) {
	    irefl = kl - soffl;
	    soffl <<= 1;
	} else {
	    irefl = (klx + kl) / 2;
	}
	irefl = max(hl,irefl);
	iwork[jx2bis + ntobis] = irefl;
	++ntobis;
	iwork[jxrbl + nrbl] = irefl;
	++nrbl;
    }

    irefu = -1;
    if (seeku) {
	if (kux == hu + 1) {
	    irefu = ku + soffu;
	    soffu <<= 1;
	} else {
	    irefu = (ku + kux) / 2;
	}
	irefu = min(hu,irefu);
	iwork[jx2bis + ntobis] = irefu;
	++ntobis;
	iwork[jxrbu + nrbu] = irefu;
	++nrbu;
    }
    dlaxrc_(n, &repr[1], &repi[1], il, iu, lgap, ugap, &ewl_ae__[1], &
	    ewl_lu__[1], &ntobis, &iwork[jx2bis], &reltol, &abstol, &rwork[
	    ixwork], &iwork[jxwork]);
    goto L90013;
L90014:
/*     Now KL and KU are correctly placed. Two possible scenarios, for */
/*     KL these are */
/*       KL = HL   meaning there is no gap between KL and JL, or */
/*       HL < KL <= JL  with RGINFO(KL-1) = GI_FULLGAP */
    xmrstats_1.xnbis_class__ += xmrstats_1.xnumfn - biscnt;
    biscnt = xmrstats_1.xnumfn;
/*   ===================================================================== */
/*     For consistency we now have to reset inner bounds for clusters */
/*     which concern us (intersect JL:JU) without being owned (not a */
/*     subset of JL:JU). The reason is that we searched for the (from our */
/*     perspective) outside end of the cluster starting from JL:JU, so */
/*     so for another call with different JL:JU, but still intersecting */
/*     the same cluster, the bounds inside the cluster will be different. */
/*     But note that the two intervals constituting the 'outer' gaps on */
/*     each side will always come out identical. */
/*     As a first step we have to determine the outermost full gaps inside */
/*     of JL:JU. Note that there may be none. */
    k = ewl_ae__[jl * 2];
L90019:
    if (k >= ju) {
	k = ku;
	goto L90020;
    }
    if (rginfo[k] == 2) {
	goto L90020;
    }
    k = ewl_ae__[(k + 1) * 2];
    goto L90019;
L90020:
    jgapl = k;
    k = ewl_ae__[(ju << 1) - 1];
L90021:
    if (k <= jgapl) {
	k = kl;
	goto L90022;
    }
    if (rginfo[k - 1] == 2) {
	goto L90022;
    }
    k = ewl_ae__[(k - 1 << 1) - 1];
    goto L90021;
L90022:
    jgapu = k;
/*     If JGAPL=KU and JGAPU=KL, then there is no gap  within JL:JU, */
/*     meaning the ews we want all belong to one cluster (which we may */
/*     or may not own). */
    nnocr = 0;
    if (jgapl == ku && jgapu == kl) {
	if (kl < jl || ku > ju) {
	    iwork[1] = kl;
	    iwork[2] = ku;
	    nnocr = 1;
	}
    } else {
	if (kl < jl) {
	    ++nnocr;
	    iwork[1] = kl;
	    iwork[2] = jgapl;
	}
	if (ku > ju) {
	    ++nnocr;
	    iwork[(nnocr << 1) - 1] = jgapu;
	    iwork[nnocr * 2] = ku;
	}
    }
    i__1 = nnocr;
    for (iocr = 1; iocr <= i__1; ++iocr) {
	i__ = ewl_ae__[iwork[(iocr << 1) - 1] * 2];
	j = ewl_ae__[(iwork[iocr * 2] << 1) - 1];
	lb = ewl_lu__[i__ * 2];
	ub = ewl_lu__[(j << 1) - 1];
	if (i__ < j) {
	    rginfo[i__] = -2;
	    ++i__;
	    --j;
	    i__2 = j;
	    for (k = i__; k <= i__2; ++k) {
		ewl_ae__[(k << 1) - 1] = i__;
		ewl_ae__[k * 2] = j;
		ewl_lu__[(k << 1) - 1] = lb;
		ewl_lu__[k * 2] = ub;
/* L90025: */
	    }
/* L90026: */
	    i__2 = j - 1;
	    for (itmp77 = i__; itmp77 <= i__2; ++itmp77) {
		rginfo[itmp77] = -1;
/* L90027: */
	    }
/* L90028: */
	    rginfo[j] = -2;
	}
/* L90023: */
    }
/* L90024: */

/*   ===================================================================== */

    return 0;
} /* dlaxrb_clssfy__ */

