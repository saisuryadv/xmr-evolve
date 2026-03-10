/* dlaxrc_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxrc_(integer *n, doublereal *repr, integer *repi, 
	integer *il, integer *iu, doublereal *lgap, doublereal *ugap, integer 
	*ewl_ae__, doublereal *ewl_lu__, integer *numew, integer *ewinds, 
	doublereal *reltol, doublereal *abstol, doublereal *rwork, integer *
	iwork)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__, j, k, m;
    doublereal lb, ub, mid;
    integer ixf77a;
    doublereal width;
    extern /* Subroutine */ int dlaxrl_refine__(integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *);
    logical chkabs, chkrel;
    doublereal absmax;
    integer jxbind;
    extern /* Subroutine */ int dlaxrm_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *);
    integer jxflgs, ntobis, jxiner;

/*     IMPLICIT NONE */




/*  Purpose */
/*  ======= */

/*     Refine eigenvalues with indices given in EWINDS until they have */
/*     relative accuracy RTOL or absolute accuracy ABSTOL. */

/*  Arguments */
/*  ========= */

/*  NUMEW   (input/output) INTEGER,  <= IU-IL+1 */
/*          On entry, the number of eigenvalues to refine. */
/*          Overwritten on exit. */

/*  EWINDS  (input/output) INTEGER array, dimension (NUMEW) */
/*          On input, the indices of the eigenvalues to refine. */
/*          Need not be sorted and may even contain duplicates, the */
/*          no unnessessary bisection steps will be done. */

/*  RELTOL  (input) DOUBLE PRECISION, >= 0 */
/*  ABSTOL  (input) DOUBLE PRECISION, >= 0 */
/*          Tolerances for refinement. An interval is converged if */
/*           (1)  |ub - lb| <= absmax(lb,ub) * RELTOL */
/*                or RELTOL = 0 */
/*          and */
/*           (2)  |ub - lb| <= ABSTOL */
/*                or  ABSTOL = 0. */

/*          Thus by setting RELTOL or ABSTOL to zero, the corresponding */
/*          test can be deactivated. */
/*          In any case an interval is also converged if it so tight that */
/*          the midpoint is computed as one of lb or ub, which normally */
/*          means there is no fp-number between them (ub = next(lb)). */


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
/*     SUBROUTINE  DLAXRM( N, REPR, REPI, M, ATAU, AXI ) */
/*     IMPLICIT NONE */

/*     INTEGER,          INTENT(IN)  ::  N, M */
/*     INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2) */
/*     DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), ATAU(M) */

/*     INTEGER,          INTENT(OUT)  ::  AXI(M) */
/*     END SUBROUTINE DLAXRM */
/*     END INTERFACE */

/*     .. Constants .. */


/*     .. Locals .. */


/*     -- Executable Statements ----------------------------------------- */

/*     ------------ */
/*      Tolerances */
/*     ------------ */
    /* Parameter adjustments */
    --repi;
    --repr;
    --ewl_ae__;
    --ewl_lu__;
    --iwork;
    --rwork;
    --ewinds;

    /* Function Body */
    chkabs = *abstol != 0.;
    chkrel = *reltol != 0.;

/*     ----------------------- */
/*      Auxiliary data fields */
/*     ----------------------- */
/*      FLGS(IL:IU) */
/*        FLGS(I) is set to 1 to indicate that the interval with left */
/*        index I is already included in the current sweep */
/*      BIND(1:NUMEW)  left indices of intervals in the current sweep */
/*      INER(1:NUMEW)  inertias of samples in the current sweep */
/*      RWORK(1:NUMEW)  midpoints of intervals in the current sweep */

    jxflgs = 1 - *il;
    jxbind = jxflgs + *iu;
    jxiner = jxbind + *numew;

    m = *numew;

L90001:

/*        ----------------------- */
/*         Check for convergence */
/*        ----------------------- */

    i__1 = jxflgs + *iu;
    for (ixf77a = jxflgs + *il; ixf77a <= i__1; ++ixf77a) {
	iwork[ixf77a] = 0;
/* L90003: */
    }
/* L90004: */
    ntobis = 0;
    k = 1;
L90005:
    if (k > m) {
	goto L90006;
    }
    i__ = ewl_ae__[(ewinds[k] << 1) - 1];
    lb = ewl_lu__[(i__ << 1) - 1];
    ub = ewl_lu__[i__ * 2];
    mid = (lb + ub) * .5;
    width = ub - lb;
/* Computing MAX */
    d__1 = abs(lb), d__2 = abs(ub);
    absmax = max(d__1,d__2);
    if (iwork[jxflgs + i__] != 0) {
/*              corresponding interval was already processed */
	++k;
    } else {
	if (mid == lb || mid == ub || (! chkabs || width <= *abstol) && (! 
		chkrel || width <= absmax * *reltol)) {
/*                 this interval is converged */
	    ewinds[k] = ewinds[m];
	    --m;
	} else {
/*                 not converged, add midpoint for one bisection step */
	    ++ntobis;
	    iwork[jxbind + ntobis] = i__;
	    rwork[ntobis] = mid;
/*                 set flag that this interval was already processed */
	    iwork[jxflgs + i__] = 1;
	    ++k;
	}
    }
    goto L90005;
L90006:

/*        ----------- */
/*         Loop Exit */
/*        ----------- */

    if (m == 0) {
	goto L90002;
    }

/*        --------------------------------------- */
/*         Bisect all unconverged intervals once */
/*        --------------------------------------- */

    dlaxrm_(n, &repr[1], &repi[1], &ntobis, &rwork[1], &iwork[jxiner + 1]);
    i__1 = ntobis;
    for (k = 1; k <= i__1; ++k) {
	i__ = iwork[jxbind + k];
	j = ewl_ae__[i__ * 2];
	dlaxrl_refine__(il, iu, lgap, ugap, &ewl_ae__[1], &ewl_lu__[1], &i__, 
		&j, &rwork[k], &iwork[jxiner + k]);
/* L90007: */
    }
/* L90008: */
    goto L90001;
L90002:
    return 0;
} /* dlaxrc_ */

