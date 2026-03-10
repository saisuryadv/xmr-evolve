/* dlaxre_initewldqds_f77.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dlaxre_initewldqds__(integer *n, doublereal *repr, 
	integer *repi, integer *il, integer *iu, doublereal *qdvals, integer *
	wil, integer *wiu, integer *ewl_ae__, doublereal *ewl_lu__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);

    /* Local variables */
    integer i__, j, k;
    doublereal gap, prec;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal absmax, reloff, septol;

/*     IMPLICIT NONE */




/*  Purpose */
/*  ======= */

/*    Init the ewlist for eigenvalues IL:IU from the bounds given */
/*    in QDVALS. Those are assumed to stem from DQDS, that is, be */
/*    aproximations with high relative accuracy. */

/*    The bounds will be relaxed a bit, but not verified. */

/* !   NOTE: At the moment the consistency requirements in DSTEXR allow */
/* !   only to use dqds for blocks where all ews are wanted. Thus this */
/* !   routine is only used with WIL=1, WIU=N. Future versions without */
/* !   consistency guarantees might use the partial feature here though. */

/*  Arguments */
/*  ========= */

/*  QDVALS  (input), DOUBLE PRECISION array, dimension(IL:IU) */
/*          The eigenvalue approximations, in ascending order. */
/*          We do not require them to have the same sign or be distinct. */
/*          However, the list should not contain more than one zero. */

/*  ====================================================================== */


/*     .. Constants .. */


/*     .. Locals .. */


/*  ===== Executable Statements ========================================== */

    /* Parameter adjustments */
    --repi;
    --repr;
    --qdvals;
    --ewl_ae__;
    --ewl_lu__;

    /* Function Body */
    prec = dlamch_("Precision", (ftnlen)9);
    septol = sqrt(prec);
    reloff = log((doublereal) (*n)) * prec;
/*     Group the values, and then relax the outer bounds of those */
/*     groups just a bit to avoid unnecessary backing off in DLAXRF. */
/*     [ Note: Carful with the backing off only to to at group boundaries, */
/*       otherwise the monotonicity of bounds in the EWL may be destroyed] */
    i__ = *il;
L90001:
    j = i__;
L90003:
    if (j == *iu) {
	goto L90004;
    }
    gap = qdvals[j + 1] - qdvals[j];
/* Computing MAX */
    d__3 = (d__1 = qdvals[j], abs(d__1)), d__4 = (d__2 = qdvals[j + 1], abs(
	    d__2));
    absmax = max(d__3,d__4);
    if (gap > absmax * septol) {
	goto L90004;
    }
    ++j;
    goto L90003;
L90004:
    i__1 = j;
    for (k = i__; k <= i__1; ++k) {
	ewl_ae__[(k << 1) - 1] = k;
	ewl_ae__[k * 2] = k;
	ewl_lu__[(k << 1) - 1] = qdvals[k];
	ewl_lu__[k * 2] = qdvals[k];
/* L90005: */
    }
/* L90006: */
/*        Note: For singletons the midpoint is taken, so the following */
/*        has no effect there, the starting value for RQI will still be */
/*        QDVALS(I). */
    ewl_lu__[(i__ << 1) - 1] = qdvals[i__] - reloff * (d__1 = qdvals[i__], 
	    abs(d__1));
    ewl_lu__[j * 2] = qdvals[j] + reloff * (d__1 = qdvals[j], abs(d__1));
    if (j == *iu) {
	goto L90002;
    }
    i__ = j + 1;
    goto L90001;
L90002:
/*      ARXFAC(1) = LOG( DBLE(N) ) * (4 * EPS) */
/*      ARXFAC(2) =   4 * N * EPS */
/*      ARXFAC(3) =  32 * N * EPS */
/*      ARXFAC(4) = 400 * N * EPS */

/*      I = IL */
/*      outer: DO */
/*         J = I */
/*         findgroup: DO */
/*            IF( J .EQ. IU )THEN */
/*               EXIT */
/*            ENDIF */
/*            GAP    = QDVALS(J+1)-QDVALS(J) */
/*            ABSMAX = MAX( ABS(QDVALS(J)), ABS(QDVALS(J+1)) ) */
/*            IF( GAP .GT. ABSMAX*GAPTOL )THEN */
/*               EXIT findgroup */
/*            ENDIF */
/*            J = J+1 */
/*         ENDDO findgroup */
/* * */
/*         IF( I.LT.J .AND. J.GE.WIL .AND. I.LE.WIU )THEN */
/* *           Cluster with wanted eigenvalues, do verify the bounds */

/* *           Lower bound */
/*            LB = QDVALS(I) */
/*            GOTBND = .FALSE. */
/*            K = 1 */
/*            getlb: DO */
/*               XI = DLAXRN( N, REPR, REPI, LB ) */
/* @STATS on */
/*               XNBIS_INIT = XNBIS_INIT + 1 */
/* @STATS ! */
/*               GOTBND = ( XI .LE. 2*I-1 ) */
/*               IF( GOTBND .OR. K.GT.NRELAX )  EXIT getlb */
/*               LB = QDVALS(I) - ABS(QDVALS(I))*ARXFAC(K) */
/*               K =  K+1 */
/*            ENDDO getlb */
/*            IF( .NOT.GOTBND )THEN */
/*               INFO = 1 */
/*               EXIT outer */
/*            ENDIF */

/* *           Upper bound */
/*            UB = QDVALS(I) */
/*            GOTBND = .FALSE. */
/*            K = 1 */
/*            getub: DO */
/*               XI = DLAXRN( N, REPR, REPI, UB ) */
/* @STATS on */
/*               XNBIS_INIT = XNBIS_INIT + 1 */
/* @STATS ! */
/*               GOTBND = ( XI .GE. 2*J-1 ) */
/*               IF( GOTBND .OR. K.GT.NRELAX )  EXIT getub */
/*               UB = QDVALS(J) + ABS(QDVALS(J))*ARXFAC(K) */
/*               K = K+1 */
/*            ENDDO getub */
/*            IF( .NOT.GOTBND )THEN */
/*               INFO = 1 */
/*               EXIT outer */
/*            ENDIF */

/*         ELSE */
/* *           Take the bounds unverified. We just inflate them a bit */
/* *           to avoid equal lower and upper bounds. */
/*            LB = QDVALS(I) - 2*EPS*ABS(QDVALS(I)) */
/*            UB = QDVALS(J) + 2*EPS*ABS(QDVALS(J)) */
/*         ENDIF */
/* * */
/*         DO K = I, J */
/*            EWL_AE(2*K-1) = I */
/*            EWL_AE(2*K)   = J */
/*            EWL_LU(2*K-1) = LB */
/*            EWL_LU(2*K)   = UB */
/*         ENDDO */
/* * */
/*         IF( J .EQ. IU )  EXIT */
/*         I = J+1 */
/*      ENDDO outer */
    return 0;
} /* dlaxre_initewldqds__ */

