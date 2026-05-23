      SUBROUTINE DLAXRE_INITEWLDQDS(
     $             N, REPR, REPI, IL, IU,
     $             QDVALS, WIL, WIU, EWL_AE, EWL_LU
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)   ::  N, IL, IU, WIL, WIU
      INTEGER,          INTENT(IN)   ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)   ::  QDVALS(IL:IU)
      DOUBLE PRECISION, INTENT(IN)   ::  REPR(4*N+3)
*
      INTEGER,          INTENT(OUT)  ::  EWL_AE(2*IL-1:2*IU)
      DOUBLE PRECISION, INTENT(OUT)  ::  EWL_LU(2*IL-1:2*IU)
*
*
*  Purpose
*  =======
*
*    Init the ewlist for eigenvalues IL:IU from the bounds given
*    in QDVALS. Those are assumed to stem from DQDS, that is, be
*    aproximations with high relative accuracy.
*
*    The bounds will be relaxed a bit, but not verified.
*
* !   NOTE: At the moment the consistency requirements in DSTEXR allow
* !   only to use dqds for blocks where all ews are wanted. Thus this
* !   routine is only used with WIL=1, WIU=N. Future versions without
* !   consistency guarantees might use the partial feature here though.
*
*  Arguments
*  =========
*
*  QDVALS  (input), DOUBLE PRECISION array, dimension(IL:IU)
*          The eigenvalue approximations, in ascending order.
*          We do not require them to have the same sign or be distinct.
*          However, the list should not contain more than one zero.
*
*  ======================================================================
*
      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH

*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
*
*     .. Locals ..
*
      DOUBLE PRECISION  PREC, SEPTOL, RELOFF
      DOUBLE PRECISION  ABSMAX, GAP, LB, UB
      INTEGER           I, J, K, XI
*
*  ===== Executable Statements ==========================================
*

      PREC = DLAMCH('Precision')
      SEPTOL = SQRT(PREC)
      RELOFF = LOG(DBLE(N))*PREC

*     Group the values, and then relax the outer bounds of those
*     groups just a bit to avoid unnecessary backing off in DLAXRF.
*     [ Note: Carful with the backing off only to to at group boundaries,
*       otherwise the monotonicity of bounds in the EWL may be destroyed]

      I = IL
      outer: DO
         J = I
         findgroup: DO
            IF( J .EQ. IU )THEN
               EXIT
            ENDIF
            GAP    = QDVALS(J+1) - QDVALS(J)
            ABSMAX = MAX( ABS(QDVALS(J)), ABS(QDVALS(J+1)) )
            IF( GAP .GT. ABSMAX*SEPTOL )THEN
               EXIT findgroup
            ENDIF
            J = J+1
         ENDDO findgroup

         DO K = I, J
            EWL_AE(2*K-1) = K
            EWL_AE(2*K)   = K
            EWL_LU(2*K-1) = QDVALS(K)
            EWL_LU(2*K)   = QDVALS(K)
         ENDDO

*        Note: For singletons the midpoint is taken, so the following
*        has no effect there, the starting value for RQI will still be
*        QDVALS(I).
         EWL_LU(2*I-1) = QDVALS(I) - RELOFF*ABS(QDVALS(I))
         EWL_LU(2*J)   = QDVALS(J) + RELOFF*ABS(QDVALS(J))

         IF( J .EQ. IU )THEN
            EXIT
         ENDIF
         I = J+1
      ENDDO outer



C      ARXFAC(1) = LOG( DBLE(N) ) * (4 * EPS)
C      ARXFAC(2) =   4 * N * EPS
C      ARXFAC(3) =  32 * N * EPS
C      ARXFAC(4) = 400 * N * EPS
C
C      I = IL
C      outer: DO
C         J = I
C         findgroup: DO
C            IF( J .EQ. IU )THEN
C               EXIT
C            ENDIF
C            GAP    = QDVALS(J+1)-QDVALS(J)
C            ABSMAX = MAX( ABS(QDVALS(J)), ABS(QDVALS(J+1)) )
C            IF( GAP .GT. ABSMAX*GAPTOL )THEN
C               EXIT findgroup
C            ENDIF
C            J = J+1
C         ENDDO findgroup
C*
C         IF( I.LT.J .AND. J.GE.WIL .AND. I.LE.WIU )THEN
C*           Cluster with wanted eigenvalues, do verify the bounds
C
C*           Lower bound
C            LB = QDVALS(I)
C            GOTBND = .FALSE.
C            K = 1
C            getlb: DO
C               XI = DLAXRN( N, REPR, REPI, LB )
C@STATS on
C               XNBIS_INIT = XNBIS_INIT + 1
C@STATS !
C               GOTBND = ( XI .LE. 2*I-1 )
C               IF( GOTBND .OR. K.GT.NRELAX )  EXIT getlb
C               LB = QDVALS(I) - ABS(QDVALS(I))*ARXFAC(K)
C               K =  K+1
C            ENDDO getlb
C            IF( .NOT.GOTBND )THEN
C               INFO = 1
C               EXIT outer
C            ENDIF
C
C*           Upper bound
C            UB = QDVALS(I)
C            GOTBND = .FALSE.
C            K = 1
C            getub: DO
C               XI = DLAXRN( N, REPR, REPI, UB )
C@STATS on
C               XNBIS_INIT = XNBIS_INIT + 1
C@STATS !
C               GOTBND = ( XI .GE. 2*J-1 )
C               IF( GOTBND .OR. K.GT.NRELAX )  EXIT getub
C               UB = QDVALS(J) + ABS(QDVALS(J))*ARXFAC(K)
C               K = K+1
C            ENDDO getub
C            IF( .NOT.GOTBND )THEN
C               INFO = 1
C               EXIT outer
C            ENDIF
C
C         ELSE
C*           Take the bounds unverified. We just inflate them a bit
C*           to avoid equal lower and upper bounds.
C            LB = QDVALS(I) - 2*EPS*ABS(QDVALS(I))
C            UB = QDVALS(J) + 2*EPS*ABS(QDVALS(J))
C         ENDIF
C*
C         DO K = I, J
C            EWL_AE(2*K-1) = I
C            EWL_AE(2*K)   = J
C            EWL_LU(2*K-1) = LB
C            EWL_LU(2*K)   = UB
C         ENDDO
C*
C         IF( J .EQ. IU )  EXIT
C         I = J+1
C      ENDDO outer
      END SUBROUTINE DLAXRE_INITEWLDQDS
*
************************************************************************
