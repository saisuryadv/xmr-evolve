      SUBROUTINE DLAXRC(
     $            N, REPR, REPI,
     $            IL, IU, LGAP, UGAP, EWL_AE, EWL_LU,
     $            NUMEW, EWINDS,
     $            RELTOL, ABSTOL, RWORK, IWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N, IL, IU, NUMEW
      INTEGER,          INTENT(IN)     ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)     ::  RELTOL, ABSTOL
      DOUBLE PRECISION, INTENT(IN)     ::  REPR(4*N+3)
*
      INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU)
      INTEGER,          INTENT(INOUT)  ::  EWINDS(NUMEW)
      INTEGER,          INTENT(INOUT)  ::  IWORK( 2*NUMEW + IU-IL+1 )
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU)
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( NUMEW )
*
*
*  Purpose
*  =======
*
*     Refine eigenvalues with indices given in EWINDS until they have
*     relative accuracy RTOL or absolute accuracy ABSTOL.
*
*  Arguments
*  =========
*
*  NUMEW   (input/output) INTEGER,  <= IU-IL+1
*          On entry, the number of eigenvalues to refine.
*          Overwritten on exit.
*
*  EWINDS  (input/output) INTEGER array, dimension (NUMEW)
*          On input, the indices of the eigenvalues to refine.
*          Need not be sorted and may even contain duplicates, the
*          no unnessessary bisection steps will be done.
*
*  RELTOL  (input) DOUBLE PRECISION, >= 0
*  ABSTOL  (input) DOUBLE PRECISION, >= 0
*          Tolerances for refinement. An interval is converged if
*           (1)  |ub - lb| <= absmax(lb,ub) * RELTOL
*                or RELTOL = 0
*          and
*           (2)  |ub - lb| <= ABSTOL
*                or  ABSTOL = 0.
*
*          Thus by setting RELTOL or ABSTOL to zero, the corresponding
*          test can be deactivated.
*          In any case an interval is also converged if it so tight that
*          the midpoint is computed as one of lb or ub, which normally
*          means there is no fp-number between them (ub = next(lb)).
*
*
*  ======================================================================
*
*     .. Declarations ..
*
      INTERFACE
      SUBROUTINE DLAXRL_REFINE(
     $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, I, J, LAMBDA, XI
     $           )
      IMPLICIT NONE
*
      INTEGER         , INTENT(IN)  ::  IL, IU, I, J, XI
      DOUBLE PRECISION, INTENT(IN)  ::  LAMBDA
*
      INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU)
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU)
      END SUBROUTINE DLAXRL_REFINE
      END INTERFACE
      INTERFACE
      SUBROUTINE  DLAXRM( N, REPR, REPI, M, ATAU, AXI )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, M
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), ATAU(M)
*
      INTEGER,          INTENT(OUT)  ::  AXI(M)
      END SUBROUTINE DLAXRM
      END INTERFACE


*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  HALF = 0.5D0
*
*     .. Locals ..
*
      DOUBLE PRECISION  MID, LB, UB, WIDTH, ABSMAX
      INTEGER           JXFLGS, JXBIND, JXINER
      INTEGER           I, J, K, M, NTOBIS
      INTEGER IXF77A, IXF77B, IXF77C
      LOGICAL           CHKABS, CHKREL
*
*     -- Executable Statements -----------------------------------------
*

*     ------------
*      Tolerances
*     ------------
      CHKABS = (ABSTOL .NE. ZERO)
      CHKREL = (RELTOL .NE. ZERO)
*
*     -----------------------
*      Auxiliary data fields
*     -----------------------
*      FLGS(IL:IU)
*        FLGS(I) is set to 1 to indicate that the interval with left
*        index I is already included in the current sweep
*      BIND(1:NUMEW)  left indices of intervals in the current sweep
*      INER(1:NUMEW)  inertias of samples in the current sweep
*      RWORK(1:NUMEW)  midpoints of intervals in the current sweep
*
      JXFLGS = 1-IL
      JXBIND = JXFLGS + IU
      JXINER = JXBIND + NUMEW
*
      M = NUMEW
*
      DO
*
*        -----------------------
*         Check for convergence
*        -----------------------
*
         DO IXF77A = JXFLGS+IL, JXFLGS+IU
            IWORK(IXF77A) = 0
         ENDDO
         NTOBIS = 0

         K = 1
         DO
            IF( K .GT. M ) EXIT

            I   = EWL_AE(2*EWINDS(K)-1)
            LB  = EWL_LU(2*I-1)
            UB  = EWL_LU(2*I)
            MID = HALF*(LB + UB)
            WIDTH  = UB - LB
            ABSMAX = MAX(ABS(LB),ABS(UB))

            IF( IWORK(JXFLGS + I) .NE. 0 )THEN
*              corresponding interval was already processed
               K = K+1
            ELSE
               IF( MID.EQ.LB .OR. MID.EQ.UB .OR.
     $             ( (.NOT.CHKABS .OR. WIDTH.LE.ABSTOL) .AND.
     $               (.NOT.CHKREL .OR. WIDTH.LE.ABSMAX*RELTOL) ) )
     $         THEN
*                 this interval is converged
                  EWINDS(K) = EWINDS(M)
                  M = M - 1
               ELSE
*                 not converged, add midpoint for one bisection step
                  NTOBIS = NTOBIS + 1
                  IWORK( JXBIND + NTOBIS ) = I
                  RWORK( NTOBIS ) = MID
*                 set flag that this interval was already processed
                  IWORK( JXFLGS + I ) = 1
                  K = K + 1
               ENDIF
            ENDIF
         ENDDO


*
*        -----------
*         Loop Exit
*        -----------
*
         IF( M .EQ. 0 )  EXIT
*
*        ---------------------------------------
*         Bisect all unconverged intervals once
*        ---------------------------------------
*
         CALL DLAXRM( N, REPR, REPI, NTOBIS, RWORK, IWORK(JXINER+1) )

         DO K = 1, NTOBIS
            I = IWORK(JXBIND + K)
            J = EWL_AE(2*I)

            CALL DLAXRL_REFINE(
     $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU,
     $             I, J, RWORK(K), IWORK(JXINER+K)
     $      )
         ENDDO

      ENDDO

      END SUBROUTINE DLAXRC
*
************************************************************************

