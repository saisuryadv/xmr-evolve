      SUBROUTINE DLAXRL_UPDATE(
     $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, LAMBDA, XI
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  IL, IU, XI
      DOUBLE PRECISION, INTENT(IN)  ::  LAMBDA
*
      INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU)
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU)
*
*  Purpose
*  =======
*
*    Incorporate the information given by sample LAMBDA with inertia
*    XI into the list.
*    Effectively this looks for an interval that contains LAMBDA and
*    refines it, but XI is employed to speed up the search. If LAMBDA
*    falls outside the intervals, the outer gaps will be updated if
*    possible.
*
*  ======================================================================
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
*
*     .. Constants ..
*
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO = 0.0D0)
*
*     .. Locals ..
*
      INTEGER           IX, JX, I, J

*
*     -- Executable Statements -----------------------------------------
*

*     Look at all intervals that might benefit from this sample
      IF( XI.GT.(2*(IL-1)-1) .AND. XI.LT.(2*(IU+1)-1) )THEN
         IX = EWL_AE( 2 * MAX( IL, XI / 2 )  -  1 )
         JX = EWL_AE( 2 * MIN( IU, (XI+1)/2 + 1 ) )

         I = IX
         DO
            J = EWL_AE(2*I)

            IF( EWL_LU(2*I-1).LT.LAMBDA .AND.
     $          LAMBDA.LT.EWL_LU(2*J) )
     $      THEN
               CALL DLAXRL_REFINE(
     $                IL, IU, LGAP, UGAP, EWL_AE, EWL_LU,
     $                I, J, LAMBDA, XI
     $              )
            ENDIF

            IF( J.EQ.JX .OR. LAMBDA.LE.EWL_LU(2*J) )THEN
               EXIT
            ENDIF

            I = J+1
         ENDDO
      ENDIF

*     Update the outer gaps
*     This must be done here, after refining the intervals, since by that
*     the outer bounds may have changed.

      IF( XI.EQ.(2*(IL-1)-1) .OR. XI.EQ.2*(IL-1) )THEN
*        is a valid left outer bound
         LGAP = MAX( LGAP, EWL_LU(2*IL-1)-LAMBDA )
      ENDIF

      IF( XI.EQ.(2*IU) .OR. XI.EQ.(2*IU+1) )THEN
*        is a valid right outer bound
         UGAP = MAX( UGAP, LAMBDA-EWL_LU(2*IU) )
      ENDIF

      END SUBROUTINE DLAXRL_UPDATE
*
************************************************************************
