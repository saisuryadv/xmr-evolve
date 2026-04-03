      SUBROUTINE DLAXRL_RESET(
     $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, L, XIL, U, XIU
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  IL, IU, XIL, XIU
      DOUBLE PRECISION, INTENT(IN)  ::  L, U
*
      INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU)
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU)
*
*  Purpose
*  =======
*
*    Reset the list. We need that L is a valid lower bound for ew IL,
*    and that U is a valid upper bound for ew IU.
*
*  ======================================================================
*
*     .. Constants ..
*
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO = 0.0D0)
*
*     .. Locals ..
*
      INTEGER           I, J

*
*  ----- Executable Statements -----------------------------------------
*

      LGAP = ZERO
      UGAP = ZERO

      DO I = IL, IU
         EWL_LU(2*I-1) = L
         EWL_LU(2*I)   = U
         EWL_AE(2*I-1) = IL
         EWL_AE(2*I)   = IU
      ENDDO

      IF( MOD(XIL,2) .EQ. 1 )THEN
         J = EWL_AE(2*IL)
         EWL_AE(2*IL) = IL
         EWL_LU(2*IL) = L
         DO I = IL+1, J
            EWL_AE(2*I-1) = IL+1
         ENDDO
      ENDIF

      IF( MOD(XIU,2) .EQ. 1 )THEN
         I = EWL_AE(2*IU-1)
         EWL_AE(2*IU-1) = IU
         EWL_LU(2*IU-1) = U
         DO J = I, IU-1
            EWL_AE(2*J) = IU-1
         ENDDO
      ENDIF
      END SUBROUTINE DLAXRL_RESET
*
************************************************************************
