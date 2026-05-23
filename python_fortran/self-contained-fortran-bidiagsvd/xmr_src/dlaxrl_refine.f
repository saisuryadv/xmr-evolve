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
*
*  Purpose
*  =======
*
*    Refine the interval [i,j] by bisecting it at lambda via the
*    inertia information xi.
*
*    Pre: [I,J] is an interval in the list, LAMBDA lies in the interior
*         of the current bounds.
*
*    Post: Eigenvalues [I,J] are covered by one or more intervals of
*         which the largest diameter is smaller than the diameter of
*         the original interval was.
*
*  ======================================================================
*
*     .. Locals ..
*
      DOUBLE PRECISION  OLDLB, OLDUB
      INTEGER           IE, JA, K


*
*     -- Executable Statements -----------------------------------------
*
*
*     Determine ew indices on the left and right of the sample
*
      IE = XI / 2
      JA = (XI + 1) / 2 + 1
*
*     Special handling of non-monotonicity. We could ignore the
*     inconsistencies, however, this might lead to an infinte
*     loop in bisection.
*
      IF( JA.LE.I .OR. IE.LT.I-1 )THEN
*        In the second condition , the inertia indicates that an ew
*        from the left jumps in.
*        In any case this means LAMBDA can be seen as lower bound
*        for [I,J].
         IE = I-1
         JA = I
      ELSEIF( IE.GE.J .OR. JA.GT.J+1 )THEN
*        Analogously, we can regard LAMBDA as upper bound for [I,J].
         JA = J+1
         IE = J
      ENDIF
*
*     Now we can split the interval i:j with bounds lb:ub into three
*     intervals
*       [lb,lambda)     containing ews  i:ie,
*       [lambda,lambda] containing ew   ie+1:ja-1, and
*       (lambda,ub]     containing ews  ja:j,
*     where the respective index ranges may be empty, but
*     nevertheless we have
*       {i:ie} union {ie+1:ja-1} union {ja:j}  =  {i:j}.
*
      OLDLB = EWL_LU(2*IL-1)
      OLDUB = EWL_LU(2*IU)

      IF( I .LE. IE )THEN
         DO K = I, IE
            EWL_AE(2*K) = IE
            EWL_LU(2*K) = LAMBDA
         ENDDO
      ENDIF

      IF( IE+1 .EQ. JA-1 )THEN
         K = IE+1
         EWL_AE(2*K-1) = K
         EWL_AE(2*K)   = K
         EWL_LU(2*K-1) = LAMBDA
         EWL_LU(2*K)   = LAMBDA
      ENDIF

      IF( JA .LE. J )THEN
         DO K = JA, J
            EWL_AE(2*K-1) = JA
            EWL_LU(2*K-1) = LAMBDA
         ENDDO
      ENDIF

      LGAP = LGAP + (EWL_LU(2*IL-1) - OLDLB)
      UGAP = UGAP + (OLDUB - EWL_LU(2*IU))

      END SUBROUTINE DLAXRL_REFINE
*
************************************************************************
