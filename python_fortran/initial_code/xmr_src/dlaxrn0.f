      FUNCTION DLAXRN0(N, K, D, OMEGAD, R, OMEGAR, GAMMA)
      IMPLICIT NONE
      INTEGER  ::  DLAXRN0
      INTEGER,          INTENT(IN)  ::  N, K
      INTEGER,          INTENT(IN)  ::  OMEGAD(1:K), OMEGAR(K:N)
      DOUBLE PRECISION, INTENT(IN)  ::  D(1:K-1), R(K+1:N), GAMMA
*
*  Purpose
*  =======
*
*    Specialised routine to determine zero inertia from base data.
*
*  ======================================================================
*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
*
*     .. Declarations ..
*
*
*     .. Locals ..
*
      INTEGER  I, NEG, XI
*
*  ----- Executable Statements ------------------------------------------
*
      NEG = 0
*
      I = 1
      DO
         IF( I .GE. K ) EXIT
         IF( OMEGAD(I+1) .NE. 0 )THEN
            NEG = NEG+1
            I = I+1
         ELSEIF( D(I) .LT. ZERO )THEN
            NEG = NEG+1
         ENDIF
         I = I+1
      ENDDO
*
      I = N
      DO
         IF( I .LE. K ) EXIT
         IF( OMEGAR(I-1) .NE. 0 )THEN
            NEG = NEG+1
            I = I-1
         ELSEIF( R(I) .LT. ZERO )THEN
            NEG = NEG+1
         ENDIF
         I = I-1
      ENDDO
*
      XI = 2 * NEG
      IF( (K.GT.1 .AND. OMEGAD(K).EQ.0) .OR.
     $    (K.LT.N .AND. OMEGAR(K).EQ.0) )
     $ THEN
         IF( GAMMA .LT. ZERO )  XI = XI + 2
         IF( GAMMA .EQ. ZERO )  XI = XI + 1
      ENDIF
*
      DLAXRN0 = XI
      END FUNCTION DLAXRN0
*
************************************************************************

