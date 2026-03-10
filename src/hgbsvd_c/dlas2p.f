*******************************************************************************
*
      SUBROUTINE DLAS2P(N, B1, R2, A2, B2, AB00, AB10,
     $  LFACT, SFACT, MU, SAWNAN, LCOUP, PCOUP)
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            N, B1, R2
      DOUBLE PRECISION   MU
      LOGICAL            SAWNAN
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A2(*), B2(*), AB10(*), AB00(*),
     $  		 LFACT(*), SFACT(*), LCOUP(*), PCOUP(*)
*
*  Purpose
*  =======
*
*  Perform the coupling transformation 
*  [LFACT, SFACT, MU] --> [LCOUP, PCOUP].
*  LFACT and SFACT are the data form the explicit factorization
*  B^T B - MU * I = L D L^T - MU * I = LFACT DFACT LFACT^T
*  using the dstqds algorithm (SFACT contains auxiliary variables).
*  The entries of the corresponding representation
*  B B^T - MU * I = LCOUP DCOUP LCOUP^T
*  as well as the auxiliary variable PCOUP are computed using couplings.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  B1      (input) INTEGER
*          First index of splitted submatrix.
*
*  R2      (input) INTEGER
*          Largest index for moderate twist elements.
*
*  A2      (input) DOUBLE PRECISION array, dimension (N)
*          Squared diagonal elements of the L D L^T representation.
*
*  B2      (input) DOUBLE PRECISION array, dimension (N-1)
*          Squared off-diagonal elements of the L D L^T representation.
*
*  AB00    (input) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of B^T B .
*
*  AB10    (input) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of B B^T .
*
*  LFACT   (input) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of the lower unit bidiagonal
*          factor computed using the dstqds algorithm.
*
*  SFACT   (input) DOUBLE PRECISION array, dimension (N)
*          Auxiliary variables from the dstqds algorithm.
*
*  MU      (input) DOUBLE PRECISION
*          The shift parameter.
*
*  SAWNAN  (input) LOGICAL
*          Indicates if a breakdown occured.
*
*  LCOUP   (output) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of the corresponding lower unit 
*          bidiagonal factor computed using coupling transformations.
*
*  PCOUP   (output) DOUBLE PRECISION array, dimension (N)
*	   Corresponding auxiliary variables by coupling transformations.
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   DCOUP, DFACT, T
      INTEGER I
*     ..
*     .. Executable Statements ..
*     ..
	IF (SAWNAN) THEN
*
*         Run slower, NaN-proof code
*
          I = B1
 10       CONTINUE
          IF (I.LE.R2-1) THEN
	    DFACT = AB00(I)/LFACT(I)
            IF (DFACT.EQ.ZERO) THEN
              DCOUP = B2(I)
              PCOUP(I) = ZERO
              LCOUP(I) = AB10(I)/DCOUP
*
*             DFACT(R2-1).EQ.ZERO not possible
*
              IF (I.LT.R2-1) THEN
                DCOUP = SFACT(I+2) - MU
                PCOUP(I+1) = -MU
                LCOUP(I+1) = AB10(I+1)/DCOUP
                I = I+2
              ELSE
                DCOUP = -MU
                PCOUP(R2) = DCOUP
                GOTO 20
              ENDIF
            ELSE
	      DFACT = AB00(I)/LFACT(I)
              T     = DFACT/(SFACT(I)-MU)
              DCOUP = (SFACT(I+1)-MU)*T
              PCOUP(I) = -MU*T
              LCOUP(I) = AB10(I)/DCOUP
              I = I+1
            ENDIF
            GOTO 10
          ENDIF
 20       CONTINUE
	ELSE 
*
*         Run normal code
*
          DO 30, I = B1,R2-1
	    DFACT    = AB00(I)/LFACT(I)
            T        = DFACT/(SFACT(I)-MU)
            DCOUP    = (SFACT(I+1)-MU)*T
            PCOUP(I) = -MU*T
            LCOUP(I) = AB10(I)/DCOUP
 30	  CONTINUE
	  IF (R2.EQ.N) THEN
	    DFACT = (SFACT(N) - MU) + A2(N)
	    PCOUP(N) = -MU*DFACT/(SFACT(N)-MU)
	  ENDIF
	ENDIF 
*
      RETURN
      END
*
*******************************************************************************
