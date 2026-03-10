*******************************************************************************
*
      SUBROUTINE DLAP2S(R1, BN, A2, AB00, AB10, UFACT, PFACT, 
     $  MU, SAWNAN, UCOUP, SCOUP)
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            R1, BN
      DOUBLE PRECISION   MU
      LOGICAL            SAWNAN
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A2(*), AB10(*), UFACT(*), PFACT(*), UCOUP(*),
     $	                 SCOUP(*), AB00(*)
*
*  Purpose
*  =======
*
*  Perform the coupling transformation 
*  [UFACT, PFACT, MU] --> [UCOUP, UCOUP].
*  UFACT and PFACT are the data form the explicit factorization
*  B^T B - MU * I = L D L^T - MU * I = UFACT RFACT UFACT^T
*  using the dqds algorithm (UFACT contains auxiliary variables).
*  The entries of the corresponding representation
*  B B^T - MU * I = UCOUP RCOUP UCOUP^T
*  as well as the auxiliary variable SCOUP are computed using couplings.
*
*  Arguments
*  =========
*
*  R1      (input) INTEGER
*          Smallest index for moderate twist elements.
*
*  BN      (input) INTEGER
*          Last index of splitted submatrix.
*
*  A2      (input) DOUBLE PRECISION array, dimension (N)
*          Squared diagonal elements of the L D L^T representation.
*
*  AB00    (input) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of B^T B .
*
*  AB10    (input) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of B B^T .
*
*  UFACT   (input) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of the upper unit bidiagonal
*          factor computed using the dqds algorithm.
*
*  PFACT   (input) DOUBLE PRECISION array, dimension (N)
*          Auxiliary variables from the dqds algorithm.
*
*  MU      (input) DOUBLE PRECISION
*          The shift parameter.
*
*  SAWNAN  (input) LOGICAL
*          Indicates if a breakdown occured.
*
*  UCOUP   (output) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of the corresponding upper unit 
*          bidiagonal factor computed using coupling transformations.
*
*  SCOUP   (output) DOUBLE PRECISION array, dimension (N)
*	   Corresponding auxiliary variables by coupling transformations.
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   RFACT, RCOUP, T
      INTEGER I
*     ..
*     .. Executable Statements ..
*     ..
	IF (SAWNAN) THEN
*
*         Run slower, NaN-proof code
*
          RCOUP = PFACT(BN)
          SCOUP(BN) = -MU
          T = ONE/PFACT(BN)
          I = BN-1
 10       CONTINUE
          IF (I.GE.R1) THEN
	    RFACT = AB00(I)/UFACT(I)
            IF (RFACT.EQ.ZERO) THEN
              UCOUP(I) = AB10(I)/RCOUP
              RCOUP    = A2(I)
              SCOUP(I) = ZERO
*
*             RFACT(2).EQ.ZERO not possible
*
              IF (I.EQ.1) GOTO 20
              I = I-1
              UCOUP(I) = AB10(I)/RCOUP
              RCOUP    = A2(I) - MU
              T        = ONE/PFACT(I)
              SCOUP(I) = -MU
              I = I-1
            ELSE
              UCOUP(I) = AB10(I)/RCOUP
	      RFACT    = AB00(I)/UFACT(I)
              RCOUP    = RFACT*PFACT(I)*T
              T        = ONE/PFACT(I)
              SCOUP(I) = -MU*RCOUP*T
              I = I-1
            ENDIF
            GOTO 10
          ENDIF
 20       CONTINUE
	ELSE 
*
*         Run normal code
*
          RCOUP = PFACT(BN)
          SCOUP(BN) = -MU
          T     = ONE/PFACT(BN)
          DO 30, I = BN-1,R1,-1
            UCOUP(I) = AB10(I)/RCOUP
	    RFACT    = AB00(I)/UFACT(I)
            RCOUP    = RFACT*PFACT(I)*T
            T        = ONE/PFACT(I)
            SCOUP(I) = -MU*RCOUP*T
 30       CONTINUE
	ENDIF 
*
      RETURN
      END
*
*******************************************************************************
