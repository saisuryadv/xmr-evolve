*******************************************************************************
*
      SUBROUTINE DLATS2( Z, ZTZ, MINGMA, L, U, LD, B1, R, BN, 
     $  		 SAWNAN, X )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER B1, R, BN
      LOGICAL SAWNAN
      DOUBLE PRECISION ZTZ, MINGMA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION Z( * ), L( * ), U( * ), LD( * ), X( * )
*
*  Purpose
*  =======
*
*  Solves N G N^T * x = z for x and copies x to z .
*  G is diagonal with G(R,R) = MINGMA, whereas N is a twisted matrix:
*  N = eye(n) + diag([L(1),...,L(R-1),0,...,0])
*             + diag([0,...,0,U(R),...,U(n-1)])
*  Thus only MINGMA, L and U is needed instead of the symbolic 
*  matrices G and N.
*
*  DLATS2 is ashortcut for DLA + T(wisted matrix) S(olver) 2 .
*
*  Arguments
*  =========
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (N)
*          On input, the solution computed with DLATS1 .
*          On output, a refined solution, typically
*          an improved eigenvector approximation.
* 
*  ZTZ     (output) DOUBLE PRECISION
*          Dot product ZTZ = Z^T Z .
*
*  MINGMA  (input) DOUBLE PRECISION
*          The twist element.
*
*  L       (input) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of the lower unit bidiagonal
*          factor computed using the dstqds algorithm.
*
*  U       (input) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of the upper unit bidiagonal
*          factor computed using the dqds algorithm.
*
*  LD      (input) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of L D L^T .
*
*  B1      (input) INTEGER
*          First index of splitted submatrix.
*
*  R       (input) INTEGER
*          The twist position.
*
*  BN      (input) INTEGER
*          Last index of splitted submatrix.
*
*  SAWNAN  (input) LOGICAL
*          Indicates if a breakdown occured.
*
*  X       (workspace) DOUBLE PRECISION array, dimension (N)
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IM1
*
*     Do not solve anything if a NaN occured.
*
      IF (SAWNAN.OR.(MINGMA.EQ.ZERO)) RETURN
*
*     Solve N*x = z .
*
      DO 10, I = B1, BN
	X(I) = ZERO
   10 CONTINUE
      IF (R.EQ.BN) THEN
        X(B1) = Z(B1)
	IM1 = B1
        DO 20, I = B1+1,BN 
          X(I) = Z(I) - L(IM1)*X(IM1)
	  IM1 = I
   20   CONTINUE
      ENDIF
      IF (R.EQ.B1) THEN
        X(BN) = Z(BN)
        DO 30, I = BN-1,R,-1 
          X(I) = Z(I) - U(I)*X(I+1)
   30   CONTINUE
      ENDIF
      IF ((B1.LT.R).AND.(R.LT.BN)) THEN
        X(B1) = Z(B1)
 	IM1 = B1
        DO 40, I = B1+1,R-1 
          X(I) = Z(I) - L(IM1)*X(IM1)
	  IM1 = I
   40   CONTINUE
        X(BN) = Z(BN)
        DO 50, I = BN-1,R+1,-1 
          X(I) = Z(I) - U(I)*X(I+1)
   50   CONTINUE
        X(R) = Z(R) - L(R-1)*X(R-1) - U(R)*X(R+1)
      ENDIF
*
*     Solve G*z = x .
*
      DO 60, I = B1,R-1 
        Z(I) = X(I)*L(I)/LD(I)
   60 CONTINUE
      Z(R) = X(R)/MINGMA
      IM1 = R
      DO 70, I = R+1,BN 
        Z(I) = X(I)*U(IM1)/LD(IM1)
	IM1 = I
   70 CONTINUE
*
*     Solve NT*x = z .
*
      X(R) = Z(R)
      DO 80, I = R-1,B1,-1 
        X(I) = Z(I) - L(I)*X(I+1)
   80 CONTINUE
      IM1 = R
      DO 90, I = R+1,BN
        X(I) = Z(I) - U(IM1)*X(IM1)
	IM1 = I
   90 CONTINUE
*
*     Copy x to z and compute Z^T Z .
*
      ZTZ = ZERO
      DO 100, I = B1,BN
	Z(I) = X(I)
	ZTZ = ZTZ + Z(I)*Z(I)
  100 CONTINUE
*
      RETURN
*
*     END OF DLATS2
*
      END
*
*******************************************************************************
