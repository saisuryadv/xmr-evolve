*******************************************************************************
*
      SUBROUTINE DLACSV(N, B1, R, BN, NDEPTH, FACTRL,
     $	A, B, A2, B2, AB00, AB10, L, UMN, S, P, 
     $  SIGMA, EVAL, GERSCH, DCOUP, LCOUP,
     $  V, U, ISUPPU, WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            N, B1, R, BN, NDEPTH
      DOUBLE PRECISION   SIGMA, EVAL
*     ..
*     .. Array Arguments ..
      INTEGER		 FACTRL(*), ISUPPU(*)
      DOUBLE PRECISION   A(*), B(*), L(*), UMN(*), S(*), P(*),
     $                   A2(*), B2(*), AB00(*), AB10(*),
     $                   DCOUP(*), LCOUP(*), U(*), V(*), GERSCH(*),
     $                   WORK(*)
*
*  Purpose
*  =======
*
*  DLACSV determines coupled singular vectors.
*  If (NDEPTH.EQ.0) the routines DLAS2P, DLAP2S and DLAG2G are called
*  to find an appropriate twisted representation for positive definite 
*  initial matrices using certain coupling transformations.
*  For (NDEPTH.GT.0) the corresponding singular vector is determined
*  directly from the coupled representation LCOUP DCOUP LCOUP^T .
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
*  R       (input) INTEGER
*          Choosen index for the optimal twist element
*          of the explicit factorization.
*
*  BN      (input) INTEGER
*          Last index of splitted submatrix.
*
*  NDEPTH  (input) INTEGER
*          Indicates if the initial matrix is positive definite 
*          (NDEPTH.EQ.0) or not.
*
*  FACTRL  (input) INTEGER array, dimension (4)
*          Contains algorithmic details from the explicit factorization,
*          is set by DLAR1V .
*
*  A       (input) DOUBLE PRECISION array, dimension (N)
*          The diagonal entries of the bidiagonal matrix B.
*
*  B       (input) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal entries of the bidiagonal matrix B.
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
*  L       (input) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of the lower unit bidiagonal
*          factor computed using the dstqds algorithm.
*
*  UMN     (input) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of the upper unit bidiagonal
*          factor computed using the dqds algorithm.
*
*  S       (input) DOUBLE PRECISION array, dimension (N)
*          Auxiliary variables from the dstqds algorithm.
*
*  P       (input) DOUBLE PRECISION array, dimension (N)
*          Auxiliary variables from the dqds algorithm.
*
*  SIGMA   (input) DOUBLE PRECISION
*          The shift parameter.
*
*  EVAL    (input) DOUBLE PRECISION
*          The refined eigenvalue of the explicitly factorized
*	   representation L D L^T .
*
*  GERSCH  (input) DOUBLE PRECISION array, dimension (N)
*          The (refined) Gerschgorin intervals.
*
*  DCOUP   (input) DOUBLE PRECISION array, dimension (N)
*          Diagonal elements of the diagonal matrix of the
*	   coupled representation LCOUP DCOUP LCOUP^T .
*
*  LCOUP   (input) DOUBLE PRECISION array, dimension (N-1)
*          Off-diagonal elements of the lower unit bidiagonal matrix 
*	   of the coupled representation LCOUP DCOUP LCOUP^T .
*
*  V       (input) DOUBLE PRECISION array, dimension (N)
*          The right singular vector (needed for arranging signs).
*
*  U       (output) DOUBLE PRECISION array, dimension (N)
*          The coupled left singular vector.
*
*  ISUPPU  (output) INTEGER array
*          Keeps the support information. The smallest and largest 
*          indices of the non-zero elements are stored in ISUPPU(1) 
*          and ISUPPU(2) resp.
*
*  U       (output) DOUBLE PRECISION array, dimension (N)
*          The coupled left singular vector.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (5*N)
*          Workspace.
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER 		 R1, R2, RCIND, I,
     $                   INDDLL, INDWRK, INDLC, INDSC, INDUC, INDPC,
     $                   UTO, UFROM
      LOGICAL            NANINS, NANINP
      DOUBLE PRECISION   MINGMA, UTU, TMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAG2G, DLAP2S, DLAS2P, DLAR1V, DLATS1, 
     $			 DLATS2, DSCAL
*     ..
*     .. Executable Statements ..
*     ..
      INDLC = 1
      INDPC = N+1
      INDUC = 2*N+1
      INDSC = 3*N+1
      INDWRK = 4*N+1
      INDDLL = INDWRK
*
*     Extract information from FACTRL
*
      R1 = FACTRL(1)
      R2 = FACTRL(2)
      NANINS = (FACTRL(3).EQ.1) 
      NANINP = (FACTRL(4).EQ.1) 
*
      IF (NDEPTH.EQ.0) THEN
*
*       Determine couplings for positive definite intitial matrices.
*
	CALL DLAS2P(N, B1, R2, A2, B2, AB00, AB10, 
     $    L, S, SIGMA, NANINS, 
     $	  WORK( INDLC ), WORK( INDPC ) )
	CALL DLAP2S(R1, BN, A2, AB00, AB10, 
     $	  UMN, P, SIGMA, NANINP, 
     $	  WORK( INDUC ), WORK( INDSC ) )
      	CALL DLAG2G(R1, R2, R, S, P, WORK( INDSC ), 
     $    SIGMA, MINGMA, RCIND )
*
*	Solve linear system with a matrix given as a twisted
*	factorization with twist index RCIND.
*
	CALL DLATS1( U, UTU, 
     $	  WORK( INDLC ), WORK( INDUC ), AB10, 
     $	  B1, RCIND, BN, ISUPPU, (NANINS.OR.NANINP) )
        CALL DLATS2( U, UTU, MINGMA, WORK(INDLC), WORK(INDUC), AB10, 
     $    ISUPPU(1), RCIND, ISUPPU(2), (NANINS.OR.NANINP), 
     $	  WORK( INDWRK ) )
      ELSE
*
*       Compute coupled singular vectors directly from the
*	coupled representation LCOUP DCOUP LCOUP^T .
*
	RCIND = 0
	DO 10, I = 1,N-1
	  WORK( INDDLL+I-1 ) = AB10(I)*LCOUP(I)
   10	CONTINUE
	CALL DLAR1V(B1, BN, SIGMA,
     $	  DCOUP, LCOUP, AB10, WORK( INDDLL ),
     $	  EVAL, GERSCH, U, UTU, MINGMA, RCIND, ISUPPU, FACTRL,
     $	  WORK( INDLC ), WORK( INDUC ),
     $	  WORK( INDSC ), WORK( INDPC ) )
        CALL DLATS2( U, UTU, MINGMA, 
     $	  WORK( INDLC ), WORK( INDUC ), AB10, 
     $    ISUPPU( 1 ), RCIND, ISUPPU( 2 ), 
     $    ((FACTRL(3).EQ.1).OR.(FACTRL(4).EQ.1)),
     $    WORK( INDWRK ) )
      ENDIF
*
*     Arrange signs for U and V .
*
*     Compute temp = r-th component of bidiagonal * V .
*
      IF (R.LT.N) THEN
        TMP = A(R)*V(R)+B(R)*V(R+1)
      ELSE
        TMP = A(N)*V(N)
      ENDIF
*
*     Compare signs of TMP and U(R) and scale appropriately.
*
      UFROM = ISUPPU( 1 )
      UTO = ISUPPU( 2 )
      IF (SIGN(ONE,TMP*U(R)).LT.ZERO) THEN
        CALL DSCAL(UTO-UFROM+1,-ONE/SQRT(UTU),U(UFROM),1)
      ELSE
        CALL DSCAL(UTO-UFROM+1,ONE/SQRT(UTU),U(UFROM),1)
      ENDIF
*
      RETURN
      END
*
*******************************************************************************
