*******************************************************************************
*
      SUBROUTINE DLAG2G(R1, R2, RF, SFACT, PFACT, SCOUP, MU, GCOUMN, RC)
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            R1, R2, RF, RC
      DOUBLE PRECISION   MU, GCOUMN
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   SFACT(*), PFACT(*), SCOUP(*) 
*
*  Purpose
*  =======
*
*  After performing the coupling transformation 
*  [LFACT, SFACT, MU] --> [LCOUP, PCOUP] (in DLAS2P) and
*  [UFACT, PFACT, MU] --> [UCOUP, UCOUP] (in DLAP2S)
*  an appropriate coupled twisted representation of B^T B - MU * I 
*  is generated yielding twist position RC and ist value GCOUMN .
*
*  Arguments
*  =========
*
*  R1      (input) INTEGER
*          Smallest index for moderate twist elements.
*
*  R2      (input) INTEGER
*          Largest index for moderate twist elements.
*
*  RF      (input) INTEGER
*          Choosen index for the optimal twist element
*          of the explicit factorization.
*
*  SFACT   (input) DOUBLE PRECISION array, dimension (N)
*          Auxiliary variables from the dstqds algorithm.
*
*  PFACT   (input) DOUBLE PRECISION array, dimension (N)
*          Auxiliary variables from the dqds algorithm.
*
*  SCOUP   (input) DOUBLE PRECISION array, dimension (N)
*	   Corresponding auxiliary variables generated
*          by coupling transformations in DLAP2S .
*
*  MU      (input) DOUBLE PRECISION
*          The shift parameter.
*
*  GCOUMN  (output) DOUBLE PRECISION
*          An acceptable twist element of the coupled twisted 
*          representation.
*
*  RC      (output) INTEGER
*          An acceptable twist position of the coupled twisted 
*          representation.
*
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TWO
      PARAMETER          ( ZERO = 0.0D0, TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER I
      DOUBLE PRECISION  EPS, GFACT, GCOUP, GFACMN
*     ..
*     .. Executable Statements ..
*     ..
      EPS = DLAMCH( 'Precision' )  
*
*     Find out if GCOUP(RF) denotes an acceptable twist position
*     by checking (ABS(GCOUP(RF)).LE.TWO*ABS(GFACT(RF)) .
*
      GFACMN = SFACT(RF) + PFACT(RF)
      IF (GFACMN.EQ.ZERO) THEN
        RC = RF
        GCOUMN = ZERO
        RETURN
      ELSE 
        GCOUMN = GFACMN*SCOUP(RF)/(SFACT(RF)-MU)
        IF (GCOUMN.EQ.ZERO) THEN
          GCOUMN = EPS*(SCOUP(RF)+MU)
        ENDIF
      ENDIF 
      RC = RF
      IF (ABS(GCOUMN).LE.TWO*ABS(GFACMN))
     $	RETURN
*
*     Determine an acceptable twist position from GFACT .
*
      DO 10, I = R1,R2
        GFACT = SFACT(I) + PFACT(I)
        GCOUP = GFACT*SCOUP(I)/(SFACT(I)-MU)
        IF (GCOUP.EQ.ZERO) THEN
          GCOUP = EPS*(SCOUP(I)+MU)
        ENDIF
        IF (ABS(GCOUP).LT.ABS(GCOUMN)) THEN
          RC = I
          GCOUMN = GCOUP
          IF (ABS(GCOUP).LE.TWO*ABS(GFACMN))
     $	    RETURN
        ENDIF
   10 CONTINUE
*
      RETURN
      END
*
*******************************************************************************
