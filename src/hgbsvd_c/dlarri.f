      SUBROUTINE DLARRI( N, WANTL, WANTR, WANTA, D, E, TOL, NSPLIT, 
     $			 ISPLIT, M, W,
     $                   GERSCH, D2, E2, DE00, DE10,
     $			 WORK, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, M, N, NSPLIT
      LOGICAL		 WANTL, WANTR, WANTA	
      DOUBLE PRECISION   TOL
*     ..
*     .. Array Arguments ..
      INTEGER            ISPLIT( * )
      DOUBLE PRECISION   D( * ), E( * ), GERSCH( * ), W( * ), WORK( * ),
     $                   D2( * ), E2( * ), DE00( * ), DE10( * )
*     ..
*     .. Common block to return timings ..
      COMMON             / TTIME / TIMNG
*     ..
*     .. Array in Common ..
      DOUBLE PRECISION   TIMNG( 12 )
*     ..
*
*  Purpose
*  =======
*
*  Given the bidiagonal matrix B, DLARRI checks for "small" elements.
*  If no splittings occur either B^T B = L D L^T (WANTL or WANTA)
*  or B B^T = L D L^T (WANTR) is formed. Note that the decomposition
*  is always positive definite. The eigenvalues are found by the dqds 
*  algorithm (subroutine DLASQ2). As an added benefit, DLARRI also outputs 
*  the n Gerschgorin intervals for the matrices L D L^T.
c
c  Reminder to Osni:
c  Note 4 : Splittings are currently not supported, flagged by INFO = 33!
c
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  WANTL   (input) LOGICAL
*          Only left singular vectors are desired.
*
*  WANTR   (input) LOGICAL
*          Only right singular vectors are desired.
*
*  WANTA   (input) LOGICAL
*          Both left and right singular vectors are desired.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the N diagonal elements of the bidiagonal
*          matrix B.
*          On exit, the N diagonal elements of the diagonal
*          matrix D.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the first (N-1) entries contain the off-diagonal
*          elements of the bidiagonal matrix B; E(N) need not be set.
*          On exit, E contains the subdiagonal elements of the unit
*          bidiagonal matrices L. 
*
*  TOL     (input) DOUBLE PRECISION
*          The threshold for splitting. If on input |D(i)| < TOL or 
*          |E(i)| < TOL, then the matrix B is split into smaller blocks.
*
*  NSPLIT  (output) INTEGER
*          The number of blocks B splits into. 1 <= NSPLIT <= N.
*
*  ISPLIT  (output) INTEGER array, dimension (2*N)
*          The splitting points, at which B breaks up into submatrices.
*          The first submatrix consists of rows/columns 1 to ISPLIT(1),
*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
*          etc., and the NSPLIT-th consists of rows/columns
*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
*
*  M       (output) INTEGER
*          The total number of eigenvalues of L D L^T found.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          The first M elements contain the eigenvalues
*          sorted in ascending order.
*
*  GERSCH  (output) DOUBLE PRECISION array, dimension (2*N)
*          The N Gerschgorin intervals.
*
*  D2      (output) DOUBLE PRECISION array, dimension (N)
*          D2(I) = D(I)*D(I)
*
*  E2      (output) DOUBLE PRECISION array, dimension (N-1)
*          E2(I) = E(I)*E(I)
*
*  DE00    (output) DOUBLE PRECISION array, dimension (N-1)
*          DE00(I) = D(I)*E(I)
*
*  DE01    (output) DOUBLE PRECISION array, dimension (N-1)
*          DE01(I) = D(I+1)*E(I)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension MAX(1,4*N)
*          Workspace.
*
*  IWORK   (workspace) INTEGER array, dimension (2*N)
*          Workspace.
*
*  INFO    (output) INTEGER
*          Output error code from DLASQ2
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Inderjit Dhillon, IBM Almaden, USA
*     Osni Marques, LBNL/NERSC, USA
*     Benedikt Grosser, BUGH Wuppertal, Germany
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IBEGIN, IEND, IN, J, JBLK
      DOUBLE PRECISION   GL, GU, OFFD, SGNDEF, T1, T2, TMP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DSECND
      EXTERNAL           DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLASQ2 
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Compute Splitting Points
*
      NSPLIT = 1
      DO 10 I = 1, N - 1
         IF( ABS( E( I ) ).LE.TOL ) THEN
            ISPLIT( NSPLIT ) = I
            NSPLIT = NSPLIT + 1
         END IF
   10 CONTINUE
      ISPLIT( NSPLIT ) = N
*
*     Check if matrix splits.
*
c
c     Reminder to Osni:
c     Cannot cope with splittings and negative entries.
c
      DO 20, I = 1,N-1
        IF ( (D(I).LE.ZERO).OR.(E(I).LE.ZERO) ) THEN
	  INFO = 33
	  RETURN
        ENDIF
        IF ( (ABS(D(I)).LE.TOL).OR.(ABS(E(I)).LE.TOL) ) THEN
	  INFO = 33
	  RETURN
        ENDIF
   20 CONTINUE
      IF ( D(N).LE.ZERO ) THEN
	INFO = 33
	RETURN
      ENDIF
*
      IBEGIN = 1
      DO 170 JBLK = 1, NSPLIT
         IEND = ISPLIT( JBLK )
         IF( IBEGIN.EQ.IEND ) THEN
            W( IBEGIN ) = D( IBEGIN )
            E( IBEGIN ) = ZERO
            IBEGIN = IEND + 1
            GO TO 170
         END IF
         IN = IEND - IBEGIN + 1
c
c        Reminder to Osni:
c	 Assumption: IBEGIN = 1, IEND = N, IN = N 
c
*
*	 Compute some auxiliary variables.
*
	 DO 30, I = IBEGIN,IEND-1
	   D2(I) = D(I)*D(I)
	   E2(I) = E(I)*E(I)
	   DE00(I) = D(I)*E(I)
	   DE10(I) = D(I+1)*E(I)
   30    CONTINUE
	 D2(IEND) = D(IEND)*D(IEND)
*
*        Form the IN Gerschgorin intervals
*	 If (WANTR) then form T = B^T B explicitly
*	   and find its Gerschgorin intervals.
*	 If (WANTL) then form T = B B^T explicitly
*	   and find its Gerschgorin intervals.
*	 If (WANTA) then unify Gerschgorin intervals 
*	   of B^T B and B B^T.
*
*	 Store input data: WORK(1:N) = D(1:N)
*	                   WORK(N+1:2*N) = E(1:N)
*
	 CALL DCOPY(IN,D(IBEGIN),1,WORK,1)
	 CALL DCOPY(IN-1,E(IBEGIN),1,WORK(N+1),1)
	 CALL DCOPY(2*IN,ZERO,0,GERSCH(2*IBEGIN-1),1)
*
	 IF (WANTR.OR.WANTA) THEN
*
*	   Compute T = B^T B and find its Gerschgorin intervals. 
*
	   D(IBEGIN) = D2(IBEGIN)
	   DO 40, I = IBEGIN,IEND-1
	     D(I+1) = D2(I+1) + E2(I)
	     E(I) = DE00(I)
   40      CONTINUE
           GL = D( IBEGIN ) - ABS( E( IBEGIN ) )
           GU = D( IBEGIN ) + ABS( E( IBEGIN ) )
           GERSCH( 2*IBEGIN-1 ) = GL
           GERSCH( 2*IBEGIN ) = GU
           GERSCH( 2*IEND-1 ) = D( IEND ) - ABS( E( IEND-1 ) )
           GERSCH( 2*IEND ) = D( IEND ) + ABS( E( IEND-1 ) )
           GL = MIN( GERSCH( 2*IEND-1 ), GL )
           GU = MAX( GERSCH( 2*IEND ), GU )
           DO 50 I = IBEGIN + 1, IEND - 1
              OFFD = ABS( E( I-1 ) ) + ABS( E( I ) )
              GERSCH( 2*I-1 ) = D( I ) - OFFD
              GL = MIN( GERSCH( 2*I-1 ), GL )
              GERSCH( 2*I ) = D( I ) + OFFD
              GU = MAX( GERSCH( 2*I ), GU )
   50      CONTINUE
	 ENDIF
	 IF (WANTL.OR.WANTA) THEN
*
*	   Compute T = B B^T and find its Gerschgorin intervals. 
*
	   DO 60, I = IBEGIN,IEND-1
	     D(I) = D2(I) + E2(I)
	     E(I) = DE10(I)
   60      CONTINUE
	   D(IEND) = D2(IEND)
           GL = MIN(GL, D( IBEGIN ) - ABS( E( IBEGIN ) ) )
           GU = MAX(GU, D( IBEGIN ) + ABS( E( IBEGIN ) ) )
           GERSCH( 2*IBEGIN-1 ) = MIN( GL, GERSCH( 2*IBEGIN-1 ) )
           GERSCH( 2*IBEGIN ) = MAX( GU, GERSCH( 2*IBEGIN ) )
           GERSCH( 2*IEND-1 ) = MIN( GERSCH( 2*IEND-1 ), 
     $	     D( IEND ) - ABS( E( IEND-1 ) ) )
           GERSCH( 2*IEND ) = MAX( GERSCH( 2*IEND ),
     $	     D( IEND ) + ABS( E( IEND-1 ) ) )
           GL = MIN( GERSCH( 2*IEND-1 ), GL )
           GU = MAX( GERSCH( 2*IEND ), GU )
           DO 70 I = IBEGIN + 1, IEND - 1
              OFFD = ABS( E( I-1 ) ) + ABS( E( I ) )
              GERSCH( 2*I-1 ) = MIN( GERSCH( 2*I-1 ), D( I ) - OFFD )
              GL = MIN( GERSCH( 2*I-1 ), GL )
              GERSCH( 2*I ) = MAX( GERSCH( 2*I ), D( I ) + OFFD )
              GU = MAX( GERSCH( 2*I ), GU )
   70      CONTINUE
	 ENDIF 
*
*	 Recover input data:
*
	 CALL DCOPY(IN,WORK,1,D(IBEGIN),1)
	 CALL DCOPY(IN-1,WORK(N+1),1,E(IBEGIN),1)
*
*	 Compute  B^T B = L D L^T or  B B^T = L D L^T
*	 Since B^T B and B B^T are already definite, there is no
*	 need to find an initial shift parameter as in DLARRE.
*
	 IF ( WANTL ) THEN
*
*          Compute B B^T = L D L^T .
*	   Input:  
*		D(I), I=1,N - diagonal entries of B
*	        E(1:N-1), I=1,N-1 - offdiagonal entries of B
*	   Output:  
*		WORK(    I) = d_i , I=1,N - diagonal entries of D
*		WORK(2*N+I) = 1/d_i , I=1,N-1, WORK(3*N) = ONE
*		WORK(  N+I) = l_i , I=1,N-1 - offdiagonal entries of L
*
	   TMP = D2(IBEGIN)
	   DO 80, I = IBEGIN,IEND-1
	     WORK(I) = TMP + E2(I)
             WORK(2*N+I) = ONE/WORK(I)
             WORK(N+I) = DE10(I)*WORK(2*N+I)
	     TMP = D2(I+1)*WORK(2*N+I)*TMP
   80      CONTINUE
	   WORK(IEND) = TMP
	   WORK(3*N) = ONE
*
	 ELSE
*
*          Compute B^T B = L D L^T .
*	   Input:  
*		D(I), I=1,N - diagonal entries of B
*	        E(1:N-1), I=1,N-1 - offdiagonal entries of B
*	   Output:  
*		WORK(    I) = d_i , I=1,N - diagonal entries of D
*		WORK(2*N+I) = 1/d_i , I=1,N-1, WORK(3*N) = ONE
*		WORK(  N+I) = l_i , I=1,N-1 - offdiagonal entries of L
*
	   DO 90, I = IBEGIN,IEND-1
	     WORK(I) = D2(I)
             WORK(2*N+I) = ONE/WORK(I)
             WORK(N+I) = E(I)/D(I)
   90      CONTINUE
	   WORK(IN) = D2(IN)
	   WORK(3*N) = ONE
	 ENDIF 
         CALL DCOPY( IN, WORK(IBEGIN), 1, D( IBEGIN ), 1 )
         CALL DCOPY( IN-1, WORK( N+IBEGIN ), 1, E( IBEGIN ), 1 )
*
*	 Initialize some relicts of DLARRE :
*
	 SGNDEF = ONE
         E( IEND ) = ZERO
*
*        Compute the eigenvalues of L D L^T.
*
         J = IBEGIN
         DO 140 I = 1, IN - 1
            WORK( 2*I-1 ) = ABS( D( J ) )
            WORK( 2*I ) = E( J )*E( J )*WORK( 2*I-1 )
            J = J + 1
  140    CONTINUE
         WORK( 2*IN-1 ) = ABS( D( IEND ) )
*
         T1 = DSECND( )
         CALL DLASQ2( IN, WORK, INFO )
         T2 = DSECND( )
         TIMNG( 1 ) = TIMNG( 1 ) + ( T2-T1 )
*
         J = IBEGIN
         IF( SGNDEF.GT.ZERO ) THEN
            DO 150 I = 1, IN
               W( J ) = WORK( IN-I+1 )
               J = J + 1
  150       CONTINUE
         ELSE
            DO 160 I = 1, IN
               W( J ) = -WORK( I )
               J = J + 1
  160       CONTINUE
         END IF
         IBEGIN = IEND + 1
  170 CONTINUE
      M = N
*   
      RETURN
*
*     End of DLARRI
*
      END
