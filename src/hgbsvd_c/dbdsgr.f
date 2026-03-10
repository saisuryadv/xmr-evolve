      SUBROUTINE DBDSGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL,
     $                   M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
     $                   LIWORK, INFO )
*
*  -- LAPACK computational routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE
      INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            ISUPPZ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )
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
*  DBDSGR computes selected singular values and, optionally, singular vectors
*  of a real symmetric bidiagonal matrix B. Singular values and
*  singular vectors can be selected by specifying either a range of values
*  or a range of indices for the desired singular values. The singular values
*  are computed by the dqds algorithm, while orthogonal eigenvectors are
*  computed from various ``good'' L D L^T representations (also known as
*  Relatively Robust Representations). Refer to DSTEGR for more information.
*  
*  For more details, see "A new O(n^2) algorithm for the symmetric
*  tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,
*  Computer Science Division Technical Report No. UCB/CSD-97-971,
*  UC Berkeley, May 1997.
*
*  The extension to the SVD is performed using a set of coupling
*  transformations relating the L D L^T representations of B^T B, B B^T 
*  and the so-called Golub-Kahan matrix respectively.
*  The technique is described in "An O(n^2) algorithm for the bidiagonal 
*  SVD", by Benedikt Grosser and Bruno Lang, to appear in LAA, 2001.
*
*  Note 1 : Currently DBDSGR is only set up to find ALL the n
*  singular values and corresponding singular vectors of B in O(n^2) time.
*  Note 2 : DSTEGR works only on machines which follow ieee-754
*  floating-point standard in their handling of infinities and NaNs.
*  Normal execution of DSTEGR may create NaNs and infinities and hence
*  may abort due to a floating point exception in environments which
*  do not conform to the ieee standard.
*  Note 3 : Only couplings for so-called positive definite initial matrices
*  are implemented. Thus the routine may return no results for very
*  tight clusters of singular values. This is flagged by INFO=3, 4 or 5.
c
c  Reminder to Osni:
c  Note 4 : Splittings are currently not supported, flagged by INFO = 33!
c
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute singular values only;
*          = 'L':  Compute singular values and left singular vectors.
*          = 'R':  Compute singular values and right singular vectors.
*          = 'A':  Compute singular values and both left and right  
*		   singular vectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all singular values will be found.
*          = 'V': all singular values in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th singular values will be found.
********** Only RANGE = 'A' is currently supported *********************
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the N diagonal elements of the bidiagonal matrix
*          B. On exit, D is overwritten.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the (N-1) offdiagonal elements of the bidiagonal
*          matrix B in elements 1 to N-1 of E. E(N) need not be set on
*          input, but is used internally as workspace.
*          On exit, E is overwritten.
*
*  VL      (input) DOUBLE PRECISION
*  VU      (input) DOUBLE PRECISION
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for singular values. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in decending order) of the
*          smallest and largest singular values to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (input) DOUBLE PRECISION
*          The absolute error tolerance for the
*          singular values/singular vectors. IF JOBZ = 'A', the singular 
*	   values and singular vectors output have residual norms bounded  
*	   ABSTOL, and the dot products between different singular vectors 
*          are bounded by ABSTOL. For JOBZ = 'L' or JOBZ = 'R' the dot 
*	   products between different singular vectors are bounded by ABSTOL,
*          whereas residual norms cannot be computed.
*	   Let T denote the Golub-Kahan matrix
*	   T =   diag([D(1),E(1),...,E(N-1),D(N)],-1)
*              + diag([D(1),E(1),...,E(N-1),D(N)],1)
*                                                          [  0  B ]
*          resulting from a perfect shuffle permutation of [       ] . 
*                                                          [ B^T 0 ]
*	   If ABSTOL is less than N*EPS*|T|, then
*          N*EPS*|T| will be used in its place, where EPS is the
*          machine precision and |T| is the 1-norm of the tridiagonal
*          matrix. The singular values are computed to an accuracy of
*          EPS*|T| irrespective of ABSTOL. If high relative accuracy
*          is important, set ABSTOL to DLAMCH( 'Safe minimum' ).
*          See Barlow and Demmel "Computing Accurate Eigensystems of
*          Scaled Diagonally Dominant Matrices", LAPACK Working Note #7
*          for a discussion of which matrices define their singular values
*          to high relative accuracy.
*
*  M       (output) INTEGER
*          The total number of singular values found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          The first M elements contain the selected singular values in
*          descending order.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) )
*          If JOBZ = 'L' or JOBZ = 'R', then if INFO = 0, the first M 
*	   columns of Z contain the orthonormal 
*          left (JOBZ = 'L') or right (JOBZ = 'R') singular vectors 
*          of the matrix  corresponding to the selected singular values, 
*          with the i-th column of Z holding the singular vector 
*          associated with W(i).
*          If JOBZ = 'A', then if INFO = 0, the first M 
*	   columns of Z contain the orthonormal 
*          right singular vectors (internal matrix V) of the matrix B 
*          in rows 1 to N,
*	   whereas rows LDZ/2+1 to LDZ/2+N contain the orthonormal 
*	   left singular vectors (internal matrix U).
*          If JOBZ = 'N', then Z is not referenced.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and an upper bound must be used.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'L' or JOBZ = 'R', LDZ >= max(1,N).
*	   If JOBZ = 'A', LDZ >= max(1,2*N)
*          Storage convention for JOBZ = 'A', LDZ = 16, M = 4, N = 6:
*		row 1       VVVV::
*		            VVVV::
*		            VVVV::
*		            VVVV::
*		            VVVV::
*		row N       VVVV::
*		            ::::::
*		row LDZ/2   ::::::
*		row LDZ/2+1 UUUU::
*		            UUUU::
*		            UUUU::
*		            UUUU::
*		            UUUU::
*		row LDZ/2+N UUUU::
*		            ::::::
*		row LDZ     ::::::
*
*  ISUPPZ  (output) INTEGER ARRAY, 
*	   If JOBZ = 'L' or 'R', dimension ( 2*max(1,LDZ) )
*	   If JOBZ = 'A', dimension ( 2*max(1,LDZ) )
*          The support of the singular vectors in Z, i.e., the indices
*          indicating the nonzero elements in Z. If JOBZ = 'L' or 'R', 
*	   the i-th singular vector is nonzero only in elements
*          ISUPPZ( 2*i-1 ) through ISUPPZ( 2*i ).
*          If JOBZ = 'A', the i-th right singular vector is nonzero only
*	   in elements ISUPPZ( 2*i-1 ) through ISUPPZ( 2*i ), and the
*	   i-th left singular vector is nonzero onlyin elements 
*	   ISUPPZ( LDZ + 2*i-1 ) through ISUPPZ( LDZ + 2*i ).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal
*          (and minimal) LWORK.
*
c
c  Reminder to Osni:
c  DBDSGR requires an workspace array of dimension 26*N .
c
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,26*N)
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.  LIWORK >= max(1,10*N)
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = 1, internal error in DLARRI,
*                if INFO > 2, internal error in DLARRV.
c
c	   Reminder to Osni: INFO = 33 - Splittings are not supported.
c
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
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LQUERY, VALEIG, WANTZ
c     new LOGICAL variables
      LOGICAL            WANTA, WANTL, WANTR
      INTEGER            I, IBEGIN, IEND, IINDBL, IINDWK, IINFO, IINSPL, 
     $                   INDGRS, INDWRK, ITMP, J, JJ, LIWMIN, LWMIN, 
     $                   NSPLIT
c     new INTEGER variables
      INTEGER 		 IOFFSU, IOFFSV
      INTEGER 		 INDAOR, INDBOR, INDD2, INDE2, INDDE0, INDDE1
      INTEGER 		 INDDC, INDLC, INDSGL, INDSGR
c     changes for DOUBLE PRECISION variables
      DOUBLE PRECISION   EPS, SCALE,
     $                   THRESH, TMP, TNRM, TOL, T1, T2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DSECND
      EXTERNAL           DLAMCH, DLANST, DSECND, LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARRI, DLARRV, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTL = LSAME( JOBZ, 'L' )
      WANTR = LSAME( JOBZ, 'R' )
      WANTA = LSAME( JOBZ, 'A' )
      WANTZ = ( WANTL .OR. WANTR .OR. WANTA )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
*
      LQUERY = ( ( LWORK.EQ.-1 ) .OR. ( LIWORK.EQ.-1 ) )
      LWMIN = 26*N
      LIWMIN = 10*N
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
*
*     The following two lines need to be removed once the
*     RANGE = 'V' and RANGE = 'I' options are provided.
*
      ELSE IF( VALEIG .OR. INDEIG ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) THEN
         INFO = -7
      ELSE IF( INDEIG .AND. ( IL.LT.1 .OR. IL.GT.N ) ) THEN
         INFO = -8
      ELSE IF( INDEIG .AND. ( IU.LT.IL .OR. IU.GT.N ) ) THEN
         INFO = -9
      ELSE IF( LDZ.LT.1 .OR. ( ( WANTZ .AND. LDZ.LT.N ) 
     $   .OR. (WANTA. AND. LDZ.LT.2*N) ) ) THEN
         INFO = -14
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -17
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -19
      END IF
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DBDSGR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      IF ( WANTA ) THEN
	IOFFSV = 0
	IOFFSU = LDZ/2
      ELSE IF ( WANTL ) THEN
	IOFFSU = 0
	IOFFSV = 0
      ELSE
	IOFFSV = 0
	IOFFSU = 0
      ENDIF
*
*     Quick return if possible
*
      M = 0
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = D( 1 )
         ELSE
            IF( VL.LT.D( 1 ) .AND. VU.GE.D( 1 ) ) THEN
               M = 1
               W( 1 ) = D( 1 )
            END IF
         END IF
         IF( WANTL.OR.WANTA ) THEN
           Z( IOFFSU+1, 1 ) = ONE
	 ENDIF
         IF( WANTR.OR.WANTA ) THEN
           Z( IOFFSV+1, 1 ) = ONE
	 ENDIF
         RETURN
      END IF
*
*     Get machine constants.
*
      EPS = DLAMCH( 'Precision' )
*
*     Scale matrix.
*
      TNRM = DLANST( 'M', N, D, E )
      SCALE  = ONE/TNRM
      IF( SCALE.NE.ONE ) THEN
         CALL DSCAL( N, SCALE, D, 1 )
         CALL DSCAL( N-1, SCALE, E, 1 )
         TNRM = TNRM*SCALE
      END IF
*
      T1 = DSECND( )
*
*     Set pointers.
*
      INDGRS =  1
      INDSGL =  2*N + 1
      INDSGR =  3*N + 1
      INDAOR =  4*N + 1
      INDBOR =  5*N + 1
      INDD2  =  6*N + 1
      INDE2  =  7*N + 1
      INDDE0 =  8*N + 1
      INDDE1 =  9*N + 1
      INDDC  = 10*N + 1
      INDLC  = 11*N + 1
      INDWRK = 12*N + 1
*
      IINSPL = 1
      IINDBL = N + 1
      IINDWK = 2*N + 1
*
*     Compute signs of B's elements.
*
      WORK( INDSGR ) = ONE
      DO 5, I = 1,N-1
	WORK( INDSGL-1+I ) = SIGN(ONE,D(I))
	D(I) = WORK( INDSGL-1+I )*D(I)
	E(I) = WORK( INDSGL-1+I )*E(I)
	WORK( INDSGR+I ) = SIGN(ONE,E(I))
	E(I) = WORK( INDSGR+I )*E(I)
	D(I+1) = WORK( INDSGR+I )*D(I+1)
  5   CONTINUE
      WORK( INDSGL-1+N ) = SIGN(ONE,D(N))
      D(N) = WORK( INDSGL-1+I )*D(N)
*
      CALL DCOPY(N,D,1,WORK(INDAOR),1)
      CALL DCOPY(N-1,E,1,WORK(INDBOR),1)
*
*     Compute the desired singular values of the bidiagonal B.
*     Form correponding Gerschgorin intervals  WORK( INDGRS ),
*     and some auxiliary variables WORK( INDD2 ), WORK( INDE2 ),  
*     WORK( INDDE0 ), WORK( INDDE1 ) needed later.
*
c
c     Reminder to Osni:
c     New routine DLARRI instead of DLARRE. 
c     Since splittings are not supported, this routine may return IINFO=33
c
      THRESH = EPS*TNRM
      CALL DLARRI( N, WANTL, WANTR, WANTA, D, E, THRESH, NSPLIT, 
     $		     IWORK( IINSPL ), M, W,
     $               WORK( INDGRS ), 
     $               WORK( INDD2 ), WORK( INDE2 ),  
     $		     WORK( INDDE0 ), WORK( INDDE1 ),
     $		     WORK( INDWRK ), 
     $               IINFO )
      T2 = DSECND( )
      TIMNG( 7 ) = TIMNG( 7 ) + ( T2-T1 )
	write(*,*) TIMNG( 7 )
      IF( IINFO.NE.0 ) THEN
         INFO = IINFO
         RETURN
      END IF
*
      IF( WANTZ ) THEN
*
*        Compute the desired singular vectors corresponding to the computed
*        singular values
*
         IBEGIN = 1
         DO 20 I = 1, NSPLIT
            IEND = IWORK( IINSPL+I-1 )
            DO 10 J = IBEGIN, IEND
               IWORK( IINDBL+J-1 ) = I
   10       CONTINUE
            IBEGIN = IEND + 1
   20    CONTINUE
*
         TOL = MAX( ABSTOL, DBLE( N )*EPS )
         T1 = DSECND( )
c
c        Reminder to Osni: Interface changed.
c        Since DLARRV needs WORK( INDWRK:INDWRK+14*N-1 ) as workspace
c	 and INDWRK = 12*N + 1, the overall workspace requirement
c	 sums up to 26*N .
c
         CALL DLARRV( N, D, E, IWORK( IINSPL ), M, W, IWORK( IINDBL ),
     $                WORK( INDGRS ), TOL, WANTA, Z, LDZ, ISUPPZ,
     $                WORK( INDAOR ), WORK( INDBOR ),  
     $                WORK( INDD2 ), WORK( INDE2 ),  
     $		      WORK( INDDE0 ), WORK( INDDE1 ),
     $		      WORK( INDDC ), WORK( INDLC ),
     $                WORK( INDWRK ), IWORK( IINDWK ), IINFO )
         T2 = DSECND( )
         TIMNG( 8 ) = TIMNG( 8 ) + ( T2-T1 )
         IF( IINFO.NE.0 ) THEN
            INFO = 2
            RETURN
         END IF
*
      END IF
*
      IBEGIN = 1
      DO 40 I = 1, NSPLIT
         IEND = IWORK( IINSPL+I-1 )
         DO 30 J = IBEGIN, IEND
	    W( J ) = SQRT( W(J) )
   30    CONTINUE
         IBEGIN = IEND + 1
   40 CONTINUE
*
*     If matrix was scaled, then rescale singular values appropriately.
*
      IF( SCALE.NE.ONE ) THEN
         CALL DSCAL( M, ONE / SCALE, W, 1 )
      END IF
*
*     Arrange signs.
*
      DO 45, I = 1,N
        IF ((WORK( INDSGL-1+I ).NE.ONE).AND.(WANTL.OR.WANTA)) THEN
	  CALL DSCAL(N,-ONE,Z(IOFFSU+I,1),LDZ)
	ENDIF
        IF ((WORK( INDSGR-1+I ).NE.ONE).AND.(WANTR.OR.WANTA)) THEN
	  CALL DSCAL(N,-ONE,Z(IOFFSV+I,1),LDZ)
	ENDIF
  45    CONTINUE
*
*     Sort singular values to descending order, along with
*     singular vectors.
*
         DO 60 J = 1, M - 1
            I = 0
            TMP = W( J )
            DO 50 JJ = J + 1, M
               IF( W( JJ ).GT.TMP ) THEN
                  I = JJ
                  TMP = W( JJ )
               END IF
   50       CONTINUE
            IF( I.NE.0 ) THEN
               W( I ) = W( J )
               W( J ) = TMP
               IF( WANTZ ) THEN
                  CALL DSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
                  ITMP = ISUPPZ( 2*I-1 )
                  ISUPPZ( 2*I-1 ) = ISUPPZ( 2*J-1 )
                  ISUPPZ( 2*J-1 ) = ITMP
                  ITMP = ISUPPZ( 2*I )
                  ISUPPZ( 2*I ) = ISUPPZ( 2*J )
                  ISUPPZ( 2*J ) = ITMP
		  IF ( WANTA ) THEN
*
*		    Assumes that left singular vectors
*		    are stored in Z(IOFFSU+1:IOFFSU+N,1:N)
*
                    CALL DSWAP( N, Z( IOFFSU+1, I ), 1, 
     $		      Z( IOFFSU+1, J ), 1 )
                    ITMP = ISUPPZ( IOFFSU+2*I-1 )
                    ISUPPZ( IOFFSU+2*I-1 ) = ISUPPZ( IOFFSU+2*J-1 )
                    ISUPPZ( IOFFSU+2*J-1 ) = ITMP
                    ITMP = ISUPPZ( IOFFSU+2*I )
                    ISUPPZ( IOFFSU+2*I ) = ISUPPZ( IOFFSU+2*J )
                    ISUPPZ( IOFFSU+2*J ) = ITMP
		  ENDIF
               END IF
            END IF
   60    CONTINUE
*
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of DBDSGR
*
      END
