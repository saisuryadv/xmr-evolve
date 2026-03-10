*******************************************************************************
*
      SUBROUTINE DLATS1( Z, ZTZ, L, U, LD, B1, R, BN, ISUPPZ,
     $	SAWNAN )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
c
c  Reminder to Osni:
c  extracted from former dlar1v.f
c
*     .. Scalar Arguments ..
      INTEGER B1, R, BN
      LOGICAL SAWNAN
      DOUBLE PRECISION ZTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION Z( * ), L( * ), U( * ), LD( * )
      INTEGER ISUPPZ( * )
*
*  Purpose
*  =======
*
*  Solves N G N^T * Z = G(R,R) * E_R .
*  G is diagonal, whereas N is a twisted matrix:
*    N = eye(n) + diag([L(1),...,L(R-1),0,...,0])
*               + diag([0,...,0,U(R),...,U(n-1)])
*  The right hand side is given by G(R,R) times the R-th unit vector.
*  One proceeds as follows:
*    Z(R) = 1
*    FOR I = R-1:-1:B1
*      Z(I) = -( L( I )*Z( I+1 ) )  
*    END
*    FOR I = R+1:BN
*      Z( I ) = -( U( I-1 )*Z( I-1 ) )
*    END
*  Thus only L and U are used instead of the symbolic matrices
*  G and N.
*
*  DLATS1 is ashortcut for DLA + T(wisted matrix) S(olver) 1 .
*
*  Arguments
*  =========
*
*  Z       (output) DOUBLE PRECISION array, dimension (N)
*          The solution of the linear system, typically
*          a very good eigenvector approximation.
* 
*  ZTZ     (output) DOUBLE PRECISION
*          Dot product ZTZ = Z^T Z .
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
*  ISUPPZ  (output) INTEGER array
*          Keeps the support information. The smallest and largest 
*          indices of the non-zero elements of Z are stored in ISUPPZ(1) 
*          and ISUPPZ(2) resp.
*
*  SAWNAN  (input) LOGICAL
*          Indicates if a breakdown occured.
*
*     .. Parameters ..
      INTEGER            BLKSIZ
      PARAMETER          ( BLKSIZ = 32 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. External functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Local Scalars ..
      INTEGER            FROM, I, TO
      DOUBLE PRECISION   EPS
*
      EPS = DLAMCH( 'Precision' )
      ISUPPZ( 1 ) = B1
      ISUPPZ( 2 ) = BN
      Z( R ) = ONE
      ZTZ = ONE
      IF( .NOT.SAWNAN ) THEN
         FROM = R - 1
         TO = MAX( R-BLKSIZ, B1 )
  120    CONTINUE
         IF( FROM.GE.B1 ) THEN
            DO 130 I = FROM, TO, -1
               Z( I ) = -( L( I )*Z( I+1 ) )
               ZTZ = ZTZ + Z( I )*Z( I )
  130       CONTINUE
            IF( ABS( Z( TO ) ).LE.EPS .AND. ABS( Z( TO+1 ) ).LE.EPS )
     $           THEN
               ISUPPZ( 1 ) = TO + 2
            ELSE
               FROM = TO - 1
               TO = MAX( TO-BLKSIZ, B1 )
               GO TO 120
            END IF
         END IF
         FROM = R + 1
         TO = MIN( R+BLKSIZ, BN )
  140    CONTINUE
         IF( FROM.LE.BN ) THEN
            DO 150 I = FROM, TO
               Z( I ) = -( U( I-1 )*Z( I-1 ) )
               ZTZ = ZTZ + Z( I )*Z( I )
  150       CONTINUE
            IF( ABS( Z( TO ) ).LE.EPS .AND. ABS( Z( TO-1 ) ).LE.EPS )
     $           THEN
               ISUPPZ( 2 ) = TO - 2
            ELSE
               FROM = TO + 1
               TO = MIN( TO+BLKSIZ, BN )
               GO TO 140
            END IF
         END IF
      ELSE
         DO 160 I = R - 1, B1, -1
            IF( Z( I+1 ).EQ.ZERO ) THEN
               Z( I ) = -( LD( I+1 ) / LD( I ) )*Z( I+2 )
            ELSE
               Z( I ) = -( L( I )*Z( I+1 ) )
            END IF
            IF( ABS( Z( I ) ).LE.EPS .AND. ABS( Z( I+1 ) ).LE.
     $               EPS ) THEN
               ISUPPZ( 1 ) = I + 2
               GO TO 170
            END IF
            ZTZ = ZTZ + Z( I )*Z( I )
  160    CONTINUE
  170    CONTINUE
         DO 180 I = R, BN - 1
            IF( Z( I ).EQ.ZERO ) THEN
               Z( I+1 ) = -( LD( I-1 ) / LD( I ) )*Z( I-1 )
            ELSE
               Z( I+1 ) = -( U( I )*Z( I ) )
            END IF
            IF( ABS( Z( I ) ).LE.EPS .AND. ABS( Z( I+1 ) ).LE.EPS )
     $                THEN
               ISUPPZ( 2 ) = I - 1
               GO TO 190
            END IF
            ZTZ = ZTZ + Z( I+1 )*Z( I+1 )
  180    CONTINUE
  190    CONTINUE
      END IF
*
      RETURN
*
*     END OF DLATS1
*
      END
*
*******************************************************************************
