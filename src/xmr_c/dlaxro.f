      SUBROUTINE DLAXRO( N, M, W, Z, LDZ, ISUPPZ, REVORD )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N, M, LDZ
*
      INTEGER,          INTENT(INOUT)  ::  ISUPPZ(2*M)
      DOUBLE PRECISION, INTENT(INOUT)  ::  W(M), Z(LDZ,M)
*
      INTEGER,          INTENT(OUT)    ::  REVORD(M)
*
*
*  Purpose
*  =======
*
*    Sort the triplets ( W(i), Z(:,i), ISUPPZ(2*i-1:2*i) ) into
*    ascending order with respect to W.
*
*    The array REVORD is set to the reverse permutation: The eigenpairs
*    i in the output order was the REVORD(i)'th in the input order.
*
*  ======================================================================
*
*     .. Delcarations ..
*
      EXTERNAL          DSWAP
*
*     .. Locals ..
*
      DOUBLE PRECISION  RTMP
      INTEGER           I, J, K, ITMP
*
*  ===== Executable Statements ==========================================
*
      DO I = 1, M
         REVORD(I) = I
      ENDDO
*
      DO I = 1, M-1
         RTMP = W(I)
         J = I
         DO K = I+1, M
            IF( W(K) .LT. RTMP )THEN
               RTMP = W(K)
               J = K
            ENDIF
         ENDDO
*
         IF( J .NE. I )THEN

            CALL DSWAP( N, Z(1,I), 1, Z(1,J), 1 )

            RTMP = W(I)
            W(I) = W(J)
            W(J) = RTMP

            ITMP = ISUPPZ(2*I-1)
            ISUPPZ(2*I-1) = ISUPPZ(2*J-1)
            ISUPPZ(2*J-1) = ITMP

            ITMP = ISUPPZ(2*I)
            ISUPPZ(2*I) = ISUPPZ(2*J)
            ISUPPZ(2*J) = ITMP

            ITMP = REVORD(I)
            REVORD(I) = REVORD(J)
            REVORD(J) = ITMP

         ENDIF
      ENDDO
      END SUBROUTINE DLAXRO
*
*************************************************************************


