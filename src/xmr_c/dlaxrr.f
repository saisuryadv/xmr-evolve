      SUBROUTINE DLAXRR( N, K, TYPE, E, PIVBASE, REPR, REPI )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, TYPE
      DOUBLE PRECISION, INTENT(IN)  ::  PIVBASE
      DOUBLE PRECISION, INTENT(IN)  ::  E(1:N-1)
*
      INTEGER,          INTENT(INOUT)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(INOUT)  ::  REPR(4*N+3)
*
*  Purpose
*  =======
*
*    Set up representation data structure.
*    Is called within dlarrf and for breadth-first traversal also
*    within dlaxrv.
*    Assumes that the following primary data elements have been set:
*      G(1:n), Omega(1:n).
*
*  ======================================================================
*
*     .. Declarations ..
*
      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH
*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0
*
*     .. Locals ..
*
      DOUBLE PRECISION  UNDEF, BDET
      INTEGER           IXG, IXBDET, IXNGN, IXGNSQ, IYOMGA, IYLBBK
      INTEGER           I, IBB, NB


*
*  ----- Executable Statements -----------------------------------------
*
      IXG    = (0)
      IXBDET = (N)
      IXNGN  = (2*N+1)
      IXGNSQ = (3*N+2)
      IYOMGA = (4)
      IYLBBK = (N+6)

      UNDEF = DLAMCH('Overflow')

*
*     Analyse block structure
*     set squares of the offdiagonals
*
      IBB = 0
      DO I = 1, K-1
         REPR(IXGNSQ + I) = E(I)**2
         IF( REPI(IYOMGA + I) .NE. 0 )THEN
            REPI(IYLBBK + IBB) = I-1
            IBB = IBB+1
         ENDIF
      ENDDO
      IF( K .EQ. N )THEN
         IF( REPI(IYOMGA + K) .NE. 0 )THEN
            REPI(IYLBBK + IBB) = I-1
            IBB = IBB+1
         ENDIF
      ENDIF

      REPI(IYLBBK + IBB) = K
      REPR(IXGNSQ + K)   = UNDEF
      IBB = IBB+1

      IF( K .EQ. 1 )THEN
         IF( REPI(IYOMGA + K) .NE. 0 )THEN
            REPI(IYLBBK + IBB) = I+1
            IBB = IBB+1
         ENDIF
      ENDIF
      DO I = K+1, N
         REPR(IXGNSQ + I) = E(I-1)**2
         IF( REPI(IYOMGA + I) .NE. 0 )THEN
            REPI(IYLBBK + IBB) = I+1
            IBB = IBB+1
         ENDIF
      ENDDO
      NB = IBB-1

      DO I = IBB, N/2
         REPI(IYLBBK + I) = 0
      ENDDO
*
*     Upper part (analogous to below)
*
      I = 1
      IBB = 0
      DO
         DO
            IF( I .EQ. REPI(IYLBBK + IBB) )THEN
               EXIT
            ENDIF
            REPR(IXNGN + I)  = REPR(IXGNSQ + I) / REPR(IXG + I)
            REPR(IXBDET + I) = ZERO
            I = I+1
         ENDDO
         IF( I .EQ. K )THEN
            EXIT
         ENDIF

         IBB = IBB+1

         BDET = REPR(IXG + I) * REPR(IXG + I+1)  -  REPR(IXGNSQ + I)
         REPR(IXNGN  + I) = ZERO
         REPR(IXBDET + I) = BDET

         IF( I .EQ. N-1 )THEN
*           I = N-1 implies K = N
            EXIT
         ENDIF

         REPR(IXNGN  + I+1) =
     $      (REPR(IXGNSQ + I+1) * REPR(IXG + I))  /  BDET
         REPR(IXBDET + I+1) = ZERO

         I = I+2
      ENDDO
*
*     The Twist
*
      REPR(IXNGN + K)  = UNDEF
      REPR(IXBDET + K) = ZERO
*
*     Lower part (analogous to above)
*
      I = N
      IBB = NB
      DO
         DO
            IF( I .EQ. REPI(IYLBBK + IBB) )THEN
               EXIT
            ENDIF
            REPR(IXNGN + I)  = REPR(IXGNSQ + I) / REPR(IXG + I)
            REPR(IXBDET + I) = ZERO
            I = I-1
         ENDDO
         IF( I .EQ. K )THEN
            EXIT
         ENDIF

         IBB = IBB-1

         BDET = REPR(IXG + I) * REPR(IXG + I-1)  -  REPR(IXGNSQ + I)
         REPR(IXNGN  + I) = ZERO
         REPR(IXBDET + I) = BDET

         IF( I .EQ. 2 )THEN
*           I = 2 implies K = 1
            EXIT
         ENDIF

         REPR(IXNGN  + I-1) =
     $      (REPR(IXGNSQ + I-1) * REPR(IXG + I))  /  BDET
         REPR(IXBDET + I-1) = ZERO

         I = I-2
      ENDDO

      REPI(1)  = TYPE
      REPI(3)    = NB
      REPI(2)     = K
      REPI(IYOMGA)     = 0
      REPI(IYOMGA+N+1) = 0
      REPR(IXNGN)      = ZERO
      REPR(IXNGN+N+1)  = ZERO
      REPR(4*N+3)   = PIVBASE



      END SUBROUTINE DLAXRR
*
************************************************************************
