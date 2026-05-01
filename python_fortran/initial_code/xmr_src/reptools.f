*  ======================================================================
*   Convenience tools for working with our compressed representation
*   data structure, in particular to facilitate easy integration in
*   C++ Grailframe.
*
*   XMR_INITREP
*     to set up the representation given the primary data
*     explicitly (instead of assuming it is already set as
*     DLAXRR does it)
*
*   XMR_REPSIZE_REAL
*     given n, returns int giving absolute size of real part
*
*   XMR_REPSIZE_INT
*     given n, returns int giving absolute size of real part
*
*  ======================================================================
      SUBROUTINE XMR_INITREP( N, K, TYPE, G, OMEGA,
     $                        E, PIVBASE, REPR, REPI )
      IMPLICIT NONE
*     ..
*     .. Scalar arguments ..
*     ..
      DOUBLE PRECISION  PIVBASE
      INTEGER           N, K, TYPE
*     ..
*     .. Array arguments ..
*     ..
      DOUBLE PRECISION  G(N), E(1:N-1), REPR(4*N+3)
      INTEGER           OMEGA(N), REPI(6+N+N/2)
*
*  Purpose
*  =======
*
*    Set up representation data structure from the data given by
*    N, K, G, Omega and E.
*
*    This is a convenience routine, where no part of REPR or REPI
*    is assumed to ba initialized already.
*
*    Note that OMEGA is indexed 1:n, no padding here.
*
*  ======================================================================
*
*     .. Declarations ..
*
      EXTERNAL DLAXRR
*
*     .. Locals ..
*
      INTEGER IXG, IYOMGA, I
*
*  ----- Executable Statements -----------------------------------------
*
      IXG    = (0)
      IYOMGA = (4)

      DO I = 1, N
         REPR(IXG + I)    = G(I)
         REPI(IYOMGA + I) = OMEGA(I)
      ENDDO

      CALL DLAXRR( N, K, TYPE, E, PIVBASE, REPR, REPI )

      END
*
*     End of subroutine XMR_INITREP
*
************************************************************************
*
      FUNCTION XMR_REPSIZE_REAL( N )
      IMPLICIT NONE
      INTEGER N, XMR_REPSIZE_REAL
      XMR_REPSIZE_REAL = (4*N+3)
      END
*
************************************************************************
*
      FUNCTION XMR_REPSIZE_INT( N )
      IMPLICIT NONE
      INTEGER N, XMR_REPSIZE_INT
      XMR_REPSIZE_INT = (6+N+N/2)
      END
*
************************************************************************


