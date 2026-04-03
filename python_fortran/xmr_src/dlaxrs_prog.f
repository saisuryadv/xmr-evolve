      SUBROUTINE DLAXRS_PROG(
     $             N, J1, J2, I0, I1, G, OMEGA, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, SHFPRT, OFF, GPLUS, OPLUS, P
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, J1, J2, I0, I1
      INTEGER,          INTENT(IN)  ::  OMEGA(0:N+1)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU, OFF
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
      DOUBLE PRECISION, INTENT(IN)  ::  SHFPRT(N)
*
      INTEGER,          INTENT(OUT)  ::  OPLUS(N)
      DOUBLE PRECISION, INTENT(OUT)  ::  GPLUS(I0:I1), P(N)
*
*  Purpose    !!! Progressive at the moment only non-blocked !!!
*  =======
*
*     Computes a progressive factorization from J1 to J2.
*
*     J1 < J2 and I0 = 1 and I1 = N-1
*        produces GPLUS(J1:J2-1), OPLUS(J1+1:J2) and P(J1+1:J2)
*
*     J1 > J2 and I0 = 2 and I1 = N
*        produces GPLUS(J2+1:J1), OPLUS(J2:J1-1) and P(J2:J1-1)
*
*     It is assumed the start index J1 is twist in the source, but not
*     in the target. Likewise, the end index J2 is assumed to be twist
*     in the target, but not the source.
*     This has consequences for OPLUS:
*     - OPLUS(J1) may be 1, ie a block may end there.
*       Then the auxiliary has been computed accordingly.
*       Will hit us as soon as we use twisted BFs, but should not
*       require much change.
*
*  ======================================================================


*
*     .. Locals ..
*
      DOUBLE PRECISION  AUX
      INTEGER           I, INEXT, J, DIR

*
*     -- Executable Statements -----------------------------------------
*

      I = J1
      J = J2

      IF( J1.LT.J2 .AND. I0.EQ.1 .AND. I1.EQ.N-1 )THEN
         DIR = +1
      ELSEIF( J2.LT.J1 .AND. I0.EQ.2 .AND. I1.EQ.N )THEN
         DIR = -1
      ELSE
         I   = -1
         DIR = 0
      ENDIF
      INEXT = I+DIR


      AUX = OFF + (G(I) - TAU*SHFPRT(I))
      DO
         GPLUS(I) = NGN(INEXT) + AUX
         IF( ABS(GPLUS(I)).LT.PIVMIN )  GPLUS(I) = -PIVMIN
         P(INEXT) = (G(INEXT)/GPLUS(I)) * AUX
         OPLUS(INEXT) = 0
         IF( INEXT.EQ.J )THEN
            EXIT
         ENDIF
         I     = INEXT
         INEXT = INEXT+DIR
         AUX = P(I) - TAU*SHFPRT(I)
      ENDDO

      END SUBROUTINE DLAXRS_PROG
*
************************************************************************
