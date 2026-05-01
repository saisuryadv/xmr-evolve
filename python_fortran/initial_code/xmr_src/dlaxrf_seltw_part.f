      SUBROUTINE DLAXRF_SELTW_PART(
     $             N, IA, IE, DIR, G, GN, OMEGA, TWISTOK, S,
     $             ARCSQ, ABETA, AGMAX, AGSMAX, ABRCFAC
     $           )
      IMPLICIT NONE
*
      INTEGER, INTENT(IN)  ::  N, IA, IE, DIR
      INTEGER, INTENT(IN)  ::
     $   TWISTOK( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ),
     $   OMEGA(   MIN(IA,IE+DIR) : MAX(IA,IE+DIR) )
      DOUBLE PRECISION, INTENT(IN) ::
     $   G(  MIN(IA,IE) : MAX(IA,IE) ),
     $   GN( MIN(IA,IE) : MAX(IA,IE) ),
     $   S( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) )
*
      DOUBLE PRECISION, INTENT(INOUT)  ::
     $   ARCSQ ( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ),
     $   ABETA ( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ),
     $   AGMAX ( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ),
     $   AGSMAX( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ),
     $   ABRCFAC( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) )
*
*  Purpose
*  =======
*
*     Compute auxiliary quantities to evaluate element growth and
*     relative condition for a bidiagonal factorization, depending
*     on IA and IE either
*      (L)  a top-to-bottom one if 1 = IA <= IE < N, DIR=+1
*      (U)  a bottom-to-top one if N = IA >= IE > 1, DIR=-1.
*     For the following assume case (L), that is, G = D and GN = LD for
*     an LDL (partial) factor.
*
*     The quantity of main interest is an estimate for the relative
*     condition number
*       RC := || X L' s || * || invX invL s ||
*     for a nonsingular diagonal scaling matrix X which minimizes it.
*     The latter is equivalent to saying
*       RC := ||r||^2 with r(i) := sqrt[ (L'*s)(i) * (invL*s)(i) ]
*
*     The result arrays are changed as follows by this routine (again
*     only for the LDL case):
*        ARCSQ(i) += ||r(1:i-1)||^2,     i = 2:IE+1.
*        ABETA(i) += (invL*s)(i) - s(i), i = 2:IE+1.
*        AGMAX(i)  is maxed to Max{ |D(i)| : 1 <= j < i }
*        AGSMAX(i) is maxed to Max{ |D(i)*s(i)| : 1 <= j < i }.
*        ABRCFAC(i) gives an additional factor to multiply the final
*          RC estimate for twist i with. This is to account for the
*          fact that the standard RC-formula needs to be modified
*          if blocks are present. The entries are maxed as well.
*
*       CC   ARCSQ_NOS is analogous, except that we bound RC by
*       CC      || X *L' || * || invX * invL ||.
*       CC    cancelled for now: need separate beta return as well
*
*     For all cases, the entry 1 is unchanged, as the corresponding
*     quantity would be zero. Besides that, only entries where a twist
*     is allowed are touched (TWISTOK(i)!=0), this should include all
*     indices where a no (non-final) block ends.
*
*     Notes.
*     a)
*         (L'*s)(i) = s(i) +
*            k(i)*s(i+2), if i starts a block and i != n-1
*            l(i)*s(i+1), otherwise and i < n,
*            0, i=n or i=n-1 and i starts a block.
*     b)
*         (invL*s)(i) = s(i) +
*            0,   if i=1 or i ends a block,
*           -l(i-1) * (invL*s)(i-1),   if i-1 does not end a block,
*           -k(i-2) * (invL*s)(i-2) - l(i-1)*s(i-1),   otherwise.
*
*  ======================================================================
*
*     .. Declarations ..
*
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
*
*     .. Local Variables ..
*
      DOUBLE PRECISION  BETA, RCSQ, GMAX, GSMAX, BRCFAC
      DOUBLE PRECISION  INVNS, TRPNS, F, Z
      INTEGER           I, J, INEXT, IPREV
*
*
*  ===== Executable Statements ==========================================
*
*
      I = IA
      J = IE+DIR
*
      BETA = ZERO
      RCSQ = ZERO
      GMAX = ZERO
      GSMAX = ZERO
      DO
*        Here we have:
*          I != J, result quantities for I have been set and are
*          given by the corresponding locals.
*
         INEXT = I+DIR
*
         INVNS = BETA + S(I)
         GMAX  = MAX( GMAX, ABS(G(I)) )
         GSMAX = MAX( GSMAX, ABS(G(I)*S(I)) )
*
         IF( OMEGA(INEXT).EQ.0 )THEN
            F = GN(I) / G(I)
            TRPNS = S(I) + F*S(INEXT)

            RCSQ  = RCSQ + ABS(INVNS*TRPNS)
            BETA  = - F*INVNS

         ELSEIF( INEXT.EQ.1 .OR. INEXT.EQ.N )THEN

            RCSQ = RCSQ + ABS(INVNS*S(I))
            BETA = ZERO

         ELSE
            Z = GN(INEXT) / (G(I)*G(INEXT) - GN(I)**2)
            F = -GN(I) * Z
*           now F=k(i)
            TRPNS = S(I) + F*S(INEXT+DIR)
            RCSQ  = RCSQ + ABS(INVNS*TRPNS)
            GMAX  = MAX( GMAX, ABS(G(I)) )
            GSMAX = MAX( GSMAX, ABS(G(I)*S(I)) )
            BETA  = - F*INVNS

            IPREV = I
            I = INEXT
            INEXT = INEXT+DIR

            F = G(IPREV)*Z
*           now F = N(iprev)
            TRPNS = S(I) + F*S(INEXT)
            RCSQ  = RCSQ + ABS(S(I)*TRPNS)
            BETA  = BETA - F*S(I)

         ENDIF

         IF( TWISTOK(INEXT) .NE. 0 )THEN
            ARCSQ(INEXT) = ARCSQ(INEXT) + RCSQ
            ABETA(INEXT) = ABETA(INEXT) + BETA
            AGMAX(INEXT) = MAX( AGMAX(INEXT),  GMAX )
            AGSMAX(INEXT) = MAX( AGSMAX(INEXT), GSMAX )
         ENDIF

         IF( INEXT .EQ. J )THEN
            EXIT
         ENDIF

         I = INEXT
         INEXT = INEXT+DIR

      ENDDO
*
      BRCFAC = ONE
      I = IA
      DO
         ABRCFAC(I) = MAX( ABRCFAC(I), BRCFAC )
         IF( I .EQ. IE+DIR ) EXIT
         IF( OMEGA(I) .NE. 0 )THEN
            BRCFAC = MAX( BRCFAC,
     $                    ABS(G(I-DIR) / GN(I-DIR)),
     $                    ABS(G(I) / GN(I-DIR)) )
         ENDIF
         I = I + DIR
      ENDDO
*
      END SUBROUTINE DLAXRF_SELTW_PART
*
*************************************************************************
