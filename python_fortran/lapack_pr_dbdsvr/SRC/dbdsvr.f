*> \brief \b DBDSVR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DBDSVR( UPLO, JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
*                          NS, S, Z, LDZ, WORK, LWORK, IWORK, LIWORK,
*                          INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, RANGE, UPLO
*       INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, N, NS
*       DOUBLE PRECISION   VL, VU
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   D( * ), E( * ), S( * ), WORK( * ),
*                          Z( LDZ, * )
*       ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>  DBDSVR computes the singular value decomposition (SVD) of a real
*>  N-by-N (upper or lower) bidiagonal matrix B, B = U * S * VT, using
*>  the Tridiagonal Golub-Kahan (TGK) MR^3 driver DBDSVDMR3.  Its
*>  argument list mirrors DBDSVDX so that DBDSVR is a drop-in
*>  alternative alongside DBDSVDX for callers that want the MR^3
*>  variant.
*>
*>  Given an upper bidiagonal B with diagonal D and superdiagonal E,
*>  DBDSVR builds the 2N-by-2N TGK tridiagonal, computes its
*>  eigendecomposition by the MR^3 algorithm rooted at the TGK
*>  (see P. Willems and B. Lang, SIAM J. Sci. Comput., 35:740-766,
*>  2013), and returns the singular values and (optionally) singular
*>  vectors of B in the same packing convention as DBDSVDX.
*>
*>  RANGE = 'V' or 'I' is realized by post-filtering the full spectrum
*>  from DBDSVDMR3 (which always returns all N singular values).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  B is upper bidiagonal;
*>          = 'L':  B is lower bidiagonal.
*> \endverbatim
*>
*> \param[in] JOBZ
*> \verbatim
*>          JOBZ is CHARACTER*1
*>          = 'N':  Compute singular values only;
*>          = 'V':  Compute singular values and singular vectors.
*> \endverbatim
*>
*> \param[in] RANGE
*> \verbatim
*>          RANGE is CHARACTER*1
*>          = 'A': all singular values will be found.
*>          = 'V': all singular values in the half-open interval [VL,VU)
*>                 will be found.
*>          = 'I': the IL-th through IU-th singular values will be found.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the bidiagonal matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N)
*>          The n diagonal elements of the bidiagonal matrix B.
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (max(1,N-1))
*>          The (n-1) off-diagonal elements of the bidiagonal matrix
*>          B in elements 1 to N-1.
*> \endverbatim
*>
*> \param[in] VL
*> \verbatim
*>         VL is DOUBLE PRECISION
*>          If RANGE='V', the lower bound of the interval to
*>          be searched for singular values. VU > VL.
*>          Not referenced if RANGE = 'A' or 'I'.
*> \endverbatim
*>
*> \param[in] VU
*> \verbatim
*>         VU is DOUBLE PRECISION
*>          If RANGE='V', the upper bound of the interval to
*>          be searched for singular values. VU > VL.
*>          Not referenced if RANGE = 'A' or 'I'.
*> \endverbatim
*>
*> \param[in] IL
*> \verbatim
*>          IL is INTEGER
*>          If RANGE='I', the index of the smallest singular value
*>          to be returned.  1 <= IL <= IU <= N, if N > 0.
*>          Not referenced if RANGE = 'A' or 'V'.
*> \endverbatim
*>
*> \param[in] IU
*> \verbatim
*>          IU is INTEGER
*>          If RANGE='I', the index of the largest singular value
*>          to be returned.  1 <= IL <= IU <= N, if N > 0.
*>          Not referenced if RANGE = 'A' or 'V'.
*> \endverbatim
*>
*> \param[out] NS
*> \verbatim
*>          NS is INTEGER
*>          The total number of singular values found.  0 <= NS <= N.
*> \endverbatim
*>
*> \param[out] S
*> \verbatim
*>          S is DOUBLE PRECISION array, dimension (N)
*>          The first NS elements contain the selected singular values
*>          in decreasing order (matching DBDSVDX's ordering).
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, K)
*>          If JOBZ = 'V', the first NS columns of Z contain the
*>          singular vectors of B corresponding to the selected
*>          singular values, with U stored in rows 1..N and V in
*>          rows N+1..2N (Z = [U; V]).  If JOBZ = 'N', Z is not
*>          referenced.  K >= NS+1 (the caller must supply an upper
*>          bound on NS when RANGE = 'V').
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          Leading dimension of Z.  LDZ >= 1; if JOBZ = 'V', LDZ >= 2N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          LWORK >= max( 1, 2*N*N + 34*N ) (matches DBDSVDMR3).  If
*>          LWORK = -1, a workspace query is performed: the optimal
*>          size is returned in WORK(1) and no work is done.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (LIWORK)
*> \endverbatim
*>
*> \param[in] LIWORK
*> \verbatim
*>          LIWORK is INTEGER
*>          LIWORK >= max( 1, 20*N ).  If LIWORK = -1, a workspace
*>          query is performed.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  DBDSVDMR3 reported an internal error
*>                (INFO = 1: DLASQ1 failure; INFO = 2: DLARRV_TGK
*>                failure).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup bdsvr
*
*  =====================================================================
      SUBROUTINE DBDSVR( UPLO, JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
     $                   NS, S, Z, LDZ, WORK, LWORK, IWORK, LIWORK,
     $                   INFO )
      IMPLICIT NONE
*
*  -- LAPACK-style driver (offline PR) --
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, N, NS
      DOUBLE PRECISION   VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), S( * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLSV, INDSV, LOWER, LQUERY, VALSV, WANTZ
      INTEGER            I, IDCOPY, IECOPY, ISFULL, IUFULL, IVTFULL,
     $                   IWRK, J, LWMIN, LIWMIN, MFOUND, NEED,
     $                   IWKENG, IIWENG, ILO, IHI, NCOL
      DOUBLE PRECISION   SVAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DBDSVDMR3, DCOPY, DLASET, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      ALLSV  = LSAME( RANGE, 'A' )
      VALSV  = LSAME( RANGE, 'V' )
      INDSV  = LSAME( RANGE, 'I' )
      WANTZ  = LSAME( JOBZ,  'V' )
      LOWER  = LSAME( UPLO,  'L' )
      LQUERY = ( LWORK.EQ.-1 ) .OR. ( LIWORK.EQ.-1 )
*
*     Test the input parameters (mirror DBDSVDX).
*
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LOWER ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( ALLSV .OR. VALSV .OR. INDSV ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.GT.0 ) THEN
         IF( VALSV ) THEN
            IF( VL.LT.ZERO ) THEN
               INFO = -7
            ELSE IF( VU.LE.VL ) THEN
               INFO = -8
            END IF
         ELSE IF( INDSV ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) THEN
               INFO = -9
            ELSE IF( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) THEN
               INFO = -10
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.2*N ) ) INFO = -14
      END IF
*
*     Workspace requirements.  The wrapper needs (in doubles):
*       * N   for a private copy of D
*       * N   for a private copy of E
*       * N   for the full unfiltered spectrum
*       * LDU*N = N*N for a private copy of U (when WANTZ)
*       * LDVT*N = N*N for a private copy of VT (when WANTZ)
*       * 2*N*N + 34*N for the DBDSVDMR3 engine itself
*     Total (WANTZ): 4*N*N + 37*N.  For JOBZ='N' we still budget the
*     larger of DBDSVDMR3's requirement and this to keep sizing
*     branch-free; the caller supplied buffer is not modified beyond
*     the used slice.
*
      IF( N.EQ.0 ) THEN
         LWMIN  = 1
         LIWMIN = 1
      ELSE IF( WANTZ ) THEN
         LWMIN  = MAX( 1, 4*N*N + 37*N )
         LIWMIN = MAX( 1, 20*N )
      ELSE
         LWMIN  = MAX( 1, 2*N*N + 37*N )
         LIWMIN = MAX( 1, 20*N )
      END IF
*
      IF( INFO.EQ.0 .AND. .NOT.LQUERY ) THEN
         IF( LWORK.LT.LWMIN ) THEN
            INFO = -16
         ELSE IF( LIWORK.LT.LIWMIN ) THEN
            INFO = -18
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DBDSVR', -INFO )
         RETURN
      END IF
*
      IF( LQUERY ) THEN
         WORK( 1 ) = DBLE( LWMIN )
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF
*
*     Quick return if possible.
*
      NS = 0
      IF( N.EQ.0 ) RETURN
*
      IF( N.EQ.1 ) THEN
         SVAL = ABS( D( 1 ) )
         IF( ALLSV .OR. INDSV ) THEN
            NS = 1
            S( 1 ) = SVAL
         ELSE
            IF( VL.LT.SVAL .AND. VU.GE.SVAL ) THEN
               NS = 1
               S( 1 ) = SVAL
            END IF
         END IF
         IF( WANTZ .AND. NS.GT.0 ) THEN
            Z( 1, 1 ) = SIGN( ONE, D( 1 ) )
            Z( 2, 1 ) = ONE
         END IF
         RETURN
      END IF
*
*     Slice the WORK buffer.
*
      IDCOPY  = 1
      IECOPY  = IDCOPY + N
      ISFULL  = IECOPY + N
      IF( WANTZ ) THEN
         IUFULL  = ISFULL + N
         IVTFULL = IUFULL + N*N
         IWRK    = IVTFULL + N*N
      ELSE
         IUFULL  = ISFULL + N
         IVTFULL = IUFULL
         IWRK    = ISFULL + N
      END IF
*
*     Stage inputs (DBDSVDMR3 overwrites D and E; we must not
*     modify the caller's arrays).
*
      CALL DCOPY( N, D, 1, WORK( IDCOPY ), 1 )
      IF( N.GT.1 ) THEN
         CALL DCOPY( N-1, E, 1, WORK( IECOPY ), 1 )
      END IF
      WORK( IECOPY + N - 1 ) = ZERO
*
*     Compute the full spectrum via the MR^3 engine.
*
      IWKENG = LWORK - IWRK + 1
      IIWENG = LIWORK
      IF( WANTZ ) THEN
         CALL DBDSVDMR3( 'V', UPLO, N,
     $                   WORK( IDCOPY ), WORK( IECOPY ),
     $                   WORK( ISFULL ),
     $                   WORK( IUFULL  ), N,
     $                   WORK( IVTFULL ), N,
     $                   MFOUND,
     $                   WORK( IWRK ), IWKENG,
     $                   IWORK, IIWENG, INFO )
      ELSE
         CALL DBDSVDMR3( 'N', UPLO, N,
     $                   WORK( IDCOPY ), WORK( IECOPY ),
     $                   WORK( ISFULL ),
     $                   WORK( IUFULL  ), MAX( 1, N ),
     $                   WORK( IVTFULL ), MAX( 1, N ),
     $                   MFOUND,
     $                   WORK( IWRK ), IWKENG,
     $                   IWORK, IIWENG, INFO )
      END IF
      IF( INFO.NE.0 ) RETURN
*
*     Filter and pack the requested slice into S and Z.
*
      IF( ALLSV ) THEN
         ILO = 1
         IHI = N
      ELSE IF( INDSV ) THEN
         ILO = IL
         IHI = IU
      ELSE
*        RANGE='V'; scan ascending spectrum for [VL, VU).
         ILO = N + 1
         IHI = 0
         DO 10 I = 1, N
            SVAL = WORK( ISFULL + I - 1 )
            IF( SVAL.GE.VL .AND. SVAL.LT.VU ) THEN
               IF( I.LT.ILO ) ILO = I
               IF( I.GT.IHI ) IHI = I
            END IF
   10    CONTINUE
      END IF
*
      IF( IHI.GE.ILO ) THEN
         NS = IHI - ILO + 1
*        Reverse ascending -> descending so S(1) is the largest,
*        matching DBDSVDX's convention (see /tmp/lapack-ref/SRC/dbdsvdx.f
*        line 747 sort loop).
         DO 20 I = 1, NS
            S( I ) = WORK( ISFULL + IHI - I )
   20    CONTINUE
         IF( WANTZ ) THEN
*           Column J of Z holds the vector for S(J).  In the engine's
*           spectrum that is column (IHI - J + 1).
*           Rows 1..N of Z(:, J) hold U(:, IHI-J+1).
*           Rows N+1..2N of Z(:, J) hold VT(IHI-J+1, :)^T.
            NCOL = NS
            DO 40 J = 1, NCOL
               DO 30 I = 1, N
                  Z( I,     J ) = WORK( IUFULL  + (IHI-J)*N + I - 1 )
                  Z( N + I, J ) = WORK( IVTFULL + (I-1)*N + (IHI-J) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE
         NS = 0
      END IF
*
      RETURN
*
*     End of DBDSVR
*
      END
