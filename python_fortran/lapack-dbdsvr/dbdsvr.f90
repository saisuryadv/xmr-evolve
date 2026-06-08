!> \brief \b DBDSVR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!     SUBROUTINE DBDSVR( UPLO, JOBZ, RANGE, N, D, E, VL, VU, IL, IU, &
!                        NS, S, Z, LDZ, WORK, LWORK, IWORK, INFO )
!
!     .. Scalar Arguments ..
!     CHARACTER          JOBZ, RANGE, UPLO
!     INTEGER            IL, INFO, IU, LDZ, LWORK, N, NS
!     DOUBLE PRECISION   VL, VU
!     ..
!     .. Array Arguments ..
!     INTEGER            IWORK( * )
!     DOUBLE PRECISION   D( * ), E( * ), S( * ), WORK( * ), Z( LDZ, * )
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DBDSVR computes the singular value decomposition (SVD) of a real
!> N-by-N (upper or lower) bidiagonal matrix B, B = U * S * VT,
!> where S is a diagonal matrix with non-negative diagonal elements
!> (the singular values of B), and U and VT are orthogonal matrices
!> of left and right singular vectors, respectively.
!>
!> DBDSVR mirrors the API of LAPACK's DBDSVDX but replaces the
!> Golub-Kahan tridiagonal eigensolver with two specialised kernels:
!>
!>   1.  Singular values: DLASQ1 (dqds) — high relative accuracy on
!>       every singular value, regardless of cluster structure.
!>
!>   2.  Singular vectors: Willems MR^3-GK (DMR3GK_SVD) — a multiple
!>       relative robust representation algorithm specialised to the
!>       Golub-Kahan tridiagonal of the bidiagonal matrix.  See
!>       P. Willems and B. Lang, "A framework for the MR^3 algorithm:
!>       theory and implementation", SIAM J. Sci. Comput., 35:740-766,
!>       2013.
!>
!> Sigmas from MR^3-GK are paired positionally with DLASQ1's sorted
!> singular values: sort the MR^3-GK output descending, then pair the
!> k-th MR^3-GK column with S(k).  This is robust to clustered
!> singular values, where greedy nearest-match would mis-pair vectors.
!>
!> Only RANGE='A' and JOBZ='V' are supported.  Other inputs are
!> rejected via XERBLA, keeping VL, VU, IL, IU in the signature for
!> drop-in compatibility with DBDSVDX (they are not referenced).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  B is upper bidiagonal;
!>          = 'L':  B is lower bidiagonal.
!> \endverbatim
!>
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER*1
!>          = 'V':  Compute singular values and singular vectors.
!>          DBDSVR rejects JOBZ='N' with XERBLA.
!> \endverbatim
!>
!> \param[in] RANGE
!> \verbatim
!>          RANGE is CHARACTER*1
!>          = 'A':  all singular values will be found.
!>          DBDSVR rejects RANGE='V' and 'I' with XERBLA.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the bidiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The N diagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (max(1,N-1))
!>          The (N-1) off-diagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] VL, VU, IL, IU
!> \verbatim
!>          Unused.  Reserved for DBDSVDX compatibility.
!> \endverbatim
!>
!> \param[out] NS
!> \verbatim
!>          NS is INTEGER.  On exit NS = N (RANGE='A').
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (N).
!>          The N singular values in DECREASING order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ, N).
!>          On exit, if INFO = 0:
!>            Z(1:N,   k) = k-th left  singular vector u_k.
!>            Z(N+1:2N,k) = k-th right singular vector v_k.
!>          For UPLO='L', U and V are swapped (mirroring DBDSVDX).
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER.  LDZ >= max(1, 2*N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)).
!>          On exit, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER.
!>          The dimension of WORK.  LWORK >= 2*N*N + 7*N.
!>          If LWORK = -1, a workspace query is assumed; the routine
!>          returns the optimal LWORK in WORK(1) and the optimal LIWORK
!>          in IWORK(1), and does no other work.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N).
!>          Used internally for the sigma-vector positional pairing.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER.
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = N+1, a singular value reported by DLASQ1
!>                  and the corresponding sigma from MR^3-GK differ by
!>                  more than 16*N*EPS*||B||.  Output is still produced
!>                  but should be inspected for numerical mismatch.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Sai Surya Duvvuri
!> \author Based on LAPACK DBDSVDX template (Univ. of Tennessee,
!>         Univ. of California Berkeley, Univ. of Colorado Denver,
!>         NAG Ltd.)
!> \date June 2026
!
!> \ingroup bdsvdx
!  =====================================================================
subroutine dbdsvr(uplo, jobz, range, n, d, e, vl, vu, il, iu, &
                  ns, s, z, ldz, work, lwork, iwork, info)
   use mr3gk_consts, only : dp
   use mr3gk,        only : dmr3gk_svd
   implicit none

!  -- LAPACK driver routine --

!     .. Scalar Arguments ..
   character          :: jobz, range, uplo
   integer            :: il, info, iu, ldz, lwork, n, ns
   double precision   :: vl, vu
!     .. Array Arguments ..
   integer            :: iwork(*)
   double precision   :: d(*), e(*), s(*), work(*), z(ldz, *)
!     ..
!     .. Parameters ..
   double precision, parameter :: zero = 0.0d0, one = 1.0d0
!     ..
!     .. Local Scalars ..
   logical            :: allsv, indsv, lower, valsv, wantz, lquery
   integer            :: i, j, k, jbest, info_loc, lwkopt
   double precision   :: eps, smax, dtemp, dbest, mismatch_tol
!     ..
!     .. Local Allocatables (F90 stylistic — pre-LAPACK-PR port to F77
!        will replace these with slices of WORK; thread-safe in either
!        form because allocatables are stack/heap-local per call) ..
   double precision, allocatable :: dcopy_loc(:), ecopy_loc(:)
   double precision, allocatable :: wlasq(:), sigma_mr3(:)
   double precision, allocatable :: u_tmp(:,:), v_tmp(:,:)
!     ..
!     .. External Functions ..
   logical, external            :: lsame
   integer, external            :: idamax
   double precision, external   :: dlamch
!     ..
!     .. External Subroutines ..
   external           :: dcopy, dlasq1, dswap, xerbla
!     ..
!     .. Intrinsic Functions ..
   intrinsic abs, dble, max, sign
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
   allsv  = lsame( range, 'A' )
   valsv  = lsame( range, 'V' )
   indsv  = lsame( range, 'I' )
   wantz  = lsame( jobz,  'V' )
   lower  = lsame( uplo,  'L' )
   lquery = ( lwork.eq.-1 )

   info   = 0
   lwkopt = max( 1, 2*n*n + 7*n )

   if( .not.lsame( uplo, 'U' ) .and. .not.lower ) then
      info = -1
   else if( .not.wantz ) then
      info = -2
   else if( .not.allsv ) then
      info = -3
   else if( n.lt.0 ) then
      info = -4
   else if( ldz.lt.max( 1, 2*n ) ) then
      info = -14
   else if( lwork.lt.lwkopt .and. .not.lquery ) then
      info = -16
   end if

   if( info.ne.0 ) then
      call xerbla( 'DBDSVR', -info )
      return
   end if

!  Workspace query: return required sizes and exit.
   if( lquery ) then
      work(1)  = dble( lwkopt )
      iwork(1) = n
      return
   end if
!
!  Quick return if possible.
!
   ns = 0
   if( n.eq.0 ) return

   if( n.eq.1 ) then
      ns     = 1
      s( 1 ) = abs( d( 1 ) )
      z( 1, 1 ) = sign( one, d( 1 ) )
      z( 2, 1 ) = one
      if( lower ) then
         dtemp = z( 1, 1 )
         z( 1, 1 ) = z( 2, 1 )
         z( 2, 1 ) = dtemp
      end if
      work(1) = dble( lwkopt )
      return
   end if
!
!  Approximate ||B|| for the mismatch tolerance.
!
   i    = idamax( n,   d, 1 ); smax = abs( d( i ) )
   i    = idamax( n-1, e, 1 ); smax = max( smax, abs( e( i ) ) )
   eps  = dlamch( 'Epsilon' )
!
!  Phase 1: high-accuracy singular values via DLASQ1.
!  DLASQ1 overwrites its D argument with sigmas in DECREASING order.
!
   allocate( dcopy_loc(n), ecopy_loc(max(n-1,1)), wlasq(4*n) )
   call dcopy( n,   d, 1, dcopy_loc, 1 )
   call dcopy( n-1, e, 1, ecopy_loc, 1 )
   call dlasq1( n, dcopy_loc, ecopy_loc, wlasq, info_loc )
   if( info_loc.ne.0 ) then
      info = info_loc
      deallocate( dcopy_loc, ecopy_loc, wlasq )
      return
   end if
   ns = n
   call dcopy( n, dcopy_loc, 1, s, 1 )
   deallocate( dcopy_loc, ecopy_loc, wlasq )
!
!  Phase 2: singular vectors via Willems MR^3-GK pipeline.
!  DMR3GK_SVD produces N sigmas (in block order, NOT sorted) plus
!  N-by-N column-major U_tmp, V_tmp.
!
   allocate( sigma_mr3(n), u_tmp(n,n), v_tmp(n,n) )
   call dmr3gk_svd( n, d, e, sigma_mr3, u_tmp, v_tmp, n, info_loc )
   if( info_loc.ne.0 ) then
      info = info_loc
      deallocate( sigma_mr3, u_tmp, v_tmp )
      return
   end if
!
!  Phase 3: positional pairing — sort MR^3-GK output by its own
!  sigma (decreasing) and pair with DLASQ1's sigmas (also decreasing).
!  Positional pairing keeps clustered triples intact, where greedy
!  nearest-match can put the wrong vector pair under a sigma.
!
   do j = 1, n
      iwork( j ) = 0
   end do

   mismatch_tol = 16.0d0 * dble( n ) * eps * smax

   do k = 1, n
      jbest = 0
      dbest = zero
      do j = 1, n
         if( iwork( j ).eq.0 .and. &
             ( jbest.eq.0 .or. sigma_mr3( j ).gt.dbest ) ) then
            jbest = j
            dbest = sigma_mr3( j )
         end if
      end do
      iwork( jbest ) = 1
      dtemp = abs( dbest - s( k ) )
      if( dtemp.gt.mismatch_tol .and. info.eq.0 ) then
         info = n + 1
      end if
!
!     Copy MR^3-GK column JBEST into Z column K; UPLO='L' swaps U and V.
!
      if( lower ) then
         call dcopy( n, v_tmp( 1, jbest ), 1, z( 1,   k ), 1 )
         call dcopy( n, u_tmp( 1, jbest ), 1, z( n+1, k ), 1 )
      else
         call dcopy( n, u_tmp( 1, jbest ), 1, z( 1,   k ), 1 )
         call dcopy( n, v_tmp( 1, jbest ), 1, z( n+1, k ), 1 )
      end if
   end do

   deallocate( sigma_mr3, u_tmp, v_tmp )

   work(1) = dble( lwkopt )

!  Suppress unused-argument warnings for VL, VU, IL, IU
!  (kept in signature for DBDSVDX compatibility).
   if( valsv .or. indsv ) then
      dtemp = vl + vu + dble( il ) + dble( iu )
   end if

   return
end subroutine dbdsvr
