! ===========================================================
! Module: mr3gk_utils
! Wrappers for LAPACK DSTEBZ (bisection eigenvalue) used by Python's
! bisect_evals and refine_eval, plus simple helpers.
! ===========================================================
module mr3gk_utils
   use mr3gk_consts
   implicit none
contains

   ! Compute eigenvalues of symmetric tridiagonal (c, e_off) of size n,
   ! returning indices il..iu (1-based, inclusive). Output stored in w(1..nfound).
   ! Equivalent to scipy.linalg.lapack._flapack.dstebz with range='I'.
   subroutine bisect_evals(c, e_off, n, il, iu, tol, w, nfound, info)
      integer, intent(in) :: n, il, iu
      real(dp), intent(in) :: c(*), e_off(*), tol
      real(dp), intent(out) :: w(*)
      integer, intent(out) :: nfound, info

      real(dp) :: vl_d, vu_d
      integer :: m_o, nsplit
      integer, allocatable :: iblock(:), isplit(:), iwork(:)
      real(dp), allocatable :: work(:)
      external :: dstebz

      nfound = 0
      info = 0
      if (n <= 0 .or. il > iu) return
      if (il < 1 .or. iu > n) then
         info = -1
         return
      end if

      vl_d = 0.0_dp
      vu_d = 0.0_dp
      allocate(iblock(n), isplit(n), iwork(3*n), work(4*n))
      call dstebz('I', 'E', n, vl_d, vu_d, il, iu, tol, &
                  c, e_off, m_o, nsplit, w, iblock, isplit, &
                  work, iwork, info)
      if (info == 0) then
         nfound = m_o
      end if
      deallocate(iblock, isplit, iwork, work)
   end subroutine bisect_evals

   ! Refine a single eigenvalue near `approx` using DSTEBZ with RANGE='V'.
   ! Mirrors refine_eval in mr3_gk.py (returns approx if no eigenvalue found).
   real(dp) function refine_eval(c, e_off, n, approx, tol)
      integer, intent(in) :: n
      real(dp), intent(in) :: c(*), e_off(*), approx, tol

      real(dp) :: gap, vl_d, vu_d
      integer  :: m_o, nsplit, info
      integer, allocatable :: iblock(:), isplit(:), iwork(:)
      real(dp), allocatable :: work(:), w_loc(:)
      integer :: jbest, j
      real(dp) :: best_diff, diff
      external :: dstebz

      refine_eval = approx
      gap = max(abs(approx)*1.0d-3, 10.0d0*EPS*abs(approx), SAFMIN)
      vl_d = approx - gap
      vu_d = approx + gap
      allocate(iblock(n), isplit(n), iwork(3*n), work(4*n), w_loc(n))
      call dstebz('V', 'E', n, vl_d, vu_d, 0, 0, tol, &
                  c, e_off, m_o, nsplit, w_loc, iblock, isplit, &
                  work, iwork, info)
      if (info == 0 .and. m_o > 0) then
         jbest = 1
         best_diff = abs(w_loc(1) - approx)
         do j = 2, m_o
            diff = abs(w_loc(j) - approx)
            if (diff < best_diff) then
               best_diff = diff
               jbest = j
            end if
         end do
         refine_eval = w_loc(jbest)
      end if
      deallocate(iblock, isplit, iwork, work, w_loc)
   end function refine_eval

   ! ||x||_2 via DNRM2
   real(dp) function vec_norm(n, x)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(*)
      real(dp), external :: dnrm2
      vec_norm = dnrm2(n, x, 1)
   end function vec_norm

end module mr3gk_utils
