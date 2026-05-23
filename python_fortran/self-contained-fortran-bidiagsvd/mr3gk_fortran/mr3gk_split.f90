! ===========================================================
! Module: mr3gk_split
! Two-phase splitting of bidiagonal matrix (mirrors split_bidiag in mr3_gk.py).
! ===========================================================
module mr3gk_split
   use mr3gk_consts
   implicit none
contains

   ! Two-phase split of bidiagonal (d, e) of size n.
   ! Output: nb = number of blocks, b_beg(1:nb), b_end(1:nb) (0-based to match Python).
   subroutine dmr3gk_split(d, e, n, nb, b_beg, b_end)
      integer, intent(in) :: n
      real(dp), intent(in) :: d(0:n-1), e(0:max(n-2,0))
      integer, intent(out) :: nb
      integer, intent(out) :: b_beg(:), b_end(:)

      integer :: i, j, beg, rb_beg, rb_end, sub_beg, rb_k
      integer :: nrel
      integer, allocatable :: rb_b(:), rb_e(:)
      real(dp) :: rel_thresh, abs_thresh, bnorm, s

      nb = 0
      if (n <= 1) then
         nb = 1
         b_beg(1) = 0
         b_end(1) = n - 1
         return
      end if

      allocate(rb_b(n), rb_e(n))
      nrel = 0
      beg = 0
      do i = 0, n - 2
         rel_thresh = EPS * (abs(d(i)) + abs(d(i+1)))
         if (abs(e(i)) <= max(rel_thresh, SAFMIN)) then
            nrel = nrel + 1
            rb_b(nrel) = beg
            rb_e(nrel) = i
            beg = i + 1
         end if
      end do
      nrel = nrel + 1
      rb_b(nrel) = beg
      rb_e(nrel) = n - 1

      do i = 1, nrel
         rb_beg = rb_b(i)
         rb_end = rb_e(i)
         rb_k = rb_end - rb_beg + 1
         if (rb_k <= 1) then
            nb = nb + 1
            b_beg(nb) = rb_beg
            b_end(nb) = rb_end
            cycle
         end if
         ! ||B_sub||_inf
         bnorm = 0.0_dp
         do j = rb_beg, rb_end
            s = abs(d(j))
            if (j < rb_end) s = s + abs(e(j))
            if (s > bnorm) bnorm = s
         end do
         abs_thresh = real(rb_k, dp) * EPS * bnorm

         sub_beg = rb_beg
         do j = rb_beg, rb_end - 1
            if (abs(e(j)) <= max(abs_thresh, SAFMIN)) then
               nb = nb + 1
               b_beg(nb) = sub_beg
               b_end(nb) = j
               sub_beg = j + 1
            end if
         end do
         nb = nb + 1
         b_beg(nb) = sub_beg
         b_end(nb) = rb_end
      end do

      deallocate(rb_b, rb_e)
   end subroutine dmr3gk_split

end module mr3gk_split
