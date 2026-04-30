! ===========================================================
! Module: mr3gk_qrsweep
! Zero-shift QR sweep (Demmel-Kahan 1990) and Givens application,
! mirroring zero_shift_qr_sweep and _apply_givens_to_rows in mr3_gk.py.
! ===========================================================
module mr3gk_qrsweep
   use mr3gk_consts
   implicit none
contains

   ! Generates a Givens plane rotation. Calls LAPACK DLARTG.
   ! [cs sn; -sn cs] * [f;g] = [r;0]
   subroutine my_dlartg(f, g, cs, sn, r)
      real(dp), intent(in) :: f, g
      real(dp), intent(out) :: cs, sn, r
      external :: dlartg
      call dlartg(f, g, cs, sn, r)
   end subroutine my_dlartg

   ! Output: rotr_cs(1:n-1), rotr_sn(1:n-1), rotr_idx(1:n-1) - right rotations
   !         rotl_cs(1:n-1), rotl_sn(1:n-1), rotl_idx(1:n-1) - left rotations
   ! In Python rotations are stored 0-based in `i`. Here we use i in [0, n-2].
   ! Modifies d(0:n-1), e(0:n-2) in place.
   subroutine dmr3gk_qrsweep(n, d, e, rotr_cs, rotr_sn, rotr_idx, &
                                       rotl_cs, rotl_sn, rotl_idx, nrots)
      integer, intent(in) :: n
      real(dp), intent(inout) :: d(0:n-1), e(0:max(n-2,0))
      real(dp), intent(out) :: rotr_cs(:), rotr_sn(:), rotl_cs(:), rotl_sn(:)
      integer, intent(out) :: rotr_idx(:), rotl_idx(:), nrots

      integer :: i
      real(dp) :: cs, sn, r, oldcs, oldsn, d1, d2, h

      nrots = 0
      if (n <= 1) return

      cs = 1.0_dp
      oldcs = 1.0_dp
      oldsn = 0.0_dp
      do i = 0, n - 2
         d1 = d(i) * cs
         call my_dlartg(d1, e(i), cs, sn, r)
         nrots = nrots + 1
         rotr_cs(nrots) = cs
         rotr_sn(nrots) = sn
         rotr_idx(nrots) = i
         if (i > 0) then
            e(i-1) = oldsn * r
         end if
         d1 = oldcs * r
         d2 = d(i+1) * sn
         call my_dlartg(d1, d2, oldcs, oldsn, d(i))
         rotl_cs(nrots) = oldcs
         rotl_sn(nrots) = oldsn
         rotl_idx(nrots) = i
      end do
      h = d(n-1) * cs
      d(n-1) = h * oldcs
      e(n-2) = h * oldsn
   end subroutine dmr3gk_qrsweep

   ! Apply G_0 @ G_1 @ ... @ G_{m-1} to M from the left (REVERSED order).
   ! G_i acts on rows (row_offset+idx, row_offset+idx+1):
   !   M[ri,:]   = cs*r0 - sn*r1
   !   M[ri+1,:] = sn*r0 + cs*r1
   ! Mirrors _apply_givens_to_rows in mr3_gk.py.
   subroutine apply_givens_rows(M, ldm, ncols, rcs, rsn, ridx, nrots, row_offset)
      integer, intent(in) :: ldm, ncols, nrots, row_offset
      real(dp), intent(inout) :: M(ldm, ncols)
      real(dp), intent(in) :: rcs(:), rsn(:)
      integer, intent(in) :: ridx(:)
      integer :: t, ri, j
      real(dp) :: cs, sn, r0, r1
      do t = nrots, 1, -1
         cs = rcs(t)
         sn = rsn(t)
         ri = row_offset + ridx(t) + 1   ! +1 because Fortran arrays are 1-based
         do j = 1, ncols
            r0 = M(ri, j)
            r1 = M(ri+1, j)
            M(ri, j)   = cs*r0 - sn*r1
            M(ri+1, j) = sn*r0 + cs*r1
         end do
      end do
   end subroutine apply_givens_rows

end module mr3gk_qrsweep
