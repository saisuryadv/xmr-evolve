! ===========================================================
! Module: mr3gk
! Public entry point: dmr3gk_svd. Mirrors bidiag_svd in mr3_gk.py.
! ===========================================================
module mr3gk
   use mr3gk_consts
   use mr3gk_split
   use mr3gk_qrsweep
   use mr3gk_tgk
   use mr3gk_postproc
   use mr3gk_utils, only: vec_norm
   implicit none

contains

   ! Bidiagonal SVD via MR³ on Golub-Kahan tridiagonal.
   !   d, e   : input bidiagonal (length n, n-1). NOT modified.
   !   sigma  : (out) singular values (length n)
   !   U, V   : (out) left/right singular vectors, n×n stored column-major
   !   info_out : 0 on success
   recursive subroutine dmr3gk_svd(n, d_in, e_in, sigma, U, V, ldUV, info_out)
      integer, intent(in) :: n, ldUV
      real(dp), intent(in) :: d_in(0:n-1), e_in(0:max(n-2,0))
      real(dp), intent(out) :: sigma(*), U(ldUV, *), V(ldUV, *)
      integer, intent(out) :: info_out

      real(dp), allocatable :: d(:), e(:)
      real(dp), allocatable :: d1_arr(:), d2_arr(:)
      real(dp) :: bnorm_scale, scale_factor, s_e, s_d
      integer  :: i, j, ii

      ! Splitting outputs
      integer, allocatable :: b_beg(:), b_end(:)
      integer :: nb

      ! Sub-block storage
      integer, allocatable :: sing_rows(:), sing_cols(:)
      integer :: n_sing
      integer, allocatable :: mblk_bbeg(:), mblk_bend(:), mblk_bk(:), mblk_col(:)
      integer :: n_mblk
      integer, allocatable :: qr_def_cols(:)
      integer :: n_qr_def

      ! QR sweep buffers
      real(dp), allocatable :: d_sw(:), e_sw(:)
      real(dp), allocatable :: rotr_cs(:), rotr_sn(:), rotl_cs(:), rotl_sn(:)
      integer,  allocatable :: rotr_idx(:), rotl_idx(:)
      integer :: nrots
      logical :: clean_split
      real(dp) :: zero_sv_thresh

      ! For zero-sv recursive call
      real(dp), allocatable :: d_sub(:), e_sub(:), sigma_sub(:)
      real(dp), allocatable :: U_sub(:,:), V_sub(:,:)
      real(dp), allocatable :: U_blk(:,:), V_blk(:,:)
      integer :: bk, bk_sub, info_sub

      ! T_GK solver buffers
      real(dp), allocatable :: w_blk(:), Z_blk(:,:)

      ! Multi-col tracking
      integer, allocatable :: mc_list(:)
      integer :: n_mc

      integer :: bbeg, bend, col, c0, c1, ncols

      info_out = 0
      if (n <= 0) return
      if (n == 1) then
         sigma(1) = abs(d_in(0))
         if (d_in(0) >= 0.0_dp) then
            U(1, 1) = 1.0_dp
         else
            U(1, 1) = -1.0_dp
         end if
         V(1, 1) = 1.0_dp
         return
      end if

      ! Working copies
      allocate(d(0:n-1), e(0:max(n-2,0)))
      do i = 0, n - 1
         d(i) = d_in(i)
      end do
      do i = 0, n - 2
         e(i) = e_in(i)
      end do

      ! Pre-scaling
      bnorm_scale = 0.0_dp
      do i = 0, n - 1
         if (abs(d(i)) > bnorm_scale) bnorm_scale = abs(d(i))
      end do
      do i = 0, n - 2
         if (abs(e(i)) > bnorm_scale) bnorm_scale = abs(e(i))
      end do
      if (bnorm_scale > 0.0_dp) then
         scale_factor = 1.0_dp / bnorm_scale
         do i = 0, n - 1
            d(i) = d(i) * scale_factor
         end do
         do i = 0, n - 2
            e(i) = e(i) * scale_factor
         end do
      else
         scale_factor = 1.0_dp
      end if

      ! Split bidiagonal
      allocate(b_beg(n+1), b_end(n+1))
      call dmr3gk_split(d, e, n, nb, b_beg, b_end)
      ! Zero out e at split boundaries (mirroring Python)
      do i = 1, nb
         if (b_end(i) < n - 1) then
            e(b_end(i)) = 0.0_dp
         end if
      end do

      ! Sign matrices D1, D2
      allocate(d1_arr(0:n-1), d2_arr(0:n-1))
      do i = 0, n - 1
         d1_arr(i) = 1.0_dp
         d2_arr(i) = 1.0_dp
      end do
      d1_arr(0) = 1.0_dp
      if (d(0) > 0.0_dp) then
         d2_arr(0) = 1.0_dp
      else if (d(0) < 0.0_dp) then
         d2_arr(0) = -1.0_dp
      else
         d2_arr(0) = 1.0_dp
      end if
      do i = 0, n - 2
         if (e(i) > 0.0_dp) then
            s_e = 1.0_dp
         else if (e(i) < 0.0_dp) then
            s_e = -1.0_dp
         else
            s_e = 1.0_dp
         end if
         if (abs(d1_arr(i)) > 0.0_dp) then
            d2_arr(i+1) = s_e / d1_arr(i)
         else
            d2_arr(i+1) = s_e
         end if
         if (d(i+1) > 0.0_dp) then
            s_d = 1.0_dp
         else if (d(i+1) < 0.0_dp) then
            s_d = -1.0_dp
         else
            s_d = 1.0_dp
         end if
         if (abs(d2_arr(i+1)) > 0.0_dp) then
            d1_arr(i+1) = s_d / d2_arr(i+1)
         else
            d1_arr(i+1) = s_d
         end if
      end do

      ! Initialize outputs
      do j = 1, n
         sigma(j) = 0.0_dp
         do i = 1, n
            U(i, j) = 0.0_dp
            V(i, j) = 0.0_dp
         end do
      end do

      ! Block dispatch with zero-shift QR deflation
      allocate(sing_rows(n), sing_cols(n))
      allocate(mblk_bbeg(n), mblk_bend(n), mblk_bk(n), mblk_col(n))
      allocate(qr_def_cols(n))
      n_sing = 0
      n_mblk = 0
      n_qr_def = 0
      zero_sv_thresh = real(n, dp) * EPS
      col = 0
      do i = 1, nb
         bbeg = b_beg(i)
         bend = b_end(i)
         bk = bend - bbeg + 1
         if (bk <= 0) cycle
         if (bk == 1) then
            n_sing = n_sing + 1
            sing_rows(n_sing) = bbeg
            sing_cols(n_sing) = col
            col = col + 1
         else
            ! One zero-shift QR sweep on block
            allocate(d_sw(0:bk-1), e_sw(0:max(bk-2,0)))
            do j = 0, bk - 1
               d_sw(j) = d(bbeg + j)
            end do
            do j = 0, bk - 2
               e_sw(j) = e(bbeg + j)
            end do
            allocate(rotr_cs(bk), rotr_sn(bk), rotr_idx(bk))
            allocate(rotl_cs(bk), rotl_sn(bk), rotl_idx(bk))
            call dmr3gk_qrsweep(bk, d_sw, e_sw, rotr_cs, rotr_sn, rotr_idx, &
                                rotl_cs, rotl_sn, rotl_idx, nrots)
            clean_split = (abs(d_sw(bk-1)) < zero_sv_thresh) .and. &
                          (abs(e_sw(bk-2)) < zero_sv_thresh)
            if (clean_split) then
               ! Solve (bk-1) sub-problem (unscale for recursive call)
               bk_sub = bk - 1
               allocate(d_sub(0:bk_sub-1), e_sub(0:max(bk_sub-2,0)))
               allocate(sigma_sub(bk_sub))
               allocate(U_sub(bk_sub, bk_sub), V_sub(bk_sub, bk_sub))
               do j = 0, bk_sub - 1
                  d_sub(j) = d_sw(j) / scale_factor
               end do
               do j = 0, bk_sub - 2
                  e_sub(j) = e_sw(j) / scale_factor
               end do
               call dmr3gk_svd(bk_sub, d_sub, e_sub, sigma_sub, U_sub, V_sub, bk_sub, info_sub)

               ! Embed into bk×bk: top-left = sub, bottom-right = 1
               allocate(U_blk(bk, bk), V_blk(bk, bk))
               do j = 1, bk
                  do ii = 1, bk
                     if (ii == j) then
                        U_blk(ii, j) = 1.0_dp
                        V_blk(ii, j) = 1.0_dp
                     else
                        U_blk(ii, j) = 0.0_dp
                        V_blk(ii, j) = 0.0_dp
                     end if
                  end do
               end do
               do j = 1, bk_sub
                  do ii = 1, bk_sub
                     U_blk(ii, j) = U_sub(ii, j)
                     V_blk(ii, j) = V_sub(ii, j)
                  end do
               end do
               ! Apply QR rotations
               call apply_givens_rows(V_blk, bk, bk, rotr_cs, rotr_sn, rotr_idx, nrots, 0)
               call apply_givens_rows(U_blk, bk, bk, rotl_cs, rotl_sn, rotl_idx, nrots, 0)

               ! Store
               do j = 1, bk_sub
                  sigma(col + j) = sigma_sub(j) * scale_factor
               end do
               sigma(col + bk) = abs(d_sw(bk-1))
               do j = 1, bk
                  do ii = 1, bk
                     U(bbeg + ii, col + j) = U_blk(ii, j)
                     V(bbeg + ii, col + j) = V_blk(ii, j)
                  end do
               end do
               do j = 0, bk - 1
                  n_qr_def = n_qr_def + 1
                  qr_def_cols(n_qr_def) = col + j
               end do
               col = col + bk
               deallocate(d_sub, e_sub, sigma_sub, U_sub, V_sub, U_blk, V_blk)
            else
               n_mblk = n_mblk + 1
               mblk_bbeg(n_mblk) = bbeg
               mblk_bend(n_mblk) = bend
               mblk_bk(n_mblk) = bk
               mblk_col(n_mblk) = col
               col = col + bk
            end if
            deallocate(d_sw, e_sw, rotr_cs, rotr_sn, rotr_idx, rotl_cs, rotl_sn, rotl_idx)
         end if
      end do

      ! Singletons
      do j = 1, n_sing
         sigma(sing_cols(j) + 1) = abs(d(sing_rows(j)))
         V(sing_rows(j) + 1, sing_cols(j) + 1) = d2_arr(sing_rows(j))
         U(sing_rows(j) + 1, sing_cols(j) + 1) = d1_arr(sing_rows(j))
      end do

      ! Multi-blocks
      allocate(mc_list(n))
      n_mc = 0
      do i = 1, n_mblk
         bbeg = mblk_bbeg(i)
         bend = mblk_bend(i)
         bk = mblk_bk(i)
         c0 = mblk_col(i)
         allocate(w_blk(bk), Z_blk(2*bk, bk))
         call mr3_tgk_multiblock(d, e, n, bbeg, bend, w_blk, Z_blk, 2*bk, info_sub)
         ncols = min(bk, n - c0)
         do j = 1, ncols
            sigma(c0 + j) = abs(w_blk(j))
            do c1 = 0, bk - 1
               V(bbeg + 1 + c1, c0 + j) = Z_blk(2*c1 + 1, j) * d2_arr(bbeg + c1)
               U(bbeg + 1 + c1, c0 + j) = Z_blk(2*c1 + 2, j) * d1_arr(bbeg + c1)
            end do
            n_mc = n_mc + 1
            mc_list(n_mc) = c0 + j   ! 1-based
         end do
         deallocate(w_blk, Z_blk)
      end do

      ! Post-processing on multi-block columns. Use ORIGINAL (un-prescaled) (d_in, e_in).
      ! NOTE: Python uses (d, e) which is the SCALED bidiagonal. To match Python exactly,
      ! we use the same scaled (d, e) in post-processing then unscale sigma at the end.
      if (n_mc > 0) then
         call post_process(n, d, e, sigma, U, V, ldUV, mc_list, n_mc)
      end if

      ! Scale sigma back if pre-scaled
      if (scale_factor /= 1.0_dp) then
         do j = 1, n
            sigma(j) = sigma(j) / scale_factor
         end do
      end if

      ! Cleanup
      deallocate(d, e, d1_arr, d2_arr)
      deallocate(b_beg, b_end)
      deallocate(sing_rows, sing_cols)
      deallocate(mblk_bbeg, mblk_bend, mblk_bk, mblk_col, qr_def_cols)
      deallocate(mc_list)
   end subroutine dmr3gk_svd

end module mr3gk
