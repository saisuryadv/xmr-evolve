! ===========================================================
! Module: mr3gk_postproc
! Multi-block post-processing: normalization, Bv recovery, sign fix,
! Gram-Schmidt completion. Mirrors lines 807-936 of mr3_gk.py.
! ===========================================================
module mr3gk_postproc
   use mr3gk_consts
   use mr3gk_utils, only: vec_norm
   implicit none

   external :: dlarnv

contains

   ! Compute B @ M(:, idx_list) where B = diag(d) + superdiag(e), n×k_in.
   ! Inplace into BM(:, 1..k_in).
   subroutine bidiag_matvec_batch(n, d, e, M, ldm, ncols, BM, ldb)
      integer, intent(in) :: n, ldm, ncols, ldb
      real(dp), intent(in) :: d(0:n-1), e(0:max(n-2,0))
      real(dp), intent(in) :: M(ldm, ncols)
      real(dp), intent(out) :: BM(ldb, ncols)
      integer :: i, j
      do j = 1, ncols
         do i = 1, n
            BM(i, j) = d(i-1) * M(i, j)
         end do
      end do
      if (n > 1) then
         do j = 1, ncols
            do i = 0, n - 2
               BM(i+1, j) = BM(i+1, j) + e(i) * M(i+2, j)
            end do
         end do
      end if
   end subroutine bidiag_matvec_batch

   ! Compute B^T @ M
   subroutine bidiagT_matvec_batch(n, d, e, M, ldm, ncols, BM, ldb)
      integer, intent(in) :: n, ldm, ncols, ldb
      real(dp), intent(in) :: d(0:n-1), e(0:max(n-2,0))
      real(dp), intent(in) :: M(ldm, ncols)
      real(dp), intent(out) :: BM(ldb, ncols)
      integer :: i, j
      do j = 1, ncols
         do i = 1, n
            BM(i, j) = d(i-1) * M(i, j)
         end do
      end do
      if (n > 1) then
         do j = 1, ncols
            do i = 0, n - 2
               BM(i+2, j) = BM(i+2, j) + e(i) * M(i+1, j)
            end do
         end do
      end if
   end subroutine bidiagT_matvec_batch

   ! Apply post-processing on multi-block columns.
   !   n     : matrix size
   !   d, e  : original (un-prescaled) bidiagonal
   !   sigma : (length n) singular values
   !   U, V  : (n×n) singular vectors, modified in place
   !   mc(1:n_mc) : 1-based indices of multi-block columns
   subroutine post_process(n, d, e, sigma, U, V, ldUV, mc, n_mc)
      integer, intent(in) :: n, ldUV, n_mc
      real(dp), intent(in) :: d(0:n-1), e(0:max(n-2,0))
      real(dp), intent(inout) :: sigma(*), U(ldUV, n), V(ldUV, n)
      integer, intent(in) :: mc(*)

      real(dp), allocatable :: nv_mc(:), nu_mc(:), max_norms(:), norm_thresh(:)
      logical, allocatable :: v_good(:), u_good(:)
      logical, allocatable :: need_recovery(:), recover_u(:), recover_v(:), both_good(:)
      logical, allocatable :: fix_v(:), fix_u(:)
      real(dp), allocatable :: scale_v_arr(:), scale_u_arr(:)
      real(dp), allocatable :: BV_buf(:,:), BTU_buf(:,:), tmp_n(:), tmp_n2(:)
      real(dp), allocatable :: res_v_arr(:), res_u_arr(:), thresh_arr(:)
      real(dp) :: sigma_max, sig_thresh, bnorm_val
      real(dp) :: nrm
      integer :: j, jj, i, idx_col, n_fix_v, n_fix_u
      integer, allocatable :: rec_u_idx(:), rec_v_idx(:), both_idx(:), nz_idx(:)
      integer, allocatable :: fix_v_idx(:), fix_u_idx(:)
      integer :: n_rec_u, n_rec_v, n_both, n_nz, n_fv, n_fu
      real(dp), allocatable :: BV_idx(:,:), BTU_idx(:,:), U_new(:,:), V_new(:,:)
      real(dp) :: dot_val
      ! GS completion
      real(dp), allocatable :: all_v_norms(:), all_u_norms(:)
      logical, allocatable :: good_mask(:)
      real(dp), allocatable :: V_gd(:,:), U_gd(:,:), vc(:), uc(:)
      integer, allocatable :: zero_v_mc(:), zero_u_mc(:)
      integer :: n_zero_v, n_zero_u, ngood, kk, isd
      integer :: iseed(4)

      if (n_mc <= 0) return

      allocate(nv_mc(n_mc), nu_mc(n_mc))
      allocate(max_norms(n_mc), norm_thresh(n_mc))
      allocate(v_good(n_mc), u_good(n_mc))
      allocate(scale_v_arr(n_mc), scale_u_arr(n_mc))

      ! Norms of V[:, mc] and U[:, mc]
      do j = 1, n_mc
         idx_col = mc(j)
         nv_mc(j) = vec_norm(n, V(1, idx_col))
         nu_mc(j) = vec_norm(n, U(1, idx_col))
      end do

      do j = 1, n_mc
         max_norms(j) = max(nv_mc(j), nu_mc(j))
         norm_thresh(j) = max_norms(j) * 1.0d-4
         v_good(j) = (nv_mc(j) > norm_thresh(j)) .and. (nv_mc(j) > SAFMIN)
         u_good(j) = (nu_mc(j) > norm_thresh(j)) .and. (nu_mc(j) > SAFMIN)
         if (v_good(j)) then
            scale_v_arr(j) = 1.0_dp / max(nv_mc(j), SAFMIN)
         else
            scale_v_arr(j) = 0.0_dp
         end if
         if (u_good(j)) then
            scale_u_arr(j) = 1.0_dp / max(nu_mc(j), SAFMIN)
         else
            scale_u_arr(j) = 0.0_dp
         end if
      end do

      ! Apply scaling to V[:, mc], U[:, mc]
      do j = 1, n_mc
         idx_col = mc(j)
         do i = 1, n
            V(i, idx_col) = V(i, idx_col) * scale_v_arr(j)
            U(i, idx_col) = U(i, idx_col) * scale_u_arr(j)
         end do
      end do

      ! ----- Bv recovery -----
      sigma_max = 0.0_dp
      do j = 1, n
         if (sigma(j) > sigma_max) sigma_max = sigma(j)
      end do
      sig_thresh = real(n, dp) * EPS * sigma_max

      allocate(need_recovery(n_mc), recover_u(n_mc), recover_v(n_mc), both_good(n_mc))
      do j = 1, n_mc
         need_recovery(j) = sigma(mc(j)) > sig_thresh
         recover_u(j) = need_recovery(j) .and. (.not. u_good(j)) .and. v_good(j)
         recover_v(j) = need_recovery(j) .and. (.not. v_good(j)) .and. u_good(j)
         both_good(j) = need_recovery(j) .and. v_good(j) .and. u_good(j)
      end do

      ! Recover U
      n_rec_u = 0
      allocate(rec_u_idx(n_mc))
      do j = 1, n_mc
         if (recover_u(j)) then
            n_rec_u = n_rec_u + 1
            rec_u_idx(n_rec_u) = mc(j)
         end if
      end do
      if (n_rec_u > 0) then
         allocate(BV_idx(n, n_rec_u), U_new(n, n_rec_u))
         ! Build V[:, idx] subarray then multiply
         do jj = 1, n_rec_u
            idx_col = rec_u_idx(jj)
            do i = 1, n
               U_new(i, jj) = V(i, idx_col)   ! reuse U_new as scratch for V columns
            end do
         end do
         call bidiag_matvec_batch(n, d, e, U_new, n, n_rec_u, BV_idx, n)
         do jj = 1, n_rec_u
            idx_col = rec_u_idx(jj)
            do i = 1, n
               U_new(i, jj) = BV_idx(i, jj) / sigma(idx_col)
            end do
            nrm = max(vec_norm(n, U_new(1, jj)), SAFMIN)
            do i = 1, n
               U(i, idx_col) = U_new(i, jj) / nrm
            end do
         end do
         deallocate(BV_idx, U_new)
      end if
      deallocate(rec_u_idx)

      ! Recover V
      n_rec_v = 0
      allocate(rec_v_idx(n_mc))
      do j = 1, n_mc
         if (recover_v(j)) then
            n_rec_v = n_rec_v + 1
            rec_v_idx(n_rec_v) = mc(j)
         end if
      end do
      if (n_rec_v > 0) then
         allocate(BTU_idx(n, n_rec_v), V_new(n, n_rec_v))
         do jj = 1, n_rec_v
            idx_col = rec_v_idx(jj)
            do i = 1, n
               V_new(i, jj) = U(i, idx_col)
            end do
         end do
         call bidiagT_matvec_batch(n, d, e, V_new, n, n_rec_v, BTU_idx, n)
         do jj = 1, n_rec_v
            idx_col = rec_v_idx(jj)
            do i = 1, n
               V_new(i, jj) = BTU_idx(i, jj) / sigma(idx_col)
            end do
            nrm = max(vec_norm(n, V_new(1, jj)), SAFMIN)
            do i = 1, n
               V(i, idx_col) = V_new(i, jj) / nrm
            end do
         end do
         deallocate(BTU_idx, V_new)
      end if
      deallocate(rec_v_idx)

      ! Both-good check
      n_both = 0
      allocate(both_idx(n_mc))
      do j = 1, n_mc
         if (both_good(j)) then
            n_both = n_both + 1
            both_idx(n_both) = mc(j)
         end if
      end do
      if (n_both > 0) then
         allocate(BV_idx(n, n_both), BTU_idx(n, n_both))
         allocate(V_new(n, n_both), U_new(n, n_both))
         allocate(res_v_arr(n_both), res_u_arr(n_both), thresh_arr(n_both))
         allocate(fix_v(n_both), fix_u(n_both))
         do jj = 1, n_both
            idx_col = both_idx(jj)
            do i = 1, n
               U_new(i, jj) = V(i, idx_col)
               V_new(i, jj) = U(i, idx_col)
            end do
         end do
         call bidiag_matvec_batch (n, d, e, U_new, n, n_both, BV_idx, n)
         call bidiagT_matvec_batch(n, d, e, V_new, n, n_both, BTU_idx, n)
         bnorm_val = 0.0_dp
         do i = 0, n - 1
            if (abs(d(i)) > bnorm_val) bnorm_val = abs(d(i))
         end do
         if (n > 1) then
            do i = 0, n - 2
               if (abs(e(i)) > bnorm_val) bnorm_val = abs(e(i))
            end do
         end if
         allocate(tmp_n(n))
         do jj = 1, n_both
            idx_col = both_idx(jj)
            ! res_v = ||BV[:,jj] - sigma * U[:,idx]||
            do i = 1, n
               tmp_n(i) = BV_idx(i, jj) - sigma(idx_col) * U(i, idx_col)
            end do
            res_v_arr(jj) = vec_norm(n, tmp_n)
            do i = 1, n
               tmp_n(i) = BTU_idx(i, jj) - sigma(idx_col) * V(i, idx_col)
            end do
            res_u_arr(jj) = vec_norm(n, tmp_n)
            thresh_arr(jj) = 10.0_dp * real(n, dp) * EPS * max(sigma(idx_col), bnorm_val)
            fix_v(jj) = (res_v_arr(jj) > thresh_arr(jj)) .and. (res_v_arr(jj) > res_u_arr(jj))
            fix_u(jj) = (res_u_arr(jj) > thresh_arr(jj)) .and. (res_u_arr(jj) > res_v_arr(jj)) .and. (.not. fix_v(jj))
         end do
         deallocate(tmp_n)
         ! Apply fix_v (recompute V from U)
         do jj = 1, n_both
            idx_col = both_idx(jj)
            if (fix_v(jj)) then
               do i = 1, n
                  V_new(i, jj) = BTU_idx(i, jj) / sigma(idx_col)
               end do
               nrm = max(vec_norm(n, V_new(1, jj)), SAFMIN)
               do i = 1, n
                  V(i, idx_col) = V_new(i, jj) / nrm
               end do
            end if
            if (fix_u(jj)) then
               do i = 1, n
                  U_new(i, jj) = BV_idx(i, jj) / sigma(idx_col)
               end do
               nrm = max(vec_norm(n, U_new(1, jj)), SAFMIN)
               do i = 1, n
                  U(i, idx_col) = U_new(i, jj) / nrm
               end do
            end if
         end do
         deallocate(BV_idx, BTU_idx, V_new, U_new)
         deallocate(res_v_arr, res_u_arr, thresh_arr, fix_v, fix_u)
      end if
      deallocate(both_idx)

      ! ----- Sign fix -----
      n_nz = 0
      allocate(nz_idx(n_mc))
      do j = 1, n_mc
         if (sigma(mc(j)) >= SAFMIN) then
            n_nz = n_nz + 1
            nz_idx(n_nz) = mc(j)
         end if
      end do
      if (n_nz > 0) then
         allocate(BV_idx(n, n_nz))
         allocate(U_new(n, n_nz))
         do jj = 1, n_nz
            idx_col = nz_idx(jj)
            do i = 1, n
               U_new(i, jj) = V(i, idx_col)
            end do
         end do
         call bidiag_matvec_batch(n, d, e, U_new, n, n_nz, BV_idx, n)
         do jj = 1, n_nz
            idx_col = nz_idx(jj)
            dot_val = 0.0_dp
            do i = 1, n
               dot_val = dot_val + BV_idx(i, jj) * U(i, idx_col)
            end do
            if (dot_val < 0.0_dp) then
               do i = 1, n
                  U(i, idx_col) = -U(i, idx_col)
               end do
            end if
         end do
         deallocate(BV_idx, U_new)
      end if
      deallocate(nz_idx)

      ! ----- Gram-Schmidt completion -----
      ! Recompute all column norms (because Bv recovery may have changed them)
      allocate(all_v_norms(n), all_u_norms(n))
      do j = 1, n
         all_v_norms(j) = vec_norm(n, V(1, j))
         all_u_norms(j) = vec_norm(n, U(1, j))
      end do
      ! Find zero columns within multi-block set
      n_zero_v = 0
      n_zero_u = 0
      allocate(zero_v_mc(n_mc), zero_u_mc(n_mc))
      do j = 1, n_mc
         idx_col = mc(j)
         if (all_v_norms(idx_col) < 0.5_dp) then
            n_zero_v = n_zero_v + 1
            zero_v_mc(n_zero_v) = idx_col
         end if
         if (all_u_norms(idx_col) < 0.5_dp) then
            n_zero_u = n_zero_u + 1
            zero_u_mc(n_zero_u) = idx_col
         end if
      end do

      if (n_zero_v > 0 .or. n_zero_u > 0) then
         allocate(good_mask(n))
         allocate(vc(n), uc(n))
      end if
      ! V completion
      do jj = 1, n_zero_v
         idx_col = zero_v_mc(jj)
         do j = 1, n
            good_mask(j) = (j /= idx_col) .and. (all_v_norms(j) > 0.5_dp)
         end do
         ngood = 0
         do j = 1, n
            if (good_mask(j)) ngood = ngood + 1
         end do
         allocate(V_gd(n, ngood))
         kk = 0
         do j = 1, n
            if (good_mask(j)) then
               kk = kk + 1
               do i = 1, n
                  V_gd(i, kk) = V(i, j)
               end do
            end if
         end do
         ! Use DLARNV with same seed as Python: (1, 3, 5, (2*j+1)&4095 |1)
         ! Note: idx_col is 1-based here (Fortran), Python uses 0-based.
         isd = (2 * (idx_col - 1) + 1)
         isd = iand(isd, 4095)
         isd = ior(isd, 1)
         iseed = (/ 1, 3, 5, isd /)
         call dlarnv(3, iseed, n, vc)
         ! Two-pass MGS
         allocate(tmp_n(ngood), tmp_n2(n))
         do kk = 1, ngood
            tmp_n(kk) = 0.0_dp
            do i = 1, n
               tmp_n(kk) = tmp_n(kk) + V_gd(i, kk) * vc(i)
            end do
         end do
         do i = 1, n
            tmp_n2(i) = 0.0_dp
            do kk = 1, ngood
               tmp_n2(i) = tmp_n2(i) + V_gd(i, kk) * tmp_n(kk)
            end do
         end do
         do i = 1, n
            vc(i) = vc(i) - tmp_n2(i)
         end do
         do kk = 1, ngood
            tmp_n(kk) = 0.0_dp
            do i = 1, n
               tmp_n(kk) = tmp_n(kk) + V_gd(i, kk) * vc(i)
            end do
         end do
         do i = 1, n
            tmp_n2(i) = 0.0_dp
            do kk = 1, ngood
               tmp_n2(i) = tmp_n2(i) + V_gd(i, kk) * tmp_n(kk)
            end do
         end do
         do i = 1, n
            vc(i) = vc(i) - tmp_n2(i)
         end do
         deallocate(tmp_n, tmp_n2)
         nrm = vec_norm(n, vc)
         if (nrm > 1.0d-14) then
            do i = 1, n
               V(i, idx_col) = vc(i) / nrm
            end do
         end if
         deallocate(V_gd)
      end do
      ! U completion
      do jj = 1, n_zero_u
         idx_col = zero_u_mc(jj)
         do j = 1, n
            good_mask(j) = (j /= idx_col) .and. (all_u_norms(j) > 0.5_dp)
         end do
         ngood = 0
         do j = 1, n
            if (good_mask(j)) ngood = ngood + 1
         end do
         allocate(U_gd(n, ngood))
         kk = 0
         do j = 1, n
            if (good_mask(j)) then
               kk = kk + 1
               do i = 1, n
                  U_gd(i, kk) = U(i, j)
               end do
            end if
         end do
         isd = (2 * (idx_col - 1) + 3)
         isd = iand(isd, 4095)
         isd = ior(isd, 1)
         iseed = (/ 2, 4, 6, isd /)
         call dlarnv(3, iseed, n, uc)
         allocate(tmp_n(ngood), tmp_n2(n))
         do kk = 1, ngood
            tmp_n(kk) = 0.0_dp
            do i = 1, n
               tmp_n(kk) = tmp_n(kk) + U_gd(i, kk) * uc(i)
            end do
         end do
         do i = 1, n
            tmp_n2(i) = 0.0_dp
            do kk = 1, ngood
               tmp_n2(i) = tmp_n2(i) + U_gd(i, kk) * tmp_n(kk)
            end do
         end do
         do i = 1, n
            uc(i) = uc(i) - tmp_n2(i)
         end do
         do kk = 1, ngood
            tmp_n(kk) = 0.0_dp
            do i = 1, n
               tmp_n(kk) = tmp_n(kk) + U_gd(i, kk) * uc(i)
            end do
         end do
         do i = 1, n
            tmp_n2(i) = 0.0_dp
            do kk = 1, ngood
               tmp_n2(i) = tmp_n2(i) + U_gd(i, kk) * tmp_n(kk)
            end do
         end do
         do i = 1, n
            uc(i) = uc(i) - tmp_n2(i)
         end do
         deallocate(tmp_n, tmp_n2)
         nrm = vec_norm(n, uc)
         if (nrm > 1.0d-14) then
            do i = 1, n
               U(i, idx_col) = uc(i) / nrm
            end do
         end if
         deallocate(U_gd)
      end do
      if (allocated(good_mask)) deallocate(good_mask)
      if (allocated(vc))   deallocate(vc)
      if (allocated(uc))   deallocate(uc)
      deallocate(zero_v_mc, zero_u_mc, all_v_norms, all_u_norms)
      deallocate(nv_mc, nu_mc, max_norms, norm_thresh, v_good, u_good)
      deallocate(scale_v_arr, scale_u_arr)
      deallocate(need_recovery, recover_u, recover_v, both_good)
   end subroutine post_process

end module mr3gk_postproc
