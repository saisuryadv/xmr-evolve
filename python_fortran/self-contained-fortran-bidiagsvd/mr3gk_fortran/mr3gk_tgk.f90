! ===========================================================
! Module: mr3gk_tgk
! Solve T_GK eigenproblem by calling DLAXRE + DLAXRV (existing XMR Fortran).
! Mirrors _solve_tgk_block_xmr and mr3_tgk_multiblock from mr3_gk.py.
! ===========================================================
module mr3gk_tgk
   use mr3gk_consts
   use mr3gk_utils
   implicit none

   external :: dlaxre, dlaxrv, wsreq_xrv

contains

   ! Solve symmetric tridiagonal eigenproblem T = (bc, bte) of size m,
   ! returning the n_nonneg = (m+1)/2 (or bk_hint, if > 0) largest
   ! eigenvalues w(1:n_nonneg) and their eigenvectors Z(1:m, 1:n_nonneg).
   !
   ! This wraps the same dlaxre + dlaxrv path that the C wrapper uses
   ! in xmr_eigenvectors(). Bidirectional bit-equivalence with Python
   ! is guaranteed because we reuse the very same Fortran subroutines.
   subroutine solve_tgk_block(bc, bte, m, bk_hint, w, Z, ldz, info_out)
      integer, intent(in) :: m, bk_hint, ldz
      real(dp), intent(in) :: bc(0:m-1)         ! T_GK diagonal (zeros)
      real(dp), intent(in) :: bte(0:m-2)        ! T_GK off-diagonals
      real(dp), intent(out) :: w(*), Z(ldz, *)
      integer, intent(out) :: info_out

      integer :: n_nonneg, wil, wiu
      real(dp), allocatable :: all_evals(:)
      real(dp) :: spdiam, gaptol_loc, gl_g, gu_g, abserr_g, tau_re
      integer  :: nfound, info_be
      ! XMR scratch arrays
      real(dp), allocatable :: d_arr(:), e_arr(:)
      real(dp), allocatable :: repr_buf(:), rwork_re(:)
      integer,  allocatable :: repi_buf(:), ewl_ae(:)
      real(dp), allocatable :: ewl_lu(:)
      real(dp), allocatable :: rwork_v(:)
      integer,  allocatable :: iwork_v(:), isuppz(:)
      real(dp), allocatable :: w_loc(:), z_loc(:)
      integer  :: reqr, reqi, j, i
      character :: emode

      info_out = 0
      n_nonneg = (m + 1) / 2
      if (bk_hint > 0) n_nonneg = bk_hint
      wil = m - n_nonneg + 1
      wiu = m

      ! Compute all m eigenvalues via DSTEBZ (matches Python's bisect_evals 1..m).
      allocate(all_evals(m))
      call bisect_evals(bc, bte, m, 1, m, 0.0_dp, all_evals, nfound, info_be)
      if (info_be /= 0 .or. nfound /= m) then
         info_out = -100 - info_be
         deallocate(all_evals)
         return
      end if
      ! Sort ascending (DSTEBZ with ORDER='E' already returns sorted).
      spdiam = all_evals(m) - all_evals(1)

      ! ----- dlaxre setup: identical to xmr_wrapper.c -----
      allocate(d_arr(m), e_arr(m))
      do i = 1, m
         d_arr(i) = 0.0_dp
      end do
      do i = 1, m - 1
         e_arr(i) = bte(i-1)
      end do
      e_arr(m) = 0.0_dp

      gl_g = all_evals(1) - 0.01_dp * abs(all_evals(1)) - 1.0d-20
      gu_g = all_evals(m) + 0.01_dp * abs(all_evals(m)) + 1.0d-20
      abserr_g = 0.0_dp

      gaptol_loc = max(GAPTOL, 0.02_dp / real(n_nonneg, dp))

      allocate(repr_buf(4*m + 3), repi_buf(6 + m + m/2))
      allocate(ewl_ae(2*m), ewl_lu(2*m))
      allocate(rwork_re(6*m + 10))
      tau_re = 0.0_dp

      emode = 'o'
      call dlaxre(m, d_arr, e_arr, gl_g, gu_g, abserr_g, gaptol_loc, &
                  repr_buf, repi_buf, tau_re, &
                  ewl_ae, ewl_lu, emode, &
                  rwork_re, info_out)
      deallocate(rwork_re, d_arr)

      if (info_out /= 0) then
         deallocate(e_arr, repr_buf, repi_buf, ewl_ae, ewl_lu, all_evals)
         info_out = -10 - info_out
         return
      end if

      ! ----- dlaxrv with proper workspace -----
      reqr = -1
      reqi = -1
      call wsreq_xrv(m, reqr, reqi)
      allocate(rwork_v(reqr + 1), iwork_v(reqi + 1))
      allocate(isuppz(2*n_nonneg))
      allocate(w_loc(n_nonneg), z_loc(m * n_nonneg))

      call dlaxrv(m, e_arr, repr_buf, repi_buf, ewl_ae, ewl_lu, &
                  wil, wiu, spdiam, gaptol_loc, &
                  w_loc, z_loc, m, isuppz, &
                  rwork_v, reqr, iwork_v, reqi, info_out)

      ! Add back the dlaxre shift and copy results
      do j = 1, n_nonneg
         w(j) = w_loc(j) + tau_re
         do i = 1, m
            Z(i, j) = z_loc(i + m*(j-1))
         end do
      end do

      deallocate(e_arr, repr_buf, repi_buf, ewl_ae, ewl_lu)
      deallocate(rwork_v, iwork_v, isuppz, w_loc, z_loc, all_evals)
   end subroutine solve_tgk_block

   ! Run MR³ on a single bidiagonal block [bbeg..bend] of (d_bidiag, e_bidiag).
   ! Returns w(1:bk) eigenvalues (singular values) and Z(1:2*bk, 1:bk) eigenvectors.
   ! Mirrors mr3_tgk_multiblock in mr3_gk.py.
   subroutine mr3_tgk_multiblock(d_bidiag, e_bidiag, n, bbeg, bend, w, Z, ldz, info_out)
      integer, intent(in) :: n, bbeg, bend, ldz
      real(dp), intent(in) :: d_bidiag(0:n-1), e_bidiag(0:max(n-2,0))
      real(dp), intent(out) :: w(*), Z(ldz, *)
      integer, intent(out) :: info_out

      integer :: bk, m2k, i
      real(dp), allocatable :: bc(:), bte(:)

      bk = bend - bbeg + 1
      m2k = 2 * bk
      allocate(bc(0:m2k-1), bte(0:m2k-2))
      do i = 0, m2k - 1
         bc(i) = 0.0_dp
      end do
      ! bte[0::2] = |d|, bte[1::2] = |e|
      do i = 0, bk - 1
         bte(2*i) = abs(d_bidiag(bbeg + i))
      end do
      do i = 0, bk - 2
         bte(2*i + 1) = abs(e_bidiag(bbeg + i))
      end do

      call solve_tgk_block(bc, bte, m2k, bk, w, Z, ldz, info_out)

      deallocate(bc, bte)
   end subroutine mr3_tgk_multiblock

end module mr3gk_tgk
