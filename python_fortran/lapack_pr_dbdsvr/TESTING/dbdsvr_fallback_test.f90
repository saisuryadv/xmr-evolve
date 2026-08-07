program dbdsvr_fallback_test
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: dp = real64, n = 4
  integer, parameter :: lwork = 5*n*n + 37*n, liwork = 20*n
  integer :: fault_mode
  common /dbdsvr_fault/ fault_mode

  call check_fallback(1, 2)
  call check_fallback(2, 3)
  write(*,'(A)') 'DBDSVR fallback tests passed.'

contains

  subroutine check_fallback(mode, expected_path)
    integer, intent(in) :: mode, expected_path
    integer :: info, iwork(liwork), ns
    real(dp) :: d(n), e(n), s(n), work(lwork), z(2*n,n)

    d = [4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
    e = [0.25_dp, 0.20_dp, 0.10_dp, 0.0_dp]
    iwork = 0
    fault_mode = mode
    call dbdsvr('U', 'V', 'A', n, d, e, 0.0_dp, 0.0_dp, 0, 0, &
         ns, s, z, 2*n, work, lwork, iwork, liwork, info)

    if (info /= 0) error stop 'fallback returned nonzero INFO'
    if (ns /= n) error stop 'fallback returned incomplete spectrum'
    if (iwork(2) /= expected_path) error stop 'wrong fallback path'
    if (any(.not. ieee_is_finite(s))) error stop 'nonfinite fallback values'
    if (any(.not. ieee_is_finite(z))) error stop 'nonfinite fallback vectors'
  end subroutine check_fallback

end program dbdsvr_fallback_test

subroutine dbdsvdmr3(jobz, uplo, n, d, e, s, u, ldu, vt, ldvt, &
     m, work, lwork, iwork, liwork, info)
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  implicit none

  integer, parameter :: dp = real64
  character, intent(in) :: jobz, uplo
  integer, intent(in) :: n, ldu, ldvt, lwork, liwork
  integer, intent(out) :: m, info
  integer, intent(inout) :: iwork(*)
  real(dp), intent(inout) :: d(*), e(*), work(*)
  real(dp), intent(out) :: s(*), u(ldu,*), vt(ldvt,*)
  integer :: col, fault_mode, row
  common /dbdsvr_fault/ fault_mode

  if (fault_mode == 1) then
    m = 0
    info = 1
    return
  end if

  m = n
  info = 0
  do col = 1, n
    s(col) = real(col, dp)
    do row = 1, n
      u(row,col) = 0.0_dp
      vt(row,col) = 0.0_dp
    end do
    u(col,col) = 1.0_dp
    vt(col,col) = 1.0_dp
  end do
  u(1,1) = ieee_value(0.0_dp, ieee_quiet_nan)
end subroutine dbdsvdmr3
