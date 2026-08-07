program dbdsv_solver_probe
  use, intrinsic :: iso_fortran_env, only: int64, real64
  implicit none

  character(len=4096) :: input_path
  character(len=16) :: method
  integer :: i, idx, info, ios, liwork, lwork, m, n
  integer(int64) :: count_begin, count_end, count_rate
  real(real64) :: elapsed, file_d, file_e
  real(real64), allocatable :: d(:), e(:), s(:), u(:,:), vt(:,:), work(:)
  integer, allocatable :: iwork(:)
  real(real64) :: dummy_q(1)
  integer :: dummy_iq(1)

  call get_command_argument(1, input_path)
  call get_command_argument(2, method)
  if (len_trim(input_path) == 0 .or. len_trim(method) == 0) then
    write(*, '(A)') 'usage: dbdsv_solver_probe MATRIX.dat mr3|dbdsdc'
    stop 2
  end if

  open(unit=10, file=trim(input_path), status='old', action='read', &
       iostat=ios)
  if (ios /= 0) then
    write(*, '(A,A)') 'cannot open ', trim(input_path)
    stop 2
  end if
  read(10, *, iostat=ios) n
  if (ios /= 0 .or. n < 1) then
    write(*, '(A,I0)') 'bad matrix order: ', n
    stop 2
  end if

  if (trim(method) == 'mr3') then
    lwork = max(1, 2*n*n + 100*n)
    liwork = max(1, 30*n)
  else if (trim(method) == 'dbdsdc') then
    lwork = max(1, 3*n*n + 4*n)
    liwork = max(1, 8*n)
  else
    write(*, '(A,A)') 'unknown method: ', trim(method)
    stop 2
  end if

  allocate(d(n), e(n), s(n), u(n,n), vt(n,n))
  allocate(work(lwork), iwork(liwork))
  do i = 1, n
    read(10, *, iostat=ios) idx, file_d, file_e
    if (ios /= 0) then
      write(*, '(A,I0)') 'bad matrix row: ', i
      stop 2
    end if
    d(i) = file_d
    e(i) = file_e
  end do
  close(10)
  e(n) = 0.0_real64

  call system_clock(count_begin, count_rate)
  if (trim(method) == 'mr3') then
    call dbdsvdmr3('V', 'U', n, d, e, s, u, n, vt, n, m, &
                   work, lwork, iwork, liwork, info)
  else
    call dbdsdc('U', 'I', n, d, e, u, n, vt, n, dummy_q, dummy_iq, &
                work, iwork, info)
    m = n
    if (info == 0) s = d
  end if
  call system_clock(count_end)
  elapsed = real(count_end - count_begin, real64) / real(count_rate, real64)

  write(*, '(A,A,A,I0,A,I0,A,ES14.6)') 'method=', trim(method), &
       ' n=', n, ' info=', info, ' solver_seconds=', elapsed
  if (info == 0) then
    write(*, '(A,ES24.16,A,ES24.16)') 'sigma_first=', s(1), &
         ' sigma_last=', s(m)
  end if
end program dbdsv_solver_probe
