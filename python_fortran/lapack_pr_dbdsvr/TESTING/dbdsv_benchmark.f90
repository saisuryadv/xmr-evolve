program dbdsv_benchmark
  use iso_fortran_env, only: int64
  implicit none

  integer, parameter :: dp = kind(1.0d0)
#ifndef BENCH_MAXN
#define BENCH_MAXN 400
#endif
#ifndef BENCH_NTYPES
#define BENCH_NTYPES 16
#endif
  integer, parameter :: maxn = BENCH_MAXN, ntypes = BENCH_NTYPES, nsizes = 11
  integer, parameter :: email_nsizes = 5
  integer, parameter :: lwork = 5*maxn*maxn + 40*maxn + 4096
  integer, parameter :: liwork = 20*maxn + 4096
  integer, parameter :: default_max_reps = 1024
  real(dp), parameter :: default_timing_seconds = 0.020_dp
  integer, parameter :: sizes(nsizes) = &
       [2, 10, 30, 60, 100, 200, 400, 500, 1000, 2000, 3000]
  character(len=7), parameter :: methods(4) = &
       [character(len=7) :: 'DBDSDC ', 'DBDSQR ', 'DBDSVDX', 'DBDSVR ']
  character(len=1), parameter :: ranges(3) = ['A', 'I', 'V']

  integer :: active_nsizes, i, info, isize, itype, meth, irange, il, iu, ns
  integer :: iseed(4), reps, max_reps, target_type, target_n, ios
  integer :: backend_path
  character(len=1) :: uplo, target_range
  character(len=16) :: target_method, arg
  character(len=4096) :: dump_path
  logical :: dump_written, timing_only
  real(dp) :: anorm, vl, vu, residual, u_orth, vt_orth, seconds
  real(dp) :: min_timing_seconds
  real(dp) :: d(maxn), e(maxn), dtmp(maxn)
  real(dp) :: s(maxn), sref(maxn), z(2*maxn,maxn)
  real(dp) :: u(maxn,maxn), vt(maxn,maxn), a(maxn,maxn)
  real(dp) :: work(lwork)
  integer :: iwork(liwork)

  real(dp), external :: dlamch
  external :: dbdt04, dbsvrmg, dort01

  iseed = [1, 2, 3, 4]
  target_type = 0
  target_n = 0
  target_range = '*'
  target_method = '*'
  dump_path = ''
  dump_written = .false.
  timing_only = .false.
  min_timing_seconds = default_timing_seconds
  max_reps = default_max_reps
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    read(arg, *, iostat=ios) target_type
    if (ios /= 0) error stop 'invalid type filter'
  end if
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg)
    read(arg, *, iostat=ios) target_n
    if (ios /= 0) error stop 'invalid N filter'
  end if
  if (command_argument_count() >= 3) call get_command_argument(3, target_range)
  if (command_argument_count() >= 4) call get_command_argument(4, target_method)
  if (command_argument_count() >= 5) then
    call get_command_argument(5, arg)
    read(arg, *, iostat=ios) min_timing_seconds
    if (ios /= 0 .or. min_timing_seconds <= 0.0_dp) &
         error stop 'invalid timing target'
  end if
  if (command_argument_count() >= 6) then
    call get_command_argument(6, arg)
    read(arg, *, iostat=ios) max_reps
    if (ios /= 0 .or. max_reps < 1) error stop 'invalid repetition cap'
  end if
  if (command_argument_count() >= 7) then
    call get_command_argument(7, dump_path)
  end if
  if (command_argument_count() >= 8) then
    call get_command_argument(8, arg)
    if (trim(arg) == 'TIMING_ONLY') then
      timing_only = .true.
    else
      error stop 'eighth argument must be TIMING_ONLY'
    end if
  end if
  if (target_n > maxn) &
       error stop 'target N exceeds compile-time BENCH_MAXN'
  active_nsizes = email_nsizes
  if (target_n > sizes(email_nsizes)) then
    active_nsizes = 0
    do isize = email_nsizes + 1, nsizes
      if (target_n == sizes(isize)) active_nsizes = isize
    end do
    if (active_nsizes == 0) &
         error stop 'unsupported large-N filter'
  end if
  write(*,'(A)') 'method,type,n,range,ns,residual,u_orthogonality,' // &
       'vt_orthogonality,orthogonality,time_seconds,info,repetitions,' // &
       'backend_path'

  do isize = 1, active_nsizes
    do itype = 1, ntypes
      call dbsvrmg(itype, sizes(isize), iseed, d, e, uplo, a, maxn, &
           dtmp, work, info)
      if (info /= 0) then
        do meth = 1, size(methods)
          do irange = 1, size(ranges)
            call emit_row(methods(meth), itype, sizes(isize), &
                 ranges(irange), 0, huge(1.0_dp), huge(1.0_dp), &
                 huge(1.0_dp), huge(1.0_dp), 0.0_dp, info, 0, 0)
          end do
        end do
        cycle
      end if
      if (target_n > 0 .and. sizes(isize) /= target_n) cycle
      if (target_type > 0 .and. itype /= target_type) cycle

      anorm = abs(d(1))
      do i = 2, sizes(isize)
        anorm = max(anorm, abs(d(i)) + abs(e(i-1)))
      end do
      if (anorm <= 0.0_dp) anorm = 1.0_dp

      il = max(1, sizes(isize)/4 + 1)
      iu = min(sizes(isize), (3*sizes(isize))/4)
      if (iu < il) then
        il = 1
        iu = sizes(isize)
      end if
      call reference_spectrum(uplo, sizes(isize), d, e, sref, info)
      if (info /= 0) then
        write(0,'(A,I0,A,I0,A,I0)') 'reference DBDSDC failed: N=', &
             sizes(isize), ' type=', itype, ' INFO=', info
        stop 2
      end if
      call value_bounds(sizes(isize), sref, il, iu, anorm, vl, vu)

      do meth = 1, size(methods)
        if (trim(target_method) /= '*' .and. &
             trim(methods(meth)) /= trim(target_method)) cycle
        do irange = 1, size(ranges)
          if (target_range /= '*' .and. ranges(irange) /= target_range) cycle
          call benchmark_case(methods(meth), uplo, ranges(irange), &
               sizes(isize), d, e, vl, vu, il, iu, ns, s, z, seconds, &
               reps, info)
          backend_path = 0
          if (info == 0 .and. trim(methods(meth)) == 'DBDSVR') &
               backend_path = iwork(2)
          if (info == 0 .and. ns > 0 .and. .not. timing_only) then
            call unpack_vectors(sizes(isize), ns, z, u, vt)
            call dbdt04(uplo, sizes(isize), d, e, s, ns, u, maxn, vt, &
                 maxn, work, residual)
            call dort01('Columns', sizes(isize), ns, u, maxn, work, &
                 lwork, u_orth)
            call dort01('Rows', ns, sizes(isize), vt, maxn, work, &
                 lwork, vt_orth)
#ifndef PAPER29
            if (trim(methods(meth)) == 'DBDSVR' .and. itype == 2) then
              if (residual /= 0.0_dp .or. u_orth /= 0.0_dp .or. &
                   vt_orth /= 0.0_dp) &
                   error stop 'DBDSVR identity metrics must be exact zero'
            end if
#endif
          else if (info /= 0 .or. ns <= 0) then
            residual = huge(1.0_dp)
            u_orth = huge(1.0_dp)
            vt_orth = huge(1.0_dp)
          else
            residual = -1.0_dp
            u_orth = -1.0_dp
            vt_orth = -1.0_dp
          end if
          call emit_row(methods(meth), itype, sizes(isize), &
               ranges(irange), ns, residual, u_orth, vt_orth, &
               max(u_orth, vt_orth), seconds, info, reps, backend_path)
          if (len_trim(dump_path) > 0) then
            call dump_solution(trim(dump_path), itype, sizes(isize), &
                 ranges(irange), methods(meth), ns, info, s, z, &
                 .not. dump_written)
            dump_written = .true.
          end if
        end do
      end do
    end do
  end do

contains

  subroutine dump_solution(path, itype, n, range, method, nfound, ierr, &
       sigma, packed, first)
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: method
    character(len=1), intent(in) :: range
    integer, intent(in) :: itype, n, nfound, ierr
    real(dp), intent(in) :: sigma(maxn), packed(2*maxn,maxn)
    logical, intent(in) :: first
    integer :: col, unit

    if (first) then
      open(newunit=unit, file=path, access='stream', form='unformatted', &
           status='replace', action='write')
    else
      open(newunit=unit, file=path, access='stream', form='unformatted', &
           status='old', position='append', action='write')
    end if
    write(unit) itype, n, nfound, ierr, range, method
    if (nfound > 0) then
      write(unit) sigma(1:nfound)
      do col = 1, nfound
        write(unit) packed(1:2*n,col)
      end do
    end if
    close(unit)
  end subroutine dump_solution

  subroutine benchmark_case(method, ul, range, n, bd, be, vlow, vhigh, &
       idxlo, idxhi, nfound, sout, zout, elapsed, nreps, ierr)
    character(len=*), intent(in) :: method
    character(len=1), intent(in) :: ul, range
    integer, intent(in) :: n, idxlo, idxhi
    real(dp), intent(in) :: bd(maxn), be(maxn), vlow, vhigh
    integer, intent(out) :: nfound, nreps, ierr
    real(dp), intent(out) :: sout(maxn), zout(2*maxn,maxn), elapsed
    integer :: rep, trial_info, trial_ns
    integer(int64) :: clock0, clock1, clock_rate
    real(dp) :: trial_s(maxn), trial_z(2*maxn,maxn), trial_time, total
    real(dp) :: solver_total

    call solve_case(method, ul, range, n, bd, be, vlow, vhigh, idxlo, &
         idxhi, nfound, sout, zout, elapsed, ierr)
    nreps = 1
    if (ierr /= 0) return

    do
      solver_total = 0.0_dp
      call system_clock(clock0, clock_rate)
      do rep = 1, nreps
        call solve_case(method, ul, range, n, bd, be, vlow, vhigh, &
             idxlo, idxhi, trial_ns, trial_s, trial_z, trial_time, &
             trial_info)
        if (trial_info /= 0) then
          ierr = trial_info
          return
        end if
        solver_total = solver_total + trial_time
      end do
      call system_clock(clock1)
      total = real(clock1-clock0, dp) / real(clock_rate, dp)
      if (total >= min_timing_seconds .or. nreps >= max_reps) exit
      nreps = min(2*nreps, max_reps)
    end do
    if (timing_only) then
      elapsed = solver_total / real(nreps, dp)
    else
      elapsed = total / real(nreps, dp)
    end if
  end subroutine benchmark_case

  subroutine solve_case(method, ul, range, n, bd, be, vlow, vhigh, &
       idxlo, idxhi, nfound, sout, zout, elapsed, ierr)
    character(len=*), intent(in) :: method
    character(len=1), intent(in) :: ul, range
    integer, intent(in) :: n, idxlo, idxhi
    real(dp), intent(in) :: bd(maxn), be(maxn), vlow, vhigh
    integer, intent(out) :: nfound, ierr
    real(dp), intent(out) :: sout(maxn), zout(2*maxn,maxn), elapsed
    integer(int64) :: clock0, clock1, clock_rate
    real(dp) :: dc(maxn), ec(maxn)
    real(dp) :: ufull(maxn,maxn), vtfull(maxn,maxn)
    real(dp) :: q(maxn*maxn)
    integer :: iq(maxn*maxn)

    dc(1:n) = bd(1:n)
    if (n > 1) ec(1:n-1) = be(1:n-1)
    ec(n) = 0.0_dp
    sout(1:n) = 0.0_dp
    zout(1:2*n,1:n) = 0.0_dp
    ierr = 0
    nfound = 0
    call system_clock(clock0, clock_rate)

    select case (trim(method))
    case ('DBDSDC')
      call dbdsdc(ul, 'I', n, dc, ec, ufull, maxn, vtfull, maxn, q, iq, &
           work, iwork, ierr)
      call system_clock(clock1)
      if (ierr == 0) call select_full(range, n, dc, ufull, vtfull, &
           vlow, vhigh, idxlo, idxhi, nfound, sout, zout)
    case ('DBDSQR')
      call dlaset('Full', n, n, 0.0_dp, 1.0_dp, ufull, maxn)
      call dlaset('Full', n, n, 0.0_dp, 1.0_dp, vtfull, maxn)
      call dbdsqr(ul, n, n, n, 0, dc, ec, vtfull, maxn, ufull, maxn, &
           q, 1, work, ierr)
      call system_clock(clock1)
      if (ierr == 0) call select_full(range, n, dc, ufull, vtfull, &
           vlow, vhigh, idxlo, idxhi, nfound, sout, zout)
    case ('DBDSVDX')
      call dbdsvdx(ul, 'V', range, n, dc, ec, vlow, vhigh, idxlo, &
           idxhi, nfound, sout, zout, 2*maxn, work, iwork, ierr)
      call system_clock(clock1)
    case ('DBDSVR')
      call dbdsvr(ul, 'V', range, n, dc, ec, vlow, vhigh, idxlo, &
           idxhi, nfound, sout, zout, 2*maxn, work, lwork, iwork, &
           liwork, ierr)
      call system_clock(clock1)
    case default
      call system_clock(clock1)
      ierr = -999
    end select
    elapsed = real(clock1-clock0, dp) / real(clock_rate, dp)
  end subroutine solve_case

  subroutine select_full(range, n, sfull, ufull, vtfull, vlow, vhigh, &
       idxlo, idxhi, nfound, sout, zout)
    character(len=1), intent(in) :: range
    integer, intent(in) :: n, idxlo, idxhi
    real(dp), intent(in) :: sfull(maxn), ufull(maxn,maxn)
    real(dp), intent(in) :: vtfull(maxn,maxn), vlow, vhigh
    integer, intent(out) :: nfound
    real(dp), intent(out) :: sout(maxn), zout(2*maxn,maxn)
    integer :: first, last, j, k

    nfound = 0
    first = 1
    last = n
    if (range == 'I') then
      first = idxlo
      last = idxhi
    end if
    do j = first, last
      if (range == 'V') then
        if (sfull(j) < vlow .or. sfull(j) >= vhigh) cycle
      end if
      nfound = nfound + 1
      sout(nfound) = sfull(j)
      do k = 1, n
        zout(k,nfound) = ufull(k,j)
        zout(n+k,nfound) = vtfull(j,k)
      end do
    end do
  end subroutine select_full

  subroutine reference_spectrum(ul, n, bd, be, sref, ierr)
    character(len=1), intent(in) :: ul
    integer, intent(in) :: n
    real(dp), intent(in) :: bd(maxn), be(maxn)
    real(dp), intent(out) :: sref(maxn)
    integer, intent(out) :: ierr
    real(dp) :: dc(maxn), ec(maxn), dummy_u(1,1), dummy_vt(1,1)
    real(dp) :: q(1), w(4*maxn)
    integer :: iq(1), iw(8*maxn)

    dc(1:n) = bd(1:n)
    if (n > 1) ec(1:n-1) = be(1:n-1)
    call dbdsdc(ul, 'N', n, dc, ec, dummy_u, 1, dummy_vt, 1, q, iq, &
         w, iw, ierr)
    sref(1:n) = dc(1:n)
  end subroutine reference_spectrum

  subroutine value_bounds(n, sv, idxlo, idxhi, matrix_norm, vlow, vhigh)
    integer, intent(in) :: n, idxlo, idxhi
    real(dp), intent(in) :: sv(maxn), matrix_norm
    real(dp), intent(out) :: vlow, vhigh
    real(dp) :: half, two, ulp, rtunfl, temp

    half = 0.5_dp
    two = 2.0_dp
    ulp = dlamch('Precision')
    rtunfl = sqrt(dlamch('Safe minimum'))
    if (idxlo /= 1) then
      temp = half*abs(sv(idxlo)-sv(idxlo-1))
      vhigh = sv(idxlo) + max(temp, ulp*matrix_norm, two*rtunfl)
    else
      temp = half*abs(sv(n)-sv(1))
      vhigh = sv(1) + max(temp, ulp*matrix_norm, two*rtunfl)
    end if
    if (idxhi /= n) then
      temp = half*abs(sv(idxhi+1)-sv(idxhi))
      vlow = sv(idxhi) - max(ulp*matrix_norm, two*rtunfl, temp)
    else
      temp = half*abs(sv(n)-sv(1))
      vlow = sv(n) - max(ulp*matrix_norm, two*rtunfl, temp)
    end if
    vlow = max(vlow, 0.0_dp)
    vhigh = max(vhigh, 0.0_dp)
    if (vlow >= vhigh) vhigh = max(two*vhigh, vhigh+vlow+half)
  end subroutine value_bounds

  subroutine unpack_vectors(n, nfound, packed, left, right_t)
    integer, intent(in) :: n, nfound
    real(dp), intent(in) :: packed(2*maxn,maxn)
    real(dp), intent(out) :: left(maxn,maxn), right_t(maxn,maxn)
    integer :: j, k

    do j = 1, nfound
      do k = 1, n
        left(k,j) = packed(k,j)
        right_t(j,k) = packed(n+k,j)
      end do
    end do
  end subroutine unpack_vectors

  subroutine emit_row(method, itype, n, range, nfound, residual, &
       left_orth, right_orth, orth, elapsed, ierr, nreps, backend_path)
    character(len=*), intent(in) :: method
    character(len=1), intent(in) :: range
    integer, intent(in) :: itype, n, nfound, ierr, nreps, backend_path
    real(dp), intent(in) :: residual, left_orth, right_orth, orth, elapsed

    write(*,'(A,",",I0,",",I0,",",A,",",I0,5(",",ES24.16E3),3(",",I0))') &
         trim(method), itype, n, range, nfound, residual, left_orth, &
         right_orth, orth, elapsed, ierr, nreps, backend_path
  end subroutine emit_row

end program dbdsv_benchmark
