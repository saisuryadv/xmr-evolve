program dbdsvr_case_metrics
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use, intrinsic :: iso_fortran_env, only: int64, real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: dp = real64
  integer :: i, idx, info, ios, liwork, lwork, n, ns
  integer(int64) :: count0, count1, count_rate
  character(len=4096) :: input_path
  character(len=32) :: route, orth_mode
  real(dp) :: email_orth_u, email_orth_v, email_residual
  real(dp) :: file_d, file_e, metric_seconds, orth_u, orth_v
  real(dp) :: residual, solver_seconds
  real(dp), allocatable :: d(:), e(:), s(:), z(:,:), work(:)
  integer, allocatable :: iwork(:)
  integer(int64) :: orth_pairs

  interface
    subroutine cblas_dgemm(layout, trans_a, trans_b, m, n, k, alpha, &
         a, lda, b, ldb, beta, c, ldc) bind(C, name='cblas_dgemm')
      import :: c_double, c_int
      integer(c_int), value :: layout, trans_a, trans_b, m, n, k
      integer(c_int), value :: lda, ldb, ldc
      real(c_double), value :: alpha, beta
      real(c_double), intent(in) :: a(*), b(*)
      real(c_double), intent(inout) :: c(*)
    end subroutine cblas_dgemm
  end interface

  call get_command_argument(1, input_path)
  if (len_trim(input_path) == 0) then
    write(*,'(A)') 'usage: dbdsvr_case_metrics MATRIX.dat'
    stop 2
  end if
  open(unit=10, file=trim(input_path), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(*,'(A,A)') 'cannot open ', trim(input_path)
    stop 2
  end if
  read(10,*,iostat=ios) n
  if (ios /= 0 .or. n < 2) then
    write(*,'(A,I0)') 'bad matrix order: ', n
    stop 2
  end if

  lwork = max(1, 5*n*n + 37*n)
  liwork = max(1, 20*n)
  allocate(d(n), e(n), s(n), z(2*n,n), work(lwork), iwork(liwork))
  iwork = 0
  do i = 1, n
    read(10,*,iostat=ios) idx, file_d, file_e
    if (ios /= 0) then
      write(*,'(A,I0)') 'bad matrix row: ', i
      stop 2
    end if
    d(i) = file_d
    e(i) = file_e
  end do
  close(10)
  e(n) = 0.0_dp

  call system_clock(count0, count_rate)
  call dbdsvr('U', 'V', 'A', n, d, e, 0.0_dp, 0.0_dp, 0, 0, ns, &
       s, z, 2*n, work, lwork, iwork, liwork, info)
  call system_clock(count1)
  solver_seconds = real(count1-count0, dp) / real(count_rate, dp)

  residual = huge(1.0_dp)
  orth_u = huge(1.0_dp)
  orth_v = huge(1.0_dp)
  email_residual = huge(1.0_dp)
  email_orth_u = huge(1.0_dp)
  email_orth_v = huge(1.0_dp)
  orth_pairs = 0_int64
  orth_mode = 'not-run'
  metric_seconds = 0.0_dp
  if (info == 0 .and. finite_solution(n, s, z)) then
    call system_clock(count0)
    call triplet_residual(n, d, e, s, z, residual)
    call orthogonality_metric(n, z, orth_u, orth_v, orth_mode, orth_pairs)
    call email_metrics(n, d, e, s, z, work, email_residual, &
         email_orth_u, email_orth_v)
    call system_clock(count1)
    metric_seconds = real(count1-count0, dp) / real(count_rate, dp)
  else if (info == 0) then
    orth_mode = 'nonfinite-output'
  end if

  if (info /= 0) then
    route = 'all-attempts-failed'
  else
    select case (iwork(2))
    case (1)
      route = 'mr3'
    case (2)
      route = 'dbdsdc-after-mr3-info'
    case (3)
      route = 'dbdsdc-after-mr3-output'
    case default
      route = 'unknown-success-path'
    end select
  end if
  write(*,'(A)') 'n,info,route,residual,orth_u,orth_v,email_residual,' // &
       'email_orth_u,email_orth_v,orth_mode,orth_pairs,solver_seconds,' // &
       'metric_seconds'
  write(*,'(I0,",",I0,",",A,6(",",ES24.16E3),",",A,",",I0,2(",",ES24.16E3))') &
       n, info, trim(route), residual, orth_u, orth_v, email_residual, &
       email_orth_u, email_orth_v, trim(orth_mode), orth_pairs, &
       solver_seconds, metric_seconds

contains

  logical function finite_solution(order, sigma, packed)
    integer, intent(in) :: order
    real(dp), intent(in) :: sigma(order), packed(2*order,order)
    integer :: col, row

    finite_solution = .false.
    do col = 1, order
      if (.not. ieee_is_finite(sigma(col))) return
      do row = 1, 2*order
        if (.not. ieee_is_finite(packed(row,col))) return
      end do
    end do
    finite_solution = .true.
  end function finite_solution

  subroutine triplet_residual(order, bd, be, sigma, packed, result)
    integer, intent(in) :: order
    real(dp), intent(in) :: bd(order), be(order), sigma(order)
    real(dp), intent(in) :: packed(2*order,order)
    real(dp), intent(out) :: result
    integer :: col, row
    real(dp) :: bnorm, eps, rnorm, scale, value
    real(dp), allocatable :: rleft(:), rright(:)
    real(dp), external :: dlamch, dnrm2

    allocate(rleft(order), rright(order))
    bnorm = max(abs(sigma(1)), tiny(1.0_dp))
    eps = dlamch('Precision')
    rnorm = 0.0_dp
    do col = 1, order
      do row = 1, order
        value = bd(row)*packed(order+row,col) - &
             sigma(col)*packed(row,col)
        if (row < order) value = value + be(row)*packed(order+row+1,col)
        rleft(row) = value

        value = bd(row)*packed(row,col) - &
             sigma(col)*packed(order+row,col)
        rright(row) = value
      end do
      do row = 2, order
        rright(row) = rright(row) + be(row-1)*packed(row-1,col)
      end do
      if (any(.not. ieee_is_finite(rleft)) .or. &
          any(.not. ieee_is_finite(rright))) then
        result = huge(1.0_dp)
        deallocate(rleft, rright)
        return
      end if
      rnorm = max(rnorm, dnrm2(order, rleft, 1), dnrm2(order, rright, 1))
    end do
    scale = max(real(order,dp)*eps*bnorm, tiny(1.0_dp))
    result = rnorm / scale
    deallocate(rleft, rright)
  end subroutine triplet_residual

  subroutine email_metrics(order, bd, be, sigma, packed, scratch, &
       result_residual, result_u, result_v)
    ! Reproduce the max-norm metric contract used by the 2026-06-12
    ! email sweep.  BLAS performs the same U*S*V^T and Gram products as
    ! the original triply nested loops, while keeping the full practical
    ! suite tractable.  WORK is solver scratch after DBDSVR returns; its
    ! documented size leaves at least two order-by-order matrices here.
    integer, intent(in) :: order
    real(dp), intent(in) :: bd(order), be(order), sigma(order)
    real(dp), intent(in) :: packed(2*order,order)
    real(dp), intent(inout) :: scratch(*)
    real(dp), intent(out) :: result_residual, result_u, result_v
    integer :: col, index, offset, row
    real(dp) :: bnorm, eps, expected, scale, value
    real(dp), external :: dlamch

    offset = order*order
    do col = 1, order
      do row = 1, order
        index = row + (col-1)*order
        scratch(index) = packed(row,col)*sigma(col)
      end do
    end do
    call cblas_dgemm(102_c_int, 111_c_int, 112_c_int, order, order, &
         order, 1.0_dp, scratch(1), order, packed(order+1,1), &
         2*order, 0.0_dp, scratch(offset+1), order)

    bnorm = max(maxval(abs(sigma)), tiny(1.0_dp))
    eps = dlamch('Precision')
    result_residual = 0.0_dp
    do col = 1, order
      do row = 1, order
        expected = 0.0_dp
        if (row == col) expected = bd(row)
        if (col == row+1) expected = be(row)
        value = scratch(offset + row + (col-1)*order)
        if (.not. ieee_is_finite(value)) then
          result_residual = huge(1.0_dp)
          result_u = huge(1.0_dp)
          result_v = huge(1.0_dp)
          return
        end if
        result_residual = max(result_residual, abs(value-expected))
      end do
    end do
    result_residual = result_residual / &
         max(real(order,dp)*eps*bnorm, tiny(1.0_dp))

    call cblas_dgemm(102_c_int, 112_c_int, 111_c_int, order, order, &
         order, 1.0_dp, packed(1,1), 2*order, packed(1,1), 2*order, &
         0.0_dp, scratch(1), order)
    result_u = 0.0_dp
    do col = 1, order
      do row = 1, order
        value = scratch(row + (col-1)*order)
        if (row == col) value = value - 1.0_dp
        if (.not. ieee_is_finite(value)) then
          result_u = huge(1.0_dp)
          result_v = huge(1.0_dp)
          return
        end if
        result_u = max(result_u, abs(value))
      end do
    end do

    call cblas_dgemm(102_c_int, 112_c_int, 111_c_int, order, order, &
         order, 1.0_dp, packed(order+1,1), 2*order, &
         packed(order+1,1), 2*order, 0.0_dp, scratch(1), order)
    result_v = 0.0_dp
    do col = 1, order
      do row = 1, order
        value = scratch(row + (col-1)*order)
        if (row == col) value = value - 1.0_dp
        if (.not. ieee_is_finite(value)) then
          result_v = huge(1.0_dp)
          return
        end if
        result_v = max(result_v, abs(value))
      end do
    end do
    scale = max(real(order,dp)*eps, tiny(1.0_dp))
    result_u = result_u / scale
    result_v = result_v / scale
  end subroutine email_metrics

  subroutine orthogonality_metric(order, packed, result_u, result_v, &
       mode, pair_count)
    integer, intent(in) :: order
    real(dp), intent(in) :: packed(2*order,order)
    real(dp), intent(out) :: result_u, result_v
    character(len=*), intent(out) :: mode
    integer(int64), intent(out) :: pair_count
    integer :: available, col, i1, i2, offset, sample, samples
    real(dp) :: dot_u, dot_v, eps, target, worst_u, worst_v
    real(dp), external :: dlamch, ddot

    eps = dlamch('Precision')
    worst_u = 0.0_dp
    worst_v = 0.0_dp
    pair_count = 0_int64
    if (order <= 512) then
      mode = 'exact-all-pairs'
      do i2 = 1, order
        do i1 = 1, i2
          target = merge(1.0_dp, 0.0_dp, i1 == i2)
          dot_u = ddot(order, packed(1,i1), 1, packed(1,i2), 1)
          dot_v = ddot(order, packed(order+1,i1), 1, &
               packed(order+1,i2), 1)
          if (.not. ieee_is_finite(dot_u) .or. &
              .not. ieee_is_finite(dot_v)) then
            result_u = huge(1.0_dp)
            result_v = huge(1.0_dp)
            return
          end if
          worst_u = max(worst_u, abs(target-dot_u))
          worst_v = max(worst_v, abs(target-dot_v))
          pair_count = pair_count + 1_int64
        end do
      end do
    else
      mode = 'sampled-norm-adjacent-powers'
      do col = 1, order
        dot_u = ddot(order, packed(1,col), 1, packed(1,col), 1)
        dot_v = ddot(order, packed(order+1,col), 1, &
             packed(order+1,col), 1)
        if (.not. ieee_is_finite(dot_u) .or. &
            .not. ieee_is_finite(dot_v)) then
          result_u = huge(1.0_dp)
          result_v = huge(1.0_dp)
          return
        end if
        worst_u = max(worst_u, abs(1.0_dp-dot_u))
        worst_v = max(worst_v, abs(1.0_dp-dot_v))
        pair_count = pair_count + 1_int64
      end do
      do col = 1, order-1
        dot_u = ddot(order, packed(1,col), 1, packed(1,col+1), 1)
        dot_v = ddot(order, packed(order+1,col), 1, &
             packed(order+1,col+1), 1)
        if (.not. ieee_is_finite(dot_u) .or. &
            .not. ieee_is_finite(dot_v)) then
          result_u = huge(1.0_dp)
          result_v = huge(1.0_dp)
          return
        end if
        worst_u = max(worst_u, abs(dot_u))
        worst_v = max(worst_v, abs(dot_v))
        pair_count = pair_count + 1_int64
      end do
      offset = 2
      do while (offset < order)
        available = order - offset
        samples = min(256, available)
        do sample = 1, samples
          i1 = 1 + ((sample-1)*available)/samples
          i2 = i1 + offset
          dot_u = ddot(order, packed(1,i1), 1, packed(1,i2), 1)
          dot_v = ddot(order, packed(order+1,i1), 1, &
               packed(order+1,i2), 1)
          if (.not. ieee_is_finite(dot_u) .or. &
              .not. ieee_is_finite(dot_v)) then
            result_u = huge(1.0_dp)
            result_v = huge(1.0_dp)
            return
          end if
          worst_u = max(worst_u, abs(dot_u))
          worst_v = max(worst_v, abs(dot_v))
          pair_count = pair_count + 1_int64
        end do
        offset = 2*offset
      end do
    end if
    result_u = worst_u / (real(order,dp)*eps)
    result_v = worst_v / (real(order,dp)*eps)
  end subroutine orthogonality_metric

end program dbdsvr_case_metrics
