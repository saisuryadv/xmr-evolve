! ===========================================================
! Program: mr3gk_run_traced
!
! Instrumented variant of mr3gk_run that opens FORTRAN unit 99 to
! the trace-log path passed as argv[3] BEFORE invoking dmr3gk_svd.
!
! This unit is what dlaxrf_traced.f writes one log line per tree-node
! visit to (NODE depth=NN n=NN icbeg=NN icend=NN tau=...
!           pivot_dmin=... pivot_dmax=... eig_lo=... eig_hi=...
!           ncd_eigs=... ncd_pivots=...).
!
! Usage: mr3gk_run_traced <input.bin> <output.bin> <trace.log>
!
! Input/output binary formats are identical to mr3gk_run.
! ===========================================================
program mr3gk_run_traced
   use mr3gk_consts
   use mr3gk
   implicit none

   character(len=512) :: in_path, out_path, trace_path
   integer :: nargs, ios, n, i, j, info_out
   integer(kind=4) :: n32
   real(dp), allocatable :: d(:), e(:), sigma(:), U(:,:), V(:,:)
   integer, parameter :: u_in = 10, u_out = 11, u_trace = 99

   nargs = command_argument_count()
   if (nargs < 3) then
      write(*,*) 'Usage: mr3gk_run_traced <input.bin> <output.bin> <trace.log>'
      stop 1
   end if
   call get_command_argument(1, in_path)
   call get_command_argument(2, out_path)
   call get_command_argument(3, trace_path)

   open(u_in, file=trim(in_path), form='unformatted', access='stream', &
        status='old', iostat=ios)
   if (ios /= 0) then
      write(*,*) 'ERROR: cannot open ', trim(in_path)
      stop 2
   end if
   read(u_in) n32
   n = n32
   allocate(d(0:n-1))
   if (n > 1) allocate(e(0:n-2))
   if (n == 1) allocate(e(0:0))
   read(u_in) (d(i), i = 0, n-1)
   if (n > 1) then
      read(u_in) (e(i), i = 0, n-2)
   end if
   close(u_in)

   ! Open the trace-log unit BEFORE calling the SVD so that
   ! every WRITE(99, ...) inside dlaxrf_traced.f is captured.
   open(u_trace, file=trim(trace_path), form='formatted', &
        status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      write(*,*) 'ERROR: cannot open trace ', trim(trace_path)
      stop 4
   end if
   write(u_trace, '(A,I0)') 'TRACE_BEGIN n=', n

   allocate(sigma(n), U(n, n), V(n, n))
   call dmr3gk_svd(n, d, e, sigma, U, V, n, info_out)

   write(u_trace, '(A,I0)') 'TRACE_END info=', info_out
   close(u_trace)

   open(u_out, file=trim(out_path), form='unformatted', access='stream', &
        status='replace', iostat=ios)
   if (ios /= 0) then
      write(*,*) 'ERROR: cannot open output ', trim(out_path)
      stop 3
   end if
   write(u_out) int(info_out, kind=4)
   write(u_out) (sigma(i), i=1, n)
   do j = 1, n
      do i = 1, n
         write(u_out) U(i, j)
      end do
   end do
   do j = 1, n
      do i = 1, n
         write(u_out) V(i, j)
      end do
   end do
   close(u_out)

   deallocate(d, e, sigma, U, V)
end program mr3gk_run_traced
