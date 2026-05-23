! ===========================================================
! Program: mr3gk_run
! CLI driver: read bidiagonal from binary input file, run dmr3gk_svd,
! write outputs (sigma, U, V) to binary output file.
!
! Input  format: int32 n, n doubles d, n-1 doubles e
! Output format: int32 info, n doubles sigma, n*n doubles U (col-major), n*n doubles V
! ===========================================================
program mr3gk_run
   use mr3gk_consts
   use mr3gk
   implicit none

   character(len=512) :: in_path, out_path
   integer :: nargs, ios, n, i, j, info_out
   integer(kind=4) :: n32
   real(dp), allocatable :: d(:), e(:), sigma(:), U(:,:), V(:,:)
   integer, parameter :: u_in = 10, u_out = 11

   nargs = command_argument_count()
   if (nargs < 2) then
      write(*,*) 'Usage: mr3gk_run <input.bin> <output.bin>'
      stop 1
   end if
   call get_command_argument(1, in_path)
   call get_command_argument(2, out_path)

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

   allocate(sigma(n), U(n, n), V(n, n))
   call dmr3gk_svd(n, d, e, sigma, U, V, n, info_out)

   open(u_out, file=trim(out_path), form='unformatted', access='stream', &
        status='replace', iostat=ios)
   if (ios /= 0) then
      write(*,*) 'ERROR: cannot open output ', trim(out_path)
      stop 3
   end if
   write(u_out) int(info_out, kind=4)
   write(u_out) (sigma(i), i=1, n)
   ! Column-major: write column by column
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
end program mr3gk_run
