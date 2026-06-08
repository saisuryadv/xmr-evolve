! ===========================================================
! Program: mr3gk_run
! CLI driver: read bidiagonal from binary input file, run DBDSVR
! (LAPACK-style: DLASQ1 + Willems MR^3-GK), write outputs
! (sigma, U, V) to binary output file.
!
! Input  format: int32 n, n doubles d, n-1 doubles e
! Output format: int32 info, n doubles sigma, n*n doubles U (col-major),
!                n*n doubles V (col-major)
! ===========================================================
program mr3gk_run
   use mr3gk_consts
   implicit none

   character(len=512) :: in_path, out_path
   integer :: nargs, ios, n, i, j, info_out, ns, lwork
   integer(kind=4) :: n32
   real(dp), allocatable :: d(:), e(:), sigma(:), U(:,:), V(:,:)
   real(dp), allocatable :: z(:,:), work(:)
   integer, allocatable  :: iwork(:)
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
   allocate(d(n))
   allocate(e(max(n-1,1)))
   read(u_in) (d(i), i = 1, n)
   if (n > 1) then
      read(u_in) (e(i), i = 1, n-1)
   end if
   close(u_in)

   allocate(sigma(n), U(n, n), V(n, n))
   allocate(z(max(1,2*n), max(1,n)))
   lwork = max(1, 2*n*n + 7*n)
   allocate(work(lwork), iwork(max(1,n)))

   call dbdsvr('U', 'V', 'A', n, d, e, 0.0d0, 0.0d0, 0, 0, &
               ns, sigma, z, max(1,2*n), work, lwork, iwork, info_out)

   ! Unpack Z = [U ; V] (column-major) into U and V.
   if (n >= 1) then
      do j = 1, n
         do i = 1, n
            U(i, j) = z(i,     j)
            V(i, j) = z(n + i, j)
         end do
      end do
   end if

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

   deallocate(d, e, sigma, U, V, z, work, iwork)
end program mr3gk_run
