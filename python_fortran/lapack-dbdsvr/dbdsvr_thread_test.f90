! ===========================================================
! Program: dbdsvr_thread_test
!
! Thread-safety smoke test for DBDSVR.  Runs four DBDSVR calls in
! parallel under OpenMP on four different bidiagonal matrices, then
! compares each thread's output bit-for-bit against a serial reference
! computed on the same matrix.  Any divergence is reported.
!
! Build:  gfortran -fopenmp -O2 dbdsvr_thread_test.f90 \
!                 mr3gk_fortran/*.o libxmr.so -llapack -lblas -lgfortran -lm
! Run:    OMP_NUM_THREADS=4 ./dbdsvr_thread_test
! ===========================================================
program dbdsvr_thread_test
   use mr3gk_consts, only : dp
   use omp_lib
   implicit none

   integer, parameter :: ncases = 4
   integer, parameter :: sizes(ncases) = (/ 32, 48, 64, 80 /)
   real(dp), allocatable :: d_all(:,:), e_all(:,:)
   real(dp), allocatable :: s_serial(:,:), z_serial(:,:,:)
   real(dp), allocatable :: s_par(:,:),    z_par(:,:,:)
   integer  :: nmax, n, k, i, j, ns, info, lwork
   real(dp), allocatable :: work(:)
   integer,  allocatable :: iwork(:)
   integer :: nfail, divergent_cases

   nmax = maxval(sizes)
   allocate(d_all(nmax, ncases), e_all(max(nmax-1,1), ncases))
   allocate(s_serial(nmax, ncases), z_serial(2*nmax, nmax, ncases))
   allocate(s_par(nmax, ncases),    z_par(2*nmax, nmax, ncases))

!  Reproducible inputs: glued-Wilkinson-ish saw tooth + small alternating off-diag.
   do k = 1, ncases
      n = sizes(k)
      do i = 1, n
         d_all(i,k) = 1.0_dp + 0.1_dp * real(mod(i,5), dp) + 0.01_dp * real(k, dp)
      end do
      do i = 1, n-1
         e_all(i,k) = 0.5_dp + 0.05_dp * real(mod(i,3), dp)
      end do
   end do

!  --- Serial reference ---
   do k = 1, ncases
      n = sizes(k)
      lwork = 2*n*n + 7*n
      allocate(work(lwork), iwork(n))
      call dbdsvr('U', 'V', 'A', n, d_all(:,k), e_all(:,k), 0.0_dp, 0.0_dp, &
                  0, 0, ns, s_serial(:,k), z_serial(:,:,k), 2*n, &
                  work, lwork, iwork, info)
      if (info /= 0) then
         write(*,'(A,I0,A,I0)') 'serial case ', k, ' INFO=', info
      end if
      deallocate(work, iwork)
   end do

!  --- Parallel run ---
   call omp_set_num_threads(min(ncases, max(1, omp_get_max_threads())))
   divergent_cases = 0

!$omp parallel do schedule(static) private(k, n, lwork, work, iwork, ns, info)
   do k = 1, ncases
      n = sizes(k)
      lwork = 2*n*n + 7*n
      allocate(work(lwork), iwork(n))
      call dbdsvr('U', 'V', 'A', n, d_all(:,k), e_all(:,k), 0.0_dp, 0.0_dp, &
                  0, 0, ns, s_par(:,k), z_par(:,:,k), 2*n, &
                  work, lwork, iwork, info)
      deallocate(work, iwork)
   end do
!$omp end parallel do

!  --- Bit-compare ---
   nfail = 0
   do k = 1, ncases
      n = sizes(k)
      do i = 1, n
         if (s_serial(i,k) /= s_par(i,k)) nfail = nfail + 1
      end do
      do j = 1, n
         do i = 1, 2*n
            if (z_serial(i,j,k) /= z_par(i,j,k)) nfail = nfail + 1
         end do
      end do
      if (nfail > 0) then
         write(*,'(A,I0,A,I0,A,I0)') '  case ', k, ' (n=', n, '): DIVERGENT, mismatches=', nfail
         divergent_cases = divergent_cases + 1
         nfail = 0
      else
         write(*,'(A,I0,A,I0,A)') '  case ', k, ' (n=', n, '): OK'
      end if
   end do

   if (divergent_cases == 0) then
      write(*,'(A,I0,A)') 'PASS: ', ncases, ' parallel cases bit-identical to serial'
   else
      write(*,'(A,I0,A)') 'FAIL: ', divergent_cases, ' parallel cases diverged'
      stop 1
   end if

   deallocate(d_all, e_all, s_serial, z_serial, s_par, z_par)
end program dbdsvr_thread_test
