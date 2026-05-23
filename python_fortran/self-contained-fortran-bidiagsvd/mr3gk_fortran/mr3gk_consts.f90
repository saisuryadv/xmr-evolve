! ===========================================================
! Module: mr3gk_consts
! Constants for MR3-GK SVD (mirrors mr3_gk.py top-level constants).
! ===========================================================
module mr3gk_consts
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   real(dp), parameter :: EPS    = 2.2204460492503131d-16    ! np.finfo(np.float64).eps
   real(dp), parameter :: SAFMIN = 2.2250738585072014d-308   ! np.finfo(np.float64).tiny
   real(dp), parameter :: GAPTOL = 1.0d-3
   integer,  parameter :: MAX_DEPTH  = 50
   integer,  parameter :: ELG_THRESH = 8
   integer,  parameter :: MAXRELCOND = 10
end module mr3gk_consts
