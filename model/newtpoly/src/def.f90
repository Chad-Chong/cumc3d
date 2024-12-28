!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER files for the a newtonian polytrope star       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE PARAMETER
USE DEFINITION
IMPLICIT NONE
SAVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Floor !
REAL*8, ALLOCATABLE, DIMENSION (:) :: prim_a
REAL*8 :: eps_a
INTEGER :: jNS
! Polytope !
REAL*8, PARAMETER :: npoly = 1.0D0
REAL*8, PARAMETER :: kpoly = 100.0D0
REAL*8, PARAMETER :: rhocen = 1.28D-3
REAL*8, PARAMETER :: alpha = SQRT((npoly + 1.0D0)*kpoly/4.0D0/pi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gravity !
REAL*8, ALLOCATABLE :: dphidx(:,:,:)

! Floors !
REAL*8 :: rho_floor
REAL*8 :: eps_floor

END MODULE
