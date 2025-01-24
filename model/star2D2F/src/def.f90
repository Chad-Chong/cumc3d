!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom modular files to contain all custom arrays 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE CUSTOM_DEF
USE DEFINITION
IMPLICIT NONE
SAVE
INCLUDE "param.h"

! Floors !
REAL*8, ALLOCATABLE, DIMENSION (:) :: primNM_a
REAL*8 :: epsNM_a
REAL*8, ALLOCATABLE, DIMENSION (:) :: primDM_a
REAL*8 :: epsDM_a


! gravitational potential !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: rho_grav
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: phi, phi_old
REAL*8, ALLOCATABLE, DIMENSION (:) :: phi_gr

! for poisson solver !
REAL*8, ALLOCATABLE, DIMENSION (:) :: ajp1, ajm1
REAL*8, ALLOCATABLE, DIMENSION (:,:) :: bkp1, bkm1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: clp1, clm1, epsc
REAL*8, PARAMETER :: omega_weight = 1.9d0

! Lapse function !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: a_lapse, apL, apR

END MODULE
