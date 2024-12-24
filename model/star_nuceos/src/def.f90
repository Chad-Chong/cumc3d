!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Custom modular files to contain all custom arrays 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE CUSTOM_DEF
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! for electron fraction !
INTEGER :: iye

! Floors !
REAL*8, ALLOCATABLE, DIMENSION (:) :: prim_a

! gravitational potential !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: phi, phi_old

! for poisson solver !
REAL*8, ALLOCATABLE, DIMENSION (:) :: ajp1, ajm1
REAL*8, ALLOCATABLE, DIMENSION (:,:) :: bkp1, bkm1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: clp1, clm1, epsc

! temperature for eos !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: temperature
REAL*8, PARAMETER :: temp_a = 1.0D-2
CHARACTER(LEN=100) :: eos_file

! for poisson solver !
REAL*8, PARAMETER :: omega_weight = 1.9d0

! epsilon at atmosphere !
REAL*8 :: eps_a

END MODULE
