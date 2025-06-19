!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Custom modular files to contain all custom arrays 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE CUSTOM_DEF
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Floors !
REAL*8, ALLOCATABLE, DIMENSION (:) :: prim_a

! gravitational potential !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: phi, phi_old

! for poisson solver !
REAL*8, ALLOCATABLE, DIMENSION (:) :: ajp1, ajm1
REAL*8, ALLOCATABLE, DIMENSION (:,:) :: bkp1, bkm1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: clp1, clm1, epsc

! for poisson solver !
REAL*8, PARAMETER :: omega_weight = 1.9d0

! temperature (HelmEOS) and epsilon at atmosphere !
REAL*8 :: temp2_a = 1.0D-1
REAL*8 :: eps_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For Chemical Composition (HelmEOS) !

! Number of isotope
INTEGER, PARAMETER :: totalion = 7 !19 !204

! The last global time for burning
REAL*8 :: last_burntime

! Density limit for nuclear burning
real*8, PARAMETER :: rho2_burn_max = 1.62D-11
real*8, PARAMETER :: rho2_burn_min = 8.10D-12

! Density limit for deflagration
real*8, PARAMETER :: rho2_flame_max = 1.62D-8
real*8, PARAMETER :: rho2_flame_min = 1.62D-11

! Density limit for detonation
real*8, PARAMETER :: rho2_deton_max = 1.62D-9
REAL*8, PARAMETER :: rho2_deton_min = 1.62D-13

! Some global quantities
! 1. Central temperature
! 2. Total thermal neutrino loss
! 3. Total non-thermal neutrino loss
! 4. Luminsoity
! 5. Luminosity by burning, def. and det.

REAL*8 :: centraltemperature
REAL*8 :: total_nu_qdot
REAL*8 :: total_ecap_nu_qdot
REAL*8 :: lumino
REAL*8 :: lumino_burn, lumino_flame, lumino_deton

! MAss burned by deflgration/detonation
REAL*8 :: burn_mass

! Atmosphere chemical composition
REAL*8 :: abar2_a
REAL*8 :: zbar2_a
REAL*8 :: ye2_a
REAL*8 :: yiso_a, qash_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! My iso7 prescription

integer :: ihe4, ic12, io16, ine20, img24, isi28, ini56
integer :: che4, cc12, co16, cne20, cmg24, csi28, cni56
character(len=5) :: ionam(1:totalion)
REAL*8 :: mev2erg,mev2gr
REAL*8, PARAMETER :: ev2erg = 1.60217648740d-12
REAL*8, PARAMETER :: avo     = 6.0221367d23
PARAMETER        (mev2erg = ev2erg*1.0d6, &
                    mev2gr  = mev2erg/clight**2)


real (selected_real_kind(15,307)), dimension(1:totalion) :: aion, zion, bion, nion, mion, wion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section of NSE table 1
! Table size
INTEGER, PARAMETER :: temp_rowno_nse = 70
INTEGER, PARAMETER :: den_rowno_nse = 30

! Binding energy and composition of NSE
REAL*8, DIMENSION(0:den_rowno_nse, 0:temp_rowno_nse):: nsetable_binde
REAL*8, DIMENSION(0:den_rowno_nse, 0:temp_rowno_nse, 1:totalion):: nsetable_xiso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section of NSE table 2
! Table size         
INTEGER, PARAMETER :: ent_rowno_nse2 = 50
INTEGER, PARAMETER :: ye_rowno_nse2 = 26
INTEGER, PARAMETER :: den_rowno_nse2 = 50

! New binding energy and composition of NSE with a larger network
REAL*8, DIMENSION(0:den_rowno_nse2, 0:ye_rowno_nse2+1, 3):: nsetable2_head
REAL*8, DIMENSION(0:den_rowno_nse2, 0:ye_rowno_nse2+1, 0:ent_rowno_nse2+1, 6):: nsetable2_binde

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Section of NSE table 3                       
! Table size
INTEGER, PARAMETER :: temp_rowno_nse3 = 48     
INTEGER, PARAMETER :: ye_rowno_nse3 = 122      
INTEGER, PARAMETER :: den_rowno_nse3 = 23

! New binding energy and composition of NSE with a larger network
REAL*8, DIMENSION(0:den_rowno_nse3+1, 0:ye_rowno_nse3+1, 0:temp_rowno_nse3+1, 6):: nsetable3_binde

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocatable objects

! Flag for being in NSE state
! nse_flag = 0 means in C-burning 
! nse_flag = 1 means in O- and Mg- burning
! nse_flag = 2 means in NSE burning

INTEGER  , allocatable, DIMENSION (:,:,:) :: nse_flag

! Energy loss by neutrino
REAL*8, allocatable, DIMENSION (:,:,:) :: nu_qdot

! Energy input by burning
real*8, allocatable, DIMENSION (:,:,:) :: burn_qdot

! Energy input by deflagration
real*8, allocatable, DIMENSION (:,:,:) :: flame_qdot

! Energy input by detonation
real*8, allocatable, DIMENSION (:,:,:) :: deton_qdot

! Mean atomic mass
real*8, allocatable, DIMENSION (:,:,:) :: abar2

! Mean atomic number
real*8, allocatable, DIMENSION (:,:,:) :: zbar2

! Temperature
REAL*8, allocatable, DIMENSION (:,:,:) :: temp2
REAL*8, allocatable, DIMENSION (:,:,:) :: temp2_old

! Dummy variables for time-evolution
REAL*8, allocatable, DIMENSION (:,:,:,:) :: xiso1, delta_xiso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For Flame !

! New set of table
! Size of deflagration table
INTEGER, PARAMETER :: temp_rowno1 = 0
INTEGER, PARAMETER :: den_rowno1 = 30

! Deflagration quantities
REAL*8, DIMENSION(0:den_rowno1, 0:temp_rowno1):: ashtable_eps
REAL*8, DIMENSION(0:den_rowno1, 0:temp_rowno1, 9):: ashtable_state
REAL*8, DIMENSION(0:den_rowno1, 0:temp_rowno1, 1:totalion):: ashtable_xiso

! Size of detonation table
INTEGER, PARAMETER :: temp_rowno2 = 0
INTEGER, PARAMETER :: den_rowno2 = 30

! Detonation quantities
REAL*8, DIMENSION(0:den_rowno2, 0:temp_rowno2):: dettable_eps, dettable_vel
REAL*8, DIMENSION(0:den_rowno2, 0:temp_rowno2, 9):: dettable_state
REAL*8, DIMENSION(0:den_rowno2, 0:temp_rowno2, 1:totalion):: dettable_xiso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For Electron Capture !

INTEGER :: iye2

! Refer to the module itself; this is to avoid redeclaration !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! For Turbulence !

INTEGER :: iturbq

REAL*8 :: turb_q_a = 1.0D-10
REAL*8 :: turb_qtotal
REAL*8 :: turb_q_max

! Output arrays
REAL*8, allocatable, dimension(:,:,:) :: turb_source
REAL*8, allocatable, dimension(:,:,:) :: turb_diff

! The k-eps component
REAL*8, allocatable, dimension(:,:,:) :: turb_eps

END MODULE