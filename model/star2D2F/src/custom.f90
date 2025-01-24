!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_EQN
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE CUSTOM_EQN

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom arrays !
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_HYDRO
USE CUSTOM_DEF
IMPLICIT NONE

! atmospheric values !
ALLOCATE (primNM_a(imin:imax))
ALLOCATE (primDM_a(imin:imax))

! gravitational potential energy !
ALLOCATE (rho_grav(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE (phi(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE (phi_old(-2:nx+3,-2:ny+3,-2:nz+3))

! for poisson equation !
ALLOCATE (ajp1(1:nx))
ALLOCATE (ajm1(1:nx))
ALLOCATE (bkp1(1:nx,1:ny))
ALLOCATE (bkm1(1:nx,1:ny))
ALLOCATE (clp1(1:nx,1:ny,1:nz))
ALLOCATE (clm1(1:nx,1:ny,1:nz))
ALLOCATE (epsc(1:nx,1:ny,1:nz))

! for gr potential !
IF (gr_potential == 1) THEN
  ALLOCATE (phi_gr(-2:nx+3))
END IF

! Lapse function !
IF (lapse_flag == 1) THEN
  ALLOCATE (a_lapse(-2:nx+3,-2:ny+3,-2:nz+3))
  ALLOCATE (apL(-2:nx+3,-2:ny+3,-2:nz+3))
  ALLOCATE (apR(-2:nx+3,-2:ny+3,-2:nz+3))
END IF

END SUBROUTINE CUSTOM_HYDRO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Populate custom arrays to GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_POPULATE
USE CUSTOM_DEF 
IMPLICIT NONE

! Now populate all necessary, and reuseable arrays to the graphic cards !
!$ACC enter DATA COPYIN(phi, phi_old, ajp1, ajm1, bkp1, bkm1, clp1, clm1, epsc)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Clear custom arrays from GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CLEAR
USE CUSTOM_DEF 
IMPLICIT NONE

! Now we clear memory in the GPU device !
!$ACC exit DATA DELETE(phi, phi_old, ajp1, ajm1, bkp1, bkm1, clp1, clm1, epsc)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_GRID
USE CUSTOM_DEF
IMPLICIT NONE

! Integer !
INTEGER :: i, j, nlines

! Read the number of lines in the file !
nlines = 0 
OPEN (970, file = './profile/hydro_x1_fgrid.dat') 
DO 
  READ (970,*, END=10) 
  nlines = nlines + 1 
END DO 
10 CLOSE (970) 

! Error message !
IF(nlines .ne. nx+7) THEN
  WRITE (*,*) 'number of grid faces from files', nlines
  WRITE (*,*) 'number of x grid faces in the program', nx+6
  STOP 'inconsistent number of grid faces, exit'
END IF

! Read !
OPEN(UNIT=970, FILE = './profile/hydro_x1_fgrid.dat', ACTION='READ')
DO i = -3, nx+3
	READ(970,*) xF(i)
ENDDO
CLOSE(970)
! xF(0) = xF(0)*10.0D0
xF = xF*lencgs2code

! nlines = 0 
! OPEN (970, file = './profile/hydro_x2_fgrid.dat')
! DO 
!   READ (970,*, END=11) 
!   nlines = nlines + 1 
! END DO 
! 11 CLOSE (970) 

! ! Error message !
! IF(nlines .ne. ny+7) THEN
!   WRITE (*,*) 'number of grid faces from files', nlines
!   WRITE (*,*) 'number of y grid faces in the program', ny+6
!   STOP 'inconsistent number of grid faces, exit'
! END IF

! ! Read !
! OPEN(UNIT=970, FILE = './profile/hydro_x2_fgrid.dat', ACTION='READ')
! DO i = -3, ny+3
! 	READ(970,*) yF(i)
! ENDDO
! CLOSE(970)

END SUBROUTINE CUSTOM_GRID

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY_X
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE CUSTOM_BOUNDARY_X

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY_Y
USE DEFINITION

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Back up fluxes from riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLUX_X
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Back up fluxes from riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLUX_Y
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Back up fluxes from riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLUX_Z
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom variable floor !
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CHECKRHO
USE CUSTOM_DEF
IMPLICIT NONE

INTEGER :: j, k, l

! Threshold for atmosphere density
REAL*8 :: factor, diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP PARALLEL DO PRIVATE(diff, factor) COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(diff, factor)
DO l = iNM, iDM
  DO k = 1, ny
    DO j = 1, nx
      ! Standard !
      IF (prim(irho,j,k,l) <= rhoatm) THEN
        prim(irho,j,k,l) = MAX(prim(irho,j,k,l), rhomin)
        prim(ivx:ivz,j,k,l) = 0.0D0
        epsilon(j,k,l) = MAX(epsilon(j,k,l), 0.0D0)
      END IF
    END DO
  END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE CUSTOM_CHECKRHO

!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom variable floor !
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOMFLOOR
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom build states   !
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BUILDSTATE(dir_in)
USE CUSTOM_DEF
USE PPMC_MODULE
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER, INTENT(INOUT) :: dir_in
INTEGER :: i, j, k, l

IF (lapse_flag == 1) THEN
  SELECT CASE (dir_in)
  CASE (x_dir)
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
    DO l = iNM, iDM
      DO k = 1, ny
        DO j = 0, nx + 1 
          CALL PPMC (wx(j,1:14), a_lapse(j-2,k,l), a_lapse(j-1,k,l), a_lapse(j,k,l), a_lapse(j+1,k,l), a_lapse(j+2,k,l), apR(j-1,k,l), apL(j,k,l))
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
  CASE (y_dir)
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
    DO l = iNM, iDM
      DO k = 0, ny + 1 
        DO j = 1, nx
          CALL PPMC (wy(k,1:14), a_lapse(j,k-2,l), a_lapse(j,k-1,l), a_lapse(j,k,l), a_lapse(j,k+1,l), a_lapse(j,k+2,l), apR(j,k-1,l), apL(j,k,l))
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
  END SELECT
  
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
  DO l = iNM, iDM
    DO k = 0, ny + 1
      DO j = 0, nx + 1
        fluxL(:,j,k,l) = fluxL(:,j,k,l)*apL(j,k,l)
        fluxR(:,j,k,l) = fluxR(:,j,k,l)*apR(j,k,l)
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO
END IF
END SUBROUTINE CUSTOM_BUILDSTATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_SOURCE
USE CUSTOM_DEF 
IMPLICIT NONE

INTEGER :: j, k, l
REAL*8 :: dphi_grdx, dphidx, dphidy
REAL*8 :: factor, diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(dphi_grdx, dphidx, dphidy, diff, factor)
DO l = iNM, iDM
  DO k = 1, ny
    DO j = 1, nx
      ! Include only non-atmosphere !
			diff = prim(irho,j,k,l) - primNM_a(irho)
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)
      ! Gravitational potential of the matter !
      dphidx = first_derivative (x(j-1), x(j), x(j+1), phi(j-1,k,l), phi(j,k,l), phi(j+1,k,l))
      dphidy = first_derivative (y(k-1), y(k), y(k+1), phi(j,k-1,l), phi(j,k,l), phi(j,k+1,l))

      IF (gr_potential == 1) THEN
        dphi_grdx = first_derivative (x(j-1), x(j), x(j+1), phi_gr(j-1), phi_gr(j), phi_gr(j+1))
        dphidx = dphidx + dphi_grdx
      END IF

      IF (lapse_flag == 1) then
        ! Geometric source for momentum equation!
        sc(ivx:ivz,j,k,l) = sc(ivx:ivz,j,k,l)*a_lapse(j,k,l)

        ! Geometric and gravitational source terms !
        sc(ivx,j,k,l) = sc(ivx,j,k,l) - a_lapse(j,k,l)*factor*(prim(irho,j,k,l) - prim(itau,j,k,l))*dphidx
        sc(ivy,j,k,l) = sc(ivy,j,k,l) - a_lapse(j,k,l)*factor*(prim(irho,j,k,l) - prim(itau,j,k,l))*dphidy/x(j)
        sc(itau,j,k,l) = sc(itau,j,k,l) - a_lapse(j,k,l)*factor*prim(irho,j,k,l)*(prim(ivx,j,k,l)*dphidx + prim(ivy,j,k,l)*dphidy/x(j))
      ELSE
      ! Add them to the source term !
      ! sc(ivx,j,k,l) = sc(ivx,j,k,l) - factor*prim(irho,j,k,l)*dphidx
      ! sc(ivy,j,k,l) = sc(ivy,j,k,l) - factor*prim(irho,j,k,l)*dphidy/x(j)
      ! sc(itau,j,k,l) = sc(itau,j,k,l) - factor*prim(irho,j,k,l)*(prim(ivx,j,k,l)*dphidx + prim(ivy,j,k,l)*dphidy/x(j))
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

	REAL*8 function first_derivative (xm1, xc, xp1, fm1, fc, fp1)
	!$acc routine seq
	implicit none
	REAL*8 :: xm1, xc, xp1, fm1, fc, fp1, h1, h2
  h2 = xp1 - xc
  h1 = xc - xm1
	first_derivative = ((fp1-fc)*h1*h1+(fc-fm1)*h2*h2)/(h1*h2*(h1+h2))
	end function

END SUBROUTINE CUSTOM_SOURCE

!!!!!!!!!!!!!!!!!!!!!
! Do custom updates !
!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_UPDATE (p_in)
USE CUSTOM_DEF 
IMPLICIT NONE

INTEGER, INTENT (IN) :: p_in
INTEGER :: j, k, l

! Update gravitational potentials !
IF (MOD(n_step, n_pot) == 0) THEN
  CALL POISSON_SOLVER

  IF (lapse_flag == 1) THEN
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
    DO j = 1, nx
      DO k = 1, ny
        DO l = iNM, iDM
          a_lapse(j,k,l) = exp(phi(j,k,l) + phi_gr(j))
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO

  CALL BOUNDARY1D_NM (a_lapse, even, even, even, even, even, even)
  END IF
END IF

END SUBROUTINE CUSTOM_UPDATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPERATOR_SPLIT
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE OPERATOR_SPLIT

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENFILE_CUSTOM
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE OPENFILE_CUSTOM

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_ANALYSIS
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE CUSTOM_ANALYSIS
