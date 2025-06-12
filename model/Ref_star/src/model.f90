!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model                                                           !
! Simulation for rotating magnetised white dwarf AIC                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE
INCLUDE "param.h"

! Integer !
INTEGER :: i, j, k, l

! Magnetic field !
REAL*8 :: maxdb
REAL*8 :: div_b

REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: a_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preperation !

! Allocate
Allocate(a_phi(-3:nx+3,-3:ny+3,-3:nz+3))

! Poisson interpolation coefficient !
call get_poisson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read and assign density !
OPEN(UNIT=970, FILE = './profile/hydro_rho.dat', ACTION='READ')
READ(970,*) ((prim(irho,j,k,1), j = 1, nx), k = 1, ny)
CLOSE(970)

prim(irho,:,:,:) = prim(irho,:,:,:)*rhocgs2code
PRINT *, "Finished reading rho"

! Assign velocity !
prim(ivx,:,:,1) = 0.0d0
prim(ivy,:,:,1) = 0.0d0

! Read and assign velocity_z !
OPEN(UNIT=970, FILE = './profile/hydro_vphi.dat', ACTION='READ')
READ(970,*) ((prim(ivz,j,k,1), j = 1, nx), k = 1, ny)
CLOSE(970)

prim(ivz,:,:,:) = prim(ivz,:,:,:)*lencgs2code/tcgs2code
PRINT *, "Finished reading vphi"

! Read for magnetic vector potential !
OPEN(UNIT=970, FILE = './profile/hydro_Aphi.dat', ACTION='READ')
READ(970,*) ((a_phi(j,k,1), j = 0, nx), k = 0, ny)
CLOSE(970)

PRINT *, "Finished reading Aphi"

! Calculate magnetic field !
! Coordinate here are in code unit but aphi is in gaussian unit !
! Unit conversion below !
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
      prim(ibx,j,k,l) = (sinf(k)*a_phi(j,k,l) - sinf(k-1)*a_phi(j,k-1,l))/(xF(j)*dcose(k)+small_num)
      prim(iby,j,k,l) = - (xF(j)*a_phi(j,k,l) - xF(j-1)*a_phi(j-1,k,l))/(x(j)*dx(j))
    END DO
  END DO
END DO
prim(ibx:iby,:,:,:) = prim(ibx:iby,:,:,:)*gauss2code*lencgs2code  ! length conversion for curl !

PRINT *, "Finished calculating poloidal field"

! OPEN(UNIT=970, FILE = './profile/hydro_bphi.dat', ACTION='READ')
! READ(970,*) ((prim(ibz,j,k,1), j = 1, nx), k = 1, ny)
! CLOSE(970)
! prim(ibz,:,:,:) = prim(ibz,:,:,:)*gauss2code

PRINT *, "Finished reading torodial field "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate
Deallocate(a_phi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL FINDPRESSURE
PRINT *, "Finished calculating pressure and epsilon"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find DIVB !
div_b = 0.0D0
maxdb = 0.0D0
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      div_b = (xF(j)*xF(j)*prim(ibx,j,k,l) - xF(j-1)*xF(j-1)*prim(ibx,j-1,k,l))/(dx_cb(j)/3.0d0) &
            + (sinf(k)*prim(iby,j,k,l) - sinf(k-1)*prim(iby,j,k-1,l))*(x(j)*dx(j))/(dx_cb(j)*dcose(k)/3.0d0) &
            + (prim(ibz,j,k,l) - prim(ibz,j,k,l-1))*(x(j)*dx(j)*dy(k))/(dx_cb(j)*dcose(k)*dz(l)/3.0d0)
      maxdb = MAX(maxdb, ABS(div_b))
    END DO
  END DO
END DO
WRITE (*,*)
WRITE (*,*) 'Maximum initial divergence B', maxdb
WRITE (*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
prim_a(:) = 0.0D0
prim_a(irho) = atmosphere
prim_a(itau) = prim(itau,nx,1,1)
eps_a = epsilon(nx,1,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign floor density !
DO j = 1, nx
  DO k = 1, ny
    IF(prim(irho,j,k,1) < atmosphere) THEN
      prim(irho,j,k,1) = atmosphere
      prim(ivx:ivz,j,k,1) = 0.0d0
      epsilon(j,k,1) = eps_a
    END IF
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set output profile interval !
total_time = 100.0D0 !0.05d0*tcgs2code
output_profiletime = total_time/100.0d0

END SUBROUTINE
