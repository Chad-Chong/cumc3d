!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model                                                           !
! Simulation for rotating magnetised white dwarf type Ia                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE CUSTOM_DEF
USE FlameTable_module
USE Ecaptable_module
USE Helmeos_module
USE TURB_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l, m, flag_eostable
REAL*8 :: dummy
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

PRINT *, "Finished reading torodial field"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate
Deallocate(a_phi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Chemical Composition !

IF (xisotran_flag == 1) THEN

  CALL read_helm_table()
  CALL initialize_network()
  WRITE(*,*), 'Finished reading helmeos, flame and nse table'

  DO i = ihe4, ini56
    m = i+1-ihe4
    OPEN(UNIT=970, FILE = trim('./profile/Helm_X'//ionam(m))//'.dat', ACTION='READ')
    READ(970,*) ((prim(i,j,k,1), j = 1, nx), k = 1, ny)
    CLOSE(970)
  ENDDO
  IF (helmeos_flag == 1) THEN
    CALL FIND_AZBAR()
  ENDIF
  PRINT *, "Finished reading chemical composition"

ENDIF

IF (helmeos_flag == 1) THEN
  OPEN(UNIT=970, FILE = trim('./profile/Helm_temp.dat'), ACTION='READ')
  READ(970,*) ((temp2(j,k,1), j = 1, nx), k = 1, ny)
  CLOSE(970)
  temp2_old = temp2
  PRINT *, "Finished reading temperature for Helmholtz EOS"

  OPEN(UNIT=970, FILE = trim('./profile/Helm_Ye.dat'), ACTION='READ')
  READ(970,*) ((prim(iye2,j,k,1), j = 1, nx), k = 1, ny)
  CLOSE(970)
  PRINT *, "Finished reading Ye for Helmholtz EOS (electron capture)"
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
prim_a(:) = 0.0D0
atmosphere = 1.0D-8*MAXVAL(prim(irho,:,:,:))
prim_a(irho) = atmosphere

IF (helmeos_flag == 1) THEN
  CALL HELM_EOSPRESSURE(atmosphere, temp2_a, abar2(nx,1,1),  zbar2(nx,1,1), prim(iye2, nx, 1, 1), prim_a(itau), dummy, dummy, flag_eostable)
  CALL HELM_EOSEPSILON(atmosphere, temp2_a, abar2(nx,1,1), zbar2(nx,1,1), prim(iye2, nx, 1, 1), eps_a)
ELSE
  prim_a(itau) = prim(itau,nx,1,1)
  eps_a = epsilon(nx,1,1)
ENDIF

IF (helmcheck_flag == 1) THEN
    WRITE(*,*) 'Atmosphere rho is', atmosphere, 'epsilon is', eps_a, 'pressure is', prim_a(itau)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign floor density !
DO j = 1, nx
  DO k = 1, ny
    IF(prim(irho,j,k,1) < atmosphere) THEN
      prim(irho,j,k,1) = atmosphere
      prim(ivx:ivz,j,k,1) = 0.0d0
      epsilon(j,k,1) = eps_a
      IF (helmeos_flag == 1) THEN
        temp2(j,k,1) = temp2_a
        temp2_old(j,k,1) = temp2_a
      ENDIF
    END IF
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sub-grid Turbulence !

IF (turb_flag == 1) THEN
  CALL GetTurb
  CALL FINDTURBULENCE
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL FINDPRESSURE
IF (helmeos_flag == 0) THEN
  DO j = 1, nx
    DO k = 1, ny
      prim(irho,j,k,1) = max(prim(irho,j,k,1), atmosphere)
      CALL EOSEPSILON_NM(prim(irho,j,k,1), prim(itau,j,k,1), epsilon(j,k,1))
    END DO
  END DO
ENDIF
PRINT *, 'Finish calculating pressure, sound speed and epsilon'

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

! Set output profile interval !
total_time = 0.0005d0*tcgs2code
output_profiletime = total_time/100.0d0
END SUBROUTINE
