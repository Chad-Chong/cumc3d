!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model                                                           !
! Simulation for rotating magnetised white dwarf type Ia                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE CUSTOM_DEF
USE MHD_MODULE
USE IEEE_ARITHMETIC
USE Ecaptable_module
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l, m, flag_eostable, found_atmo, flag_notfindtemp
REAL*8 :: dummy
! Magnetic field !
REAL*8 :: maxdb
REAL*8 :: div_b

REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: a_phi

! For Atmosphere
REAL*8 :: diff, factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flag Checking !
IF (axissym_flag == 1 .and. coordinate_flag == 1) THEN

  IF (ny > 1) THEN
    WRITE(*,*) 'For axis symmetry there can only be one grid in the azimuth direction'
    STOP
  ENDIF

  IF (n_dim /= 3) THEN
    WRITE(*,*) 'For axis symmetry in cylindrical coordinates n_dim needs to be 3'
    STOP
  ENDIF
ELSE
  WRITE(*,*) 'axissym_flag only imposes axis symmetry for cylindrical coordinates'
  STOP  
ENDIF

IF (turb_flag == 1 .and. coordinate_flag /= 1) THEN
  WRITE(*,*) 'SGS Turbulence is only implemented for cylindrical coordinates at the moment'
  STOP 
ENDIF

IF (helmeos_flag == 1 .and. xisotran_flag /= 1) THEN
  WRITE(*,*) 'Helm eos must be initiated with isotope tracking'
  STOP
ENDIF

IF (levelset_flag == 1 .and. coordinate_flag /= 1) THEN
  WRITE(*,*) 'The level set method (2D) for flame tracking is only implemented for cylindrical coordinates at the moment'
  STOP 
ENDIF

IF (flame_flag == 1 .and. coordinate_flag /= 1) THEN
  WRITE(*,*) 'The level set method (2D) for flame tracking is only implemented for cylindrical coordinates at the moment'
  STOP 
ENDIF

IF (flame_flag == 1 .and. xisotran_flag /= 1) THEN
  WRITE(*,*) 'The flame must be initiated with the helm eos with isotope tracking'
  STOP
ENDIF

IF (coordinate_flag == 0) THEN
  WRITE(*,*) 'This code is not written for Cartesian coordinates'
  STOP
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preperation !

! Allocate
Allocate(a_phi(-3:nx+3,-3:ny+3,-3:nz+3))

! Poisson interpolation coefficient !
call get_poisson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read and assign density !
OPEN(UNIT=970, FILE = './profile/hydro_rho.dat', ACTION='READ')
IF (coordinate_flag == 2) THEN
  READ(970,*) ((prim(irho,j,k,1), j = 1, nx), k = 1, ny)
ELSEIF (coordinate_flag == 1) THEN
  READ(970,*) ((prim(irho,j,1,l), j = 1, nx), l = 1, nz)
ENDIF
CLOSE(970)
prim(irho,:,:,:) = prim(irho,:,:,:)*rhocgs2code
PRINT *, "Finished reading rho"
WRITE(*,*) 'Maximum rho is', MAXVAL(prim(irho,:,:,:))/rhocgs2code

! Assign velocity !
IF (coordinate_flag == 2) THEN
  prim(ivx,:,:,1) = 0.0d0
  prim(ivy,:,:,1) = 0.0d0
ELSEIF (coordinate_flag == 1) THEN
  prim(ivx,:,1,:) = 0.0d0
  prim(ivz,:,1,:) = 0.0d0
ENDIF

IF (rotate_flag == 1) THEN
  ! Read and assign azimuth direction velocity !
  OPEN(UNIT=970, FILE = './profile/hydro_vphi.dat', ACTION='READ')
  IF (coordinate_flag == 2) THEN
    READ(970,*) ((prim(ivz,j,k,1), j = 1, nx), k = 1, ny)
    prim(ivz,:,:,:) = prim(ivz,:,:,:)*lencgs2code/tcgs2code
  ELSEIF (coordinate_flag == 1) THEN
    READ(970,*) ((prim(ivy,j,1,l), j = 1, nx), l = 1, nz)
    prim(ivy,:,:,:) = prim(ivy,:,:,:)*lencgs2code/tcgs2code
  ENDIF
  CLOSE(970)
  PRINT *, "Finished reading vphi"
ELSE
  IF (coordinate_flag == 2) THEN
    prim(ivz,:,:,:) = 0.0D0
  ELSEIF (coordinate_flag == 1) THEN
    prim(ivy,:,:,:) = 0.0D0
  ENDIF
  PRINT *, "No rotation in this run"
ENDIF

IF (mhd_flag == 1) THEN
  ! Read for magnetic vector potential !
  OPEN(UNIT=970, FILE = './profile/hydro_Aphi.dat', ACTION='READ')
  IF (coordinate_flag == 2) THEN
    READ(970,*) ((a_phi(j,k,1), j = 0, nx), k = 0, ny)
  ELSEIF (coordinate_flag == 1) THEN
    READ(970,*) ((a_phi(j,1,l), j = 0, nx), l = 0, nz)
    a_phi(j,:,l) = a_phi(j,1,l)
  ENDIF

  CLOSE(970)

  PRINT *, "Finished reading Aphi"

  ! In the direction of symmetry, cell centered is the same as face centered
  OPEN(UNIT=970, FILE = './profile/hydro_bphi.dat', ACTION='READ')
  IF (coordinate_flag == 2) THEN
    READ(970,*) ((prim(ibz,j,k,1), j = 1, nx), k = 1, ny)
    prim(ibz,:,:,0) = prim(ibz,:,:,1)
  ELSEIF (coordinate_flag == 1) THEN
    READ(970,*) ((prim(iby,j,1,l), j = 1, nx), l = 1, nz)
    prim(iby,:,0,:) = prim(iby,:,1,:)
    prim(iby,:,:,:) = prim(iby,:,:,:)*gauss2code*lencgs2code
  ENDIF

  CLOSE(970)

  PRINT *, "Finished reading Bphi"

  ! Coordinate here are in code unit but aphi is in gaussian unit !
  ! Unit conversion below !
  IF (coordinate_flag == 2) THEN
    DO l = 1, nz
      DO k = 1, ny
        DO j = 0, nx
          prim(ibx,j,k,l) = (sinf(k)*a_phi(j,k,l) - sinf(k-1)*a_phi(j,k-1,l))/(xF(j)*dcose(k)+small_num)
        END DO
      END DO
    END DO

    DO l = 1, nz
      DO k = 0, ny
        DO j = 1, nx
          prim(iby,j,k,l) = - (xF(j)*a_phi(j,k,l) - xF(j-1)*a_phi(j-1,k,l))/(x(j)*dx(j))
        END DO
      END DO
    END DO

    prim(ibx:ibz,:,:,:) = prim(ibx:iby,:,:,:)*gauss2code*lencgs2code  ! length conversion for curl !
  ELSEIF(coordinate_flag == 1) THEN
    DO l = 1, nz
      DO k = 1, ny
        DO j = 0, nx
          prim(ibx,j,k,l) = - (a_phi(j,k,l) - a_phi(j,k,l-1))/(dz(l))
        END DO
      END DO
    END DO

    DO l = 0, nz
      DO k = 1, ny
        DO j = 1, nx
          prim(ibz,j,k,l) = (xF(j)*a_phi(j,k,l) - xF(j-1)*a_phi(j-1,k,l))/(x(j)*dx(j))
        END DO
      END DO
    END DO

    prim(ibx,:,:,:) = prim(ibx,:,:,:)*gauss2code*lencgs2code  ! length conversion for curl !
    prim(ibz,:,:,:) = prim(ibz,:,:,:)*gauss2code*lencgs2code  ! length conversion for curl !
  ENDIF

  PRINT *, "Finished calculating magnetic field"
ELSE
  prim(ibx:ibz,:,:,:) = 0.0D0
  PRINT *, "No magnetic field in this run"
ENDIF
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
    IF (coordinate_flag == 2) THEN
      READ(970,*) ((prim(i,j,k,1), j = 1, nx), k = 1, ny)
    ELSEIF (coordinate_flag == 1) THEN
      READ(970,*) ((prim(i,j,1,l), j = 1, nx), l = 1, nz)
    ENDIF
    CLOSE(970)
  ENDDO
  IF (helmeos_flag == 1) THEN
    CALL FIND_AZBAR()
  ENDIF
  PRINT *, "Finished reading chemical composition"

ENDIF

IF (helmeos_flag == 1) THEN
  OPEN(UNIT=970, FILE = trim('./profile/Helm_temp.dat'), ACTION='READ')
  IF (coordinate_flag == 2) THEN
    READ(970,*) ((temp2(j,k,1), j = 1, nx), k = 1, ny)
  ELSEIF (coordinate_flag == 1) THEN
    READ(970,*) ((temp2(j,1,l), j = 1, nx), l = 1, nz)
  ENDIF
  CLOSE(970)
  temp2_old = temp2
  PRINT *, "Finished reading temperature for Helmholtz EOS"

  OPEN(UNIT=970, FILE = trim('./profile/Helm_Ye.dat'), ACTION='READ')
  IF (coordinate_flag == 2) THEN
    READ(970,*) ((prim(iye2,j,k,1), j = 1, nx), k = 1, ny)
  ELSEIF (coordinate_flag == 1) THEN
    READ(970,*) ((prim(iye2,j,1,l), j = 1, nx), l = 1, nz)
  ENDIF
  CLOSE(970)
  PRINT *, "Finished reading Ye for Helmholtz EOS (electron capture)"
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL FINDPRESSURE

IF (helmeos_flag == 1) THEN
	DO l = 1, nz
		DO k = 1, ny
			DO j = 1, nx

				CALL HELM_EOSEPSILON(prim(irho,j,k,l), temp2(j,k,l), abar2(j,k,l), zbar2(j,k,l), prim(iye2,j,k,l), epsilon(j,k,l))
				IF (ieee_is_nan(epsilon(j,k,l))) THEN
					WRITE(*,*) 'Global time', global_time, 'Input Rho', prim(irho,j,k,l), 'Input temp', temp2(j,k,l), 'abar2', abar2(j,k,l), 'zbar2', zbar2(j,k,l), 'ye', prim(iye2, j, k, l), 'at j,k,l', j,k,l
					STOP
				ENDIF
				
			END DO
		END DO
	END DO
	CALL BOUNDARY1D_NM(epsilon, even, even, even, even, even, even)
ENDIF

IF (helmeos_flag == 0 .and. coordinate_flag == 2) THEN
  DO j = 1, nx
    DO k = 1, ny
      prim(irho,j,k,1) = max(prim(irho,j,k,1), atmosphere)
      CALL EOSEPSILON_NM(prim(irho,j,k,1), prim(itau,j,k,1), epsilon(j,k,1))
    END DO
  END DO
ELSEIF (helmeos_flag == 0 .and. coordinate_flag == 1) THEN 
  DO j = 1, nx
    DO l = 1, nz
      prim(irho,j,k,1) = max(prim(irho,j,1,l), atmosphere)
      CALL EOSEPSILON_NM(prim(irho,j,1,l), prim(itau,j,1,l), epsilon(j,1,l))
    END DO
  END DO
ENDIF

PRINT *, 'Finish calculating pressure, sound speed and epsilon'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sub-grid Turbulence !

IF (turb_flag == 1) THEN
  IF (ABS(dx(1) - dz(1)) > 1e-3*lencgs2code) THEN
    WRITE(*,*) 'dx(1)=', dx(1), 'dz(1)=', dz(1)
    WRITE(*,*) 'Grid is not uniform for cylindrical sub-grid scale turbulence'
    STOP
  ENDIF
  CALL GetTurb
  CALL FINDTURBULENCE
ENDIF
CALL BOUNDARY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
prim_a(:) = 0.0D0
IF (turb_flag == 1) THEN
  prim_a(iturbq) = turb_q_a
ENDIF
atmosphere = atmospheric*MAXVAL(prim(irho,:,:,:))
prim_a(irho) = atmosphere

IF (helmeos_flag == 1) THEN
  prim_a(ihe4) = xiso_ahe4
  prim_a(ic12) = xiso_ac12
  prim_a(io16) = xiso_ao16
  CALL PRIVATE_HELMEOS_AZBAR(prim_a(ihe4:ini56), abar2_a, zbar2_a, prim_a(iye2))

  found_atmo = 0 ! Find atmosphere flag

  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        ! Standard !
        diff = prim(irho,j,k,l) - prim_a(irho)
        factor = MAX(SIGN(1.0D0, diff), 0.0D0)

        IF (factor == 0.0D0 .and. found_atmo == 0) THEN
          found_atmo = 1
          eps_a = epsilon(j,k,l)
          WRITE(*,*) 'Atmosphere epsilon is', eps_a
          CALL private_invert_helm_ed(epsilon(j,k,l), &
                            prim_a(irho), abar2_a, &
                            zbar2_a, prim_a(iye2), &
                                            temp2_old(j,k,l), temp2_a, flag_notfindtemp)
          CALL HELM_EOSPRESSURE(atmosphere, temp2_a, abar2_a, zbar2_a, prim_a(iye2), prim_a(itau), dummy, dummy, flag_eostable)
          CALL HELM_EOSEPSILON(atmosphere, temp2_a, abar2_a, zbar2_a, prim_a(iye2), eps_a)
        ENDIF

        prim(imin:ibx-1,j,k,l) = factor*prim(imin:ibx-1,j,k,l) + (1.0D0 - factor)*prim_a(:)
        temp2(j,k,l) = factor*temp2(j,k,l) + (1.0D0 - factor)*temp2_a
        abar2(j,k,l) = factor*abar2(j,k,l) + (1.0D0 - factor)*abar2_a
        zbar2(j,k,l) = factor*zbar2(j,k,l) + (1.0D0 - factor)*zbar2_a
        epsilon(j,k,l) = factor*epsilon(j,k,l) + (1.0D0 - factor)*eps_a
      END DO
    END DO
  END DO

  ELSE
    prim_a(itau) = prim(itau,nx,1,1)
    eps_a = epsilon(nx,1,1)
ENDIF

IF (helmcheck_flag == 1) THEN
    WRITE(*,*) 'Atmosphere rho is', atmosphere, 'epsilon is', eps_a, 'pressure is', prim_a(itau), 'abar is', abar2_a, 'zbar is', zbar2_a, 'Ye is',  prim_a(iye2), 'temp is', temp2_a  
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign floor variables !
CALL CUSTOM_CHECKRHO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (xisotran_flag == 1 .and. helmeos_flag == 1 .and. levelset_flag == 1 .and. flame_flag == 1) THEN
  CALL GetFlame
  CALL FLAME_INI()
  CALL FIND_AZBAR()	
  CALL FINDHELMTEMP()
  WRITE(*,*) 'Finished initializing flame'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! OPEN (UNIT = 123, FILE = './rho_temp_bfcr.dat', STATUS = 'REPLACE')
!     DO l = 1, nz, 1 
!         DO j = 1, nx, 1
!             WRITE(123, *) prim(irho,j,1,l), temp2(j,1,l)
!         ENDDO
!     ENDDO
! CLOSE(123)

! Assign floor variables !
CALL CUSTOM_CHECKRHO

! OPEN (UNIT = 123, FILE = './rho_temp_afcr.dat', STATUS = 'REPLACE')
!     DO l = 1, nz, 1 
!         DO j = 1, nx, 1
!             WRITE(123, *) prim(irho,j,1,l), temp2(j,1,l)
!         ENDDO
!     ENDDO
! CLOSE(123)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find DIVB !
CALL FIND_DIVB
WRITE (*,*)
WRITE (*,*) 'Maximum initial divergence B', maxDivB
WRITE (*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set output profile interval !
total_time = 2.0D0*tcgs2code ! cgs
output_profiletime = total_time/100.0d0
END SUBROUTINE
