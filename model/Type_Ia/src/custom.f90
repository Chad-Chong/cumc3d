!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_EQN(no_of_eq)
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

INTEGER, INTENT(INOUT) :: no_of_eq

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Isotopes !
! Set up the code for chemical composition !

IF (turb_flag == 1) THEN
  no_of_eq = no_of_eq + 1
  imax = no_of_eq
  iturbq = no_of_eq
  WRITE(*,*) 'Make iturbq = ', no_of_eq
END IF

IF (xisotran_flag == 1) THEN
  ! Helium-4
  no_of_eq = no_of_eq + 1
  imax = no_of_eq
  ihe4 = no_of_eq
  WRITE(*,*) 'Make ihe4 = ', no_of_eq

  ! Carbon-12
  no_of_eq = no_of_eq + 1
  imax = no_of_eq
  ic12 = no_of_eq
  WRITE(*,*) 'Make ic12 = ', no_of_eq
  
  ! Oxygen-16
  no_of_eq = no_of_eq + 1
  imax = no_of_eq
  io16 = no_of_eq
  WRITE(*,*) 'Make io16 = ', no_of_eq

  ! Neon-20
  no_of_eq = no_of_eq + 1
  imax = no_of_eq
  ine20 = no_of_eq
  WRITE(*,*) 'Make ine20 = ', no_of_eq

  ! Magnesium-24
  no_of_eq = no_of_eq + 1
  imax = no_of_eq
  img24 = no_of_eq
  WRITE(*,*) 'Make img24 = ', no_of_eq

  ! Silicon-28
  no_of_eq = no_of_eq + 1
  imax = no_of_eq
  isi28 = no_of_eq
  WRITE(*,*) 'Make isi28 = ', no_of_eq

  ! Nickel-56
  no_of_eq = no_of_eq + 1
  imax = no_of_eq
  ini56 = no_of_eq
  WRITE(*,*) 'Make ini56 = ', no_of_eq

  ! Electron Fraction
  no_of_eq = no_of_eq + 1
  imax = no_of_eq
  iye2 = no_of_eq
  WRITE(*,*) 'Make iye2 = ', no_of_eq

  IF (levelset_flag == 1) THEN
    no_of_eq = no_of_eq + 1 ! Deflagration level set
    imax = no_of_eq
    iscaG1 = no_of_eq
    WRITE(*,*) 'Make iscaG1 = ', no_of_eq

    no_of_eq = no_of_eq + 1 ! Deteonation level set
    imax = no_of_eq
    iscaG2 = no_of_eq
    WRITE(*,*) 'Make iscaG2 = ', no_of_eq
  ENDIF
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom arrays !
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_HYDRO
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

! atmospheric values !
ALLOCATE (prim_a(imin:ibx-1))

! gravitational potential energy !
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

IF (turb_flag == 1) THEN
  WRITE(*,*) 'Build sub-grid turbulence variables'
  CALL buildTurb
  WRITE(*,*) 'Done building sub-grid turbulence variables'
  WRITE(*,*)
ENDIF

IF (xisotran_flag == 1) THEN
  WRITE(*,*) 'Build chemical composition & nuclear variables'
  CALL buildHelm
  WRITE(*,*) 'Done building chemical composition & nuclear variables'
  WRITE(*,*)
ENDIF

IF (levelset_flag == 1) THEN
  WRITE(*,*) 'Build level set variables'
  CALL buildLevelSet
  WRITE(*,*) 'Done building level set variables'
  WRITE(*,*)
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Populate custom arrays to GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_POPULATE
USE CUSTOM_DEF 
USE DEFINITION
IMPLICIT NONE

! Now populate all necessary, and reuseable arrays to the graphic cards !
!$ACC enter DATA COPYIN(phi, phi_old, ajp1, ajm1, bkp1, bkm1, clp1, clm1, epsc):q

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Clear custom arrays from GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CLEAR
USE CUSTOM_DEF 
USE DEFINITION
IMPLICIT NONE

! Now we clear memory in the GPU device !
!$ACC exit DATA DELETE(phi, phi_old, ajp1, ajm1, bkp1, bkm1, clp1, clm1, epsc)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_GRID
USE DEFINITION
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
  IF (xF(i) == 0.0D0) THEN
    xF(i) = 1.0D-50 
  ENDIF
ENDDO
CLOSE(970)
xF = xF*lencgs2code
IF (coordinate_flag == 2) THEN
  nlines = 0 
  OPEN (970, file = './profile/hydro_x2_fgrid.dat')
  DO 
    READ (970,*, END=11) 
    nlines = nlines + 1 
  END DO 
  11 CLOSE (970) 

  ! Error message !
  IF(nlines .ne. ny+7) THEN
    WRITE (*,*) 'number of grid faces from files', nlines
    WRITE (*,*) 'number of y grid faces in the program', ny+6
    STOP 'inconsistent number of grid faces, exit'
  END IF

  ! Read !
  OPEN(UNIT=970, FILE = './profile/hydro_x2_fgrid.dat', ACTION='READ')
  DO i = -3, ny+3
    READ(970,*) yF(i)
  ENDDO
  CLOSE(970)
ENDIF

IF (coordinate_flag == 1 .and. axissym_flag == 1) THEN
  nlines = 0 
  OPEN (970, file = './profile/hydro_x3_fgrid.dat')
  DO 
    READ (970,*, END=12) 
    nlines = nlines + 1 
  END DO 
  12 CLOSE (970) 

  ! Error message !
  IF(nlines .ne. nz+7) THEN
    WRITE (*,*) 'number of grid faces from files', nlines
    WRITE (*,*) 'number of z grid faces in the program', ny+6
    STOP 'inconsistent number of grid faces, exit'
  END IF

  ! Read !
  OPEN(UNIT=970, FILE = './profile/hydro_x3_fgrid.dat', ACTION='READ')
  DO i = -3, nz+3
    READ(970,*) zF(i)
    IF (zF(i) == 0.0D0) THEN
    zF(i) = 1.0D-50 
  ENDIF
  ENDDO
  CLOSE(970)
  zF = zF*lencgs2code
ENDIF

END SUBROUTINE
 
!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY_X
USE DEFINITION
USE CUSTOM_DEF
USE MHD_MODULE
IMPLICIT NONE

INTEGER :: i, j, k, l

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

! x boundary !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)   
DO l = - 2, nz + 3
  DO k = -2, ny + 3
    DO j = 1, 3
      prim(ivx,1-j,k,l) = MIN(prim(ivx,1-j,k,l), 0.0D0)
      prim(ivx,nx+j,k,l) = MAX(prim(ivx,nx+j,k,l), 0.0D0)
    END DO
  END DO               
ENDDO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'custom boundary = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY_Y
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY_Z
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Back up fluxes from riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLUX_X
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Back up fluxes from riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLUX_Y
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Back up fluxes from riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLUX_Z
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom variable floor !
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CHECKRHO
USE CUSTOM_DEF
USE MHD_MODULE
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l, flag_eostable
REAL*8 :: dummy

! Threshold for atmosphere density
REAL*8 :: diff, factor, bfield, alven, rho_old, m_local

! Minimum / maximum internal energy density
REAL*8 :: epsilon_temp_min, epsilon_temp_max

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

IF (helmeos_flag == 1) THEN
  atmosphere = MIN(atmosphere, 1.0D-4 * MAXVAL(prim(irho,:,:,:)))
  prim_a(irho) = atmosphere
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx

      ! Standard !
      diff = prim(irho,j,k,l) - prim_a(irho)
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)

      IF (diff <  0.0D0 ) THEN
        prim(irho,j,k,l) = prim_a(irho) ! Change the thermodynamic / hydro properties, not the composition
        ! prim(ivx:ivz,j,k,l) = 0.0D0
        IF (helmeos_flag == 1) THEN
          CALL HELM_EOSPRESSURE(prim(irho,j,k,l), temp2(j,k,l), abar2(j,k,l), zbar2(j,k,l), prim(iye2,j,k,l), prim(itau,j,k,l), dummy, dummy, flag_eostable)
          CALL HELM_EOSEPSILON(prim(irho,j,k,l), temp2(j,k,l), abar2(j,k,l), zbar2(j,k,l), prim(iye2,j,k,l), epsilon(j,k,l))
        ELSE
          epsilon(j,k,l) = factor*epsilon(j,k,l) + (1.0D0 - factor)*epsilon(nx,k,1)
        ENDIF
        IF (turb_flag == 1) THEN
          prim(iturbq,j,k,l) = prim_a(iturbq)
        ENDIF
      ENDIF

      ! This code segment gives a bug where it assigns epsilon_temp_min in the interior of the star (higher density)
      ! IF (helmeos_flag == 1) THEN 
      !   CALL HELM_EOSEPSILON(prim(irho,j,k,l), temp_min, abar2(j,k,l), zbar2(j,k,l), prim(iye2,j,k,l), epsilon_temp_min)
      !   CALL HELM_EOSEPSILON(prim(irho,j,k,l), temp_max, abar2(j,k,l), zbar2(j,k,l), prim(iye2,j,k,l), epsilon_temp_max)
      !     IF (epsilon(j,k,l) < epsilon_temp_min) THEN
      !       epsilon(j,k,l) = epsilon_temp_min
      !       temp2(j,k,l) = temp_min
      !     ENDIF
      !     IF (epsilon(j,k,l) > epsilon_temp_max) THEN
      !       epsilon(j,k,l) = epsilon_temp_max
      !       temp2(j,k,l) = temp_max
      !     ENDIF
      ! ENDIF


    END DO
  END DO
END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom variable floor !
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOMFLOOR
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_SOURCE
USE MHD_MODULE
USE CUSTOM_DEF 
USE DEFINITION
USE IEEE_ARITHMETIC
IMPLICIT NONE

 ! Integers !
INTEGER :: i, j, k,l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For Gravity !

! Derivatives of gravitational potential
  REAL*8 :: dphidx, dphidy, dphidz

  ! Threshold for atmosphere density
  REAL*8 :: factor, diff

  ! Check timing with or without openmp
  #ifdef DEBUG
  INTEGER :: time_start, time_end
  INTEGER :: cr
  REAL*8 :: rate
  CALL system_clock(count_rate=cr)
  rate = REAL(cr)
  CALL system_clock(time_start)
  #endif

IF (coordinate_flag == 2)THEN
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(dphidx, dphidy, dphidz, factor, diff) 
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(dphidx, dphidy, dphidz, factor, diff) 
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        ! Include only non-atmosphere !
        diff = prim(irho,j,k,l) - prim_a(irho)
        factor = MAX(SIGN(1.0D0, diff), 0.0D0)

        ! Gravitational potential of the matter !
        dphidx = first_derivative (x(j-1), x(j), x(j+1), phi(j-1,k,l), phi(j,k,l), phi(j+1,k,l))
        dphidy = first_derivative (y(k-1), y(k), y(k+1), phi(j,k-1,l), phi(j,k,l), phi(j,k+1,l))
        dphidz = first_derivative (z(l-1), z(l), z(l+1), phi(j,k,l-1), phi(j,k,l), phi(j,k,l+1))
            
        ! Add them to the source term !
        sc(ivx,j,k,l) = sc(ivx,j,k,l) - factor*prim(irho,j,k,l)*dphidx
        sc(ivy,j,k,l) = sc(ivy,j,k,l) - factor*prim(irho,j,k,l)*dphidy/x(j)
        sc(ivz,j,k,l) = sc(ivz,j,k,l) - factor*prim(irho,j,k,l)*dphidz/x(j)/sine(k)
        sc(itau,j,k,l) = sc(itau,j,k,l) - factor*prim(irho,j,k,l)* &
                          (prim(ivx,j,k,l)*dphidx + prim(ivy,j,k,l)*dphidy/x(j) + &
                          prim(ivz,j,k,l)*dphidz/x(j)/sine(k))
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END PARALLEL DO
  
ELSEIF(coordinate_flag == 1) THEN
  IF (gravity_flag == 1) THEN
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(dphidx, dphidy, dphidz, factor, diff) 
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(dphidx, dphidy, dphidz, factor, diff) 
    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx
          ! Include only non-atmosphere !
          diff = prim(irho,j,k,l) - prim_a(irho)
          factor = MAX(SIGN(1.0D0, diff), 0.0D0)

          ! Gravitational potential of the matter !
          dphidx = first_derivative (x(j-1), x(j), x(j+1), phi(j-1,k,l), phi(j,k,l), phi(j+1,k,l))
          dphidy = first_derivative (y(k-1), y(k), y(k+1), phi(j,k-1,l), phi(j,k,l), phi(j,k+1,l))
          dphidz = first_derivative (z(l-1), z(l), z(l+1), phi(j,k,l-1), phi(j,k,l), phi(j,k,l+1))
              
          ! Add them to the source term !
          sc(ivx,j,k,l) = sc(ivx,j,k,l) - factor*prim(irho,j,k,l)*dphidx
          sc(ivy,j,k,l) = sc(ivy,j,k,l) - factor*prim(irho,j,k,l)*dphidy/x(j)
          sc(ivz,j,k,l) = sc(ivz,j,k,l) - factor*prim(irho,j,k,l)*dphidz
          sc(itau,j,k,l) = sc(itau,j,k,l) - factor*prim(irho,j,k,l)*(prim(ivx,j,k,l)*dphidx + prim(ivy,j,k,l)*dphidy/x(j)+prim(ivz,j,k,l)*dphidz)

          IF (ieee_is_nan(factor*prim(irho,j,k,l)*dphidx)) THEN
            WRITE(*,*) prim(irho,j,k,l)
            WRITE(*,*) dphidx
            WRITE(*,*) 'dphidx term is nan'
            STOP 
          ENDIF

          IF (ieee_is_nan( factor*prim(irho,j,k,l)*dphidy/x(j))) THEN
            WRITE(*,*) dphidy
            WRITE(*,*) 'dphidy term is nan'
            STOP 
          ENDIF

          IF (ieee_is_nan(factor*prim(irho,j,k,l)*dphidz)) THEN
              WRITE(*,*) dphidz
            WRITE(*,*) 'dphidz term is nan'
            STOP 
          ENDIF
        
        END DO
      END DO
    END DO
    !$ACC END PARALLEL
    !$OMP END PARALLEL DO
  ENDIF
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For Turbulence !
IF (turb_flag == 1) THEN
  DO k = 1, nz, 1
    DO j = 1, ny, 1
      DO i = 1, nx, 1
        sc(itau,i,j,k) = sc(itau,i,j,k) - turb_source(i,j,k) 
        sc(iturbq,i,j,k) = sc(iturbq,i,j,k) + turb_source(i,j,k) + turb_diff(i,j,k)
      ENDDO
    ENDDO
  ENDDO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'custom source = ', REAL(time_end - time_start) / rate
#endif

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

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!
! Do custom updates !
!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_UPDATE (p_in)
USE DEFINITION
USE CUSTOM_DEF
USE MHD_MODULE
IMPLICIT NONE

! Integer !
INTEGER, INTENT (IN) :: p_in

! Integer !
INTEGER :: j, k, l, n

! For poisson solver !
REAL*8 :: abserror, rhs
real*8 :: mono
real*8, dimension(3) :: posit, dipo
real*8, dimension(3,3) :: quad
character(len=99) :: charac_n

! Density threshold !
REAL*8 :: rho_in, factor, diff

! nuceos !
REAL*8 :: eps_out, p_out, eosdummy
! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
WRITE(*,*) 'START custom update'
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For SGS Turbulence !

IF(turb_flag == 1) CALL FINDTURBULENCE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For Gravity !

IF (gravity_flag == 1) THEN
  ! Update gravitational potentials !
  IF (p_in == 0 .OR. MOD(n_step, n_pot) == 0) THEN
    ! special treatment for initial model !
    IF(p_in == 0) THEN
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! First, give a guessing potential !
      IF (coordinate_flag == 2) THEN
      !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
      !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
      DO l = 0, nz+1 
        DO k = 0, ny+1
          DO j = 0, nx+1
            phi(j,k,l) = 0.0d0
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
      !$OMP END PARALLEL DO
      ENDIF

      IF (coordinate_flag == 1) THEN
        CALL multipole_expansion(mono, dipo, quad)
        WRITE(*,*) 'Mass is', mono
        !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
        !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
        DO l = 0, nz+1 
          DO k = 0, ny+1
            DO j = 0, nx+1
              posit(1) = abs(x(j)*DCOS(y(k)))
              posit(2) = abs(x(j)*DSIN(y(k)))
              posit(3) = z(l)
              phi(j,k,l) = -mono/DSQRT(x(j)**2+z(l)**2) - dot_product(dipo, posit)/DSQRT(x(j)**2+z(l)**2)**3 - quadSum(quad, posit, posit)/(2*DSQRT(x(j)**2+z(l)**2)**5)
            END DO
          END DO
        END DO
        !$ACC END PARALLEL
        !$OMP END PARALLEL DO    
      ENDIF

      IF (phitest_flag  == 1) THEN
        OPEN(UNIT = 123, FILE = './BCphi.dat', STATUS = 'REPLACE')
        DO l = 1, nz, 1 
            DO j = 1, nx, 1
                WRITE(123, *) prim(irho, j, 1, l)
            ENDDO
        ENDDO
        CLOSE(123)

        OPEN(UNIT = 123, FILE = './vol.dat', STATUS = 'REPLACE')
        DO l = 1, nz, 1 
            DO j = 1, nx, 1
                WRITE(123, *) vol(j, 1, l)
            ENDDO
        ENDDO
        CLOSE(123)

        OPEN(UNIT = 123, FILE = './z.dat', STATUS = 'REPLACE')
        DO l = 0, nz, 1 
          WRITE(123,*) zF(l)
        ENDDO
        CLOSE(123)

        OPEN(UNIT = 123, FILE = './x.dat', STATUS = 'REPLACE')
          DO j = 0, nx, 1
              WRITE(123, *) xF(j)
          ENDDO
        CLOSE(123)

        OPEN(UNIT = 123, FILE = './y.dat', STATUS = 'REPLACE')
        DO k = 0, ny, 1 
          WRITE(123,*) yF(k)
        ENDDO
        CLOSE(123)

        OPEN(UNIT = 123, FILE = './dy.dat', STATUS = 'REPLACE')
        DO k = 1, ny, 1 
          WRITE(123,*) dy(k)
        ENDDO
        CLOSE(123)

        WRITE(*,*) 'ny is', ny

        STOP
      ENDIF
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calucaltes potential by RBSOR
    DO n = 1, relax_max
      !$OMP PARALLEL PRIVATE(diff,rhs,rho_in,factor)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Back up potential !
      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
      !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
      DO l = 1, nz
        DO k = 1, ny
          DO j = 1, nx
            phi_old(j,k,l) = phi(j,k,l)
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
      !$OMP END DO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set error !
      !$OMP SINGLE
      !$ACC SERIAL
      abserror = 1.0D-50
      !$ACC END SERIAL
      !$OMP END SINGLE

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Red chess !
      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
      !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(diff,rhs,rho_in,factor)
      DO l = 1, nz
        DO k = 1, ny
          DO j = 1, nx
            IF ((-1)**(j+k+l)>0) THEN	
              diff = prim(irho,j,k,l) - prim_a(irho)
              factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
              rho_in = factor*prim(irho,j,k,l)
              rhs = (4.0d0*pi*rho_in - & 
                  (ajp1(j)*phi(j+1,k,l) + ajm1(j)*phi(j-1,k,l) + & 
                  bkp1(j,k)*phi(j,k+1,l) + bkm1(j,k)*phi(j,k-1,l) + & 
                  clp1(j,k,l)*phi(j,k,l+1) + clm1(j,k,l)*phi(j,k,l-1)))/epsc(j,k,l)
              phi(j,k,l) = (1.0d0 - omega_weight)*phi(j,k,l) + omega_weight*rhs
            ELSE 
              CYCLE
            END IF
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
      !$OMP END DO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Black chess !
      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
      !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(diff,rhs,rho_in,factor)
      DO l = 1, nz
        DO k = 1, ny
          DO j = 1, nx
            IF ((-1)**(j+k+l)<0) THEN
              diff = prim(irho,j,k,l) - prim_a(irho)
              factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
              rho_in = factor*prim(irho,j,k,l)
              rhs = (4.0d0*pi*rho_in - & 
                  (ajp1(j)*phi(j+1,k,l) + ajm1(j)*phi(j-1,k,l) + & 
                  bkp1(j,k)*phi(j,k+1,l) + bkm1(j,k)*phi(j,k-1,l) + & 
                  clp1(j,k,l)*phi(j,k,l+1) + clm1(j,k,l)*phi(j,k,l-1)))/epsc(j,k,l)
              phi(j,k,l) = (1.0d0 - omega_weight)*phi(j,k,l) + omega_weight*rhs
            ELSE 
              CYCLE
            END IF
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
      !$OMP END DO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Look for maximum abserror !
      !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
      !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) Reduction(MAX:abserror)
      DO l = 1, nz
        DO k = 1, ny
          DO j = 1, nx
            abserror = max(abserror, abs((phi(j,k,l) - phi_old(j,k,l)) / phi_old(j,k,l)))
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
      !$OMP END DO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Boundary conditions !
      IF (coordinate_flag == 2) THEN

        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
        !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
        DO l = 1, nz
          DO k = 1, ny
            phi(0,k,l) = phi(1,k,l)
            phi(nx+1,k,l) = 0.0d0
          END DO
        END DO
        !$ACC END PARALLEL
        !$OMP END DO

        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
        !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
        DO l = 1, nz
          DO j = 1, nx
            phi(j,0,l) = phi(j,1,l)
            phi(j,ny+1,l) = phi(j,ny,l)
          END DO
        END DO
        !$ACC END PARALLEL
        !$OMP END DO

        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
        !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
        DO k = 1, ny
          DO j = 1, nx
            phi(j,k,0) = phi(j,k,1)
            phi(j,k,nz+1) = phi(j,k,nz)
          END DO
        END DO
        !$ACC END PARALLEL
        !$OMP END DO

      ELSEIF(coordinate_flag == 1) THEN
        CALL multipole_expansion(mono, dipo, quad)

        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
        !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
        DO l = 1, nz
          DO k = 1, ny
            phi(0,k,l) = phi(1,k,l)
          END DO
        END DO
        !$ACC END PARALLEL
        !$OMP END DO

        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
        !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
        DO l = 1, nz
          DO j = 1, nx
            phi(j,0,l) = phi(j,1,l)
            phi(j,ny+1,l) = phi(j,ny,l)
          END DO
        END DO
        !$ACC END PARALLEL
        !$OMP END DO
      ENDIF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !$OMP END PARALLEL
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Debug and exit !
      !WRITE (*,*) n, abserror
      IF(abserror <= tolerance) EXIT 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Stop condition !
      IF(n == relax_max) THEN
        WRITE (*,*) n, relax_max
        STOP 'Convergence error in poisson solver'
      END IF
    END DO
  ENDIF
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'custom update = ', REAL(time_end - time_start) / rate
#endif

contains
  REAL*8 FUNCTION quadSum(vec_ab, vec_i, vec_j)
    !$ACC routine seq
    IMPLICIT NONE
    REAL*8, DIMENSION (3,3) :: vec_ab
    REAL*8, DIMENSION (3) :: vec_i, vec_j
    INTEGER :: i, j
    quadSum = 0.0d0
    DO i = 1, 3
      DO j = 1, 3
        quadSum = quadSum + vec_ab(i,j)*vec_i(i)*vec_j(j)
      END DO
    END DO
  END FUNCTION quadSum

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPERATOR_SPLIT
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

INTEGER :: j, k, l

! output_file = .true.

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENFILE_CUSTOM
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

! OPEN !
OPEN (UNIT = 600, FILE = './outfile/density_analysis.dat')
WRITE(600,'(a)') 'Time, Max density, Central density'
FLUSH(600)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_ANALYSIS
USE DEFINITION
USE CUSTOM_DEF
USE CUSTOM_DEF
IMPLICIT NONE

REAL*8 :: graviE, intE, polMagE, torMagE, totMagE, rotE, totKE, totE
INTEGER :: j, k, l

REAL*8 :: maxrho, cenrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (energy_analysis_flag == 1) THEN
  ! Gravitational potential energy !
  graviE = 0.0D0
  ! Internal energy !
  intE = 0.0D0
  ! Poloidal magnetic energy !
  polMagE = 0.0D0
  ! Toroidal magnetic energy !
  torMagE = 0.0D0
  ! Total magnetic energy !
  totMagE = 0.0D0
  ! Rotational energy !
  rotE = 0.0D0
  ! Total kinetic energy !
  totKE = 0.0D0
  ! Total energy !
  totE = 0.0D0

  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
  DO l = 1, nz
    DO k = 1, ny
      DO j = 1, nx
        graviE  = graviE  + 0.5d0*prim(irho,j,k,l)*vol(j,k,l)* phi(j,k,l)
        intE    = intE    +                       vol(j,k,l)* prim(itau,j,k,l)
        rotE    = rotE    + 0.5d0*prim(irho,j,k,l)*vol(j,k,l)* prim(ivz,j,k,l)**2
        totKE   = totKE   + 0.5d0*prim(irho,j,k,l)*vol(j,k,l)*(prim(ivx,j,k,l)**2 + prim(ivy,j,k,l)**2 + prim(ivz,j,k,l)**2)
        polMagE = polMagE + 0.5D0                 *vol(j,k,l)*(prim(ibx,j,k,l)**2  + prim(iby,j,k,l)**2)
        torMagE = torMagE + 0.5D0                 *vol(j,k,l)* prim(ibz,j,k,l)**2
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
  !$OMP END DO

  polMagE = polMagE/(4.0D0*pi)
  torMagE = torMagE/(4.0D0*pi)
  totMagE = polMagE + torMagE
  totE = graviE + totMagE + totKE + intE

  WRITE(601,'(9es30.20e2)') global_time, graviE, intE, polMagE, torMagE, totMagE, rotE, totKE, totE
  FLUSH(601)
ENDIF
maxrho = MAXVAL(prim(irho,:,:,1))
cenrho = prim(irho,0,0,1)

! output !
WRITE (600,*) global_time/tcgs2code, maxrho/rhocgs2code, cenrho/rhocgs2code
! FLUSH(600)

END SUBROUTINE
