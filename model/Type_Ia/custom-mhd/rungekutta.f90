!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine does one single Runge-Kutta full step
! It uses the opeator splitting and separate
! all non-gravitational source term to be done 
! after the hydro step.
! Written by Leung Shing Chi in 2016
! Updated by Leung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RUNGEKUTTA
USE CUSTOM_DEF
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l
! Dummy !
REAL*8 :: rhoaold, dummy

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end, time1, time0
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time0)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

#ifdef DEBUG
CALL system_clock(time_start)
#endif

! Backup old arrays !
!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
			DO i = imin, imax 
				u_old (i,j,k,l) = cons (i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'backup = ', REAL(time_end - time_start) / rate
#endif


IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 1st backup'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1st iteration

! Discretize !
CALL SPATIAL

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 1st spatial'
ENDIF


#ifdef DEBUG
CALL system_clock(time_start)
#endif
! NM sector !
!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
			DO i = imin, imax 
				cons (i,j,k,l) = u_old (i,j,k,l) + dt * l_rk (i,j,k,l) 
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'rk1 = ', REAL(time_end - time_start) / rate
#endif


IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 1st opearation'
ENDIF

! Convert from conservative to primitive !
CALL FROMUTORVE

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 1st FROMUTORVE'
ENDIF

! Check quantities !
CALL CUSTOM_CHECKRHO

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 1st CUSTOM_CHECKRHO'
ENDIF

! set boundary conditions !
call BOUNDARY

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 1st Boundary'
ENDIF

! Update 
CALL UPDATE (1)

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 1st Updat'
ENDIF

! set boundary conditions !
call BOUNDARY

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 1st Boundary'
ENDIF

! Convert from primitive to conservative !
CALL FROMRVETOU

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 1st FROMRVETOU'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2nd iteration

! Discretize !
CALL SPATIAL

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 2nd spatial'
ENDIF

#ifdef DEBUG
CALL system_clock(time_start)
#endif

! NM sector !
!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
			DO i = imin, imax 
				cons (i,j,k,l) = rk20 * u_old(i,j,k,l) + rk21 * cons (i,j,k,l) + rk22 * dt * l_rk (i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'rk2 = ', REAL(time_end - time_start) / rate
#endif

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 2nd operation'
ENDIF

! Convert from conservative to primitive
CALL FROMUTORVE

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 2nd FROMUTORVE'
ENDIF

! Check quantities !
CALL CUSTOM_CHECKRHO

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 2nd CUSTOM_CHECKRHO'
ENDIF

! set boundary conditions !
call BOUNDARY

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 2nd BOUNDARY'
ENDIF

! Update 
CALL UPDATE (2)

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 2nd UPDATE'
ENDIF

! set boundary conditions !
call BOUNDARY

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 2nd Boundary'
ENDIF

! Convert from primitive to conservative !
CALL FROMRVETOU

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 2nd FROMRVETOU'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare for next step

CALL SPATIAL

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 3rd SPATIAL'
ENDIF

#ifdef DEBUG
CALL system_clock(time_start)
#endif

! NM sector !
!$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(STATIC)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
			DO i = imin, imax 
				cons (i,j,k,l) = rk30 * u_old(i,j,k,l) + rk31 * cons (i,j,k,l) + rk32 * dt * l_rk (i,j,k,l)
			END DO
		END DO
	END DO
END DO 
!$ACC END PARALLEL
!$OMP END PARALLEL DO

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'rk3 = ', REAL(time_end - time_start) / rate
#endif

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 3rd opeartion'
ENDIF

! Convert from conservative to primitive
CALL FROMUTORVE 

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 3rd FROMUTORVE'
ENDIF

! Check quantities !
CALL CUSTOM_CHECKRHO

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 3rd CUSTOM_CHECKRHO'
ENDIF

! set boundary conditions !
call BOUNDARY

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 3rd BOUNDARY'
ENDIF

! Convert from primitive to conservative !
CALL FROMRVETOU

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 3rd FROMRVETOU'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for operator splitting

CALL OPERATOR_SPLIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! uncomment the below if conservative variables are changed with opeartor split

! ! Convert from conservative to primitive
! CALL FROMUTORVE 

! ! Check quantities !
! CALL CUSTOM_CHECKRHO

! ! set boundary conditions !
! call BOUNDARY

! ! Update again !
! CALL FROMRVETOU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for burning (burns only at the last step of RK)
  
IF(levelset_flag == 1 .and. xisotran_flag == 1) THEN

  ! If there is level-set, update it

  CALL UPDATE_FLAME_RADIUS

  	IF (say_flag == 1) THEN
		WRITE(*,*) 'RK: Finished UPDATE_FLAME_RADIUS'
	ENDIF
  
  ! This trigger the burning package proposed by
  ! Reinecke 1999b
  IF(burn_flag == 1) THEN

    ! This does the Carbon burning
  
    IF(carburn_flag == 1) CALL BURN_PHASE1B

	IF (say_flag == 1) THEN
		WRITE(*,*) 'RK: Finished BURN_PHASE1B'
	ENDIF


    ! This do the O- and Si- burning

    IF(advburn_flag == 1) CALL BURN_PHASE2B

	IF (say_flag == 1) THEN
		WRITE(*,*) 'RK: Finished BURN_PHASE2B'
	ENDIF

    ! Update the AZbar and temperature accordingly
    CALL FIND_AZBAR
    CALL FINDHELMTEMP

	IF (say_flag == 1) THEN
		WRITE(*,*) 'RK: Finished EOS find azbar and helmtemp'
	ENDIF

    ! For completely burnt zone, check if NSE applies

    IF(convert_nse_flag == 1) CALL NSE2

	IF (say_flag == 1) THEN
		WRITE(*,*) 'RK: Finished NSE2'
	ENDIF

    ! Copy the new Xiso and epsilon to ghost cells
    CALL BOUNDARY

	IF (say_flag == 1) THEN
		WRITE(*,*) 'RK: Finished final boundary'
	ENDIF

    ! Check if the change of isotope perserve the sum
    !CALL system_clock(time_start2)

    CALL CHECKXISOTOPE

	IF (say_flag == 1) THEN
		WRITE(*,*) 'RK: Finished checkxisotope'
	ENDIF

    ! Update the burntime
    last_burntime = global_time

    ! Update Abar and Zbar and temperature again
    CALL FIND_AZBAR
    CALL FINDHELMTEMP

	IF (say_flag == 1) THEN
		WRITE(*,*) 'RK: Finished EOS find azbar and helmtemp'
	ENDIF

	ENDIF
ENDIF

! Check quantities !
CALL CUSTOM_CHECKRHO

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished CUSTOM_CHECKRHO'
ENDIF

! set boundary conditions !
call BOUNDARY

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished BOUNDARY'
ENDIF

! Update physical quantities
CALL UPDATE (3)

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished 3rd update'
ENDIF

! set boundary conditions !
call BOUNDARY

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished BOUNDARY'
ENDIF

! Do conversion again !
CALL FROMRVETOU

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished FROMRVETOU'
ENDIF

IF (nuspec_flag == 1) THEN
	CALL FINDNUSPEC
ENDIF

IF (say_flag == 1) THEN
	WRITE(*,*) 'RK: Finished FINDNUSPEC'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time1)
WRITE(*,*) 'rk total = ', REAL(time1 - time0) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
! This subroutine calculates the maximum time step
! which satisfies the Courant condition 
! Written by Leung Shing Chi in 2016   
! If you modify the Euler equation, make sure you change this 
! part to include the new effective sound speed
! Limiters are posed based on output time and running time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDDT
USE DEFINITION
USE MHD_MODULE
USE CUSTOM_DEF
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

! For MHD speed !
REAL*8 :: a2_mhd, b2_mhd
REAL*8 :: a4_mhd, b4_mhd
REAL*8 :: b2x_mhd, b2y_mhd, b2z_mhd
REAL*8 :: cfx_mhd, cfy_mhd, cfz_mhd

! Local maximum effective speed
REAL*8 :: lambda, lambda1, lambda2, lambda3

! Local minimum dt for DM, NM and 1st overlayer
REAL*8 :: dt_temp1, dt_temp2

! Local minimum dt for DM, NM and 1st overlayer
REAL*8 :: dt_out1, dt_out2

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set !
dt_out1 = 1.0D5
dt_out2 = 1.0D5
lambdas = -1.0D0

! Now we find the minimum time constrained by NM sector
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) & 
!$OMP PRIVATE(a2_mhd, b2_mhd, a4_mhd, b4_mhd, b2x_mhd, b2y_mhd, b2z_mhd, &
!$OMP cfx_mhd, cfy_mhd, cfz_mhd, lambda, lambda1, lambda2, lambda3, dt_temp2) REDUCTION(MIN:dt_out2)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) &
!$ACC PRIVATE(a2_mhd, b2_mhd, a4_mhd, b4_mhd, b2x_mhd, b2y_mhd, b2z_mhd, & 
!$ACC cfx_mhd, cfy_mhd, cfz_mhd, lambda, lambda1, lambda2, lambda3, dt_temp2) REDUCTION(MIN:dt_out2)
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx

			! Only grid with density above threshold density is counted
			IF (prim(irho,j,k,l) > 1.0D1*prim_a(irho)) THEN
				a2_mhd = cs(j,k,l)*cs(j,k,l)
				a4_mhd = a2_mhd*a2_mhd
				b2x_mhd = (bcell(ibx,j,k,l)*bcell(ibx,j,k,l)/prim(irho,j,k,l))
				b2y_mhd = (bcell(iby,j,k,l)*bcell(iby,j,k,l)/prim(irho,j,k,l))
				b2z_mhd = (bcell(ibz,j,k,l)*bcell(ibz,j,k,l)/prim(irho,j,k,l))
				b2_mhd = b2x_mhd + b2y_mhd + b2z_mhd
				b4_mhd = b2_mhd*b2_mhd
				cfx_mhd = DSQRT(0.5D0*(a2_mhd + b2_mhd + DSQRT((a4_mhd + 2.0d0*a2_mhd*b2_mhd + b4_mhd) - 4.0D0*a2_mhd*b2x_mhd)))
				cfy_mhd = DSQRT(0.5D0*(a2_mhd + b2_mhd + DSQRT((a4_mhd + 2.0d0*a2_mhd*b2_mhd + b4_mhd) - 4.0D0*a2_mhd*b2y_mhd)))
				cfz_mhd = DSQRT(0.5D0*(a2_mhd + b2_mhd + DSQRT((a4_mhd + 2.0d0*a2_mhd*b2_mhd + b4_mhd) - 4.0D0*a2_mhd*b2z_mhd)))
				lambda1 = ABS(prim(ivx,j,k,l)) + cfx_mhd
				lambda2 = ABS(prim(ivy,j,k,l)) + cfy_mhd
				lambda3 = ABS(prim(ivz,j,k,l)) + cfz_mhd
				lambda = MAX(lambda1, lambda2, lambda3)
							
				! Look for minimum grid size !
				dt_temp2 = dx(j)
				IF(coordinate_flag == 0) THEN
					IF(n_dim > 1) THEN
						dt_temp2 = MIN(dt_temp2, dy(k))
					END IF
					IF(n_dim > 2) THEN
						dt_temp2 = MIN(dt_temp2, dz(l))
					END IF
				ELSEIF(coordinate_flag == 1) THEN
					IF(n_dim > 1 .AND. ny > 1) THEN
						dt_temp2 = MIN(dt_temp2, x(j)*dy(k))
					END IF
					IF(n_dim > 2) THEN
						dt_temp2 = MIN(dt_temp2, dz(l))
					END IF
				ELSEIF(coordinate_flag == 2) THEN
					IF(n_dim > 1) THEN
						dt_temp2 = MIN(dt_temp2, x(j)*dy(k))
					END IF
					IF(n_dim > 2) THEN
						dt_temp2 = MIN(dt_temp2, x(j)*sine(k)*dz(l))
					END IF
				END IF
				dt_temp2 = dt_temp2*cfl/lambda
				lambdas(j,k,l) = lambda
				dt_out2 = MIN(dt_out2, dt_temp2)
			END IF
		END DO
	ENDDO
ENDDO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Only the minimum one is chosen
dt = MIN(dt_out1, dt_out2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'finddt = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE FindDt
