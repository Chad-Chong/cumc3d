!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine contains the essential for the calculation of      !
! pressure, sound speed and other thermodynamic quantities.          !
! The EOSEPSILON and EOSSOUNDSPEED are required in reconstruction    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDPRESSURE
USE DEFINITION
USE CUSTOM_DEF
USE HELMEOS_MODULE
USE IEEE_ARITHMETIC
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l, flag_eostable
! Real !
REAL*8 :: xe, dummy1, dummy2


IF (helmeos_flag == 1) THEN
	DO l = 1, nz
		DO k = 1, ny
			DO j = 1, nx
				flag_eostable = 1
				CALL HELM_EOSPRESSURE(prim(irho,j,k,l), temp2(j,k,l), abar2(j,k,l), zbar2(j,k,l), prim(iye2, j, k, l), prim(itau,j,k,l), dummy1, dummy2, flag_eostable)
				IF (ieee_is_nan(prim(itau,j,k,l))) THEN
					WRITE(*,*) 'Global time', global_time, 'Input Rho', prim(irho,j,k,l), 'Input temp', temp2(j,k,l), 'abar2', abar2(j,k,l), 'zbar2', zbar2(j,k,l), 'ye', prim(iye2, j, k, l), 'at j,k,l', j,k,l
					STOP
				ENDIF
				IF (flag_eostable == 0) THEN
					WRITE(*,*) 'EOS Failure: Pressure'
					STOP
				ENDIF
				CALL HELM_EOSEPSILON(prim(irho,j,k,l), temp2(j,k,l), abar2(j,k,l), zbar2(j,k,l), prim(iye2,j,k,l), epsilon(j,k,l))
				IF (ieee_is_nan(epsilon(j,k,l))) THEN
					WRITE(*,*) 'Global time', global_time, 'Input Rho', prim(irho,j,k,l), 'Input temp', temp2(j,k,l), 'abar2', abar2(j,k,l), 'zbar2', zbar2(j,k,l), 'ye', prim(iye2, j, k, l), 'at j,k,l', j,k,l
					STOP
				ENDIF
				CALL HELM_EOSSOUNDSPEED(prim(irho,j,k,l), temp2(j,k,l), abar2(j,k,l), zbar2(j,k,l), cs(j,k,l))
				IF (ieee_is_nan(cs(j,k,l))) THEN
					WRITE(*,*) 'Global time', global_time, 'Input Rho', prim(irho,j,k,l), 'Input temp', temp2(j,k,l), 'abar2', abar2(j,k,l), 'zbar2', zbar2(j,k,l), 'ye', prim(iye2, j, k, l), 'at j,k,l', j,k,l
					STOP
				ENDIF
			END DO
		END DO
	END DO
	CALL BOUNDARY1D_NM(epsilon, even, even, even, even, even, even)
	CALL BOUNDARY1D_NM(prim(itau,:,:,:), even, even, even, even, even, even)
	CALL BOUNDARY1D_NM(cs, even, even, even, even, even, even)

ELSE
	! Check timing with or without openmp
	#ifdef DEBUG
	INTEGER :: time_start, time_end
	INTEGER :: cr
	REAL*8 :: rate
	CALL system_clock(count_rate=cr)
	rate = REAL(cr)
	CALL system_clock(time_start)
	#endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(xe)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
	DO l = -2, nz + 3
		DO k = -2, ny + 3
			DO j = -2, nx + 3
				xe = (prim(irho,j,k,l)/bmax)**(1.0D0/3.0D0)
				IF (xe > 1.0D-2) THEN
					prim(itau,j,k,l) = amax*large_pressure(xe)
				ELSE
					prim(itau,j,k,l) = amax*small_pressure(xe)
				END IF
				cs(j,k,l) = DSQRT(amax*dpdx(xe)/3.0D0/(prim(irho,j,k,l)**2*bmax)**(1.0D0/3.0D0))
			END DO
		END DO
	END DO
	!$ACC END PARALLEL
	!$OMP END PARALLEL DO

	#ifdef DEBUG
	CALL system_clock(time_end)
	WRITE(*,*) 'findpressure = ', REAL(time_end - time_start) / rate
	#endif

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

  REAL*8 function large_pressure(x)
  !$acc routine seq
	implicit none
	REAL*8 :: x
	large_pressure = x*DSQRT(x**2 + 1.0D0)*(2.0D0*x**2 - 3.0D0) + 3.0D0*log(x + DSQRT(x**2 + 1.0D0))
	end function

	REAL*8 function small_pressure(x)
  !$acc routine seq
	implicit none
	REAL*8 :: x
	small_pressure = 1.6D0*x**5 - (4.0D0/7.0D0)*x**7 + (1.0D0/3.0D0)*x**9 - (5.0D0/2.2D1)*x**11 & 
			+ (3.5D1/2.08D2)*x**13 - (2.1D1/1.6D2)*x**15 + (2.31D2/2.176D3)*x**17
	end function

  REAL*8 function dpdx(x)
	!$acc routine seq
	implicit none
	REAL*8 :: x
	dpdx = 8.0D0*x**4/DSQRT(x**2 + 1.0D0)
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EOSEPSILON_NM (rho_in, p_in, eps_out)
!$acc routine seq
USE DEFINITION
IMPLICIT NONE
include 'param.h'

REAL*8 :: xe
REAL*8, INTENT (IN) :: rho_in, p_in
REAL*8, INTENT (OUT) :: eps_out

xe = (rho_in/bmax)**(1.0D0/3.0D0)
IF(xe > 1.0D-2) THEN
	eps_out = amax*large_energy(xe)/rho_in
ELSE
	eps_out = amax*small_energy(xe)/rho_in
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	REAL*8 function large_energy(x)
	!$acc routine seq
	implicit none
	REAL*8 :: x
	large_energy = 3.0D0*x*DSQRT(x**2 + 1.0D0)*(1.0D0 + 2.0D0*x**2) - 3.0D0*log(x + DSQRT(x**2 + 1.0D0)) - 8.0D0*x**3
	end function

	REAL*8 function small_energy(x)
	!$acc routine seq
	implicit none
	REAL*8 :: x
	small_energy = 8.0D0*x**3 + (1.2D1/5.0D0)*x**5 - (3.0D0/7.0D0)*x**7 + (1.0D0/6.0D0)*x**9 - (1.5D1/1.76D2)*x**11 & 
							 + (2.1D1/4.16D2)*x**13 - (2.1D1/6.40D2)*x**15 + (9.9D1/4.352D3)*x**17 - 8.0D0*x**3
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EOSSOUNDSPEED(p_in, rho_in, eps_in, cs_out)
!$acc routine seq
USE DEFINITION
IMPLICIT NONE
include 'param.h'

REAL*8 :: xe
REAL*8, INTENT (IN) :: p_in, rho_in, eps_in
REAL*8, INTENT (OUT) :: cs_out

xe = (rho_in/bmax)**(1.0D0/3.0D0)
! cs_out = DSQRT(amax*dpdx(xe)*xe/rho_in/3.0D0)
cs_out = DSQRT(amax*dpdx(xe)/3.0D0/(rho_in**2*bmax)**(1.0D0/3.0D0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	REAL*8 function dpdx(x)
	!$acc routine seq
	implicit none
	REAL*8 :: x
	dpdx = 8.0D0*x**4/DSQRT(x**2 + 1.0D0)
	end function

END SUBROUTINE
