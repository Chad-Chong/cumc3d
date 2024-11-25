!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine contains the essential for the calculation of      !
! pressure, sound speed and other thermodynamic quantities.          !
! The EOSEPSILON and EOSSOUNDSPEED are required in reconstruction    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDPRESSURE
USE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

! Integer !
INTEGER :: j, k, l
! Real !
REAL*8 :: xe

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
