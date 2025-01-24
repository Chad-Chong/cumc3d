!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine contains the essential for the calculation of      !
! pressure, sound speed and other thermodynamic quantities.          !
! The EOSEPSILON and EOSSOUNDSPEED are required in reconstruction    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDPRESSURE
USE CUSTOM_DEF
IMPLICIT NONE

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
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(xe)
DO k = -2, ny + 3
	DO j = -2, nx + 3
		xe = (prim(irho,j,k,iNM)/bmaxNM)**(1.0D0/3.0D0)
		IF (xe > 1.0D-2) THEN
			prim(itau,j,k,iNM) = amaxNM*large_pressure(xe)
		ELSE
			prim(itau,j,k,iNM) = amaxNM*small_pressure(xe)
		END IF
		CALL EOSEPSILON(prim(irho,j,k,iNM), prim(itau,j,k,iNM), epsilon(j,k,iNM),iNM)
		CALL EOSSNDSPEED(prim(itau,j,k,iNM), prim(irho,j,k,iNM), epsilon(j,k,iNM), cs(j,k,iNM),iNM)
	END DO
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(xe)
DO k = -2, ny + 3
	DO j = -2, nx + 3
		xe = (prim(irho,j,k,iDM)/bmaxDM)**(1.0D0/3.0D0)
		IF (xe > 1.0D-2) THEN
			prim(itau,j,k,iDM) = amaxDM*large_pressure(xe)
		ELSE
			prim(itau,j,k,iDM) = amaxDM*small_pressure(xe)
		END IF
		CALL EOSEPSILON(prim(irho,j,k,iDM), prim(itau,j,k,iDM), epsilon(j,k,iDM),iDM)
		CALL EOSSNDSPEED(prim(itau,j,k,iDM), prim(irho,j,k,iDM), epsilon(j,k,iDM), cs(j,k,iDM),iDM)
	END DO
END DO
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
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EOSEPSILON (rho_in, p_in, eps_out, itype)
!$acc routine seq
USE CUSTOM_DEF
IMPLICIT NONE

REAL*8 :: xe
INTEGER, INTENT (IN) :: itype
REAL*8, INTENT (IN) :: rho_in, p_in
REAL*8, INTENT (OUT) :: eps_out

SELECT CASE(itype)
	CASE (iNM)
		xe = (rho_in/bmaxNM)**(1.0D0/3.0D0)
		IF(xe > 1.0D-2) THEN
			eps_out = amaxNM*large_energy(xe)/rho_in
		ELSE
			eps_out = amaxNM*small_energy(xe)/rho_in
		END IF
	CASE (iDM)
		xe = (rho_in/bmaxDM)**(1.0D0/3.0D0)
		IF(xe > 1.0D-2) THEN
			eps_out = amaxDM*large_energy(xe)/rho_in
		ELSE
			eps_out = amaxDM*small_energy(xe)/rho_in
		END IF
END SELECT
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

SUBROUTINE EOSSNDSPEED(p_in, rho_in, eps_in, cs_out, itype)
!$acc routine seq
USE CUSTOM_DEF
IMPLICIT NONE

REAL*8 :: xe
INTEGER, INTENT (IN) :: itype
REAL*8, INTENT (IN) :: p_in, rho_in, eps_in
REAL*8, INTENT (OUT) :: cs_out

SELECT CASE(itype)
	CASE (iNM)
		xe = (rho_in/bmaxNM)**(1.0D0/3.0D0)
		cs_out = DSQRT(amaxNM*dpdx(xe)/3.0D0/(rho_in**2*bmaxNM)**(1.0D0/3.0D0))
	CASE (iDM)
		xe = (rho_in/bmaxDM)**(1.0D0/3.0D0)
		cs_out = DSQRT(amaxDM*dpdx(xe)/3.0D0/(rho_in**2*bmaxDM)**(1.0D0/3.0D0))
END SELECT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

	REAL*8 function dpdx(x)
	!$acc routine seq
	implicit none
	REAL*8 :: x
	dpdx = 8.0D0*x**4/DSQRT(x**2 + 1.0D0)
	end function

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dummy Subroutines for other non-PPMC reconstruction, this model only use PPMC !
SUBROUTINE EOSSOUNDSPEED(p_in, rho_in, eps_in, cs_out)
IMPLICIT NONE
REAL*8, INTENT (IN) :: p_in, rho_in, eps_in
REAL*8, INTENT (OUT) :: cs_out
cs_out = 0.0D0
END SUBROUTINE

SUBROUTINE EOSEPSILON_NM (rho_in, p_in, eps_out)
IMPLICIT NONE
REAL*8, INTENT (IN) :: rho_in, p_in
REAL*8, INTENT (OUT) :: eps_out
eps_out = 0.0D0
END SUBROUTINE
