!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine contains the essential for the calculation of      !
! pressure, sound speed and other thermodynamic quantities.          !
! The EOSEPSILON and EOSSOUNDSPEED are required in reconstruction    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDPRESSURE
USE nuceos_module
USE CUSTOM_DEF
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l
! Real !
REAL*8 :: rho_in,temp_in,ye_in
REAL*8 :: p_out,cs2_out

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

CALL FINDNUCTEMP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(rho_in,temp_in,ye_in,p_out,cs2_out)
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
DO l = -2, nz + 3
	DO k = -2, ny + 3
		DO j = -2, nx + 3
			rho_in = prim(irho,j,k,l)/rhocgs2code
			temp_in = temperature(j,k,l)
			ye_in = prim(iye,j,k,l)
			CALL nuc_eos_custom(rho_in,temp_in,ye_in,p_out,cs2_out)
            prim(itau,j,k,l) = p_out*taucgs2code
            cs(j,k,l) = DSQRT(cs2_out)/clight
   		 END DO
  	END DO
END DO
!$ACC END PARALLEL
!$OMP END PARALLEL DO

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'findpressure = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDNUCTEMP
USE nuceos_module
USE CUSTOM_DEF
IMPLICIT NONE

INTEGER :: i, j, k, l
INTEGER :: keyerr

! dummies for EOS call
REAL*8 :: rho_in, eps_in, temp_out, eps_min, eps_max

! Check timing with or without openmp
#ifdef DEBUG
INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
WRITE(*,*) 'start findnuctemp'
CALL system_clock(count_rate=cr)
rate = REAL(cr)
CALL system_clock(time_start)
#endif

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(rho_in, eps_min, eps_max, eps_in, temp_out, keyerr)
DO l = 1, nz
   DO k = 1, ny
      DO j = 1, nx
         IF(prim(irho,j,k,l) > prim_a(irho)) THEN
            rho_in = (prim(irho,j,k,l)/rhocgs2code)

            CALL nuc_eos_one(rho_in, eos_tempmin, prim(iye,j,k,l),eps_min,2)
            CALL nuc_eos_one(rho_in, eos_tempmax, prim(iye,j,k,l),eps_max,2)
            
            eps_min = (eps_min*energycgs2code)
            eps_max = (eps_max*energycgs2code)
            IF(epsilon(j,k,l) >= eps_min .and. epsilon(j,k,l) <= eps_max) THEN
               eps_in = (epsilon(j,k,l)/energycgs2code)
               temp_out = temperature(j,k,l)
               CALL nuc_eos_findtemp(rho_in,temp_out,prim(iye,j,k,l),eps_in,keyerr,1.0D-11)
               temperature(j,k,l) = MAX(temp_out, eos_tempmin)
            ELSEIF(epsilon(j,k,l) < eps_min) THEN
               epsilon(j,k,l) = eps_min
               temperature(j,k,l) = eos_tempmin
            ELSEIF(epsilon(j,k,l) > eps_max) THEN
               WRITE(*,*) "NUC temp ceiling", epsilon(j,k,l), eps_max
               epsilon(j,k,l) = eps_max
               temperature(j,k,l) = eos_tempmax
            ENDIF
         ELSE
            prim(irho,j,k,l) = prim_a(irho)
            temperature(j,k,l) = temp_a
            epsilon(j,k,l) = eps_a
         ENDIF
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO
CALL BOUNDARY1D_NM (temperature, even, even, even, even, even, even)
CALL BOUNDARY1D_NM (epsilon, even, even, even, even, even, even)

#ifdef DEBUG
CALL system_clock(time_end)
WRITE(*,*) 'findnuctemp = ', REAL(time_end - time_start) / rate
#endif

END SUBROUTINE findnuctemp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EOSEPSILON_NM (rho_in, p_in, eps_out)
!$acc routine seq
USE DEFINITION
IMPLICIT NONE

REAL*8, INTENT (IN) :: rho_in, p_in
REAL*8, INTENT (OUT) :: eps_out

eps_out = 0.0D0

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EOSSOUNDSPEED(p_in, rho_in, eps_in, cs_out)
!$acc routine seq
USE DEFINITION
IMPLICIT NONE

REAL*8, INTENT (IN) :: rho_in, p_in, eps_in
REAL*8, INTENT (OUT) :: cs_out

cs_out = 0.0D0

END SUBROUTINE
