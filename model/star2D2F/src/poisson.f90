!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the poisson equation coefficient for the relaxation method !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_poisson
USE CUSTOM_DEF 
IMPLICIT NONE

INTEGER :: i, j, k
DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nz
      CALL poisson_coef(x(i-1), x(i), x(i+1), y(j-1), y(j), y(j+1), z(k-1), z(k), z(k+1), &
              ajp1(i), ajm1(i), bkp1(i,j), bkm1(i,j), clp1(i,j,k), clm1(i,j,k), epsc(i,j,k))
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE get_poisson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get left hand side of the discrete poisson equation ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE poisson_coef(xm1, xc, xp1, ym1, yc, yp1, zm1, zc, zp1, &
    alphajp1, alphajm1, betakp1, betakm1, gammalp1, gammalm1, epsc)
IMPLICIT NONE

REAL*8, INTENT(IN) :: xm1, xc, xp1, ym1, yc, yp1, zm1, zc, zp1
REAL*8, INTENT(OUT) :: epsc, alphajp1, alphajm1, betakp1, betakm1, gammalp1, gammalm1

epsc = 2.0d0*(xp1+xm1-3.0d0*xc)/xc/(xp1-xc)/(xc-xm1) &
+ (yp1+ym1-2.0d0*yc-2.0d0*DTAN(yc))/xc**2/DTAN(yc)/(yp1-yc)/(yc-ym1) &
- 2.0d0/xc**2/DSIN(yc)**2/(zp1-zc)/(zc-zm1)

alphajp1 = 2.0d0*(2.0d0*xc-xm1)/(xp1-xm1)/(xp1-xc)/xc
alphajm1 = 2.0d0*(2.0d0*xc-xp1)/(xp1-xm1)/(xc-xm1)/xc

betakp1 = (2.0d0*DTAN(yc) + yc - ym1)/xc**2/DTAN(yc)/(yp1 - yc)/(yp1 - ym1)
betakm1 = (2.0d0*DTAN(yc) + yc - yp1)/xc**2/DTAN(yc)/(yp1 - ym1)/(yc - ym1)

gammalp1 = 2.0d0/xc**2/DSIN(yc)**2/(zp1-zm1)/(zp1-zc)
gammalm1 = 2.0d0/xc**2/DSIN(yc)**2/(zp1-zm1)/(zc-zm1) 

END SUBROUTINE poisson_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve the poisson equation by RBSOR !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE POISSON_SOLVER
USE CUSTOM_DEF
IMPLICIT NONE

INTEGER :: j, k, l, n
REAL*8 :: abserror, rhs
REAL*8 :: rho_in, factor, diff

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
      rho_grav(j,k,l) = prim(irho,j,k,iNM) + prim(irho,j,k,iDM)
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
        diff = rho_grav(j,k,l) - primNM_a(irho)
		factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
		rho_in = factor*rho_grav(j,k,l)
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
      diff = rho_grav(j,k,l) - primNM_a(irho)
      factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
      rho_in = factor*rho_grav(j,k,l)
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

! Boundary conditions !
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
!$OMP END PARALLEL

IF(abserror <= tolerance) EXIT 

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stop condition !
IF(n == relax_max) THEN
  WRITE (*,*) n, relax_max
  STOP 'Convergence error in poisson solver'
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (gr_potential == 1) THEN
  STOP 'GR potential not implemented for DM'
  CALL GRPOTENTIAL
END IF

END SUBROUTINE POISSON_SOLVER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GR modified potential     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GRPOTENTIAL
USE CUSTOM_DEF
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l

! Real !
REAL*8, DIMENSION(-2:nx+3) :: dx_total,phi_avg,rho_avg,tau_avg,eps_avg,vr_avg,mtov_avg,lorentz2,integrand,sum_temp
REAL*8 :: check_value, int_old, total_old, int_new, total_new
REAL*8, PARAMETER :: compact = 0.1D0
INTEGER, PARAMETER :: corr_max = 5

! Initialize !
mtov_avg = 0.0D0
dx_total(0) = 0.5D0*dx(1)
dx_total(1:nx+1) = 0.5D0*(dx(1:nx+1) + dx(2:nx+2))
phi_avg(1:nx+1) = SUM(phi(1:nx+1, 1:ny, 1), DIM=2) / ny
rho_avg(1:nx+1) = SUM(prim(irho, 1:nx+1, 1:ny, 1), DIM=2) / ny
tau_avg(1:nx+1) = SUM(prim(itau, 1:nx+1, 1:ny, 1), DIM=2) / ny
eps_avg(1:nx+1) = SUM(epsilon(1:nx+1, 1:ny, 1), DIM=2) / ny
vr_avg(1:nx+1)  = SUM(prim(ivx, 1:nx+1, 1:ny, 1), DIM=2) / ny

integrand(1:nx+1) = 4.0D0*pi*x(1:nx+1)**2*rho_avg(1:nx+1)*(1.0D0 + eps_avg(1:nx+1))

DO j = 1, nx + 1
  IF (j == 1) THEN
    int_old = 0.0D0
    total_old = 0.0D0
  ELSE
    check_value = 1.0D0 + vr_avg(j-1)**2 - 2.0D0*mtov_avg(j-1)/x(j-1)
    IF (check_value < 0.0D0) THEN
      check_value = 1.0D0 - compact
    END IF
    int_old = SQRT(check_value)
    total_old = int_old*integrand(j-1)
  END IF

  ! Update !
  mtov_avg(j) = mtov_avg(j-1) + dx_total(j-1)*total_old

  DO k = 1, corr_max
    check_value = 1.0D0 + vr_avg(j)**2 - 2.0D0*mtov_avg(j)/x(j)
    IF (check_value < 0.0D0) THEN
      check_value = 1.0D0 - compact
    END IF
    int_new = SQRT(check_value)
    total_new = int_new*integrand(j)
    mtov_avg(j) = mtov_avg(j-1) + 0.5D0*dx_total(j-1)*(total_old + total_new)
  END DO
END DO

! Initialise !
integrand(:) = 0.0D0

! Get the integrand !
DO j = 1, nx + 1
  lorentz2(j) = 1.0D0 + vr_avg(j)*vr_avg(j) - 2.0D0*mtov_avg(j)/x(j)
  IF (x(j) > 0.0D0) THEN
    check_value = 1.0D0 + eps_avg(j) + tau_avg(j)/x(j)
  ELSE
    check_value = 1.0D0
  END IF
  integrand(j) = ((mtov_avg(j)/4.0D0/pi) + x(j)**3*tau_avg(j))*check_value/(x(j)**2)/(lorentz2(j))
END DO

! Initialise !
phi_gr(:) = 0.0D0

! Solve for the TOV potential !
DO j = 2, nx + 1
	phi_gr(1) = phi_gr(1) + 0.5D0*(integrand(j) + integrand(j-1))*dx_total(j-1)
END DO
DO j = 2, nx + 1
	phi_gr(j) = phi_gr(j-1) - 0.5D0*(integrand(j) + integrand(j-1))*dx_total(j-1)
END DO

! Add the outer integral term !
phi_gr(:) = phi_gr(:) + log(x(nx+1)/(x(nx+1) - 2.0D0*mtov_avg(nx+1)))/(8.0D0*pi)
phi_gr(:) = - phi_avg(:) - 4.0D0*pi*phi_gr(:)

! Boundary conditions !
phi_gr(0) = phi_gr(1)
phi_gr(nx+1) = 0.0D0

END SUBROUTINE GRPOTENTIAL
