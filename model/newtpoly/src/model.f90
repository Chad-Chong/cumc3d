!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model                                                           !
! Simulation for rotating magnetised white dwarf AIC                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

INTEGER :: j, k, l
REAL*8 :: xi, xl

ggas = 1.0D0 + 1.0D0/npoly
rho_floor = rhocen*1.0D-6
eps_floor = kpoly*rho_floor**(ggas - 1.0)/(ggas - 1.0)

prim(:,:,:,:) = 0.0d0
epsilon(:,:,:) = 0.0d0

DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      xl = xF(j)
      xi = xl/alpha
      IF (xi < pi*0.9D0) THEN
        prim(irho,j,k,l) = rhocen*DSIN(xi)/xi
        jNS = j
      ELSE
        prim(irho,j,k,l) = prim(irho,jNS,k,l)*(xF(jNS)/xl)**3.1d0
      END IF
      prim(itau,j,k,l) = kpoly*prim(irho,j,k,l)**2.0D0
      epsilon(j,k,l) = prim(itau,j,k,l)/prim(irho,j,k,l)               ! p=rho*eps*(gamma-1)
      dphidx(j,k,l) = -(ggas*prim(itau,j,k,l)/prim(irho,j,k,l)**2)*&   ! dphidx = -(1/rho)*(dp/drho)*(drho/dr) 
                      rhocen*(DCOS(xi)/xl - alpha*DSIN(xi)/xl**2)
    END DO
  END DO
END DO
PRINT *, "Finished generating star"

prim_a(:) = 0.0D0
prim_a(irho) = rho_floor
prim_a(itau) = prim(itau,nx,1,1)
eps_a = epsilon(nx,1,1)

total_time = 120.0D0
output_profiletime = total_time/10.0D0


END SUBROUTINE
