SUBROUTINE FINDPRESSURE
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

INTEGER :: j, k, l

DO l = -2, nz + 3
  DO k = -2, ny + 3
    DO j = -2, nx + 3
      prim(itau,j,k,l) = kpoly*prim(irho,j,k,l)**2.0D0
      epsilon(j,k,l) = prim(itau,j,k,l)/prim(irho,j,k,l)
      cs(j,k,l) = DSQRT(ggas*prim(itau,j,k,l)/prim(irho,j,k,l))
    END DO
  END DO
END DO

END SUBROUTINE

SUBROUTINE EOSSOUNDSPEED(p_in, rho_in, eps_in, cs_out)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

REAL*8, INTENT (IN) :: p_in, rho_in, eps_in
REAL*8, INTENT (OUT) :: cs_out

cs_out = DSQRT(ggas*p_in/rho_in)

END SUBROUTINE

SUBROUTINE EOSEPSILON_NM(p_in, rho_in, eps_out)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

REAL*8, INTENT (IN) :: p_in, rho_in
REAL*8, INTENT (OUT) :: eps_out

eps_out = p_in/rho_in

END SUBROUTINE