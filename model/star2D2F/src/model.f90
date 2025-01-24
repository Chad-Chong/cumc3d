!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preperation !
! Poisson interpolation coefficient !
call get_poisson

! Read and assign density !
OPEN(UNIT=970, FILE = './profile/hydro_rhoNM.dat', ACTION='READ')
DO j = 1, nx
  READ(970,*) prim(irho,j,1,iNM)
END DO
CLOSE(970)
PRINT *, "Finished reading rho"


OPEN(UNIT=970, FILE = './profile/hydro_rhoDM.dat', ACTION='READ')
DO j = 1, nx
  READ(970,*) prim(irho,j,1,iDM)
END DO
CLOSE(970)

DO k = 1, ny
  DO l = iNM, iDM
    prim(irho,:,k,l) = prim(irho,:,1,l)
  END DO
END DO

prim(irho,:,:,:) = prim(irho,:,:,:)*rhocgs2code
PRINT *, "Finished reading rho"

! Assign velocity !
prim(ivx:ivz,:,:,:) = 0.0d0
prim(ibx:ibz,:,:,:) = 0.0d0

DO j = 1, nx
  DO k = 1, ny
    DO l = iNM, iDM
      prim(ivz,j,k,l) = 2.0D0*x(j)*DSIN(y(k))/tcgs2code
    END DO
  END DO
END DO


CALL FINDPRESSURE
PRINT *, "Finished calculating pressure and epsilon"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
primNM_a(:) = 0.0D0
primNM_a(irho) = rhoatm
primNM_a(itau) = prim(itau,nx,1,iNM)
epsNM_a = epsilon(nx,1,iDM)

primDM_a(:) = 0.0D0
primDM_a(irho) = rhoatm
primDM_a(itau) = prim(itau,nx,1,iNM)
epsDM_a = epsilon(nx,1,iDM)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign floor density !
DO j = 1, nx
  DO k = 1, ny
    DO l = iNM, iDM
      IF(prim(irho,j,k,l) < rhoatm) THEN
        prim(irho,j,k,l) = rhoatm
        prim(ivx:ivz,j,k,l) = 0.0d0
        epsilon(j,k,l) = epsNM_a
      END IF
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set output profile interval !
total_time = 0.5d0*tcgs2code
output_profiletime = total_time/50.0d0

END SUBROUTINE
