!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This subroutine prepare the data necessary for constructing
! the flux for the spatial discretization.
! It takes input/output of the U array and the 
! mode p which decides whether or not to update
! the gravitational potential
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE UPDATE (p_in)
USE DEFINITION
USE CUSTOM_DEF
USE HELMEOS_MODULE
USE TURB_MODULE
IMPLICIT NONE

! Integer !
INTEGER, INTENT (IN) :: p_in
character(len=99) :: charac_p
INTEGER :: j,k, dummy
REAL*8 :: temp_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Updates to hydrodynamic variables !

IF (helmcheck_flag == 1) THEN
write(charac_p,'(I)') count
OPEN (UNIT = 123, FILE = './rho_eps_bfup'//trim(adjustl(charac_p))//'.dat', STATUS = 'REPLACE')
    DO k = 1, ny, 1 
        DO j = 1, nx, 1
            WRITE(123, *) prim(irho,j,k,1), epsilon(j,k,1)
        ENDDO
    ENDDO
ENDIF

! Update Abar and Zbar !
IF(xisotran_flag == 1) CALL FIND_AZBAR()

! Find Temperature if using Helmholtz EOS !
if(helmeos_flag == 1) CALL FINDHELMTEMP()

IF (helmcheck_flag == 1) THEN
write(charac_p,'(I)') count
OPEN (UNIT = 123, FILE = './rho_temp'//trim(adjustl(charac_p))//'.dat', STATUS = 'REPLACE')
    DO k = 1, ny, 1 
        DO j = 1, nx, 1
            WRITE(123, *) prim(irho,j,k,1), temp2(j,k,1)
        ENDDO
    ENDDO
ENDIF

! Find pressure and speed of sound !
CALL FINDPRESSURE

IF (helmcheck_flag == 1) THEN
    OPEN (UNIT = 123, FILE = './rho_cs'//trim(adjustl(charac_p))//'.dat', STATUS = 'REPLACE')
    DO k = 1, ny, 1 
        DO j = 1, nx, 1
            WRITE(123, *) prim(irho,j,k,1), cs(j,k,1)
        ENDDO
    ENDDO
    CLOSE(123)

    OPEN (UNIT = 123, FILE = './rho_press'//trim(adjustl(charac_p))//'.dat', STATUS = 'REPLACE')
    DO k = 1, ny, 1 
        DO j = 1, nx, 1
            WRITE(123, *) prim(irho,j,k,1), prim(itau,j,k,1)
        ENDDO
    ENDDO
    CLOSE(123)

    OPEN (UNIT = 123, FILE = './rho_eps'//trim(adjustl(charac_p))//'.dat', STATUS = 'REPLACE')
    DO k = 1, ny, 1 
        DO j = 1, nx, 1
            WRITE(123, *) prim(irho,j,k,1), epsilon(j,k,1)
        ENDDO
    ENDDO
    CLOSE(123)
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom updates !

! Do custom updates !
CALL CUSTOM_UPDATE (p_in)

count = count + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE