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
IMPLICIT NONE

! Integer !
INTEGER, INTENT (IN) :: p_in
character(len=99) :: charac_p
INTEGER :: j,k,l, dummy
REAL*8 :: temp_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Updates to hydrodynamic variables !

IF (helmcheck_flag == 1) THEN
    write(charac_p,'(I)') count
    OPEN (UNIT = 123, FILE = './rho_eps_bfup'//trim(adjustl(charac_p))//'.dat', STATUS = 'REPLACE')
        DO l = 1, nz, 1 
            DO j = 1, nx, 1
                WRITE(123, *) prim(irho,j,1,l), epsilon(j,1,l)
            ENDDO
        ENDDO
ENDIF

IF (helmcheck_flag == 1) THEN
    write(charac_p,'(I)') count
    OPEN (UNIT = 123, FILE = './rho_temp_bfup'//trim(adjustl(charac_p))//'.dat', STATUS = 'REPLACE')
        DO l = 1, nz, 1 
            DO j = 1, nx, 1
                WRITE(123, *) prim(irho,j,1,l), temp2(j,1,l)
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
        DO l = 1, nz, 1 
            DO j = 1, nx, 1
                WRITE(123, *) prim(irho,j,1,l), temp2(j,1,l)
            ENDDO
        ENDDO
ENDIF

! Find pressure and speed of sound !
CALL FINDPRESSURE

IF (helmcheck_flag == 1) THEN
    OPEN (UNIT = 123, FILE = './rho_cs'//trim(adjustl(charac_p))//'.dat', STATUS = 'REPLACE')
    DO l = 1, nz, 1 
        DO j = 1, nx, 1
            WRITE(123, *) prim(irho,j,1,l), cs(j,1,l)
        ENDDO
    ENDDO
    CLOSE(123)

    OPEN (UNIT = 123, FILE = './rho_press'//trim(adjustl(charac_p))//'.dat', STATUS = 'REPLACE')
    DO l = 1, nz, 1 
        DO j = 1, nx, 1
            WRITE(123, *) prim(irho,j,1,l), prim(itau,j,1,l)
        ENDDO
    ENDDO
    CLOSE(123)

    OPEN (UNIT = 123, FILE = './rho_eps'//trim(adjustl(charac_p))//'.dat', STATUS = 'REPLACE')
    DO l = 1, nz, 1 
        DO j = 1, nx, 1
            WRITE(123, *) prim(irho,j,1,l), epsilon(j,1,l)
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