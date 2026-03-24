!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model                                                           !
! Simulation for rotating magnetised white dwarf type Ia                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE CUSTOM_DEF
USE MHD_MODULE
USE IEEE_ARITHMETIC
USE Ecaptable_module
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l, m, n, flag_eostable, flag_notfindtemp, found_atmo
REAL*8 :: dummy
! Magnetic field !
REAL*8 :: maxdb
REAL*8 :: div_b

REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: a_phi

! For Atmosphere
REAL*8 :: diff, factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flag Checking !
IF (axissym_flag == 1 .and. coordinate_flag == 1) THEN

  IF (ny > 1) THEN
    WRITE(*,*) 'For axis symmetry there can only be one grid in the azimuth direction'
    STOP
  ENDIF

  IF (n_dim /= 3) THEN
    WRITE(*,*) 'For axis symmetry in cylindrical coordinates n_dim needs to be 3'
    STOP
  ENDIF
ELSE
  WRITE(*,*) 'axissym_flag only imposes axis symmetry for cylindrical coordinates'
  STOP  
ENDIF

IF (turb_flag == 1 .and. coordinate_flag /= 1) THEN
  WRITE(*,*) 'SGS Turbulence is only implemented for cylindrical coordinates at the moment'
  STOP 
ENDIF

IF (helmeos_flag == 1 .and. xisotran_flag /= 1) THEN
  WRITE(*,*) 'Helm eos must be initiated with isotope tracking'
  STOP
ENDIF

IF (levelset_flag == 1 .and. coordinate_flag /= 1) THEN
  WRITE(*,*) 'The level set method (2D) for flame tracking is only implemented for cylindrical coordinates at the moment'
  STOP 
ENDIF

IF (flame_flag == 1 .and. coordinate_flag /= 1) THEN
  WRITE(*,*) 'The level set method (2D) for flame tracking is only implemented for cylindrical coordinates at the moment'
  STOP 
ENDIF

IF (flame_flag == 1 .and. xisotran_flag /= 1) THEN
  WRITE(*,*) 'The flame must be initiated with the helm eos with isotope tracking'
  STOP
ENDIF

IF (coordinate_flag == 0) THEN
  WRITE(*,*) 'This code is not written for Cartesian coordinates'
  STOP
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preperation !

! Allocate
Allocate(a_phi(-3:nx+3,-3:ny+3,-3:nz+3))

! Poisson interpolation coefficient !
IF (gravity_flag == 1) THEN
  call get_poisson
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (restart_flag == 0) THEN
  ! Read and assign density !
  OPEN(UNIT=970, FILE = './profile/hydro_rho.dat', ACTION='READ')
  IF (coordinate_flag == 2) THEN
    READ(970,*) ((prim(irho,j,k,1), j = 1, nx), k = 1, ny)
  ELSEIF (coordinate_flag == 1) THEN
    READ(970,*) ((prim(irho,j,1,l), j = 1, nx), l = 1, nz)
  ENDIF
  CLOSE(970)
  prim(irho,:,:,:) = prim(irho,:,:,:)*rhocgs2code
  atmosphere = atmospheric*maxval(prim(irho,:,:,:))
  DO l = 1, nz
    DO j = 1, nx
      IF (prim(irho,j,1,l) < atmosphere) THEN
        prim(irho,j,1,l) = atmosphere
      ENDIF
    ENDDO
  ENDDO
  PRINT *, "Finished reading rho"
  WRITE(*,*) 'Maximum rho is', atmosphere/atmospheric/rhocgs2code

  ! Assign velocity !
  IF (coordinate_flag == 2) THEN
    prim(ivx,:,:,1) = 0.0d0
    prim(ivy,:,:,1) = 0.0d0
  ELSEIF (coordinate_flag == 1) THEN
    prim(ivx,:,1,:) = 0.0d0
    prim(ivz,:,1,:) = 0.0d0
  ENDIF

  IF (rotate_flag == 1) THEN
    ! Read and assign azimuth direction velocity !
    OPEN(UNIT=970, FILE = './profile/hydro_vphi.dat', ACTION='READ')
    IF (coordinate_flag == 2) THEN
      READ(970,*) ((prim(ivz,j,k,1), j = 1, nx), k = 1, ny)
      prim(ivz,:,:,:) = prim(ivz,:,:,:)*lencgs2code/tcgs2code
    ELSEIF (coordinate_flag == 1) THEN
      READ(970,*) ((prim(ivy,j,1,l), j = 1, nx), l = 1, nz)
      prim(ivy,:,:,:) = prim(ivy,:,:,:)*lencgs2code/tcgs2code
    ENDIF
    CLOSE(970)
    PRINT *, "Finished reading vphi"
  ELSE
    IF (coordinate_flag == 2) THEN
      prim(ivz,:,:,:) = 0.0D0
    ELSEIF (coordinate_flag == 1) THEN
      prim(ivy,:,:,:) = 0.0D0
    ENDIF
    PRINT *, "No rotation in this run"
  ENDIF

  IF (mhd_flag == 1) THEN
    ! Read for magnetic vector potential !
    OPEN(UNIT=970, FILE = './profile/hydro_Aphi.dat', ACTION='READ')
    IF (coordinate_flag == 2) THEN
      READ(970,*) ((a_phi(j,k,1), j = 0, nx), k = 0, ny)
    ELSEIF (coordinate_flag == 1) THEN
      READ(970,*) ((a_phi(j,1,l), j = 0, nx), l = 0, nz)
      a_phi(j,:,l) = a_phi(j,1,l)
    ENDIF

    CLOSE(970)

    PRINT *, "Finished reading Aphi"

    ! In the direction of symmetry, cell centered is the same as face centered
    OPEN(UNIT=970, FILE = './profile/hydro_bphi.dat', ACTION='READ')
    IF (coordinate_flag == 2) THEN
      READ(970,*) ((prim(ibz,j,k,1), j = 1, nx), k = 1, ny)
      prim(ibz,:,:,0) = prim(ibz,:,:,1)
    ELSEIF (coordinate_flag == 1) THEN
      READ(970,*) ((prim(iby,j,1,l), j = 1, nx), l = 1, nz)
      DO k = -2, ny+3, 1 
        prim(iby,:,k,:) = prim(iby,:,1,:)
      ENDDO
      prim(iby,:,:,:) = prim(iby,:,:,:)*gauss2code
    ENDIF

    CLOSE(970)

    PRINT *, "Finished reading Bphi"

    ! WRITE(*,*) ABS(xF(0)/xF(1))

    ! Coordinate here are in code unit but aphi is in gaussian unit !
    ! Unit conversion below !
    IF (coordinate_flag == 2) THEN
      DO l = 1, nz
        DO k = 1, ny
          DO j = 0, nx
            prim(ibx,j,k,l) = (sinf(k)*a_phi(j,k,l) - sinf(k-1)*a_phi(j,k-1,l))/(xF(j)*dcose(k)+small_num)
          END DO
        END DO
      END DO

      DO l = 1, nz
        DO k = 0, ny
          DO j = 1, nx
            prim(iby,j,k,l) = - (xF(j)*a_phi(j,k,l) - xF(j-1)*a_phi(j-1,k,l))/(x(j)*dx(j))
          END DO
        END DO
      END DO

      prim(ibx:ibz,:,:,:) = prim(ibx:iby,:,:,:)*gauss2code*lencgs2code  ! length conversion for curl !
    ELSEIF(coordinate_flag == 1) THEN
      DO l = 0, nz
        DO k = 0, ny
          DO j = 0, nx
            prim(ibx,j,k,l) = - (a_phi(j,k,l) - a_phi(j,k,l-1))/(dz(l))

            IF (j==0) THEN
              IF ( ABS(xF(0)/xF(1)) < 1e-10 ) THEN
                prim(ibx,j,k,l) = 0.0D0
              ENDIF
            ENDIF

          END DO
        END DO
      END DO

      DO l = 0, nz
        DO k = 0, ny
          DO j = 0, nx
            prim(ibz,j,k,l) = (xF(j)*a_phi(j,k,l) - xF(j-1)*a_phi(j-1,k,l))/(x(j)*dx(j))
          END DO
        END DO
      END DO

      prim(ibx,:,:,:) = prim(ibx,:,:,:)*gauss2code*lencgs2code  ! length conversion for curl !
      prim(ibz,:,:,:) = prim(ibz,:,:,:)*gauss2code*lencgs2code  ! length conversion for curl !
    ENDIF

    PRINT *, "Finished calculating magnetic field"
  ELSE
    prim(ibx:ibz,:,:,:) = 0.0D0
    PRINT *, "No magnetic field in this run"
  ENDIF

  ! Deallocate
  Deallocate(a_phi)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define the Kronecker Delta

do m = 1, 3, 1
   do n = 1, 3, 1
      if (m==n) then
         eye(m,n) = 1.0D0
      else
         eye(m,n) = 0.0D0
      endif
   enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Chemical Composition !

IF (xisotran_flag == 1) THEN

  CALL read_helm_table()
  CALL initialize_network()
  WRITE(*,*), 'Finished reading helmeos, flame and nse table'

  IF (restart_flag == 0) THEN
    DO i = ihe4, ini56
      m = i+1-ihe4
      OPEN(UNIT=970, FILE = trim('./profile/Helm_X'//ionam(m))//'.dat', ACTION='READ')
      IF (coordinate_flag == 2) THEN
        READ(970,*) ((prim(i,j,k,1), j = 1, nx), k = 1, ny)
      ELSEIF (coordinate_flag == 1) THEN
        READ(970,*) ((prim(i,j,1,l), j = 1, nx), l = 1, nz)
      ENDIF
      CLOSE(970)
    ENDDO
  
    PRINT *, "Finished reading chemical composition"

  ENDIF

  IF (helmeos_flag == 1) THEN
    CALL FIND_AZBAR()
  ENDIF
  

ENDIF

IF (helmeos_flag == 1) THEN

  IF (restart_flag == 0) THEN

    OPEN(UNIT=970, FILE = trim('./profile/Helm_temp.dat'), ACTION='READ')
    IF (coordinate_flag == 2) THEN
      READ(970,*) ((temp2(j,k,1), j = 1, nx), k = 1, ny)
    ELSEIF (coordinate_flag == 1) THEN
      READ(970,*) ((temp2(j,1,l), j = 1, nx), l = 1, nz)
    ENDIF
    CLOSE(970)
    temp2 = temp2/1.0D9
    temp2_old = temp2
    PRINT *, "Finished reading temperature for Helmholtz EOS"

    OPEN(UNIT=970, FILE = trim('./profile/Helm_Ye.dat'), ACTION='READ')
    IF (coordinate_flag == 2) THEN
      READ(970,*) ((prim(iye2,j,k,1), j = 1, nx), k = 1, ny)
    ELSEIF (coordinate_flag == 1) THEN
      READ(970,*) ((prim(iye2,j,1,l), j = 1, nx), l = 1, nz)
    ENDIF
    CLOSE(970)
    PRINT *, "Finished reading Ye for Helmholtz EOS (electron capture)"

  ENDIF


ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (nuspec_flag == 1) THEN
  CALL read_nutable
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (restart_flag == 0) THEN

  CALL FINDPRESSURE

  IF (helmeos_flag == 1) THEN
  
    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx

          CALL HELM_EOSEPSILON(prim(irho,j,k,l), temp2(j,k,l), abar2(j,k,l), zbar2(j,k,l), prim(iye2,j,k,l), epsilon(j,k,l))
          
        END DO
      END DO
    END DO
    CALL BOUNDARY1D_NM(epsilon, even, even, even, even, even, even)

    ! OPEN (UNIT = 125, FILE = './eps.dat', STATUS = 'REPLACE')
    !       DO k = 1, nz, 1
    !          DO i = 1, nx, 1
    !             WRITE(125,*) epsilon(i,1,k)
    !          ENDDO
    !       ENDDO
    ! CLOSE(125)

    IF (helmeos_flag == 0 .and. coordinate_flag == 2) THEN
      DO j = 1, nx
        DO k = 1, ny
          prim(irho,j,k,1) = max(prim(irho,j,k,1), atmosphere)
          CALL EOSEPSILON_NM(prim(irho,j,k,1), prim(itau,j,k,1), epsilon(j,k,1))
        END DO
      END DO
    ELSEIF (helmeos_flag == 0 .and. coordinate_flag == 1) THEN 
      DO j = 1, nx
        DO l = 1, nz
          prim(irho,j,k,1) = max(prim(irho,j,1,l), atmosphere)
          CALL EOSEPSILON_NM(prim(irho,j,1,l), prim(itau,j,1,l), epsilon(j,1,l))
        END DO
      END DO
    ENDIF

    PRINT *, 'Finished calculating pressure, sound speed and (epsilon)'

  ENDIF

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sub-grid Turbulence !

IF (turb_flag == 1) THEN
  IF (ABS(dx(1) - dz(1)) > 1e-3*lencgs2code) THEN
    WRITE(*,*) 'dx(1)=', dx(1), 'dz(1)=', dz(1)
    WRITE(*,*) 'Grid is not uniform for cylindrical sub-grid scale turbulence'
    STOP
  ENDIF
  IF (restart_flag == 0) THEN
    CALL GetTurb
  ENDIF
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (restart_flag == 0) THEN
  ! Set atmospheric primitive variables!
  prim_a(:) = 0.0D0
  ! prim_a(irho) = atmospheric*maxval(prim(irho,:,:,:))
  IF (turb_flag == 1) THEN
    prim_a(iturbq) = turb_q_a
  ENDIF

  ! Ensure epsilon do not oscillate. Somehow the code dislikes osciallating energy very much

  IF (helmeos_flag == 1) THEN
    ! prim_a(ihe4) = xiso_ahe4
    prim_a(ic12) = xiso_ac12
    prim_a(io16) = xiso_ao16
    CALL PRIVATE_HELMEOS_AZBAR(prim_a(ihe4:ini56), abar2_a, zbar2_a, prim_a(iye2))

    eps_a = MINVAL(epsilon(1:nx, 1, 1:nz))

    found_atmo = 0

    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx
          IF (helmeos_flag == 1) THEN
            IF (epsilon(j,k,l) == eps_a .and. found_atmo == 0) THEN

              eps_a = epsilon(j,k,l)
              prim_a(irho) = prim(irho,j,k,l)

              
              WRITE(*,*) 'Atmosphere epsilon is', eps_a

              CALL private_invert_helm_ed(epsilon(j,k,l), &
                                prim_a(irho), abar2_a, &
                                zbar2_a, prim_a(iye2), &
                                                temp2_old(j,k,l), temp2_a, flag_notfindtemp)
              WRITE(*,*) 'Found temp_a is', temp2_a
              CALL HELM_EOSPRESSURE(prim_a(irho), temp2_a, abar2_a, zbar2_a, prim_a(iye2), prim_a(itau), dummy, dummy, flag_eostable)

              found_atmo = 1
            
            ENDIF
          ELSE

            prim_a(itau) = prim(itau,nx,1,1)
            eps_a = epsilon(nx,1,1)

            EXIT

          ENDIF

        END DO
      END DO
    END DO

    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx

          diff = prim(irho,j,k,l) - prim_a(irho)
          factor = MAX(SIGN(1.0D0, diff), 0.0D0)
          prim(imin:ibx-1,j,k,l) = factor*prim(imin:ibx-1,j,k,l) + (1.0D0 - factor)*prim_a(:)
          IF (helmeos_flag == 1) THEN
            temp2(j,k,l) = factor*temp2(j,k,l) + (1.0D0 - factor)*temp2_a
            abar2(j,k,l) = factor*abar2(j,k,l) + (1.0D0 - factor)*abar2_a
            zbar2(j,k,l) = factor*zbar2(j,k,l) + (1.0D0 - factor)*zbar2_a
            epsilon(j,k,l) = factor*epsilon(j,k,l) + (1.0D0 - factor)*eps_a
          ENDIF
          
        END DO
      END DO
    END DO
  ENDIF

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (restart_flag == 0) THEN
  ! Assign floor variables !
  CALL CUSTOM_CHECKRHO
  WRITE(*,*) 'Finished atmosphere handling'
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (restart_flag == 0) THEN

  IF (xisotran_flag == 1 .and. helmeos_flag == 1 .and. levelset_flag == 1 .and. flame_flag == 1) THEN
    CALL GetFlame
    WRITE(*,*) 'Finished assigning deflagration level set'
    CALL FLAME_INI()
    WRITE(*,*) 'Finished injecting flame energy'
    CALL FIND_AZBAR()	
    CALL FINDHELMTEMP()
    WRITE(*,*) 'Finished initializing flame'

    ! Assign floor variables !
    CALL CUSTOM_CHECKRHO
  ENDIF

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (restart_flag == 1) THEN
  CALL READ_HDF5
  CALL BOUNDARY
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find DIVB !
CALL FIND_DIVB
WRITE (*,*)
WRITE (*,*) 'Maximum initial divergence B', maxDivB
WRITE (*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set output profile interval !
total_time = 1.5D0*tcgs2code + global_time ! cgs; global time is in code unit already
output_profiletime = (total_time-global_time)/100.0d0
WRITE(*,*) 'The code outputs every', output_profiletime, 'code unit'

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE READ_HDF5
USE HDF5
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

INTEGER :: error
INTEGER(HID_T) :: file_id,dset_id,mem_id,dspace_id
INTEGER(HSIZE_T) :: prim_dims(4), dims(3), bfield_dims(4), prim_adim(1), start(4), start_s(1),stride(4), block(4), stride_s(1), block_s(1), ddims(4), tmp(4), start3(3), stride3(3), block3(3)
character(len=99) :: restart_file, index_time

WRITE(*,*) 'Enter the the restart index for the code.'
READ(*,*) start_index

write(index_time,'(I)') start_index

! WRITE(*,*) 'Restarting the code from output index', index_time

restart_file = './outfile/rkiter-'// trim(adjustl(index_time)) //'-nm.hdf5'

call h5open_f(error)
call h5fopen_f(restart_file,H5F_ACC_RDONLY_F,file_id,error)

! Read primitive
start = 0
stride = 1
block = 1
prim_dims(1) = iye2-imin+1
prim_dims(2) = nx+2
prim_dims(3) = ny+2
prim_dims(4) = nz+2

call h5dopen_f(file_id,"primitive",dset_id,error)
call h5dget_space_f(dset_id, dspace_id, error) 
call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, prim_dims, error, stride, block)
call h5screate_simple_f(4, prim_dims, mem_id, error)
call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, prim(imin:iye2,0:nx+1,0:ny+1,0:nz+1), prim_dims, error, mem_id, dspace_id)
CALL h5dclose_f(dset_id, error)
call h5sclose_f(mem_id, error)
call h5sclose_f(dspace_id, error)

start_s(1) = 0
prim_adim(1) = iye2-imin+1
stride_s = 1
block_s = 1

call h5dopen_f(file_id,"prim_a",dset_id,error)
call h5dget_space_f(dset_id, dspace_id, error) 
call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start_s, prim_adim, error, stride_s, block_s)
call h5screate_simple_f(1, prim_adim, mem_id, error)
call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, prim_a(imin:iye2), prim_adim, error, mem_id, dspace_id)
CALL h5dclose_f(dset_id, error)
call h5sclose_f(mem_id, error)
call h5sclose_f(dspace_id, error)

prim_adim(1) = 1
call h5dopen_f(file_id,"time",dset_id,error) !For atmosphere
call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,global_time,prim_adim,error)
CALL h5dclose_f(dset_id, error)

IF (error == 0) THEN
  PRINT *, "Finished reading primitive (hydro+composition) and time from hdf5 file"
ELSE
  PRINT *, "Error reading primitive (hydro+composition) and time from hdf5 file"
  STOP
END IF

! Read epsilon
dims(1) = nx+2
dims(2) = ny+2
dims(3) = nz+2
error = 0
call h5dopen_f(file_id,"epsilon",dset_id,error)
call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,epsilon(0:nx+1,0:ny+1,0:nz+1),dims,error)
CALL h5dclose_f(dset_id, error)

prim_adim(1) = 1

call h5dopen_f(file_id,"eps_a",dset_id,error)
call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,eps_a,prim_adim,error)
CALL h5dclose_f(dset_id, error)

IF (error == 0) THEN
  PRINT *, "Finished reading epsilon from hdf5 file"
ELSE
  PRINT *, "Error reading epsilon hdf5 file"
  STOP
END IF

! Read temperature
dims(1) = nx+2
dims(2) = ny+2
dims(3) = nz+2
error = 0
call h5dopen_f(file_id,"temp",dset_id,error)
call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,temp2(0:nx+1,0:ny+1,0:nz+1),dims,error)
CALL h5dclose_f(dset_id, error)

temp2_old = temp2
prim_adim(1) = 1

call h5dopen_f(file_id,"temp2_a",dset_id,error)
call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,temp2_a,prim_adim,error)
CALL h5dclose_f(dset_id, error)

IF (error == 0) THEN
  PRINT *, "Finished reading temperature from hdf5 file"
ELSE
  PRINT *, "Error reading temperature hdf5 file"
  STOP
END IF

! Read magnetic field
bfield_dims(1) = (ibz-ibx) + 1
bfield_dims(2) = nx+6
bfield_dims(3) = ny+6
bfield_dims(4) = nz+6

error = 0
call h5dopen_f(file_id,"bfield",dset_id,error)
call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,prim(ibx:ibz,-2:nx+3,-2:ny+3,-2:nz+3),bfield_dims,error)
CALL h5dclose_f(dset_id, error)

IF (error == 0) THEN
  PRINT *, "Finished reading magnetic field from hdf5 file"
ELSE
  PRINT *, "Error reading magnetic field hdf5 file"
  STOP
END IF

! Check if turbulence is present, check only when turbulence is calculated in restart
IF (turb_flag == 1) THEN

  prim_dims(1) = 1
  prim_dims(2) = nx+2
  prim_dims(3) = ny+2
  prim_dims(4) = nz+2
  start = (/iturbq-1,0,0,0/)
  error = 0

  CALL h5dopen_f(file_id,"primitive",dset_id,error)
  call h5dget_space_f(dset_id, dspace_id, error)
  call h5sget_simple_extent_dims_f(dspace_id, ddims, tmp, error)

  IF ((ddims(1)-2) >= (iturbq-imin+1) ) THEN !There is turbulence in the previous run
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, prim_dims, error, stride, block)
    call h5screate_simple_f(4, prim_dims, mem_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, prim(imin:iye2,0:nx+1,0:ny+1,0:nz+1), prim_dims, error, mem_id, dspace_id)
    CALL h5dclose_f(dset_id, error)
    call h5sclose_f(mem_id, error)
    call h5sclose_f(dspace_id, error)

    start_s(1) = iturbq-1
    prim_adim(1) = 1
    stride_s = 1
    block_s = 1

    call h5dopen_f(file_id,"prim_a",dset_id,error)
    call h5dget_space_f(dset_id, dspace_id, error) 
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start_s, prim_adim, error, stride_s, block_s)
    call h5screate_simple_f(1, prim_adim, mem_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, prim_a(iturbq), prim_adim, error, mem_id, dspace_id)
    CALL h5dclose_f(dset_id, error)
    call h5sclose_f(mem_id, error)
    call h5sclose_f(dspace_id, error)

    IF (error == 0) THEN
      PRINT *, "Finished reading turbulence from hdf5 file"
    ELSE
      PRINT *, "Error reading turbulence hdf5 file"
      STOP
    END IF

  ELSE !There is no turbulence in the previous run
    CALL h5dclose_f(dset_id, error)
    CALL GetTurb
    prim_a(iturbq) = turb_q_a

    IF (error == 0) THEN
      PRINT *, "Finished initialising turbulence"
    ELSE
      PRINT *, "Error initialising turbulence"
      STOP
    END IF

  ENDIF

ENDIF

IF (rs_burn == 1) THEN !What to do when burning resumes

  start = (/iscaG1-1,0,0,0/)
  prim_dims(1) = 2
  prim_dims(2) = nx+2
  prim_dims(3) = ny+2
  prim_dims(4) = nz+2

  CALL h5dopen_f(file_id,"primitive",dset_id,error)
  call h5dget_space_f(dset_id, dspace_id, error)
  call h5sget_simple_extent_dims_f(dspace_id, ddims, tmp, error)

  IF (ddims(1) == (iscaG2-imin+1) ) THEN !There is burning in the previous run

    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, prim_dims, error, stride, block)
    call h5screate_simple_f(4, prim_dims, mem_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, prim(iscaG1:iscaG2,0:nx+1,0:ny+1,0:nz+1), prim_dims, error, mem_id, dspace_id)
    CALL h5dclose_f(dset_id, error)
    call h5sclose_f(mem_id, error)
    call h5sclose_f(dspace_id, error)

    prim_adim(1) = 1
    call h5dopen_f(file_id,"found_deton_flag",dset_id,error) !For detonation
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,found_deton_flag,prim_adim,error)
    CALL h5dclose_f(dset_id, error)

    start3 = (/0,0,0/) !h5 indexing
    dims(1) = nx
    dims(2) = 1
    dims(3) = nz
    stride3 = 1
    block3 = 1
    error = 0

    ! The only quantity that is not redefined for level-set on the fly (depends on previous run)
    call h5dopen_f(file_id,"nse_flag",dset_id,error)
    call h5dget_space_f(dset_id, dspace_id, error)
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start3, dims, error, stride3, block3)
    call h5screate_simple_f(3, dims, mem_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nse_flag(1:nx,1:1,1:nz), dims, error, mem_id, dspace_id)
    CALL h5dclose_f(dset_id, error)
    call h5sclose_f(mem_id, error)
    call h5sclose_f(dspace_id, error)

    IF(found_deton_flag == 0 .OR. (found_deton_flag == 1 .and. ABS(global_time - found_deton_time) <= 25000.0D0)) THEN
      IF (update_flag == 1) CALL reinitialization2()
      CALL identify_flamegrid()
      CALL compute_flameratio()
    ENDIF

    ! Do the reinitilization and 
    ! compute the geometry for the 2nd level-set
    IF(found_deton_flag == 1) THEN
      IF (update_flag == 1) CALL reinitialization3()
      CALL identify_detongrid()
      CALL compute_detonratio()
    ENDIF

    ! Conmpute the total local fraction
    burn_ratio(:,:,:) = flame_ratio(:,:,:) + deton_ratio(:,:,:)

    ! Backup the data
    flame_ratio_old = flame_ratio
    deton_ratio_old = deton_ratio

    IF (error == 0) THEN
      PRINT *, "Finished reading level sets from hdf5 file"
    ELSE
      PRINT *, "Error reading level sets hdf5 file"
      STOP
    END IF

  ELSE !There is no burning in the previous run

    CALL h5dclose_f(dset_id, error)
    CALL GetFlame
    WRITE(*,*) 'Finished assigning deflagration level set'
    CALL FLAME_INI()
    WRITE(*,*) 'Finished injecting flame energy'
    CALL FIND_AZBAR()	
    CALL FINDHELMTEMP()
    WRITE(*,*) 'Finished initializing flame'

    ! Assign floor variables !
    CALL CUSTOM_CHECKRHO

  ENDIF

ENDIF

ENDSUBROUTINE