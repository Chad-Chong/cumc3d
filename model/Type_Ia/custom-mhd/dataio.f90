!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output the hydrodynamic variable profiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE print_hydroprofile
USE HDF5
USE DEFINITION
USE CUSTOM_DEF
USE MHD_MODULE
IMPLICIT NONE

! for HDF5 !
integer :: error, space_rank
character(len=99) :: globalt
character(len=99) :: filename
integer(HID_T) :: file_id, dspace_id, dset_id1
integer(HSIZE_T) :: sup_dims(3), data_dims(4), dist_dims(1),df_dims(4)

! integer !
INTEGER :: j  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for GPU !

#ifdef GPU
!$ACC UPDATE HOST(prim(imin:imax,:,:,:), epsilon)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! write to character !
write(globalt,'(I)') n_iter

! assign !
filename = './outfile/rkiter-'// trim(adjustl(globalt)) //'-nm.hdf5'

! create interface !
call h5open_f(error)

! open the file !
call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"temp2_a",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,temp2_a,dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"eps_a",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,eps_a,dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"found_deton_flag",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,found_deton_flag,dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"time",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,global_time,dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (burn_flag == 1) THEN
    ! define DIMENSION !
    space_rank = 1
    dist_dims(1) = 1

    ! open dataspace !
    call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

    ! create dataset !
    call h5dcreate_f(file_id,"Enuc",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

    ! write dataset !
    call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,Enuc,dist_dims,error)

    ! close dataset !
    call h5dclose_f(dset_id1,error)

    ! close data space !
    call h5sclose_f(dspace_id,error)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
! define DIMENSION !
space_rank = 1
dist_dims(1) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"dimension",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,DBLE(n_dim),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"coordinate",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,DBLE(coordinate_flag+1),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = ibx-imin

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"prim_a",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,prim_a,dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = nx + 3

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"x-interface",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,xF(-1:nx+1),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = ny + 3

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"y-interface",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,yF(-1:ny+1),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
dist_dims(1) = nz + 3

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"z-interface",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,zF(-1:nz+1),dist_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
IF (levelset_flag == 1 .and. det_count /= 0 .and. deton_flag == 1) THEN
    space_rank = 1
    dist_dims(1) = det_count

    ! open dataspace !
    call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

    ! create dataset !
    call h5dcreate_f(file_id,"detonations",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

    ! write dataset !
    call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,det_times,dist_dims,error)

    ! close dataset !
    call h5dclose_f(dset_id1,error)

    ! close data space !
    call h5sclose_f(dspace_id,error)

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 4
data_dims(1) = (ibx - 1 - imin) + 1
data_dims(2) = nx + 2
data_dims(3) = ny + 2
data_dims(4) = nz + 2

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"primitive",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,prim(imin:ibx-1,0:nx+1,0:ny+1,0:nz+1),data_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 3
! sup_dims(1) = nx
! sup_dims(2) = ny 
! sup_dims(3) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"check",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE, debug_flux(1,ivx,1:nx,1:ny,1:nz)+debug_flux(2,ivx,1:nx,1:ny,1:nz)+debug_flux(3,ivx,1:nx,1:ny,1:nz)+sc(ivx,1:nx,1:ny,1:nz),sup_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 4
df_dims(1) = 3
df_dims(2) = nx
df_dims(3) = ny
df_dims(4) = nz

! open dataspace !
call h5screate_simple_f(space_rank,df_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"l_rk",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,l_rk(ivx:ivz,1:nx,1:ny,1:nz),df_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 4
! df_dims(1) = 3
! df_dims(2) = nx
! df_dims(3) = ny
! df_dims(4) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,df_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"debug_flux_x",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,debug_flux(1,ivx:ivz,1:nx,1:ny,1:nz),df_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 4
! df_dims(1) = 3
! df_dims(2) = nx
! df_dims(3) = ny
! df_dims(4) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,df_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"gravity_sc",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,gravity_sc(ivx:ivz,1:nx,1:ny,1:nz),df_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 3
! sup_dims(1) = nx
! sup_dims(2) = ny 
! sup_dims(3) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"sc(ivx)",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE, sc(ivx,1:nx,1:ny,1:nz),sup_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 3
! sup_dims(1) = nx
! sup_dims(2) = ny 
! sup_dims(3) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"sc(ivy)",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE, sc(ivy,1:nx,1:ny,1:nz),sup_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 3
! sup_dims(1) = nx
! sup_dims(2) = ny 
! sup_dims(3) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"sc(ivz)",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE, sc(ivz,1:nx,1:ny,1:nz),sup_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 4
! df_dims(1) = 3
! df_dims(2) = nx
! df_dims(3) = ny
! df_dims(4) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,df_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"debug_flux_y",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,debug_flux(2,ivx:ivz,1:nx,1:ny,1:nz),df_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 4
! df_dims(1) = 3
! df_dims(2) = nx
! df_dims(3) = ny
! df_dims(4) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,df_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"debug_flux_z",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,debug_flux(3,ivx:ivz,1:nx,1:ny,1:nz),df_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 3
! sup_dims(1) = nx 
! sup_dims(2) = ny 
! sup_dims(3) = nz

! CALL FIND_DIVB

! ! open dataspace !
! call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"divb",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,divb_arr(1:nx,1:ny,1:nz),sup_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
sup_dims(1) = nx 
sup_dims(2) = ny 
sup_dims(3) = nz


! open dataspace !
call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"efield_x",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,efield_x(1:nx,1:ny,1:nz),sup_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
sup_dims(1) = nx 
sup_dims(2) = ny 
sup_dims(3) = nz


! open dataspace !
call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"efield_y",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,efield_y(1:nx,1:ny,1:nz),sup_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
sup_dims(1) = nx 
sup_dims(2) = ny 
sup_dims(3) = nz


! open dataspace !
call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"efield_z",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,efield_z(1:nx,1:ny,1:nz),sup_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
sup_dims(1) = nx+2 
sup_dims(2) = ny+2
sup_dims(3) = nz+2

! open dataspace !
call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"eface(iyx)",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,eface(iyx,0:nx+1,0:ny+1,0:nz+1),sup_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
sup_dims(1) = nx+2 
sup_dims(2) = ny+2
sup_dims(3) = nz+2

! open dataspace !
call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"eface(iyz)",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,eface(iyz,0:nx+1,0:ny+1,0:nz+1),sup_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 3
! sup_dims(1) = nx+1 
! sup_dims(2) = ny+1 
! sup_dims(3) = nz+1

! ! open dataspace !
! call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"pz_efield_x-2*px_efield_z",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,pz_efield_x(0:nx,0:ny,0:nz)-2*px_efield_z(0:nx,0:ny,0:nz),sup_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 3
! sup_dims(1) = nx
! sup_dims(2) = ny 
! sup_dims(3) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"vx Bphi",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE, prim(ivx,1:nx,1:ny,1:nz)*bcell(iby,1:nx,1:ny,1:nz),sup_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 3
! sup_dims(1) = nx
! sup_dims(2) = ny 
! sup_dims(3) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"vphi Bx",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE, prim(ivy,1:nx,1:ny,1:nz)*bcell(ibx,1:nx,1:ny,1:nz),sup_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 3
! sup_dims(1) = nx
! sup_dims(2) = ny 
! sup_dims(3) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"vz Bphi",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE, prim(ivz,1:nx,1:ny,1:nz)*bcell(iby,1:nx,1:ny,1:nz),sup_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! space_rank = 3
! sup_dims(1) = nx
! sup_dims(2) = ny 
! sup_dims(3) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"vphi Bz",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE, prim(ivy,1:nx,1:ny,1:nz)*bcell(ibz,1:nx,1:ny,1:nz),sup_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! ! define DIMENSION !
! space_rank = 3
! sup_dims(1) = nx
! sup_dims(2) = ny 
! sup_dims(3) = nz

! ! open dataspace !
! call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"sc(ivx)",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE, sc(ivx,1:nx,1:ny,1:nz),sup_dims,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
sup_dims(1) = nx + 2
sup_dims(2) = ny + 2
sup_dims(3) = nz + 2

! open dataspace !
call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"epsilon",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,epsilon(0:nx+1,0:ny+1,0:nz+1),sup_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
sup_dims(1) = nx + 2
sup_dims(2) = ny + 2
sup_dims(3) = nz + 2

! open dataspace !
call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"cs",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,cs(0:nx+1,0:ny+1,0:nz+1),sup_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
sup_dims(1) = nx
sup_dims(2) = ny
sup_dims(3) = nz

! open dataspace !
call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"lambdas",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,lambdas(1:nx,1:ny,1:nz),sup_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
IF (helmeos_flag == 1) THEN
    space_rank = 3
    sup_dims(1) = nx + 2
    sup_dims(2) = ny + 2
    sup_dims(3) = nz + 2

    ! open dataspace !
    call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

    ! create dataset !
    call h5dcreate_f(file_id,"temp",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

    ! write dataset !
    call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,temp2(0:nx+1,0:ny+1,0:nz+1),sup_dims,error)

    ! close dataset !
    call h5dclose_f(dset_id1,error)

    ! close data space !
    call h5sclose_f(dspace_id,error)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! IF (helmeos_flag == 1) THEN
!     space_rank = 3
!     sup_dims(1) = nx + 2
!     sup_dims(2) = ny + 2
!     sup_dims(3) = nz + 2

!     ! open dataspace !
!     call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

!     ! create dataset !
!     call h5dcreate_f(file_id,"abar",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

!     ! write dataset !
!     call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,abar2(0:nx+1,0:ny+1,0:nz+1),sup_dims,error)

!     ! close dataset !
!     call h5dclose_f(dset_id1,error)

!     ! close data space !
!     call h5sclose_f(dspace_id,error)
! ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! define DIMENSION !
! IF (helmeos_flag == 1) THEN
!     space_rank = 3
!     sup_dims(1) = nx + 2
!     sup_dims(2) = ny + 2
!     sup_dims(3) = nz + 2

!     ! open dataspace !
!     call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

!     ! create dataset !
!     call h5dcreate_f(file_id,"zbar",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

!     ! write dataset !
!     call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,zbar2(0:nx+1,0:ny+1,0:nz+1),sup_dims,error)

!     ! close dataset !
!     call h5dclose_f(dset_id1,error)

!     ! close data space !
!     call h5sclose_f(dspace_id,error)
! ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (n_iter == 0) THEN
    space_rank = 3
    sup_dims(1) = nx
    sup_dims(2) = ny
    sup_dims(3) = nz

    ! open dataspace !
    call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

    ! create dataset !
    call h5dcreate_f(file_id,"vol",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

    ! write dataset !
    call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,vol(1:nx,1:ny,1:nz),sup_dims,error)

    ! close dataset !
    call h5dclose_f(dset_id1,error)

    ! close data space !
    call h5sclose_f(dspace_id,error)
ENDIF

! define DIMENSION !
IF (levelset_flag == 1) THEN
    space_rank = 3
    sup_dims(1) = nx + 2
    sup_dims(2) = ny + 2
    sup_dims(3) = nz + 2

    ! open dataspace !
    call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

    ! create dataset !
    call h5dcreate_f(file_id,"flame_state",H5T_NATIVE_INTEGER,dspace_id,dset_id1,error)

    ! write dataset !
    call h5dwrite_f(dset_id1,H5T_NATIVE_INTEGER,flamegrid_flag(0:nx+1,0:ny+1,0:nz+1),sup_dims,error)

    ! close dataset !
    call h5dclose_f(dset_id1,error)

    ! close data space !
    call h5sclose_f(dspace_id,error)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
if (found_deton_flag == 1) THEN
    space_rank = 3
    sup_dims(1) = nx + 2
    sup_dims(2) = ny + 2
    sup_dims(3) = nz + 2

    ! open dataspace !
    call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

    ! create dataset !
    call h5dcreate_f(file_id,"deton_state",H5T_NATIVE_INTEGER,dspace_id,dset_id1,error)

    ! write dataset !
    call h5dwrite_f(dset_id1,H5T_NATIVE_INTEGER,detongrid_flag(0:nx+1,0:ny+1,0:nz+1),sup_dims,error)

    ! close dataset !
    call h5dclose_f(dset_id1,error)

    ! close data space !
    call h5sclose_f(dspace_id,error)

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
IF (gravity_flag == 1) THEN
    space_rank = 3
    sup_dims(1) = nx + 2
    sup_dims(2) = ny + 2
    sup_dims(3) = nz + 2

    ! open dataspace !
    call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

    ! create dataset !
    call h5dcreate_f(file_id,"phi",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

    ! write dataset !
    call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,phi(0:nx+1,0:ny+1,0:nz+1),sup_dims,error)

    ! close dataset !
    call h5dclose_f(dset_id1,error)

    ! close data space !
    call h5sclose_f(dspace_id,error)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
IF (burn_flag == 1) THEN

    space_rank = 3
    sup_dims(1) = nx
    sup_dims(2) = 1
    sup_dims(3) = nz

    ! open dataspace !
    call h5screate_simple_f(space_rank,sup_dims,dspace_id,error)

    ! create dataset !
    call h5dcreate_f(file_id,"nse_flag",H5T_NATIVE_INTEGER,dspace_id,dset_id1,error)

    ! write dataset !
    call h5dwrite_f(dset_id1,H5T_NATIVE_INTEGER,nse_flag(1:nx,1,1:nz),sup_dims,error)

    ! close dataset !
    call h5dclose_f(dset_id1,error)

    ! close data space !
    call h5sclose_f(dspace_id,error)

ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
IF (nuspec_flag == 1) THEN

    space_rank = 1
    dist_dims(1) = 10

    ! open dataspace !
    call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

    ! create dataset !
    call h5dcreate_f(file_id,"nu_phi",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

    ! write dataset !
    call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,nu_phi(:),dist_dims,error)

    ! close dataset !
    call h5dclose_f(dset_id1,error)

    ! close data space !
    call h5sclose_f(dspace_id,error)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 4
data_dims(1) = (ibz - ibx) + 1
data_dims(2) = nx + 1
data_dims(3) = ny + 1
data_dims(4) = nz + 1

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"bcell",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,bcell(ibx:ibz,0:nx,0:ny,0:nz),data_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 4
data_dims(1) = (ibz - ibx) + 1
data_dims(2) = nx + 6
data_dims(3) = ny + 6
data_dims(4) = nz + 6

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"dbdt",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,l_rk(ibx:ibz,-2:nx+3,-2:ny+3,-2:nz+3),data_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 4
data_dims(1) = (ibz - ibx) + 1
data_dims(2) = nx + 6
data_dims(3) = ny + 6
data_dims(4) = nz + 6

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"bfield",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,prim(ibx:ibz,-2:nx+3,-2:ny+3,-2:nz+3),data_dims,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! close the file !
call h5fclose_f(file_id,error)

! close interface !
call h5close_f(error)

END SUBROUTINE
