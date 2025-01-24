!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output the hydrodynamic variable profiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE print_hydroprofile
USE HDF5
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

! for HDF5 !
integer :: error, space_rank
character(len=99) :: globalt
character(len=99) :: filename
integer(HID_T) :: file_id, dspace_id, dset_id1
integer(HSIZE_T) :: dims3(3), dims4(4), dist_dims(1)

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
call h5dcreate_f(file_id,"time",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,global_time,dist_dims,error)

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
call h5dcreate_f(file_id,"dimension",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,DBLE(2),dist_dims,error)

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
space_rank = 4
dims4(1) = (ibx - 1 - imin) + 1
dims4(2) = nx
dims4(3) = ny
dims4(4) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dims4,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"primitiveNM",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,prim(imin:ibx-1,1:nx,1:ny,iNM),dims4,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
dims3(1) = nx
dims3(2) = ny
dims3(3) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dims3,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"epsilonNM",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,epsilon(1:nx,1:ny,iNM),dims3,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 4
dims4(1) = (ibx - 1 - imin) + 1
dims4(2) = nx
dims4(3) = ny
dims4(4) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dims4,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"primitiveDM",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,prim(imin:ibx-1,1:nx,1:ny,iDM),dims4,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
dims3(1) = nx
dims3(2) = ny
dims3(3) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dims3,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"epsilonDM",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,epsilon(1:nx,1:ny,iDM),dims3,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Output potential !
dims3(1) = nx + 1

! open dataspace !
call h5screate_simple_f(space_rank,dims3,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"phi",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! write dataset !
call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,phi(0:nx,1:ny,1:nz),dims3,error)

! close dataset !
call h5dclose_f(dset_id1,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
! space_rank = 4
! dims4(1) = (ibz - ibx) + 1
! dims4(2) = nx + 1
! dims4(3) = ny + 1
! dims4(4) = 1

! ! open dataspace !
! call h5screate_simple_f(space_rank,dims4,dspace_id,error)

! ! create dataset !
! call h5dcreate_f(file_id,"bfield",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)

! ! write dataset !
! call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,prim(ibx:ibz,0:nx,0:ny,iNM),dims4,error)

! ! close dataset !
! call h5dclose_f(dset_id1,error)

! ! close data space !
! call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! close the file !
call h5fclose_f(file_id,error)

! close interface !
call h5close_f(error)

END SUBROUTINE
