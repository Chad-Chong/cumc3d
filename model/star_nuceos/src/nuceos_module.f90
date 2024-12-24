!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! This module contains all the tools to access the 
! NUCEOS, for nuclear density matter
! More instruction about NUCEOS can be found from GR1D 
! This part is extracted from GR1D for reading the table
! Written by Leung Shing Chi in 2016   
! Modified a bit by Cheung Siu Hei in 2024 for CUMC3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module nuceos_module
implicit none

integer,save :: nrho,ntemp,nye
real*8,save  :: energy_shift = 0.0d0

! min-max values:
real*8,save  :: eos_rhomin,eos_rhomax
real*8,save  :: eos_yemin,eos_yemax
real*8,save  :: eos_tempmin,eos_tempmax

real*8  :: t_max_hack = 240.0d0

! basics
integer, parameter :: nvars = 19
real*8, allocatable,save :: alltables(:,:,:,:)
real*8, allocatable,save :: custom_table(:,:,:,:)

real*8,allocatable,save :: logrho(:)
real*8,allocatable,save :: logtemp(:)
real*8,allocatable,save :: ye(:)

   contains

      !====================================================!
      !== This subroutine reads the .h5 format EOS table ==!
      !== And initializes all the variable in the module ==!
      !====================================================!
      subroutine readtable(eos_filename)
      use hdf5 
      implicit none

      character(*) eos_filename
      character(len=100) message
  
      integer(HID_T) file_id,dset_id,dspace_id
      integer(HSIZE_T) dims1(1), dims3(3)
      integer error,rank,accerr
      integer i,j,k
      INTEGER(HID_T), DIMENSION(nvars) :: dset_ids
      CHARACTER(LEN=50), DIMENSION(nvars) :: dataset_names

      accerr=0
  
      write(*,*) "Reading Nuclear EOS Table"
  
      call h5open_f(error)
      call h5fopen_f (trim(adjustl(eos_filename)), H5F_ACC_RDONLY_F, file_id, error)
      write(6,*) trim(adjustl(eos_filename))
  
!!!   read scalars
      dims1(1)=1
      call h5dopen_f(file_id, "pointsrho", dset_id, error)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nrho, dims1, error)
      call h5dclose_f(dset_id,error)
  
      if(error.ne.0) then
         stop "Could not read EOS table file"
      endif
  
      dims1(1)=1
      call h5dopen_f(file_id, "pointstemp", dset_id, error)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntemp, dims1, error)
      call h5dclose_f(dset_id,error)
  
      if(error.ne.0) then
         stop "Could not read EOS table file"
      endif
  
      dims1(1)=1
      call h5dopen_f(file_id, "pointsye", dset_id, error)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nye, dims1, error)
      call h5dclose_f(dset_id,error)
  
      if(error.ne.0) then
         stop "Could not read EOS table file"
      endif
  
      write(message,"(a25,1P,3i5)") "We have nrho ntemp nye: ", nrho,ntemp,nye
      write(*,*) message
  
      allocate(alltables(nrho,ntemp,nye,nvars))
      ! index variable mapping in .h5 file:
      !  1 -> logpress
      !  2 -> logenergy
      !  3 -> entropy
      !  4 -> munu
      !  5 -> cs2
      !  6 -> dedT
      !  7 -> dpdrhoe
      !  8 -> dpderho
      !  9 -> muhat
      ! 10 -> mu_e
      ! 11 -> mu_p
      ! 12 -> mu_n
      ! 13 -> xa
      ! 14 -> xh
      ! 15 -> xn
      ! 16 -> xp
      ! 17 -> abar
      ! 18 -> zbar
      ! 19 -> gamma
  
      dims3(1)=nrho
      dims3(2)=ntemp
      dims3(3)=nye

      ! Initialize variables
      accerr = 0

      ! List of 3d dataset names
      ! Cheung Siu Hei changed the reading of same shaped 3d datasets to 3 seperated loops
      ! Weirdly the compiler only look at the first element to determine the number of char
      ! If you remove the space at the end of "logpress ", it will only read the first 8 characters
      ! for all the other names: "logenergy" -> "logenerg"
      dataset_names = ["logpress ", "logenergy", "entropy", "munu", "cs2", &
         "dedt", "dpdrhoe", "dpderho", "muhat", "mu_e", &
         "mu_p", "mu_n", "Xa", "Xh", "Xn", "Xp", "Abar", &
         "Zbar", "gamma"]

      ! Open all 3d datasets
      DO i = 1, nvars
         CALL h5dopen_f(file_id, TRIM(dataset_names(i)), dset_ids(i), error)
         CALL h5dread_f(dset_ids(i), H5T_NATIVE_DOUBLE, alltables(:,:,:,i), dims3, error)
         CALL h5dclose_f(dset_ids(i), error)
         accerr = accerr + error
      END DO

      allocate(logrho(nrho))
      dims1(1)=nrho
      call h5dopen_f(file_id, "logrho", dset_id, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logrho, dims1, error)
      call h5dclose_f(dset_id,error)
      accerr=accerr+error
  
      allocate(logtemp(ntemp))
      dims1(1)=ntemp
      call h5dopen_f(file_id, "logtemp", dset_id, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logtemp, dims1, error)
      call h5dclose_f(dset_id,error)
      accerr=accerr+error
  
      allocate(ye(nye))
      dims1(1)=nye
      call h5dopen_f(file_id, "ye", dset_id, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ye, dims1, error)
      call h5dclose_f(dset_id,error)
      accerr=accerr+error
  
  
      call h5dopen_f(file_id, "energy_shift", dset_id, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, energy_shift, dims1, error)
      call h5dclose_f(dset_id,error)
      accerr=accerr+error
  
      if(accerr.ne.0) then
        stop "Problem reading EOS table file"
      endif
  
      call h5fclose_f (file_id,error)
      call h5close_f (error)
    
      eos_rhomin = 10.0d0**logrho(1)
      eos_rhomax = 10.0d0**logrho(nrho)
      write (*,*) 'eos_rhomin', eos_rhomin
      write (*,*) 'eos_rhomax', eos_rhomax

      eos_yemin = ye(1)
      eos_yemax = ye(nye)
      write (*,*) 'eos_yemin', eos_yemin
      write (*,*) 'eos_yemax', eos_yemax

      eos_tempmin = 10.0d0**logtemp(1)
      eos_tempmax = 10.0d0**logtemp(ntemp)
      write (*,*) 'eos_tempmin', eos_tempmin
      write (*,*) 'eos_tempmax', eos_tempmax

      allocate(custom_table(nrho,ntemp,nye,2))
      ! Assign required variables to the temporary array
      custom_table(:,:,:,1)  = alltables(:,:,:,1)    ! logpress
      custom_table(:,:,:,2)  = alltables(:,:,:,5)    ! cs squared

      write(6,*) "Done reading eos tables"

      end subroutine readtable

      subroutine nuc_eos_custom(xrho,xtemp,xye,xprs,xcs2)
      implicit none

      real*8, intent(in)    :: xrho,xye,xtemp
      real*8, intent(out)   :: xprs,xcs2
      real*8,parameter :: idealK1 =  1.2435d15 * (0.5d0**(4.d0/3.d0))
      real*8,parameter :: idealgamma = 1.41d0
      real*8 ffx(2,1)

      ! local variables
      real*8 :: lr,lt,y
      
      if (xrho > eos_rhomax)   stop "nuc_eos_custom: rho > rhomax"
      if (xtemp > eos_tempmax) stop "nuc_eos_custom: temp > tempmax"
      if (xye > eos_yemax)     stop "nuc_eos_custom: ye > yemax"
      if (xye < eos_yemin)     stop "nuc_eos_custom: ye < yemin"

      if(xrho.lt.eos_rhomin*1.2d0 .or. xtemp.lt.eos_tempmin) then
         xprs = idealK1*xrho**(idealgamma)
         xcs2 = idealgamma*xprs/xrho
         return
      endif

      lr = log10(xrho)
      lt = log10(xtemp)
      y = xye
  
      call intp3d_many(lr,lt,y,ffx,1,custom_table,nrho,ntemp,nye,2,logrho,logtemp,ye)

      xprs = 10.0d0**ffx(1,1)
      xcs2 = ffx(2,1)
      end subroutine nuc_eos_custom

      subroutine nuc_eos_one(xrho,xtemp,xye,xvalue,index)
      implicit none

      real*8, intent(in)    :: xye,xrho,xtemp
      real*8, intent(out)   :: xvalue
      integer, intent(in)    :: index

      ! local variables
      real*8 :: lr,lt,y
      real*8 :: ff(1)
      real*8 :: ffx(1,1)

      if(xrho.gt.eos_rhomax)   THEN; WRITE (*,*) "nuc_eos_one: rho > rhomax"  ; CALL TRACEBACKQQ; ENDIF
      if(xrho.lt.eos_rhomin)   THEN; WRITE (*,*) "nuc_eos_one: rho < rhomin"  ; CALL TRACEBACKQQ; ENDIF
      if(xye.gt.eos_yemax)     THEN; WRITE (*,*) "nuc_eos_one: ye > yemax"    ; CALL TRACEBACKQQ; ENDIF
      if(xye.lt.eos_yemin)     THEN; WRITE (*,*) "nuc_eos_one: ye < yemin"    ; CALL TRACEBACKQQ; ENDIF
      if(xtemp.gt.eos_tempmax) THEN; WRITE (*,*) "nuc_eos_one: temp > tempmax"; CALL TRACEBACKQQ; ENDIF
      if(xtemp.lt.eos_tempmin) THEN; WRITE (*,*) "nuc_eos_one: temp < tempmin"; CALL TRACEBACKQQ; ENDIF

      lr = log10(xrho)
      lt = log10(xtemp)
      y = xye

      call intp3d_many(lr,lt,y,ffx,1,alltables(:,:,:,index), &
            nrho,ntemp,nye,1,logrho,logtemp,ye)
      ff(:) = ffx(:,1)

      xvalue = ff(1)
      if (index.eq.1) xvalue = 10.0d0**xvalue                  ! pressure 
      if (index.eq.2) xvalue = 10.0d0**xvalue - energy_shift   ! energy

      end subroutine nuc_eos_one

      ! Written by Cheung Siu Hei for findnuctemp only !
      RECURSIVE subroutine nuc_eos_findtemp(xrho,xtemp,xye,xenr,keyerr,rfeps)
      implicit none
  
      real*8, intent(in)    :: xrho,xye,xenr
      real*8, intent(out)   :: xtemp
      real*8, intent(in)    :: rfeps
      integer, intent(out)  :: keyerr
  
      ! local variables
      real*8 :: lr,lt,y,xx,xeps,leps
      real*8 :: d1,d2,d3,ff(8)
      integer :: keyerrt = 0
      
      if (xrho > eos_rhomax)  stop "nuc_eos_findtemp: rho > rhomax"
      if (xye > eos_yemax)    stop "nuc_eos_findtemp: ye > yemax"
      if (xye < eos_yemin)    stop "nuc_eos_findtemp: ye < yemin"

      lr = log10(xrho)
      lt = log10(xtemp)
      y = xye
      xeps = xenr + energy_shift
      leps = log10(max(xeps,1.0d0))
  
      keyerr = 0
      call findtemp(lr,lt,y,leps,keyerrt,rfeps)
      if(keyerrt.ne.0) then
         keyerr = keyerrt
         return
      endif
      xtemp = 10.0d0**lt
      end subroutine nuc_eos_findtemp


      !SVN:EOSdriver Revision #7
      RECURSIVE subroutine findtemp(lr,lt0,y,epsin,keyerrt,rfeps)
      implicit none

      real*8 lr,lt0,y,epsin
      real*8 eps,lt,ldt
      real*8 tol
      real*8 d1,d2,d3
      real*8 eps0,eps1,lt1
    
      real*8 ltn,ltmax,ltmin
      real*8 tinput,rfeps
    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !my check
      real*8 eps_min, eps_max, eosdummy
      integer keyerr, keytemp 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
      integer :: rl = 0
      integer itmax,i,keyerrt
      integer ii,jj,kk
      
      keyerrt=0
    
      tol=rfeps  ! need to find energy to less than 1 in 10^-6
      itmax=100 ! use at most 25 iterations, then bomb
    
      lt=lt0
      lt1=lt 
    
      eps0=epsin
      eps1=eps0
      
      ltmax=logtemp(ntemp)
      ltmin=logtemp(1)
    
      ! Note: We are using Ewald's Lagrangian interpolator here!
    
      !preconditioning 1: do we already have the right temperature?
      call intp3d(lr,lt,y,eps,1,alltables(:,:,:,2),nrho,ntemp,nye,logrho,logtemp,ye,d1,d2,d3)

      if (abs(eps-eps0).lt.tol*abs(eps0)) then
         return
      endif
      lt1=lt
      eps1=eps

      do i=1,itmax
         !d2 is the derivative deps/dlogtemp;
         ldt = -(eps - eps0)/d2 
         ltn = lt+ldt
         ltn = min(ltn,ltmax)
         ltn = max(ltn,ltmin)
         lt1=lt
         lt=ltn
         eps1=eps
         call intp3d(lr,lt,y,eps,1,alltables(:,:,:,2),nrho,ntemp,nye,logrho,logtemp,ye,d1,d2,d3)
         if (abs(eps - eps0).lt.tol*abs(eps0)) then
            exit
         endif
         if(lt <= log10(eos_tempmin)) exit

         if(abs(eps-eps0).lt.1.0d-3*abs(eps0)) then
            d2 = (eps-eps1)/(lt-lt1)
         endif
     enddo
    
     if(i.ge.itmax) then
        keyerrt=667
        call bisection(lr,lt0,y,eps0,lt,alltables(:,:,:,2),keyerrt,1)
        if(keyerrt.eq.667) then
           if(lt.ge.log10(t_max_hack)) then
              ! handling for too high temperatures
              lt = log10(t_max_hack)
              keyerrt=0
              goto 12
           else if(abs(lt-log10(t_max_hack))/log10(t_max_hack).lt.0.025d0) then
              lt0 = min(lt,log10(t_max_hack))
              keyerrt=0
              goto 12
           else
           endif
        endif
        
        lt0=min(lt,log10(t_max_hack))
        return
      endif
    
      12 continue
    
      lt0=min(lt,log10(t_max_hack))
    
      end subroutine findtemp

      ! SVN:EOSdriver Revision #7
      RECURSIVE SUBROUTINE intp3d ( x, y, z, f, kt, ft, nx, ny, nz, xt, yt, zt, &
                        d1, d2, d3 )     
      implicit none

      integer kt,nx,ny,nz,ktx
      real*8 x,y,z,f
      real*8 xt(nx),yt(ny),zt(nz)
      real*8 ft(nx,ny,nz)
      real*8 d1,d2,d3

      PARAMETER   (ktx = 400)
      real*8  fh(ktx,8), delx(ktx), dely(ktx), delz(ktx), &
                 a1(ktx), a2(ktx), a3(ktx), a4(ktx),                &
                 a5(ktx), a6(ktx), a7(ktx), a8(ktx)

      real*8 dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi
      integer n,ix,iy,iz

      IF (kt .GT. ktx)  STOP '***KTX**'

      dx    = (xt(nx) - xt(1)) / FLOAT(nx-1)
      dy    = (yt(ny) - yt(1)) / FLOAT(ny-1)
      dz    = (zt(nz) - zt(1)) / FLOAT(nz-1)

      dxi   = 1. / dx
      dyi   = 1. / dy
      dzi   = 1. / dz

      dxyi  = dxi * dyi
      dxzi  = dxi * dzi
      dyzi  = dyi * dzi

      dxyzi = dxi * dyi * dzi

      DO  n = 1, kt
         ix = 2 + INT( (x - xt(1) - 1.e-10) * dxi )
         iy = 2 + INT( (y - yt(1) - 1.e-10) * dyi )
         iz = 2 + INT( (z - zt(1) - 1.e-10) * dzi )
                                                     
         ix = MAX( 2, MIN( ix, nx ) )
         iy = MAX( 2, MIN( iy, ny ) )
         iz = MAX( 2, MIN( iz, nz ) )
               
         delx(n) = xt(ix) - x
         dely(n) = yt(iy) - y
         delz(n) = zt(iz) - z
                                                                   
         fh(n,1) = ft(ix  , iy  , iz  )                             
         fh(n,2) = ft(ix-1, iy  , iz  )                             
         fh(n,3) = ft(ix  , iy-1, iz  )                             
         fh(n,4) = ft(ix  , iy  , iz-1)                             
         fh(n,5) = ft(ix-1, iy-1, iz  )                             
         fh(n,6) = ft(ix-1, iy  , iz-1)                             
         fh(n,7) = ft(ix  , iy-1, iz-1)                             
         fh(n,8) = ft(ix-1, iy-1, iz-1)                             

         a1(n) = fh(n,1)                             
         a2(n) = dxi   * ( fh(n,2) - fh(n,1) )       
         a3(n) = dyi   * ( fh(n,3) - fh(n,1) )       
         a4(n) = dzi   * ( fh(n,4) - fh(n,1) )       
         a5(n) = dxyi  * ( fh(n,5) - fh(n,2) - fh(n,3) + fh(n,1) )
         a6(n) = dxzi  * ( fh(n,6) - fh(n,2) - fh(n,4) + fh(n,1) )
         a7(n) = dyzi  * ( fh(n,7) - fh(n,3) - fh(n,4) + fh(n,1) )
         a8(n) = dxyzi * ( fh(n,8) - fh(n,1) + fh(n,2) + fh(n,3) + &
                          fh(n,4) - fh(n,5) - fh(n,6) - fh(n,7) )

         d1 = -a2(n)
         d2 = -a3(n)
         d3 = -a4(n)
         f  = a1(n) +  a2(n) * delx(n)                          & 
                    +  a3(n) * dely(n)                          &
                    +  a4(n) * delz(n)                          &
                    +  a5(n) * delx(n) * dely(n)                &
                    +  a6(n) * delx(n) * delz(n)                &
                    +  a7(n) * dely(n) * delz(n)                &
                    +  a8(n) * delx(n) * dely(n) * delz(n)     

      ENDDO
                                                                    
      RETURN                                                         
      END SUBROUTINE

      RECURSIVE SUBROUTINE intp3d_many ( x, y, z, f, kt, ft, nx, ny, nz, nvars, xt, yt, zt)
      implicit none

      integer kt,nx,ny,nz,iv,nvars
      real*8 :: ft(nx,ny,nz,nvars)

      real*8 x,y,z,f(nvars)
      real*8 xt(nx),yt(ny),zt(nz)
      real*8 d1,d2,d3

      integer,parameter :: ktx = 1
      real*8  fh(ktx,8,nvars), delx(ktx), dely(ktx), delz(ktx), &
           a1(ktx,nvars), a2(ktx,nvars), a3(ktx,nvars), a4(ktx,nvars), &
           a5(ktx,nvars), a6(ktx,nvars), a7(ktx,nvars), a8(ktx,nvars)

      real*8 dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi
      integer n,ix,iy,iz

      IF (kt .GT. ktx)  STOP '***KTX**'

      dx    = (xt(nx) - xt(1)) / FLOAT(nx-1)
      dy    = (yt(ny) - yt(1)) / FLOAT(ny-1)
      dz    = (zt(nz) - zt(1)) / FLOAT(nz-1)

      dxi   = 1. / dx
      dyi   = 1. / dy
      dzi   = 1. / dz

      dxyi  = dxi * dyi
      dxzi  = dxi * dzi
      dyzi  = dyi * dzi

      dxyzi = dxi * dyi * dzi

      dO  n = 1, kt                                            
         ix = 2 + INT( (x - xt(1) - 1.e-10) * dxi )
         iy = 2 + INT( (y - yt(1) - 1.e-10) * dyi )
         iz = 2 + INT( (z - zt(1) - 1.e-10) * dzi )
                                                     
         ix = MAX( 2, MIN( ix, nx ) )
         iy = MAX( 2, MIN( iy, ny ) )
         iz = MAX( 2, MIN( iz, nz ) )
                                                               
         delx(n) = xt(ix) - x
         dely(n) = yt(iy) - y
         delz(n) = zt(iz) - z

         do iv = 1, nvars
            fh(n,1,iv) = ft(ix  , iy  , iz, iv  )                             
            fh(n,2,iv) = ft(ix-1, iy  , iz, iv  )                             
            fh(n,3,iv) = ft(ix  , iy-1, iz, iv  )                             
            fh(n,4,iv) = ft(ix  , iy  , iz-1, iv)                             
            fh(n,5,iv) = ft(ix-1, iy-1, iz, iv  )                             
            fh(n,6,iv) = ft(ix-1, iy  , iz-1, iv)                             
            fh(n,7,iv) = ft(ix  , iy-1, iz-1, iv)                             
            fh(n,8,iv) = ft(ix-1, iy-1, iz-1, iv)                             

            a1(n,iv) = fh(n,1,iv)                             
            a2(n,iv) = dxi   * ( fh(n,2,iv) - fh(n,1,iv) )       
            a3(n,iv) = dyi   * ( fh(n,3,iv) - fh(n,1,iv) )       
            a4(n,iv) = dzi   * ( fh(n,4,iv) - fh(n,1,iv) )       
            a5(n,iv) = dxyi  * ( fh(n,5,iv) - fh(n,2,iv) - fh(n,3,iv) + fh(n,1,iv) )
            a6(n,iv) = dxzi  * ( fh(n,6,iv) - fh(n,2,iv) - fh(n,4,iv) + fh(n,1,iv) )
            a7(n,iv) = dyzi  * ( fh(n,7,iv) - fh(n,3,iv) - fh(n,4,iv) + fh(n,1,iv) )
            a8(n,iv) = dxyzi * ( fh(n,8,iv) - fh(n,1,iv) + fh(n,2,iv) + fh(n,3,iv) + &
                 fh(n,4,iv) - fh(n,5,iv) - fh(n,6,iv) - fh(n,7,iv) )

            f(iv)  = a1(n,iv) +  a2(n,iv) * delx(n)            &
                 +  a3(n,iv) * dely(n)                         &
                 +  a4(n,iv) * delz(n)                         &
                 +  a5(n,iv) * delx(n) * dely(n)               &
                 +  a6(n,iv) * delx(n) * delz(n)               &
                 +  a7(n,iv) * dely(n) * delz(n)               &
                 +  a8(n,iv) * delx(n) * dely(n) * delz(n)     
         enddo
      enddo
      
      end SUBROUTINE intp3d_many

      RECURSIVE subroutine bisection(lr,lt0,y,eps0,lt,bivar,keyerrt,keybisect)
      implicit none

      real*8 lr,lt0,y,eps0,lt
      integer keyerrt
  
      integer keybisect

      !temporary vars
      real*8 lt1,lt2,ltmin,ltmax
      real*8 f1,f2,fmid,dlt,ltmid
      real*8 d1,d2,d3,tol
      real*8 f1a,f2a
  
      integer bcount,i,itmax,maxbcount
  
      real*8 bivar(*)
  
      bcount = 0
      maxbcount = 80
  
      tol=1.d-9 ! need to find energy to less than 1 in 10^-9
      itmax=50
  
      ltmax=logtemp(ntemp)
      ltmin=logtemp(1)
  
      lt = lt0
      lt1 = dlog10(min(10.0d0**ltmax,1.10d0*(10.0d0**lt0)))
      lt2 = dlog10(max(10.0d0**ltmin,0.90d0*(10.0d0**lt0)))

      call intp3d(lr,lt1,y,f1a,1,bivar,nrho,ntemp,nye,logrho,logtemp,ye,d1,d2,d3)
      call intp3d(lr,lt2,y,f2a,1,bivar,nrho,ntemp,nye,logrho,logtemp,ye,d1,d2,d3)

      f1=f1a-eps0
      f2=f2a-eps0
  
      keyerrt=0
      
      do while(f1*f2.ge.0.0d0)
         bcount=bcount+1
         lt1=dlog10(min(10.0d0**ltmax,1.2d0*(10.0d0**lt1)))
         lt2=dlog10(max(10.0d0**ltmin,0.8d0*(10.0d0**lt2)))
         call intp3d(lr,lt1,y,f1a,1,bivar,nrho,ntemp,nye,logrho,logtemp,ye,d1,d2,d3)
         call intp3d(lr,lt2,y,f2a,1,bivar,nrho,ntemp,nye,logrho,logtemp,ye,d1,d2,d3)
         f1=f1a-eps0
         f2=f2a-eps0
         if(bcount.ge.maxbcount) then
            keyerrt=667
            return
         endif
      enddo
   
      if(f1.lt.0.0d0) then
         lt=lt1
         dlt=dlog10((10.0D0**lt2)-(10.0d0**lt1))
      else
         lt=lt2
         dlt=dlog10((10.0D0**lt1)-(10.0d0**lt2))
      endif
  
      do i=1,itmax
         dlt=dlog10((10.0d0**dlt)*0.5d0)
         ltmid=dlog10(10.0d0**lt+10.0d0**dlt)
         call intp3d(lr,ltmid,y,f2a,1,bivar,nrho,ntemp,nye,logrho,logtemp,ye,d1,d2,d3)
         fmid=f2a-eps0
         if(fmid.le.0.0d0) lt=ltmid
         if(abs(1.0d0-f2a/eps0).lt.tol) then
            lt=ltmid
            return
         endif
      enddo

      end subroutine bisection
end module nuceos_module

