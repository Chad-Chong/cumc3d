!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! New variables for turbulent flame used when turb_flag == 1
!
! This module contains the components for calculating
! sub-grid scale (SGS) turbulence using the 
! one-equation form (See Niemeyer1995b). 
! This assumes the turbulence can be characterized by the 
! local velocity flucation v', and the associated 
! turbulent kinetic energy q.
! By solving the evolution of q, one can approximate
! the local eddy motion, which can be used for 
! 
! Written by Leung Shing Chi in 2015
! Documented by Leung Shing Chi in 2016
! Updated by Leung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE buildTurb
USE definition
USE CUSTOM_DEF
IMPLICIT NONE

ALLOCATE(turb_source(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(turb_diff(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(turb_eps(-2:nx+3,-2:ny+3,-2:nz+3))

END SUBROUTINE buildTurb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE destroyTurb
USE definition
USE CUSTOM_DEF
IMPLICIT NONE

DEALLOCATE(turb_source)
DEALLOCATE(turb_eps)

END SUBROUTINE destroyTurb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetTurb
USE definition
USE CUSTOM_DEF
IMPLICIT NONE

integer :: m,n

!This is to initialize the sub-grid turbulence-dissipation field

!First give the turbulent kinetic energy field a very small number
prim(iturbq,:,:,:) = turb_q_a

END SUBROUTINE GetTurb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds all component for the sub-grid turbulence

subroutine findturbulence
use DEFINITION
use CUSTOM_DEF
use ieee_arithmetic
implicit none

integer :: i, j, k, m, n, neighbour_flag
real*8 :: sum
real*8 :: turb_C, turb_D, turb_F
real*8 :: at_no, g_eff, phi_r, phi_z
real*8 :: inv_sqrt2 = 1.0D0 / DSQRT(2.0D0)
real*8 :: c_lambda = 2.0D0/3.0D0

! Character for debug
character(len=99) :: globalt

! For Archimedis Production
real*8, allocatable, dimension(:,:,:) :: turb_RT

! For turbulence compression
real*8, allocatable, dimension(:,:,:) :: turb_comp

! For turbulence production
REAL*8 :: dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz
real*8, allocatable, dimension(:,:,:) :: turb_div_v
real*8, allocatable, dimension(:,:,:,:,:) :: turb_str_tensor
real*8, allocatable, dimension(:,:,:) :: turb_str
real*8, allocatable, dimension(:,:,:,:) :: turb_eta
real*8, allocatable, dimension(:,:,:) :: gm
real*8, allocatable, dimension(:,:,:,:,:) :: turb_velgrad

! For turbulence diffusion
REAL*8 :: dqdx, dqdy, dqdz, dqdevdx, dqdevdy, dqdevdz
real*8, allocatable, dimension(:,:,:,:) :: turb_qdev
! if(debug_flag == 1) write(*,*) 'In Find turbulence'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(turb_RT(-2:nx+3,-2:ny+3,-2:nz+3))
allocate(turb_comp(-2:nx+3,-2:ny+3,-2:nz+3))

allocate(turb_div_v(-2:nx+3,-2:ny+3,-2:nz+3))
allocate(turb_str_tensor(-2:nx+3,-2:ny+3,-2:nz+3,3,3))
allocate(turb_str(-2:nx+3,-2:ny+3,-2:nz+3))
allocate(turb_eta(-2:nx+3,-2:ny+3,-2:nz+3,3))
allocate(turb_velgrad(-2:nx+3,-2:ny+3,-2:nz+3,3,3))
allocate(turb_qdev(-2:nx+3,-2:ny+3,-2:nz+3,3))
allocate(gm(-2:nx+3,-2:ny+3,-2:nz+3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (coordinate_flag == 1) THEN

   do k = 0, nz+1, 1
      do j = 0, ny+1, 1
         do i=0, nx+1, 1
            gm(i,j,k) = vol(i,j,k)**(1.0D0/3.0D0)
         enddo
      enddo
   enddo

   IF (turbcheck_flag == 1) THEN
      write(globalt,'(I)') count
      OPEN (UNIT = 127, FILE = './turb_F'// trim(adjustl(globalt)) //'.dat', STATUS = 'REPLACE')
   ENDIF

   do k = 0, nz+1, 1
      do j = 0, ny+1, 1
         do i=0, nx+1, 1

            turb_F = MIN(100.0D0, MAX(0.1D0, 1.0D-4 * epsilon(i,j,k)/prim(iturbq,i,j,k)))
            turb_C = 0.1D0 * turb_F
            turb_D = 0.5D0 / turb_F

            IF (turbcheck_flag == 1) THEN
               IF (i /= 0 .and. j == 1 .and. k /=0 .and. i /= nx+1 .and. k /= nz+1) THEN
                  WRITE(127,*) 1.0D-4 * epsilon(i,j,k)/prim(iturbq,i,j,k)
               ENDIF
            ENDIF

            ! Different viscosity at different direction for diffusion (note that dx(i) = dz(k)) !
            turb_eta(i,j,k,1) = turb_C * 1.6D0 * dx(i)**2/gm(i,j,k) * DSQRT(DSQRT(prim(iturbq,i,j,k)**2))
            turb_eta(i,j,k,2) = turb_C * 1.6D0 * x(i)*dy(j)**2/gm(i,j,k) * DSQRT(DSQRT(prim(iturbq,i,j,k)**2))
            turb_eta(i,j,k,3) = turb_C * 1.6D0 * dz(k)**2/gm(i,j,k) * DSQRT(DSQRT(prim(iturbq,i,j,k)**2))
            ! Same viscosity for dissipation !
            turb_eps(i,j,k) = prim(irho,i,j,k) * turb_D * DSQRT(DSQRT(prim(iturbq,i,j,k)**6)) / (1.6D0 * gm(i,j,k))
            turb_RT(i,j,k) = 0.0D0
            
            IF (gravity_flag == 1 .and. update_flag == 1) THEN
               if(flamegrid_flag(i,j,k) /= 0 .and. flamegrid_flag(i,j,k) /= 1) then
                  at_no = MAX(0.5D0 * (0.0522D0 + 0.145D0 / DSQRT(prim(irho,i,j,k)/1.62D-9) - 0.01D0 / (prim(irho,i,j,k)/1.62D-9)), 0.0D0)
                  phi_r = first_derivative(x(i-1), x(i), x(i+1), phi(i-1,j,k), phi(i,j,k), phi(i+1,j,k))
                  phi_z = first_derivative(z(k-1), z(k), z(k+1), phi(i,j,k-1), phi(i,j,k), phi(i,j,k+1))
                  g_eff = at_no * DSQRT(phi_r**2 + phi_z**2)
	               turb_RT(i,j,k) = 0.625D0 * prim(irho,i,j,k) * DSQRT(DSQRT(prim(iturbq,i,j,k)**2)) * g_eff
               endif
            ENDIF

         enddo
      enddo
   enddo


   IF (turbcheck_flag == 1) THEN
      CLOSE(127)
      OPEN (UNIT = 125, FILE = './turb_eps'// trim(adjustl(globalt)) //'.dat', STATUS = 'REPLACE')
      DO k = 1, nz, 1
         DO i = 1, nx, 1
            WRITE(125,*) turb_eps(i,1,k)
         ENDDO
      ENDDO
      CLOSE(125)
      
      OPEN (UNIT = 126, FILE = './turb_RT'// trim(adjustl(globalt)) //'.dat', STATUS = 'REPLACE')
      DO k = 1, nz, 1
         DO i = 1, nx, 1
            WRITE(126,*) turb_RT(i,1,k)
         ENDDO
      ENDDO
      CLOSE(126)
   ENDIF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do k = 1, nz, 1
      do j = 1, ny, 1
         do i = 1, nx, 1

            dvxdx = first_derivative (x(i-1), x(i), x(i+1), prim(ivx,i-1,j,k), prim(ivx,i,j,k), prim(ivx,i+1,j,k))
            dvxdy = 0.0D0 !first_derivative (y(j-1), y(j), y(j+1), prim(ivx,i,j-1,k), prim(ivx,i,j,k), prim(ivx,i,j+1,k))
            dvxdz = first_derivative (z(k-1), z(k), z(k+1), prim(ivx,i,j,k-1), prim(ivx,i,j,k), prim(ivx,i,j,k+1))
            dvydx = first_derivative (x(i-1), x(i), x(i+1), prim(ivy,i-1,j,k), prim(ivy,i,j,k), prim(ivy,i+1,j,k))
            dvydy = 0.0D0 !first_derivative (y(j-1), y(j), y(j+1), prim(ivy,i,j-1,k), prim(ivy,i,j,k), prim(ivy,i,j+1,k))
            dvydz = first_derivative (z(k-1), z(k), z(k+1), prim(ivy,i,j,k-1), prim(ivy,i,j,k), prim(ivy,i,j,k+1))
            dvzdx = first_derivative (x(i-1), x(i), x(i+1), prim(ivz,i-1,j,k), prim(ivz,i,j,k), prim(ivz,i+1,j,k))
            dvzdy = 0.0D0 !first_derivative (y(j-1), y(j), y(j+1), prim(ivz,i,j-1,k), prim(ivz,i,j,k), prim(ivz,i,j+1,k))
            dvzdz = first_derivative (z(k-1), z(k), z(k+1), prim(ivz,i,j,k-1), prim(ivz,i,j,k), prim(ivz,i,j,k+1))

            turb_velgrad(i,j,k,1,1) = dvxdx
            turb_velgrad(i,j,k,1,2) = dvydx
            turb_velgrad(i,j,k,1,3) = dvzdx
            turb_velgrad(i,j,k,2,1) = 1.0D0/x(i)*(dvxdy-prim(ivy,i,j,k))
            turb_velgrad(i,j,k,2,2) = 1.0D0/x(i)*(dvydy+prim(ivx,i,j,k))
            turb_velgrad(i,j,k,2,3) = 1.0D0/x(i)*dvzdy 
            turb_velgrad(i,j,k,3,1) = dvxdz 
            turb_velgrad(i,j,k,3,2) = dvydz
            turb_velgrad(i,j,k,3,3) = dvzdz

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            turb_div_v(i,j,k) = prim(ivx,i,j,k)/x(i)+dvxdx + 1.0D0/x(i)*dvydy + dvzdz 

            turb_comp(i,j,k) = c_lambda * prim(irho,i,j,k) * prim(iturbq,i,j,k) * turb_div_v(i,j,k)

         enddo
      enddo
   enddo

   IF (turbcheck_flag == 1) THEN

      OPEN (UNIT = 123, FILE = './turb_comp'// trim(adjustl(globalt)) //'.dat', STATUS = 'REPLACE')
      DO k = 1, nz, 1
         DO i = 1, nx, 1
            WRITE(123,*) turb_comp(i,1,k)
         ENDDO
      ENDDO
      CLOSE(123)

   ENDIF

   do k = 1, nz, 1
      do j = 1, ny, 1
         do i = 1, nx, 1
            turb_str_tensor(i,j,k,:,:) = turb_velgrad(i,j,k,:,:) + TRANSPOSE(turb_velgrad(i,j,k,:,:)) - c_lambda*turb_div_v(i,j,k)*eye
         enddo
      enddo
   enddo

   do k = 1, nz, 1
      do j = 1, ny, 1
         do i = 1, nx, 1
            turb_str_tensor(i,j,k,:,1) = prim(irho,i,j,k)*turb_eta(i,j,k,1)*turb_str_tensor(i,j,k,:,1)
            turb_str_tensor(i,j,k,:,2) = prim(irho,i,j,k)*turb_eta(i,j,k,2)*turb_str_tensor(i,j,k,:,2)
            turb_str_tensor(i,j,k,:,3) = prim(irho,i,j,k)*turb_eta(i,j,k,3)*turb_str_tensor(i,j,k,:,3)
         enddo
      enddo
   enddo

   do k = 1, nz, 1
      do j = 1, ny, 1
         do i = 1, nx, 1
            sum = 0.0D0
            do m = 1, 3, 1
               do n = 1, 3, 1
                     sum = sum + turb_str_tensor(i,j,k,m,n) * turb_velgrad(i,j,k,m,n)
               enddo
            enddo
            turb_str(i,j,k) = sum
         enddo
      enddo
   enddo

   IF (turbcheck_flag == 1) THEN

      OPEN (UNIT = 124, FILE = './turb_str'// trim(adjustl(globalt)) //'.dat', STATUS = 'REPLACE')
      DO k = 1, nz, 1
         DO i = 1, nx, 1
            WRITE(124,*) turb_str(i,1,k)
         ENDDO
      ENDDO
      CLOSE(124)
   ENDIF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Find diffusion

   do k = 0, nz+1, 1
      do j = 0, ny+1, 1
         do i = 0, nx+1, 1
            dqdx = first_derivative (x(i-1), x(i), x(i+1), prim(iturbq,i-1,j,k), prim(iturbq,i,j,k), prim(iturbq,i+1,j,k))
            dqdy = 0.0D0 !first_derivative (y(j-1), y(j), y(j+1), prim(iturbq,i,j-1,k), prim(iturbq,i,j,k), prim(iturbq,i,j+1,k))
            dqdz = first_derivative (z(k-1), z(k), z(k+1), prim(iturbq,i,j,k-1), prim(iturbq,i,j,k), prim(iturbq,i,j,k+1))

            turb_qdev(i,j,k,1) = prim(irho,i,j,k)*turb_eta(i,j,k,1)*dqdx
            turb_qdev(i,j,k,2) = prim(irho,i,j,k)*turb_eta(i,j,k,2)*1.0D0/x(i)*dqdy
            turb_qdev(i,j,k,3) = prim(irho,i,j,k)*turb_eta(i,j,k,3)*dqdz

         enddo
      enddo
   enddo

   do k = 1, nz, 1
      do j = 1, ny, 1
         do i = 1, nx, 1

               dqdevdx = first_derivative (x(i-1), x(i), x(i+1), turb_qdev(i-1,j,k,1), turb_qdev(i,j,k,1), turb_qdev(i+1,j,k,1))
               dqdevdy = 0.0D0 !first_derivative (y(j-1), y(j), y(j+1), turb_qdev(i,j-1,k,2), turb_qdev(i,j,k,2), turb_qdev(i,j+1,k,2))
               dqdevdz = first_derivative (z(k-1), z(k), z(k+1), turb_qdev(i,j,k-1,3), turb_qdev(i,j,k,3), turb_qdev(i,j,k+1,3))

         enddo
      enddo
   enddo
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find source

do k = 1, nz, 1
   do j = 1, ny, 1
      do i = 1, nx, 1

         neighbour_flag = 1         
         turb_source(i,j,k) = 0.0D0
         turb_diff(i,j,k) = 0.0D0

         IF (turb_comp(i,j,k)*dt/(prim(irho,i,j,k)*prim_a(iturbq)) < 1e-2) THEN
            turb_comp(i,j,k) = 0.0D0
         ENDIF

         IF (turb_str(i,j,k)*dt/(prim(irho,i,j,k)*prim_a(iturbq)) < 1e-2) THEN
            turb_str(i,j,k) = 0.0D0
         ENDIF

         IF (turb_eps(i,j,k)*dt/(prim(irho,i,j,k)*prim_a(iturbq)) < 1e-2) THEN
            turb_eps(i,j,k) = 0.0D0
         ENDIF

         IF (turb_RT(i,j,k)*dt/(prim(irho,i,j,k)*prim_a(iturbq)) < 1e-2) THEN
            turb_RT(i,j,k) = 0.0D0
         ENDIF

         do m = -turb_neighbour, turb_neighbour, 1
            do n = -turb_neighbour, turb_neighbour, 1

               IF ( (prim(irho,i+m,j,k+n)-prim_a(irho))/prim_a(irho) < 1e-2 ) THEN
                  neighbour_flag = 0
                  EXIT
               ENDIF

            enddo
         enddo
         
         IF (neighbour_flag == 1) THEN
            ! Only calculate SGS turbulence for non-atmosphere terms whose neighbour also not atmosphere !
            turb_source(i,j,k) = -turb_comp(i,j,k) + turb_str(i,j,k) - turb_eps(i,j,k) + turb_RT(i,j,k)
            turb_diff(i,j,k) = turb_qdev(i,j,k,1)/x(i)+dqdevdx + 1.0D0/x(i)*dqdevdy + dqdevdz
         ENDIF

      enddo
   enddo
enddo

IF (turbcheck_flag == 1) THEN

   OPEN (UNIT = 128, FILE = './turb_diff'// trim(adjustl(globalt)) //'.dat', STATUS = 'REPLACE')
   DO k = 1, nz, 1
      DO i = 1, nx, 1
         WRITE(128,*) turb_diff(i,1,k)
      ENDDO
   ENDDO
   CLOSE(128)
ENDIF

CALL BOUNDARY1D_NM(turb_source, even, even, even, even, even, even)
CALL BOUNDARY1D_NM(turb_diff, even, even, even, even, even, even)
CALL BOUNDARY1D_NM(turb_eps, even, even, even, even, even, even)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DEALLOCATE(turb_RT)
DEALLOCATE(turb_comp)

DEALLOCATE(turb_div_v)
DEALLOCATE(turb_str)
DEALLOCATE(turb_eta)
DEALLOCATE(turb_velgrad)

DEALLOCATE(turb_qdev)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
   REAL*8 function first_derivative (xm1, xc, xp1, fm1, fc, fp1)
   !$acc routine seq
   implicit none
   REAL*8 :: xm1, xc, xp1, fm1, fc, fp1, h1, h2
   h2 = xp1 - xc
   h1 = xc - xm1
   first_derivative = ((fp1-fc)*h1*h1+(fc-fm1)*h2*h2)/(h1*h2*(h1+h2))
   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE findturbulence
