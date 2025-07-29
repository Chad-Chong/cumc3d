!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine allocates the neccessary arrays for 
! modeling level-sets.
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine buildLevelSet
use definition
USE CUSTOM_DEF
implicit none

! First allocate integer-type array
allocate(flamegrid_flag(-2:nx+3, -2:ny+3, -2:nz+3))
allocate(detongrid_flag(-2:nx+3, -2:ny+3, -2:nz+3))
allocate(flamecorn_flag(-2:nx+3, -2:ny+3, -2:nz+3))
allocate(detoncorn_flag(-2:nx+3, -2:ny+3, -2:nz+3))

! Then, allocate real-type array
allocate(flame_ratio(-2:nx+3, -2:ny+3, -2:nz+3))
allocate(deton_ratio(-2:nx+3, -2:ny+3, -2:nz+3))
allocate(flame_ratio_old(-2:nx+3, -2:ny+3, -2:nz+3))
allocate(deton_ratio_old(-2:nx+3, -2:ny+3, -2:nz+3))
allocate(burn_ratio(-2:nx+3, -2:ny+3, -2:nz+3))

! Extra array for NSE in HELMEOS_MODULE
allocate(flame_loc_ratio(-2:nx+3, -2:ny+3, -2:nz+3))
allocate(deton_loc_ratio(-2:nx+3, -2:ny+3, -2:nz+3))

! For Burning
ALLOCATE(nse_flag(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(nu_qdot(-2:nx+3,-2:ny+3,-2:nz+3)) 
ALLOCATE(burn_qdot(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(flame_qdot(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(deton_qdot(-2:nx+3,-2:ny+3,-2:nz+3))

end subroutine buildLevelSet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetFlame
USE definition
USE CUSTOM_DEF
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization of the scalar G

! Dummy variables
integer :: j , k

! Dist from the front
real*8 :: dist

! Initilization
flame_ratio_old = 0.0D0
flame_ratio = 0.0D0

! Assign initial flame
do j = 1, nx, 1   
    do k  = 1, nz, 1 

    ! big c3 flame
    prim(iscaG1,j,1,k) = 148.90D0 - DSQRT(x(j)**2 + z(k)**2) + &
                   60.92D0 * ABS(DSIN(ASIN(z(k) / DSQRT(x(j)**2 + z(k)**2)) * 6.0D0))

    ! c3 flame
        ! prim(iscaG1,j,1,k) = 74.45D0 - DSQRT(x(j)**2 + z(k)**2) + & 
        !             30.46D0 * ABS(DSIN(ASIN(z(k) / DSQRT(x(j)**2 + z(k)**2)) * 6.0D0))

    ! small c3 flame
    !prim(iscaG1,j,1,k) = 37.23D0 - DSQRT(x(j)**2 + z(k)**2) + &
        !            15.23D0 * ABS(DSIN(ASIN(z(k) / DSQRT(x(j)**2 + z(k)**2)) * 6.0D0))

    ! b1 flame 23.93 -> 50 km, 28.72 -> 60 km
    !prim(iscaG1,j,1,k) = -DSQRT((x(j) - 47.86D0)**2 + (z(k) - 47.86D0)**2) + 7.0D0

    ! b1 flame 23.93 -> 50 km, 28.72 -> 60 km
        !prim(iscaG1,j,1,k) = -DSQRT((x(j))**2 + (z(k) - 101.52D0)**2) + 15.0D0

    ! A bubble along z-axis in the helium sphere
    !prim(iscaG1,j,1,k) = -DSQRT((x(j))**2 + (z(k) - 1.0D3)**2) + 15.0D0

    ! Two bubbles along two axis
    

    ! b5 flame Type A
        !prim(iscaG1,j,1,k) = MAX(10.0D0 - DSQRT((DBLE(j-1)*dx(j) - 20.0D0)**2 + (DBLE(k-1)*dx(j) - 40.0D0)**2), &      
        !               10.0D0 - DSQRT((DBLE(j-1)*dx(j) - 40.0D0)**2 + (DBLE(k-1)*dx(j) - 20.0D0)**2), &
        !               10.0D0 - DSQRT((DBLE(j-1)*dx(j) - 25.0D0)**2 + (DBLE(k-1)*dx(j) - 75.0D0)**2), &
        !               10.0D0 - DSQRT((DBLE(j-1)*dx(j) - 50.0D0)**2 + (DBLE(k-1)*dx(j) - 50.0D0)**2), &
        !               10.0D0 - DSQRT((DBLE(j-1)*dx(j) - 75.0D0)**2 + (DBLE(k-1)*dx(j) - 25.0D0)**2))

    ! b5 flame Type B
    !prim(iscaG1,j,1,k) = MAX(10.0D0 - DSQRT((DBLE(j-1)*dx(j) - 33.84D0)**2 + (DBLE(k-1)*dx(j) - 67.68D0)**2), &      
        !               10.0D0 - DSQRT((DBLE(j-1)*dx(j) - 67.68D0)**2 + (DBLE(k-1)*dx(j) - 33.84D0)**2), &
        !               10.0D0 - DSQRT((DBLE(j-1)*dx(j) - 51.80D0)**2 + (DBLE(k-1)*dx(j) - 125.07D0)**2), &
        !               10.0D0 - DSQRT((DBLE(j-1)*dx(j) - 95.72D0)**2 + (DBLE(k-1)*dx(j) - 95.72D0)**2), &
        !               10.0D0 - DSQRT((DBLE(j-1)*dx(j) - 125.07D0)**2 + (DBLE(k-1)*dx(j) - 51.80D0)**2))

    enddo
enddo 

! Copy the data to ghost cells
CALL boundary1D_NM(prim(iscaG1,:,:,:), even, even, even, even, even, even)

! Readjust the initial data
call reinitialization2()               !You need this for aspherical initial flame
call identify_flamegrid()
call compute_flameratio()

! Backup initial result
flame_ratio_old = flame_ratio

! Set the deton part
deton_ratio(:,:,:) = 0.0D0
deton_ratio_old(:,:,:) = 0.0D0

! Sum all level-set ratio up
burn_ratio(:,:,:) = flame_ratio(:,:,:) + deton_ratio(:,:,:) 

end subroutine GetFlame

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine do the reinitialization in order to 
! maintain the distance property of the 2nd level set
! Based on the level set paper Reinecke (1999a)
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!

subroutine reinitialization2
use DEFINITION
USE CUSTOM_DEF
implicit none

! Dummy variables
integer :: j, k

! Delta and d0
real*8 :: delta, d0

! The number of point intersected by the grid and the front
integer :: flame_count

! The scaled distance
real*8 :: H_d

! The signed distance
real*8 :: sgn_scaG

! Array storing the minimal distance of the grid to the nearest front
real*8, dimension(-2:nx+3, -2:nz+3) :: flame_distance

! Character for debug
character(len=99) :: globalt

! Array storing all intersection points
real*8, dimension (nx * nz, 2):: flame_position

! First find all intersections
call locate_flame_grid(1, flame_count, flame_position)

! write(globalt,'(I)') n_step
! IF (MOD(n_step, 10) == 0) THEN
!     OPEN (UNIT = 123, FILE = './flame_posit'// trim(adjustl(globalt)) //'.dat', STATUS = 'REPLACE')
! ENDIF

!     DO k = 1, flame_count, 1

!         IF (MOD(n_step, 10) == 0) THEN
!             WRITE(123, *) flame_position(k,1), flame_position(k,2)
!         ENDIF

!     ENDDO
! IF (MOD(n_step, 10) == 0) THEN
!     CLOSE(123)
! ENDIF

! Then find the minimal distance
call locate_min_flame_distance(flame_count, flame_position, flame_distance)

! write(globalt,'(I)') n_step
! IF (MOD(n_step, 10) == 0) THEN
!     OPEN (UNIT = 123, FILE = './flame_distance'// trim(adjustl(globalt)) //'.dat', STATUS = 'REPLACE')
!     OPEN (UNIT = 124, FILE = './scaG1_before'// trim(adjustl(globalt)) //'.dat', STATUS = 'REPLACE')
! ENDIF
!     DO k = 1, nz, 1 
!         DO j = 1, nx, 1

!             IF (MOD(n_step, 10) == 0) THEN
!                 WRITE(123, *) flame_distance(j,k)
!                 WRITE(124,*) prim(iscaG1,j,1,k)
!             ENDIF

!         ENDDO
!     ENDDO
! IF (MOD(n_step, 10) == 0) THEN
!     CLOSE(123)
!     CLOSE(124)
! ENDIF


d0 = 3.0D0*dx(1)
delta = 1.0D0*dx(1)

! write(globalt,'(I)') n_step
! IF (MOD(n_step, 10) == 0) THEN
!     OPEN (UNIT = 123, FILE = './H_d'// trim(adjustl(globalt)) //'.dat', STATUS = 'REPLACE')
!     OPEN (UNIT = 124, FILE = './scaG1_after'// trim(adjustl(globalt)) //'.dat', STATUS = 'REPLACE')
!     OPEN (UNIT = 125, FILE = './prim'// trim(adjustl(globalt)) //'.dat', STATUS = 'REPLACE')
! ENDIF

do k = 1, nz, 1
    do j = 1, nx, 1

    ! Correct the level set
        H_d = (1.0D0 - TANH(3.0D0*(flame_distance(j,k) - d0)/delta)) / & 
            (1.0D0 - TANH(-3.0D0*d0/delta))

        if(prim(iscaG1,j,1,k) >= 0.0D0) then
        sgn_scaG = 1.0D0
        else
        sgn_scaG = -1.0D0
        endif
        prim(iscaG1,j,1,k) = H_d * prim(iscaG1,j,1,k) + (1.0D0 - H_d) * sgn_scaG * flame_distance(j,k)
        
        ! IF (MOD(n_step, 10) == 0) THEN
        !     WRITE(123, *) flame_distance(j,k), H_d
        !     WRITE(124,*) prim(iscaG1,j,1,k)
        !     WRITE(125,*) prim(ic12,j,1,k)
        ! ENDIF

    enddo
enddo

! IF (MOD(n_step, 10) == 0) THEN
!     CLOSE(123)
!     CLOSE(124)
!     CLOSE(125)
! ENDIF

! Copy the results to ghost cells
CALL BOUNDARY1D_NM(prim(iscaG1,:,:,:), even, even, even, even, even, even)  

100 FORMAT (6E13.6)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine do the reinitialization in order to 
! maintain the distance property of the 2nd level set
! Based on the level set paper Reinecke (1999a)
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine set the reinitialization, using grid as counter

subroutine reinitialization3
use definition
use custom_def
implicit none

! Dummy variables                 
integer :: j, k          

! Delta and d0
real*8 :: delta, d0

! The number of point intersected by the grid and the front
integer :: flame_count

! The scaled distance
real*8 :: H_d

! The signed distance
real*8 :: sgn_scaG            

! Array storing the minimal distance of the grid to the nearest front
real*8, dimension(-2:nx+3, -2:nz+3) :: flame_distance

! Array storing all intersection points
real*8, dimension(nx * nz, 2):: flame_position
            
! First find all intersections         
call locate_flame_grid(2, flame_count, flame_position)

! Then find the minimal distance
call locate_min_flame_distance(flame_count, flame_position, flame_distance)

d0 = 3.0D0*dx(1)
delta = 1.0D0*dx(1)


do k = 1, nz, 1
    do j = 1, nx, 1  
    
    ! Correct the level set
        H_d = (1.0D0 - TANH(3.0D0*(flame_distance(j,k) - d0)/delta)) / & 
            (1.0D0 - TANH(-3.0D0*d0/delta))
        if(prim(iscaG2,j,1,k) >= 0.0D0) then
            sgn_scaG = 1.0D0
        else
            sgn_scaG = -1.0D0
        endif
            prim(iscaG2,j,1,k) = H_d * prim(iscaG2,j,1,k) + (1.0D0 - H_d) * sgn_scaG * flame_distance(j,k)

    enddo
enddo


! Copy the results to the ghost cells
CALL BOUNDARY1D_NM(prim(iscaG2,:,:,:), even, even, even, even, even, even)  

100 FORMAT (6E13.6)

end subroutine reinitialization3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine find all points where the 
! deflagration/detonation front cuts the grid
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutines find all the intersections points of the 
! front and the grid. The subroutine takes input of the mode
! Mode 1 = first level set
! Mode 2 = second level set
! and gives the number of intersection point and their
! corresponding positions. 
!
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine locate_flame_grid(mode, flame_count, flame_position)
use definition
USE CUSTOM_DEF
implicit none

!dummy variables
integer :: j, k

! Output intersection point number
integer :: flame_count

! Input mode
integer :: mode

! distance of grid to the intersection
real*8 :: flame_dx

! Output array of intersection point
real*8, dimension(nx*nz, 2) :: flame_position

! Initialization
flame_count = 0 
flame_rad = 0.0D0

!write(*,*) 'In locate flame grid'

if(mode == 1) then

    ! For the first level set
    do k = 1, nz, 1
        do j = 1, nx, 1

        if(prim(iscaG1,j,1,k) * prim(iscaG1,j+1,1,k) < 0.0D0) then
            flame_count = flame_count + 1
            flame_dx = ABS(prim(iscaG1,j,1,k)) / ABS(prim(iscaG1,j+1,1,k) - prim(iscaG1,j,1,k))
            flame_position(flame_count, 1) = x(j) + flame_dx * dx(j)
            flame_position(flame_count, 2) = z(k)
            if(flame_position(flame_count, 1) > flame_rad) flame_rad = flame_position(flame_count, 1)
        endif

        if(prim(iscaG1,j,1,k) * prim(iscaG1,j,1,k+1) < 0.0D0) then             
            flame_count = flame_count + 1
            flame_dx = ABS(prim(iscaG1,j,1,k)) / ABS(prim(iscaG1,j,1,k+1) - prim(iscaG1,j,1,k))
            flame_position(flame_count, 1) = x(j)
            flame_position(flame_count, 2) = z(k) + flame_dx * dz(k)
            if(flame_position(flame_count, 1) > flame_rad) flame_rad = flame_position(flame_count, 1) 
        endif

        enddo
    enddo

elseif(mode == 2) then

    ! For the second level set
    do k = 1, nz, 1
        do j = 1, nx, 1

            if(prim(iscaG2,j,1,k) * prim(iscaG2,j+1,1,k) < 0.0D0) then
                flame_count = flame_count + 1
                flame_dx = ABS(prim(iscaG2,j,1,k)) / ABS(prim(iscaG2,j+1,1,k) - prim(iscaG2,j,1,k))
                flame_position(flame_count, 1) = x(j) + flame_dx * dx(j)
                flame_position(flame_count, 2) = z(k)
                if(flame_position(flame_count, 1) > flame_rad) flame_rad = flame_position(flame_count, 1)
            endif

            if(prim(iscaG2,j,1,k) * prim(iscaG2,j,1,k+1) < 0.0D0) then
                flame_count = flame_count + 1
                flame_dx = ABS(prim(iscaG2,j,1,k)) / ABS(prim(iscaG2,j,1,k+1) - prim(iscaG2,j,1,k))
                flame_position(flame_count, 1) = x(j) 
                flame_position(flame_count, 2) = z(k) + flame_dx * dz(k)
                if(flame_position(flame_count, 1) > flame_rad) flame_rad = flame_position(flame_count, 1)
            endif

        enddo
    enddo

endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the minimum distance from that
! particular grid point to the flame surface.
! The subroutines take input of the intersection
! point number, their positions. Then it gives
! the minimum distances from all grid point to 
! the surface
!
! Written by Leung Shing Chi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine locate_min_flame_distance(flame_count, flame_position, flame_distance)
use definition
USE CUSTOM_DEF
implicit none

! dummy variables
integer :: j, k, k2

! Input intersection point number
integer :: flame_count

! The local minimal distance
real (selected_real_kind(15,307)):: distance, last_distance

! Output array for the minimal distance
real (selected_real_kind(15,307)), dimension(-2:nx+3, -2:nz+3):: flame_distance

! Input array of intersection point positions
real (selected_real_kind(15,307)), dimension(nx * nz, 2) :: flame_position


do j = 1, nx, 1
    do k = 1, nz, 1

        last_distance = 10.0D0 * DBLE(MAX(nx,nz)) * MAX(dx(j), dz(k))

        ! Search for the minimal distance
        do k2 = 1, flame_count, 1

            distance = DSQRT((flame_position(k2,1) - x(j)) ** 2 + & 
                        (flame_position(k2,2) - z(k)) ** 2) 

            if((flame_count > 1 .and. distance <= last_distance) .or. flame_count == 1) then
                flame_distance(j,k) = distance 
                last_distance = distance
                endif

            if(k2 == 1) last_distance = distance
        enddo

    enddo
enddo


!WRITE(*,*) flame_distance(1,1), flame_count

100 FORMAT (6E13.6)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine classifies the type of grid
! (filled or unfilled?) / (partially or completely?)
!
! The classification is defined as
! 
! 0 = completely unburnt
! 1 = completely burnt
! 2 = 
!
! The data transfer includes:
! In: scaG_in (the level-set field)
! Out: 
!  
!
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine identify_flamegrid()
use definition
USE CUSTOM_DEF
implicit none

! Dummy variables 
integer :: j, k

! Initilization
flamecorn_flag = 0

do j = 0, nx, 1
    do k = 0, nz, 1

        ! The technique is as follows
        ! If all corners are positive/negative, then it is either completely occpied/emptied
        ! If any one or above corners have different signs, that means the front cut this grid
        ! which means you need to calculate what fraction this front occupies the grid 

        if(prim(iscaG1,j,1,k) > 0.0D0 .or. prim(iscaG1,j+1,1,k) > 0.0D0 .or. prim(iscaG1,j,1,k+1) > 0.0D0 .or. prim(iscaG1,j+1,1,k+1) > 0.0D0) then
        if(prim(iscaG1,j,1,k) > 0.0D0 .and. prim(iscaG1,j+1,1,k) > 0.0D0 .and. prim(iscaG1,j,1,k+1) > 0.0D0 .and. prim(iscaG1,j+1,1,k+1) > 0.0D0) then

        ! Completely occupied
        flamecorn_flag(j,1,k) = 1

        elseif(prim(iscaG1,j,1,k) * prim(iscaG1,j+1,1,k) * prim(iscaG1,j,1,k+1) * prim(iscaG1,j+1,1,k+1) > 0.0D0) then

        if((prim(iscaG1,j,1,k) * prim(iscaG1,j,1,k+1) <= 0.0D0) .and. (prim(iscaG1,j,1,k) * prim(iscaG1,j+1,1,k) <= 0.0D0)) then
        ! The ambiguous case
            flamecorn_flag(j,1,k) = 4
        else
        ! The parallelogram case
            flamecorn_flag(j,1,k) = 3
        endif

        else

            ! The triangle case
            flamecorn_flag(j,1,k) = 2

        endif

            elseif(prim(iscaG1,j,1,k) <= 0.0D0 .and. prim(iscaG1,j+1,1,k) <= 0.0D0 .and. &
                prim(iscaG1,j,1,k+1) <= 0.0D0 .and. prim(iscaG1,j+1,1,k+1) <= 0.0D0) then

        ! Empty case
            flamecorn_flag(j,1,k) = 0

            endif

    enddo
enddo

100 FORMAT (6I5)
101 FORMAT (6E13.6)


end subroutine identify_flamegrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine classifies the type of grid 
! (filled or unfilled?) / (partially or completely?)
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine identify_detongrid()
use definition
use custom_def
implicit none

! dummy variables
integer :: j, k  

! initialization
detoncorn_flag = 0

do k = 0, nz, 1   
    do j = 0, nx, 1

        if(prim(iscaG2,j,1,k) > 0.0D0 .or. prim(iscaG2,j+1,1,k) > 0.0D0 .or. prim(iscaG2,j,1,k+1) > 0.0D0 .or. prim(iscaG2,j+1,1,k+1) > 0.0D0) then

        if(prim(iscaG2,j,1,k) > 0.0D0 .and. prim(iscaG2,j+1,1,k) > 0.0D0 .and. prim(iscaG2,j,1,k+1) > 0.0D0 .and. prim(iscaG2,j+1,1,k+1) > 0.0D0) then
        
        ! Completely filled
            detoncorn_flag(j,1,k) = 1

        elseif(prim(iscaG2,j,1,k) * prim(iscaG2,j+1,1,k) * prim(iscaG2,j,1,k+1) * prim(iscaG2,j+1,1,k+1) > 0.0D0) then

            if((prim(iscaG2,j,1,k) * prim(iscaG2,j,1,k+1) <= 0.0D0) .and. (prim(iscaG2,j,1,k) * prim(iscaG2,j+1,1,k) <= 0.0D0)) then

        ! The ambigous case
                detoncorn_flag(j,1,k) = 4

            else

        ! The parallelogram case
                detoncorn_flag(j,1,k) = 3

            endif

        else

        ! The triangle case
            detoncorn_flag(j,1,k) = 2

        endif

        elseif(prim(iscaG2,j,1,k) <= 0.0D0 .and. prim(iscaG2,j+1,1,k) <= 0.0D0 .and. &
            prim(iscaG2,j,1,k+1) <= 0.0D0 .and. prim(iscaG2,j+1,1,k+1) <= 0.0D0) then

    ! The empty case
        detoncorn_flag(j,1,k) = 0

        endif
    enddo
enddo

!CALL BOUNDARY1D_INT(detongrid_flag, even)

100 FORMAT (6I5)
101 FORMAT (6E13.6)

end subroutine identify_detongrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine calculates the fraction occupied by one of
! the phase (ash if in the context of supernova)
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_flameratio()
use definition
use custom_def
implicit none

! Dummy variables
integer :: j, k

! Local storage
real*8 :: loc_flame_ratio 

!I need to include the opposite case as well!

! Initialization
!write(*,*) 'In compute norm snd flame ratio'
flame_ratio = 0.0D0

do k = 0, nz, 1
    do j = 0, nx, 1

        if(flamecorn_flag(j,1,k) /= 0) then

    if(flamecorn_flag(j,1,k) == 1) then

        ! Completely filled is one	       
        loc_flame_ratio = 1.0D0

        elseif(flamecorn_flag(j,1,k) == 2) then

        ! Incomplete filled needs to be studied case by case

        if(prim(iscaG1,j,1,k) >= 0.0D0 .and. prim(iscaG1,j,1,k+1) < 0.0D0 .and. &
                prim(iscaG1,j+1,1,k) < 0.0D0 .and. prim(iscaG1,j+1,1,k+1) < 0.0D0) then
                loc_flame_ratio = ((prim(iscaG1,j,1,k)/(prim(iscaG1,j,1,k)-prim(iscaG1,j+1,1,k))) * & 
                                    (prim(iscaG1,j,1,k)/(prim(iscaG1,j,1,k)-prim(iscaG1,j,1,k+1)))) * 0.5D0
            elseif(prim(iscaG1,j,1,k) < 0.0D0 .and. prim(iscaG1,j,1,k+1) >= 0.0D0 .and. &
                prim(iscaG1,j+1,1,k) >= 0.0D0 .and. prim(iscaG1,j+1,1,k+1) >= 0.0D0) then
                loc_flame_ratio = 1.0D0 - ((prim(iscaG1,j,1,k)/(prim(iscaG1,j,1,k)-prim(iscaG1,j+1,1,k))) * &
                                    (prim(iscaG1,j,1,k)/(prim(iscaG1,j,1,k)-prim(iscaG1,j,1,k+1)))) * 0.5D0

            elseif(prim(iscaG1,j,1,k) < 0.0D0 .and. prim(iscaG1,j,1,k+1) >= 0.0D0 .and. &
                    prim(iscaG1,j+1,1,k) < 0.0D0 .and. prim(iscaG1,j+1,1,k+1) < 0.0D0) then
                loc_flame_ratio = ((prim(iscaG1,j,1,k+1)/(prim(iscaG1,j,1,k+1)-prim(iscaG1,j+1,1,k+1))) * &
                                    (prim(iscaG1,j,1,k+1)/(prim(iscaG1,j,1,k+1)-prim(iscaG1,j,1,k)))) * 0.5D0
            elseif(prim(iscaG1,j,1,k) >= 0.0D0 .and. prim(iscaG1,j,1,k+1) < 0.0D0 .and. &
                    prim(iscaG1,j+1,1,k) >= 0.0D0 .and. prim(iscaG1,j+1,1,k+1) >= 0.0D0) then
                loc_flame_ratio = 1.0D0 - ((prim(iscaG1,j,1,k+1)/(prim(iscaG1,j,1,k+1)-prim(iscaG1,j+1,1,k+1))) * &
                                    (prim(iscaG1,j,1,k+1)/(prim(iscaG1,j,1,k+1)-prim(iscaG1,j,1,k)))) * 0.5D0

            elseif(prim(iscaG1,j,1,k) < 0.0D0 .and. prim(iscaG1,j,1,k+1) < 0.0D0 .and. &
                    prim(iscaG1,j+1,1,k) >= 0.0D0 .and. prim(iscaG1,j+1,1,k+1) < 0.0D0) then
                loc_flame_ratio = ((prim(iscaG1,j+1,1,k)/(prim(iscaG1,j+1,1,k)-prim(iscaG1,j+1,1,k+1))) * &
                                    (prim(iscaG1,j+1,1,k)/(prim(iscaG1,j+1,1,k)-prim(iscaG1,j,1,k)))) * 0.5D0    
            elseif(prim(iscaG1,j,1,k) >= 0.0D0 .and. prim(iscaG1,j,1,k+1) > 0.0D0 .and. &
                    prim(iscaG1,j+1,1,k) < 0.0D0 .and. prim(iscaG1,j+1,1,k+1) > 0.0D0) then
                loc_flame_ratio = 1.0D0 - ((prim(iscaG1,j+1,1,k)/(prim(iscaG1,j+1,1,k)-prim(iscaG1,j+1,1,k+1))) * &
                                    (prim(iscaG1,j+1,1,k)/(prim(iscaG1,j+1,1,k)-prim(iscaG1,j,1,k)))) * 0.5D0

            elseif(prim(iscaG1,j,1,k) < 0.0D0 .and. prim(iscaG1,j,1,k+1) < 0.0D0 .and. &
                    prim(iscaG1,j+1,1,k) < 0.0D0 .and. prim(iscaG1,j+1,1,k+1) >= 0.0D0) then
                loc_flame_ratio = ((prim(iscaG1,j+1,1,k+1)/(prim(iscaG1,j+1,1,k+1)-prim(iscaG1,j+1,1,k))) * &
                                    (prim(iscaG1,j+1,1,k+1)/(prim(iscaG1,j+1,1,k+1)-prim(iscaG1,j,1,k+1)))) * 0.5D0
            elseif(prim(iscaG1,j,1,k) >= 0.0D0 .and. prim(iscaG1,j,1,k+1) >= 0.0D0 .and. &
                    prim(iscaG1,j+1,1,k) >= 0.0D0 .and. prim(iscaG1,j+1,1,k+1) < 0.0D0) then
                loc_flame_ratio = 1.0D0 - ((prim(iscaG1,j+1,1,k+1)/(prim(iscaG1,j+1,1,k+1)-prim(iscaG1,j+1,1,k))) * &
                                    (prim(iscaG1,j+1,1,k+1)/(prim(iscaG1,j+1,1,k+1)-prim(iscaG1,j,1,k+1)))) * 0.5D0
            endif

    elseif(flamecorn_flag(j,1,k) == 3) then

            ! Parallelogram case

            if(prim(iscaG1,j,1,k) >= 0.0D0 .and. prim(iscaG1,j,1,k+1) < 0.0D0) then
                loc_flame_ratio = ((prim(iscaG1,j,1,k)/(prim(iscaG1,j,1,k)-prim(iscaG1,j,1,k+1))) + &
                                    (prim(iscaG1,j+1,1,k)/(prim(iscaG1,j+1,1,k)-prim(iscaG1,j+1,1,k+1)))) * 0.5D0  
            elseif(prim(iscaG1,j,1,k) < 0.0D0 .and. prim(iscaG1,j,1,k+1) >= 0.0D0) then      
                loc_flame_ratio = ((prim(iscaG1,j,1,k+1)/(prim(iscaG1,j,1,k+1)-prim(iscaG1,j,1,k))) + &
                                    (prim(iscaG1,j+1,1,k+1)/(prim(iscaG1,j+1,1,k+1)-prim(iscaG1,j+1,1,k)))) * 0.5D0

            elseif(prim(iscaG1,j,1,k) >= 0.0D0 .and. prim(iscaG1,j+1,1,k) < 0.0D0) then
                loc_flame_ratio = ((prim(iscaG1,j,1,k)/(prim(iscaG1,j,1,k)-prim(iscaG1,j+1,1,k))) + &
                                    (prim(iscaG1,j,1,k+1)/(prim(iscaG1,j,1,k+1)-prim(iscaG1,j+1,1,k+1)))) * 0.5D0
            elseif(prim(iscaG1,j,1,k) < 0.0D0 .and. prim(iscaG1,j+1,1,k) >= 0.0D0) then
                loc_flame_ratio = ((prim(iscaG1,j+1,1,k)/(prim(iscaG1,j+1,1,k)-prim(iscaG1,j,1,k))) + &  
                                    (prim(iscaG1,j+1,1,k+1)/(prim(iscaG1,j+1,1,k+1)-prim(iscaG1,j,1,k+1)))) * 0.5D0

            endif
        
        elseif(flamecorn_flag(j,1,k) == 4) then

        ! The ambiguous case averages all possibilities

        if(prim(iscaG1,j,1,k) > 0.0D0) then
                loc_flame_ratio = 0.5D0 * &
                                    ((((prim(iscaG1,j,1,k)/(prim(iscaG1,j,1,k) - prim(iscaG1,j,1,k+1))) * &
                                    (prim(iscaG1,j,1,k)/(prim(iscaG1,j,1,k) - prim(iscaG1,j+1,1,k)))) * 0.5D0 + &
                                    ((prim(iscaG1,j+1,1,k+1)/(prim(iscaG1,j+1,1,k+1)-prim(iscaG1,j+1,1,k))) * &
                                    (prim(iscaG1,j+1,1,k+1)/(prim(iscaG1,j+1,1,k+1)-prim(iscaG1,j,1,k+1)))) * 0.5D0) + &
                                    (1.0D0 - ((prim(iscaG1,j,1,k+1)/(prim(iscaG1,j,1,k+1) - prim(iscaG1,j,1,k))) * &
                                    (prim(iscaG1,j,1,k+1)/(prim(iscaG1,j,1,k+1) - prim(iscaG1,j+1,1,k+1)))) * 0.5D0 - &  
                                    ((prim(iscaG1,j+1,1,k)/(prim(iscaG1,j+1,1,k)-prim(iscaG1,j,1,k))) * &
                                    (prim(iscaG1,j+1,1,k)/(prim(iscaG1,j+1,1,k)-prim(iscaG1,j+1,1,k+1)))) * 0.5D0))
            else
                loc_flame_ratio = 0.5D0 * &
                                    ((1.0D0 - ((prim(iscaG1,j,1,k)/(prim(iscaG1,j,1,k) - prim(iscaG1,j,1,k+1))) * &
                                    (prim(iscaG1,j,1,k)/(prim(iscaG1,j,1,k) - prim(iscaG1,j+1,1,k)))) * 0.5D0 - &
                                    ((prim(iscaG1,j+1,1,k+1)/(prim(iscaG1,j+1,1,k+1)-prim(iscaG1,j+1,1,k))) * &
                                    (prim(iscaG1,j+1,1,k+1)/(prim(iscaG1,j+1,1,k+1)-prim(iscaG1,j,1,k+1)))) * 0.5D0) + &
                                    (((prim(iscaG1,j,1,k+1)/(prim(iscaG1,j,1,k+1) - prim(iscaG1,j,1,k))) * &  
                                    (prim(iscaG1,j,1,k+1)/(prim(iscaG1,j,1,k+1) - prim(iscaG1,j+1,1,k+1)))) * 0.5D0 + & 
                                    ((prim(iscaG1,j+1,1,k)/(prim(iscaG1,j+1,1,k)-prim(iscaG1,j,1,k))) * &             
                                    (prim(iscaG1,j+1,1,k)/(prim(iscaG1,j+1,1,k)-prim(iscaG1,j+1,1,k+1)))) * 0.5D0))
            endif
        
    endif

        else

    ! Completely empty is zero
        loc_flame_ratio = 0.0D0

        endif

    ! Notice that the ratio is defined with grid corner
    ! as center, so the contribution of the fraction is 
    ! applied to all the overlapped grids
    loc_flame_ratio = 0.25D0 * loc_flame_ratio

    flame_ratio(j,1,k) = flame_ratio(j,1,k) + loc_flame_ratio
    flame_ratio(j+1,1,k) = flame_ratio(j+1,1,k) + loc_flame_ratio
    flame_ratio(j,1,k+1) = flame_ratio(j,1,k+1) + loc_flame_ratio
    flame_ratio(j+1,1,k+1) = flame_ratio(j+1,1,k+1) + loc_flame_ratio

    enddo
enddo

! Copy the result to the current fraction
flame_loc_ratio(:,:,:) = flame_ratio(:,:,:)

! According to the summed fraction, find the 
! grid type according to the geometry
do k = 1, nz, 1
    do j = 1, nx, 1

        ! Note: For this case, only three types
        ! of grid are needed because we have 
        ! already found the fraction	 

        ! Get the effective flame ratio
            !flame_ratio(j,1,k) = MIN(MAX(flame_ratio(j,1,k), flame_ratio_old(j,1,k)), 1.0D0 - deton_ratio(j,1,k))
        flame_ratio(j,1,k) = MAX(flame_ratio(j,1,k), flame_ratio_old(j,1,k))
        !flame_ratio(j,1,k) = MIN(flame_ratio(j,1,k), 1.0D0 - deton_ratio(j,1,k))

        if(flame_ratio(j,1,k) == 0.0D0) then

            ! Completely empty
            flamegrid_flag(j,1,k) = 0

            elseif(flame_ratio(j,1,k) == 1.0D0) then

            ! Completely filled
            flamegrid_flag(j,1,k) = 1

            else

            ! Partially filled	
            flamegrid_flag(j,1,k) = 2

        endif

    enddo
enddo

! Copy the results to ghost cells
CALL BOUNDARY1D_NM(flame_ratio, even, even, even, even, even, even)
!CALL BOUNDARY1D_INT(flamegrid_flag,even)

100 FORMAT (6E13.6)

end subroutine compute_flameratio

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine calculates the fraction occupied by one of
! the phase (ash if in the context of supernova)
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_detonratio()
use definition
use custom_def
implicit none

! Dummy variables         
integer :: j, k

! Local stroage
real*8 :: loc_deton_ratio
        
!I need to include the opposite case as well!
!write(*,*) 'In compute norm snd flame ratio'
deton_ratio(:,:,:) = 0.0D0
    
do k = 0, nz, 1
    do j = 0, nx, 1

        if(detoncorn_flag(j,1,k) /= 0) then

        if(detoncorn_flag(j,1,k) == 1) then

        ! Completely filled is one
            loc_deton_ratio = 1.0D0

        elseif(detoncorn_flag(j,1,k) == 2) then

        ! For triangle case, need to study case by case

        if(prim(iscaG2,j,1,k) >= 0.0D0 .and. prim(iscaG2,j,1,k+1) < 0.0D0 .and. &
                prim(iscaG2,j+1,1,k) < 0.0D0 .and. prim(iscaG2,j+1,1,k+1) < 0.0D0) then
                loc_deton_ratio = ((prim(iscaG2,j,1,k)/(prim(iscaG2,j,1,k)-prim(iscaG2,j+1,1,k))) * &
                                    (prim(iscaG2,j,1,k)/(prim(iscaG2,j,1,k)-prim(iscaG2,j,1,k+1)))) * 0.5D0
            elseif(prim(iscaG2,j,1,k) < 0.0D0 .and. prim(iscaG2,j,1,k+1) >= 0.0D0 .and. &     
                prim(iscaG2,j+1,1,k) >= 0.0D0 .and. prim(iscaG2,j+1,1,k+1) >= 0.0D0) then
                loc_deton_ratio = 1.0D0 - ((prim(iscaG2,j,1,k)/(prim(iscaG2,j,1,k)-prim(iscaG2,j+1,1,k))) * &
                                    (prim(iscaG2,j,1,k)/(prim(iscaG2,j,1,k)-prim(iscaG2,j,1,k+1)))) * 0.5D0
                                
            elseif(prim(iscaG2,j,1,k) < 0.0D0 .and. prim(iscaG2,j,1,k+1) >= 0.0D0 .and. &
                    prim(iscaG2,j+1,1,k) < 0.0D0 .and. prim(iscaG2,j+1,1,k+1) < 0.0D0) then
                loc_deton_ratio = ((prim(iscaG2,j,1,k+1)/(prim(iscaG2,j,1,k+1)-prim(iscaG2,j+1,1,k+1))) * &
                                    (prim(iscaG2,j,1,k+1)/(prim(iscaG2,j,1,k+1)-prim(iscaG2,j,1,k)))) * 0.5D0
            elseif(prim(iscaG2,j,1,k) >= 0.0D0 .and. prim(iscaG2,j,1,k+1) < 0.0D0 .and. &   
                    prim(iscaG2,j+1,1,k) >= 0.0D0 .and. prim(iscaG2,j+1,1,k+1) >= 0.0D0) then
                loc_deton_ratio = 1.0D0 - ((prim(iscaG2,j,1,k+1)/(prim(iscaG2,j,1,k+1)-prim(iscaG2,j+1,1,k+1))) * &
                                    (prim(iscaG2,j,1,k+1)/(prim(iscaG2,j,1,k+1)-prim(iscaG2,j,1,k)))) * 0.5D0
                                
            elseif(prim(iscaG2,j,1,k) < 0.0D0 .and. prim(iscaG2,j,1,k+1) < 0.0D0 .and. &
                    prim(iscaG2,j+1,1,k) >= 0.0D0 .and. prim(iscaG2,j+1,1,k+1) < 0.0D0) then
                loc_deton_ratio = ((prim(iscaG2,j+1,1,k)/(prim(iscaG2,j+1,1,k)-prim(iscaG2,j+1,1,k+1))) * &
                                    (prim(iscaG2,j+1,1,k)/(prim(iscaG2,j+1,1,k)-prim(iscaG2,j,1,k)))) * 0.5D0
            elseif(prim(iscaG2,j,1,k) >= 0.0D0 .and. prim(iscaG2,j,1,k+1) > 0.0D0 .and. &
                    prim(iscaG2,j+1,1,k) < 0.0D0 .and. prim(iscaG2,j+1,1,k+1) > 0.0D0) then
                loc_deton_ratio = 1.0D0 - ((prim(iscaG2,j+1,1,k)/(prim(iscaG2,j+1,1,k)-prim(iscaG2,j+1,1,k+1))) * &
                                    (prim(iscaG2,j+1,1,k)/(prim(iscaG2,j+1,1,k)-prim(iscaG2,j,1,k)))) * 0.5D0

            elseif(prim(iscaG2,j,1,k) < 0.0D0 .and. prim(iscaG2,j,1,k+1) < 0.0D0 .and. &
                    prim(iscaG2,j+1,1,k) < 0.0D0 .and. prim(iscaG2,j+1,1,k+1) >= 0.0D0) then
                loc_deton_ratio = ((prim(iscaG2,j+1,1,k+1)/(prim(iscaG2,j+1,1,k+1)-prim(iscaG2,j+1,1,k))) * &
                                    (prim(iscaG2,j+1,1,k+1)/(prim(iscaG2,j+1,1,k+1)-prim(iscaG2,j,1,k+1)))) * 0.5D0
            elseif(prim(iscaG2,j,1,k) >= 0.0D0 .and. prim(iscaG2,j,1,k+1) >= 0.0D0 .and. &
                    prim(iscaG2,j+1,1,k) >= 0.0D0 .and. prim(iscaG2,j+1,1,k+1) < 0.0D0) then
                loc_deton_ratio = 1.0D0 - ((prim(iscaG2,j+1,1,k+1)/(prim(iscaG2,j+1,1,k+1)-prim(iscaG2,j+1,1,k))) * &
                                    (prim(iscaG2,j+1,1,k+1)/(prim(iscaG2,j+1,1,k+1)-prim(iscaG2,j,1,k+1)))) * 0.5D0
            endif
                                            
        elseif(detoncorn_flag(j,1,k) == 3) then

        ! The parallelogram case

        if(prim(iscaG2,j,1,k) >= 0.0D0 .and. prim(iscaG2,j,1,k+1) < 0.0D0) then
                loc_deton_ratio = ((prim(iscaG2,j,1,k)/(prim(iscaG2,j,1,k)-prim(iscaG2,j,1,k+1))) + &
                                    (prim(iscaG2,j+1,1,k)/(prim(iscaG2,j+1,1,k)-prim(iscaG2,j+1,1,k+1)))) * 0.5D0
            elseif(prim(iscaG2,j,1,k) < 0.0D0 .and. prim(iscaG2,j,1,k+1) >= 0.0D0) then
                loc_deton_ratio = ((prim(iscaG2,j,1,k+1)/(prim(iscaG2,j,1,k+1)-prim(iscaG2,j,1,k))) + &
                                (prim(iscaG2,j+1,1,k+1)/(prim(iscaG2,j+1,1,k+1)-prim(iscaG2,j+1,1,k)))) * 0.5D0
        
            elseif(prim(iscaG2,j,1,k) >= 0.0D0 .and. prim(iscaG2,j+1,1,k) < 0.0D0) then
                loc_deton_ratio = ((prim(iscaG2,j,1,k)/(prim(iscaG2,j,1,k)-prim(iscaG2,j+1,1,k))) + &
                                    (prim(iscaG2,j,1,k+1)/(prim(iscaG2,j,1,k+1)-prim(iscaG2,j+1,1,k+1)))) * 0.5D0
            elseif(prim(iscaG2,j,1,k) < 0.0D0 .and. prim(iscaG2,j+1,1,k) >= 0.0D0) then  
                loc_deton_ratio = ((prim(iscaG2,j+1,1,k)/(prim(iscaG2,j+1,1,k)-prim(iscaG2,j,1,k))) + &
                                    (prim(iscaG2,j+1,1,k+1)/(prim(iscaG2,j+1,1,k+1)-prim(iscaG2,j,1,k+1)))) * 0.5D0

            endif
            
        elseif(detoncorn_flag(j,1,k) == 4) then

        ! The ambigous case, average of all possible configuration

        if(prim(iscaG2,j,1,k) > 0.0D0) then
                loc_deton_ratio = 0.5D0 * &
                                    ((((prim(iscaG2,j,1,k)/(prim(iscaG2,j,1,k) - prim(iscaG2,j,1,k+1))) * &
                                    (prim(iscaG2,j,1,k)/(prim(iscaG2,j,1,k) - prim(iscaG2,j+1,1,k)))) * 0.5D0 + &  
                                    ((prim(iscaG2,j+1,1,k+1)/(prim(iscaG2,j+1,1,k+1)-prim(iscaG2,j+1,1,k))) * &
                                    (prim(iscaG2,j+1,1,k+1)/(prim(iscaG2,j+1,1,k+1)-prim(iscaG2,j,1,k+1)))) * 0.5D0) + &
                                    (1.0D0 - ((prim(iscaG2,j,1,k+1)/(prim(iscaG2,j,1,k+1) - prim(iscaG2,j,1,k))) * &
                                    (prim(iscaG2,j,1,k+1)/(prim(iscaG2,j,1,k+1) - prim(iscaG2,j+1,1,k+1)))) * 0.5D0 - &
                                    ((prim(iscaG2,j+1,1,k)/(prim(iscaG2,j+1,1,k)-prim(iscaG2,j,1,k))) * &
                                    (prim(iscaG2,j+1,1,k)/(prim(iscaG2,j+1,1,k)-prim(iscaG2,j+1,1,k+1)))) * 0.5D0))
            else
                loc_deton_ratio = 0.5D0 * &
                                    ((1.0D0 - ((prim(iscaG2,j,1,k)/(prim(iscaG2,j,1,k) - prim(iscaG2,j,1,k+1))) * &
                                    (prim(iscaG2,j,1,k)/(prim(iscaG2,j,1,k) - prim(iscaG2,j+1,1,k)))) * 0.5D0 - &
                                    ((prim(iscaG2,j+1,1,k+1)/(prim(iscaG2,j+1,1,k+1)-prim(iscaG2,j+1,1,k))) * &
                                    (prim(iscaG2,j+1,1,k+1)/(prim(iscaG2,j+1,1,k+1)-prim(iscaG2,j,1,k+1)))) * 0.5D0) + &
                                    (((prim(iscaG2,j,1,k+1)/(prim(iscaG2,j,1,k+1) - prim(iscaG2,j,1,k))) * &
                                    (prim(iscaG2,j,1,k+1)/(prim(iscaG2,j,1,k+1) - prim(iscaG2,j+1,1,k+1)))) * 0.5D0 + &
                                    ((prim(iscaG2,j+1,1,k)/(prim(iscaG2,j+1,1,k)-prim(iscaG2,j,1,k))) * &
                                    (prim(iscaG2,j+1,1,k)/(prim(iscaG2,j+1,1,k)-prim(iscaG2,j+1,1,k+1)))) * 0.5D0))
            endif
                                
        endif

        else   

    ! Completely empty mean zero
        loc_deton_ratio = 0.0D0

        endif

! Notice that the ratio is defined with grid corner
    ! as center, so the contribution of the fraction is
    ! applied to all the overlapped grids		
loc_deton_ratio = 0.25D0 * loc_deton_ratio
    deton_ratio(j,1,k) = deton_ratio(j,1,k) + loc_deton_ratio
    deton_ratio(j+1,1,k) = deton_ratio(j+1,1,k) + loc_deton_ratio
    deton_ratio(j,1,k+1) = deton_ratio(j,1,k+1) + loc_deton_ratio
    deton_ratio(j+1,1,k+1) = deton_ratio(j+1,1,k+1) + loc_deton_ratio

    enddo                          
enddo

! Copy the result to the current fraction
deton_loc_ratio(:,:,:) = deton_ratio(:,:,:)

    
! Now classify the grid according to the summed fraction    
do k = 1, nz, 1
    do j = 1, nx, 1

    ! Find out the effective ratio
        deton_ratio(j,1,k) = MIN(MAX(deton_ratio(j,1,k), deton_ratio_old(j,1,k)), 1.0D0 - flame_ratio(j,1,k))
    !deton_ratio(j,1,k) = MIN(deton_ratio(j,1,k), 1.0D0 - flame_ratio(j,1,k))

    if(deton_ratio(j,1,k) == 0.0D0) then

    ! Completely empty
    detongrid_flag(j,1,k) = 0

    elseif(deton_ratio(j,1,k) == 1.0D0) then

    ! Completely filled
    detongrid_flag(j,1,k) = 1

    else

    ! Partially filled
    detongrid_flag(j,1,k) = 2

    endif

    enddo
enddo

! Copy the results to ghost cells
CALL BOUNDARY1D_NM(deton_ratio,even, even, even, even, even, even)
!CALL BOUNDARY1D_INT(detongrid_flag, even)
            
100 FORMAT (6E13.6)
                                
end subroutine compute_detonratio

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine finds the initial input
! of energy due to deflagration/detonation
!
! Written by Leung Shing Chi in 2016
! Updated by Leung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE flame_ini()
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

! dummy variables
INTEGER :: j, k

! Local variables
REAL*8 :: epsilon_temp, rho_mid, temp_mid, rho_ash, eps_ash

! Local dummy
REAL*8 :: dummy
REAL*8 :: flame_mid

! Ash composition
REAL*8, dimension(totalion) :: x_burn

! Initilzation
WRITE(*,*) 'Set initial flame_qdot'
flame_qdot(:,:,:) = 0.0D0

! Find the energy input
!WRITE(*,*) 'Give the initial def/det energy'
!WRITE(*,*) length_step_z_min_part, length_step_r_part, length_step_z_part
!WRITE(*,"(10ES15.7)") rho2(1,1), temp2(1,1), flame_ratiO(1,1), epsilon2(1,1)

DO k = 1, nz, 1
    DO j = 1, nx, 1

    ! Store to local variables
    rho_mid = prim(irho,j,1,k) 
    temp_mid = temp2(j,1,k) 

    ! Energy injection and change in chemical composition are done
    ! to grids which are first assumed to be burnt
    IF(rho_mid > rho2_flame_min .and. flame_ratio(j,1,k) > 0.0D0) THEN

        CALL readtable_flameenergy(temp_mid, rho_mid, flame_mid, x_burn(:))
        CALL readtable_flamestate(temp_mid, rho_mid, rho_ash, dummy, eps_ash)

        !epsilon2(j,1,k) = epsilon2(j,1,k) + flame_mid * flame_ratio(j,1,k)
        epsilon(j,1,k) = (1.0D0 - flame_ratio(j,1,k)) * epsilon(j,1,k) + flame_ratio(j,1,k) * eps_ash
        prim(irho,j,1,k) = (1.0D0 - flame_ratio(j,1,k)) * rho_mid + flame_ratio(j,1,k) * rho_ash
        prim(ihe4:ini56,j,1,k) = (1.0D0 - flame_ratio(j,1,k)) * prim(ihe4:ini56,j,1,k) + flame_ratio(j,1,k) * x_burn(:) 
        burn_mass = burn_mass + rho_mid * 2.0D0 * pi * x(j) * dx(j) * dz(k)
        
    ENDIF

    ENDDO
ENDDO


! Copy to new results to ghost grids
CALL BOUNDARY()

END SUBROUTINE flame_ini

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine collects all the subroutines
! for updating the level-set
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_flame_radius()
use definition
use custom_def
implicit none

! Check timing with or without openmp
INTEGER :: time_start, time_end

INTEGER :: cr, cm
REAL :: rate            

!CALL system_clock(count_rate=cr)
!CALL system_clock(count_max=cm)
!rate = REAL(cr)
!WRITE(*,*) "system_clock rate ",rate
!CALL system_clock(time_start)

! Update the first level-set
! The deflagration is neglected as long as detonation has
! takes place long enough to surround the whole deflagration
IF(found_deton_flag == 0 .or. (found_deton_flag == 1 .AND. ABS(global_time - found_deton_time) <= 25000.0D0)) THEN
    CALL identify_flamegrid()
    CALL compute_flameratio()
ENDIF                

! Update the second level-set
IF(found_deton_flag == 1) THEN
    CALL identify_detongrid()
    CALL compute_detonratio()
ENDIF

! Backup the data
flame_ratio_old = flame_ratio
deton_ratio_old = deton_ratio

! Let the first level-set self-Propagate
IF(found_deton_flag == 0 .OR. (found_deton_flag == 1 .AND. ABS(global_time - found_deton_time) <= 25000.0D0)) THEN
    ! CALL update_scaG()
ENDIF

! Let the second level-set self-propagate
IF(found_deton_flag == 1) THEN
    ! call update_scaG2()
ENDIF

! Do the reinitilization and 
! comput the geometry for the 1st level-set
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

!CALL system_clock(time_end)   
!WRITE(*,*) 'Level set = ', REAL(time_end - time_start) / rate

end subroutine update_flame_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine do the propagation of level-set due 
! to its normal-flow (propagation of flame in the 
! supernovae context)
! Written by Leung SHing CHi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_scaG()
use definition
use custom_def
implicit none

! dummy variables
integer :: j, k  

! local geometry information
real*8 :: norm_r, norm_z

! local sound speed
real*8 :: cs_grid

! Local propagation
real*8, allocatable, dimension(:,:,:) :: scaG_flux1
allocate(scaG_flux1(-2:nx+3, -2:ny+3, -2:nz+3))

! initilization
cs_grid = 0.0D0
scaG_flux1 = 0.0D0

do k = 1, nz, 1
    do j = 1, nx, 1

    if(flamegrid_flag(j,1,k) > 1) then
        call update_flame_velocity(j, k, cs_grid)
        scaG_flux1(j,1,k) = cs_grid
    endif

    enddo
enddo  

! Copy the results to ghost cells
call boundary1D_NM(scaG_flux1, even, even, even, even, even, even)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now update the level set
do j = 1, nx, 1
    do k = 1, nz, 1
        prim(iscaG1,j,1,k) = prim(iscaG1,j,1,k) + dt * scaG_flux1(j,1,k)
    enddo
enddo

! Fill the level set ghost grid
CALL BOUNDARY1D_NM(prim(iscaG1,:,:,:), even, even, even, even, even, even)

deallocate(scaG_flux1)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine do the propagation of level-set due
! to its normal-flow (propagation of flame in the
! supernovae context)
! Written by Leung SHing CHi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_scaG2()
use definition
use custom_def
implicit none
    
! dummy variables
integer :: j, k

! Geometry information
real*8 :: norm_r, norm_z

! local sound speed
real*8 :: cs_grid

! local propagation
real*8, dimension(-2:nx+3,-2:nz+3) :: scaG_flux1
    
! initialization
cs_grid = 0.0D0
scaG_flux1 = 0.0D0

do k = 1, nz, 1
    do j = 1, nx, 1

    ! Only the front propagate
        if(detongrid_flag(j,1,k) > 1) then

        if(prim(irho,j,1,k) > 3.2D-11) then        
            
            ! Pathological substained detonation
            call readtable_detonvel(prim(irho,j,1,k), cs_grid)

        elseif(prim(irho,j,1,k) > 1.6D-12 .and. prim(irho,j,1,k) < 3.2D-11) then  ! between 1x10^6 and 2x10^7

            ! Chapman-Jouguet detonation
            call findsoundspeed(j, k, cs_grid)

        else

            ! No detonation otherwise
            cs_grid = 0.0D0

        endif

        scaG_flux1(j,k) = cs_grid       

        endif

    enddo
enddo
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, I calculate the flux by the flame
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do k = 1, nz, 1
    do j = 1, nx, 1
        prim(iscaG2,j,1,k) = prim(iscaG2,j,1,k) + dt * scaG_flux1(j,k)
    enddo
enddo
        
CALL BOUNDARY1D_NM(prim(iscaG2,:,:,:), even, even, even, even, even, even)

end subroutine update_scaG2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine finds the flame propagation velocity by finding 
! the maximum of the three terms
! 1. laminar flame (Timmes 1992)
! 2. Rayleigh-taylor instability (Schmidty2006a)
! 3. turbulent flame (Niemeyer1995)
!
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_flame_velocity(j_in, k_in, flame_localvel)
use definition
use custom_def
implicit none

! Critical number for DDT transition
real*8 :: Ka_crit = 0.5D0

! The Gibson length scale
real*8 :: Gibson_length

! The width of flame
real*8 :: flame_length

! The various flame speed
real*8 :: flame_vel_lam, flame_vel_helm, flame_vel_turb

! Output effective flame speed
real*8 :: flame_localvel

! Check if really the level set is changed
integer :: modify_scaG_flag

! Input grid number
integer :: j_in, k_in

! Dummy variables
integer :: j, k

! Local dummy variables
real*8 :: rho_fuel, dist

! First store in local variables
rho_fuel = prim(irho,j_in,1,k_in)

! Then find the laminar flame speed
flame_vel_lam = 3.06666667D-4 * (rho_fuel / 3.24D-9) ** 0.805D0
!flame_vel_lam = 4.29D-4 * (rho_fuel / 3.24D-9) ** 0.805D0

! When there is turbulence, then map to the effective turblent flame speed
if(turb_flag == 1) then

    ! The traditional one
    flame_vel_turb = flame_vel_lam * DSQRT(1.0D0 + 8.0D0 / 3.0D0 * prim(iturbq,j_in,1,k_in) / flame_vel_lam ** 2)  

    ! The updated one (but slower)
    !flame_vel_turb = flame_vel_lam * DSQRT(1.0D0 + 1.228D0 * prim(iturbq,j_in,1,k_in) / flame_vel_lam ** 2)  !1.228, 200.0

    ! The best fitted one
    !flame_vel_turb = flame_vel_lam * (1.0D0 + 0.654D0 * (prim(iturbq,j_in,1,k_in) / flame_vel_lam**2) ** 0.5985)

    ! The pre-Schmidt one
    !flame_vel_turb = MAX(flame_vel_lam, DSQRT(2.0D0 * prim(iturbq,j_in,1,k_in)))

endif

! Only flame above density threshold can propagate
if(rho_fuel > rho2_flame_min .and. rho_fuel < rho2_flame_max) then

    flame_localvel = flame_vel_turb
    !flame_localvel = flame_vel_lam

else

    flame_localvel = 0.0D0

endif

! If deton_flag is on, then check together whether the
! flame is thick enough to cause DDT
if(deton_flag == 1) then

    ! Typing DDT required density
    if(rho_fuel <= 3.24D-11 .and. rho_fuel > rho2_deton_min) then

        ! Characteristic length scale
        flame_length = 5.1381D-9 * (rho_fuel/1.62D-9)**(-1.9375D0)
        gibson_length = dx(j) * (flame_vel_lam**2 / 2.0D0 / prim(iturbq,j_in,1,k_in))**1.5D0

        ! Flame width > turbulence eddy-overturn scale is needed
        ! But the exact value is not yet discovered
        if(flame_length > Ka_crit * gibson_length) then
            !if(flame_vel_turb >= 0.2268D0 * (rho_fuel/1.62D-9) ** 1.4508) then
            !if(prim(iturbq,j_in,1,k_in) >= 1.0D0**(2.0D0/3.0D0) * 0.01982D0 * (rho_fuel/1.62D-9) ** 8.705D0/3) then
        
        ! The grid must be only partially burnt
            if(rho_fuel <= rho2_deton_max .and. flamegrid_flag(j_in,1,k_in) > 1) then

            ! For the first flame, then switch to 
            ! the output timescale and initialize
            ! the related detonation level set
                IF(found_deton_flag == 0) then

                    write(*,*) 'Occured at = ', global_time
                    write(*,*) 'Found deton at ', j_in, k_in
                    write(*,*) 'Flame length = ', flame_length
                    write(*,*) 'Gibson length = ', gibson_length
                    write(*,*)

                    !cfl = 0.2D0  
                    ! output_profiletime = 0.5D4
                    ! output_flametime = 0.5D4
                    ! output_turbtime = 0.5D4
                    ! output_Helmtime = 0.5D4
                    found_deton_flag = 1
                    found_deton_time = global_time

                    deton_ratio_old(:,:,:) = 0.0D0

                    ! Plant the level set
                    do j = 1, nx, 1     
                        do k = 1, nz, 1
                        dist = DSQRT((x(j) - x(j_in))**2 + (z(k) - z(k_in))**2) - 10.0D0
                        prim(iscaG2,j,1,k) = -dist
                        enddo
                    enddo
        
                    CALL BOUNDARY1D_NM(prim(iscaG2,:,:,:), even, even, even, even, even, even)
                
                    ! Make sure the distance property of the level set
                    ! needs to be preserved always
                    call reinitialization3()               !You need this for aspherical initial flame
                    !call update_scaG2corn()
                    call identify_detongrid()
                    call compute_detonratio()

                    !deton_ratio_old(:,:) = deton_ratio(:,:)

                ELSEIF(found_deton_flag == 1 .and. prim(iscaG2,j_in,1,k_in) < -100.0D0) THEN

                    ! If this is the 2nd time (or above) for the
                    ! detonation to start, then simply plant
                    ! the detonation seed 
                    modify_scaG_flag = 0

                    write(*,*) 'Occured at = ', global_time
                    write(*,*) 'Found deton at ', j_in, k_in
                    write(*,*) 'Flame length = ', flame_length
                    write(*,*) 'Gibson length = ', gibson_length
                    write(*,*)

                        do j = j_in - 3, j_in + 3, 1
                            do k = k_in - 3, k_in + 3, 1

                            dist = DSQRT((x(j) - x(j_in))**2 + (z(k) - z(k_in))**2) - 10.0D0 
                            prim(iscaG2,j,1,k) = -dist
                            modify_scaG_flag = 1
                            
                            enddo
                        enddo   

                    if(modify_scaG_flag == 1) then

                        CALL BOUNDARY1D_NM(prim(iscaG2,:,:,:), even, even, even, even, even, even)
                
                        ! Preserve the distance property
                        call reinitialization3()               !You need this for aspherical initial flame
                        !call update_scaG2corn()
                        call identify_detongrid()
                        call compute_detonratio()

                    endif

                endif
        
            endif
    
        endif

    endif

endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine calculates the energy release 
! By assuming 12C --> 24Mg
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE burn_phase1b
USE definition
use custom_def
IMPLICIT NONE

! Dummy variables
INTEGER :: j, k

! Flag for finding temp
integer :: flag_notfindtemp

! Local variables
REAL*8 :: rho_mid, temp_mid, vol_mid

! Change of fraction for both level sets
REAL*8 :: x1, x2

! Local Atwood Number
REAL*8 :: at_no

! Ne-22 binding energy                       
!real*8, parameter :: qbar_He4 = 7.57791D-3
!real*8, parameter :: qbar_C12 = 8.23266D-3
!real*8, parameter :: qbar_O16 = 8.55261D-3
!real*8, parameter :: qbar_Ne20 = 8.61316D-3
!real*8, parameter :: qbar_Mg24 = 8.86028D-3
!real*8, parameter :: qbar_Si28 = 9.06265D-3
!real*8, parameter :: qbar_Ne22 = 8.66501D-3
!real*8, parameter :: qbar_Ni56 =  9.27328D-3

! local flame energy, defined by Dqbar_C12 * xiso_ini(c12) + Dqbar_Ne22 * metallicity
! i.e. 6.2762E-4 * (xiso_ini(c12) - Z/2) + 1.9527E-4 * Z
! Example: Z = 0.00, flame_ene = 3.13810E-4
! Example: Z = 0.02, flame_ene = 3.11439E-4
! Example: Z = 0.04, flame_ene = 3.09068E-4
! Example: Z = 0.06, flame_ene = 3.06698E-4
! Example: Z = 0.08, flame_ene = 3.04327E-4
! Example: Z = 0.10, flame_ene = 3.01956E-4
real*8, parameter :: flame_ene = 3.11439D-4

! Ash composition
real*8, dimension(totalion) :: x_fuel1 = (/0.0D0, 0.51D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
real*8, dimension(totalion) :: x_ash1 = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.51D0, 0.0D0, 0.0D0/)

! local deton energy, defined by Dqbar_C12 * xiso_ini(c12) + Dqbar_Ne22 * metallicity
! i.e. 6.2762E-4 * (xiso_ini(c12) - Z/2) + 1.9527E-4 * Z
real*8, parameter :: deton_ene = 3.11439D-4  ! He to Ne


! Ash composition
! real*8, dimension(totalion) :: x_fuel2 = (/1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
! real*8, dimension(totalion) :: x_ash2 = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0/)

! Other variables
real*8 :: esum, nqse_burntime2, ratio

! if(debug_flag == 1) write(*,*) 'In Flame'

! First, initialize the grid, to erase all past data
flame_qdot(:,:,:) = 0.0D0
deton_qdot(:,:,:) = 0.0D0

!call findenergy
!write(*,*) 'Before: ', energy2
!esum = 0.0D0

! Then, find the energy release due to nuclear fusion

do k = 1, nz, 1
    do j = 1, nx, 1

        rho_mid = prim(irho,j,1,k)
        temp_mid = temp2(j,1,k)
        vol_mid = vol(j,1,k)

        ! Just compute the change of fraction occupies
        ! by deflagration or detonation
        !x1 = flame_ratio(j,1,k) - flame_ratio_old(j,1,k)
        x1 = MIN(flame_ratio(j,1,k) - flame_ratio_old(j,1,k), prim(ic12,j,1,k) / 0.51D0)
        x2 = MIN(deton_ratio(j,1,k) - deton_ratio_old(j,1,k), prim(ic12,j,1,k) / 0.51D0)
        !x2 = MIN(deton_ratio(j,1,k) - deton_ratio_old(j,k), xiso(j,k,che4))

        ! Remember to switch flame_ini as well
        ! When there is a change in area fraction
        ! Change the chemical composition accordingly
        ! And inject temperature
        if(x1 > 0.0D0 .and. rho_mid > rho2_flame_min .and. rho_mid < rho2_flame_max) then
        !if(x1 > 0.0D0 .and. rho_mid > rho2_deton_min) then
        flame_qdot(j,1,k) = flame_ene * x1
        epsilon(j,1,k) = epsilon(j,1,k) + flame_qdot(j,1,k) 
        prim(ihe4:ini56,j,1,k) = prim(ihe4:ini56,j,1,k) - x1 * x_fuel1(:) + x1 * x_ash1(:)     
        burn_mass = burn_mass + x1 * rho_mid * vol_mid
        endif
    
        !Repeat the same procedure for detonation
        if(x2 > 0.0D0 .and. rho_mid > rho2_deton_min .and. rho_mid < rho2_deton_max) then
        !if(x2 > 0.0D0 .and. rho_mid > 1.62D-13) then

        !nqse_burntime2 = 3.4838D0 * (rho_mid / 1.62D-10)**(-2)
        !ratio = MIN(1.0D0, dt / nqse_burntime2)

        deton_qdot(j,1,k) = deton_ene * x2 
        epsilon(j,1,k) = epsilon(j,1,k) + deton_qdot(j,1,k) 
        prim(ihe4:ini56,j,1,k) = prim(ihe4:ini56,j,1,k) - x2 *  x_fuel1(:) + x2 * x_ash1(:)
        burn_mass = burn_mass + x2 * rho_mid * vol_mid

        endif

    enddo
enddo

! Copy the results to ghost cells
!call boundary1D(epsilon2, even)

100 format(10ES15.6)   
101 format(A6, 2I5, 10ES15.6)

end subroutine burn_phase1b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine calculates the energy release
! By assuming 16O + 24Mg --> 28Si
!
! Written by Leung Shing Chi in 2016
! Updated by LEung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine burn_phase2b
use definition
use custom_def
implicit none

! Dummy variables
integer :: j, k

! Flag for finding temp
integer :: flag_notfindtemp               

! Local variables
real*8 :: rho_mid, temp_mid       

! Change of fraction for both level sets
real*8 :: x1, x2

! Local Atwood Number
real*8 :: at_no

! Ne-22 binding energy
!real*8, parameter :: qbar_C12 = 8.23266D-3
!real*8, parameter :: qbar_O16 = 8.55261D-3
!real*8, parameter :: qbar_Ne20 = 8.61316D-3
!real*8, parameter :: qbar_Mg24 = 8.86028D-3
!real*8, parameter :: qbar_Si28 = 9.06265D-3
!real*8, parameter :: qbar_Ne22 = 8.66501D-3
!real*8, parameter :: qbar_Ni56 = 9.27328D-3

! local flame energy, defined by 
! Energy release form 16O  --> 28Si: Dqbar_O16 * (xiso_ini(o16) - Z/2))
! Energy release from 24Mg --> 28Si: Dqbar_Mg24 * (xiso_ini(c12) + Z/2)
! i.e. 5.1004E-4 * (xiso_ini(o16) - Z/2)
! i.e. 2.2037E-4 * (xiso_ini(c12) + Z/2)
! Example: 
real*8 :: burn_ene1 = 5.1004D-4	! Energy release form 16O  --> 28Si
real*8 :: burn_ene2 = 2.0237D-4   ! Energy release from 24Mg --> 28Si  
real*8 :: burn_ene3 = 1.69537D-3   ! Energy release from 4He  --> 56Ni
!real*8 :: burn_ene3 = 6.6012D-4   ! Energy release from 20Ne  --> 56Ni
!real*8 :: burn_ene4 = 2.0163D-4   ! Energy release from 28Si  --> 56Ni

! local deton energy            
real*8 :: deton_ene 

! Ash composition
real*8, dimension(totalion) :: x_mid
real*8, dimension(totalion) :: x_fuel1 = (/0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
real*8, dimension(totalion) :: x_fuel2 = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0/)
real*8, dimension(totalion) :: x_ash = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0/)

real*8, dimension(totalion) :: x_fuelb = (/1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
real*8, dimension(totalion) :: x_ashb = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0/)

! NQSE timescale
real*8 :: nqse_burntime, nqse_burntime2, temp_nse

! Other variables
real*8 :: esum

! if(debug_flag == 1) write(*,*) 'In Flame'

! First, initialize the grid, to erase all past data
burn_qdot(:,:,:) = 0.0D0

!call findenergy
!write(*,*) 'Before: ', energy2
!esum = 0.0D0          

! Then, find the energy release due to nuclear fusion
do k = 1, nz, 1
    do j = 1, nx, 1

        ! Pass the data to local variables
        rho_mid = prim(irho,j,1,k)
        temp_mid = temp2(j,1,k)
        x_mid(:) = prim(ihe4:ini56,j,1,k)

    !if(rho_mid < 2.43D-11) cycle  !~3.6 standardized by the rho3, Z2, Ka1 case
        !if(rho_mid < 1.62D-11) cycle  !~3.6 standardized by the rho3, Z2, Ka3.2 (SQRT(10)) case
    !if(ABS(p2(j+1,k) - p2(j-1,k)) / p2(j,k) >= 0.5D0 .or. &
    !   ABS(p2(j,k+1) - p2(j,k-1)) / p2(j,k) >= 0.5D0) cycle


    if(nse_flag(j,1,k) == 2) cycle

    ! Only the flame_ratio or deton_ratio = 1 can burn
    ! Meaning that completely C-burnt cells can start
    ! Notice the MAX nature of flame_ratio or deton_ratio
    !if(flame_loc_ratio(j,1,k) == 1.0D0 .and. rho_mid > rho2_flame_min) then
        if(burn_ratio(j,1,k) == 1.0D0 .and. rho_mid > 2.43D-11) then

            ! When there is a change in area fraction
            ! Change the chemical composition accordingly
            ! And inject temperature

            ! Calculate the energy timescale
            !nqse_burntime = EXP(182.06D0 / temp_mid - 46.054D0) / 4.9282D-6
            temp_nse = -8.6803D0 + 1.8787D0 * LOG10(rho_mid * 6.171D17)
            nqse_burntime = EXP(182.06D0 / temp_nse - 46.054D0) / 4.9282D-6

            ! Just compute the change of fraction by considering
            ! the maximum energy release equals to the amount of fuel remained
            x1 = x_mid(co16)
            x2 = x_mid(cmg24)

            if(dt > nqse_burntime) then
        
                ! Complete burning occurs
                burn_qdot(j,1,k) = burn_ene1 * x1 + burn_ene2 * x2
                epsilon(j,1,k) = epsilon(j,1,k) + burn_qdot(j,1,k)
                prim(ihe4:ini56,j,1,k) = x_mid(:) - x1 * x_fuel1(:) - x2 * x_fuel2(:) + (x1 + x2) * x_ash(:)

                ! If it is completely burnt, then go to the next burning phase
                nse_flag(j,1,k) = 2

                else

                ! InComplete burning occurs
                burn_qdot(j,1,k) = (burn_ene1 * x1 + burn_ene2 * x2) * (dt / nqse_burntime)
                epsilon(j,1,k) = epsilon(j,1,k) + burn_qdot(j,1,k)
                prim(ihe4:ini56,j,1,k) = x_mid(:) + (-x1 * x_fuel1(:) - x2 * x_fuel2(:) + &
                    (x1 + x2) * x_ash(:)) * (dt / nqse_burntime)

                nse_flag(j,1,k) = 1

            endif

        endif
    enddo
enddo                  

end subroutine burn_phase2b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine compute the binding energy for 
! a given composition. The subroutine takes 
! chemical composition as input and return the
! binding energy.
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_binde(x_in, binde_out)
use definition
use custom_def
implicit none

! Input chemical composition
real*8 :: x_in(1:totalion)

! Output binding energy 
real*8 :: binde_out

! Dummy variable
integer :: i

! Binding energy
real*8 :: binde_sum

! Initialiation
binde_sum = 0.0D0

! Sum all the isotopes
do i = 1, totalion, 1
    binde_sum = binde_sum + x_in(i) / mion(i) * bion(i)
enddo

! Change to code unit
binde_out = 1.78D-27 * binde_sum

end subroutine compute_binde

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This subroutine calculate the sound speed 
! assuming you are using the HELMHOLTZ EOS
! Switch it accordingly if you have other EOS
! Written by Leung Shing Chi
! For more information about HELMHOLTZ, refer to Timmes website 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine findsoundspeed(j_in, k_in, cs_out)

use definition
use custom_def
implicit none

! Input grid number
integer :: j_in, k_in

! Output sound speed
real*8 :: cs_out

! Pass the global hydro variables to the local variables
real*8 :: temp_in, rho_in, abar_in, zbar_in

! Do the pass
rho_in = prim(irho,j_in,1,k_in)
temp_in = temp2(j_in,1,k_in)
abar_in = abar2(j_in,1,k_in)
zbar_in = zbar2(j_in,1,k_in)

! Pass to the HELMHOLTZ
call HELM_EOSSOUNDSPEED(rho_in, temp_in, abar_in, zbar_in, cs_out)
cs_out = cs_out ! * xiso(j_in,k_in,cc12) / xc12_ini

end subroutine findsoundspeed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!