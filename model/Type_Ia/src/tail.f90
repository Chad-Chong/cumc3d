subroutine locate_tail_grid(dens, tail_count, tail_position)
use definition
USE CUSTOM_DEF
implicit none

!dummy variables
integer :: j, k

!input density
real*8, dimension(-2:nx+3, 1, -2:nz+3) :: dens

! Output intersection point number
integer :: tail_count

REAL*8 :: tail_rad

! Output array of intersection point
real*8, dimension(nx*nz, 2) :: tail_position

! Initialization
tail_count = 0 
tail_rad = 0.0D0

!write(*,*) 'In locate tail grid'

! For the first level set
do k = 1, nz, 1
    do j = 1, nx, 1

    if(dens(j,1,k) * dens(j+1,1,k) < 0.0D0) then
        tail_count = tail_count + 1
        tail_position(tail_count, 1) = x(j)
        tail_position(tail_count, 2) = z(k)
        if(tail_position(tail_count, 1) > tail_rad) tail_rad = tail_position(tail_count, 1)
    endif

    if(dens(j,1,k) * dens(j,1,k+1) < 0.0D0) then             
        tail_count = tail_count + 1
        tail_position(tail_count, 1) = x(j)
        tail_position(tail_count, 2) = z(k)
        if(tail_position(tail_count, 1) > tail_rad) tail_rad = tail_position(tail_count, 1) 
    endif

    enddo
enddo


end subroutine

subroutine locate_min_tail_distance(tail_count, tail_position, tail_distance)
use definition
USE CUSTOM_DEF
implicit none

! dummy variables
integer :: j, k, k2

! Input intersection point number
integer :: tail_count

! The local minimal distance
real (selected_real_kind(15,307)):: distance, last_distance

! Output array for the minimal distance
real (selected_real_kind(15,307)), dimension(-2:nx+3, -2:nz+3), intent(out):: tail_distance

! Input array of intersection point positions
real (selected_real_kind(15,307)), dimension(nx * nz, 2), intent(in) :: tail_position


do j = 1, nx, 1
    do k = 1, nz, 1

        last_distance = 10.0D0 * DBLE(MAX(nx,nz)) * MAX(dx(j), dz(k))

        ! Search for the minimal distance
        do k2 = 1, tail_count, 1

            distance = DSQRT((tail_position(k2,1) - x(j)) ** 2 + & 
                        (tail_position(k2,2) - z(k)) ** 2)

            if ((x(j)**2+z(k)**2) < (tail_position(k2,1)**2+tail_position(k2,2)**2)) then
                distance = -distance
            endif

            if((tail_count > 1 .and. ABS(distance) <= ABS(last_distance)) .or. tail_count == 1) then
                tail_distance(j,k) = distance
                last_distance = distance
            endif

            if(k2 == 1) last_distance = distance
        enddo

    enddo
enddo


!WRITE(*,*) tail_distance(1,1), tail_count

100 FORMAT (6E13.6)

end subroutine