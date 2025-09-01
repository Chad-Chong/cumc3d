!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine interpolates the neutrino table to 
! Get the effective mass and the emissivities
!
! Written by Leung Shing Chi in 2016 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine findnuspec
USE DEFINITION
USE CUSTOM_DEF
implicit none

! dummy variables
integer :: i, j, k

! Target grid in the table
integer :: itemp, irho2

! Derivation of the actual point from the grid boundary
real*8 :: dtemp, drho

! transverse mass
real*8 :: mass_t

! Neutrino emissivities
real*8 :: emiss

! Neutrino frequence
real*8 :: nu_freq

! Initialization
nu_phi = 0.0D0

do k = 1, nz
    do j = 1, nx

        if(prim(irho,j,1,k) > 1.62D-11 .and. temp2(j,1,k) > 1.0D0) then

    ! Do the standard bilinear interpolation 
        irho2 = INT((LOG10(prim(irho,j,1,k) * 6.171D17) - 7.0D0) / 0.1D0) + 1
        itemp = INT(temp2(j,1,k)) 

        drho = ((LOG10(prim(irho,j,1,k) * 6.171D17) - 7.0D0) / 0.1D0 + 1.0D0) - DBLE(irho2)
        dtemp = temp2(j,1,k) - DBLE(itemp) 

        mass_t = nutable_mass(itemp, irho2) + &
                    dtemp * (1.0D0 - drho) * (nutable_mass(itemp + 1, irho2) - nutable_mass(itemp, irho2)) + &
                    drho * (1.0D0 - dtemp) * (nutable_mass(itemp, irho2 + 1) - nutable_mass(itemp, irho2)) + &
                    dtemp * drho * (nutable_mass(itemp + 1, irho2 + 1) - nutable_mass(itemp, irho2))

    emiss = nutable_emiss(itemp, irho2) + &
                dtemp * (1.0D0 - drho) * (nutable_emiss(itemp + 1, irho2) - nutable_emiss(itemp, irho2)) + &
                drho * (1.0D0 - dtemp) * (nutable_emiss(itemp, irho2 + 1) - nutable_emiss(itemp, irho2)) + &
                dtemp * drho * (nutable_emiss(itemp + 1, irho2 + 1) - nutable_emiss(itemp, irho2))

        do i = 1, 10, 1

        ! Fix the frequence
            !nu_freq = 10.0D0 ** (-1.0D0 + 0.1D0 * DBLE(i))
        nu_freq = DBLE(i)

        ! Use the analytic form given in the article
        ! 3.22474D15 = length unit^3 = (1.4774 km)^3
            nu_phi(i) = nu_phi(i) + (4.5215D25 * emiss * 0.1425776426D0 / (0.08614D0 * temp2(j,1,k)) * &
            (nu_freq / 0.08614D0 / temp2(j,1,k)) ** 3.180657028D0 * &
                EXP(-1.108192299D0 * nu_freq / 0.08614D0 / temp2(j,1,k)) + &
            2.115D30 * (0.08614D0 * temp2(j,1,k)) * (mass_t * 0.511D0 ** 2) ** 3 * &
            EXP(-nu_freq / 0.08614D0 / temp2(j,1,k))) * vol(j,1,k) * 3.22474D15  
    
        enddo

        endif

    enddo
enddo

end subroutine findnuspec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This subroutine reads in the nu_mass.dat and 
! nu_emiss.dat table. 
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_nutable
USE DEFINITION
USE CUSTOM_DEF
implicit none

! dummy variables
integer :: i, j

! open the files
open(unit=898, file='./src/lib/nu_mass.dat', action='read')
open(unit=899, file='./src/lib/nu_emiss.dat', action='read')

! read the table
do i = 1, temp_rowno3, 1
    read(898,*) (nutable_mass(i,j), j=1,den_rowno3,1)
    read(899,*) (nutable_emiss(i,j), j=1,den_rowno3,1)
enddo

close(898)
close(899)

end subroutine read_nutable