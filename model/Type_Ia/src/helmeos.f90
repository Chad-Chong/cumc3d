!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This SUBROUTINE allocates the necessary variables 
! that will be used in this module.
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE buildHelm
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

ALLOCATE(abar2 (-2:nx+3,-2:ny+3,-2:nz+3)) 
ALLOCATE(zbar2(-2:nx+3,-2:ny+3,-2:nz+3))
! ALLOCATE(xiso(-2:nx+3,-2:ny+3,-2:nz+3, totalion))
ALLOCATE(temp2(-2:nx+3,-2:ny+3,-2:nz+3))  

!    if(AMR_level1_flag == 1) then
!       ALLOCATE(abar3 (-2:nx+3,-2:ny+3,-2:nz+3))
!       ALLOCATE(zbar3 (-2:nx+3,-2:ny+3,-2:nz+3))

!       ALLOCATE(xiso3(-2:nx+3,-2:ny+3,-2:nz+3, totalion))
!    endif



! ALLOCATE(xiso_old(-2:nx+3,-2:ny+3,-2:nz+3, totalion))
! ALLOCATE(xiso1(-2:nx+3,-2:ny+3,-2:nz+3, totalion))
! ALLOCATE(delta_xiso(-2:nx+3,-2:ny+3,-2:nz+3, totalion))

END SUBROUTINE buildHelm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE read_helm_table()
INCLUDE 'implno.dek'
INCLUDE 'helm_table_storage.dek'

! This routine reads the helmholtz eos file, and
! must be called once before the helmeos routine is invoked.

! declare local variables
integer          i,j
double precision tsav,dsav,dth,dt2,dti,dt2i,dt3i, &
                    dd,dd2,ddi,dd2i,dd3i


! open the file (use softlinks to input the desired table)

OPEN(unit=19,file='/home/cnchong/Codes/cumc3d/model/Type_Ia/src/lib/helm_table.dat',status='old')


! for standard table limits
tlo   = 3.0d0
thi   = 13.0d0
tstp  = (thi - tlo)/float(jmax-1)
tstpi = 1.0d0/tstp
dlo   = -12.0d0
dhi   = 15.0d0
dstp  = (dhi - dlo)/float(imax-1)
dstpi = 1.0d0/dstp

! read the helmholtz free energy and its derivatives
do j=1,jmax
    tsav = tlo + (j-1)*tstp
    t(j) = 10.0d0**(tsav)
    do i=1,imax
        dsav = dlo + (i-1)*dstp
        d(i) = 10.0d0**(dsav)
        read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                fddt(i,j),fdtt(i,j),fddtt(i,j)
    enddo
enddo


! read the pressure derivative with density table
do j=1,jmax
    do i=1,imax
        read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
    enddo
enddo

! read the electron chemical potential table
do j=1,jmax
    do i=1,imax
        read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
    enddo
enddo

! read the number density table
do j=1,jmax
    do i=1,imax
        read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
    enddo
enddo

! close the file and write a summary message
close(unit=19)


! construct the temperature and density deltas and their inverses
do j=1,jmax-1
    dth          = t(j+1) - t(j)
    dt2         = dth * dth
    dti         = 1.0d0/dth
    dt2i        = 1.0d0/dt2
    dt3i        = dt2i*dti
    dt_sav(j)   = dth
    dt2_sav(j)  = dt2
    dti_sav(j)  = dti
    dt2i_sav(j) = dt2i
    dt3i_sav(j) = dt3i
end do
do i=1,imax-1
    dd          = d(i+1) - d(i)
    dd2         = dd * dd
    ddi         = 1.0d0/dd
    dd2i        = 1.0d0/dd2
    dd3i        = dd2i*ddi
    dd_sav(i)   = dd
    dd2_sav(i)  = dd2
    ddi_sav(i)  = ddi
    dd2i_sav(i) = dd2i
    dd3i_sav(i) = dd3i
enddo



!      write(6,*)
!      write(6,*) 'finished reading eos table'
!      write(6,04) 'imax=',imax,' jmax=',jmax
!04    format(1x,4(a,i4))
!      write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
!      write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
!03    format(1x,4(a,1pe11.3))
!      write(6,*)

RETURN
END SUBROUTINE read_helm_table

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE initialize_network()
USE ecaptable_module
IMPLICIT NONE

CALL init_my_iso7
CALL read_nse_table
CALL read_nse_table2
CALL read_flame_table
CALL read_deton_table
CALL readEcapRate

END SUBROUTINE initialize_network

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE init_my_iso7
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

INTEGER :: i

burn_mass = 0.0D0

! set the id numbers of the elements
che4  = 1
cc12  = 2
co16  = 3
cne20 = 4
cmg24 = 5
csi28 = 6
cni56 = 7

! set the names of the elements
ionam(che4)  = 'he4'
ionam(cc12)  = 'c12'
ionam(co16)  = 'o16'
ionam(cne20) = 'ne20'
ionam(cmg24) = 'mg24'
ionam(csi28) = 'si28'
ionam(cni56) = 'ni56'

! set the number of nucleons in the element
aion(che4)  = 4.0d0
aion(cc12)  = 12.0d0
aion(co16)  = 16.0d0
aion(cne20) = 20.0d0
aion(cmg24) = 24.0d0
aion(csi28) = 28.0d0
aion(cni56) = 56.0d0

! set the number of protons in the element
zion(che4)  = 2.0d0
zion(cc12)  = 6.0d0
zion(co16)  = 8.0d0
zion(cne20) = 10.0d0
zion(cmg24) = 12.0d0
zion(csi28) = 14.0d0
zion(cni56) = 28.0d0

! set the binding energy of the element
bion(che4)  =  28.29603d0
bion(cc12)  =  92.16294d0
bion(co16)  = 127.62093d0
bion(cne20) = 160.64788d0
bion(cmg24) = 198.25790d0
bion(csi28) = 236.53790d0
bion(cni56) = 484.00300d0

! set the number of neutrons and mass
DO i = 1, totalion, 1
    nion(i) = aion(i) - zion(i)
ENDDO

! mass of each isotope
DO i = 1, totalion, 1
    mion(i) = nion(i)*1.67492721184d-24 + zion(i)*1.67262163783d-24 - bion(i)*mev2gr
ENDDO

! molar mass      
DO i = 1, totalion, 1
wion(i) = avo * mion(i)             
ENDDO

! a common approximation   
DO i = 1, totalion, 1
    wion(i) = aion(i)       
ENDDO

END SUBROUTINE init_my_iso7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE read_nse_table
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

INTEGER :: i, j, k
REAL*8 :: dummy

OPEN(unit=500, file='/home/cnchong/Codes/cumc3d/model/Type_Ia/src/lib/nse_table_7iso.dat',action='read')
DO i = 0, den_rowno_nse, 1
    DO j = 0, temp_rowno_nse, 1
        READ(500,*) dummy, dummy, nsetable_binde(i,j), (nsetable_xiso(i,j,k), k = 1, totalion)
        nsetable_binde(i,j) = nsetable_binde(i,j) / 9.0D20
    ENDDO
ENDDO

CLOSE(500)

END SUBROUTINE read_nse_table

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE read_nse_table2
USE DEFINITION
USE CUSTOM_DEF
IMPLICIT NONE

INTEGER :: i, j, k, k2
REAL*8 :: dummy

OPEN(unit=500, file='/home/cnchong/Codes/cumc3d/model/Type_Ia/src/lib/nse_table_495iso.dat',action='read')

DO i = 0, den_rowno_nse2, 1
    DO j = 0, ye_rowno_nse2, 1
    READ(500,*) (nsetable2_head(i,j,k2), k2 = 1, 3)
    DO k = 0, ent_rowno_nse2, 1
        READ(500,*) (nsetable2_binde(i,j,k,k2), k2 = 1, 6)
    nsetable2_binde(i,j,k,4) = nsetable2_binde(i,j,k,4) / 1.0D9
    ENDDO
    READ(500,*)
    ENDDO
ENDDO

!write(*,*) 'Header'
!write(*,*) (nsetable2_head(0,0,k2), k2 = 1, 3)
!write(*,*) (nsetable2_binde(0,0,0,k2), k2 = 1, 6)
!write(*,*) (nsetable2_head(0,1,k2), k2 = 1, 3)
!write(*,*) (nsetable2_binde(0,1,0,k2), k2 = 1, 6)
!write(*,*) (nsetable2_head(1,0,k2), k2 = 1, 3)
!write(*,*) (nsetable2_binde(1,0,0,k2), k2 = 1, 6)

CLOSE(500)

END SUBROUTINE read_nse_table2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
!
! This subroutine finds the mean atomic mass and mean atomic
! number of all grids. 
! Written by Leung Shing Chi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutines finds abar and zbar in all grid

subroutine find_AZbar()
use definition
use custom_def
implicit none

! Dummy variables
integer :: i, j, k

! Dummy variable
REAL*8 :: dummy

! if(debug_flag == 1) write(*,*) 'In FindAZBar'

! Simply call the subroutine to help you
! find the abar and zbar directly

do k = 1, nz, 1
    do j = 1, ny, 1
        do i = 1, nx, 1 
            CALL PRIVATE_HELMEOS_AZBAR(prim(ihe4:ini56,i,j,k), abar2(i,j,k), zbar2(i,j,k), dummy)
        enddo
    enddo
enddo

CALL BOUNDARY1D_NM(abar2, even, even, even, even, even, even)
CALL BOUNDARY1D_NM(zbar2, even, even, even, even, even, even)
CALL BOUNDARY1D_NM(prim(iye2,:,:,:), even, even, even, even, even, even)

END SUBROUTINE find_AZbar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
!
! This subroutine bridges the NUCLEAR_NET subroutine which 
! aims at finding the abar and zbar.
! Written by Leung Shin Chi in 2016.
! For more details about nuclear net, refer Timmes (2000b)
! or his website.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine private_helmeos_azbar(xiso_in, abar_out, zbar_out, electron_frac_out)
USE CUSTOM_DEF
implicit none

!For 7 isotopes only, expand the function when more isotopes are included
!The sequence is as follows: He4, C12, O16, Ne20, Mg24, Si28, Ni56

! Input quantities
real (selected_real_kind(15,307)), dimension(totalion):: xiso_in

! Output quantitites
real (selected_real_kind(15,307)) :: abar_out, zbar_out, electron_frac_out

! Dummy variables
real (selected_real_kind(15,307)) :: dummy
real (selected_real_kind(15,307)), dimension(totalion) :: dummy_x

! The traditional way
!ymass(1:ionmax) = xmass(1:ionmax)/wion(1:ionmax)      
!wbar  = 1.0d0/sum(ymass(1:ionmax))
!sum1  = sum(aion(1:ionmax)*ymass(1:ionmax))
!abar2  = wbar * sum1
!electron_frac2  = sum(zion(1:ionmax)*ymass(1:ionmax))
!zbar2  = wbar * ye

! Call the preset subroutine
call azbar(xiso_in,aion,zion,wion,totalion,dummy_x,abar_out,zbar_out,dummy,electron_frac_out,dummy)

end subroutine private_helmeos_azbar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
!
! This subroutine is obtained from Timmes' nuclear network
! program to calculate the abar and zbar for a given 
! chemical composition in a flexible structure
! Merged by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine azbar(xmass,aion,zion,wion,ionmax, &
                    ymass,abar,zbar,wbar,ye,nxcess)
use definition
USE CUSTOM_DEF
implicit none

! this routine calculates composition variables

! input:
! mass fractions               = xmass(1:ionmax)  dimensionless
! number of nucleons           = aion(1:ionmax)   dimensionless
! charge of nucleus            = zion(1:ionmax)   dimensionless
! atomic weight or molar mass  = wion(1:ionmax)    g/mole
! number of isotopes           = ionmax
!
! output:
! molar abundances        = ymass(1:ionmax)   mole/g
! mean number of nucleons = abar              dimensionless
! mean nucleon charge     = zbar              dimensionless
! mean weight             = wbar              g/mole
! electron fraction       = ye                mole/g
! neutron excess          = xcess


! declare the pass
integer          ionmax
real (selected_real_kind(15,307)), dimension(1:totalion) :: xmass,aion,zion,wion,ymass
real (selected_real_kind(15,307)) :: abar,zbar,wbar,ye,nxcess


! local variables
real (selected_real_kind(15,307)) asum,sum1

! molar abundances
ymass(1:totalion) = xmass(1:totalion)/wion(1:totalion)

! mean molar mass
wbar  = 1.0d0/sum(ymass(1:totalion))

! mean number of nucleons
sum1  = sum(aion(1:totalion)*ymass(1:totalion))
abar  = wbar * sum1

! mean charge
ye  = sum(zion(1:totalion)*ymass(1:totalion))
zbar  = wbar * ye

! neutron excess
nxcess = sum1 - 2.0d0 * ye

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine findhelmtemp()
use definition
use CUSTOM_DEF
implicit none

! The maximum and minimum allowed epsilon
REAL*8 :: epsilon_temp_max, epsilon_temp_min

! Dummy variables
REAL*8 :: dummy
REAL*8 :: abar_dum, zbar_dum

! Dummy variables
INTEGER :: i, j, k

! Helmholtz flag
INTEGER :: flag_notfindtemp

! reply when unsuccessful temp-finding occurs
! But I do not use it anymore
INTEGER :: flag_ans

DO k = 1, nz, 1
    DO j = 1, ny, 1
        DO i = 1, nx, 1
            IF(prim(irho,i,j,k) >= prim_a(irho)) THEN
                ! reset the flag
                flag_notfindtemp = 0

                ! Only grid with density above threshold density threshold is counted
                call private_invert_helm_ed(epsilon(i,j,k), &
                            prim(irho,i,j,k), abar2(i,j,k), &
                            zbar2(i,j,k), prim(iye2,i,j,k), &
                                            temp2_old(i,j,k), temp2(i,j,k), flag_notfindtemp)
                
                IF (invert == 0 .and. flag_notfindtemp == 1) THEN
                    invert = 1
                    ! STOP
                ENDIF
            
                IF(temp2(i,j,k) < temp_min) THEN
                    call HELM_EOSEPSILON(prim(irho,i,j,k), temp_min, abar2(i,j,k), zbar2(i,j,k), prim(iye2,i,j,k), epsilon_temp_min)
                    temp2(i,j,k) = temp_min
                    epsilon(i,j,k) = epsilon_temp_min
                ELSEIF(temp2(i,j,k) > temp_max) THEN
                    call HELM_EOSEPSILON(prim(irho,i,j,k), temp_max, abar2(i,j,k), zbar2(i,j,k), prim(iye2,i,j,k), epsilon_temp_max)
                        temp2(i,j,k) = temp_max
                        epsilon(i,j,k) = epsilon_temp_max
                ENDIF
            ! ELSE
            !     temp2(i,j,k) = temp2_a
            !     epsilon(i,j,k) = eps_a
            ENDIF
        ENDDO
    ENDDO
ENDDO
temp2_old = temp2
CALL BOUNDARY1D_NM(temp2, even, even, even, even, even, even)
CALL BOUNDARY1D_NM(epsilon, even, even, even, even, even, even)
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HELM_EOSEPSILON(rho_in, temp_in, abar_in, zbar_in, ye_in, epsilon_out)
USE IEEE_ARITHMETIC
INCLUDE 'implno.dek'
INCLUDE 'vector_eos.dek'


! Input quantities
REAL (selected_real_kind(15,307)) :: rho_in, temp_in, abar_in, zbar_in, ye_in

! Output quantity
REAL (selected_real_kind(15,307)) :: epsilon_out


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert the data from code unit to cgs and input the blanks

den_row(1) = rho_in * 6.1710D17
temp_row(1) = temp_in * 1.0D9
abar_row(1) = abar_in
! zbar_row(1) = zbar_in
zbar_row(1) = MIN(zbar_in, ye_in * abar_in)
jlo_eos = 1 
jhi_eos = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Call the EOS

call helmeos


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the results back to code unit

epsilon_out = etot_row(1) * 1.0D-4 / 9.0D16

IF (ieee_is_nan(epsilon_out)) THEN
    WRITE(*,*) 'Helmeos: epsilon is nan, input in cgs is'
    WRITE(*,*) 'Density', den_row(1)
    WRITE(*,*) 'Temperature', temp_row(1)
    WRITE(*,*) 'Abar', abar_row(1)
    WRITE(*,*) 'Zbar', zbar_row(1)
    STOP
ENDIF

END SUBROUTINE HELM_EOSEPSILON

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HELM_EOSSOUNDSPEED(rho_in, temp_in, abar_in, zbar_in, cs_out)
USE CUSTOM_DEF
USE DEFINITION
USE IEEE_ARITHMETIC
include 'implno.dek'
include 'vector_eos.dek'

! Input quantitites           
real (selected_real_kind(15,307)):: rho_in, temp_in, abar_in, zbar_in

! Output quantity
real (selected_real_kind(15,307)):: cs_out

! Signal Flag
integer :: flag_eostable
                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert data from code unit to cgs and input the necessary blanks
        
temp_row(1) = temp_in * 1.0D9
den_row(1) = rho_in * 6.1710D17   
abar_row(1) = abar_in
zbar_row(1) = zbar_in
! zbar_row(1) = MIN(zbar_in, ye_in * abar_in)
jlo_eos = 1
jhi_eos = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Call the EOS

call helmeos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output to code unit

cs_out = cs_row(1) / 3.0D10

IF (ieee_is_nan(cs_out)) THEN
    WRITE(*,*) 'Speed of sound is nan'
    STOP
ENDIF

END SUBROUTINE HELM_EOSSOUNDSPEED

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HELM_EOSPRESSURE(rho_in, temp_in, abar_in, zbar_in, ye_in, p_out, dpdrho_out, dpdeps_out, flag_eostable)
USE IEEE_ARITHMETIC
IMPLICIT NONE
INCLUDE 'vector_eos.dek'

REAL (selected_real_kind(15,307)) :: rho_in, temp_in, abar_in, zbar_in, ye_in
REAL (selected_real_kind(15,307)) :: p_out, dpdrho_out, dpdeps_out
INTEGER :: flag_eostable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert data from code unit to cgs and input the necessary blanks

den_row(1) = rho_in * 6.1710D17
temp_row(1) = temp_in * 1.0D9
abar_row(1) = abar_in
zbar_row(1) = MIN(zbar_in, ye_in * abar_in)
jlo_eos = 1 
jhi_eos = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Call the EOS      

CALL helmeos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert the results back to code unit

p_out = ptot_row(1) / 5.5539D38
dpdrho_out = dpd_row(1) / 9.0D20
dpdeps_out = (dpt_row(1) / det_row(1))/6.1710D17   !den_row(1) * (gam1_row(1) - 1.0D0) / 6.1710D17
IF(eosfail .eqv. .true.) flag_eostable = 0

IF (ieee_is_nan(p_out)) THEN
    WRITE(*,*) 'den_row(1)', den_row(1)
    WRITE(*,*) 'temp_row(1)', temp_row(1)
    WRITE(*,*) 'abar_row(1)', abar_row(1)
    WRITE(*,*) 'zbar_row(1)', zbar_row(1)
    WRITE(*,*) 'Pressure is nan'
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine private_invert_helm_ed(epsilon_input, rho_input, abar_input, zbar_input, ye_input, temp_input, temp_output, &
                                flag_notfindtemp)
USE IEEE_ARITHMETIC
include 'implno.dek'
include 'const.dek'
include 'vector_eos.dek'


! given the specific internal energy density, density, and composition
! find everything else

! it is assumed that etot_row(j), den_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input temp_row(j) conatins a guess for the temperature,
! on output temp_row(j) contains the converged temperature.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.

! extra variables for merging with the outside subroutine 

    integer          			   :: flag_notfindtemp
    double precision			   :: epsilon_input, rho_input, abar_input, &
                        zbar_input, ye_input, temp_input, temp_output



! local variables
    integer          i,j,jlo_save,jhi_save
    double precision tmpold,tmp,f,df,tmpnew,eostol,fpmin
    parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)



! initialize

    jlo_eos = 1
    jhi_eos = 1
    etot_row(1) 	= epsilon_input * 9.0D20
    temp_row(1) 	= temp_input * 1.0D9
    den_row(1) 	= rho_input * 6.1710D17 
    abar_row(1) 	= abar_input
    zbar_row(1) 	= zbar_input
    !zbar_row(1) 	= MIN(zbar_input(1), ye_input(1) * abar_input(1))
    
    jlo_save = jlo_eos
    jhi_save = jhi_eos
    do j=jlo_eos, jhi_eos
    eoswrk01(j) = 0.0d0
    eoswrk02(j) = 0.0d0
    eoswrk03(j) = etot_row(j)
    eoswrk04(j) = temp_row(j)
    end do


! do the first newton loop with all elements in the pipe

    call helmeos

    do j = jlo_eos, jhi_eos

    f     = etot_row(j)/eoswrk03(j) - 1.0d0
    df    = det_row(j)/eoswrk03(j)

    eoswrk02(j) = f/df

! limit excursions to factor of two changes
    tmp    = temp_row(j)
    tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
    eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new temperature, keep it within the table limits
    temp_row(j)  = min(1.0d14,max(tmpnew,1.0d-11))


    enddo



! now loop over each element of the pipe individually

    do j = jlo_save, jhi_save

    do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

        call helmeos

        f     = etot_row(j)/eoswrk03(j) - 1.0d0
        df    = det_row(j)/eoswrk03(j)

        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        tmp    = temp_row(j)
        tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
        eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new temperature, keep it within the table limits
        temp_row(j)  = min(1.0d13,max(tmpnew,1.0d3))


! end of netwon loop
    end do


! we did not converge if we land here
!      write(6,*)
!      write(6,*) 'newton-raphson failed in routine invert_helm_ed'
!      write(6,*) 'pipeline element',j
!      write(6,01) 'ewant  =',eoswrk03(j)
! 01   format(1x,5(a,1pe16.8))
!      write(6,01) 'error =',eoswrk01(j), &
!                  '  eostol=',eostol,'  fpmin =',fpmin
!      write(6,01) 'tmp   =',temp_row(j),'  tmpold=',eoswrk04(j)
!      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
!      write(6,*)

!     ! stop 'could not find a temperature in routine invert_helm_ed'
!     write(*,*) 'attnetion! cannot find temperature'
    flag_notfindtemp = 1


! land here if newton loop converged, back for another pipe element
20    continue
    end do



! call eos one more time with the converged value of the density

    jlo_eos = jlo_save
    jhi_eos = jhi_save

    call helmeos

    temp_output = temp_row(1) * 1.0D-9

    IF (ieee_is_nan(temp_output)) THEN
        WRITE(*,*) 'Temperature is nan'
        STOP
    ENDIF
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine helmeos
include 'implno.dek'
include 'helm_const.dek'
include 'vector_eos.dek'
include 'helm_table_storage.dek'


! given a temperature temp [K], density den [g/cm**3], and a composition
! characterized by abar and zbar, this routine returns most of the other
! thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
! specific thermal energy [erg/gr], the entropy [erg/g/K], along with
! their derivatives with respect to temperature, density, abar, and zbar.
! other quantites such the normalized chemical potential eta (plus its
! derivatives), number density of electrons and positron pair (along
! with their derivatives), adiabatic indices, specific heats, and
! relativistically correct sound speed are also returned.
!
! this routine assumes planckian photons, an ideal gas of ions,
! and an electron-positron gas with an arbitrary degree of relativity
! and degeneracy. interpolation in a table of the helmholtz free energy
! is used to return the electron-positron thermodynamic quantities.
! all other derivatives are analytic.
!
! references: cox & giuli chapter 24 ; timmes & swesty apj 1999


! declare
    integer          i,j
    double precision temp,den,abar,zbar,ytot1,ye, &
                    x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
                    dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
                    dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
                    deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
                    dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
                    sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
                    dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
                    gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
                    detadt,detadd,xnefer,dxnedt,dxnedd,s

    double precision pgas,dpgasdd,dpgasdt,dpgasda,dpgasdz, &
                    egas,degasdd,degasdt,degasda,degasdz, &
                    sgas,dsgasdd,dsgasdt,dsgasda,dsgasdz, &
                    cv_gas,cp_gas,gam1_gas,gam2_gas,gam3_gas, &
                    chit_gas,chid_gas,nabad_gas,sound_gas


    double precision sioncon,forth,forpi,kergavo,ikavo,asoli3,light2
    parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h), &
                    forth   = 4.0d0/3.0d0, &
                    forpi   = 4.0d0 * pi, &
                    kergavo = kerg * avo, &
                    ikavo   = 1.0d0/kergavo, &
                    asoli3  = asol/3.0d0, &
                    light2  = clight * clight)

! for the abar derivatives
    double precision dpradda,deradda,dsradda, &
                    dpionda,deionda,dsionda, &
                    dpepda,deepda,dsepda, &
                    dpresda,denerda,dentrda, &
                    detada,dxneda

! for the zbar derivatives
    double precision dpraddz,deraddz,dsraddz, &
                    dpiondz,deiondz,dsiondz, &
                    dpepdz,deepdz,dsepdz, &
                    dpresdz,denerdz,dentrdz, &
                    detadz,dxnedz

! for the interpolations
    integer          iat,jat
    double precision free,df_d,df_t,df_dd,df_tt,df_dt
    double precision xt,xd,mxt,mxd, &
                    si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                    si0d,si1d,si2d,si0md,si1md,si2md, &
                    dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                    dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                    ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
                    ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md, &
                    z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, &
                    dpsi2,ddpsi2,din,h5,fi(36), &
                    xpsi0,xdpsi0,xpsi1,xdpsi1,h3, &
                    w0t,w1t,w2t,w0mt,w1mt,w2mt, &
                    w0d,w1d,w2d,w0md,w1md,w2md


! for the uniform background coulomb correction
    double precision dsdd,dsda,lami,inv_lami,lamida,lamidd, &
                    plasg,plasgdd,plasgdt,plasgda,plasgdz, &
                    ecoul,decouldd,decouldt,decoulda,decouldz, &
                    pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                    scoul,dscouldd,dscouldt,dscoulda,dscouldz, &
                    a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
    parameter        (a1    = -0.898004d0, &
                    b1    =  0.96786d0, &
                    c1    =  0.220703d0, &
                    d1    = -0.86097d0, &
                    e1    =  2.5269d0, &
                    a2    =  0.29561d0, &
                    b2    =  1.9885d0, &
                    c2    =  0.288675d0, &
                    third =  1.0d0/3.0d0, &
                    esqu  =  qe * qe)


! quintic hermite polynomial statement functions
! psi0 and its derivatives
    psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
    dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
    ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)


! psi1 and its derivatives
    psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
    dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
    ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)


! psi2  and its derivatives
    psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
    dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
    ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)


! biquintic hermite polynomial statement function
    h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= &
            fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
        + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
        + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
        + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
        + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
        + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
        + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
        + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
        + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
        + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
        + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
        + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
        + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
        + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
        + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
        + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
        + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
        + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



! cubic hermite polynomial statement functions
! psi0 & derivatives
    xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
    xdpsi0(z) = z * (6.0d0*z - 6.0d0)


! psi1 & derivatives
    xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
    xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


! bicubic hermite polynomial statement function
    h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = &
            fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
        + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
        + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
        + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
        + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
        + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
        + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
        + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



! popular format statements
01    format(1x,5(a,1pe11.3))
02    format(1x,a,1p4e16.8)
03    format(1x,4(a,1pe11.3))
04    format(1x,4(a,i4))



! start of pipeline loop, normal execution starts here
    eosfail = .false.
    do j=jlo_eos,jhi_eos

!       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
!       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

    temp  = temp_row(j)
    den   = den_row(j)
    abar  = abar_row(j)
    zbar  = zbar_row(j)
    ytot1 = 1.0d0/abar
    ye    = max(1.0d-16,ytot1 * zbar)


! initialize
    deni    = 1.0d0/den
    tempi   = 1.0d0/temp
    kt      = kerg * temp
    ktinv   = 1.0d0/kt


! radiation section:
    prad    = asoli3 * temp * temp * temp * temp
    dpraddd = 0.0d0
    dpraddt = 4.0d0 * prad*tempi
    dpradda = 0.0d0
    dpraddz = 0.0d0

    erad    = 3.0d0 * prad*deni
    deraddd = -erad*deni
    deraddt = 3.0d0 * dpraddt*deni
    deradda = 0.0d0
    deraddz = 0.0d0

    srad    = (prad*deni + erad)*tempi
    dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
    dsraddt = (dpraddt*deni + deraddt - srad)*tempi
    dsradda = 0.0d0
    dsraddz = 0.0d0


! ion section:
    xni     = avo * ytot1 * den
    dxnidd  = avo * ytot1
    dxnida  = -xni * ytot1

    pion    = xni * kt
    dpiondd = dxnidd * kt
    dpiondt = xni * kerg
    dpionda = dxnida * kt
    dpiondz = 0.0d0

    eion    = 1.5d0 * pion*deni
    deiondd = (1.5d0 * dpiondd - eion)*deni
    deiondt = 1.5d0 * dpiondt*deni
    deionda = 1.5d0 * dpionda*deni
    deiondz = 0.0d0


! sackur-tetrode equation for the ion entropy of
! a single ideal gas characterized by abar
    x       = abar*abar*sqrt(abar) * deni/avo
    s       = sioncon * temp
    z       = x * s * sqrt(s)
    y       = log(z)

!        y       = 1.0d0/(abar*kt)
!        yy      = y * sqrt(y)
!        z       = xni * sifac * yy
!        etaion  = log(z)


    sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
    dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
                - kergavo * deni * ytot1
    dsiondt = (dpiondt*deni + deiondt)*tempi - &
                (pion*deni + eion) * tempi*tempi &
                + 1.5d0 * kergavo * tempi*ytot1
    x       = avo*kerg/abar
    dsionda = (dpionda*deni + deionda)*tempi &
                + kergavo*ytot1*ytot1* (2.5d0 - y)
    dsiondz = 0.0d0



! electron-positron section:


! assume complete ionization
    xnem    = xni * zbar


! enter the table with ye*den
    din = ye*den


! bomb proof the input
    ! if (temp .gt. t(jmax)) then
    !  write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
    !  write(6,*) 'temp too hot, off grid'
    !  write(6,*) 'setting eosfail to true and returning'
    !  eosfail = .true.
    !  !return
    ! end if
    ! if (temp .lt. t(1)) then
    !  write(6,01) 'temp=',temp,' t(1)=',t(1)
    !  write(6,*) 'temp too cold, off grid'
    !  write(6,*) 'setting eosfail to true and returning'
    !  eosfail = .true.
    !  !return
    ! end if
    ! if (din  .gt. d(imax)) then
    !  write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
    !  write(6,*) 'ye*den too big, off grid'
    !  write(6,*) 'setting eosfail to true and returning'
    !  eosfail = .true.
    !  !return
    ! end if
    ! if (din  .lt. d(1)) then
    !  write(6,01) 'ye*den=',din,' d(1)=',d(1)
    !  write(6,*) 'ye*den too small, off grid'
    !  write(6,*) 'setting eosfail to true and returning'
    !  eosfail = .true.
    !  !return
    ! end if

! hash locate this temperature and density
    jat = int((log10(temp) - tlo)*tstpi) + 1
    jat = max(1,min(jat,jmax-1))
    iat = int((log10(din) - dlo)*dstpi) + 1
    iat = max(1,min(iat,imax-1))


! access the table locations only once
    fi(1)  = f(iat,jat)
    fi(2)  = f(iat+1,jat)
    fi(3)  = f(iat,jat+1)
    fi(4)  = f(iat+1,jat+1)
    fi(5)  = ft(iat,jat)
    fi(6)  = ft(iat+1,jat)
    fi(7)  = ft(iat,jat+1)
    fi(8)  = ft(iat+1,jat+1)
    fi(9)  = ftt(iat,jat)
    fi(10) = ftt(iat+1,jat)
    fi(11) = ftt(iat,jat+1)
    fi(12) = ftt(iat+1,jat+1)
    fi(13) = fd(iat,jat)
    fi(14) = fd(iat+1,jat)
    fi(15) = fd(iat,jat+1)
    fi(16) = fd(iat+1,jat+1)
    fi(17) = fdd(iat,jat)
    fi(18) = fdd(iat+1,jat)
    fi(19) = fdd(iat,jat+1)
    fi(20) = fdd(iat+1,jat+1)
    fi(21) = fdt(iat,jat)
    fi(22) = fdt(iat+1,jat)
    fi(23) = fdt(iat,jat+1)
    fi(24) = fdt(iat+1,jat+1)
    fi(25) = fddt(iat,jat)
    fi(26) = fddt(iat+1,jat)
    fi(27) = fddt(iat,jat+1)
    fi(28) = fddt(iat+1,jat+1)
    fi(29) = fdtt(iat,jat)
    fi(30) = fdtt(iat+1,jat)
    fi(31) = fdtt(iat,jat+1)
    fi(32) = fdtt(iat+1,jat+1)
    fi(33) = fddtt(iat,jat)
    fi(34) = fddtt(iat+1,jat)
    fi(35) = fddtt(iat,jat+1)
    fi(36) = fddtt(iat+1,jat+1)


! various differences
    xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
    xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
    mxt = 1.0d0 - xt
    mxd = 1.0d0 - xd

! the six density and six temperature basis functions
    si0t =   psi0(xt)
    si1t =   psi1(xt)*dt_sav(jat)
    si2t =   psi2(xt)*dt2_sav(jat)

    si0mt =  psi0(mxt)
    si1mt = -psi1(mxt)*dt_sav(jat)
    si2mt =  psi2(mxt)*dt2_sav(jat)

    si0d =   psi0(xd)
    si1d =   psi1(xd)*dd_sav(iat)
    si2d =   psi2(xd)*dd2_sav(iat)

    si0md =  psi0(mxd)
    si1md = -psi1(mxd)*dd_sav(iat)
    si2md =  psi2(mxd)*dd2_sav(iat)

! derivatives of the weight functions
    dsi0t =   dpsi0(xt)*dti_sav(jat)
    dsi1t =   dpsi1(xt)
    dsi2t =   dpsi2(xt)*dt_sav(jat)

    dsi0mt = -dpsi0(mxt)*dti_sav(jat)
    dsi1mt =  dpsi1(mxt)
    dsi2mt = -dpsi2(mxt)*dt_sav(jat)

    dsi0d =   dpsi0(xd)*ddi_sav(iat)
    dsi1d =   dpsi1(xd)
    dsi2d =   dpsi2(xd)*dd_sav(iat)

    dsi0md = -dpsi0(mxd)*ddi_sav(iat)
    dsi1md =  dpsi1(mxd)
    dsi2md = -dpsi2(mxd)*dd_sav(iat)

! second derivatives of the weight functions
    ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
    ddsi1t =   ddpsi1(xt)*dti_sav(jat)
    ddsi2t =   ddpsi2(xt)

    ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
    ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
    ddsi2mt =  ddpsi2(mxt)

!        ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
!        ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
!        ddsi2d =   ddpsi2(xd)

!        ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
!        ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
!        ddsi2md =  ddpsi2(mxd)


! the free energy
    free  = h5(iat,jat, &
            si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
            si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to density
    df_d  = h5(iat,jat, &
            si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
            dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)


! derivative with respect to temperature
    df_t = h5(iat,jat, &
            dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
            si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to density**2
!        df_dd = h5(iat,jat,
!     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
!     2          ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

! derivative with respect to temperature**2
    df_tt = h5(iat,jat, &
            ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
            si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to temperature and density
    df_dt = h5(iat,jat, &
            dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
            dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



! now get the pressure derivative with density, chemical potential, and
! electron positron number densities
! get the interpolation weight functions
    si0t   =  xpsi0(xt)
    si1t   =  xpsi1(xt)*dt_sav(jat)

    si0mt  =  xpsi0(mxt)
    si1mt  =  -xpsi1(mxt)*dt_sav(jat)

    si0d   =  xpsi0(xd)
    si1d   =  xpsi1(xd)*dd_sav(iat)

    si0md  =  xpsi0(mxd)
    si1md  =  -xpsi1(mxd)*dd_sav(iat)


! derivatives of weight functions
    dsi0t  = xdpsi0(xt)*dti_sav(jat)
    dsi1t  = xdpsi1(xt)

    dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
    dsi1mt = xdpsi1(mxt)

    dsi0d  = xdpsi0(xd)*ddi_sav(iat)
    dsi1d  = xdpsi1(xd)

    dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
    dsi1md = xdpsi1(mxd)


! look in the pressure derivative only once
    fi(1)  = dpdf(iat,jat)
    fi(2)  = dpdf(iat+1,jat)
    fi(3)  = dpdf(iat,jat+1)
    fi(4)  = dpdf(iat+1,jat+1)
    fi(5)  = dpdft(iat,jat)
    fi(6)  = dpdft(iat+1,jat)
    fi(7)  = dpdft(iat,jat+1)
    fi(8)  = dpdft(iat+1,jat+1)
    fi(9)  = dpdfd(iat,jat)
    fi(10) = dpdfd(iat+1,jat)
    fi(11) = dpdfd(iat,jat+1)
    fi(12) = dpdfd(iat+1,jat+1)
    fi(13) = dpdfdt(iat,jat)
    fi(14) = dpdfdt(iat+1,jat)
    fi(15) = dpdfdt(iat,jat+1)
    fi(16) = dpdfdt(iat+1,jat+1)

! pressure derivative with density
    dpepdd  = h3(iat,jat, &
                    si0t,   si1t,   si0mt,   si1mt, &
                    si0d,   si1d,   si0md,   si1md)
    dpepdd  = max(ye * dpepdd,1.0d-30)



! look in the electron chemical potential table only once
    fi(1)  = ef(iat,jat)
    fi(2)  = ef(iat+1,jat)
    fi(3)  = ef(iat,jat+1)
    fi(4)  = ef(iat+1,jat+1)
    fi(5)  = eft(iat,jat)
    fi(6)  = eft(iat+1,jat)
    fi(7)  = eft(iat,jat+1)
    fi(8)  = eft(iat+1,jat+1)
    fi(9)  = efd(iat,jat)
    fi(10) = efd(iat+1,jat)
    fi(11) = efd(iat,jat+1)
    fi(12) = efd(iat+1,jat+1)
    fi(13) = efdt(iat,jat)
    fi(14) = efdt(iat+1,jat)
    fi(15) = efdt(iat,jat+1)
    fi(16) = efdt(iat+1,jat+1)


! electron chemical potential etaele
    etaele  = h3(iat,jat, &
                    si0t,   si1t,   si0mt,   si1mt, &
                    si0d,   si1d,   si0md,   si1md)


! derivative with respect to density
    x       = h3(iat,jat, &
                    si0t,   si1t,   si0mt,   si1mt, &
                dsi0d,  dsi1d,  dsi0md,  dsi1md)
    detadd  = ye * x

! derivative with respect to temperature
    detadt  = h3(iat,jat, &
                dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                    si0d,   si1d,   si0md,   si1md)

! derivative with respect to abar and zbar
    detada = -x * din * ytot1
    detadz =  x * den * ytot1



! look in the number density table only once
    fi(1)  = xf(iat,jat)
    fi(2)  = xf(iat+1,jat)
    fi(3)  = xf(iat,jat+1)
    fi(4)  = xf(iat+1,jat+1)
    fi(5)  = xft(iat,jat)
    fi(6)  = xft(iat+1,jat)
    fi(7)  = xft(iat,jat+1)
    fi(8)  = xft(iat+1,jat+1)
    fi(9)  = xfd(iat,jat)
    fi(10) = xfd(iat+1,jat)
    fi(11) = xfd(iat,jat+1)
    fi(12) = xfd(iat+1,jat+1)
    fi(13) = xfdt(iat,jat)
    fi(14) = xfdt(iat+1,jat)
    fi(15) = xfdt(iat,jat+1)
    fi(16) = xfdt(iat+1,jat+1)

! electron + positron number densities
    xnefer   = h3(iat,jat, &
                    si0t,   si1t,   si0mt,   si1mt, &
                    si0d,   si1d,   si0md,   si1md)

! derivative with respect to density
    x        = h3(iat,jat, &
                    si0t,   si1t,   si0mt,   si1mt, &
                dsi0d,  dsi1d,  dsi0md,  dsi1md)
    x = max(x,1.0d-30)
    dxnedd   = ye * x

! derivative with respect to temperature
    dxnedt   = h3(iat,jat, &
                dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                    si0d,   si1d,   si0md,   si1md)

! derivative with respect to abar and zbar
    dxneda = -x * din * ytot1
    dxnedz =  x  * den * ytot1


! the desired electron-positron thermodynamic quantities

! dpepdd at high temperatures and low densities is below the
! floating point limit of the subtraction of two large terms.
! since dpresdd doesn't enter the maxwell relations at all, use the
! bicubic interpolation done above instead of the formally correct expression
    x       = din * din
    pele    = x * df_d
    dpepdt  = x * df_dt
!        dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
    s       = dpepdd/ye - 2.0d0 * din * df_d
    dpepda  = -ytot1 * (2.0d0 * pele + s * din)
    dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)


    x       = ye * ye
    sele    = -df_t * ye
    dsepdt  = -df_tt * ye
    dsepdd  = -df_dt * x
    dsepda  = ytot1 * (ye * df_dt * din - sele)
    dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)


    eele    = ye*free + temp * sele
    deepdt  = temp * dsepdt
    deepdd  = x * df_d + temp * dsepdd
    deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
    deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz




! coulomb section:

! uniform background corrections only
! from yakovlev & shalybkov 1989
! lami is the average ion seperation
! plasg is the plasma coupling parameter

    z        = forth * pi
    s        = z * xni
    dsdd     = z * dxnidd
    dsda     = z * dxnida

    lami     = 1.0d0/s**third
    inv_lami = 1.0d0/lami
    z        = -third * lami
    lamidd   = z * dsdd/s
    lamida   = z * dsda/s

    plasg    = zbar*zbar*esqu*ktinv*inv_lami
    z        = -plasg * inv_lami
    plasgdd  = z * lamidd
    plasgda  = z * lamida
    plasgdt  = -plasg*ktinv * kerg
    plasgdz  = 2.0d0 * plasg/zbar


! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
    if (plasg .ge. 1.0) then
        x        = plasg**(0.25d0)
        y        = avo * ytot1 * kerg
        ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
        pcoul    = third * den * ecoul
        scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
                + d1 * (log(plasg) - 1.0d0) - e1)

        y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
        decouldd = y * plasgdd
        decouldt = y * plasgdt + ecoul/temp
        decoulda = y * plasgda - ecoul/abar
        decouldz = y * plasgdz

        y        = third * den
        dpcouldd = third * ecoul + y*decouldd
        dpcouldt = y * decouldt
        dpcoulda = y * decoulda
        dpcouldz = y * decouldz


        y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
        dscouldd = y * plasgdd
        dscouldt = y * plasgdt
        dscoulda = y * plasgda - scoul/abar
        dscouldz = y * plasgdz


! yakovlev & shalybkov 1989 equations 102, 103, 104
    else if (plasg .lt. 1.0) then
        x        = plasg*sqrt(plasg)
        y        = plasg**b2
        z        = c2 * x - third * a2 * y
        pcoul    = -pion * z
        ecoul    = 3.0d0 * pcoul/den
        scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

        s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
        dpcouldd = -dpiondd*z - pion*s*plasgdd
        dpcouldt = -dpiondt*z - pion*s*plasgdt
        dpcoulda = -dpionda*z - pion*s*plasgda
        dpcouldz = -dpiondz*z - pion*s*plasgdz

        s        = 3.0d0/den
        decouldd = s * dpcouldd - ecoul/den
        decouldt = s * dpcouldt
        decoulda = s * dpcoulda
        decouldz = s * dpcouldz

        s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
        dscouldd = s * plasgdd
        dscouldt = s * plasgdt
        dscoulda = s * plasgda - scoul/abar
        dscouldz = s * plasgdz
    end if


! bomb proof
    x   = prad + pion + pele + pcoul
    y   = erad + eion + eele + ecoul
    z   = srad + sion + sele + scoul

!        write(6,*) x,y,z
    if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then
!        if (x .le. 0.0 .or. y .le. 0.0) then
!        if (x .le. 0.0) then

!         write(6,*)
!         write(6,*) 'coulomb corrections are causing a negative pressure'
!         write(6,*) 'setting all coulomb corrections to zero'
!         write(6,*)

        pcoul    = 0.0d0
        dpcouldd = 0.0d0
        dpcouldt = 0.0d0
        dpcoulda = 0.0d0
        dpcouldz = 0.0d0
        ecoul    = 0.0d0
        decouldd = 0.0d0
        decouldt = 0.0d0
        decoulda = 0.0d0
        decouldz = 0.0d0
        scoul    = 0.0d0
        dscouldd = 0.0d0
        dscouldt = 0.0d0
        dscoulda = 0.0d0
        dscouldz = 0.0d0
    end if


! sum all the gas components
    pgas    = pion + pele + pcoul
    egas    = eion + eele + ecoul
    sgas    = sion + sele + scoul

    dpgasdd = dpiondd + dpepdd + dpcouldd
    dpgasdt = dpiondt + dpepdt + dpcouldt
    dpgasda = dpionda + dpepda + dpcoulda
    dpgasdz = dpiondz + dpepdz + dpcouldz

    degasdd = deiondd + deepdd + decouldd
    degasdt = deiondt + deepdt + decouldt
    degasda = deionda + deepda + decoulda
    degasdz = deiondz + deepdz + decouldz

    dsgasdd = dsiondd + dsepdd + dscouldd
    dsgasdt = dsiondt + dsepdt + dscouldt
    dsgasda = dsionda + dsepda + dscoulda
    dsgasdz = dsiondz + dsepdz + dscouldz




! add in radiation to get the total
    pres    = prad + pgas
    ener    = erad + egas
    entr    = srad + sgas

    dpresdd = dpraddd + dpgasdd
    dpresdt = dpraddt + dpgasdt
    dpresda = dpradda + dpgasda
    dpresdz = dpraddz + dpgasdz

    denerdd = deraddd + degasdd
    denerdt = deraddt + degasdt
    denerda = deradda + degasda
    denerdz = deraddz + degasdz

    dentrdd = dsraddd + dsgasdd
    dentrdt = dsraddt + dsgasdt
    dentrda = dsradda + dsgasda
    dentrdz = dsraddz + dsgasdz


! for the gas
! the temperature and density exponents (c&g 9.81 9.82)
! the specific heat at constant volume (c&g 9.92)
! the third adiabatic exponent (c&g 9.93)
! the first adiabatic exponent (c&g 9.97)
! the second adiabatic exponent (c&g 9.105)
! the specific heat at constant pressure (c&g 9.98)
! and relativistic formula for the sound speed (c&g 14.29)

    zz        = pgas*deni
    zzi       = den/pgas
    chit_gas  = temp/pgas * dpgasdt
    chid_gas  = dpgasdd*zzi
    cv_gas    = degasdt
    x         = zz * chit_gas/(temp * cv_gas)
    gam3_gas  = x + 1.0d0
    gam1_gas  = chit_gas*x + chid_gas
    nabad_gas = x/gam1_gas
    gam2_gas  = 1.0d0/(1.0d0 - nabad_gas)
    cp_gas    = cv_gas * gam1_gas/chid_gas
    z         = 1.0d0 + (egas + light2)*zzi
    sound_gas = clight * sqrt(gam1_gas/z)



! for the totals
    zz    = pres*deni
    zzi   = den/pres
    chit  = temp/pres * dpresdt
    chid  = dpresdd*zzi
    cv    = denerdt
    x     = zz * chit/(temp * cv)
    gam3  = x + 1.0d0
    gam1  = chit*x + chid
    nabad = x/gam1
    gam2  = 1.0d0/(1.0d0 - nabad)
    cp    = cv * gam1/chid
    z     = 1.0d0 + (ener + light2)*zzi
    sound = clight * sqrt(gam1/z)



! maxwell relations; each is zero if the consistency is perfect
    x   = den * den

    dse = temp*dentrdt/denerdt - 1.0d0

    dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0

    dsp = -dentrdd*x/dpresdt - 1.0d0


! store this row
    ptot_row(j)   = pres
    dpt_row(j)    = dpresdt
    dpd_row(j)    = dpresdd
    dpa_row(j)    = dpresda
    dpz_row(j)    = dpresdz

    etot_row(j)   = ener
    det_row(j)    = denerdt
    ded_row(j)    = denerdd
    dea_row(j)    = denerda
    dez_row(j)    = denerdz

    stot_row(j)   = entr
    dst_row(j)    = dentrdt
    dsd_row(j)    = dentrdd
    dsa_row(j)    = dentrda
    dsz_row(j)    = dentrdz


    pgas_row(j)   = pgas
    dpgast_row(j) = dpgasdt
    dpgasd_row(j) = dpgasdd
    dpgasa_row(j) = dpgasda
    dpgasz_row(j) = dpgasdz

    egas_row(j)   = egas
    degast_row(j) = degasdt
    degasd_row(j) = degasdd
    degasa_row(j) = degasda
    degasz_row(j) = degasdz

    sgas_row(j)   = sgas
    dsgast_row(j) = dsgasdt
    dsgasd_row(j) = dsgasdd
    dsgasa_row(j) = dsgasda
    dsgasz_row(j) = dsgasdz


    prad_row(j)   = prad
    dpradt_row(j) = dpraddt
    dpradd_row(j) = dpraddd
    dprada_row(j) = dpradda
    dpradz_row(j) = dpraddz

    erad_row(j)   = erad
    deradt_row(j) = deraddt
    deradd_row(j) = deraddd
    derada_row(j) = deradda
    deradz_row(j) = deraddz

    srad_row(j)   = srad
    dsradt_row(j) = dsraddt
    dsradd_row(j) = dsraddd
    dsrada_row(j) = dsradda
    dsradz_row(j) = dsraddz


    pion_row(j)   = pion
    dpiont_row(j) = dpiondt
    dpiond_row(j) = dpiondd
    dpiona_row(j) = dpionda
    dpionz_row(j) = dpiondz

    eion_row(j)   = eion
    deiont_row(j) = deiondt
    deiond_row(j) = deiondd
    deiona_row(j) = deionda
    deionz_row(j) = deiondz

    sion_row(j)   = sion
    dsiont_row(j) = dsiondt
    dsiond_row(j) = dsiondd
    dsiona_row(j) = dsionda
    dsionz_row(j) = dsiondz

    xni_row(j)    = xni

    pele_row(j)   = pele
    ppos_row(j)   = 0.0d0
    dpept_row(j)  = dpepdt
    dpepd_row(j)  = dpepdd
    dpepa_row(j)  = dpepda
    dpepz_row(j)  = dpepdz

    eele_row(j)   = eele
    epos_row(j)   = 0.0d0
    deept_row(j)  = deepdt
    deepd_row(j)  = deepdd
    deepa_row(j)  = deepda
    deepz_row(j)  = deepdz

    sele_row(j)   = sele
    spos_row(j)   = 0.0d0
    dsept_row(j)  = dsepdt
    dsepd_row(j)  = dsepdd
    dsepa_row(j)  = dsepda
    dsepz_row(j)  = dsepdz

    xnem_row(j)   = xnem
    xne_row(j)    = xnefer
    dxnet_row(j)  = dxnedt
    dxned_row(j)  = dxnedd
    dxnea_row(j)  = dxneda
    dxnez_row(j)  = dxnedz
    xnp_row(j)    = 0.0d0
    zeff_row(j)   = zbar

    etaele_row(j) = etaele
    detat_row(j)  = detadt
    detad_row(j)  = detadd
    detaa_row(j)  = detada
    detaz_row(j)  = detadz
    etapos_row(j) = 0.0d0

    pcou_row(j)   = pcoul
    dpcout_row(j) = dpcouldt
    dpcoud_row(j) = dpcouldd
    dpcoua_row(j) = dpcoulda
    dpcouz_row(j) = dpcouldz

    ecou_row(j)   = ecoul
    decout_row(j) = decouldt
    decoud_row(j) = decouldd
    decoua_row(j) = decoulda
    decouz_row(j) = decouldz

    scou_row(j)   = scoul
    dscout_row(j) = dscouldt
    dscoud_row(j) = dscouldd
    dscoua_row(j) = dscoulda
    dscouz_row(j) = dscouldz

    plasg_row(j)  = plasg

    dse_row(j)    = dse
    dpe_row(j)    = dpe
    dsp_row(j)    = dsp

    cv_gas_row(j)    = cv_gas
    cp_gas_row(j)    = cp_gas
    gam1_gas_row(j)  = gam1_gas
    gam2_gas_row(j)  = gam2_gas
    gam3_gas_row(j)  = gam3_gas
    nabad_gas_row(j) = nabad_gas
    cs_gas_row(j)    = sound_gas

    cv_row(j)     = cv
    cp_row(j)     = cp
    gam1_row(j)   = gam1
    gam2_row(j)   = gam2
    gam3_row(j)   = gam3
    nabad_row(j)  = nabad
    cs_row(j)     = sound

! end of pipeline loop
    enddo
    !return
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
!
! This subroutine checks if the local isotope mass fraction sum
! to one. If not, normalize it. 
! Written by Leung Shing Chi in 2016.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkxisotope()
use definition
use custom_def
implicit none

! dummy variables
integer :: i, j, k, l

! local variables
real*8:: xiso_sum

! if(debug_flag == 1) write(*,*) 'In CheckXIsotope'

! First check if any grid has any unusual mass fraction

do i = ihe4, ini56, 1
    do l = 1, nz, 1
        do k = 1, ny, 1
            do j = 1, nx, 1
                if (prim(i,j,k,l) < 1.0D-30) prim(i,j,k,l) = 1.0D-30
            enddo
        enddo
    enddo
enddo


! Then do the sum and normalization


do l = 1, nz, 1
    do k = 1, ny, 1
        do j = 1, nx, 1
            
            xiso_sum = 0.0D0
            do i = ihe4, ini56, 1
                xiso_sum = xiso_sum + prim(i,j,k,l)
            enddo
            prim(ihe4:ini56,j,k,l) = prim(ihe4:ini56,j,k,l) / xiso_sum

        enddo
    enddo
enddo

CALL BOUNDARYP_NM

end subroutine checkxisotope

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine finds the chemical composition
! due to nuclear statistical equilibrium (NSE)
! Written by Leung Shing Chi in 2016, but 
! assuming the binding energy change due to 
! the composition difference causes gain/loss in 
! local energy density
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine NSE2
use definition                               
use custom_def
use ecaptable_module
implicit none                               

! Flag for finding temperature
integer :: flag_notfindtemp             

! Dummy variables
integer :: i, j, k, k2

! Dummy variables
real*8 :: ye_sample = 0.5D0, dummy

! Local variables
real*8 :: abar_mid, zbar_mid, ye_mid

! Input variables
real*8 :: temp_beg    
real*8 :: eps_beg

! Trial variables
real*8 :: temp_mid                        
real*8 :: eps_mid         
real*8 :: rho_mid

! Change in temperature
real*8 :: dtemp_mid 

! Number of successful trial                
integer  :: count_digit

! Ecap rate and neutrino energy loss
real*8 :: ecaprate_mid
real*8 :: eneurate_mid

! Initial and Expected chemical composition
real*8, dimension(totalion) :: x_mid
real*8, dimension(totalion) :: x_burn

! Binding energy
real*8 :: binde_bf, binde_af, deps_nuc

! Energy balance equation
real*8 :: check_e, check_e_last

! Timescale for burning to NSE
real*8 :: nse_burntime, temp_nse

! if(debug_flag == 1) write(*,*) 'In NSE2'

! Initialization
! Already done in burn_phase2b

do k = 1, nz, 1              
    do j = 1, nx, 1                    

        ! Do the checking first               
        if(temp2(j,1,k) < 5.0D0) cycle
        if(temp2(j,1,k) > 11.00) cycle
        if(prim(iye2,j,1,k) < 0.40D0) cycle
        if(prim(iye2,j,1,k) > 0.50D0) cycle
        if(prim(irho,j,1,k) < 1.62D-11) cycle
        if(burn_ratio(j,1,k) /= 1.0D0) cycle

        ! Special test if multi-stage burning is used
        !if(advburn_flag == 1) then
        if(nse_flag(j,1,k) == 1) cycle
        !endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! If they pass then start the iteration
        nse_flag(j,1,k) = 2                

        ! The density is not changed, so 
        ! the trial value is also the 
        ! final value
        rho_mid = prim(irho,j,1,k)

        ! Give a trial temperature
        temp_beg = temp2(j,1,k)                  
        eps_beg = epsilon(j,1,k)

        ! Also give the trial composition
        abar_mid = abar2(j,1,k)
        zbar_mid = zbar2(j,1,k)
        ye_mid = prim(iye2,j,1,k)                
        x_mid(:) = prim(ihe4:ini56,j,1,k)

        !if(flame_loc_ratio(j,1,k) == 1.0D0 .or. deton_loc_ratio(j,1,k) == 1.0D0) then
        if(burn_ratio(j,1,k) == 1.0D0) then

            ! Scheme for electron capture
            ! Get the Ecap rate and energy loss
            call getecaprate(rho_mid, temp_beg, prim(iye2,j,1,k), ecaprate_mid, eneurate_mid)
            !if(ye_mid < 0.46D0) then
            !   ecaprate_mid = 0.0D0; eneurate_mid = 0.0D0
            !endif


            ! Update the Ye rate
            ! Note: No iteration is done here because
            ! the temperature sensitivity of Ecap
            ! rate is much smaller than NSE composition
            ye_mid = prim(iye2,j,1,k) + ecaprate_mid * dt
            !ye_mid = prim(iye2,j,1,k)                 

            ! Compute the binding energy of the 
            ! initial composition
            call compute_binde(x_mid, binde_bf)

            !write(*,*) 'Temp Before: ', temp_beg
            !write(*,*) 'Abar before: ', abar2(j,1,k)
            !write(*,*) 'Zbar before: ', zbar2(j,1,k)         
            !write(*,*) 'Eps before: ', eps_beg
            !write(*,*) 'Binde before: ', binde_bf
            !write(*,*)                 

            ! Prepare for the search of temperature
            count_digit = 0
            temp_mid = temp_beg
            dtemp_mid = 0.01D0 * temp_beg          

            ! Some debug stuff
            !if(j == 1 .and. k == 1) write(*,*) 'Now in grid', j, k
            !if(j == 1 .and. k == 1) write(*,*) 'Input density = ', rho_mid
            !if(j == 1 .and. k == 1) write(*,*) 'Old temp = ', temp_beg
            !if(j == 1 .and. k == 1) write(*,101) (xiso(j,k,k2),k2=1,totalion)

                do i = 0, 400, 1   

                ! Now do the search

                    !call nse_interface(temp_mid ,rho_mid, ye_mid, x_burn)

                ! Get the trial NSE state by the trial temp
                    call getnsestate(rho_mid, temp_mid, x_burn)

                ! Compute the trial binding energy
                    call compute_binde(x_burn, binde_af)

                ! Calculate the trial abar and zbar
                    call private_helmeos_azbar(x_burn, abar_mid, zbar_mid, dummy)

                ! Get the trial epsilon
                    call HELM_EOSEPSILON(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid)

                ! Calculate the binding energy change
                    deps_nuc = binde_af - binde_bf

                ! Check if the energy is balanced
                    check_e = eps_mid - eps_beg - deps_nuc - &  !Original +
                        8.32696D-4 * ecaprate_mid * dt + eneurate_mid * dt

                ! Make sure you go to the right direction of dtemp
                    if(i == 0) then
                        if(check_e > 0.0D0) dtemp_mid = -dtemp_mid
                    endif

                    ! Use shooting method    
                    if(check_e_last * check_e < 0.0D0 .and. i /= 0) then
                        temp_mid = temp_mid - dtemp_mid
                        dtemp_mid = dtemp_mid * 0.1D0
                        temp_mid = temp_mid + dtemp_mid  
                        count_digit = count_digit + 1
                    else
                        temp_mid = temp_mid + dtemp_mid
                        check_e_last = check_e
                    endif

                    if(count_digit == 4) exit
                    if(i == 400 .or. temp_mid < 5.0D0) then
                        !stop 'Check NSE solver'
                        temp_mid = temp_beg        
                        eps_mid = eps_beg                  
                        !call nse_interface(temp_mid ,rho_mid, ye_mid, x_burn)
                        call getnsestate(rho_mid, temp_mid, x_burn)
                        exit
                    endif            

                    !if(j == 1 .and. k == 1) write(*,100) i, temp_mid, check_e, ecaprate_mid, eneurate_mid, (x_burn(k2), k2 = 1, totalion)
                    !if(mod(i,100) == 0 .and. j == 1 .and. k == 1) read(*,*)

                enddo

            !nse_burntime = EXP(196.02D0 / temp_mid - 41.646D0) / 4.9282D-6
            temp_nse = -8.6803D0 + 1.8787D0 * LOG10(rho_mid * 6.171D17)
            nse_burntime = EXP(196.02D0 / temp_nse - 41.646D0) / 4.9282D-6

            if(dt > nse_burntime) then

            ! When things are done, move the refined
            ! trial results to the outpuit
                temp2(j,1,k) = temp_mid       
                epsilon(j,1,k) = eps_mid
                prim(iye2,j,1,k) = ye_mid                
                prim(ihe4:ini56,j,1,k) = x_burn(:) 
                ! burn_qdot(j,1,k) = burn_qdot(j,1,k) + binde_af - binde_bf

            else

            ! When things are partially done, use
            ! linear interpolation

                temp_mid = temp_beg + (temp_mid - temp_beg) * dt / nse_burntime
                x_burn(:) = x_mid(:) + (x_burn(:) - x_mid(:)) * dt / nse_burntime
                call compute_binde(x_burn, binde_af)
                call private_helmeos_azbar(x_burn, abar_mid, zbar_mid, dummy)	
                call HELM_EOSEPSILON(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid)	   

            ! Now patch the result
                temp2(j,1,k) = temp_mid
                epsilon(j,1,k) = eps_mid
                prim(iye2,j,1,k) = ye_mid    
                prim(ihe4:ini56,j,1,k) = x_burn(:)
                ! burn_qdot(j,1,k) = burn_qdot(j,1,k) + binde_af - binde_bf

            endif

        else !(burn_ratio =/ 1)       

            ! For detonation, no NSE energy is accounted
                temp_mid = temp2(j,1,k)

            ! Electron capture for deflagration zone
            !if(flame_loc_ratio(j,1,k) == 1.0D0) then
            !   call getecaprate(rho_mid, temp_mid, ye2(j,1,k), ecaprate_mid, eneurate_mid)
                !   ye2(j,1,k) = ye2(j,1,k) + ecaprate_mid * dt
            !endif

            ! Get the NSE state
                call getnsestate(rho_mid, temp_mid, x_burn) 

            ! Copy the result to output
            prim(ihe4:ini56,j,1,k) = x_burn(:)
            ! burn_qdot(j,1,k) = 0.0D0  
        endif

        Enuc = (epsilon(j,1,k)-eps_beg)*prim(irho,j,1,k)*vol(j,1,k)

    enddo         
enddo

! Copy the results to ghost cells
call BOUNDARY()        

100 format(I5, 20ES15.7)
101 format(10ES15.7)

end subroutine nse2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetNSEState(rho_in, temp_in, xiso_nse_out)
USE DEFINITiON
USE CUSTOM_DEF
IMPLICIT NONE

REAL*8 :: rho_in, temp_in
REAL*8 :: xiso_nse_out(1:totalion)

REAL*8 :: log10rho

INTEGER :: rho_grid, temp_grid, ye_grid
REAL*8 :: rho_dgrid, temp_dgrid, ye_dgrid

log10rho = LOG10(rho_in * 6.171D17)                

IF(log10rho > 7.0D0) THEN            ! Minimum = 5.3
    if(temp_in > 5.0D0) THEN          ! Minimum = 3.0 and maximum 11.0

        rho_grid = INT((log10rho - 7.0D0) / 0.1D0)
        temp_grid = INT((temp_in - 4.0D0) / 0.1D0)

        rho_dgrid = (log10rho - (7.0D0 + (DBLE(rho_grid) * 0.1D0))) / 0.1D0
        temp_dgrid = (temp_in - (4.0D0 + (DBLE(temp_grid) * 0.1D0))) / 0.1D0

    IF(rho_grid >= 30) THEN
        rho_grid = 30         
        rho_dgrid = 0.0D0
        ENDIF

        IF(temp_grid >= 70) THEN
        temp_grid = 70
        temp_dgrid = 0.0D0                    
        ENDIF

        xiso_nse_out(:) =  nsetable_xiso(rho_grid, temp_grid, :) + &
                        rho_dgrid * temp_dgrid * &
        (nsetable_xiso(rho_grid+1, temp_grid+1, :) - nsetable_xiso(rho_grid, temp_grid, :)) + & 
                        rho_dgrid * (1.0D0 - temp_dgrid) * &
        (nsetable_xiso(rho_grid+1, temp_grid, :) - nsetable_xiso(rho_grid, temp_grid, :)) + &
                        (1.0D0 - rho_dgrid) * temp_dgrid * &
        (nsetable_xiso(rho_grid, temp_grid+1, :) - nsetable_xiso(rho_grid, temp_grid, :)) 

    ENDIF
ENDIF

END SUBROUTINE GetNSEState