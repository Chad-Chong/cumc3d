!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute the poisson equation coefficient for the relaxation method 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_poisson
USE CUSTOM_DEF 
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Get poisson coefficient !
Do i = 1, nx
	Do j = 1, ny
		Do k = 1, nz
			CALL poisson_coef(x(i-1), x(i), x(i+1), &
							y(j-1), y(j), y(j+1), &
							z(k-1), z(k), z(k+1), &
							ajp1(i), ajm1(i), bkp1(i,j), bkm1(i,j), &
							clp1(i,j,k), clm1(i,j,k), epsc(i,j,k))
		end do
	end do
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get left hand side of the discrete poisson equation 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE poisson_coef(xm1, xc, xp1, ym1, yc, yp1, zm1, zc, zp1, &
  alphajp1, alphajm1, betakp1, betakm1, gammalp1, gammalm1, epsc)
USE CUSTOM_DEF 
USE DEFINITION
implicit none

! input !
real*8, INTENT(IN) :: xm1, xc, xp1, ym1, yc, yp1, zm1, zc, zp1

! dummy variables for calculation !
real*8 :: dxp1, dxm1, dyp1, dym1, dzp1, dzm1, fxp1, fxm1, fyp1, fym1, fzp1, fzm1 

! output !
real*8, INTENT(OUT) :: epsc
real*8, INTENT(OUT) :: alphajp1, alphajm1
real*8, INTENT(OUT) :: betakp1, betakm1
real*8, INTENT(OUT) :: gammalp1, gammalm1

IF (coordinate_flag == 2) THEN
	! assign !
	epsc = 2.0d0*(xp1+xm1-3.0d0*xc)/xc/(xp1-xc)/(xc-xm1) &
	+ (yp1+ym1-2.0d0*yc-2.0d0*DTAN(yc))/xc**2/DTAN(yc)/(yp1-yc)/(yc-ym1) &
	- 2.0d0/xc**2/DSIN(yc)**2/(zp1-zc)/(zc-zm1)

	! assign !
	alphajp1 = 2.0d0*(2.0d0*xc-xm1)/(xp1-xm1)/(xp1-xc)/xc
	alphajm1 = 2.0d0*(2.0d0*xc-xp1)/(xp1-xm1)/(xc-xm1)/xc

	! assign !
	betakp1 = (2.0d0*DTAN(yc) + yc - ym1)/xc**2/DTAN(yc)/(yp1 - yc)/(yp1 - ym1)
	betakm1 = (2.0d0*DTAN(yc) + yc - yp1)/xc**2/DTAN(yc)/(yp1 - ym1)/(yc - ym1)

	! assign !
	gammalp1 = 2.0d0/xc**2/DSIN(yc)**2/(zp1-zm1)/(zp1-zc)
	gammalm1 = 2.0d0/xc**2/DSIN(yc)**2/(zp1-zm1)/(zc-zm1)

ELSEIF (coordinate_flag == 1) THEN

	dxp1 = xp1 - xc
	dxm1 = xm1 - xc
	dyp1 = yp1 - yc
	dym1 = ym1 - yc
	dzp1 = zp1 - zc
	dzm1 = zm1 - zc

	fxp1 = 1/(dxp1*(dxm1-dxp1))
	fxm1 = 1/(dxm1*(dxm1-dxp1))
	fyp1 = 1/(dyp1*(dym1-dyp1))
	fym1 = 1/(dym1*(dym1-dyp1))
	fzp1 = 1/(dzp1*(dzm1-dzp1))
	fzm1 = 1/(dzm1*(dzm1-dzp1))

	epsc = 2.0D0 *(fxp1-fxm1 + (fyp1-fym1)/xc**2 + fzp1-fzm1) + 1.0D0/xc*(dxp1*fxm1 - dxm1*fxp1)

	alphajp1 = 1.0D0/xc*dxm1*fxp1 - 2.0D0*fxp1
	alphajm1 = -1.0D0/xc*dxp1*fxm1 + 2.0D0*fxm1

	betakp1 = -2.0D0*fyp1/xc**2
	betakm1 = 2.0D0*fym1/xc**2

	gammalp1 = -2.0D0*fzp1
	gammalm1 = 2.0D0*fzm1

ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute the multipole expansion
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE multipole_expansion(mono, dipo, quad)
USE CUSTOM_DEF 
USE DEFINITION
IMPLICIT NONE

real*8 :: diff, factor, rho_in, carte_x, carte_y, carte_z
integer :: j,k,l
real*8, intent(out) :: mono
real*8, dimension(3) :: posit
real*8, dimension(3), intent(out) :: dipo
real*8, dimension(3,3) :: qposit
real*8, dimension(3,3), intent(out) :: quad

mono = 0.0D0
dipo(:) = 0.0D0
quad(:,:) = 0.0D0

IF (n_pole == 0) THEN
	DO j = 1, nx
		DO k = 1, ny
			DO l = 1, nz
				diff = prim(irho,j,k,l) - prim_a(irho)
				factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
				rho_in = factor*prim(irho,j,k,l)
				mono = mono + rho_in*vol(j,k,l)
			ENDDO
		ENDDO
	ENDDO

ELSEIF (n_pole == 1) THEN

DO j = 1, nx
	DO k = 1, ny
		DO l = 1, nz
			diff = prim(irho,j,k,l) - prim_a(irho)
			factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
			rho_in = factor*prim(irho,j,k,l)

			mono = mono + rho_in*vol(j,k,l)

			posit(1) = x(j)*DCOS(y(k))
			posit(2) = x(j)*DSIN(y(k))
			posit(3) = z(l)
			dipo(:) = dipo(:) + posit(:)*rho_in*vol(j,k,l)

		ENDDO
	ENDDO
ENDDO

ELSEIF (n_pole == 2) THEN

DO j = 1, nx
	DO k = 1, ny
		DO l = 1, nz
			diff = prim(irho,j,k,l) - prim_a(irho)
			factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
			rho_in = factor*prim(irho,j,k,l)
			mono = mono + rho_in*vol(j,k,l)

			posit(1) = x(j)*DCOS(y(k))
			posit(2) = x(j)*DSIN(y(k))
			posit(3) = z(l)
			dipo(:) = dipo(:) + posit(:)*rho_in*vol(j,k,l)

			qposit(1,1) = posit(1)*posit(1)
			qposit(1,2) = posit(1)*posit(2)
			qposit(1,3) = posit(1)*posit(3)
			qposit(2,1) = posit(2)*posit(1)
			qposit(2,2) = posit(2)*posit(2)
			qposit(2,3) = posit(2)*posit(3)
			qposit(3,1) = posit(3)*posit(1)
			qposit(3,2) = posit(3)*posit(2)
			qposit(3,3) = posit(3)*posit(3)
			quad(:,:) = quad(:,:) + (3.0D0*qposit(:,:)-(x(j)**2+z(l)**2)*eye(:,:))*rho_in*vol(j,k,l)

		ENDDO
	ENDDO
ENDDO

ELSE

	WRITE(*,*), 'n_pole not within range'
	STOP

ENDIF
END SUBROUTINE