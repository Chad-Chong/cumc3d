!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER files for 2-fluid star !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! mode !
INTEGER, PARAMETER :: iNM = 1
INTEGER, PARAMETER :: iDM = 2

! Unit constants !

! Physical constants to be as one !
REAL*8, PARAMETER :: gconst = 6.67430D-8
REAL*8, PARAMETER :: clight = 2.99792458D10
REAL*8, PARAMETER :: solar = 1.98847D33
REAL*8, PARAMETER :: GeV2gram = 1.78266191D-24

! Here, mu_0 is not in cgs !
REAL*8, PARAMETER :: mu_0 = 1.25663706212D-6 ! unit in kg m s-2 A-2 or H/m !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unit conversions !
! Conversion between units !
REAL*8, PARAMETER :: lencgs2code = (clight**2)/(solar*gconst)
REAL*8, PARAMETER :: masscgs2code = (1.0D0/solar)
REAL*8, PARAMETER :: tcgs2code = (clight**3)/(solar*gconst)

! Derived conversion !
REAL*8, PARAMETER :: rhocgs2code = (masscgs2code/lencgs2code**3)
REAL*8, PARAMETER :: energycgs2code = (1.0D0/clight**2)
REAL*8, PARAMETER :: taucgs2code = (rhocgs2code*energycgs2code)

REAL*8, PARAMETER :: amp2code = (mu_0*1.0D5*masscgs2code*lencgs2code)**(0.5D0)/tcgs2code
REAL*8, PARAMETER :: gauss2code = 1.0D-1*masscgs2code/amp2code/tcgs2code**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8, PARAMETER :: me2 = 9.1093837015D-28*masscgs2code
REAL*8, PARAMETER :: mb2 = 1.66053906660D-24*masscgs2code
REAL*8, PARAMETER :: mx = 1.0D-1*GeV2gram*masscgs2code

REAL*8, PARAMETER :: ye = 0.5D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters !
REAL*8, PARAMETER :: rhomax = 1.0D6*rhocgs2code
REAL*8, PARAMETER :: rhoatm = 1.0D3*rhocgs2code
REAL*8, PARAMETER :: rhomin = 1.0D3*rhocgs2code
REAL*8, PARAMETER :: h_bar = (1.054571817D-27)*(lencgs2code**2*masscgs2code/tcgs2code)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constant for fermi equation of state  !
! Note that the speed of light is unity !
REAL*8, PARAMETER :: amaxNM = (me2**4)/(2.4D1*pi**2*h_bar**3)
REAL*8, PARAMETER :: bmaxNM = (mb2*me2**3)/(3.0D0*pi**2*h_bar**3*ye)

REAL*8, PARAMETER :: amaxDM = (mx**4)/(2.4D1*pi**2*h_bar**3)
REAL*8, PARAMETER :: bmaxDM = (mx**4)/(3.0D0*pi**2*h_bar**3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for solving gravity !

! Solve the potential per how many steps
INTEGER, PARAMETER :: n_pot = 20

! maximum number of relaxation !
INTEGER, PARAMETER :: relax_max = 1000000

! Tolerance in relaxation of the potential			
REAL*8, PARAMETER :: tolerance = 1.0D-8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Physics flags !
! GR potential !
! We use the modified case A potential (Muller et al. 2008)!
INTEGER, PARAMETER :: gr_potential = 0

! Lapse function !
INTEGER, PARAMETER :: lapse_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
