!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER files for the rotating magnetised white dwarf  !
! The cold fermi gas is used as the equation of state      !
! Normal Newtonian gravity is assumed                      !
! Electron fraction is assumed to be 0.5 always everywhere !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Feature Flags !

!Flag for Rotation (set initial vphi to 0 if 0)!
INTEGER, PARAMETER :: rotate_flag = 1

!Flag for MHD (set initial magnetic field to 0 if 0)!
INTEGER, PARAMETER :: mhd_flag = 1

! Flag for gravity
INTEGER, PARAMETER :: gravity_flag = 1

! Flag for testing phi (for this to work gravity flag needs to be 1)!
INTEGER, PARAMETER :: phitest_flag = 0

! Chemical Composition and Transport (ISO) Flag
INTEGER, PARAMETER :: xisotran_flag = 1

! Helmholtz EOS Flag and its Checking Flag
INTEGER, PARAMETER :: helmeos_flag = 1
INTEGER, PARAMETER :: helmcheck_flag = 0
INTEGER :: count = 0

! Turbulence Flag (only for cylindrical coordinates (uniform r z grid); spherical coordinates are not implemented correctly, the stress tensor viscosity definition needs to be modified)

INTEGER, PARAMETER :: turb_flag = 1
INTEGER, PARAMETER :: turbcheck_flag = 0
INTEGER, PARAMETER :: turb_neighbour = 3 ! Number of cells that is not atmosphere for turbulence

! Flame Flags ! (the level set is written with uniform grid in mind)
INTEGER, PARAMETER :: levelset_flag = 0
INTEGER, PARAMETER :: flame_flag = 0 ! Flag to initialize flame
INTEGER, PARAMETER :: update_flag = 0 ! Flag to change the level set by reinitialization
INTEGER, PARAMETER :: burn_flag = 0

! Flag for allowing the interaction between level set 1 and energy input
! 1 = Allow energy input by detonation + finding detonation

INTEGER, PARAMETER :: deton_flag = 1

! Flag for allowing final burning input for level set 1 & 2
! 1 = Allow energy input by NSE evolution

INTEGER, PARAMETER :: convert_nse_flag = 1

! Flag for allowing 1st step burning input for level set 1 & 2
! 1 = Allow energy input by carbon burning

INTEGER, PARAMETER 	:: carburn_flag = 1

! Flag for allowing 2nd step burning input for level set 1 & 2
! 1 = Allow energy input by advanced burning

INTEGER, PARAMETER 	:: advburn_flag = 1

! Flag for thermal neutrino !
INTEGER, PARAMETER 	:: thermal_flag = 1

! Flag for neutrino spectrum 
INTEGER, PARAMETER :: nuspec_flag = 1

! Flag for output
INTEGER, PARAMETER :: say_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Unit constants !

! Physical constants to be as one !
REAL*8, PARAMETER :: gconst = 6.67430D-8
REAL*8, PARAMETER :: clight = 2.99792458D10
REAL*8, PARAMETER :: solar = 1.98847D33

! Here, mu_0 is not in cgs !
REAL*8, PARAMETER :: mu_0 = 1.25663706212D-6 ! unit in kg m s-2 A-2 or H/m !

! Solar Radius !
REAL*8, PARAMETER :: rsolar = 6.96342D10

! 1 GeV !
REAL*8, PARAMETER :: GeV2gram = 1.78266191D-24

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unit conversions !
! Conversion between units !
REAL*8, PARAMETER :: lencgs2code = (clight**2)/(solar*gconst)
REAL*8, PARAMETER :: masscgs2code = (1.0D0/solar)
REAL*8, PARAMETER :: tcgs2code = (clight**3)/(solar*gconst)

! Derived conversion !
REAL*8, PARAMETER :: rhocgs2code = (masscgs2code/lencgs2code**3)
REAL*8, PARAMETER :: energycgs2code = (1.0D0/clight**2)
REAL*8, PARAMETER :: taucgs2code = (rhocgs2code*energycgs2code)        ! (masscgs2code*lencgs2code**2/tcgs2code**2) !
REAL*8, PARAMETER :: h_bar = (1.054571817D-27)*(lencgs2code**2*masscgs2code/tcgs2code)

! Current conversion !
REAL*8, PARAMETER :: amp2code = (mu_0*1.0D5*masscgs2code*lencgs2code)**(0.5D0)/tcgs2code

! Magnetic field !
REAL*8, PARAMETER :: gauss2code = 1.0D-1*masscgs2code/amp2code/tcgs2code**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8, PARAMETER :: me2 = 9.1093837015D-28*masscgs2code
REAL*8, PARAMETER :: mb2 = 1.66053906660D-24*masscgs2code
REAL*8, PARAMETER :: ye = 0.5D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters !
REAL*8, PARAMETER :: rhomax = 1.0D9*rhocgs2code
REAL*8 :: atmosphere 
REAL*8, PARAMETER :: atmospheric = 1.0D-7

! Constant for fermi equation of state !
! Note that the speed of light is unity !
REAL*8, PARAMETER :: amax = (me2**4)/(2.4D1*pi**2*h_bar**3)
REAL*8, PARAMETER :: bmax = (mb2*me2**3)/(3.0D0*pi**2*h_bar**3*ye)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for solving gravity !

! Solve the potential per how many steps
INTEGER, PARAMETER :: n_pot = 20

! maximum number of relaxation !
INTEGER, PARAMETER :: relax_max = 100000

! Tolerance in relaxation of the potential			
REAL*8, PARAMETER :: tolerance = 1.0D-8

! Number of multipoles at the boundary, maximum is quadrupole !
INTEGER, PARAMETER :: n_pole = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Helmholtz EOS Parameters !

! Maximum and Minimum Temperature Allowed, 1 unit = 1 billion kelvin !
REAL*8, PARAMETER	:: temp_max = 7.0D1		
REAL*8, PARAMETER 	:: temp_min = 1.0D-4

! Atmosphere Parameters !
REAL*8, PARAMETER    :: xiso_ahe4 = 0.0D0
REAL*8, PARAMETER    :: xiso_ac12 = 0.49D0
REAL*8, PARAMETER    :: xiso_ao16 = 0.51D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Analysis flags !
INTEGER, PARAMETER :: energy_analysis_flag = 0
