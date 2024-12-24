!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER files for the rotating magnetised white dwarf  !
! The cold fermi gas is used as the equation of state      !
! Normal Newtonian gravity is assumed                      !
! Electron fraction is assumed to be 0.5 always everywhere !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unit constants !

! Physical constants to be as one !
REAL*8, PARAMETER :: gconst = 6.67430D-8
REAL*8, PARAMETER :: clight = 2.99792458D10
REAL*8, PARAMETER :: solar = 1.98847D33

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
REAL*8, PARAMETER :: taucgs2code = (rhocgs2code*energycgs2code)        ! (masscgs2code*lencgs2code**2/tcgs2code**2) !

! Current conversion !
REAL*8, PARAMETER :: amp2code = (mu_0*1.0D5*masscgs2code*lencgs2code)**(0.5D0)/tcgs2code

! Magnetic field !
REAL*8, PARAMETER :: gauss2code = 1.0D-1*masscgs2code/amp2code/tcgs2code**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters !
REAL*8, PARAMETER :: rhomax = 5.0D10*rhocgs2code
REAL*8, PARAMETER :: atmosphere = rhomax*1.0D-7
REAL*8, PARAMETER :: atmospheric = 1.0D-7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for solving gravity !

! Solve the potential per how many steps
INTEGER, PARAMETER :: n_pot = 20

! maximum number of relaxation !
INTEGER, PARAMETER :: relax_max = 1000000

! Tolerance in relaxation of the potential			
REAL*8, PARAMETER :: tolerance = 1.0D-8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8, PARAMETER :: ye_ini = 0.5D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NucEOS !
! EOS file path can be assigned from param.dat !
CHARACTER(LEN=100), PARAMETER :: default_eos = '/project/MCChu/eos/STOS.h5'
! CHARACTER(LEN=100), PARAMETER :: default_eos = '/ds/cmyip/eos_h5/STOS.h5'
! available: B139.h5  B145.h5  B155.h5  B165.h5  B165Z.h5  DD2F1.2.h5  LNS.h5  LS180.h5  LS220.h5  NRAPR.h5  SFHo.h5  SFHx.h5  STOS.h5 <-- this is HShen !
