!----------------------------------------------------------------------
!
!   module name        - problem
!
!----------------------------------------------------------------------
!
!   latest revision    - April 17
!
!   purpose            - problem dependent data
!
!----------------------------------------------------------------------
!
module problem
!
  use parametersDPG
  implicit none
#if C_MODE
#define V_TYPE  complex*16
#else
#define V_TYPE double precision
#endif
!
!  ...control flag defining a problem and geometry
  integer :: NO_PROBLEM, GEOM_NO
!
!  ...control flag for the test inner product
  integer :: INNER_PRODUCT
!  ...number of steps
  integer :: NSTEPS
!
!  ...toggle laser mode, nonlinear problem
  integer :: LASER_MODE, NONLINEAR_FLAG
!
!  ...kappa, deltaT, Nsteps, HEAT_LOAD for the heat equation
  double precision :: KAPPA, DELTAT, TMAX, HEAT_LOAD

! ...Maxwell equation parameters
  double precision :: OMEGA, BEAM_WAIST, GAMMA_IMP
!  ...material constants
  double precision :: MU, NTHETA
  double precision :: EPSILON
  double precision :: SIGMA
!  PI and sqrt(-1)
  complex*16, parameter :: ZI = (0.d0,1.d0)
  REAL*8, parameter :: PI             = 3.14159265358979D0
!
!............NON-DIMENSIONALIZATION....................
!  ...Parameters for the gain function
  double precision, parameter :: REF_INDEX_CORE = 1.45001d0
  double precision, parameter :: REF_INDEX_CLAD = 1.45d0
  REAL*8 TAU
  PARAMETER (TAU = 8.014D-4)
  REAL*8 MU0
  PARAMETER (MU0 = 1.2566370614D-6)
  REAL*8 LIGHT_SPEED
  PARAMETER (LIGHT_SPEED = 2.99792458D8)
  REAL*8 CORE_RAD
  PARAMETER (CORE_RAD = 12.32155D-6)
  REAL*8 N0
  PARAMETER (N0 = 8.75D29)
  REAL*8 H_BAR
  PARAMETER (H_BAR = 1.05457266D-34)
  REAL*8, parameter :: OMEGA0 = 1.7704D15
  double precision, parameter :: OMEGA_RATIO = 1064.d0/976.d0
  double precision, parameter :: INTENSITY_RATIO = sqrt(25.d0/ &
                                             (1000.d0))
  double precision, parameter :: I_INC = 1000.d0/(3.14159*CORE_RAD**2)
  double precision, parameter :: E_INC = sqrt(I_INC*MU0*LIGHT_SPEED &
                                                  /REF_INDEX_CORE)
  double precision, parameter :: E_0 = INTENSITY_RATIO*E_INC
  REAL*8, parameter :: SIGMA_S_ABS    = 6.0D-27
  REAL*8, parameter :: SIGMA_S_EM     = 3.58D-25
  REAL*8, parameter :: SIGMA_P_ABS    = 1.429D-24
  REAL*8, parameter :: SIGMA_P_EM     = 1.776D-24

!  ...IBCFlag for the impedance BC
  integer     :: IBCFlag
!
! parameters for polynomial solution, exact solution number, components etc.
  integer :: ORDER_APPROX
  integer :: NPX, NPY, NPZ,ICOMP_EXACT
  integer :: ISOL
  integer,parameter :: MY_NR_RHS=1
  integer :: VARIABLE, COMPONENT

!  ...graphics options
  integer :: IEXACT_DISP
  integer :: ICHOOSE_DISP
  integer :: ICHOOSE_COMP
!
!.... paraview parameters
  integer, parameter ::   iParAdap=0;
  integer, parameter ::   iParAttr(6) = (/1,0,1,0,0,  1/)
  integer, dimension(10) :: ISEL_PARAVIEW = (/0,0,0,0,0,1,1,1,1,1/)
!
end module problem


