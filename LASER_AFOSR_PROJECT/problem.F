c----------------------------------------------------------------------
c
c   module name        - problem
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 17
c
c   purpose            - problem dependent data
c
c----------------------------------------------------------------------
c
      module problem
c
      use parametersDPG
#if C_MODE
#define V_TYPE  complex*16
#else
#define V_TYPE double precision
#endif
c
c  ...control flag defining a problem and geometry
      integer :: NO_PROBLEM, GEOM_NO
c
c  ...control flag for the test inner product
      integer :: INNER_PRODUCT
c
c  ...number of steps
      integer :: NSTEPS
c
c  ...toggle laser mode
      integer :: LASER_MODE, NONLINEAR_FLAG

c  ...kappa, deltaT, Nsteps, HEAT_LOAD for the heat equation
      double precision :: KAPPA, DELTAT, TMAX, HEAT_LOAD

c  ...Maxwell equation parameters
      double precision :: OMEGA, BEAM_WAIST
c
c  ...Parameters for the gain function
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
      double precision, parameter :: INTENSITY_RATIO = sqrt(25.d0/
     .                                           (1000.d0))
      double precision, parameter :: I_INC = 1000/(3.14159*CORE_RAD**2)
      double precision, parameter :: E_INC = sqrt(I_INC*MU0*LIGHT_SPEED
     .                                             /REF_INDEX_CORE)
      double precision, parameter :: E_0 = INTENSITY_RATIO*E_INC
      REAL*8, parameter :: SIGMA_S_ABS    = 6.0D-27
      REAL*8, parameter :: SIGMA_S_EM     = 3.58D-25
      REAL*8, parameter :: SIGMA_P_ABS    = 1.429D-24
      REAL*8, parameter :: SIGMA_P_EM     = 1.776D-24
      REAL*8, parameter :: PI             = 3.14159265358979D0


c  ...IBCFlag for the impedance BC
      integer     :: IBCFlag

c  ...material constants
      double precision :: MU, NTHETA
      double precision :: EPSILON
      double precision :: SIGMA
      complex*16, parameter :: ZI = (0.d0,1.d0)
c
c.... paraview parameters
      integer, parameter ::   iParAdap=0;
      integer, parameter ::   iParAttr(6) = (/1,0,1,0,0,  1/)
      integer, dimension(10) :: ISEL_PARAVIEW = (/0,0,0,0,0,1,1,1,1,1/)
c
      end module problem


