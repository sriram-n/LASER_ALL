c----------------------------------------------------------------------
c
c   module name        - DPGLaser
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 17
c
c   purpose            - work space for the primal DPG method related
c                        computations of optimal test functions and
c                        energy error for Laser problem
c
c----------------------------------------------------------------------
c
      module DPGLaser
c
      use parametersDPG
#if C_MODE
#define V_TYPE  complex*16
#else
#define V_TYPE double precision
#endif
c

c ........Heat equation variable allocation

c  ...H1 discontinuous shape functions
      double precision shapHH(MAXbrickHH),gradHH(3,MAXbrickHH)
c$OMP THREADPRIVATE (shapHH,gradHH)
c  ...Riesz matrix size
      integer, parameter :: MAXtestH = MAXbrickHH

c  ...stiffnes matrix for the local Riesz H1 matrix in Lapack format
      V_TYPE AP_Heat(MAXtestH*(MAXtestH+1)/2)
c$OMP THREADPRIVATE (AP_Heat)
c
c  ...load vector for the enriched space
      V_TYPE BLOADH(MAXtestH)
c$OMP THREADPRIVATE (BLOADH)
c
c  ...copy of load vector for the enriched space
      V_TYPE BLOADHc(MAXtestH)
c$OMP THREADPRIVATE (BLOADHc)
c
c  ...stiffnes matrices for the enriched test space
      V_TYPE STIFFHH(MAXtestH,MAXbrickH)
c$OMP THREADPRIVATE(STIFFHH)
      V_TYPE STIFFHV(MAXtestH,MAXbrickV)
c$OMP THREADPRIVATE(STIFFHV)
c
c  ...copies of the same matrices
      V_TYPE STIFFHHc(MAXtestH,MAXbrickH)
c$OMP THREADPRIVATE(STIFFHHc)
      V_TYPE STIFFHVc(MAXtestH,MAXbrickV)
c$OMP THREADPRIVATE(STIFFHVc)
c
c  ....STIFF_ALL for alternative computation of stiffness
      V_TYPE STIFF_ALLH(MAXtestH,MAXbrickH+MAXbrickV+1)
c$OMP THREADPRIVATE(STIFF_ALLH)
c ........Maxwell equation variable allocation
c  ...Hcurl discontinuous shape functions
      double precision shapEE(3,MAXbrickEE),curlEE(3,MAXbrickEE)
c$OMP THREADPRIVATE(shapEE,curlEE)
c
c  ...Riesz matrix size
      integer, parameter :: MAXtestE = 2*MAXbrickEE
c
c  ...stiffnes matrix for the local Riesz Hcurl matrix in Lapack format
      V_TYPE AP_Maxwell(MAXtestE*(MAXtestE+1)/2)
c$OMP THREADPRIVATE(AP_Maxwell)
c  ... for diagonal preconditioning
#if C_MODE
      complex*16, allocatable, save :: AP_eig(:)
c$OMP THREADPRIVATE(AP_eig)
      complex*16, allocatable, save :: DIAG_E(:)
c$OMP THREADPRIVATE(DIAG_E)
      complex*16, allocatable, save :: DIAG_H(:)
c$OMP THREADPRIVATE(DIAG_H)
#else
      real*8, allocatable, save     :: AP_eig(:)
c$OMP THREADPRIVATE(AP_eig)
      real*8, allocatable, save     :: DIAG_E(:)
c$OMP THREADPRIVATE(DIAG_E)
      real*8, allocatable, save     :: DIAG_H(:)
c$OMP THREADPRIVATE(DIAG_H)
#endif
c
c  ...load vector for the enriched space
      V_TYPE BLOADE(MAXtestE)
c$OMP THREADPRIVATE(BLOADE)
c
c  ...copy of load vector for the enriched space
      V_TYPE BLOADEc(MAXtestE)
c$OMP THREADPRIVATE(BLOADEc)
c
c  ...stiffnes matrices for the enriched test space
      V_TYPE STIFFEQ(MAXtestE,MAXbrickQ*6)
c$OMP THREADPRIVATE(STIFFEQ)
      V_TYPE STIFFEE(MAXtestE,MAXbrickE*2)
c$OMP THREADPRIVATE(STIFFEE)
c
c  ...copies of the same matrices
      V_TYPE STIFFEQc(MAXtestE,MAXbrickQ*6)
c$OMP THREADPRIVATE(STIFFEQc)
      V_TYPE STIFFEEc(MAXtestE,MAXbrickE*2)
c$OMP THREADPRIVATE(STIFFEEc)

c  ....STIFF_ALL for alternative computation of stiffness
      V_TYPE STIFF_ALLE(MAXtestE,2*MAXbrickE+6*MAXbrickQ+1)
c$OMP THREADPRIVATE(STIFF_ALLE)

c  ...additional parameters including those required by the system
      integer :: ORDER_APPROX
      integer :: NPX, NPY, NPZ,ICOMP_EXACT
      integer :: ISOL
      integer,parameter :: MY_NR_RHS=1
      integer :: VARIABLE, COMPONENT

c  ...graphics options, there should be a better place for it...
      integer :: IEXACT_DISP  ! approximate/exact solution  (0/1)
      integer :: ICHOOSE_COMP ! unknown/flux                (1/2)

c ..... dummy for elem_residual
      V_TYPE AP(MAXtestH*(MAXtestH+1)/2)
c$OMP THREADPRIVATE(AP)
c
      end module DPGLaser


