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
      complex*16, parameter :: ZI = (0.d0,1.d0)
!
!c  ...material constants
      double precision :: MU,EPSILON,SIGMA
!c
!c  ...frequency
      double precision :: OMEGA
!c
!c  ...impedance constant
      double precision :: GAMMA
!c
!c  ...PI
      double precision :: PI
!c
!c  ...additional parameters including those required by the system
      integer :: ORDER_APPROX
      integer :: NPX, NPY, NPZ, ICOMP_EXACT
      integer,parameter :: MY_NR_RHS=1
      integer :: ICHOOSE_DISP
      integer :: IEXACT_DISP
!c
!c  ...choose case for exact
      integer :: ISOL
!c.... paraview parameters
      integer, parameter ::   iParAdap=0;
      integer, parameter ::   iParAttr(6) = (/1,0,1,0,0,1/)
      integer, dimension(10) :: ISEL_PARAVIEW = (/0,0,0,0,0,1,1,1,1,1/)
!
end module problem


