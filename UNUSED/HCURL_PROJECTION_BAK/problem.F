c----------------------------------------------------------------------
c
c   module name        - problem
c
c----------------------------------------------------------------------
c
c   latest revision    - May 14
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
      complex*16, parameter :: ZI = (0.d0,1.d0)
c
c  ...control flag for the test inner product
      integer :: INNER_PRODUCT
c
c  ...material constants
      double precision :: MU,EPSILON,SIGMA
c
c  ...frequency
      double precision :: OMEGA
c
c
      end module problem


