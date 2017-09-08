!----------------------------------------------------------------------
!
!   module name        - matrices
!
!----------------------------------------------------------------------
!
!   latest revision    - May 2
!
!   purpose            - store element matrices for the first layer
!                        elements
!
!----------------------------------------------------------------------
!
module matrices
!
  use parametersDPG
  implicit none
#if C_MODE
#define VTYPE  complex*16
#else
#define VTYPE double precision
#endif
!
! max # of elements in the first layer
  integer, parameter :: MAXNRFL=64
!
! # of stored elements in the first layer
  integer :: NRFL
!!$OMP THREADPRIVATE (NRFL)
!
! xy coordinates of the first vertex node
  real*8, dimension(2,3,MAXNRFL) :: XYVERT
!!$OMP THREADPRIVATE (XYVERT)
!
!
#if C_MODE
      complex*16, allocatable, save :: ZFL_EE(:,:,:)
!!$OMP THREADPRIVATE(ZFL_EE)
      complex*16, allocatable, save :: ZFL_EQ(:,:,:)
!!$OMP THREADPRIVATE(ZFL_EQ)
      complex*16, allocatable, save :: ZFL_QQ(:,:,:)
!!$OMP THREADPRIVATE(ZFL_QQ)
#else
      real*8, allocatable, save     :: ZFL_EE(:,:,:)
!!$OMP THREADPRIVATE(ZFL_EE)
      real*8, allocatable, save     :: ZFL_EQ(:,:,:)
!!$OMP THREADPRIVATE(ZFL_EQ)
      real*8, allocatable, save     :: ZFL_QQ(:,:,:)
!!$OMP THREADPRIVATE(ZFL_QQ)
#endif
! order of elements
      integer, parameter :: MYP=7
!
! matrix dimensions
      integer, parameter :: MYE = 3*MYP*(MYP+1)**2*2
      integer, parameter :: MYQ = MYP**3*6
end module matrices



