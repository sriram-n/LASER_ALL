!------------------------------------------------------------------------------
!> Purpose : source term
!!
!! @param[in]  Mdle  - element (middle node) number
!! @param[in]  X     - physical coordinates
!! @param[out] Fval  - rhs
!------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine getf(Mdle,X, Jval)
!
      use control          , only : NEXACT
      use assembly         , only : NR_RHS 
      use data_structure3D , only : NR_COMP,ADRES,NRINDEX
      use parameters       , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
      use problem          , only : ZI,MU,EPSILON,SIGMA,OMEGA
!
      implicit none
      integer,                  intent(in)  :: Mdle
      real*8,dimension(3),      intent(in)  :: X
      VTYPE,dimension(3),       intent(out) :: Jval
!------------------------------------------------------------------------------
!      
!     exact solution
      VTYPE,dimension(  MAXEQNH    ) ::   valH
      VTYPE,dimension(  MAXEQNH,3  ) ::  dvalH
      VTYPE,dimension(  MAXEQNH,3,3) :: d2valH
      VTYPE,dimension(3,MAXEQNE    ) ::   valE
      VTYPE,dimension(3,MAXEQNE,3  ) ::  dvalE
      VTYPE,dimension(3,MAXEQNE,3,3) :: d2valE
      VTYPE,dimension(3,MAXEQNV    ) ::   valV
      VTYPE,dimension(3,MAXEQNV,3  ) ::  dvalV
      VTYPE,dimension(3,MAXEQNV,3,3) :: d2valV
      VTYPE,dimension(  MAXEQNQ    ) ::   valQ
      VTYPE,dimension(  MAXEQNQ,3  ) ::  dvalQ
      VTYPE,dimension(  MAXEQNQ,3,3) :: d2valQ
!
!     miscellaneus
      integer :: iload,ivar,ibeg,icomp,jcomp,k,l
      complex*16 :: zk2
!
!     printing flag
      integer :: iprint
!------------------------------------------------------------------------------
!
!     printing flag : 0 - silent ; 1 - verbose
      iprint=0
!      
      if (iprint.eq.1) then
        write(*,7001) Mdle,X
 7001   format(' getf: Mdle,X = ',i8,2x,3(f8.3,2x))
      endif
!      
!     initialize source terms
      Jval = ZERO; 
!
      select case(NEXACT)
!==============================================================================
!  UNKOWN EXACT SOLUTION                                                      |
!==============================================================================
      case(0)  
!!!        Jval(1) = 1.d0
! 
!==============================================================================
!  KNOWN EXACT SOLUTION, NON-ZERO RHS                                                    
!==============================================================================
      case(1)
!              
!     compute exact solution      
      call exact(X,Mdle, valH,dvalH,d2valH, valE,dvalE,d2valE, &
                         valV,dvalV,d2valV, valQ,dvalQ,d2valQ)
!
      Jval(1) = valE(1,1)
      Jval(2) = valE(2,1)
      Jval(3) = valE(3,1)
  
      

!
!==============================================================================
!  KNOWN EXACT SOLUTION , HOMOGENEOUS RHS                                     |
!==============================================================================
      case(2)
!              
      endselect        
!
      if (iprint.eq.1) then
        write(*,7010) Jval
 7010   format(' getf: Jval = ',e12.5)
      endif
!
!
endsubroutine getf
