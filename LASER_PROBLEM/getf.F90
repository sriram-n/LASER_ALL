!------------------------------------------------------------------------------
!> Purpose : source term
!!
!! @param[in]  Mdle  - element (middle node) number
!! @param[in]  X     - physical coordinates
!! @param[out] Zfval - rhs
!------------------------------------------------------------------------------
!
!
subroutine getf(Mdle,X, Zfval,zJval)
!
      use control          , only : NEXACT
      use parameters       , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
      use problem
!
#if C_MODE
#define VTYPE  complex*16
#else
#define VTYPE double precision
#endif
      implicit none
      integer,                  intent(in)  :: Mdle
      real*8,dimension(3),      intent(in)  :: X
      VTYPE,                    intent(out) :: Zfval
      VTYPE,dimension(3),       intent(out) :: ZJval
      VTYPE                                 :: zaux
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
      Zfval = ZERO
      ZJval = ZERO
!
      select case(NEXACT)
!==============================================================================
!  UNKNOWN EXACT SOLUTION                                                      |
!==============================================================================
      case(0)
!
        select case(NO_PROBLEM)
!
!  .....single time step of heat
        case(1)
        Zfval = 1.d0
        ZJval(1:3) = ZERO
!
!  .....transient heat equation
        case(2)
        Zfval = 1.d0
        ZJval(1:3) = ZERO

!  .....time harmonic Maxwell
        case(3)
        Zfval = 1.d0
        ZJval(1:3) = ZERO
        end select

!
!==============================================================================
!  KNOWN EXACT SOLUTION, NON-ZERO RHS
!==============================================================================
      case(1)
!
!       compute exact solution
        call exact(X,Mdle, valH,dvalH,d2valH, valE,dvalE,d2valE, &
                           valV,dvalV,d2valV, valQ,dvalQ,d2valQ)
!
        select case(NO_PROBLEM)
!
!  .....single time step
        case(1)
          Zfval = -KAPPA*DELTAT*(D2valH(1,1,1) + D2valH(1,2,2) + D2valH(1,3,3)) + valH(1)
!
!  .....transient heat equation
        case(2)
          write(*,*) 'getf: INCONSISTENCY'; stop 1

!  .....time harmonic Maxwell
        case(3)
        zaux = ZI*OMEGA*EPSILON + SIGMA
        zJval(1) = DvalE(3,2,2) - DvalE(2,2,3) - zaux*ValE(1,1)
        zJval(2) = DvalE(1,2,3) - DvalE(3,2,1) - zaux*ValE(2,1)
        zJval(3) = DvalE(2,2,1) - DvalE(1,2,2) - zaux*ValE(3,1)
        end select
!
!
!==============================================================================
!  KNOWN EXACT SOLUTION , HOMOGENEOUS RHS                                     |
!==============================================================================
      case(2)
!
      endselect
!
      if (iprint.eq.1) then
        write(*,7010) Zfval
 7010   format(' getf: Zfval = ',2e12.5)
        call pause
      endif
!
!
endsubroutine getf



subroutine get_bdSource(Mdle,X,Rn, Imp_val)
!
      use control          , only : NEXACT
      use assembly         , only : NR_RHS
      use data_structure3D , only : NR_COMP,ADRES,NRINDEX
      use parameters       , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
      use problem
!
      implicit none
      integer,                  intent(in)  :: Mdle
      real*8,dimension(3),      intent(in)  :: X
      real*8,dimension(3),      intent(in)  :: Rn
      !integer,                  intent(in)  :: IBCFlag
      VTYPE,dimension(3),       intent(out) :: Imp_val
!------------------------------------------------------------------------------
!
!     exact solution
      VTYPE,dimension(  MAXEQNH    ) ::   zvalH
      VTYPE,dimension(  MAXEQNH,3  ) ::  zdvalH
      VTYPE,dimension(  MAXEQNH,3,3) :: zd2valH
      VTYPE,dimension(3,MAXEQNE    ) ::   zvalE
      VTYPE,dimension(3,MAXEQNE,3  ) ::  zdvalE
      VTYPE,dimension(3,MAXEQNE,3,3) :: zd2valE
      VTYPE,dimension(3,MAXEQNV    ) ::   zvalV
      VTYPE,dimension(3,MAXEQNV,3  ) ::  zdvalV
      VTYPE,dimension(3,MAXEQNV,3,3) :: zd2valV
      VTYPE,dimension(  MAXEQNQ    ) ::   zvalQ
      VTYPE,dimension(  MAXEQNQ,3  ) ::  zdvalQ
      VTYPE,dimension(  MAXEQNQ,3,3) :: zd2valQ
!
!     miscellaneus
      integer :: iload,ivar,ibeg,icomp,jcomp,k,l
      complex*16 :: zaux
      VTYPE,dimension(3) ::   rntimesE,rn2timesE
      VTYPE,dimension(3) ::   rntimesH
      real*8                    :: impedanceConstant
      real*8                    :: E   ! vector field
      real*8, dimension(3)      :: dE  ! 1st derivative
      real*8, dimension(3,3)    :: d2E ! 2nd derivative
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
 7001   format(' get_bdSource: Mdle,X = ',i8,2x,3(f8.3,2x))
      endif
!
!     initialize source terms
      Imp_val = ZERO
      if(GEOM_NO.eq.1) then
        !impedanceConstant = sqrt(1.d0-(PI**2/OMEGA**2))
        impedanceConstant = GAMMA_IMP
      else
        impedanceConstant = GAMMA_IMP
      endif
!
      select case(NEXACT)
!  ... known HOMOGENEOUS solution: do nothing
      case(0,2)
!
!  ...exact solution known manufactured solution
      case(1)
        select case(IBCFlag)
!
!  .....impedance BC
        case(3)
          call exact(X,Mdle, zvalH,zdvalH,zd2valH, &
                        zvalE,zdvalE,zd2valE, &
                        zvalV,zdvalV,zd2valV, &
                        zvalQ,zdvalQ,zd2valQ)

          call my_cross_product(Rn,zvalE(1:3,1), rntimesE)
          call my_cross_product(Rn,rntimesE, rn2timesE)
          call my_cross_product(Rn,zvalE(1:3,2), rntimesH)

!
          Imp_val = rntimesH - ((impedanceConstant))*rn2timesE
          !write(*,*) 'Imp_val is = ', Imp_val
!
        case default
          !write(*,*) 'get_bdSource: IBCFlag = ',IBCFlag;stop
        end select
!
!  ...exact solution unknown
      case default
        write(*,*) 'get_bdSource: UNSPECIFIED NEXACT';stop
      end select
      if (iprint.eq.1) then
        write(*,7002) Imp_val
 7002   format('get_bsource: Imp_val = ',2e12.5)
      endif
!
end subroutine get_bdSource
!
!
subroutine my_cross_product(A,Zb,Zc)
!
      use problem
!
      implicit none
      real*8,dimension(3),      intent(in)  :: A
      VTYPE,dimension(3),       intent(in)  :: Zb
      VTYPE,dimension(3),       intent(out) :: Zc
!
      Zc(1) =   A(2)*Zb(3) - A(3)*Zb(2)
      Zc(2) = - A(1)*Zb(3) + A(3)*Zb(1)
      Zc(3) =   A(1)*Zb(2) - A(2)*Zb(1)
!
end subroutine my_cross_product


subroutine zcross_product(A, Zb, Zc)
#include "syscom.blk"
      dimension A(3),Zb(3), Zc(3)
!
      Zc(1) =   A(2)*Zb(3) - A(3)*Zb(2)
      Zc(2) = - A(1)*Zb(3) + A(3)*Zb(1)
      Zc(3) =   A(1)*Zb(2) - A(2)*Zb(1)
!
end subroutine zcross_product
