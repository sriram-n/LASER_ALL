!------------------------------------------------------------------------------
!> Purpose : exact (manufactured, no miracle!) solution
!!
!> @param[in]  X      - a point in physical space
!> @param[in]  Mdle   - element (middle node) number
!> @param[out] ValH   - value of the H1 solution
!> @param[out] DvalH  - corresponding first derivatives
!> @param[out] D2valH - corresponding second derivatives
!> @param[out] ValE   - value of the H(curl) solution
!> @param[out] DvalE  - corresponding first derivatives
!> @param[out] D2valE - corresponding second derivatives
!> @param[out] ValV   - value of the H(div) solution
!> @param[out] DvalV  - corresponding first derivatives
!> @param[out] D2valV - corresponding second derivatives
!> @param[out] ValQ   - value of the H(div) solution
!> @param[out] DvalQ  - corresponding first derivatives
!> @param[out] D2valQ - corresponding second derivatives
!------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine exact(Xp,Mdle, ValH,DvalH,D2valH, ValE,DvalE,D2valE, &
                          ValV,DvalV,D2valV, ValQ,DvalQ,D2valQ)
!
  use data_structure3D
  use problem
!
  implicit none
  real*8,dimension(3),            intent(in)  :: Xp
  integer                       , intent(in)  :: Mdle
  VTYPE,dimension(  MAXEQNH    ), intent(out) ::   ValH
  VTYPE,dimension(  MAXEQNH,3  ), intent(out) ::  DvalH
  VTYPE,dimension(  MAXEQNH,3,3), intent(out) :: D2valH
  VTYPE,dimension(3,MAXEQNE    ), intent(out) ::   ValE
  VTYPE,dimension(3,MAXEQNE,3  ), intent(out) ::  DvalE
  VTYPE,dimension(3,MAXEQNE,3,3), intent(out) :: D2valE
  VTYPE,dimension(3,MAXEQNV    ), intent(out) ::   ValV
  VTYPE,dimension(3,MAXEQNV,3  ), intent(out) ::  DvalV
  VTYPE,dimension(3,MAXEQNV,3,3), intent(out) :: D2valV
  VTYPE,dimension(  MAXEQNQ    ), intent(out) ::   ValQ
  VTYPE,dimension(  MAXEQNQ,3  ), intent(out) ::  DvalQ
  VTYPE,dimension(  MAXEQNQ,3,3), intent(out) :: D2valQ
!
!------------------------------------------------------------------------------
!     Space for temporary solutions
!
  real*8                   :: u
  VTYPE                    :: E
  real*8, dimension(3)     :: gradu
  VTYPE,dimension(3)       :: dE
  real*8, dimension(3,3)   :: grad2u
  VTYPE, dimension(3,3)    :: d2E
  real*8                   :: ut
  integer                  :: icomp

!
!     initialize exact solution
  ValE=ZERO ; DvalE=ZERO ; D2valE=ZERO
  ValH=ZERO ; DvalH=ZERO ; D2valH=ZERO
  ValV=ZERO ; DvalV=ZERO ; D2valV=ZERO
  ValQ=ZERO ; DvalQ=ZERO ; D2valQ=ZERO
!     set icomp
  icomp = ICOMP_EXACT

!     initialize variables
  u = ZERO; E = ZERO; ut = ZERO
  gradu = ZERO; dE = ZERO
  grad2u = ZERO; d2E = ZERO
!
  select case (NO_PROBLEM)
  case(1,2)
    call h1_solution(Xp, u,gradu,grad2u,ut)
  case(3)
    call hcurl_solution(Xp, E,dE,d2E)
  endselect
!
  ValH(1) = u                               ! H1 variable
  DvalH(1,1:3) = gradu(1:3)                 ! 1st der
  D2valH(1,1:3,1:3) = grad2u(1:3,1:3)       ! 2nd der
!
  ValV(1:3,1) = gradu(1:3)                  ! Hdiv
  DvalV(1:3,1,1:3) = grad2u(1:3,1:3)        ! 1st der

 !  Efield.....value
  ValE(icomp,1    ) =   E
  !ValE(1:3,1    ) =   E
!
!  .....1st order derivatives
  DvalE(icomp,1,1  ) =  dE(1)
  DvalE(icomp,1,2  ) =  dE(2)
  DvalE(icomp,1,3  ) =  dE(3)
  ! DvalE(1:3,1,1  ) =  dE(1)
  ! DvalE(1:3,1,2  ) =  dE(2)
  ! DvalE(1:3,1,3  ) =  dE(3)
!
!  .....2nd order derivatives
  D2valE(icomp,1,1,1) = d2E(1,1)
  D2valE(icomp,1,1,2) = d2E(1,2)
  D2valE(icomp,1,1,3) = d2E(1,3)
  D2valE(icomp,1,2,1) = d2E(2,1)
  D2valE(icomp,1,2,2) = d2E(2,2)
  D2valE(icomp,1,2,3) = d2E(2,3)
  D2valE(icomp,1,3,1) = d2E(3,1)
  D2valE(icomp,1,3,2) = d2E(3,2)
  D2valE(icomp,1,3,3) = d2E(3,3)
  ! D2valE(1:3,1,1,1) = d2E(1,1)
  ! D2valE(1:3,1,1,2) = d2E(1,2)
  ! D2valE(1:3,1,1,3) = d2E(1,3)
  ! D2valE(1:3,1,2,1) = d2E(2,1)
  ! D2valE(1:3,1,2,2) = d2E(2,2)
  ! D2valE(1:3,1,2,3) = d2E(2,3)
  ! D2valE(1:3,1,3,1) = d2E(3,1)
  ! D2valE(1:3,1,3,2) = d2E(3,2)
  ! D2valE(1:3,1,3,3) = d2E(3,3)
!
!     2nd H(curl) ATTRIBUTE = curl of the first attribute/-i omega \mu
!
!  ...value
  ValE(1,2) = DvalE(3,1,2) - DvalE(2,1,3)
  ValE(2,2) = DvalE(1,1,3) - DvalE(3,1,1)
  ValE(3,2) = DvalE(2,1,1) - DvalE(1,1,2)
  ValE(1:3,2) = ValE(1:3,2)/(-ZI*OMEGA*MU)
!
!  ...1st order derivatives
  DvalE(1,2,1) = D2valE(3,1,2,1) - D2valE(2,1,3,1)
  DvalE(1,2,2) = D2valE(3,1,2,2) - D2valE(2,1,3,2)
  DvalE(1,2,3) = D2valE(3,1,2,3) - D2valE(2,1,3,3)
!
  DvalE(2,2,1) = D2valE(1,1,3,1) - D2valE(3,1,1,1)
  DvalE(2,2,2) = D2valE(1,1,3,2) - D2valE(3,1,1,2)
  DvalE(2,2,3) = D2valE(1,1,3,3) - D2valE(3,1,1,3)
!
  DvalE(3,2,1) = D2valE(2,1,1,1) - D2valE(1,1,2,1)
  DvalE(3,2,2) = D2valE(2,1,1,2) - D2valE(1,1,2,2)
  DvalE(3,2,3) = D2valE(2,1,1,3) - D2valE(1,1,2,3)
!
  DvalE(1:3,2,1:3) = DvalE(1:3,2,1:3)/(-ZI*OMEGA*MU)
!
!  ...fake 2nd order derivatives (not needed)
  D2valE(1:3,2,1:3,1:3) = ZERO
!
!  ...L2 components, derivatives not needed
  ValQ(1:3) = ValE(1:3,1)
  ValQ(4:6) = ValE(1:3,2)
!
end subroutine exact


subroutine h1_solution(X, U,GradU,Grad2U,Ut)
  use data_structure3D
  use problem
  implicit none
!-----------------------------------------------------------------------------------
  real*8, dimension(3),     intent(in)  :: X
!
  real*8,                   intent(out) :: U
  real*8, dimension(3),     intent(out) :: GradU  ! 1st derivative
  real*8, dimension(3,3),   intent(out) :: Grad2U ! 2nd derivative
  real*8,                   intent(out) :: Ut     ! time derivative
! !-----------------------------------------------------------------------------------

  real*8  :: x1,x2,x3
  real*8 :: nn        ! for EXPONENTIAL solution
  real*8 :: cn,dn     ! for singular solution
  real*8 :: tn        ! for time step
  real*8 :: np_x,np_y,np_z,f_x,f_y,f_z,df_x,df_y,df_z,ddf_x,ddf_y,ddf_z
! !-----------------------------------------------------------------------------------

! initialize variables
  U = ZERO; GradU = ZERO; Grad2U = ZERO; Ut = ZERO
! !
! !-----------------------------------------------------------------------------------
! !      D E C L A R E    S O L U T I O N    V A R I A B L E S                       |
 ! !-----------------------------------------------------------------------------------
! !
  x1 = X(1); x2 = X(2); x3 = X(3)
!--------------- 1st prob -------------------------------------------------------
  if (ISOL .eq. 1) then
    np_x=real(NPX,8); np_y=real(NPY,8); np_z=real(NPZ,8)
!
!       value
    f_x = x1**np_x
    f_y = x2**np_y
    f_z = x3**np_z
!
!       derivatives
    select case(int(np_x))
    case(0); df_x = 0.d0; ddf_x = 0.d0
    case(1); df_x = 1.d0; ddf_x = 0.d0
    case default
    df_x = np_x * x1**(np_x-1.d0)
    ddf_x = np_x * (np_x-1.d0) * x1**(np_x-2.d0)
    end select
    select case(int(np_y))
    case(0); df_y = 0.d0; ddf_y = 0.d0
    case(1); df_y = 1.d0; ddf_y = 0.d0
    case default
    df_y = np_y * x2**(np_y-1.d0)
    ddf_y = np_y * (np_y-1.d0) * x2**(np_y-2.d0)
    end select
    select case(int(np_z))
    case(0); df_z = 0.d0; ddf_z = 0.d0
    case(1); df_z = 1.d0; ddf_z = 0.d0
    case default
    df_z = np_z * x3**(np_z-1.d0)
    ddf_z = np_z * (np_z-1.d0) * x3**(np_z-2.d0)
    end select
!  .....1st order derivatives
    GradU(1)=  df_x *   f_y *   f_z
    GradU(2) =   f_x *  df_y *   f_z
    GradU(3)=   f_x *   f_y *  df_z
!
!  .....2nd order derivatives
    Grad2U(1,1) = ddf_x *   f_y *   f_z
    Grad2U(1,2) =  df_x *  df_y *   f_z
    Grad2U(1,3) =  df_x *   f_y *  df_z
    Grad2U(2,1) =  Grad2U(1,2)
    Grad2U(2,2) =   f_x * ddf_y *   f_z
    Grad2U(2,3) =   f_x *  df_y *  df_z
    Grad2U(3,1) =  Grad2U(1,3)
    Grad2U(3,2) =  Grad2U(2,3)
    Grad2U(3,3) =   f_x *   f_y * ddf_z
!
    U=f_x*f_y*f_z
!
    Ut = -U
!----------- 2nd prob -------------------------------------------------------
  elseif (ISOL .eq. 2) then
! !   solution
    U = dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
! !   1st order derivatives
    GradU(1) = PI*dcos(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
    GradU(2) = PI*dsin(PI*x1)*dcos(PI*x2)*dsin(PI*x3)
    GradU(3) = PI*dsin(PI*x1)*dsin(PI*x2)*dcos(PI*x3)

! !   2nd order derivatives
    Grad2U(1,1) =-PI**2*dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
    Grad2U(1,2) = PI**2*dcos(PI*x1)*dcos(PI*x2)*dsin(PI*x3)
    Grad2U(1,3) = PI**2*dcos(PI*x1)*dsin(PI*x2)*dcos(PI*x3)
    Grad2U(2,1) = Grad2U(1,2)
    Grad2U(2,2) =-PI**2*dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
    Grad2U(2,3) = PI**2*dsin(PI*x1)*dcos(PI*x2)*dcos(PI*x3)
    Grad2U(3,1) = Grad2U(1,3)
    Grad2U(3,2) = Grad2U(2,3)
    Grad2U(3,3) =-PI**2*dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
! !   time derivative at tn
    Ut = -U
  elseif (ISOL .eq. 3) then
!-------------- 3rd prob -------------------------------------------------------
    nn=3.d0; ! power in exponential solution
! !   solution
    U = dexp(x1**nn)*dexp(x2**nn)*dexp(x3**nn)
        ! !   1st order derivatives
    GradU(1) = nn*x1**(nn - 1.d0)*dexp(x1**nn)*dexp(x2**nn)*dexp(x3**nn)
    GradU(2) = nn*x2**(nn - 1.d0)*dexp(x1**nn)*dexp(x2**nn)*dexp(x3**nn)
    GradU(3) = nn*x3**(nn - 1.d0)*dexp(x1**nn)*dexp(x2**nn)*dexp(x3**nn)

                ! !   2nd order derivatives
    Grad2U(1,1) = (nn*x1**nn*dexp(x1**nn + x2**nn + x3**nn)*(nn + nn*x1**nn - 1.d0))/x1**2.d0
    Grad2U(1,2) = nn**2.d0*x1**(nn - 1.d0)*x2**(nn - 1.d0)*dexp(x1**nn)*dexp(x2**nn)*dexp(x3**nn)
    Grad2U(1,3) = nn**2.d0*x1**(nn - 1.d0)*x3**(nn - 1.d0)*dexp(x1**nn)*dexp(x2**nn)*dexp(x3**nn)
    Grad2U(2,1) = nn**2.d0*x1**(nn - 1.d0)*x2**(nn - 1.d0)*dexp(x1**nn)*dexp(x2**nn)*dexp(x3**nn)
    Grad2U(2,2) = (nn*x2**nn*dexp(x1**nn + x2**nn + x3**nn)*(nn + nn*x2**nn - 1.d0))/x2**2.d0
    Grad2U(2,3) = nn**2.d0*x2**(nn - 1.d0)*x3**(nn - 1.d0)*dexp(x1**nn)*dexp(x2**nn)*dexp(x3**nn)
    Grad2U(3,1) = nn**2.d0*x1**(nn - 1.d0)*x3**(nn - 1.d0)*dexp(x1**nn)*dexp(x2**nn)*dexp(x3**nn)
    Grad2U(3,2) = nn**2.d0*x2**(nn - 1.d0)*x3**(nn - 1.d0)*dexp(x1**nn)*dexp(x2**nn)*dexp(x3**nn)
    Grad2U(3,3) = (nn*x3**nn*dexp(x1**nn + x2**nn + x3**nn)*(nn + nn*x3**nn - 1.d0))/x3**2.d0
    Ut = -U
    elseif (ISOL .eq. 4) then
  endif
end subroutine h1_solution

subroutine hcurl_solution(Xp, E,dE,d2E)
  use data_structure3D
  use problem
  implicit none
!-----------------------------------------------------------------------------------
  real*8, dimension(3),     intent(in)  :: Xp
!
    VTYPE,                   intent(out) :: E   ! vector field
    VTYPE, dimension(3),     intent(out) :: dE  ! 1st derivative
    VTYPE, dimension(3,3),   intent(out) :: d2E ! 2nd derivative
    ! !-----------------------------------------------------------------------------------

    real*8  :: x1,x2,x3
    real*8 :: nn        ! for EXPONENTIAL solution
    real*8 :: cn,dn     ! for singular solution
    real*8 :: tn        ! for time step
    real*8 :: np_x,np_y,np_z,r0,k0,w0,phase,amplitude
    real*8 :: impedanceConstant, om
    VTYPE :: f_x,f_y,f_z,df_x,df_y,df_z,ddf_x,ddf_y,ddf_z
    integer :: icomp
    VTYPE  :: c2z,uz,uz_x,uz_y,uz_z,uz_xx,uz_xy,uz_xz,uz_yy,uz_yx,uz_yz
    VTYPE  :: uz_zz,uz_zy,uz_zx
    VTYPE  :: pz,pz_x,pz_y,pz_z,pz_xx,pz_xy,pz_xz,pz_yy,pz_yx,pz_yz
    VTYPE  :: pz_zz,pz_zy,pz_zx
    ! !-----------------------------------------------------------------------------------

    ! initialize variables
    E = ZERO; dE = ZERO; d2E = ZERO;
    icomp = ICOMP_EXACT
    f_x = ZERO; f_y = ZERO; f_z = ZERO
    df_x = ZERO; df_y = ZERO; df_z = ZERO
    ddf_x = ZERO; ddf_y = ZERO; ddf_z = ZERO
    !
    impedanceConstant = GAMMA_IMP

    ! !
    ! !-----------------------------------------------------------------------------------
    ! !      D E C L A R E    S O L U T I O N    V A R I A B L E S                       |
    ! !-----------------------------------------------------------------------------------
    ! !
    x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)


    !--------------- 1st prob -------------------------------------------------------
    if (ISOL .eq. 1) then

         np_x=real(NPX,8); np_y=real(NPY,8); np_z=real(NPZ,8)
!
!       value
        f_x = x1**np_x
        f_y = x2**np_y
        f_z = x3**np_z
!
!       derivatives
        select case(int(np_x))
          case(0); df_x = 0.d0; ddf_x = 0.d0
          case(1); df_x = 1.d0; ddf_x = 0.d0
          case default
            df_x = np_x * x1**(np_x-1.d0)
            ddf_x = np_x * (np_x-1.d0) * x1**(np_x-2.d0)
        end select
        select case(int(np_y))
          case(0); df_y = 0.d0; ddf_y = 0.d0
          case(1); df_y = 1.d0; ddf_y = 0.d0
          case default
            df_y = np_y * x2**(np_y-1.d0)
            ddf_y = np_y * (np_y-1.d0) * x2**(np_y-2.d0)
        end select
        select case(int(np_z))
          case(0); df_z = 0.d0; ddf_z = 0.d0
          case(1); df_z = 1.d0; ddf_z = 0.d0
          case default
            df_z = np_z * x3**(np_z-1.d0)
            ddf_z = np_z * (np_z-1.d0) * x3**(np_z-2.d0)
        end select



!----------- 2nd prob -------------------------------------------------------
!  ...a smooth solution
      elseif (ISOL .eq. 2) then

      om = OMEGA
            
      f_x= sin(om*Xp(1))!dexp(-Xp(1)**2)/20.d0!sin(OMEGA*Xp(1))
      f_y= sin(om*Xp(2))!dexp(-Xp(2)**2)/20.d0!sin(OMEGA*Xp(2))
      f_z= sin(om*Xp(3))
!
!     1st order derivatives
      df_x=(om)*cos(om*Xp(1))!-Xp(1)*dexp(-Xp(1)**2)/10.d0!(OMEGA)*cos(OMEGA*Xp(1))
      df_y=(om)*cos(om*Xp(2))!-Xp(2)*dexp(-Xp(2)**2)/10.d0!(OMEGA)*cos(OMEGA*Xp(2))
      df_z=(om)*cos(om*Xp(3))
!
!     2nd order derivatives
      ddf_x=-om**2*f_x!dexp(-Xp(1)**2)*(2.d0*Xp(1)**2-1.d0)/10.d0!-OMEGA**2*f_x
      ddf_y=-om**2*f_y!dexp(-Xp(2)**2)*(2.d0*Xp(2)**2-1.d0)/10.d0!-OMEGA**2*f_y
      ddf_z=-om**2*f_z

        !     an exponential


      elseif (ISOL .eq. 3) then
!
!       value
        f_x=1.d0
        f_y=1.d0
        f_z=exp(Xp(3))
!
!       1st order derivatives
        df_x=0.d0
        df_y=0.d0
        df_z=exp(Xp(3))
!
!       2nd order derivatives
        ddf_x=0.d0
        ddf_y=0.d0
        ddf_z=exp(Xp(3))
!

!
      elseif (ISOL .eq. 4) then
!
!     value
      f_x=1.d0
      f_y=1.d0
      f_z=Xp(3)**2 * (1.d0-Xp(3))
!
!     1st order derivatives
      df_x=0.d0
      df_y=0.d0
      df_z=2.d0*Xp(3)-3.d0*Xp(3)**2
!
!     2nd order derivatives
      ddf_x=0.d0
      ddf_y=0.d0
      ddf_z=2.d0-6.d0*Xp(3)
!
      elseif (ISOL .eq. 5) then
!
!     value
      f_x=1.d0
      f_y=1.d0
      f_z=Xp(3)**2
!
!     1st order derivatives
      df_x=0.d0
      df_y=0.d0
      df_z=2.d0*Xp(3)
!
!     2nd order derivatives
      ddf_x=0.d0
      ddf_y=0.d0
      ddf_z=2.d0

      !  ...fundamental TE10 mode for rectangular waveguide
      elseif (ISOL .eq. 6) then
      f_x=-ZI*(OMEGA/PI)*sin(PI*Xp(1))
      f_y= 1.d0
      f_z=cdexp(-ZI*OMEGA*Xp(3)*impedanceConstant)
!
!     1st order derivatives
      df_x=-ZI*(OMEGA/PI)*PI*cos(PI*Xp(1))
      df_y=0.d0
      df_z=(-ZI*OMEGA*impedanceConstant)*f_z
!
!     2nd order derivatives
      ddf_x=-PI**2*f_x
      ddf_y=0.d0
      ddf_z=(-ZI*OMEGA*impedanceConstant)*df_z

!----------- 2nd prob -------------------------------------------------------
!  ...a smooth solution
      elseif (ISOL .eq. 50) then

      f_x= Xp(1)*(1.d0-Xp(1))!dexp(-Xp(1)**2)/20.d0!sin(OMEGA*Xp(1))
      f_y= Xp(2)*(1.d0-Xp(2))!dexp(-Xp(2)**2)/20.d0!sin(OMEGA*Xp(2))
      f_z= Xp(3)*(128.d0-Xp(3))
!
!     1st order derivatives
      df_x= 1.d0-2.d0*Xp(1)!-Xp(1)*dexp(-Xp(1)**2)/10.d0!(OMEGA)*cos(OMEGA*Xp(1))
      df_y=1.d0-2.d0*Xp(2)!-Xp(2)*dexp(-Xp(2)**2)/10.d0!(OMEGA)*cos(OMEGA*Xp(2))
      df_z=(128.d0-2.d0*Xp(3))
!
!     2nd order derivatives
      ddf_x=-2.d0!dexp(-Xp(1)**2)*(2.d0*Xp(1)**2-1.d0)/10.d0!-OMEGA**2*f_x
      ddf_y=-2.d0!dexp(-Xp(2)**2)*(2.d0*Xp(2)**2-1.d0)/10.d0!-OMEGA**2*f_y
      ddf_z=-2.d0

        !     an exponential

    endif

!     !  .....1st order derivatives
         dE(1)=  df_x *   f_y *   f_z
         dE(2) =   f_x *  df_y *   f_z
         dE(3)=   f_x *   f_y *  df_z
! !
! !  .....2nd order derivatives
         d2E(1,1) = ddf_x *   f_y *   f_z
         d2E(1,2) =  df_x *  df_y *   f_z
         d2E(1,3) =  df_x *   f_y *  df_z
         d2E(2,1) =  d2E(1,2)
         d2E(2,2) =   f_x * ddf_y *   f_z
         d2E(2,3) =   f_x *  df_y *  df_z
         d2E(3,1) =  d2E(1,3)
         d2E(3,2) =  d2E(2,3)
         d2E(3,3) =   f_x *   f_y * ddf_z

         E=f_x*f_y*f_z

!...... 3D Gaussian Beam
      if (ISOL .eq. 7) then
!
!.... beam waist
          w0 = BEAM_WAIST
          k0 = OMEGA
          c2z = w0**2
          x1 = Xp(1)-0.5d0
          x2 = Xp(2)-0.5d0
          r0 = dsqrt((x1)**2+(x2)**2)
          uz = cdexp(ZI*k0*Xp(3))*cdexp(-r0**2/c2z)/c2z
    !
          uz_x = -2.d0*x1*uz/c2z
          uz_y = -2.d0*x2*uz/c2z
          uz_z = uz*(ZI*k0)
    !
          uz_xx = -2.d0*(uz+x1*uz_x)/c2z
          uz_xy = -2.d0*x1*uz_y/c2z
          uz_xz = -2.d0*x1*uz_z/c2z

          uz_yy = -2.d0*(uz+x2*uz_y)/c2z
          uz_yx = uz_xy
          uz_yz = -2.d0*x2*uz_z/c2z

          uz_zz = uz_z*ZI*k0
          uz_zx = uz_xz
          uz_zy = uz_yz

          E=uz
    !  .....1st order derivatives
          dE(1) = uz_x
          dE(2) = uz_y
          dE(3) = uz_z
!
!  .....2nd order derivatives
          d2E(1,1) = uz_xx
          d2E(1,2) = uz_xy
          d2E(1,3) = uz_xz
          d2E(2,1) = d2E(1,2)
          d2E(2,2) = uz_yy
          d2E(2,3) = uz_yz
          d2E(3,1) = d2E(1,3)
          d2E(3,2) = d2E(2,3)
          d2E(3,3) = uz_zz
      endif
!...... 3D Plane Wave in z-direction
      if (ISOL .eq. 8) then
!
!.... plane wave parameters
      phase = Xp(3)+Xp(2)+Xp(1)
      amplitude = INTENSITY_RATIO**-2
      pz = amplitude*exp(-ZI*OMEGA*phase)
      pz_x = -pz*ZI*OMEGA
      pz_y = -pz*ZI*OMEGA
      pz_z = -pz*ZI*OMEGA
      pz_xx = -pz*OMEGA**2
      pz_xy = -pz*OMEGA**2
      pz_xz = -pz*OMEGA**2
      pz_yx = -pz*OMEGA**2
      pz_yy = -pz*OMEGA**2
      pz_yz = -pz*OMEGA**2
      pz_zx = -pz*OMEGA**2
      pz_zy = -pz*OMEGA**2
      pz_zz = -pz*OMEGA**2
!
!
      E=pz
    !  .....1st order derivatives
      dE(1) = pz_x
      dE(2) = pz_y
      dE(3) = pz_z
!
!  .....2nd order derivatives
      d2E(1,1) = pz_xx
      d2E(1,2) = pz_xy
      d2E(1,3) = pz_xz
      d2E(2,1) = d2E(1,2)
      d2E(2,2) = pz_yy
      d2E(2,3) = pz_yz
      d2E(3,1) = d2E(1,3)
      d2E(3,2) = d2E(2,3)
      d2E(3,3) = pz_zz
      endif
!...... 3D Plane Wave + Gaussian beam
      if (ISOL .eq. 9) then
!
!
!.... beam waist
          w0 = BEAM_WAIST
          k0 = OMEGA
          c2z = w0**2 + ZI*2.d0*Xp(3)/k0
          r0 = dsqrt(Xp(1)**2+Xp(2)**2)
          uz = cdexp(ZI*k0*Xp(3))*cdexp(-r0**2/c2z)/c2z
    !
          uz_x = -2.d0*Xp(1)*uz
          uz_y = -2.d0*Xp(2)*uz
          uz_z = uz*(ZI*k0+(2.d0/k0)*ZI*(r0/c2z)**2-2.d0*ZI/(k0*c2z))
    !
          uz_xx = -2.d0*(uz+Xp(1)*uz_x)
          uz_xy = -2.d0*Xp(1)*uz_y
          uz_xz = -2.d0*Xp(1)*uz_z

          uz_yy = -2.d0*(uz+Xp(2)*uz_y)
          uz_yx = uz_xy
          uz_yz = -2.d0*Xp(2)*uz_z

          uz_zz = uz_z*(ZI*k0+(2.d0/k0)*ZI*(r0/c2z)**2-2.d0*ZI/(k0*c2z)) &
                + uz*(-4.d0/((k0*c2z)**2)+6.d0*(c2z**3)*(r0/k0)**2)
          uz_zx = uz_xz
          uz_zy = uz_yz

!.... plane wave parameters
          phase = Xp(3)
          amplitude = INTENSITY_RATIO**-2
          pz = amplitude*exp(-ZI*OMEGA*phase)
          pz_x = ZERO
          pz_y = ZERO
          pz_z = -pz*ZI*OMEGA
          pz_xx = ZERO
          pz_xy = ZERO
          pz_xz = ZERO
          pz_yx = ZERO
          pz_yy = ZERO
          pz_yz = ZERO
          pz_zx = ZERO
          pz_zy = ZERO
          pz_zz = -pz*OMEGA**2
!       Check if on the launching end of fiber
        if((Xp(3).eq.0.d0)) then
!     If in the core, exact is sum of Gaussian beam with plane wave
!     else, only plane wave in the cladding
        if((Xp(1).le.1.d0).and.(Xp(2).le.1.d0)) then
        E= pz+uz
    !  .....1st order derivatives
        dE(1) = pz_x+uz_x
        dE(2) = pz_y+uz_y
        dE(3) = pz_z+uz_z
!
!  .....2nd order derivatives
        d2E(1,1) = pz_xx+uz_xx
        d2E(1,2) = pz_xy+uz_xy
        d2E(1,3) = pz_xz+uz_xz
        d2E(2,1) = d2E(1,2)
        d2E(2,2) = pz_yy+uz_yy
        d2E(2,3) = pz_yz+uz_yz
        d2E(3,1) = d2E(1,3)
        d2E(3,2) = d2E(2,3)
        d2E(3,3) = pz_zz+uz_zz
        else
        E= pz
    !  .....1st order derivatives
        dE(1) = pz_x
        dE(2) = pz_y
        dE(3) = pz_z
!
!  .....2nd order derivatives
        d2E(1,1) = pz_xx
        d2E(1,2) = pz_xy
        d2E(1,3) = pz_xz
        d2E(2,1) = d2E(1,2)
        d2E(2,2) = pz_yy
        d2E(2,3) = pz_yz
        d2E(3,1) = d2E(1,3)
        d2E(3,2) = d2E(2,3)
        d2E(3,3) = pz_zz
!.....  end if for core/cladding check
        endif
!.....  end if for z=0 check
       endif
!.....  end if for ISOL=9
      endif
!...... sum of functions exact
      if (ISOL .eq. 10) then
!
!.... pz = sum of solutions
      pz = Xp(1)+Xp(2)+dsin(OMEGA*Xp(3))!dsin(OMEGA*Xp(1))+dsin(OMEGA*Xp(2))+dsin(OMEGA*Xp(3))
      pz_x = 1.d0!OMEGA*dcos(OMEGA*Xp(1))
      pz_y = 1.d0!OMEGA*dcos(OMEGA*Xp(2))
      pz_z = OMEGA*dcos(OMEGA*Xp(3))
      pz_xx = ZERO!-OMEGA**2*dsin(OMEGA*Xp(1))
      pz_xy = ZERO
      pz_xz = ZERO
      pz_yx = ZERO
      pz_yy = ZERO!-OMEGA**2*dsin(OMEGA*Xp(2))
      pz_yz = ZERO
      pz_zx = ZERO
      pz_zy = ZERO
      pz_zz = -OMEGA**2*dsin(OMEGA*Xp(3))
!
!
      E=pz
    !  .....1st order derivatives
      dE(1) = pz_x
      dE(2) = pz_y
      dE(3) = pz_z
!
!  .....2nd order derivatives
      d2E(1,1) = pz_xx
      d2E(1,2) = pz_xy
      d2E(1,3) = pz_xz
      d2E(2,1) = d2E(1,2)
      d2E(2,2) = pz_yy
      d2E(2,3) = pz_yz
      d2E(3,1) = d2E(1,3)
      d2E(3,2) = d2E(2,3)
      d2E(3,3) = pz_zz
      endif

!...... sum of functions exact
      if (ISOL .eq. 11) then
!
!.... pz = sum + product of solutions
      pz = (Xp(1)+Xp(2))*dsin(OMEGA*Xp(3))
      pz_x = dsin(OMEGA*Xp(3))
      pz_y = dsin(OMEGA*Xp(3))
      pz_z = OMEGA*(Xp(1)+Xp(2))*dcos(OMEGA*Xp(3))
      pz_xx = ZERO
      pz_xy = ZERO
      pz_xz = OMEGA*dcos(OMEGA*Xp(3))
      pz_yx = ZERO
      pz_yy = ZERO
      pz_yz = OMEGA*dcos(OMEGA*Xp(3))
      pz_zx = OMEGA*dcos(OMEGA*Xp(3))
      pz_zy = OMEGA*dcos(OMEGA*Xp(3))
      pz_zz = -OMEGA**2*(Xp(1)+Xp(2))*dsin(OMEGA*Xp(3))
!
!
      E=pz
    !  .....1st order derivatives
      dE(1) = pz_x
      dE(2) = pz_y
      dE(3) = pz_z
!
!  .....2nd order derivatives
      d2E(1,1) = pz_xx
      d2E(1,2) = pz_xy
      d2E(1,3) = pz_xz
      d2E(2,1) = d2E(1,2)
      d2E(2,2) = pz_yy
      d2E(2,3) = pz_yz
      d2E(3,1) = d2E(1,3)
      d2E(3,2) = d2E(2,3)
      d2E(3,3) = pz_zz
      endif
!...... sum of functions exact
      if (ISOL .eq. 12) then
!
!.... pz = sum of solutions
      pz = dsin(OMEGA*(Xp(1)+Xp(2)+Xp(3)))
      pz_x = OMEGA*dcos(OMEGA*(Xp(1)+Xp(2)+Xp(3)))
      pz_y = OMEGA*dcos(OMEGA*(Xp(1)+Xp(2)+Xp(3)))
      pz_z = OMEGA*dcos(OMEGA*(Xp(1)+Xp(2)+Xp(3)))
      pz_xx = -OMEGA**2*dsin((Xp(1)+Xp(2)+Xp(3)))
      pz_xy = -OMEGA**2*dsin((Xp(1)+Xp(2)+Xp(3)))
      pz_xz = -OMEGA**2*dsin((Xp(1)+Xp(2)+Xp(3)))
      pz_yx = -OMEGA**2*dsin((Xp(1)+Xp(2)+Xp(3)))
      pz_yy = -OMEGA**2*dsin((Xp(1)+Xp(2)+Xp(3)))
      pz_yz = -OMEGA**2*dsin((Xp(1)+Xp(2)+Xp(3)))
      pz_zx = -OMEGA**2*dsin((Xp(1)+Xp(2)+Xp(3)))
      pz_zy = -OMEGA**2*dsin((Xp(1)+Xp(2)+Xp(3)))
      pz_zz = -OMEGA**2*dsin((Xp(1)+Xp(2)+Xp(3)))
!
!
      E=pz
    !  .....1st order derivatives
      dE(1) = pz_x
      dE(2) = pz_y
      dE(3) = pz_z
!
!  .....2nd order derivatives
      d2E(1,1) = pz_xx
      d2E(1,2) = pz_xy
      d2E(1,3) = pz_xz
      d2E(2,1) = d2E(1,2)
      d2E(2,2) = pz_yy
      d2E(2,3) = pz_yz
      d2E(3,1) = d2E(1,3)
      d2E(3,2) = d2E(2,3)
      d2E(3,3) = pz_zz
      endif

     end subroutine hcurl_solution
