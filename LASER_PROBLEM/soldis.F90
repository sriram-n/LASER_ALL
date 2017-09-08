!---------------------------------------------------------------------------------------
!> Purpose : display user-defined quantity
!! @param[in] Mdle   - element (middle node) number
!! @param[in] Xi     - master element coordinates
!! @param[in] X      - physical coordinates
!! @param[in] Rn     - outward normal unit vector
!! @param[in] ZsolH  - H1    sol
!! @param[in] ZgradH - H1    grad
!! @param[in] ZsolE  - Hcurl sol
!! @param[in] ZcurlE - Hcurl curl
!! @param[in] ZsolV  - Hdiv  sol
!! @param[in] ZdivV  - Hdiv  div
!! @param[in] ZsolQ  - L2    sol
!!
!! @param[out] val   - quantity to display
!---------------------------------------------------------------------------------------
!
      subroutine soldis(Mdle,Xi,X,Rn,ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ, Val)
!
      use data_structure3D
      use problem
#if C_MODE
#define V_TYPE  complex*16
#else
#define V_TYPE double precision
#endif
!---------------------------------------------------------------------------------------
      implicit none
      integer,                      intent(in)  :: Mdle
      ! real*8,dimension(3),          intent(in)  :: Xi,X,Rn
      ! real*8,dimension(  MAXEQNH  ),intent(in)  :: ZsolH
      ! real*8,dimension(  MAXEQNH,3),intent(in)  :: ZgradH
      ! real*8,dimension(3,MAXEQNE  ),intent(in)  :: ZsolE
      ! real*8,dimension(3,MAXEQNE  ),intent(in)  :: ZcurlE
      ! real*8,dimension(3,MAXEQNV  ),intent(in)  :: ZsolV
      ! real*8,dimension(  MAXEQNV  ),intent(in)  :: ZdivV
      ! real*8,dimension(  MAXEQNQ  ),intent(in)  :: ZsolQ
      ! real*8,                       intent(out) :: Val
      ! real*8                     :: rntimesE(3)

      real*8,dimension(3),          intent(in)  :: Xi,X,Rn
      V_TYPE,dimension(  MAXEQNH  ),intent(in)  :: ZsolH
      V_TYPE,dimension(  MAXEQNH,3),intent(in)  :: ZgradH
      V_TYPE,dimension(3,MAXEQNE  ),intent(in)  :: ZsolE
      V_TYPE,dimension(3,MAXEQNE  ),intent(in)  :: ZcurlE
      V_TYPE,dimension(3,MAXEQNV  ),intent(in)  :: ZsolV
      V_TYPE,dimension(  MAXEQNV  ),intent(in)  :: ZdivV
      V_TYPE,dimension(  MAXEQNQ  ),intent(in)  :: ZsolQ
      real*8,                       intent(out) :: Val
      V_TYPE                     :: rntimesE(3)
!---------------------------------------------------------------------------------------
!  ...work space for routine 'exact'
      ! real*8 :: &
      !    zvalH(MAXEQNH), &
      !    zdvalH(MAXEQNH,3),zd2valH(MAXEQNH,3,3), &
      !    zvalE(3,MAXEQNE), &
      !    zdvalE(3,MAXEQNE,3),zd2valE(3,MAXEQNE,3,3), &
      !    zvalV(3,MAXEQNV), &
      !    zdvalV(3,MAXEQNV,3),zd2valV(3,MAXEQNV,3,3), &
      !    zvalQ(MAXEQNQ), &
      !    zdvalQ(MAXEQNQ,3),zd2valQ(MAXEQNQ,3,3)
      V_TYPE :: &
         zvalH(MAXEQNH), &
         zdvalH(MAXEQNH,3),zd2valH(MAXEQNH,3,3), &
         zvalE(3,MAXEQNE), &
         zdvalE(3,MAXEQNE,3),zd2valE(3,MAXEQNE,3,3), &
         zvalV(3,MAXEQNV), &
         zdvalV(3,MAXEQNV,3),zd2valV(3,MAXEQNV,3,3), &
         zvalQ(MAXEQNQ), &
         zdvalQ(MAXEQNQ,3),zd2valQ(MAXEQNQ,3,3)
      integer :: icase = 1, iprint
!
!===========================
! redirect to system routine
!  call soldis_system(Mdle,Xi,X,Rn,ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ, Val)
!  return
!===========================
!
      iprint=0
      if (iprint.eq.1) then
        write(*,7001) Mdle,X(1:3)
 7001   format('soldis: Mdle = ',i5,'  X = ',3f8.3)
        write(*,7002) ZsolH(1)
 7002   format('ZsolH = ',e12.5)
      endif
!
      select case (IEXACT_DISP)
      case(1)
        call exact(X, icase, &
                   zvalH,zdvalH,zd2valH, &
                   zvalE,zdvalE,zd2valE, &
                   zvalV,zdvalV,zd2valV, &
                   zvalQ,zdvalQ,zd2valQ)
        select case (ICHOOSE_COMP)
!
!  .....exact solution
        case(1)
          Val = dreal(ZvalH(1))

!  .....exact (tangential) trace
        case(2)
          call cross_product(rn,ZvalE(1:3,1), rntimesE)
          Val = dreal(rntimesE(1))
          !  .....exact (tangential) trace
        case(3)
          call cross_product(rn,ZvalE(1:3,1), rntimesE)
          Val = dreal(rntimesE(2))
          !  .....exact (tangential) trace
        case(4)
          call cross_product(rn,ZvalE(1:3,1), rntimesE)
          Val = dreal(rntimesE(3))

!  .....exact (tangential) trace
        case(5)
          call cross_product(rn,ZvalE(1:3,2), rntimesE)
          Val = dreal(rntimesE(1))
          !  .....exact (tangential) trace
        case(6)
          call cross_product(rn,ZvalE(1:3,2), rntimesE)
          Val = dreal(rntimesE(2))
          !  .....exact (tangential) trace
        case(7)
          call cross_product(rn,ZvalE(1:3,2), rntimesE)
          Val = dreal(rntimesE(3))

!  .....exact (normal) flux
        case(8)
          Val = ZvalV(1,1)*Rn(1) + ZvalV(2,1)*Rn(2) + ZvalV(3,1)*Rn(3)

!  .....exact E field
        case(9)
          Val = dreal(zvalQ(1))
        case(10)
          Val = dreal(zvalQ(2))
        case(11)
          Val = dreal(zvalQ(3))

!  .....exact H Field
        case(12)
          Val = dreal(zvalQ(4))
        case(13)
          Val = dreal(zvalQ(5))
        case(14)
          Val = dreal(zvalQ(6))
        end select
!
      case(0)
        select case (ICHOOSE_COMP)
!
!  .....approximate solution
        case(1)
          Val = dreal(ZsolH(1))
!  .....approximate (tangential) trace
        case(2)
          call cross_product(rn,ZsolE(1:3,1), rntimesE)
          Val = dreal(rntimesE(1))
          !  .....exact (tangential) trace
        case(3)
          call cross_product(rn,ZsolE(1:3,1), rntimesE)
          Val = dreal(rntimesE(2))
          !  .....exact (tangential) trace
        case(4)
          call cross_product(rn,ZsolE(1:3,1), rntimesE)
          Val = dreal(rntimesE(3))

!  .....approximate (tangential) trace
        case(5)
          call cross_product(rn,ZsolE(1:3,2), rntimesE)
          Val = dreal(rntimesE(1))
          !  .....exact (tangential) trace
        case(6)
          call cross_product(rn,ZsolE(1:3,2), rntimesE)
          Val = dreal(rntimesE(2))
          !  .....exact (tangential) trace
        case(7)
          call cross_product(rn,ZsolE(1:3,2), rntimesE)
          Val = dreal(rntimesE(3))
!  .....approximate (normal) flux
        case(8)
          Val = ZsolV(1,1)*Rn(1) + ZsolV(2,1)*Rn(2) + ZsolV(3,1)*Rn(3)
!  .....approximate E field
        case(9)
          Val = dreal(zsolQ(1))
        case(10)
          Val = dreal(zsolQ(2))
        case(11)
          Val = dreal(zsolQ(3))
!  .....approximate H Field
        case(12)
          Val = dreal(zsolQ(4))
        case(13)
          Val = dreal(zsolQ(5))
        case(14)
          Val = dreal(zsolQ(6))
        end select
      end select
      if (iprint.eq.1) then
        write(*,7003) IEXACT_DISP, ICHOOSE_COMP, Val
 7003   format('IEXACT_DISP, ICHOOSE_COMP, Val = ',2i2,e12.5)
        call pause
      endif
!
      end subroutine soldis
!
!
!
!---------------------------------------------------------------------------------------
!> Purpose : select quantities to display
!---------------------------------------------------------------------------------------
      subroutine soldis_select
!
      use control,    only: NEXACT    ! exact solution flag
      use parameters, only: NSTD_OUT  ! display file
      use problem
!
      implicit none
      integer :: iprev
!---------------------------------------------------------------------------------------
!
!===========================
!     redirect to system routine if you do not want to customize the routine
!     call soldis_select_system
!     return
!===========================
!
      select case(ICHOOSE_COMP)
      case(1,2)
        write(*,*) 'USE PREVIOUS SELECTION? 1-Y,0-N'
        read(*,*) iprev
        if (iprev.eq.1) then
          call disp_soldis(NSTD_OUT) ; return
        endif
      end select
!
      IEXACT_DISP=0
      if ((NEXACT.eq.1).or.(NEXACT.eq.2)) then
   10   write(*,*) 'DISPLAY APPROXIMATE OR EXACT SOLUTION (0/1) ?'
        read(*,*) IEXACT_DISP
        if ((IEXACT_DISP.ne.0).and.(IEXACT_DISP.ne.1)) go to 10
      endif
!
   20 write(*,*) 'SET VARIABLE: tempr(1), EEhat(2-4), HHhat(5-7), hflux(8), EHfld(9-14)'
      read(*,*) ICHOOSE_COMP
      if ((ICHOOSE_COMP.lt.1).and.(ICHOOSE_COMP.gt.5)) go to 20
!
      call disp_soldis(NSTD_OUT)
!
      end subroutine soldis_select
!
!---------------------------------------------------------------------------------------
!> Purpose : show the quantities selected to display
!---------------------------------------------------------------------------------------

      subroutine disp_soldis(Nstream)
!
      use problem,      only: IEXACT_DISP,ICHOOSE_COMP
!
      implicit none
      integer, intent(in) :: Nstream

      write(Nstream,1000)
      write(Nstream,3100)
      select case (IEXACT_DISP)
      case(1); write(Nstream,1010)
      case(2); write(Nstream,1020)
      end select
!
      select case (ICHOOSE_COMP)
      case(1); write(Nstream,3010)
      case(8); write(Nstream,3020)
      end select
!
      write(Nstream,3100)

 1000 format('DISPLAY SETUP')
 1010 format('DISPLAYING EXACT SOLUTION')
 1020 format('DISPLAYING APPROXIMATE SOLUTION')
 3010 format('DISPLAYING TEMPERATURE')
 3020 format('DISPLAYING HEAT FLUX')
 3100 format('-----------------------')
!
      end subroutine disp_soldis
