c--------------------------------------------------------------------
c
c     routine name      - elem
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns unconstrained (ordinary)
c                         stiffness matrix and load vector
c                         for the Hcurl projection
c
c
c     arguments:
c
c     in:
c             Mdle        - an element middle node number, identified
c                         with the element
c     out:
c            Itest,Itrial - flags indicating
c                         presence of corresponding load vectors
c                         and stiffness matrices
c
c-----------------------------------------------------------------------
c
c  ...this is a system routine, the head cannot be changed
      subroutine elem(Mdle, Itest,Itrial)
c
      use control
      use data_structure3D
      use assembly
#include "syscom.blk"
c
      dimension Itest(NR_PHYSA),Itrial(NR_PHYSA)
c
      Itest(1:NR_PHYSA)=0; Itrial(1:NR_PHYSA)=0
c
      select case(NODES(Mdle)%case)
c
      case(1)
c
c  .....this is user-defined routine, you are in charge of it
        Itest(1)=1; Itrial(1)=1
        call elem_hcurlProj(Mdle,BLOC(1)%nrow,
     .                  BLOC(1)%array,ALOC(1,1)%array)
c
      case default
        write(*,*) 'elem: Mdle,NODES(Mdle)%case = ',
     .              Mdle,NODES(Mdle)%case
        stop1
      end select

c
c
      end subroutine elem

      subroutine elem_hcurlProj(Mdle,MdE, ZblocE,ZalocEE)
!
      use control
      use parametersDPG
      use element_data
      use data_structure3D
      use problem
!.......no implicit statements
      implicit none
#if C_MODE
#define VTYPE  complex*16
#else
#define VTYPE double precision
#endif
!
!.......declare input/output variables
      integer,                     intent(in)  :: Mdle
      integer,                     intent(in)  :: MdE
      VTYPE, dimension(MdE),       intent(out) :: ZblocE
      VTYPE, dimension(MdE,MdE),   intent(out) :: ZalocEE

!
!.......declare edge/face type varibles
      character(len=4) :: etype,ftype
!
! ...declare element order, orientation for edges and faces
      integer, dimension(19)    :: norder
      integer, dimension(12)    :: norient_edge
      integer, dimension(6)     :: norient_face
      integer, dimension(19)    :: norderc
! ...face order
      integer, dimension(5)     :: norderf
!
! ...geometry dof (work space for nodcor)
      real*8, dimension(3,MAXbrickH) :: xnod
!

! ... variables for geometry
      real*8, dimension(3)      :: xi,x,rn
      real*8, dimension(3,2)    :: dxidt,dxdt,rt
      real*8, dimension(3,3)    :: dxdxi,dxidx
      real*8, dimension(2)      :: t
! .... number of vertices,edge,faces per element type
      integer :: nrv, nre, nrf
! ...H1 shape functions
      real*8, dimension(MAXbrickH)    :: shapH
      real*8, dimension(3,MAXbrickH)  :: gradH
! ...H(curl) shape functions
      real*8, dimension(3,MAXbrickE)  :: shapE
      real*8, dimension(3,MAXbrickE)  :: curlE
! ... store nr of dofs
      integer  :: nrdofH,nrdofE
! ...3D quadrature data
      real*8, dimension(3,MAXNINT3ADD)  :: xiloc
      real*8, dimension(MAXNINT3ADD)    :: waloc
!
! ...2D quadrature data
      real*8, dimension(2,MAXNINT2ADD)  :: tloc
      real*8, dimension(MAXNINT2ADD)    :: wtloc
!
! ...BC's flags
      integer, dimension(6,NR_PHYSA)    :: ibc
!
! ...for debug printing
      VTYPE, dimension(10)  :: zaux
!..workspace for exact
      VTYPE :: ZvalH(MAXEQNH),
     .          ZdvalH(MAXEQNH,3),Zd2valH(MAXEQNH,3,3),
     .          ZvalE(3,MAXEQNE),
     .          ZdvalE(3,MAXEQNE,3),Zd2valE(3,MAXEQNE,3,3),
     .          ZvalV(3,MAXEQNV),
     .          ZdvalV(3,MAXEQNV,3),Zd2valV(3,MAXEQNV,3,3),
     .          ZvalQ(MAXEQNQ),
     .         ZdvalQ(MAXEQNQ,3),Zd2valQ(MAXEQNQ,3,3)

! ....Maxwell load and auxiliary variables
      VTYPE, dimension(3) :: zJ,zcurlJ
      real*8, dimension(3) :: E1,curlE1,E2,curlE2
!
      integer,save :: ivis=0
      integer      :: nk
!     ..... various variables for the problem
      real*8  :: rjac,weight,wa
      integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,nint
      integer :: iflag,kE,k,iprint,l
      integer :: N,nRHS,nordP,nsign,if,ndom,info,info1,info2,info3
!
!---------------------------------------------------------------------
!
      iprint=0
!
!  ...element type
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
!  ...determine order of approximation
      call find_order(Mdle, norder)
      norderc(1:nre+nrf) = norder(1:nre+nrf)
!
! ...set the enriched order of appoximation
      select case(etype)
        case('mdlb')
        nordP = NODES(Mdle)%order+NORD_ADD*111
        norderc(nre+nrf+1) = 111
        case('mdln','mdld')
        nordP = NODES(Mdle)%order+NORD_ADD
        norderc(nre+nrf+1) = 1
        case('mdlp')
        nordP = NODES(Mdle)%order+NORD_ADD*11
        norderc(nre+nrf+1) = 11
      end select
!
!  ...determine edge and face orientations
      call find_orient(Mdle, norient_edge,norient_face)
!
!  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
!
!  ...get the element boundary conditions flags
      call find_bc(Mdle, ibc)
      iprint = 0
      if (iprint.ge.1) then
        write(*,7001) Mdle
 7001   format('elem_hcurlProj: BC FLAGS FOR Mdle = ',i5)
        do i=1,NR_PHYSA
          write(*,7002) PHYSA(i), ibc(1:nrf,i)
 7002     format('          ATTRIBUTE = ',a6,' FLAGS = ',6i2)
        enddo
      endif
      iprint = 0
!
!  ...clear space for stiffness matrix and rhsv:
      ZblocE = ZERO;
      ZalocEE = ZERO;
!-----------------------------------------------------------------------
!
!  ...element integrals...
!
!  ...use the enriched order to set the quadrature
      INTEGRATION = NORD_ADD
      call set_3Dint(etype,norder, nint,xiloc,waloc)
      INTEGRATION = 0
! ... loop through integration points
      do l=1,nint
        xi(1:3) = xiloc(1:3,l)
        wa = waloc(l)
!
!  .....determine element H1 shape functions (for geometry)
        call shape3H(etype,xi,norder,norient_edge,norient_face,
     .               nrdofH,shapH,gradH)
!
!  .....determine element H(curl) shape functions
        call shape3E(etype,xi,norder,norient_edge,norient_face,
     .               nrdofE,shapE,curlE)
!
!
!  .....geometry
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,
     .              x,dxdxi,dxidx,rjac,iflag)
!
!  .....integration weight
        weight = rjac*wa
!
!
!  ...compute RHS
      call exact(x,1,
     .            zvalH,zdvalH,zd2valH,
     .            zvalE,zdvalE,zd2valE,
     .            zvalV,zdvalV,zd2valV,
     .             zvalQ,zdvalQ,zd2valQ)
!  ...exact solution and derivatives
      zJ(1:3) = zvalE(1:3,1)
      zcurlJ(1) = zDvalE(3,1,2) - zDvalE(2,1,3)
      zcurlJ(2) = zDvalE(1,1,3) - zDvalE(3,1,1)
      zcurlJ(3) = zDvalE(2,1,1) - zDvalE(1,1,2)
!  .....loop through enriched H(curl) test functions
        do k1=1,nrdofE
          E1(1:3) =   shapE(1,k1)*dxidx(1,1:3)
     .              + shapE(2,k1)*dxidx(2,1:3)
     .              + shapE(3,k1)*dxidx(3,1:3)
          curlE1(1:3) = dxdxi(1:3,1)*curlE(1,k1)
     .                  + dxdxi(1:3,2)*curlE(2,k1)
     .                  + dxdxi(1:3,3)*curlE(3,k1)
          curlE1(1:3) = curlE1(1:3)/rjac
!
!  .......compute the RHS
          ZblocE(k1) = ZblocE(k1)
     .      + (zJ(1)*E1(1)+zJ(2)*E1(2)+zJ(3)*E1(3)
     .      + zcurlJ(1)*curlE1(1)+zcurlJ(2)*curlE1(2)
     .      +zcurlJ(3)*curlE1(3))*weight
!
!  .......loop through Hcurl trial functions
          do k2=1,nrdofE
            E2(1:3) = shapE(1,k2)*dxidx(1,1:3)
     .              + shapE(2,k2)*dxidx(2,1:3)
     .              + shapE(3,k2)*dxidx(3,1:3)
            curlE2(1:3) = dxdxi(1:3,1)*curlE(1,k2)
     .                  + dxdxi(1:3,2)*curlE(2,k2)
     .                  + dxdxi(1:3,3)*curlE(3,k2)
            curlE2(1:3) = curlE2(1:3)/rjac
!
!  .........accumulate for the extended stiffness matrix
            ZalocEE(k1,k2) = ZalocEE(k1,k2)
     .      + ((curlE1(1)*curlE2(1)+curlE1(2)*curlE2(2)
     .                   +curlE1(3)*curlE2(3))
     .        +(E1(1)*E2(1)+E1(2)*E2(2)+E1(3)*E2(3)))*weight
           enddo
        enddo
      enddo

      end subroutine elem_hcurlProj

