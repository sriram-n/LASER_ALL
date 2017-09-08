!--------------------------------------------------------------------

!     routine name      - elem
!
!--------------------------------------------------------------------
!
!     latest revision:  - August 17
!
!     purpose:          - routine returns unconstrained (ordinary)
!                         stiffness matrix and load vector
!                         for the L2 projection
!
!
!     arguments:
!
!     in:
!             Mdle        - an element middle node number, identified
!                         with the element
!     out:
!
!             Ites,Itrial - flags indicating
!                         presence of corresponding load vectors
!                         and stiffness matrices
!-----------------------------------------------------------------------
!
!  ...this is a system routine, the head cannot be changed
  subroutine elem(Mdle, Itest,Itrial)
!
      use control
      use data_structure3D
      use assembly
#include "syscom.blk"
!
      dimension Itest(NR_PHYSA),Itrial(NR_PHYSA)
!
      Itest(1:NR_PHYSA)=0; Itrial(1:NR_PHYSA)=0
!
      select case(NODES(Mdle)%case)
!
!  ...optimal test problem
      case(1)
!
!  .....this is user-defined routine, you are in charge of it
        Itest(1)=1; Itrial(1)=1
        call elem_l2Proj(Mdle,BLOC(1)%nrow, &
                       BLOC(1)%array,ALOC(1,1)%array)
!
      case default
        write(*,*) 'elem: Mdle,NODES(Mdle)%case = ', &
                   Mdle,NODES(Mdle)%case
        stop 1
      end select

!
!
      end subroutine elem

!--------------------------------------------------------------------
!
!     routine name      - elem_l2Proj
!
!--------------------------------------------------------------------
!
!     latest revision:  - August 17
!
!     purpose:          - routine returns unconstrained (ordinary)
!                         stiffness matrix and load vector
!                         for the L2 Projection
!
!     arguments:
!
!     in:
!             Mdle      - an element middle node number, identified
!                         with the element
!             MdQ       - column length of ZalocQQ
!     out:
!             ZblocQ - load vector
!             ZalocQQ- stiffness matrix
!
!---------------------------------------------------------------------
!
      subroutine elem_l2Proj(Mdle,MdQ, &
                              ZblocQ,ZalocQQ)
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
  integer,                     intent(in)  :: MdQ
  VTYPE, dimension(MdQ),       intent(out) :: ZblocQ
  VTYPE, dimension(MdQ,MdQ),   intent(out) :: ZalocQQ

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
! ... L2 shape functions
      real*8, dimension(MAXbrickQ)    :: shapQ
! ... store nr of dofs
      integer  :: nrdofH,nrdofQ
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
      VTYPE :: ZvalH(MAXEQNH), &
               ZdvalH(MAXEQNH,3),Zd2valH(MAXEQNH,3,3), &
               ZvalE(3,MAXEQNE), &
               ZdvalE(3,MAXEQNE,3),Zd2valE(3,MAXEQNE,3,3), &
               ZvalV(3,MAXEQNV), &
               ZdvalV(3,MAXEQNV,3),Zd2valV(3,MAXEQNV,3,3), &
               ZvalQ(MAXEQNQ), &
              ZdvalQ(MAXEQNQ,3),Zd2valQ(MAXEQNQ,3,3)

! ....Maxwell load and auxiliary variables
      VTYPE, dimension(6) :: zJ
      real*8              :: E, F
!
      integer,save :: ivis=0
      integer      :: nk
!     ..... various variables for the problem
      real*8  :: rjac,weight,wa
      integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,nint
      integer :: iflag,kE,k,iprint,l
      integer :: N,nRHS,nordP,nsign,if,ndom,info,info1,info2,info3
      integer :: n1,n2,n3,n4,n5,n6,m1,m2,m3,m4,m5,m6
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
 7001   format('elem_l2Proj: BC FLAGS FOR Mdle = ',i5)
        do i=1,NR_PHYSA
          write(*,7002) PHYSA(i), ibc(1:nrf,i)
 7002     format('          ATTRIBUTE = ',a6,' FLAGS = ',6i2)
        enddo
      endif
      iprint = 0
!
!  ...clear space for stiffness matrix and rhsv:
      ZblocQ = ZERO;
      ZalocQQ = ZERO;
!
!
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
        call shape3H(etype,xi,norder,norient_edge,norient_face, &
                    nrdofH,shapH,gradH)
!
!  ...L2 shape functions for the trial space
        call shape3Q(etype,xi,norder, nrdofQ,shapQ)
!
!
!  .....geometry
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, &
                   x,dxdxi,dxidx,rjac,iflag)
!
!  .....integration weight
        weight = rjac*wa
!
!   ......compute load vector
        call exact(x,1,                  &
                 zvalH,zdvalH,zd2valH, &
                 zvalE,zdvalE,zd2valE, &
                 zvalV,zdvalV,zd2valV, &
                 zvalQ,zdvalQ,zd2valQ)
        zJ = zvalQ(1:6)
!  .....loop through L2 test functions
        do k1=1,nrdofQ
           n1 =(k1-1)*6+1
           n2 =(k1-1)*6+2
           n3 =(k1-1)*6+3
           n4 =(k1-1)*6+4
           n5 =(k1-1)*6+5
           n6 =(k1-1)*6+6

           E = shapQ(k1)/rjac;

           ZblocQ(n1) = ZblocQ(n1) + zJ(1)*E*weight
           ZblocQ(n2) = ZblocQ(n2) + zJ(2)*E*weight
           ZblocQ(n3) = ZblocQ(n3) + zJ(3)*E*weight
           ZblocQ(n4) = ZblocQ(n4) + zJ(4)*E*weight
           ZblocQ(n5) = ZblocQ(n5) + zJ(5)*E*weight
           ZblocQ(n6) = ZblocQ(n6) + zJ(6)*E*weight
           do k2=1,nrdofQ
             m1 =(k2-1)*6+1
             m2 =(k2-1)*6+2
             m3 =(k2-1)*6+3
             m4 =(k2-1)*6+4
             m5 =(k2-1)*6+5
             m6 =(k2-1)*6+6

              F = shapQ(k2)/rjac;

              ZalocQQ(n1,m1) = ZalocQQ(n1,m1) + E*F*weight
              ZalocQQ(n2,m2) = ZalocQQ(n2,m2) + E*F*weight
              ZalocQQ(n3,m3) = ZalocQQ(n3,m3) + E*F*weight
              ZalocQQ(n4,m4) = ZalocQQ(n4,m4) + E*F*weight
              ZalocQQ(n5,m5) = ZalocQQ(n5,m5) + E*F*weight
              ZalocQQ(n6,m6) = ZalocQQ(n6,m6) + E*F*weight
          enddo
        enddo

! ... end loop over integration points
    enddo

      if (iprint.ge.1) then
        write(*,7010)
7010   format('elem_l2Proj: ZblocQ = ')
        write(*,7011) ZblocQ(1:NrdofQ)
7011   format(10e12.5)
        write(*,7012)
7012   format('elem_l2Proj: ZalocQQ = ')
        do i=1,NrdofQ
          write(*,7013) i,ZalocQQ(i,1:NrdofQ)
7013     format('i = ',i3,10(/,10(2e12.5,2x)))
        enddo
        call pause
      endif

  end subroutine elem_l2Proj









