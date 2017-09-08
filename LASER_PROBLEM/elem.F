c--------------------------------------------------------------------
c
c     routine name      - elem
c
c--------------------------------------------------------------------
c
c     latest revision:  - Mar 17
c
c     purpose:          - routine returns unconstrained (ordinary)
c                         stiffness matrix and load vector
c                         for either:
c                         1) single step of transient heat problem
c                         2) multiple steps of transient heat problem
c                         3) steady state Maxwell problem
c
c     arguments:
c
c     in:
c             Mdle      - an element middle node number, identified
c                         with the element
c
c     out:    Itest     -
c             Itrial
c
c-----------------------------------------------------------------------
c
      subroutine elem(Mdle, Itest,Itrial)
c
      use control
      use data_structure3D
      use assembly
      use problem
#include "syscom.blk"
c
      dimension Itest(NR_PHYSA),Itrial(NR_PHYSA)
c
      Itest(1:NR_PHYSA)=0 ; Itrial(1:NR_PHYSA)=0
c
      select case(NODES(Mdle)%case)
c
c     node supports all physical attributes
c     4 physical attributes: case = 2^4-1 = 15
      case(15)
c
        select case(NO_PROBLEM)
c
c  .....single step of transient heat equation (1)
c  .....multiple steps of transient heat equation (2)
          case(1,2)
          Itest(1)=1 ; Itest(3)=1
          Itrial(1)=1; Itrial(3)=1
          call elem_dpgHeat(Mdle,
     .            BLOC(1)%nrow,BLOC(3)%nrow,
     .            BLOC(1)%array,ALOC(1,1)%array,ALOC(1,3)%array,
     .            BLOC(3)%array,ALOC(3,1)%array,ALOC(3,3)%array)
c
c  ...... time harmonic Maxwell (3)
          case(3)
          Itest(2)=1 ; Itest(4)=1;
          Itrial(2)=1; Itrial(4)=1;
          call elem_dpgMaxwell(Mdle,
     .            BLOC(2)%nrow,BLOC(4)%nrow,
     .            BLOC(2)%array,ALOC(2,2)%array,ALOC(2,4)%array,
     .            BLOC(4)%array,ALOC(4,2)%array,ALOC(4,4)%array)

        end select
c
      case default
        write(*,*) 'elem: Mdle,NODES(Mdle)%case = ',
     .              Mdle,NODES(Mdle)%case
        stop
      endselect
c
c
      end subroutine elem
c
c--------------------------------------------------------------------
c
c     routine name      - elem_dpgHeat
c
c--------------------------------------------------------------------
c
c     latest revision:  - Mar 17
c
c     purpose:          - routine returns unconstrained (ordinary)
c                         stiffness matrix and load vector
c                         for the primal H1 formulation for a single
c                         step of transient heat equation
c
c     arguments:
c
c     in:
c             Mdle      - an element middle node number, identified
c                         with the element
c             MdH       - column length of ZalocHH,ZalocHV
c             MdV       - column length of ZalocVH,ZalocVV
c     out:
c             ZblocH,ZblocV - load vectors
c             ZalocHH,ZalocHV,ZalocVH,ZalocVV - stiffness matrices
c
c---------------------------------------------------------------------
c
      subroutine elem_dpgHeat(Mdle,MdH,MdV,
     .                      ZblocH,ZalocHH,ZalocHV,
     .                      ZblocV,ZalocVH,ZalocVV)
c
      use control
      use parametersDPG
      use element_data
      use data_structure3D
      use DPGLaser
      use problem
#include "syscom.blk"
c
      dimension ZblocH(MdH),ZalocHH(MdH,MdH),ZalocHV(MdH,MdV),
     .          ZblocV(MdV),ZalocVH(MdV,MdH),ZalocVV(MdV,MdV)
c
      character(len=4) :: etype,ftype
c
c  ...element order, orientation for edges and faces
      dimension norder(19),norient_edge(12),norient_face(6)
c
c  ...face order
      dimension norderf(5)
c
c  ...geometry dof (work space for nodcor)
      dimension xnod(3,MAXbrickH)
c
c  ...solution dof (work space for solelm)
      dimension zdofH(MAXEQNH,MAXbrickH)
      dimension zdofE(MAXEQNE,MAXbrickE)
      dimension zdofV(MAXEQNV,MAXbrickV)
      dimension zdofQ(MAXEQNQ,MAXbrickQ)
c  .... current L2 solution zsolQ
      dimension zsolQ(6)
c
c  ...geometry
      dimension xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3),
     .          dxdt(3,2),rt(3,2),rn(3),t(2)
c
c  ...H1 shape functions
      dimension shapH(MAXbrickH),gradH(3,MAXbrickH)

c  ...Hcurl shape functions
      dimension shapE(3,MAXbrickE),curlE(3,MAXbrickE)
c
c  ...Hdiv shape functions
      dimension shapV(3,MAXbrickV),divV(MAXbrickV)
c  ...L2 shape functions
      dimension shapQ(MAXbrickQ)
c
c  ...3D quadrature data
      dimension xiloc(3,MAXNINT3ADD),waloc(MAXNINT3ADD)
c
c  ...2D quadrature data
      dimension tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
c
c  ...BC's flags
      dimension ibc(6,NR_PHYSA)
c
c  ...derivatives wrt physical coordinates, flux
      dimension dv1(3),dv2(3),vec(3)
c
c  ...for debug printing
      dimension zaux(10)

c  ....Maxwell load
      dimension zJ(3)
c


      integer,save :: ivis=0
      character uplo
      nk(k1,k2) = (k2-1)*k2/2+k1
c
c---------------------------------------------------------------------
c
      iprint=0
c
      ivis=1
      if (ivis.eq.0) then
        write(*,*) 'elem_dpgHeat: NO_PROBLEM    = ',NO_PROBLEM
        write(*,*) 'elem_dpgHeat: INNER_PRODUCT = ',INNER_PRODUCT
        call pause
        ivis=1
      endif
c
c  ...element type
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
c
c  ...determine order of approximation
      call find_order(Mdle, norder)
c
c  ...set the enriched order of appoximation
      select case(etype)
      case('mdlb')        ; nordP = NODES(Mdle)%order + NORD_ADD*111
      case('mdln','mdld') ; nordP = NODES(Mdle)%order + NORD_ADD*1
      case('mdlp')        ; nordP = NODES(Mdle)%order + NORD_ADD*11
      end select
c
c  ...determine edge and face orientations
      call find_orient(Mdle, norient_edge,norient_face)
c
c  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
c
c  ...determine dof for current solution if running multiple
c  ... steps of heat equation
      select case(NO_PROBLEM)
      case(2)
        call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
      end select
c
c  ...determine element size and scale correction
      h_elem = min(abs(xnod(1,2)-xnod(1,1)),abs(xnod(2,4)-xnod(2,1)),
     .             abs(xnod(3,5)-xnod(3,1)))


c -----------
c -----------  CHECK OUT INNER PRODUCT SCALING FOR HEAT STEPPING
c -----------
c      select case(INNER_PRODUCT)
c      case(2)
c        scale = min(1.d0,EPS_DIFF**3/h_elem**2)
c      case(3)
c        scale =  min(1.d0,EPS_DIFF**2/h_elem**2)
c      end select
c      scale = 1.d0
c -----------
c -----------  CHECK OUT INNER PRODUCT SCALING FOR HEAT STEPPING
c -----------
c
c  ...get the element boundary conditions flags
      call find_bc(Mdle, ibc)
      if (iprint.ge.1) then
        write(*,7001) Mdle
 7001   format('elem_dpgHeat: BC FLAGS FOR Mdle = ',i5)
        do i=1,NR_PHYSA
          write(*,7002) PHYSA(i), ibc(1:nrf,i)
 7002     format('          ATTRIBUTE = ',a6,' FLAGS = ',6i2)
        enddo
      endif
c
c  ...clear space for stiffness matrix and rhsv:
      ZblocH = ZERO; ZblocV = ZERO
      ZalocHH = ZERO; ZalocHV = ZERO; ZalocVH = ZERO; ZalocVV = ZERO
c
c     extended load vector and extended stiffness matrices
      BLOADH=ZERO ; STIFFHH=ZERO ; STIFFHV=ZERO
c
c     Gram matrix
      AP_Heat=ZERO
c
c
c-----------------------------------------------------------------------
c     E L E M E N T    I N T E G R A L S                               |
c-----------------------------------------------------------------------
c
c  ...use the enriched order to set the quadrature
      INTEGRATION = NORD_ADD
      call set_3Dint(etype,norder, nint,xiloc,waloc)
      INTEGRATION = 0
c
c  ...loop over integration points
      do l=1,nint
c
        xi(1:3)=xiloc(1:3,l) ; wa=waloc(l)
c
c  .....H1 shape functions
        call shape3H(etype,xi,norder,norient_edge,norient_face,
     .               nrdofH,shapH,gradH)
c
c ....... element L2 solution for heat Load calculation
c  .....determine element L2 shape functions in LASER_MODE
        if(LASER_MODE.eq.1) then
          call shape3Q(etype,xi,norder, nrdofQ,shapQ)
        endif
c
c  .....discontinuous H1 shape functions
        call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
c       write(*,*) 'nrdofHH is: ', nrdofHH
c        call pause
c
c  .....geometry computations
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,
     .              x,dxdxi,dxidx,rjac,iflag)
c
c..................................................................
c..................................................................
c.........COMPUTE H1 and Hcurl Solutions for Heat stepping and
c ......... Load for Heat equation
c ......... heat load = |E|^2
c..................................................................
c..................................................................
c  .....compute the H1 solution at the point
        zsolH = ZERO
        do kH=1,nrdofH
          zsolH = zsolH + zdofH(1,kH)*shapH(kH)
        enddo
c  .....compute the Hcurl solution at the point if LASER_MOD
        if(LASER_MODE.eq.1) then
          zsolQ = ZERO
          do kE=1,nrdofQ
           zsolQ(1:6) = zsolQ(1:6) + zdofQ(1:6,kE)*shapQ(kE)
          enddo
        endif
c
c  .....integration weight
        weight = rjac*wa
c
c  .....get the RHS
        call getf(Mdle,x, zfval,zJ)
c
c  .....loop through enriched H1 test functions
        do k1=1,nrdofHH
c
          v1 = shapHH(k1)
c
          dv1(1:3) = gradHH(1,k1)*dxidx(1,1:3)
     .             + gradHH(2,k1)*dxidx(2,1:3)
     .             + gradHH(3,k1)*dxidx(3,1:3)
c
c         -- EXTENDED LOAD VECTOR --
c
          select case(NO_PROBLEM)
          case(1)
            BLOADH(k1) = BLOADH(k1) + zfval*v1*weight
          case(2)
            if(LASER_MODE.eq.1) then
              zfval = DELTAT*((zsolQ(1)**2+zsolQ(2)**2
     .                       +zsolQ(3)**2)**0.5+1.d0)
            endif
            BLOADH(k1) = BLOADH(k1) + (zfval+zsolH)*v1*weight
          end select
c
c  .......loop through enriched H1 trial functions
          do k2=k1,nrdofHH
c
            v2 = shapHH(k2)
c
            dv2(1:3) = gradHH(1,k2)*dxidx(1,1:3)
     .               + gradHH(2,k2)*dxidx(2,1:3)
     .               + gradHH(3,k2)*dxidx(3,1:3)
c
c           -- GRAM MATRIX --
c
c           determine index in triagular format
            k = nk(k1,k2)
c
c  .........problem-dependent inner product
            select case(INNER_PRODUCT)
            case(1)
              AP_Heat(k) = AP_Heat(k)
     .              + (dv1(1)*dv2(1)+dv1(2)*dv2(2)+dv1(3)*dv2(3)
     .                +v1*v2)*weight
            case default; stop1
            end select
          enddo
c
c ........loop through H1 trial functions
          do k2=1,nrdofH
            v2 = shapH(k2)
            dv2(1:3) = gradH(1,k2)*dxidx(1,1:3)
     .               + gradH(2,k2)*dxidx(2,1:3)
     .               + gradH(3,k2)*dxidx(3,1:3)
c
c           -- EXTENDED HH STIFFNESS MATRIX --
              STIFFHH(k1,k2) = STIFFHH(k1,k2)
     .      +(KAPPA*DELTAT*(dv1(1)*dv2(1)+dv1(2)*dv2(2)+dv1(3)*dv2(3))
     .                +  v1*v2 )*weight
           enddo
        enddo
      enddo
c
c  ...printing
      iprint = 0
      if (iprint.eq.1) then
        write(*,*) 'elem_dpgHeat: GRAM MATRIX = '
        do i=1,10
          do j=1,i-1
            zaux(j) = AP_Heat(nk(j,i))
          enddo
          do j=i,10
            zaux(j) = AP_Heat(nk(i,j))
          enddo
          write(*,7011) zaux
        enddo
        call pause
      endif


c
c
c-----------------------------------------------------------------------
c     B O U N D A R Y    I N T E G R A L S                             |
c-----------------------------------------------------------------------
c
c  ...loop through element faces
      do if=1,nrf
c
c  .....sign factor to determine the OUTWARD normal unit vector
        nsign = nsign_param(etype,if)
c
c  .....face type
        ftype = face_type(etype,if)
c
c  .....face order of approximation
        call face_order(etype,if,norder, norderf)
c
c  .....set 2D quadrature
        INTEGRATION = NORD_ADD
        call set_2Dint(ftype,norderf, nint,tloc,wtloc)
        INTEGRATION = 0
c
c  .....loop through integration points
        do l=1,nint
c
c  .......face coordinates
          t(1:2) = tloc(1:2,l)
c
c  .......face parametrization
          call face_param(etype,if,t, xi,dxidt)
c
c  .......determine discontinuous H1 shape functions
          call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
c
c  .......determine element H1 shape functions (for geometry)
          call shape3H(etype,xi,norder,norient_edge,norient_face,
     .                 nrdofH,shapH,gradH)
c
c  .......determine element Hdiv shape functions (for fluxes)
          call shape3V(etype,xi,norder,norient_face,
     .                 nrdofV,shapV,divV)
c
c  .......geometry
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,
     .                 x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
          weight = bjac*wtloc(l)
c
c  .......loop through enriched H1 test functions
          do k1=1,nrdofHH
            v1 = shapHH(k1)
c
c  .........loop through H(div) trial functions
            do k2=1,nrdofV
c
c  ...........normal component (Piola transformation at work!)
              vec(1:3) =  dxdxi(1:3,1)*shapV(1,k2)
     .                   +dxdxi(1:3,2)*shapV(2,k2)
     .                   +dxdxi(1:3,3)*shapV(3,k2)
              vec(1:3) = vec(1:3)/rjac
              v2n = vec(1)*rn(1)+vec(2)*rn(2)+vec(3)*rn(3)
c
c             -- EXTENDED HV STIFFNESS MATRIX --
c
              STIFFHV(k1,k2) = STIFFHV(k1,k2)
     .                       - v1*v2n*weight
            enddo
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c
c  ...factorize the Gram matrix
      uplo = 'U'

      call ZPPTRF(uplo, nrdofHH, AP_Heat, info)
      if (info.ne.0) then
        write(*,*) 'elem_dpgHeat: info = ',info ; stop
      endif

c! ....................................................................
c! ...Alternative construction of normal matrix
c! ! ....................................................................
c
c!  ...total trial dof for the element
        nrTEST = nrdofHH
      call celndof(NODES(Mdle)%type,norder,
     .             nrdofH,nrdofE,nrdofV,nrdofQ)

      i1 = MAXtestH ; j1 = nrdofH ; j2 = nrdofV

      STIFF_ALLH(1:i1,1:j1) = STIFFHH(1:i1,1:j1)
      STIFF_ALLH(1:i1,j1+1:j1+j2) = STIFFHV(1:i1,1:j2)
      STIFF_ALLH(1:i1,j1+j2+1) = BLOADH(1:i1)

c
      uplo  = 'U' ; trans = 'C' ;   diag  = 'N'
      N     = nrTEST
      NRHS  = nrdofH + nrdofV + 1
c
c
      call ZTPTRS(uplo,trans,diag,N,NRHS,AP_heat,STIFF_ALLH,
     .             MAXtestH,info)
c
c
      call ZHERK(uplo,trans,NRHS,nrTEST,ZONE,STIFF_ALLH,
     .            MAXtestH,ZERO,
     .            STIFF_ALLH(1:NRHS,1:NRHS),NRHS)

c     ZHERK for complex case

      do i=1,NRHS-1
        STIFF_ALLH(i+1:NRHS,i) = conjg(STIFF_ALLH(i,i+1:NRHS))
      enddo
c
      ZblocH(1:j1) = STIFF_ALLH(1:j1,j1+j2+1)
      ZblocV(1:j2) = STIFF_ALLH(j1+1:j1+j2,j1+j2+1)
c
      ZalocHH(1:j1,1:j1) = STIFF_ALLH(1:j1,1:j1)
      ZalocHV(1:j1,1:j2) = STIFF_ALLH(1:j1,j1+1:j1+j2)
c
      ZalocVH(1:j2,1:j1) = STIFF_ALLH(j1+1:j1+j2,1:j1)
      ZalocVV(1:j2,1:j2) = STIFF_ALLH(j1+1:j1+j2,j1+1:j1+j2)
c
c
#if C_MODE
 7011   format(6(2e10.3,2x))
#else
 7011   format(10e12.5)
#endif
c

c
      end subroutine elem_dpgHeat





c--------------------------------------------------------------------
c
c     routine name      - elem_dpgMaxwell
c
c--------------------------------------------------------------------
c
c     latest revision:  - Mar 17
c
c     purpose:          - routine returns unconstrained (ordinary)
c                         stiffness matrix and load vector
c                         for the UW DPG formulation for Maxwell
c                         equations
c
c     arguments:
c
c     in:
c             Mdle      - an element middle node number, identified
c                         with the element
c             MdE       - column length of ZalocEE,ZalocEQ
c             MdQ       - column length of ZalocQE,ZalocQQ
c     out:
c             ZblocE,ZblocQ - load vectors
c             ZalocEE,ZalocEQ,ZalocQE,ZalocQQ - stiffness matrices
c
c---------------------------------------------------------------------
c
      subroutine elem_dpgMaxwell(Mdle,MdE,MdQ,
     .                         ZblocE,ZalocEE,ZalocEQ,
     .                         ZblocQ,ZalocQE,ZalocQQ)
c
      use control
      use parametersDPG
      use element_data
      use data_structure3D
      use DPGLaser
      use problem
#include "syscom.blk"
c
      dimension ZblocE(MdE),ZalocEE(MdE,MdE),ZalocEQ(MdE,MdQ),
     .          ZblocQ(MdQ),ZalocQE(MdQ,MdE),ZalocQQ(MdQ,MdQ)
c
      character(len=4) :: etype,ftype
c
c  ...element order, orientation for edges and faces
      dimension norder(19),norient_edge(12),norient_face(6),norderc(19)
c
c  ...face order
      dimension norderf(5)
c
c  ...geometry dof
      dimension xnod(3,MAXbrickH)

c  ...solution dof (work space for solelm)
      dimension zdofH(MAXEQNH,MAXbrickH)
      dimension zdofE(MAXEQNE,MAXbrickE)
      dimension zdofV(MAXEQNV,MAXbrickV)
      dimension zdofQ(MAXEQNQ,MAXbrickQ)
c  .... current L2 solution zsolQ
      dimension zsolQ(6)
c  .... current H solution zsolH
      !V_TYPE :: zsolH
c
c  ...geometry
      dimension xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3),
     .          dxdt(3,2),rt(3,2),rn(3),t(2)
c
c  ...H1 shape functions (for geometry)
      dimension shapH(MAXbrickH),gradH(3,MAXbrickH)
c
c  ...Hcurl shape functions
      dimension shapE(3,MAXbrickE),curlE(3,MAXbrickE)
c
c  ...L2 shape functions
      dimension shapQ(MAXbrickQ)
c
c  ...3D quadrature data
      dimension xiloc(3,MAXNINT3ADD),waloc(MAXNINT3ADD)
c
c  ...2D quadrature data
      dimension tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
c
c  ...BC's flags
      dimension ibc(6,NR_PHYSA)
c
c  ...impressed volume current
      dimension zJ(3), zImp(3), qq(3), p(3),rntimesp(3),rn2timesp(3)
      dimension uE(3), rntimesuE(3)
c
c  ...shape functions and their curl after Piola transform, rotated
c     tangential flux
      dimension E1(3),curlE1(3),E2(3),curlE2(3),rntimesE(3)
c
c  ...for debug printing
      dimension aux(10)
c
c  ...for debugging
      dimension zE(24),zQ(6)
c  ...for lapack eigensolve
      complex*16, allocatable :: Z(:,:), WORK(:)
      real*8, allocatable     :: W(:),   RWORK(:)
      integer, allocatable    :: IWORK(:)
      real*8  :: bg_gain, ion_gain
c

      character uplo
      nk(k1,k2) = (k2-1)*k2/2+k1
c
c---------------------------------------------------------------------
c
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'elem: Mdle = ',Mdle
      endif
c
c  ...element type
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
c
c  ...determine order of approximation
      call find_order(Mdle, norder)
      norderc(1:nre+nrf) = norder(1:nre+nrf)

c
c  ...set the enriched order of appoximation
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
c
c  ...determine edge and face orientations
      call find_orient(Mdle, norient_edge,norient_face)
c
c  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
c  ... get current solution dofs
c
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
c
      select case(NO_PROBLEM)
      case(3)
        call find_domain(Mdle, ndom)
        select case(GEOM_NO)
        case(1)
          bg_gain = 1.d0
        case(2)
ccccc          bg_gain = REF_INDEX_CORE**2
          bg_gain = 1.d0
        case(3)
          select case(ndom)
          case(1,2,3,4)
            bg_gain = REF_INDEX_CORE**2
          case(5,6,7,8)
            bg_gain = REF_INDEX_CLAD**2
          end select
c          write(*,*)'ndom, bg_gain is: ', ndom, bg_gain
        case(4)
          bg_gain = REF_INDEX_CORE**2
        case(5)
          select case(ndom)
          case(1,2,3,4,5)
            bg_gain = REF_INDEX_CORE**2
          case(6,7,8,9)
            bg_gain = REF_INDEX_CLAD**2
          end select
        end select
      end select

c
c  ...get the element boundary conditions flags
      call find_bc(Mdle, ibc)
      if (iprint.ge.1) then
        write(*,7001) Mdle
 7001   format('elem: BC FLAGS FOR Mdle = ',i5)
        do i=1,NR_PHYSA
          write(*,7002) PHYSA(i), ibc(1:nrf,i)
 7002     format('     ATTRIBUTE = ',a6,' FLAGS = ',6i2)
        enddo
      endif
c
c  ...clear space for stiffness matrix and rhsv:
      ZblocE = ZERO; ZblocQ = ZERO
      ZalocEE = ZERO; ZalocEQ = ZERO; ZalocQE = ZERO; ZalocQQ = ZERO
c
c  ...clear space for auxiliary matrices
      BLOADE = ZERO; STIFFEE = ZERO; STIFFEQ = ZERO; AP_Maxwell = ZERO
c
c  ...shortcut: adjust frequency
c
c  ...element size:
      h = min(abs(xnod(1,2)-xnod(1,1)),abs(xnod(2,4)-xnod(2,1)),
     .        abs(xnod(3,5)-xnod(3,1)))
      omeg = min(OMEGA,6.d0/h)
c
c  ...auxiliary constant1
       zaux = ZI*OMEGA*bg_gain*EPSILON + SIGMA
       za = ZI*omeg*bg_gain*EPSILON - SIGMA

c
c-----------------------------------------------------------------------
c
c  ...element integrals...
c
c  ...use the enriched order to set the quadrature
      INTEGRATION = NORD_ADD
      call set_3Dint(etype,norder, nint,xiloc,waloc)
      INTEGRATION = 0
      do l=1,nint
        xi(1:3) = xiloc(1:3,l)
        wa = waloc(l)
c
c  .....determine element H1 shape functions (for geometry)
        call shape3H(etype,xi,norder,norient_edge,norient_face,
     .               nrdofH,shapH,gradH)
c
c  .....determine element H(curl) shape functions
        call shape3E(etype,xi,norder,norient_edge,norient_face,
     .               nrdofE,shapE,curlE)
c  .....determine element L2 shape functions
        call shape3Q(etype,xi,norder, nrdofQ,shapQ)
c
c  ..... find current domain number to update ion_gain in core only
        if(GEOM_NO.ne.1) then
          call find_domain(Mdle, ndom)
          select case(ndom)
          case(1,2,3,4)
c  ... if solving nonlinear ion_gain then update ion_gain
            if(NONLINEAR_FLAG.eq.1) then

c  .... compute current H(curl) solution to update ion_gain
              zsolQ = ZERO
              do kE=1,nrdofQ
               zsolQ(1:6) = zsolQ(1:6) + zdofQ(1:6,kE)*shapQ(kE)
              enddo
c  ... evaluate gain function
              EfieldNorm = (zsolQ(1)**2+zsolQ(2)**2
     .                       +zsolQ(3)**2)**0.5
c        write(*,*) 'from elem.F: EfieldNorm = '
              call get_gainFunction(EfieldNorm, gainFunction)
cc  ...  evaluate ion_gain = REF_INDEX*gainFunction = sqrt(bg_gain)*gainFunction
              ion_gain = sqrt(bg_gain)*gainFunction
cc  ...  evaluate auxiliary constant1
              zaux = ZI*OMEGA*bg_gain*EPSILON + ion_gain*SIGMA
              za = ZI*omeg*bg_gain*EPSILON - ion_gain*SIGMA
            endif
           end select
        endif
c
c
c  .....determine discontinuous H(curl) shape functions
        call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
c
c  .....geometry
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,
     .              x,dxdxi,dxidx,rjac,iflag)
c
c..................................................................
c..................................................................
c.........COMPUTE H1 Solution for Maxwell
c ......... current heat solution u
c..................................................................
c  .....compute the H1 solution at the point
        alpha_scale = 1.d0
        if (LASER_MODE.eq.1) then
        call find_domain(Mdle, ndom)
        select case(GEOM_NO)
        case(1)
          alpha_scale = REF_INDEX_CORE
        case(2)
          alpha_scale = REF_INDEX_CORE
        case(3)
          select case(ndom)
          case(1,2,3,4)
            alpha_scale = REF_INDEX_CORE
          case(5,6,7,8)
            alpha_scale = REF_INDEX_CLAD
          end select
        case(4)
          alpha_scale = REF_INDEX_CORE
        case(5)
          select case(ndom)
          case(1,2,3,4,5)
            alpha_scale = REF_INDEX_CORE
          case(6,7,8,9)
            alpha_scale = REF_INDEX_CLAD
          end select
        end select
        !write(*,*) 'GEOM_NO is : ', GEOM_NO
        !write(*,*) 'domain number is: ', ndom
        !write(*,*) ' alpha_scale is :', alpha_scale
        endif


        if(LASER_MODE.eq.1) then
          zsolH = ZERO
          do kH=1,nrdofH
            zsolH = zsolH + zdofH(1,kH)*shapH(kH)
          enddo
          alpha_scale = alpha_scale*(alpha_scale+2.d0*NTHETA*zsolH)
        endif

c  .....integration weight
        weight = rjac*wa
c
c  .....get the RHS
        call getf(Mdle,x, zfval,zJ)
c
c  .....loop through enriched H(curl) test functions
        do k1=1,nrdofEE
c
c  .......Piola transform
          E1(1:3) = shapEE(1,k1)*dxidx(1,1:3)
     .            + shapEE(2,k1)*dxidx(2,1:3)
     .            + shapEE(3,k1)*dxidx(3,1:3)
          curlE1(1:3) = dxdxi(1:3,1)*curlEE(1,k1)
     .                + dxdxi(1:3,2)*curlEE(2,k1)
     .                + dxdxi(1:3,3)*curlEE(3,k1)
          curlE1(1:3) = curlE1(1:3)/rjac
c
c  .......compute the RHS
          k = 2*k1
          BLOADE(k) = BLOADE(k)
     .      + (zJ(1)*E1(1)+zJ(2)*E1(2)+zJ(3)*E1(3))*weight
c
c  .......loop through enriched H(curl) trial functions
          do k2=k1,nrdofEE
c
c  .........Piola transform
            E2(1:3) = shapEE(1,k2)*dxidx(1,1:3)
     .              + shapEE(2,k2)*dxidx(2,1:3)
     .              + shapEE(3,k2)*dxidx(3,1:3)
            curlE2(1:3) = dxdxi(1:3,1)*curlEE(1,k2)
     .                  + dxdxi(1:3,2)*curlEE(2,k2)
     .                  + dxdxi(1:3,3)*curlEE(3,k2)
            curlE2(1:3) = curlE2(1:3)/rjac
c
c  .........auxiliary quantities
            CC = curlE1(1)*curlE2(1) + curlE1(2)*curlE2(2)
     .         + curlE1(3)*curlE2(3)
            EE = E1(1)*E2(1) + E1(2)*E2(2) + E1(3)*E2(3)
            CE = curlE1(1)*E2(1) + curlE1(2)*E2(2) + curlE1(3)*E2(3)
            EC = E1(1)*curlE2(1) + E1(2)*curlE2(2) + E1(3)*curlE2(3)
c
c  .........accumulate for the test stiffness matrix
            select case(INNER_PRODUCT)
            case(1)
              k = nk(2*k1-1,2*k2-1)
              AP_Maxwell(k) = AP_Maxwell(k)
     .              + (CC + (abs(za)**2 + 1.d0)*EE)*weight
ccc     .              + (CC + (omeg**2*MU**2 + 1.d0)*EE)*weight
              k = nk(2*k1-1,2*k2  )
              AP_Maxwell(k) = AP_Maxwell(k)
     .              + (-ZI*omeg*MU*CE - conjg(za)*EC)*weight
ccc     .              + (za*CE+ZI*omeg*MU*EC)*weight
              if (k1.ne.k2) then
                k = nk(2*k1  ,2*k2-1)
                AP_Maxwell(k) = AP_Maxwell(k)
     .                + (-za*CE + ZI*omeg*MU*EC)*weight
ccc     .              - (za*EC + ZI*omeg*MU*CE)*weight
              endif
              k = nk(2*k1  ,2*k2  )
              AP_Maxwell(k) = AP_Maxwell(k)
     .         + (CC + ((omeg*MU)**2 + 1.d0)*EE)*weight
ccc     .         + (CC + (abs(za)**2 + 1.d0)*EE)*weight
c
            end select
          enddo
c
c  .......loop through L2 trial functions
          do k2=1,nrdofQ
            q = shapQ(k2)/rjac
c
c  .........accumulate for the extended stiffness matrix
            do icomp=1,3
              k = (k2-1)*6 + icomp
              STIFFEQ(2*k1-1,k) = STIFFEQ(2*k1-1,k)
     .            + curlE1(icomp)*q*weight
              k = (k2-1)*6 + 3+ icomp
              STIFFEQ(2*k1-1,k) = STIFFEQ(2*k1-1,k)
     .            + ZI*OMEGA*MU*E1(icomp)*q*weight
              k = (k2-1)*6 + icomp
              STIFFEQ(2*k1  ,k) = STIFFEQ(2*k1  ,k)
     .            -zaux*alpha_scale*E1(icomp)*q*weight
              k = (k2-1)*6 + 3+ icomp
              STIFFEQ(2*k1  ,k) = STIFFEQ(2*k1  ,k)
     .            + curlE1(icomp)*q*weight
            enddo
          enddo
        enddo
      enddo
      iprint = 0
      if (iprint.eq.1) then
        write(*,*) 'elem: AP_Maxwell = '
        do i=1,10
          do j=1,i-1
            aux(j) = AP_Maxwell(nk(j,i))
          enddo
          do j=i,10
            aux(j) = AP_Maxwell(nk(i,j))
          enddo
          write(*,7011) aux
        enddo
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...boundary integrals
c
        if(GEOM_NO.eq.1) then
          !impedanceConstant = 1.d0
          impedanceConstant = 1.d0/sqrt(1.d0-(PI**2/OMEGA**2))
        else
          impedanceConstant = 1.d0
        endif
c  ...loop through element faces
      do if=1,nrf
c
c  .....sign factor to determine the OUTWARD normal unit vector
        nsign = nsign_param(etype,if)
c
c  .....face type
        ftype = face_type(etype,if)
c
c  .....face order of approximation
        call face_order(etype,if,norder, norderf)
c
c  .....set 2D quadrature
        INTEGRATION = NORD_ADD
        call set_2Dint(ftype,norderf, nint,tloc,wtloc)
        INTEGRATION = 0
c
c  .....loop through integration points
        do l=1,nint
c
c  .......face coordinates
          t(1:2) = tloc(1:2,l)
c
c  .......face parametrization
          call face_param(etype,if,t, xi,dxidt)
c
c  .......determine discontinuous Hcurl shape functions
          call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
c
c  .......determine element H1 shape functions (for geometry)
          call shape3H(etype,xi,norder,norient_edge,norient_face,
     .                 nrdofH,shapH,gradH)
c
c  .......determine element H(curl) shape functions (for fluxes)
          call shape3E(etype,xi,norderc,norient_edge,norient_face,
     .                 nrdofE,shapE,curlE)
c
c  .......geometry
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,
     .                 x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
          weight = bjac*wtloc(l)
c  ......check for impedance BC
c  .......impedance BC boundary ......................................
          if (ibc(if,2).eq.9) then
            !write(*,*) 'putting impedance BCs'
c
c  .........get the boundary source
            call get_bdSource(Mdle,x,rn,IBCFlag, zImp)
c
c  .........loop through Hcurl enriched test functions
            do k1=1,nrdofEE
c
c  ...........value of the shape function at the point
              qq(1:3) = shapEE(1,k1)*dxidx(1,1:3)
     .                + shapEE(2,k1)*dxidx(2,1:3)
     .                + shapEE(3,k1)*dxidx(3,1:3)
c
c  ...........accumulate for the load vector
              k = 2*k1
              BLOADE(k) = BLOADE(k)
     .      - (zImp(1)*qq(1)+zImp(2)*qq(2)+zImp(3)*qq(3))*weight
c
c  ...........loop through Hcurl trial functions
              do k2=1,nrdofE
c
c  .............value of the shape function at the point
                p(1:3) = shapE(1,k2)*dxidx(1,1:3)
     .                 + shapE(2,k2)*dxidx(2,1:3)
     .                 + shapE(3,k2)*dxidx(3,1:3)

                call cross_product(rn,p, rntimesp)
                call cross_product(rn,rntimesp, rn2timesp)
c
c  .............accumulate for the extended stiffness matrix
                STIFFEE(2*k1,2*k2-1) = STIFFEE(2*k1,2*k2-1)
     .              + (qq(1)*rn2timesp(1)+ qq(2)*rn2timesp(2)+
     .       + qq(3)*rn2timesp(3))*(1.d0/impedanceConstant)*weight

                STIFFEE(2*k1-1,2*k2-1) = STIFFEE(2*k1-1,2*k2-1)
     .       +  (qq(1)*rntimesp(1)+ qq(2)*rntimesp(2)+
     .       +   qq(3)*rntimesp(3))*weight
              enddo
            enddo

c  .......regular boundary.............................................
          else
c  ....... check for LASER MODE, if
c  ....... LASER MODE is active, get boundary data for source
c
          if(LASER_MODE.eq.1) then
c  .........get the boundary source
            call get_bdSource(Mdle,x,rn,0, zImp)
          else
            zImp = ZERO
          end if
c  .......loop through the enriched H(curl) test functions
          do k1=1,nrdofEE
            E1(1:3) = shapEE(1,k1)*dxidx(1,1:3)
     .              + shapEE(2,k1)*dxidx(2,1:3)
     .              + shapEE(3,k1)*dxidx(3,1:3)
c
c  ...........accumulate for the load vector
              k = 2*k1
              BLOADE(k) = BLOADE(k)
     .      + (zImp(1)*E1(1)+zImp(2)*E1(2)+zImp(3)*E1(3))*weight
c  .........loop through H(curl) trial functions
            do k2=1,nrdofE
              E2(1:3) = shapE(1,k2)*dxidx(1,1:3)
     .                + shapE(2,k2)*dxidx(2,1:3)
     .                + shapE(3,k2)*dxidx(3,1:3)
              call cross_product(rn,E2, rntimesE)
c
c  ...........accumulate for the extended stiffness matrix
              STIFFEE(2*k1-1,2*k2-1) = STIFFEE(2*k1-1,2*k2-1)
     .       + (E1(1)*rntimesE(1)+E1(2)*rntimesE(2)+E1(3)*rntimesE(3))
     .          *weight
              STIFFEE(2*k1  ,2*k2  ) = STIFFEE(2*k1  ,2*k2  )
     .       + (E1(1)*rntimesE(1)+E1(2)*rntimesE(2)+E1(3)*rntimesE(3))
     .          *weight
            enddo
          enddo
c .......end if for impedance BC
        endif
c
        enddo
      enddo
c
      if (iprint.eq.5) then
        write(*,*) 'elem: STIFFEE, STIFFEQ, BLOADE = '
        do j=1,2*nrdofEE
          write(*,7200) j, BLOADE(j)
          write(*,*) 'STIFFEE for n x E = '
          write(*,7210) (STIFFEE(j,2*k-1),k=1,nrdofE)
          write(*,*) 'STIFFEE for n x H = '
          write(*,7210) (STIFFEE(j,2*k  ),k=1,nrdofE)
          write(*,*) 'STIFFEQ = '
          write(*,7210) STIFFEQ(j,1:6*nrdofQ)
 7200     format('j = ',i5,' BLOAD(j) = ',2e12.5)
 7210     format(6(2e12.5,1x))
          write(*,*) 'residual = '
          if (j/2*2.eq.j) then
            write(*,7210) BLOADE(j) - STIFFEQ(j,1)
          else
            write(*,7210) BLOADE(j) - STIFFEQ(j,1)
     .      - (STIFFEE(j,1)+STIFFEE(j,5)+STIFFEE(j,9)+STIFFEE(j,13))
          endif
          call pause
        enddo
      endif

c! ....................................................................
c! ...Alternative construction of normal matrix
c! ! ....................................................................
c
c!  ...total trial dof for the element
      nrTEST = nrdofEE*2
      call celndof(NODES(Mdle)%type,norder,
     .             nrdofH,nrdofE,nrdofV,nrdofQ)

      i1 = MAXtestE ; j1 = 2*nrdofE ; j2 = 6*nrdofQ

      STIFF_ALLE(1:i1,1:j1) = STIFFEE(1:i1,1:j1)
      STIFF_ALLE(1:i1,j1+1:j1+j2) = STIFFEQ(1:i1,1:j2)
      STIFF_ALLE(1:i1,j1+j2+1) = BLOADE(1:i1)

c
      uplo  = 'U' ; trans = 'C' ;   diag  = 'N'
      N     = nrTEST
      NRHS  = 2*nrdofE + 6*nrdofQ + 1
c
c-----------------------------------------------------------------------
c----------------------------------------------------------------------

      allocate(DIAG_E(nrTest))
c  ...preconditioning: AP_Maxwell = D^-1/2 * AP_Maxwell * D^-1/2
      do k1=1,nrTEST
          k = nk(k1,k1)
          DIAG_E(k1) = AP_Maxwell(k)
      enddo
      do k1=1,nrTEST
         do k2=k1,nrTEST
            k = nk(k1,k2)
            AP_Maxwell(k) = AP_Maxwell(k)/sqrt(DIAG_E(k1)*DIAG_E(k2))
         enddo
      enddo

c  ...preconditioning: STIFF_ALLE = D^-1/2 * STIFF_ALLE
      do k2=1,NRHS
        do k1=1,nrTEST
          STIFF_ALLE(k1,k2) = STIFF_ALLE(k1,k2)/sqrt(DIAG_E(k1))
        enddo
      enddo
      deallocate(DIAG_E)

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c.....check condition number
C       nrdof_eig = nrdofEE*2
C       kk = nrdof_eig*(nrdof_eig+1)/2
C       allocate(AP_eig(kk))
C       AP_eig(1:kk) = AP_Maxwell(1:kk)
C       allocate(W(nrdof_eig))
C       allocate(Z(1,nrdof_eig))
C       allocate(WORK(nrdof_eig))
C       allocate(RWORK(nrdof_eig))
C       allocate(IWORK(1))
C       call ZHPEVD('N','U',nrdof_eig,
C      .             AP_eig, W,Z,1,WORK,nrdof_eig,
C      .             RWORK,nrdof_eig,IWORK,1,info)
C       if (info .ne. 0) then
C          write(*,*) 'eig_solve_sc: info = ', info
C          stop 1
C       endif

C       write(*,6999) W(nrdof_eig),W(1)
C  6999 format('elem_dpgMaxwell: AP_Maxwell: max_eig, min_eig = ', 2e13.4)


C       write(*,7000)  W(nrdof_eig)/W(1)
C  7000 format('elem_dpgMaxwell: AP_Maxwell condition number = ', 1e13.4)
C       deallocate(IWORK,W,WORK)
C       deallocate(RWORK,Z)
C       deallocate(AP_eig)
C       call pause
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------

c  ...factorize the test stiffness matrix
      uplo = 'U'
      call ZPPTRF(uplo, nrTEST, AP_Maxwell, info)
      if (info.ne.0) then
        write(*,*) 'elem_dpgMaxwell: info = ',info
        stop1
      endif


      call ZTPTRS(uplo,trans,diag,N,NRHS,AP_Maxwell,STIFF_ALLE,
     .             MAXtestE,info)
c
      call ZHERK(uplo,trans,NRHS,nrTEST,ZONE,STIFF_ALLE,
     .            MAXtestE,ZERO,
     .            STIFF_ALLE(1:NRHS,1:NRHS),NRHS)

c     ZHERK for complex case

      do i=1,NRHS-1
        STIFF_ALLE(i+1:NRHS,i) = conjg(STIFF_ALLE(i,i+1:NRHS))
      enddo
c
c
      ZblocE(1:j1) = STIFF_ALLE(1:j1,j1+j2+1)
      ZblocQ(1:j2) = STIFF_ALLE(j1+1:j1+j2,j1+j2+1)
c
      ZalocEE(1:j1,1:j1) = STIFF_ALLE(1:j1,1:j1)
      ZalocEQ(1:j1,1:j2) = STIFF_ALLE(1:j1,j1+1:j1+j2)
c
      ZalocQE(1:j2,1:j1) = STIFF_ALLE(j1+1:j1+j2,1:j1)
      ZalocQQ(1:j2,1:j2) = STIFF_ALLE(j1+1:j1+j2,j1+1:j1+j2)

      if (iprint.ge.1) then
        write(*,7010)
 7010   format('elem_dpgMaxwell: ZblocE,ZblocQ = ')
        write(*,7011) ZblocE(1:2*NrdofE)
        write(*,7011) ZblocQ(1:6*NrdofQ)
 7011   format(10e12.5)
        call pause
        write(*,7012)
 7012   format('elem_dpgMaxwell: ZalocEE = ')
        do i=1,2*NrdofE
          write(*,7013) i,ZalocEE(i,1:2*NrdofE)
 7013     format('i = ',i3,10(/,5(2e12.5,2x)))
        enddo
        call pause
        write(*,7014)
 7014   format('elem_dpgMaxwell: ZalocEQ = ')
        do i=1,2*NrdofE
          write(*,7013) i,ZalocEQ(i,1:6*NrdofQ)
        enddo
        call pause
        write(*,7015)
 7015   format('elem_dpgMaxwell: ZalocQQ = ')
        do i=1,6*NrdofQ
          write(*,7013) i,ZalocQQ(i,1:6*NrdofQ)
        enddo
        call pause
      endif
c
      end subroutine elem_dpgMaxwell


      subroutine get_gainFunction(EfieldNorm, gainFunction)
      use DPGLaser
      use problem
      real*8 alpha,f_em,f_abs,wEratio,Nratio,Nex,Ngd
      real*8,                       intent(in)  :: EfieldNorm
      real*8,                       intent(out)  :: gainFunction
      sigma_s_abs_scaled = SIGMA_S_ABS*N0*CORE_RAD
      sigma_s_em_scaled = SIGMA_S_EM*N0*CORE_RAD
      sigma_p_abs_scaled = SIGMA_P_ABS*N0*CORE_RAD
      sigma_p_em_scaled = SIGMA_P_EM*N0*CORE_RAD

      alpha = (TAU*REF_INDEX_CORE*E_0**2)/(MU0*LIGHT_SPEED*CORE_RAD*
     .                                      N0*H_BAR*OMEGA0);
      f_em = sigma_s_em_scaled+((1.d0/OMEGA_RATIO)*sigma_p_em_scaled);
      f_abs = sigma_s_abs_scaled+((1.d0/OMEGA_RATIO)
     .                             *sigma_p_abs_scaled);
      wEratio = OMEGA/(alpha*EfieldNorm**2);
      Nratio= (wEratio+f_em)/f_abs;
      Nex = 1.0/(1.0+Nratio);
      Ngd = 1.0-Nex;
      gainFunction = sigma_s_em_scaled *Nex - sigma_s_abs_scaled*Ngd;
c      write(*,*) 'from get_gainFunction: EfieldNorm = '
c     .             EfieldNorm
      end subroutine get_gainFunction







