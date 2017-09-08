!--------------------------------------------------------------------
!
!     routine name      - elem
!
!--------------------------------------------------------------------
!
!     latest revision:  - April 17
!
!     purpose:          - routine returns unconstrained (ordinary)
!                         stiffness matrix and load vector
!                         for either:
!                         1) single step of transient heat problem
!                         2) multiple steps of transient heat problem
!                         3) steady state Maxwell problem
!
!     arguments:
!    in:
!            Mdle      - an element middle node number, identified
!                        with the element
!
!    out:    Itest     - index for assembly
!            Itrial    - index for assembly
!
!
!-----------------------------------------------------------------------
!
subroutine elem(Mdle, Itest,Itrial)
!
!......modules used...............
  use control
  use data_structure3D
  use assembly
  use problem
!......no implicit statements.....
  implicit none
!
!......declare input/output variables..........
!
  integer,                    intent(in)  :: Mdle
  integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!
  Itest(1:NR_PHYSA)=0 ; Itrial(1:NR_PHYSA)=0
!
  select case(NODES(Mdle)%case)
!
!    node supports all physical attributes
!    4 physical attributes: case = 2^4-1 = 15
    case(15)
!
    select case(NO_PROBLEM)
! ..... NO_PROBLEM: (1) - single step of transient heat equation
!...................(2) - multiple steps of transient heat equation
!...................(3) - time harmonic Maxwell
!.........Heat cases
      case(1,2)
      Itest(1)=1 ; Itest(3)=1
      Itrial(1)=1; Itrial(3)=1
      call elem_dpgHeat(Mdle, &
        BLOC(1)%nrow,BLOC(3)%nrow, &
        BLOC(1)%array,ALOC(1,1)%array,ALOC(1,3)%array, &
        BLOC(3)%array,ALOC(3,1)%array,ALOC(3,3)%array)
!.........Maxwell case
      case(3)
      Itest(2)=1 ; Itest(4)=1;
      Itrial(2)=1; Itrial(4)=1;
      call elem_dpgMaxwell(Mdle, &
        BLOC(2)%nrow,BLOC(4)%nrow, &
        BLOC(2)%array,ALOC(2,2)%array,ALOC(2,4)%array, &
        BLOC(4)%array,ALOC(4,2)%array,ALOC(4,4)%array)
!
!.... end select for select case(NO_PROBLEM)
    end select
!
    case default
      write(*,*) 'elem: Mdle,NODES(Mdle)%case = ', &
        Mdle,NODES(Mdle)%case
      call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
      stop
!
!.... end select for select case(NODES(Mdle)%case)
  endselect
!
!
end subroutine elem
!
!--------------------------------------------------------------------
!
!    routine name      - elem_dpgHeat
!
!--------------------------------------------------------------------
!
!    latest revision:  - April 17
!
!    purpose:          - routine returns unconstrained (ordinary)
!                        stiffness matrix and load vector
!                        for the primal H1 formulation for a single
!                        step of transient heat equation
!
!    arguments:
!
!    in:
!            Mdle      - an element middle node number, identified
!                        with the element
!            MdH       - column length of ZalocHH,ZalocHV
!            MdV       - column length of ZalocVH,ZalocVV
!    out:
!            ZblocH,ZblocV - load vectors
!            ZalocHH,ZalocHV,ZalocVH,ZalocVV - stiffness matrices
!
!---------------------------------------------------------------------
!
subroutine elem_dpgHeat(Mdle,MdH,MdV, &
  ZblocH,ZalocHH,ZalocHV,ZblocV,ZalocVH,ZalocVV)
!.......modules used
  use control
  use parametersDPG
  use element_data
  use data_structure3D
!  use DPGLaser
  use problem
!.......no implicit statements
  implicit none
#if C_MODE
#define VTYPE  complex*16
#else
#define VTYPE double precision
#endif
!.......declare input/output variables
  integer,                     intent(in)  :: Mdle
  integer,                     intent(in)  :: MdH
  integer,                     intent(in)  :: MdV
  VTYPE, dimension(MdH),       intent(out) :: ZblocH
  VTYPE, dimension(MdH,MdH),   intent(out) :: ZalocHH
  VTYPE, dimension(MdH,MdV),   intent(out) :: ZalocHV
  VTYPE, dimension(MdV),       intent(out) :: ZblocV
  VTYPE, dimension(MdV,MdH),   intent(out) :: ZalocVH
  VTYPE, dimension(MdV,MdV),   intent(out) :: ZalocVV
!
!.......declare edge/face type varibles
  character(len=4) :: etype,ftype
!
! ...declare element order, orientation for edges and faces
  integer, dimension(19)    :: norder
  integer, dimension(12)    :: norient_edge
  integer, dimension(6)     :: norient_face
! ...face order
  integer, dimension(5)     :: norderf
!
! ...geometry dof (work space for nodcor)
  real*8, dimension(3,MAXbrickH) :: xnod
!
! ...solution dof (work space for solelm)
  VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
  VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
  VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
  VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
! .... current L2 solution zsolQ
  VTYPE, dimension(6) :: zsolQ
!
! ... variables for geometry
  real*8, dimension(3)      :: xi,x,rn
  real*8, dimension(3,2)    :: dxidt,dxdt,rt
  real*8, dimension(3,3)    :: dxdxi,dxidx
  real*8, dimension(2)      :: t
!
! ...H1 shape functions
  real*8, dimension(MAXbrickH)    :: shapH
  real*8, dimension(3,MAXbrickH)  :: gradH
! ...H(curl) shape functions
  real*8, dimension(3,MAXbrickE)  :: shapE
  real*8, dimension(3,MAXbrickE)  :: curlE
! ....H(div) shape functions
  real*8, dimension(3,MAXbrickV)  :: shapV
  real*8, dimension(MAXbrickV)    :: divV
! ....L2 shape functions
  real*8, dimension(MAXbrickQ)    :: shapQ
! .... Enriched H1 shape functions
  real*8 , dimension(MAXbrickHH)    :: shapHH
  real*8 , dimension(3,MAXbrickHH)  :: gradHH
! ... nrdof for various spaces
  integer  :: nrdofH,nrdofE,nrdofV,nrdofQ,nrdofHH
!
! ....... space for DPG Computations (Gram Matrix, Stiffness etc.)
  integer, parameter   :: MAXtestH = MAXbrickHH
!  ...stiffness matrix for the local Riesz H1 matrix in Lapack format
  VTYPE, dimension(MAXtestH*(MAXtestH+1)/2) :: AP_Heat
!  ...load vector for the enriched space
  VTYPE, dimension(MAXtestH) :: BLOADH
!  ...stiffness matrices for the enriched test space
  VTYPE, dimension(MAXtestH,MAXbrickH) :: STIFFHH
  VTYPE, dimension(MAXtestH,MAXbrickV) :: STIFFHV
!  ....STIFF_ALL for alternative computation of stiffness
  VTYPE, dimension(MAXtestH,MAXbrickH+MAXbrickV+1) :: STIFF_ALLH
#if C_MODE
  complex*16, allocatable :: AP_eig(:)
  complex*16, allocatable :: DIAG_E(:)
  complex*16, allocatable :: DIAG_H(:)
#else
  real*8, allocatable     :: AP_eig(:)
  real*8, allocatable     :: DIAG_E(:)
  real*8, allocatable     :: DIAG_H(:)
#endif
! ..... dummy for elem_residual
  VTYPE, dimension(MAXtestH*(MAXtestH+1)/2) :: AP
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
! ...derivatives wrt physical coordinates, flux
  real*8, dimension(3)  :: dv1,dv2,vec
!
! ...for debug printing
  VTYPE, dimension(10)  :: zaux

! ....Maxwell load
  VTYPE, dimension(3) :: zJ
!
  integer,save :: ivis=0
  character    :: uplo, trans,diag
!
! .... number of vertices,edge,faces per element type
  integer :: nrv, nre, nrf
! .... for Gram matrix
  integer      :: nk
! ..... various variables for the problem
  real*8  :: h_elem,rjac,weight,wa,v2n,v1,v2,bjac,gainFunction,EfieldNorm
  integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,nrTEST,nint,iflag,kE,k,iprint
  integer :: N,nRHS,nordP,l,nsign,if,ndom,info
  VTYPE   :: zfval,zsolH
  nk(k1,k2) = (k2-1)*k2/2+k1

!
!------------------INITIALIZE THE ELEMENT ORDER, ORIENTATION ETC.------
!
  iprint=0
!
  ivis=1
  if (ivis.eq.0) then
    write(*,*) 'elem_dpgHeat: NO_PROBLEM    = ',NO_PROBLEM
    write(*,*) 'elem_dpgHeat: INNER_PRODUCT = ',INNER_PRODUCT
    call pause
    ivis=1
  endif
!
! ...element type
  etype = NODES(Mdle)%type
  nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
! ...determine order of approximation
  call find_order(Mdle, norder)
!
! ...set the enriched order of appoximation depending on element type
  select case(etype)
    case('mdlb')        ; nordP = NODES(Mdle)%order + NORD_ADD*111
    case('mdln','mdld') ; nordP = NODES(Mdle)%order + NORD_ADD*1
    case('mdlp')        ; nordP = NODES(Mdle)%order + NORD_ADD*11
  end select
!
! ...determine edge and face orientations
  call find_orient(Mdle, norient_edge,norient_face)
!
! ...determine nodes coordinates
  call nodcor(Mdle, xnod)
!
! ...determine dof for current solution if running multiple
! ... steps of heat equation
  select case(NO_PROBLEM)
  case(2)
    call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
  end select
!
! ...determine element size and scale correction
  h_elem = min(abs(xnod(1,2)-xnod(1,1)),abs(xnod(2,4)-xnod(2,1)), &
           abs(xnod(3,5)-xnod(3,1)))
!-----------
!-----------  CHECK OUT INNER PRODUCT SCALING FOR HEAT STEPPING
!-----------
!     select case(INNER_PRODUCT)
!     case(2)
!       scale = min(1.d0,EPS_DIFF**3/h_elem**2)
!     case(3)
!       scale =  min(1.d0,EPS_DIFF**2/h_elem**2)
!     end select
!     scale = 1.d0
!-----------
!-----------  CHECK OUT INNER PRODUCT SCALING FOR HEAT STEPPING
!-----------
!
! ...get the element boundary conditions flags
  call find_bc(Mdle, ibc)
  if (iprint.ge.1) then
    write(*,7001) Mdle
7001   format('elem_dpgHeat: B!FLAGS FOR Mdle = ',i5)
    do i=1,NR_PHYSA
      write(*,7002) PHYSA(i), ibc(1:nrf,i)
7002     format('          ATTRIBUTE = ',a6,' FLAGS = ',6i2)
    enddo
  endif
!
! ...clear space for stiffness matrix and rhsv:
  ZblocH = ZERO; ZblocV = ZERO
  ZalocHH = ZERO; ZalocHV = ZERO; ZalocVH = ZERO; ZalocVV = ZERO
!
!    extended load vector and extended stiffness matrices
  BLOADH=ZERO ; STIFFHH=ZERO ; STIFFHV=ZERO; STIFF_ALLH=ZERO
!
!    Gram matrix
  AP_Heat=ZERO
!
!
!-----------------------------------------------------------------------
!    E L E M E N T    I N T E G R A L S                               |
!-----------------------------------------------------------------------
!
! ...use the enriched order to set the quadrature
  INTEGRATION = NORD_ADD
    call set_3Dint(etype,norder, nint,xiloc,waloc)
  INTEGRATION = 0
!
! ...loop over integration points
  do l=1,nint
    xi(1:3)=xiloc(1:3,l) ; wa=waloc(l)
!
! .....H1 shape functions
    call shape3H(etype,xi,norder,norient_edge,norient_face, &
         nrdofH,shapH,gradH)
!
!....... element L2 solution for heat Load calculation
! .....determine element L2 shape functions in LASER_MODE
    if(LASER_MODE.eq.1) then
      call shape3Q(etype,xi,norder, nrdofQ,shapQ)
    endif
!
! .....discontinuous H1 shape functions
    call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!   write(*,*) 'nrdofHH is: ', nrdofHH
!   call pause
!
! .....geometry computations
    call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, &
                x,dxdxi,dxidx,rjac,iflag)
!
!..................................................................
!..................................................................
!.........COMPUTE H1 and Hcurl Solutions for Heat stepping and
!......... Load for Heat equation
!......... heat load = |E|^2
!..................................................................
!..................................................................
! .....compute the H1 solution at the point
    zsolH = ZERO
    do kH=1,nrdofH
      zsolH = zsolH + zdofH(1,kH)*shapH(kH)
    enddo
! .....compute the Hcurl solution at the point if LASER_MOD
    !if((LASER_MODE.eq.1).or.(NONLINEAR_FLAG.eq.1))then
    if(LASER_MODE.eq.1) then
      zsolQ = ZERO
        do kE=1,nrdofQ
          zsolQ(1:6) = zsolQ(1:6) + zdofQ(1:6,kE)*shapQ(kE)
        enddo
    endif
    EfieldNorm = abs(zsolQ(1))**2+abs(zsolQ(2))**2+abs(zsolQ(3))**2
    EfieldNorm = dsqrt(EfieldNorm)
    call get_gainFunction(EfieldNorm, gainFunction)

!
! .....integration weight
    weight = rjac*wa
!
! .....get the RHS
    call getf(Mdle,x, zfval,zJ)
!
! .....1st loop through enriched H1 test functions
    do k1=1,nrdofHH
      v1 = shapHH(k1)
!
      dv1(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
              + gradHH(2,k1)*dxidx(2,1:3) &
              + gradHH(3,k1)*dxidx(3,1:3)
!
!        -- EXTENDED LOAD VECTOR --
!
      select case(NO_PROBLEM)
        case(1)
          BLOADH(k1) = BLOADH(k1) + zfval*v1*weight
        case(2)
!          if((LASER_MODE.eq.1).or.(NONLINEAR_FLAG.eq.1)) then
          if(LASER_MODE.eq.1) then
            zfval = DELTAT*((zsolQ(1)**2+zsolQ(2)**2 &
                    +zsolQ(3)**2)*REF_INDEX_CORE*gainFunction)
          endif
          BLOADH(k1) = BLOADH(k1) + (zfval+zsolH)*v1*weight
      end select
!
! .......2nd loop through enriched H1 trial functions
! ....... for Gram matrix
      do k2=k1,nrdofHH
!
        v2 = shapHH(k2)
!
        dv2(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                 + gradHH(2,k2)*dxidx(2,1:3) &
                 + gradHH(3,k2)*dxidx(3,1:3)
!
! ----------- GRAM MATRIX --------------
!
!          determine index in triagular format
        k = nk(k1,k2)
!
! .........problem-dependent inner product
        select case(INNER_PRODUCT)
          case(1)
            AP_Heat(k) = AP_Heat(k) &
                   + (dv1(1)*dv2(1)+dv1(2)*dv2(2)+dv1(3)*dv2(3) &
                   +v1*v2)*weight
          case default
            stop
        end select
!
! .......enddo 2nd loop through enriched H1 trial functions
      enddo
!
!........loop through H1 trial functions
! ........ for stiffness matrix
      do k2=1,nrdofH
        v2 = shapH(k2)
        dv2(1:3) = gradH(1,k2)*dxidx(1,1:3) &
                  + gradH(2,k2)*dxidx(2,1:3) &
                  + gradH(3,k2)*dxidx(3,1:3)
!
!          -- EXTENDED HH STIFFNESS MATRIX --
        STIFFHH(k1,k2) = STIFFHH(k1,k2) &
           +(KAPPA*DELTAT*(dv1(1)*dv2(1)+dv1(2)*dv2(2)+dv1(3)*dv2(3)) &
            +  v1*v2 )*weight
!........end loop through H1 trial functions
      enddo
!........end 1st loop through enriched H1 test functions
    enddo
!........end loop through integration points
  enddo
!
! ...printing Gram matrix
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
!
!
!-----------------------------------------------------------------------
!    B O U N D A R Y    I N T E G R A L S                             |
!-----------------------------------------------------------------------
!
! ...loop through element faces
  do if=1,nrf
!
! .....sign factor to determine the OUTWARD normal unit vector
    nsign = nsign_param(etype,if)
!
! .....face type
    ftype = face_type(etype,if)
!
! .....face order of approximation
    call face_order(etype,if,norder, norderf)
!
! .....set 2D quadrature
    INTEGRATION = NORD_ADD
      call set_2Dint(ftype,norderf, nint,tloc,wtloc)
    INTEGRATION = 0
!
! .....loop through integration points
    do l=1,nint
!
! .......face coordinates
      t(1:2) = tloc(1:2,l)
!
! .......face parametrization
      call face_param(etype,if,t, xi,dxidt)
!
! .......determine discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
! .......determine element H1 shape functions (for geometry)
      call shape3H(etype,xi,norder,norient_edge,norient_face, &
                  nrdofH,shapH,gradH)
!
! .......determine element Hdiv shape functions (for fluxes)
      call shape3V(etype,xi,norder,norient_face, &
                  nrdofV,shapV,divV)
!
! .......geometry
      call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                  x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
      weight = bjac*wtloc(l)
!
! .......loop through enriched H1 test functions
      do k1=1,nrdofHH
        v1 = shapHH(k1)
!
! .........loop through H(div) trial functions
        do k2=1,nrdofV
!
! ...........normal component (Piola transformation at work!)
          vec(1:3) =  dxdxi(1:3,1)*shapV(1,k2) &
                    +dxdxi(1:3,2)*shapV(2,k2) &
                    +dxdxi(1:3,3)*shapV(3,k2)
          vec(1:3) = vec(1:3)/rjac
          v2n = vec(1)*rn(1)+vec(2)*rn(2)+vec(3)*rn(3)
!
!            -- EXTENDED HV STIFFNESS MATRIX --
!
          STIFFHV(k1,k2) = STIFFHV(k1,k2) &
                          - v1*v2n*weight
! ..... end loop through H(div) trial functions
        enddo
! ..... end loop through enriched H1 test functions
      enddo
! ..... end loop through integration points
    enddo
! ..... end loop through element faces
  enddo
!
!-----------------------------------------------------------------------
!
! ...factorize the Gram matrix
  uplo = 'U'
  call ZPPTRF(uplo, nrdofHH, AP_Heat, info)
  if (info.ne.0) then
    write(*,*) 'elem_dpgHeat: info = ',info ; stop
  endif

! ....................................................................
!  construction of DPG system
! ....................................................................
!
!  ...total trial dof for the element
  nrTEST = nrdofHH
  call celndof(NODES(Mdle)%type,norder, &
              nrdofH,nrdofE,nrdofV,nrdofQ)

  i1 = MAXtestH ; j1 = nrdofH ; j2 = nrdofV
!
  STIFF_ALLH(1:i1,1:j1) = STIFFHH(1:i1,1:j1)
  STIFF_ALLH(1:i1,j1+1:j1+j2) = STIFFHV(1:i1,1:j2)
  STIFF_ALLH(1:i1,j1+j2+1) = BLOADH(1:i1)

!
  uplo  = 'U' ; trans = 'C' ;   diag  = 'N'
  N     = nrTEST
  nRHS  = nrdofH + nrdofV + 1
!
!
  call ZTPTRS(uplo,trans,diag,N,NRHS,AP_heat,STIFF_ALLH, &
              MAXtestH,info)
!
!
  call ZHERK(uplo,trans,NRHS,nrTEST,ZONE,STIFF_ALLH, &
              MAXtestH,ZERO, &
              STIFF_ALLH(1:NRHS,1:NRHS),NRHS)
!    ZHERK for complex case
  do i=1,NRHS-1
    STIFF_ALLH(i+1:NRHS,i) = conjg(STIFF_ALLH(i,i+1:NRHS))
  enddo
!
  ZblocH(1:j1) = STIFF_ALLH(1:j1,j1+j2+1)
  ZblocV(1:j2) = STIFF_ALLH(j1+1:j1+j2,j1+j2+1)
!
  ZalocHH(1:j1,1:j1) = STIFF_ALLH(1:j1,1:j1)
  ZalocHV(1:j1,1:j2) = STIFF_ALLH(1:j1,j1+1:j1+j2)
!
  ZalocVH(1:j2,1:j1) = STIFF_ALLH(j1+1:j1+j2,1:j1)
  ZalocVV(1:j2,1:j2) = STIFF_ALLH(j1+1:j1+j2,j1+1:j1+j2)
!
!
#if C_MODE
 7011   format(6(2e10.3,2x))
#else
 7011   format(10e12.5)
#endif
!
end subroutine elem_dpgHeat

!-------------------------------------------------------------------
!
!    routine name      - elem_dpgMaxwell
!
!--------------------------------------------------------------------
!
!    latest revision:  - April 17
!
!    purpose:          - routine returns unconstrained (ordinary)
!                        stiffness matrix and load vector
!                        for the UW DPG formulation for Maxwell
!                        equations
!
!    arguments:
!
!    in:
!            Mdle      - an element middle node number, identified
!                        with the element
!            MdE       - column length of ZalocEE,ZalocEQ
!            MdQ       - column length of ZalocQE,ZalocQQ
!    out:
!            ZblocE,ZblocQ - load vectors
!            ZalocEE,ZalocEQ,ZalocQE,ZalocQQ - stiffness matrices
!
!---------------------------------------------------------------------
!
subroutine elem_dpgMaxwell(Mdle,MdE,MdQ, &
                              ZblocE,ZalocEE,ZalocEQ, &
                              ZblocQ,ZalocQE,ZalocQQ)
! ... modules used ....
  use control
  use parametersDPG
  use element_data
  use data_structure3D
  use matrices
  use problem
! ... no implicit statements
  implicit none
!
#if C_MODE
#define VTYPE  complex*16
#else
#define VTYPE double precision
#endif
!.......declare input/output variables
  integer,                     intent(in)  :: Mdle
  integer,                     intent(in)  :: MdE
  integer,                     intent(in)  :: MdQ
  VTYPE, dimension(MdE),       intent(out) :: ZblocE
  VTYPE, dimension(MdE,MdE),   intent(out) :: ZalocEE
  VTYPE, dimension(MdE,MdQ),   intent(out) :: ZalocEQ
  VTYPE, dimension(MdQ),       intent(out) :: ZblocQ
  VTYPE, dimension(MdQ,MdE),   intent(out) :: ZalocQE
  VTYPE, dimension(MdQ,MdQ),   intent(out) :: ZalocQQ
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
! ...solution dof (work space for solelm)
  VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
  VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
  VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
  VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
! .... current L2 solution zsolQ
  VTYPE, dimension(6) :: zsolQ
! .... current H1 solution
  VTYPE               :: zsolH
!
! ... variables for geometry
  real*8, dimension(3)      :: xi,x,rn
  real*8, dimension(3,2)    :: dxidt,dxdt,rt
  real*8, dimension(3,3)    :: dxdxi,dxidx
  real*8, dimension(2)      :: t
!
! ...H1 shape functions
  real*8, dimension(MAXbrickH)    :: shapH
  real*8, dimension(3,MAXbrickH)  :: gradH
! ...H(curl) shape functions
  real*8, dimension(3,MAXbrickE)  :: shapE
  real*8, dimension(3,MAXbrickE)  :: curlE
! ....H(div) shape functions
  real*8, dimension(3,MAXbrickV)  :: shapV
  real*8, dimension(MAXbrickV)    :: divV
! ....L2 shape functions
  real*8, dimension(MAXbrickQ)    :: shapQ
! .... Enriched H1 shape functions
  real*8 , dimension(3,MAXbrickEE)    :: shapEE
  real*8 , dimension(3,MAXbrickEE)    :: curlEE
! ... nrdof for various spaces
  integer  :: nrdofH,nrdofE,nrdofV,nrdofQ,nrdofEE
! ....... space for DPG Computations (Gram Matrix, Stiffness etc.)
  integer, parameter  :: MAXtestE = 2*MAXbrickEE
!  ...stiffness matrix for the local Riesz H1 matrix in Lapack format
  VTYPE, dimension(MAXtestE*(MAXtestE+1)/2) :: AP_Maxwell
!  ...load vector for the enriched space
  VTYPE, dimension(MAXtestE) :: BLOADE
!  ...stiffnes matrices for the enriched test space
  VTYPE, dimension(MAXtestE,MAXbrickQ*6) :: STIFFEQ
  VTYPE, dimension(MAXtestE,MAXbrickE*2) :: STIFFEE
!  ....STIFF_ALL for alternative computation of stiffness
  VTYPE, dimension(MAXtestE,2*MAXbrickE+6*MAXbrickQ+1) :: STIFF_ALLE
#if C_MODE
  complex*16, allocatable :: AP_eig(:)
  complex*16, allocatable :: DIAG_E(:)
  complex*16, allocatable :: DIAG_H(:)
#else
  real*8, allocatable     :: AP_eig(:)
  real*8, allocatable     :: DIAG_E(:)
  real*8, allocatable     :: DIAG_H(:)
#endif
! ..... dummy for elem_residual
  VTYPE, dimension(MAXtestE*(MAXtestE+1)/2) :: AP
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
! ...derivatives wrt physical coordinates, flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  real*8, dimension(3)  :: dv1,dv2,vec
!
! ...for debug printing
  VTYPE, dimension(10)  :: zaux
  !VTYPE, dimension(24)  :: zE
  !VTYPE, dimension(6)   :: zQ

! ....Maxwell load and auxiliary variables
  VTYPE, dimension(3) :: zJ,zImp
  real*8, dimension(3):: qq,p,rntimesp,rn2timesp
  real*8, dimension(3) :: E1,curlE1,E2,curlE2,rntimesE
!
  integer,save :: ivis=0
  character    :: uplo, trans,diag
!
! .... number of vertices,edge,faces per element type
  integer :: nrv, nre, nrf
! .... for Gram matrix
  integer      :: nk
! ..... various variables for the problem
  real*8  :: h_elem,rjac,weight,wa,v2n,CC,EE,CE,E,EC,q,h,omeg,alpha_scale
  real*8  :: bjac,impedanceConstant
  integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,nrTEST,nint,iflag,kE,k,iprint,l
  integer :: N,nRHS,nordP,nsign,if,ndom,info,icomp,nrdof_eig,idec
  VTYPE   :: zfval,za,zb,zc
!  ... for constraints
  VTYPE, allocatable  :: zConstraintMatrix(:,:)
  VTYPE, allocatable  ::  zConstraintLoad(:)
! ...for lapack eigensolve
  complex*16, allocatable :: Z(:,:), WORK(:)
  real*8, allocatable     :: W(:),   RWORK(:)
  integer, allocatable    :: IWORK(:)
! ... for gain function
  real*8  :: bg_gain, ion_gain,EfieldNorm,gainFunction,rndotE,alpha

  nk(k1,k2) = (k2-1)*k2/2+k1
! ..... allocate copy matrices
!
!---------------------------------------------------------------------
!
  iprint=0
  if (iprint.eq.1) then
    write(*,*) 'elem: Mdle = ',Mdle
  endif
!
! ...element type
  etype = NODES(Mdle)%type
  nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
! ...determine order of approximation
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
! ...determine edge and face orientations
  call find_orient(Mdle, norient_edge,norient_face)
!
! ...determine nodes coordinates
  call nodcor(Mdle, xnod)

! ...get the element boundary conditions flags
  call find_bc(Mdle, ibc)
  iprint = 0
  if (iprint.ge.1) then
    write(*,7001) Mdle
7001   format('elem: B!FLAGS FOR Mdle = ',i5)
  do i=1,NR_PHYSA
    write(*,7002) PHYSA(i), ibc(1:nrf,i)
7002     format('     ATTRIBUTE = ',a6,' FLAGS = ',6i2)
  enddo
  endif
  iprint = 0

  !do i=1,MAXbrickH
  !write(*,*) 'i = ', i
  !write(*,*) 'xnod(1,i) = ', xnod(1,i)
  !write(*,*) 'xnod(2,i) = ', xnod(2,i)
  !write(*,*) 'xnod(3,i) = ', xnod(3,i)
  !enddo
  !call pause
!
! ! !..............FOR LONG WAVEGUIDE: RE-USING STIFFNESS MATRICES
! ...check if one needs to recompute element matrices
!   idec=0
!   call copy_element_matrices(MdE,MdQ,xnod,ibc, idec, &
!                              ZalocEE,ZalocEQ, &
!                              ZalocQE,ZalocQQ,ZblocE,ZblocQ)
   if (idec.eq.1) return
!..............FOR LONG WAVEGUIDE: RE-USING STIFFNESS MATRICES
! ... get current solution dofs
  call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
  call find_domain(Mdle, ndom)
! ... select case of GEOM_NO to set
! ... refractive index according to domain
  select case(GEOM_NO)
    case(1)
    bg_gain = 1.d0
    case(2)
!   bg_gain = REF_INDEX_CORE**2
    bg_gain = 1.d0
    case(3)
    select case(ndom)
      case(1,2,3,4)
!      bg_gain = REF_INDEX_CORE**2
		bg_gain = 1.d0
      case(5,6,7,8)
!      bg_gain = REF_INDEX_CLAD**2
		bg_gain = 1.d0
    end select
!   write(*,*)'ndom, bg_gain is: ', ndom, bg_gain
    case(4)
!    bg_gain = REF_INDEX_CORE**2
	 bg_gain = 1.d0
    case(5)
    select case(ndom)
      case(1,2,3,4,5)
      bg_gain = REF_INDEX_CORE**2
      case(6,7,8,9)
      bg_gain = REF_INDEX_CLAD**2
    end select
! ... end select case of GEOM_NO
  end select
!

!
! ...clear space for stiffness matrix and rhsv:
  ZblocE = ZERO; ZblocQ = ZERO
  ZalocEE = ZERO; ZalocEQ = ZERO; ZalocQE = ZERO; ZalocQQ = ZERO
!
! ...clear space for auxiliary matrices
  BLOADE = ZERO; STIFFEE = ZERO; STIFFEQ = ZERO; AP_Maxwell = ZERO
  STIFF_ALLE = ZERO

! ...shortcut: adjust frequency
! ...element size:
  h = min(abs(xnod(1,2)-xnod(1,1)),abs(xnod(2,4)-xnod(2,1)), &
      abs(xnod(3,5)-xnod(3,1)))
  omeg = min(OMEGA,6.d0/h)


!
  omeg = OMEGA
! ...auxiliary constant1
  zb = ZI*omeg*bg_gain*EPSILON - SIGMA
  za = ZI*omeg*bg_gain*EPSILON + SIGMA
  zc = ZI*omeg*MU

!
!-----------------------------------------------------------------------
!
! ...element integrals...
!
! ...use the enriched order to set the quadrature
  INTEGRATION = NORD_ADD
  call set_3Dint(etype,norder, nint,xiloc,waloc)
  INTEGRATION = 0
!
! ... loop through integration points
  !write(*,*)'elem_dpgMaxwell: nint is: ', nint
  do l=1,nint
    xi(1:3) = xiloc(1:3,l)
    wa = waloc(l)
!
! .....determine element H1 shape functions (for geometry)
    call shape3H(etype,xi,norder,norient_edge,norient_face, &
          nrdofH,shapH,gradH)
!
! .....determine element H(curl) shape functions
    call shape3E(etype,xi,norder,norient_edge,norient_face, &
          nrdofE,shapE,curlE)
! .....determine element L2 shape functions
    call shape3Q(etype,xi,norder, nrdofQ,shapQ)
!
! ..... find current domain number to update ion_gain in core only
    if(GEOM_NO.ne.1) then
      call find_domain(Mdle, ndom)
      select case(ndom)
      case(1,2,3,4)
! ... if solving nonlinear ion_gain then update ion_gain
        if(NONLINEAR_FLAG.eq.1) then

! .... compute current H(curl) solution to update ion_gain
        zsolQ = ZERO
        do kE=1,nrdofQ
          zsolQ(1:6) = zsolQ(1:6) + zdofQ(1:6,kE)*shapQ(kE)
        enddo
! ... evaluate gain function
        EfieldNorm = (zsolQ(1)**2+zsolQ(2)**2 &
                          +zsolQ(3)**2)**0.5
!       write(*,*) 'from elem.F: EfieldNorm = '
        call get_gainFunction(EfieldNorm, gainFunction)
! ...  evaluate ion_gain = REF_INDEX*gainFunction = sqrt(bg_gain)*gainFunction
        ion_gain = sqrt(bg_gain)*gainFunction
! ...  evaluate auxiliary constant1
        zb = ZI*OMEGA*bg_gain*EPSILON - ion_gain-SIGMA
        za = ZI*OMEGA*bg_gain*EPSILON + ion_gain+SIGMA
! ...... end if for NONLINEAR FLAG check
        endif
! ...... end select for checking core domain number
      end select
! ...... end if fir GEO_NO != 1
    endif
!
! .....determine discontinuous H(curl) shape functions
    call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
!
! .....geometry
    call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, &
                   x,dxdxi,dxidx,rjac,iflag)
!..................................................................
!..................................................................
!.........COMPUTE H1 Solution for Maxwell
!......... current heat solution u
!..................................................................
! .....compute the H1 solution at the point
    alpha_scale = 1.d0
    if(LASER_MODE.eq.1) then
      zsolH = ZERO
      do kH=1,nrdofH
        zsolH = zsolH + zdofH(1,kH)*shapH(kH)
      enddo
        alpha_scale = 2.d0*zsolH+1.d0
    endif

! .....integration weight
    weight = rjac*wa
!
! .....get the RHS
    call getf(Mdle,x, zfval,zJ)
!
! .....1st loop through enriched H(curl) test functions
    do k1=1,nrdofEE
!
! .......Piola transform
      E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                 + shapEE(2,k1)*dxidx(2,1:3) &
                 + shapEE(3,k1)*dxidx(3,1:3)
      curlE1(1:3) = dxdxi(1:3,1)*curlEE(1,k1) &
                     + dxdxi(1:3,2)*curlEE(2,k1) &
                     + dxdxi(1:3,3)*curlEE(3,k1)
      curlE1(1:3) = curlE1(1:3)/rjac
!
! .......compute the RHS
      k = 2*k1
      BLOADE(k) = BLOADE(k) &
           + (zJ(1)*E1(1)+zJ(2)*E1(2)+zJ(3)*E1(3))*weight
!
! ....... 2nd loop through enriched H(curl) test functions
      do k2=k1,nrdofEE
!
! .........Piola transform
      E2(1:3) = shapEE(1,k2)*dxidx(1,1:3) &
                   + shapEE(2,k2)*dxidx(2,1:3) &
                   + shapEE(3,k2)*dxidx(3,1:3)
      curlE2(1:3) = dxdxi(1:3,1)*curlEE(1,k2) &
                   + dxdxi(1:3,2)*curlEE(2,k2) &
                   + dxdxi(1:3,3)*curlEE(3,k2)
      curlE2(1:3) = curlE2(1:3)/rjac
!
! .........auxiliary quantities
      CC  = curlE1(1)*curlE2(1) + curlE1(2)*curlE2(2) &
              + curlE1(3)*curlE2(3)
      EE = E1(1)*E2(1) + E1(2)*E2(2) + E1(3)*E2(3)
      CE = curlE1(1)*E2(1) + curlE1(2)*E2(2) + curlE1(3)*E2(3)
      EC  = E1(1)*curlE2(1) + E1(2)*curlE2(2) + E1(3)*curlE2(3)
! .........accumulate for the test Gram matrix: adjoint graph
      k = nk(2*k1-1,2*k2-1)
      AP_Maxwell(k) = AP_Maxwell(k) &
                   + (CC+ (abs(zb)**2 + 1.d0)*EE)*weight
!                   + (CC+ (abs(ZI*OMEGA*EPSILON-SIGMA)**2 + 1.d0)*EE)*weight
      k = nk(2*k1-1,2*k2  )
      AP_Maxwell(k) = AP_Maxwell(k) &
                  + (zc*CE + (zb)*EC)*weight
                  !+ (ZI*OMEGA*MU*CE + (ZI*OMEGA*EPSILON-SIGMA)*EC)*weight
      if (k1.ne.k2) then
        k = nk(2*k1  ,2*k2-1)
        AP_Maxwell(k) = AP_Maxwell(k) &
                    + (conjg(zb)*CE -zc*EC)*weight
                    !+ (conjg(ZI*OMEGA*EPSILON-SIGMA)*CE -ZI*OMEGA*MU*EC)*weight
      endif
      k = nk(2*k1  ,2*k2  )
      AP_Maxwell(k) = AP_Maxwell(k) &
                + (CC+ ((abs(zc))**2 + 1.d0)*EE)*weight
                !+ (CC+ ((abs(ZI*OMEGA*MU))**2 + 1.d0)*EE)*weight
! ! .........accumulate for the test Gram matrix: mathematicians norm
!       k = nk(2*k1-1,2*k2-1)
!       AP_Maxwell(k) = AP_Maxwell(k) &
!                    + (CC+ 1.d0*EE)*weight
!       k = nk(2*k1  ,2*k2  )
!       AP_Maxwell(k) = AP_Maxwell(k) &
!                 + (CC+ 1.d0*EE)*weight
!....... enddo 2nd loop over enriched test H(curl) function
      enddo
!
! .......loop through L2 trial functions
      do k2=1,nrdofQ
        q = shapQ(k2)/rjac
! .........accumulate for the extended stiffness matrix
        do icomp=1,3
          k = (k2-1)*6 + icomp
          STIFFEQ(2*k1-1,k) = STIFFEQ(2*k1-1,k) &
                 + curlE1(icomp)*q*weight
          k = (k2-1)*6 + 3+ icomp
          STIFFEQ(2*k1-1,k) = STIFFEQ(2*k1-1,k) &
                 +zc*E1(icomp)*q*weight
          k = (k2-1)*6 + icomp
          STIFFEQ(2*k1  ,k) = STIFFEQ(2*k1  ,k) &
               -(za)*alpha_scale*E1(icomp)*q*weight
          k = (k2-1)*6 + 3+ icomp
          STIFFEQ(2*k1  ,k) = STIFFEQ(2*k1  ,k) &
               + curlE1(icomp)*q*weight
! .......end loop through icomp
        enddo
! .......end 1st loop through L2 trial functions
      enddo
! ....... end loop over through enriched H(curl) test functions
    enddo
! ....... end loop over integration points
  enddo
! .... printing Gram matrix
  iprint = 0
  if (iprint.eq.1) then
  write(*,*) 'elem: AP_Maxwell = '
  do i=1,10
    do j=1,i-1
      zaux(j) = AP_Maxwell(nk(j,i))
    enddo
    do j=i,10
      zaux(j) = AP_Maxwell(nk(i,j))
    enddo
    write(*,7011) zaux
  enddo
  call pause
  endif
!
!-----------------------------------------------------------------------
!
! ...boundary integrals
!
  if(GEOM_NO.eq.1) then
    !impedanceConstant = sqrt(1.d0-(PI**2/OMEGA**2))
    impedanceConstant = GAMMA_IMP
  else
    impedanceConstant = GAMMA_IMP
  endif
! ...loop through element faces
  do if=1,nrf
!
! .....sign factor to determine the OUTWARD normal unit vector
  nsign = nsign_param(etype,if)
!
! .....face type
  ftype = face_type(etype,if)
!
! .....face order of approximation
  call face_order(etype,if,norder, norderf)
!
! .....set 2D quadrature
  INTEGRATION = NORD_ADD
  call set_2Dint(ftype,norderf, nint,tloc,wtloc)
  INTEGRATION = 0
!
! .....loop through integration points
  !write(*,*)'elem_dpgMaxwell: nint boundary is: ', nint
  do l=1,nint
!
! .......face coordinates
    t(1:2) = tloc(1:2,l)
!
! .......face parametrization
    call face_param(etype,if,t, xi,dxidt)
!
! .......determine discontinuous Hcurl shape functions
    call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
!
! .......determine element H1 shape functions (for geometry)
    call shape3H(etype,xi,norder,norient_edge,norient_face, &
                     nrdofH,shapH,gradH)
!
! .......determine element H(curl) shape functions (for fluxes)
    call shape3E(etype,xi,norderc,norient_edge,norient_face, &
                     nrdofE,shapE,curlE)
!
! .......geometry
    call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                     x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
    weight = bjac*wtloc(l)
! ......check for impedance BC
! .......impedance boundary ......................................
    if (ibc(if,2).eq.9) then
      !write(*,*) 'putting impedance BCs'
!
! .........get the boundary source
      call get_bdSource(Mdle,x,rn, zImp)
!
! .........loop through Hcurl enriched test functions
      do k1=1,nrdofEE
!
! ...........value of the shape function at the point
        qq(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                     + shapEE(2,k1)*dxidx(2,1:3) &
                     + shapEE(3,k1)*dxidx(3,1:3)
!
! ...........accumulate for the load vector
      k = 2*k1
      BLOADE(k) = BLOADE(k) &
          - (zImp(1)*qq(1)+zImp(2)*qq(2)+zImp(3)*qq(3))*weight
!
! ...........loop through Hcurl trial functions
        do k2=1,nrdofE
!
! .............value of the shape function at the point
          p(1:3) = shapE(1,k2)*dxidx(1,1:3) &
                      + shapE(2,k2)*dxidx(2,1:3) &
                      + shapE(3,k2)*dxidx(3,1:3)
!
          call cross_product(rn,p, rntimesp)
          call cross_product(rn,rntimesp, rn2timesp)
!
! .............accumulate for the extended stiffness matrix
          STIFFEE(2*k1,2*k2-1) = STIFFEE(2*k1,2*k2-1) &
                + (qq(1)*rn2timesp(1)+ qq(2)*rn2timesp(2)+ &
                + qq(3)*rn2timesp(3))*(impedanceConstant)*weight

          STIFFEE(2*k1-1,2*k2-1) = STIFFEE(2*k1-1,2*k2-1) &
                +  (qq(1)*rntimesp(1)+ qq(2)*rntimesp(2)+ &
                +   qq(3)*rntimesp(3))*weight
! ...........end loop through Hcurl trial functions
        enddo
! ...........end loop through Hcurl enriched test functions
      enddo
! .......regular boundary.............................................
    else
! ....... check for LASER MODE, if
! ....... LASER MODE is active, get boundary data for source
!
    if(LASER_MODE.eq.1) then
! .........get the boundary source
      !call get_bdSource(Mdle,x,rn,0, zImp)
    else
      zImp = ZERO
    end if
! .......loop through the enriched H(curl) test functions
    do k1=1,nrdofEE
      E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                   + shapEE(2,k1)*dxidx(2,1:3) &
                   + shapEE(3,k1)*dxidx(3,1:3)
! ...........accumulate for the load vector
      k = 2*k1
      BLOADE(k) = BLOADE(k) &
          + (zImp(1)*E1(1)+zImp(2)*E1(2)+zImp(3)*E1(3))*weight
! .........loop through H(curl) trial functions
      do k2=1,nrdofE
        E2(1:3) = shapE(1,k2)*dxidx(1,1:3) &
                     + shapE(2,k2)*dxidx(2,1:3) &
                     + shapE(3,k2)*dxidx(3,1:3)
      call cross_product(rn,E2, rntimesE)
!
! ...........accumulate for the extended stiffness matrix
      STIFFEE(2*k1-1,2*k2-1) = STIFFEE(2*k1-1,2*k2-1) &
            + (E1(1)*rntimesE(1)+E1(2)*rntimesE(2)+E1(3)*rntimesE(3)) &
               *weight
      STIFFEE(2*k1  ,2*k2  ) = STIFFEE(2*k1  ,2*k2  ) &
            + (E1(1)*rntimesE(1)+E1(2)*rntimesE(2)+E1(3)*rntimesE(3)) &
               *weight
! .........end loop through H(curl) trial functions
      enddo
! .......end loop through the enriched H(curl) test functions
    enddo
!.......end if for impedance BC
    endif
! ......end loop through integration points
  enddo
!......end loop through faces
  enddo
!
! printing element matrices
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
      write(*,7210) BLOADE(j) - STIFFEQ(j,1) &
          - (STIFFEE(j,1)+STIFFEE(j,5)+STIFFEE(j,9)+STIFFEE(j,13))
    endif
    call pause
  enddo
  endif

!! ....................................................................
!! ...  construction of DPG system
!! ....................................................................
!
!!  ...total trial dof for the element
  nrTEST = nrdofEE*2
  call celndof(NODES(Mdle)%type,norder, &
              nrdofH,nrdofE,nrdofV,nrdofQ)
  i1 = MAXtestE ; j1 = 2*nrdofE ; j2 = 6*nrdofQ
!
  STIFF_ALLE(1:i1,1:j1) = STIFFEE(1:i1,1:j1)
  STIFF_ALLE(1:i1,j1+1:j1+j2) = STIFFEQ(1:i1,1:j2)
  STIFF_ALLE(1:i1,j1+j2+1) = BLOADE(1:i1)
!
!
  uplo  = 'U' ; trans = 'C' ;   diag  = 'N'
  N     = nrTEST
  nRHS  = 2*nrdofE + 6*nrdofQ + 1
!
!----Preconditioning AP_Maxwell= D^-1/2 * AP_Maxwell * D^-1/2--------
!----------------------------------------------------------------------
  allocate(DIAG_E(nrTest))
! ...preconditioning: getting diagonal entries
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
! ...preconditioning: STIFF_ALLE = D^-1/2 * STIFF_ALLE
  do k2=1,NRHS
    do k1=1,nrTEST
      STIFF_ALLE(k1,k2) = STIFF_ALLE(k1,k2)/sqrt(DIAG_E(k1))
    enddo
  enddo
  deallocate(DIAG_E)

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
! !.....check condition number
!      nrdof_eig = nrdofEE*2
!      kk = nrdof_eig*(nrdof_eig+1)/2
!      allocate(AP_eig(kk))
!      AP_eig(1:kk) = AP_Maxwell(1:kk)
!      allocate(W(nrdof_eig))
!      allocate(Z(1,nrdof_eig))
!      allocate(WORK(nrdof_eig))
!      allocate(RWORK(nrdof_eig))
!      allocate(IWORK(1))
!      call ZHPEVD('N','U',nrdof_eig, &
!                  AP_eig, W,Z,1,WORK,nrdof_eig, &
!                  RWORK,nrdof_eig,IWORK,1,info)
!      if (info .ne. 0) then
!         write(*,*) 'eig_solve_sc: info = ', info
!         stop 1
!      endif

!      write(*,6999) W(nrdof_eig),W(1)
! 6999 format('elem_dpgMaxwell: AP_Maxwell: max_eig, min_eig = ', 2e13.4)


!      write(*,7000)  W(nrdof_eig)/W(1)
! 7000 format('elem_dpgMaxwell: AP_Maxwell condition number = ', 1e13.4)
!      deallocate(IWORK,W,WORK)
!      deallocate(RWORK,Z)
!      deallocate(AP_eig)
!      call pause
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------

! ...factorize the test stiffness matrix
  uplo = 'U'
  call ZPPTRF(uplo, nrTEST, AP_Maxwell, info)
  if (info.ne.0) then
    write(*,*) 'elem_dpgMaxwell: info = ',info
    stop
  endif
  call ZTPTRS(uplo,trans,diag,N,NRHS,AP_Maxwell,STIFF_ALLE, &
              MAXtestE,info)
!
  call ZHERK(uplo,trans,NRHS,nrTEST,ZONE,STIFF_ALLE, &
              MAXtestE,ZERO, &
              STIFF_ALLE(1:NRHS,1:NRHS),NRHS)
!    ZHERK for complex case
  do i=1,NRHS-1
    STIFF_ALLE(i+1:NRHS,i) = conjg(STIFF_ALLE(i,i+1:NRHS))
  enddo
!

!
  ZblocE(1:j1) = STIFF_ALLE(1:j1,j1+j2+1)
  ZblocQ(1:j2) = STIFF_ALLE(j1+1:j1+j2,j1+j2+1)
!
  ZalocEE(1:j1,1:j1) = STIFF_ALLE(1:j1,1:j1)
  ZalocEQ(1:j1,1:j2) = STIFF_ALLE(1:j1,j1+1:j1+j2)
!
  ZalocQE(1:j2,1:j1) = STIFF_ALLE(j1+1:j1+j2,1:j1)
  ZalocQQ(1:j2,1:j2) = STIFF_ALLE(j1+1:j1+j2,j1+1:j1+j2)
!

!   if (idec.ne.2) then
!    write(*,*) 'elem_dpgMaxwell: idec = ',idec
!    stop 1
!   endif
!   call copy_element_matrices(MdE,MdQ,xnod,ibc, idec, &
!                             ZalocEE,ZalocEQ, &
!                             ZalocQE,ZalocQQ,ZblocE,ZblocQ)
  !write(*,*)'ZalocEE is: ', ZalocEE(1:MdE,1:MdE)
 ! write(*,*) '   '
 ! write(*,*)'ZalocEQ is: ', ZalocEE(1:MdE,1:MdQ)
 ! write(*,*) '   '
 ! write(*,*)'ZalocQE is: ', ZalocEE(1:MdQ,1:MdE)
 ! write(*,*) '   '
 ! write(*,*)'ZalocQQ is: ', ZalocEE(1:MdQ,1:MdQ)
 ! write(*,*) '   '
 ! call pause

 ! iprint = 2
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

end subroutine elem_dpgMaxwell

!-------------------------------------------------------------------
!Routine to evaluare gain function
!Last modified : April 17
!input: EfieldNorm (norm of E field)
!output: gainFunction (value of gain function)
!-------------------------------------------------------------------
subroutine get_gainFunction(EfieldNorm, gainFunction)
!  use DPGLaser
  use problem
  implicit none
!
  real*8,                       intent(in)  :: EfieldNorm
  real*8,                       intent(out)  :: gainFunction
  real*8 :: alpha,f_em,f_abs,wEratio,Nratio,Nex,Ngd
  real*8 :: sigma_s_abs_scaled,sigma_s_em_scaled
  real*8 :: sigma_p_abs_scaled, sigma_p_em_scaled

  sigma_s_abs_scaled = SIGMA_S_ABS*N0*CORE_RAD
  sigma_s_em_scaled = SIGMA_S_EM*N0*CORE_RAD
  sigma_p_abs_scaled = SIGMA_P_ABS*N0*CORE_RAD
  sigma_p_em_scaled = SIGMA_P_EM*N0*CORE_RAD
!
  alpha = (TAU*REF_INDEX_CORE*E_0**2)/(MU0*LIGHT_SPEED*CORE_RAD* &
                                       N0*H_BAR*OMEGA0)
  f_em = sigma_s_em_scaled+((1.d0/OMEGA_RATIO)*sigma_p_em_scaled)
  f_abs = sigma_s_abs_scaled+((1.d0/OMEGA_RATIO)*sigma_p_abs_scaled)
  wEratio = OMEGA/(alpha*EfieldNorm**2)
  Nratio= (wEratio+f_em)/f_abs
  Nex = 1.0/(1.0+Nratio)
  Ngd = 1.0-Nex
  gainFunction = sigma_s_em_scaled *Nex - sigma_s_abs_scaled*Ngd
!     write(*,*) 'from get_gainFunction: EfieldNorm = '
!    .             EfieldNorm
end subroutine get_gainFunction





subroutine copy_element_matrices(MdE,MdQ,Xnod,ibc, Idec, &
                                 ZalocEE,ZalocEQ, &
                                 ZalocQE,ZalocQQ,ZblocE,ZblocQ)
!
  use control
  use parameters
  use matrices
  implicit none
!
#if C_MODE
#define VTYPE  complex*16
#else
#define VTYPE double precision
#endif
!.......declare input/output variables
  integer,                     intent(in)    :: MdE
  integer,                     intent(in)    :: MdQ
  integer, dimension(6,4), intent(in) :: ibc
  VTYPE, dimension(MdE,MdE) :: ZalocEE
  VTYPE, dimension(MdE,MdQ) :: ZalocEQ
  VTYPE, dimension(MdQ,MdE) :: ZalocQE
  VTYPE, dimension(MdQ,MdQ) :: ZalocQQ
  VTYPE, dimension(MdE) :: ZblocE
  VTYPE, dimension(MdQ) :: ZblocQ
!
! ...element geometry dof
  real*8, dimension(3,MAXbrickH),intent(in) :: Xnod
  integer, intent(inout)  :: Idec
!
  integer :: iprint,ilayer,j,i
!
  iprint=1
  ZblocE = ZERO; ZblocQ = ZERO
  !write(*,*) 'copying matrices'
!
    select case(Idec)
!
! ..check if one should recompute the element matrices
    case(0)
! ... check if IBC is imposed, if so, always compute matrices
      do j=1,6
        if(ibc(j,2).eq.9) then
          write(*,*) 'ibc is imposed'
          Idec = 2
          return
        endif
      enddo
!
! ....find the corresponding element in the first layer
      do ilayer=1,NRFL
        if ((abs(Xnod(1,1)-XYVERT(1,1,ilayer)).lt.GEOM_TOL).and. &
            (abs(Xnod(2,1)-XYVERT(2,1,ilayer)).lt.GEOM_TOL).and. &
            (abs(Xnod(1,2)-XYVERT(1,2,ilayer)).lt.GEOM_TOL).and. &
            (abs(Xnod(2,2)-XYVERT(2,2,ilayer)).lt.GEOM_TOL).and. &
            (abs(Xnod(1,3)-XYVERT(1,3,ilayer)).lt.GEOM_TOL).and. &
            (abs(Xnod(2,3)-XYVERT(2,3,ilayer)).lt.GEOM_TOL)) go to 10
      enddo
!
! ....no element with this position has been computed yet
		write(*,*) 'computing matrices'
!
      Idec=2; return
!
! ....element matrices have been computed
10    Idec=1
      write(*,*) 'not computing the stiffness matrices!'
!      write(*,*) 'ilayer: ', ilayer
!      write(*,*) 'NRFL: ', NRFL
!      write(*,*) 'xnod is'
!	  write(*,*) 'xnod(1,1) is', Xnod(1,1)
!	  write(*,*) 'xnod(2,1) is', Xnod(2,1)
!	  call pause

! ....copy the element matrices from the temporary data structure
!!$OMP CRITICAL
      ZalocEE(1:MdE,1:MdE) = ZFL_EE(1:MdE,1:MdE,ilayer)
      ZalocEQ(1:MdE,1:MdQ) = ZFL_EQ(1:MdE,1:MdQ,ilayer)
      ZalocQQ(1:MdQ,1:MdQ) = ZFL_QQ(1:MdQ,1:MdQ,ilayer)
!      ZalocQE(1:MdQ,1:MdE) = ZFL_QE(1:MdQ,1:MdE,ilayer)
!      ZblocE(1:MdE) = ZLOADFL_E(1:MdE,ilayer)
!      ZblocQ(1:MdQ) = ZLOADFL_Q(1:MdQ,ilayer)
!      write(*,*) 'ilayer,NRFL is ', ilayer, NRFL
!      write(*,*) ''
!!$OMP END CRITICAL
!      write(*,8010)
!8010   format('elem_dpgMaxwell: ZblocE,ZblocQ = ')
!    write(*,8011) ZLOADFL_E(1:MdE,ilayer)
!    write(*,8011) ZLOADFL_Q(1:MdQ,ilayer)
!8011   format(10e12.5)
!    call pause
!    write(*,8012)
!8012   format('elem_dpgMaxwell: ZalocEE = ')
!    do i=1,MdE
!    write(*,8013) i,ZFL_EE(i,1:MdE,ilayer)
!8013     format('i = ',i3,10(/,5(2e12.5,2x)))
!    enddo
!    call pause
!    write(*,8014)
!8014   format('elem_dpgMaxwell: ZalocEQ = ')
!   do i=1,MdE
!      write(*,8013) i,ZFL_EQ(i,1:MdQ,ilayer)
!    enddo
!    call pause
!    write(*,8015)
!8015   format('elem_dpgMaxwell: ZalocQQ = ')
!    do i=1,MdQ
!      write(*,8013) i,ZFL_QQ(i,1:MdQ,ilayer)
!    enddo
!    call pause
!
! ....finish generating the matrices
       do i=1,MdQ
         ZalocQE(i,1:MdE) = conjg(ZalocEQ(1:MdE,i))
       enddo
      return
!
! ..store the element matrices
    case(2)
!
      !write(*,*) 'ilayer: ', ilayer
      !write(*,*) 'NRFL: ', NRFL
! ....store first vertex node xy coordinates
!!$OMP CRITICAL
      NRFL=NRFL+1
      if (NRFL.gt.MAXNRFL) then
        write(*,*) 'copy_element_matrices: INCREASE MAXNRFL'
        stop 1
      endif
!      write(*,*) 'computing the stiffness matrices!'
!		write(*,*) 'xnod is'
!	    write(*,*) 'xnod(1,1) is', Xnod(1,1)
!	    write(*,*) 'xnod(2,1) is', Xnod(2,1)
!	    write(*,*) 'NRFL: ', NRFL
!	    call pause
      XYVERT(1:2,1:3,NRFL) = Xnod(1:2,1:3)
!
! ....copy the element matrices from the temporary data structure
      ZFL_EE(1:MdE,1:MdE,NRFL) = ZalocEE(1:MdE,1:MdE)
      ZFL_EQ(1:MdE,1:MdQ,NRFL) = ZalocEQ(1:MdE,1:MdQ)
      ZFL_QQ(1:MdQ,1:MdQ,NRFL) = ZalocQQ(1:MdQ,1:MdQ)
!     ZFL_QE(1:MdQ,1:MdE,NRFL) = ZalocQE(1:MdQ,1:MdE)
!     ZLOADFL_E(1:MdE,NRFL) =  ZblocE(1:MdE)
!     ZLOADFL_Q(1:MdQ,NRFL) = ZblocQ(1:MdQ)
!!$OMP END CRITICAL
!
    end select
!
!
  end subroutine copy_element_matrices





! !----------------------------------------------------------------------
! ! c   ..penalty constraint

!   allocate(zConstraintMatrix(2,6*nrdofQ+2*NrdofE))
!   allocate(zConstraintLoad(2))
! !
!   alpha = 1000.d0
!   zConstraintLoad = ZERO; zConstraintMatrix = ZERO
! !   ...boundary integrals...

! !   ...loop through element edges
! ! ...loop through element faces
!   do if=1,nrf
! !
! ! .....sign factor to determine the OUTWARD normal unit vector
!   nsign = nsign_param(etype,if)
! !
! ! .....face type
!   ftype = face_type(etype,if)
! !
! ! .....face order of approximation
!   call face_order(etype,if,norder, norderf)
! !
! ! .....set 2D quadrature
!   INTEGRATION = NORD_ADD
!   call set_2Dint(ftype,norderf, nint,tloc,wtloc)
!   INTEGRATION = 0
! !
! ! .....loop through integration points
!   write(*,*)'elem_dpgMaxwell: nint boundary is: ', nint
!   do l=1,nint
! !
! ! .......face coordinates
!     t(1:2) = tloc(1:2,l)
! !
! ! .......face parametrization
!     call face_param(etype,if,t, xi,dxidt)
! !
! ! .......determine element H1 shape functions (for geometry)
!     call shape3H(etype,xi,norder,norient_edge,norient_face, &
!                      nrdofH,shapH,gradH)
! !
! ! .......determine element H(curl) shape functions (for fluxes)
!     call shape3E(etype,xi,norderc,norient_edge,norient_face, &
!                      nrdofE,shapE,curlE)
! !
! ! .......geometry
!     call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
!                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
!     weight = bjac*wtloc(l)
! ! .........loop through Hcurl trial functions
!       do k1=1,nrdofE
! !
! ! ...........value of the shape function at the point
!         qq(1:3) = shapE(1,k1)*dxidx(1,1:3) &
!                      + shapE(2,k1)*dxidx(2,1:3) &
!                      + shapE(3,k1)*dxidx(3,1:3)
!         dot_product(rn,qq, rndotE)
! !
! ! ...........accumulate for the load vector
!       zConstraintMatrix = alpha*zConstraintMatrix
!       zConstraintLoad = alpha*zConstraintLoad
!       k = 2*k1-1
!       zConstraintMatrix(k) = zConstraintMatrix(k) &
!                     + rndotE*weight
!       k = 2*k1
!       zConstraintMatrix(k) = zConstraintMatrix(k) &
!                     + rndotE*weight
! !...    end loop over trial space
!       enddo
! !...    end loop over integration points
!     enddo
! !...    end loop over faces
!   enddo

!   call ZHERK(uplo,trans,NRHS-1,2,ZONE,zConstraintMatrix, &
!               NRHS-1,ZONE, &
!               STIFF_ALLE(1:NRHS-1,1:NRHS-1),NRHS-1)
!       STIFF_ALLE(NRHS,1:NRHS-1) = &
!             conjg(zConstraintMatrix(1,1:NRHS-1))*zConstraintLoad(1) &
!             +conjg(zConstraintMatrix(2,1:NRHS-1))*zConstraintLoad(2)
!   ...element integrals...
!       call set_2Dint(type,norderP, nint,xiloc,waloc)
!       do l=1,nint
!         xi(1:2) = xiloc(1:2,l)
!         wa = waloc(l)
!
!   .....determine element H1 shape functions (for geometry)
!         call shape2DH(type,xi,norder,norient, nrdofHf,shapH,gradH)
!
!   .....determine trial L2 shape functions
!         call shape2DQ(type,xi,norder, nrdofQ,shapQ)
!
!   .....geometry
!         call geom(xnod,shapH,gradH,nrdofHf, x,dxdxi,dxidx,rjac)
!         weight = rjac*wa
!
!   .....get source
!         call get_source(Mdle,xi,x, zf)
!
!           zd = zd + zf(3)*weight
!

!   .......loop through L2 trial shape functions
!           do k2=1,nrdofQ
!
!   .........Piola transformation
!             p = shapQ(k2)/rjac
!
!             m = (k2-1)*3+3
!             zCp(m) = zCp(m) + ZI*OMEGA*p*weight
!           enddo
!
!   ...end of loop through integration points
!       enddo



!       zCp = conjg(zCp)
!       zCu = conjg(zCu)


!    ...modify zalocVV for the constraints
!       call zgerc(nrdofV,nrdofV,alpha,zCu,1,zCu,1,zalocVV(1:j2,1:j2)
!      .            ,nrdofV)



!       call zgerc(3*nrdofQ,3*nrdofQ,alpha,zCp,1,zCp,1,
!      .            zalocQQ(1:j3,1:j3),3*nrdofQ)





!       call zgerc(nrdofV,3*nrdofQ,alpha,zCu,1,zCp,1,
!      .            zalocVQ(1:j2,1:j3),nrdofV)



!       call zgerc(3*nrdofQ,nrdofV,alpha,zCp,1,zCu,1,
!      .            zalocQV(1:j3,1:j2),3*nrdofQ)






!       ZblocV(1:j2) = ZblocV(1:j2) + alpha*zd*zCu(1:j2)
!       ZblocQ(1:j3) = ZblocQ(1:j3) + alpha*zd*zCp(1:j3)


! C  ...................................................................
! C  ...................................................................
! c
!       deallocate(zCu,zCp)
