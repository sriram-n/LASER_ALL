!--------------------------------------------------------------------

!     routine name      - elem
!
!--------------------------------------------------------------------
!
!     latest revision:  - August 17
!
!     purpose:          - routine returns unconstrained (ordinary)
!                         stiffness matrix and load vector
!                         for the primal DPG method for Maxwell
!                         eqn
!
!
!     arguments:
!
!     in:
!             Mdle      - an element middle node number, identified
!                         with the element
!             Maxdof    - column length of Zaloc
!     out:
!             Nrdof     - number of dof for a single component
!             Itest_e,Itest_a,Itrial_e,Itrial_a - flags indicating
!                         presence of corresponding load vectors
!                         and stiffness matrices
!             Zbloc     - load vector
!             Zaloc     - stiffness matrix
!
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
      case(3)
!
!  .....this is user-defined routine, you are in charge of it
        Itest(1:2)=1; Itrial(1:2)=1
        call elem_primalMaxwell(Mdle,BLOC(1)%nrow,BLOC(2)%nrow, &
                       BLOC(1)%array,ALOC(1,1)%array,ALOC(1,2)%array, &
                       BLOC(2)%array,ALOC(2,1)%array,ALOC(2,2)%array)
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
!     routine name      - elem_primalMaxwell
!
!--------------------------------------------------------------------
!
!     latest revision:  - September 17
!
!     purpose:          - routine returns unconstrained (ordinary)
!                         stiffness matrix and load vector
!                         for the Primal formulation for Maxwell
!                         equation
!
!     arguments:
!
!     in:
!             Mdle      - an element middle node number, identified
!                         with the element
!             MdE       - column length of ZalocEE,ZalocEF
!             MdF       - column length of ZalocFF,ZalocFE
!     out:
!             ZblocE,ZblocF - load vectors
!             ZalocEE,ZalocEF,ZalocFE,ZalocFF - stiffness matrices
!
!---------------------------------------------------------------------
!
      subroutine elem_primalMaxwell(Mdle,MdE,MdF, &
                              ZblocE,ZalocEE,ZalocEF, &
                              ZblocF,ZalocFE,ZalocFF)
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
  integer,                     intent(in)  :: MdF
  VTYPE, dimension(MdE),       intent(out) :: ZblocE
  VTYPE, dimension(MdE,MdE),   intent(out) :: ZalocEE
  VTYPE, dimension(MdE,MdF),   intent(out) :: ZalocEF
  VTYPE, dimension(MdF),       intent(out) :: ZblocF
  VTYPE, dimension(MdF,MdE),   intent(out) :: ZalocFE
  VTYPE, dimension(MdF,MdF),   intent(out) :: ZalocFF

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
!  VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
!  VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
!  VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
!  VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
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

! .... Enriched Hcurl shape functions
  real*8 , dimension(3,MAXbrickEE)    :: shapEE
  real*8 , dimension(3,MAXbrickEE)    :: curlEE
! ... nrdof for various spaces
  integer  :: nrdofH,nrdofE,nrdofV,nrdofQ,nrdofEE
! ....... space for DPG Computations (Gram Matrix, Stiffness etc.)
  integer, parameter  :: MAXtestE = MAXbrickEE
!  ...stiffness matrix for the local Riesz H1 matrix in Lapack format
  VTYPE, dimension(MAXtestE*(MAXtestE+1)/2) :: AP_Maxwell
!  ...load vector for the enriched space
  VTYPE, dimension(MAXtestE) :: BLOADE
!  ...stiffnes matrices for the enriched test space
!  VTYPE, dimension(MAXtestE,MAXbrickQ*6) :: STIFFEQ
  VTYPE, dimension(MAXtestE,MAXbrickE) :: STIFFEE
  VTYPE, dimension(MAXtestE,MAXbrickE) :: STIFFEF

!  ....STIFF_ALL for alternative computation of stiffness
  VTYPE, dimension(MAXtestE,2*MAXbrickE+1) :: STIFF_ALLE
#if C_MODE
  complex*16, allocatable :: AP_eig(:)
  complex*16, allocatable :: DIAG_E(:)
#else
  real*8, allocatable     :: AP_eig(:)
  real*8, allocatable     :: DIAG_E(:)
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
! ...for debug printing
  VTYPE, dimension(10)  :: zaux

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
  integer :: N,nRHS,nordP,nsign,if,ndom,info,info1,info2,info3,icomp,nrdof_eig,idec
  VTYPE   :: zfval,za,zb,zc,zk2
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
      !write(*,*) 'elem: Mdle = ',Mdle
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
 7001   format('elem_primalMaxwell: BC FLAGS FOR Mdle = ',i5)
        do i=1,NR_PHYSA
          write(*,7002) PHYSA(i), ibc(1:nrf,i)
 7002     format('          ATTRIBUTE = ',a6,' FLAGS = ',6i2)
        enddo
      endif
      iprint = 0
!
!  ...clear space for stiffness matrix and rhsv:
      ZblocE = ZERO; ZblocF = ZERO
      ZalocEE = ZERO; ZalocEF = ZERO; ZalocFE = ZERO; ZalocFF = ZERO
!
!  ...clear space for auxiliary matrices
      BLOADE = ZERO; STIFFEE = ZERO; STIFFEF = ZERO; AP_Maxwell = ZERO
      STIFF_ALLE = ZERO
!
!  ...complex wave number
      zk2 = OMEGA**2*EPSILON*MU - ZI*OMEGA*SIGMA*MU
!
!-----------------------------------------------------------------------
!
!  ...element integrals...
!
!  ...use the enriched order to set the quadrature
      INTEGRATION = NORD_ADD
      call set_3Dint(etype,norder, nint,xiloc,waloc)
      INTEGRATION = 0
!
! ... loop through integration points
      do l=1,nint
        xi(1:3) = xiloc(1:3,l)
        wa = waloc(l)
!
!  .....determine element H1 shape functions (for geometry)
        call shape3H(etype,xi,norder,norient_edge,norient_face, &
                    nrdofH,shapH,gradH)
!
!  .....determine element H(curl) shape functions
        call shape3E(etype,xi,norder,norient_edge,norient_face, &
                    nrdofE,shapE,curlE)
!
!  .....determine discontinuous H(curl) shape functions
        call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
!
!  .....geometry
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, &
                   x,dxdxi,dxidx,rjac,iflag)
!
!  .....integration weight
        weight = rjac*wa
!
!  .....get the RHS
        call getf(Mdle,x, zJ)
!
!  .....loop through enriched H(curl) test functions
        do k1=1,nrdofEE
          E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                 + shapEE(2,k1)*dxidx(2,1:3) &
                 + shapEE(3,k1)*dxidx(3,1:3)
          curlE1(1:3) = dxdxi(1:3,1)*curlEE(1,k1) &
                     + dxdxi(1:3,2)*curlEE(2,k1) &
                     + dxdxi(1:3,3)*curlEE(3,k1)
          curlE1(1:3) = curlE1(1:3)/rjac
!
!  .......compute the RHS
          BLOADE(k1) = BLOADE(k1) &
           + (zJ(1)*E1(1)+zJ(2)*E1(2)+zJ(3)*E1(3))*weight
!
!  .......loop through enriched H(curl) trial functions
          do k2=k1,nrdofEE
            E2(1:3) = shapEE(1,k2)*dxidx(1,1:3) &
                   + shapEE(2,k2)*dxidx(2,1:3) &
                   + shapEE(3,k2)*dxidx(3,1:3)
            curlE2(1:3) = dxdxi(1:3,1)*curlEE(1,k2) &
                       + dxdxi(1:3,2)*curlEE(2,k2) &
                       + dxdxi(1:3,3)*curlEE(3,k2)
            curlE2(1:3) = curlE2(1:3)/rjac
!
!  .........accumulate for the test stiffness matrix
            k = nk(k1,k2)
              AP_Maxwell(k) = AP_Maxwell(k) &
                   + ((E1(1)*E2(1)+E1(2)*E2(2)+E1(3)*E2(3)) &
                     + curlE1(1)*curlE2(1)+curlE1(2)*curlE2(2) &
                     + curlE1(3)*curlE2(3))*weight
          enddo
!
!  .......loop through Hcurl trial functions
          do k2=1,nrdofE
            E2(1:3) = shapE(1,k2)*dxidx(1,1:3) &
                   + shapE(2,k2)*dxidx(2,1:3) &
                   + shapE(3,k2)*dxidx(3,1:3)
            curlE2(1:3) = dxdxi(1:3,1)*curlE(1,k2) &
                       + dxdxi(1:3,2)*curlE(2,k2) &
                       + dxdxi(1:3,3)*curlE(3,k2)
            curlE2(1:3) = curlE2(1:3)/rjac
!
!  .........accumulate for the extended stiffness matrix
            STIFFEE(k1,k2) = STIFFEE(k1,k2) &
           + ((curlE1(1)*curlE2(1)+curlE1(2)*curlE2(2) &
                        +curlE1(3)*curlE2(3)) &
             -zk2*(E1(1)*E2(1)+E1(2)*E2(2)+E1(3)*E2(3)))*weight
           enddo
        enddo
      enddo

      iprint = 0
      if (iprint.eq.1) then
        write(*,*) 'elem_primalMaxwell: AP_Maxwell = '
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
      iprint = 0
!

!-----------------------------------------------------------------------
!
!  ...boundary integrals
!
!  ...loop through element faces
      do if=1,nrf
!
!  ..... bc of faces:
       ! write(*,*) 'ibc(if,1) is : ', ibc(if,1)
       !write(*,*) 'ibc(if,2) is : ', ibc(if,2)
!  .....sign factor to determine the OUTWARD normal unit vector
        nsign = nsign_param(etype,if)
!
!  .....face type
        ftype = face_type(etype,if)
!
!  .....face order of approximation
        call face_order(etype,if,norder, norderf)
!
!  .....set 2D quadrature
        INTEGRATION = NORD_ADD
        call set_2Dint(ftype,norderf, nint,tloc,wtloc)
        INTEGRATION = 0
!
!  .....loop through integration points
        do l=1,nint
!
!  .......face coordinates
          t(1:2) = tloc(1:2,l)
!
!  .......face parametrization
          call face_param(etype,if,t, xi,dxidt)
!
!  .......determine discontinuous Hcurl shape functions
          call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
!
!  .......determine element H1 shape functions (for geometry)
          call shape3H(etype,xi,norder,norient_edge,norient_face, &
                       nrdofH,shapH,gradH)
!
!  .......determine element H(curl) shape functions (for fluxes)
          call shape3E(etype,xi,norder,norient_edge,norient_face, &
                       nrdofE,shapE,curlE)
!
!  .......geometry
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                       x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
          weight = bjac*wtloc(l)
!
! ........impedance boundary ......................................
          if (ibc(if,2).eq.9) then
            !write(*,*) 'putting impedance BCs'
!
! ..........get the boundary source
            call get_bdSource(Mdle,x,rn, zImp)
            !write(*,*) 'elem_primalMaxwell rn is ', rn

!
! ..........loop through Hcurl enriched test functions
            do k1=1,nrdofEE
!
! ............value of the shape function at the point
              qq(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                      + shapEE(2,k1)*dxidx(2,1:3) &
                      + shapEE(3,k1)*dxidx(3,1:3)
!
! ............accumulate for the load vector
              BLOADE(k1) = BLOADE(k1) &
                         + ZI*OMEGA*MU*(zImp(1)*qq(1)+zImp(2)*qq(2)+zImp(3)*qq(3))*weight
!
! ............loop through Hcurl trial functions
              do k2=1,nrdofE
!
! ..............value of the shape function at the point
                p(1:3) = shapE(1,k2)*dxidx(1,1:3) &
                       + shapE(2,k2)*dxidx(2,1:3) &
                       + shapE(3,k2)*dxidx(3,1:3)
!
                call cross_product(rn,p, rntimesp)
                call cross_product(rn,rntimesp, rn2timesp)
                !write(*,*) 'rn2timesp is : ', rn2timesp
                !call pause
!
! ..............accumulate for the extended stiffness matrix
                STIFFEE(k1,k2) = STIFFEE(k1,k2) &
                      - ZI*OMEGA*MU*(GAMMA)*  &
                      (qq(1)*rn2timesp(1)+qq(2)*rn2timesp(2)+qq(3)*rn2timesp(3))*weight
!
! ............end loop through Hcurl trial functions
              enddo
!
! ..........end loop through Hcurl enriched test functions
            enddo
!
! ........regular boundary.............................................
          else
!  .........loop through the enriched H(curl) test functions
            do k1=1,nrdofEE
              E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                      + shapEE(2,k1)*dxidx(2,1:3) &
                      + shapEE(3,k1)*dxidx(3,1:3)
!
!  .........  loop through H(curl) trial functions
              do k2=1,nrdofE
                E2(1:3) = shapE(1,k2)*dxidx(1,1:3) &
                        + shapE(2,k2)*dxidx(2,1:3) &
                        + shapE(3,k2)*dxidx(3,1:3)
                call cross_product(rn,E2, rntimesE)
!
!  .............accumulate for the extended stiffness matrix
                STIFFEF(k1,k2) = STIFFEF(k1,k2) &
                + (E1(1)*rntimesE(1)+E1(2)*rntimesE(2)+E1(3)*rntimesE(3))*weight
!
! ............end loop through H(curl) trial functions
              enddo
!
! ..........end loop through the enriched H(curl) test functions
            enddo
!
! ........end if for impedance BC
          endif
!
! ......end loop through integration points
        enddo
!
! ....end loop through faces
      enddo
!
!-----------------------------------------------------------------------
!! ....................................................................
!! ...  construction of DPG system
!! ....................................................................
!!  ...total trial dof for the element
  nrTEST = nrdofEE
  call celndof(NODES(Mdle)%type,norder, &
              nrdofH,nrdofE,nrdofV,nrdofQ)
  i1 = MAXtestE ; j1 = nrdofE ; j2 = nrdofE
!
  STIFF_ALLE(1:i1,1:j1) = STIFFEE(1:i1,1:j1)
  STIFF_ALLE(1:i1,j1+1:j1+j2) = STIFFEF(1:i1,1:j2)
  STIFF_ALLE(1:i1,j1+j2+1) = BLOADE(1:i1)
!
!
  uplo  = 'U' ; trans = 'C' ;   diag  = 'N'
  N     = nrTEST
  nRHS  = 2*nrdofE + 1
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
  call ZPPTRF(uplo, nrdofEE, AP_Maxwell, info)
  if (info.ne.0) then
    write(*,*) 'elem_primalMaxwell: info = ',info
    write(*,*) 'elem_primalMaxwell: AP_Maxwell third pass '
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
  ZblocF(1:j2) = STIFF_ALLE(j1+1:j1+j2,j1+j2+1)
!
  ZalocEE(1:j1,1:j1) = STIFF_ALLE(1:j1,1:j1)
  ZalocEF(1:j1,1:j2) = STIFF_ALLE(1:j1,j1+1:j1+j2)
!
  ZalocFE(1:j2,1:j1) = STIFF_ALLE(j1+1:j1+j2,1:j1)
  ZalocFF(1:j2,1:j2) = STIFF_ALLE(j1+1:j1+j2,j1+1:j1+j2)
!

  iprint = 0
  if (iprint.ge.1) then
    write(*,7010)
7010   format('elem_primalMaxwell: ZblocE,ZblocF = ')
    write(*,7011) ZblocE(1:NrdofE)
    write(*,7011) ZblocF(1:NrdofE)
7011   format(10e12.5)
    call pause
    write(*,7012)
7012   format('elem_primalMaxwell: ZalocEE = ')
    do i=1,2*NrdofE
    write(*,7013) i,ZalocEE(i,1:NrdofE)
7013     format('i = ',i3,10(/,5(2e12.5,2x)))
    enddo
    call pause
    write(*,7014)
7014   format('elem_primalMaxwell: ZalocEF = ')
    do i=1,2*NrdofE
      write(*,7013) i,ZalocEF(i,1:NrdofE)
    enddo
    call pause
    write(*,7015)
7015   format('elem_primalMaxwell: ZalocFF = ')
    do i=1,6*NrdofQ
      write(*,7013) i,ZalocFF(i,1:NrdofE)
    enddo
    call pause
    endif


      end subroutine elem_primalMaxwell









