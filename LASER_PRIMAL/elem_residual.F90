!--------------------------------------------------------------------
!
!     routine name      - elem_residual
!
!--------------------------------------------------------------------
!
!     latest revision:  - August 17
!
!     purpose:          - routine returns element residual (squared)
!                         for the Primal formulation for Maxwell
!                         equations
!
!     arguments:
!
!     in:
!             Mdle      - an element middle node number, identified
!                         with the element
!     out:
!             Resid     - element residual (squared)
!             Nref_flag - suggested h-refinement flag
!
!---------------------------------------------------------------------
!
      subroutine elem_residual(Mdle, Resid,Nref_flag)
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
!
!
!.......declare input/output variables
  integer,                     intent(in)  :: Mdle
  integer,                     intent(out)  :: Nref_flag
  real*8,                      intent(out)  :: Resid
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
  VTYPE, dimension(MAXtestE) :: BLOADE,BLOADEc
!  ...stiffnes matrices for the enriched test space
!  VTYPE, dimension(MAXtestE,MAXbrickQ*6) :: STIFFEQ
  VTYPE, dimension(MAXtestE,MAXbrickE) :: STIFFEE
  VTYPE, dimension(MAXtestE,MAXbrickE) :: STIFFEF
!  ....STIFF_ALL for alternative computation of stiffness
  VTYPE, dimension(MAXtestE,2*MAXbrickE+1) :: STIFF_ALLE
#if C_MODE
  complex*16, allocatable :: DIAG_E(:)
#else
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

!  ...approximate solution
  VTYPE, dimension(3) :: zsolExi(3),zcurlExi(3),zsolE(3),zcurlE(3)
  VTYPE, dimension(3) :: zsolFxi(3),zsolF(3),zflux(3),zflux2(3)
  VTYPE, dimension(3) :: zcurl_exact(3),zflux_exact(3)
  VTYPE, dimension(3) :: zpsi(3),zcurl_xi_psi(3)
!  ...exact solution (for debugging)
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

      VTYPE :: zresid,Residual
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
  integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,nrTEST,nint,nint3,iflag,kE,k,iprint,l
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

!
!  ...directional contributions to element residual
!      dimension residd(0:4),zresidd(0:4),nref(3)
!
!  ...error representation function
!      dimension zpsi(3),zcurl_xi_psi(3)
!
!      dimension aux(10)
!
!---------------------------------------------------------------------
!
      select case(Mdle)
      case(1)
        iprint=0
      case default
        iprint=0
      end select
!
!  ...element type
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
!  ...determine order of approximation
      call find_order(Mdle, norder)
      norderc(1:nre+nrf) = norder(1:nre+nrf)
!
!  ...set the enriched order of appoximation
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
      call find_orient( Mdle, norient_edge,norient_face)
!
!  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
!
!  ...determine solution dof
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
      if (iprint.eq.1) then
        write(*,7020) xnod(1,1:nrv)
 7020   format('elem_residual: xnod  = ',8(f8.3,2x))
        write(*,7025) xnod(2,1:nrv)
        write(*,7025) xnod(3,1:nrv)
 7025   format('                       ',8(f8.3,2x))
        write(*,7030) 1,zdofE(1,1:nre)
        write(*,7030) 2,zdofE(2,1:nre)
 7030   format('elem_residual: zdofE(',i1',*) = ',2(/,6(2e12.5,2x)))
        call pause
      endif
!
!  ...clear space for auxiliary matrices
      BLOADE = ZERO; AP = ZERO
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
      call set_3Dint(etype,norder, nint3,xiloc,waloc)
      INTEGRATION = 0
      do l=1,nint3
        xi(1:3) = xiloc(1:3,l)
        wa = waloc(l)
!
!  .....determine element H1 shape functions
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
!  .....compute the approximate solution
        zsolExi = ZERO; zcurlExi(1:3) = ZERO
        do k=1,nrdofE
          zsolExi(1:3)  = zsolExi(1:3)  + zdofE(1,k)*shapE(1:3,k)
          zcurlExi(1:3) = zcurlExi(1:3) + zdofE(1,k)*curlE(1:3,k)
        enddo
        zsolE(1:3) = zsolExi(1)*dxidx(1,1:3) &
                  + zsolExi(2)*dxidx(2,1:3) &
                  + zsolExi(3)*dxidx(3,1:3)
        zcurlE(1:3) = dxdxi(1:3,1)*zcurlExi(1) &
                   + dxdxi(1:3,2)*zcurlExi(2) &
                   + dxdxi(1:3,3)*zcurlExi(3)
        zcurlE(1:3) = zcurlE(1:3)/rjac
!
        call exact(x,Mdle, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, &
                          zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
        zcurl_exact(1) = zdvalE(3,1,2)-zdvalE(2,1,3)
        zcurl_exact(2) = zdvalE(1,1,3)-zdvalE(3,1,1)
        zcurl_exact(3) = zdvalE(2,1,1)-zdvalE(1,1,2)

        if (iprint.eq.1) then
          write(*,8010) Mdle,x(1:3)
 8010     format('elem_residual: Mdle = ',i6,',  x = ',3f8.3)
          write(*,8020) zsolE(1:3), zcurlE(1:3)
 8020     format('  approximate solution and curl = ', &
                3(2e12.5,2x),3x,3(2e12.5,2x))
          write(*,8030) zvalE(1:3,1),zcurl_exact(1:3)
 8030     format('  exact solution and curl       = ', &
                3(2e12.5,2x),3x,3(2e12.5,2x))

        endif
!
!  .....debugging with the exact solution...
        ! if (iprint.eq.10) then
        !   zsolE(1:3)  = zvalE(1:3,1)
        !   zcurlE(1:3) = zcurl_exact(1:3)
        ! endif
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
!  .......accumulate for the load
          BLOADE(k1) = BLOADE(k1) &
           + ((curlE1(1)*zcurlE(1)+curlE1(2)*zcurlE(2)+curlE1(3)*zcurlE(3)) &
          -zk2*(E1(1)*zsolE(1)+E1(2)*zsolE(2)+E1(3)*zsolE(3)) &
         - (zJ(1)*E1(1)+zJ(2)*E1(2)+zJ(3)*E1(3)))*weight
!          write(*,*) 'l,k1,BLOADE(k1) = ', l,k1,BLOADE(k1)
!          if (k1.eq.nrdofEE) call pause
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
            select case(INNER_PRODUCT)
            case(1)
              AP(k) = AP(k) &
                   + (E1(1)*E2(1)+E1(2)*E2(2)+E1(3)*E2(3) &
                     + curlE1(1)*curlE2(1)+curlE1(2)*curlE2(2) &
                     + curlE1(3)*curlE2(3))*weight
            end select
          enddo
        enddo
      enddo
      if (iprint.eq.1) call pause
      if (iprint.eq.2) then
        do i=1,10
          do j=1,i-1
            zaux(j) = AP(nk(j,i))
          enddo
          do j=i,10
            zaux(j) = AP(nk(i,j))
          enddo
          write(*,7011) zaux
 7011   format(10e12.5)
        enddo
        call pause
      endif
      if (iprint.ge.1) then
        write(*,7014) BLOADE(1:nrdofEE)
 7014   format('elem_residual: BLOADE AFTER VOL INT = ', &
              10(/,6(2e12.5,2x)))
        call pause
      endif
!
!-----------------------------------------------------------------------
!
!  ...boundary integrals
!
!  ...loop through element faces
      do if=1,nrf
!
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
!  .......determine discontinuous H(curl) shape functions
          call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
!
!  .......determine element H1 shape functions (for geometry)
          call shape3H(etype,xi,norder,norient_edge,norient_face, &
                      nrdofH,shapH,gradH)
!
!  .......determine element H(curl) shape functions (for fluxes)
          call shape3E(etype,xi,norderc,norient_edge,norient_face, &
                      nrdofE,shapE,curlE)
!
!  .......geometry
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
          weight = bjac*wtloc(l)
!

!  .......compute approximate flux at the point
          zsolFxi = ZERO
          do k=1,nrdofE
!  .......if ibc compute E
            if(ibc(if,2).eq.9) then
              zsolFxi(1:3)  = zsolFxi(1:3)  + zdofE(1,k)*shapE(1:3,k)
            else
              zsolFxi(1:3)  = zsolFxi(1:3)  + zdofE(2,k)*shapE(1:3,k)
            endif
          enddo
          zsolF(1:3) = zsolFxi(1)*dxidx(1,1:3) &
                    + zsolFxi(2)*dxidx(2,1:3) &
                    + zsolFxi(3)*dxidx(3,1:3)
          call zcross_product(rn,zsolF, zflux)
!
!  .......if ibc compute n x n x E
          if(ibc(if,2).eq.9) then
            call zcross_product(rn,zflux, zflux2)
! ..........get the boundary source for ibc
            call get_bdSource(Mdle,x,rn, zImp)
          end if

          call exact(x,Mdle, &
                   zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, &
                   zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
          call zcross_product(rn,zvalE(1:3,2), zflux_exact)
          if (iprint.eq.1) then
            write(*,8040) if,x(1:3)
 8040       format('elem_residualPrimalMaxwell: if = ',i6,' x = ',3f8.3)
            write(*,8050) zflux
 8050       format('  approximate flux = ',3(2e12.5,2x))
            write(*,8060) zflux_exact
 8060       format('  exact flux       = ',3(2e12.5,2x))
          endif
!
!  .......debugging with the exact solution...
          ! if (iprint.eq.10) then
          !   zflux(1:3) = zflux_exact(1:3)
          ! endif
!
!  .......loop through enriched test functions
          do k1=1,nrdofEE
            E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                   + shapEE(2,k1)*dxidx(2,1:3) &
                   + shapEE(3,k1)*dxidx(3,1:3)
!
!  .........accumulate for the load vector
!  .......if ibc compute according to boundary load
          if(ibc(if,2).eq.9) then
            BLOADE(k1) = BLOADE(k1) &
           + (E1(1)*zflux2(1)+E1(2)*zflux2(2)+E1(3)*zflux2(3)) &
            *weight
          else
            BLOADE(k1) = BLOADE(k1) &
           + (E1(1)*zflux(1)+E1(2)*zflux(2)+E1(3)*zflux(3)) &
            *weight
!  ..... end ibc check
          end if
          enddo
        enddo
        if (iprint.eq.1) call pause
      enddo
      if (iprint.ge.1) then
        write(*,7015) BLOADE(1:nrdofEE)
 7015   format('elem_residualPrimalMaxwell: FINAL BLOADE = ',10(/,6(2e12.5,2x)))
        call pause
      endif
!
!-----------------------------------------------------------------------
!
!   allocate(DIAG_E(nrTest))
! ! ...preconditioning: getting diagonal entries
!   do k1=1,nrTEST
!     k = nk(k1,k1)
!     DIAG_E(k1) = AP_Maxwell(k)
!   enddo
!   do k1=1,nrTEST
!     do k2=k1,nrTEST
!       k = nk(k1,k2)
!       AP_Maxwell(k) = AP_Maxwell(k)/sqrt(DIAG_E(k1)*DIAG_E(k2))
!     enddo
!   enddo
!   deallocate(DIAG_E)
!  ...factorize the test stiffness matrix
      uplo = 'U'
      call ZPPTRF(uplo, nrdofEE, AP, info)
      if (info.ne.0) then
        write(*,*) 'elem_residualPrimalMaxwell: info = ',info
        stop
      endif
!
!  ...save copies of the RHS to compute later the residual
      BLOADEc = BLOADE
!
!  ...compute the product of inverted test Gram matrix with RHS,
!     BLOADE is overwritten with the solution
      call ZPPTRS(uplo, nrdofEE, 1, AP, BLOADE, MAXbrickEE, info1 )
      if (info1.ne.0) then
        write(*,*) 'elem_residualPrimalMaxwell: info1 = ',info1
        stop
      endif
!
!  ...compute the residual
      zresid = ZERO
      do k=1,NRDOFEE
        zresid = zresid + BLOADEc(k)*conjg(BLOADE(k))
      enddo
      Resid = zresid
      Nref_flag = 111
!
!-----------------------------------------------------------------------
! !
! !  ...recompute the element residual through direct integration to
! !     establish anisotropy flags
!       zresidd(0:3) = ZERO
!       do l=1,nint3
!         xi(1:3) = xiloc(1:3,l)
!         wa = waloc(l)
! !
! !  .....determine element H1 shape functions (for geometry)
!         call shape3H(etype,xi,norder,norient_edge,norient_face,
!      .               nrdofH,shapH,gradH)
! !
! !  .....determine discontinuous H(curl) shape functions
!         call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
! !
! !  .....geometry
!         call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,
!      .              x,dxdxi,dxidx,rjac,iflag)
! !
! !  .....integration weight
!         weight = rjac*wa
! !
! !  .....compute the error representation function
!         zpsi = ZERO; zcurl_xi_psi(1:3) = ZERO
!         do k=1,nrdofEE
!           zpsi(1:3) = zpsi(1:3) + BLOADE(k)*shapEE(1:3,k)
!           zcurl_xi_psi(1:3) = zcurl_xi_psi(1:3)
!      .                      + BLOADE(k)*curlEE(1:3,k)
!         enddo
! !
!         select case(INNER_PRODUCT)
!         case(1)
!           do i=1,3; do j=1,3
!             axi = 0.d0; bxi = 0.d0
!             do k=1,3
!               axi = axi + dxidx(i,k)*dxidx(j,k)
!               bxi = bxi + dxdxi(k,i)*dxdxi(k,j)
!             enddo
!             bxi = bxi/rjac**2
! !
! !  .........only the vector contributes to directional residuals
!             if (i.eq.j) then
!               zresidd(i) = zresidd(i) + axi*abs(zpsi(i))**2*weight
!             else
!               zresidd(0) = zresidd(0)
!      .                   + axi*zpsi(i)*conjg(zpsi(j))*weight
!             endif
!             zresidd(0) = zresidd(0)
!      .      + bxi*zcurl_xi_psi(i)*conjg(zcurl_xi_psi(j))*weight
!           enddo; enddo
!         end select
! !
! !  ...end of loop through integration points
!       enddo
!       residd(0:3) = zresidd(0:3)
!       diff = residd(0)+residd(1)+residd(2)+residd(3) - Resid
!       if (abs(diff).gt.1.d-8*abs(Resid)) then
!         write(*,*) 'Resid = ',Resid,
!      .              residd(0)+residd(1)+residd(2)+residd(3)
!         write(*,*) 'residd = ',residd(0:3)
!         call pause
!       endif
! !
! !  ...determine the refinement flag
!       select case(etype)
!       case('mdlb')
!         if (residd(0).lt..1d0*Resid) then
!           nref(1:3) = 1
!           do i=1,3
!             if (residd(i).lt..1d0*Resid) nref(i)=0
!           enddo
!           Nref_flag = nref(1)*100+nref(2)*10+nref(3)
!         else
!           Nref_flag = 111
!         endif
!       case('mdln','mdld')
!         Nref_flag = 1
!       case('mdlp')
!         if (residd(0).lt..1d0*Resid) then
!           nref(1:2) = 1
!           if (residd(1)+residd(2).lt..2d0*Resid) nref(1)=0
!           if (residd(3).lt..1d0*Resid) nref(2)=0
!           Nref_flag = nref(1)*10+nref(2)
!         else
!           Nref_flag = 111
!         endif
!       end select
! ccc      write(*,*) 'residd = ',residd(0:3)
! ccc      write(*,*) 'Mdle,Nref_flag = ',Mdle,Nref_flag
! ccc      call pause
!
!-----------------------------------------------------------------------
!
      if (iprint.eq.1) then
        write(*,7010) Mdle, Resid
 7010   format('elem_residualPrimalMaxwell: Mdle, Resid = ',i5,3x,e12.5)
        call pause
      endif
!
      end subroutine elem_residual

!-----------------------------------------------------------------------
!
!     routine name      - compute_residual
!
!--------------------------------------------------------------------
!
!     latest revision:  - August 17
!
!     purpose:          - routine returns global residual
!                         for the Primal formulation of Maxwell
!
!---------------------------------------------------------------------
!
      subroutine compute_residual
!
      use data_structure3D
      use environment      , only : QUIET_MODE
      implicit none
#if C_MODE
#define VTYPE  complex*16
#else
#define VTYPE double precision
#endif
!
!  ...visitation flag
      integer, save :: ivis = 0
!
!  ...total number of dof for the old mesh
      integer, save :: nrdof_total_old
!
!  ...residual for the old mesh
      real*8,  save :: residual_old
!
!  ...residuals and rates to display
      real*8 , dimension(2,10), save :: rwork
!
!  ...number of dof to display
      integer , dimension(10), save :: iwork
      integer :: iprint,nrdof_total,nref_flag,mdle,iel,ndofH,ndofV,ndofE,ndofQ,i
      real*8 :: residual,resid,rate

!
      iprint=0
!
!  ...compute total residual and number of dof
      nrdof_total = 2*NRDOFSE; residual = 0.d0
      mdle = 0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_residual(mdle, resid,nref_flag)
        if (iprint.eq.1) then
          write(*,7010) iel, mdle, resid
 7010     format('compute_residual: iel, mdle = ',2i8, &
                ' element residual = ',e12.5)
        endif
        residual = residual + resid
        call find_ndof(mdle, ndofH,ndofE,ndofV,ndofQ)
!
!  .....subtract the middle node H(curl) dof from the global count
        nrdof_total = nrdof_total - ndofE
      enddo
      residual = sqrt(residual)
!
!  ...compute rate
      rate = 0.d0
      if (ivis.ne.0) then
        if (nrdof_total.gt.nrdof_total_old) then
          rate = log(residual_old/residual) &
               /log(float(nrdof_total_old)/float(nrdof_total))
        endif
      endif
!
!  ...save current data
      ivis = ivis+1
      nrdof_total_old = nrdof_total
      residual_old = residual
!
!  ...store data to display
      iwork(ivis) = nrdof_total
      rwork(1,ivis) = residual; rwork(2,ivis) = rate
!
!  ...display the convergence history
      if (.NOT. QUIET_MODE) then
        write(*,*)''
        write(*,*)'         -- Error Report --'
        write(*,7100)
 7100   format(' Mesh  ','  Nrdof  ', ' Residual   ','     Rate ')
        do i=1,ivis
          write(*,7110)i,iwork(i),rwork(1:2,i)
 7110     format(i3,4x,i7,2x,e12.5,2x,f8.3)
        enddo
        write(*,*)''
      endif
!
!
      end subroutine compute_residual






