c--------------------------------------------------------------------
c
c     routine name      - elem_residual
c
c--------------------------------------------------------------------
c
c     latest revision:  - Mar 17
c
c     purpose:          - routine returns element residual (squared)
c                         for the Primal Poisson and UW Time Harmonic
c                         Maxwell equation
c
c     arguments:
c
c     in:
c             Mdle      - an element middle node number, identified
c                         with the element
c     out:
c             Resid     - element residual (squared)
c             Nref_flag - suggested h-refinement flag
c
c---------------------------------------------------------------------
c
      subroutine elem_residual(Mdle, Resid,Nref_flag)
c
      use control
      use parametersDPG
      use element_data
      use data_structure3D
      use problem
#include "syscom.blk"
c
c
      character(len=4) :: etype,ftype
c
c  ...element order, orientation for edges and faces
      dimension norder(19),norient_edge(12),norient_face(6)
      dimension norderc(19)
c
c  ...face order
      dimension norderf(5)
c
c  ...geometry dof
      dimension xnod(3,MAXbrickH)
c
c  ...solution dof
      dimension zdofH(MAXEQNH,MAXbrickH),
     .          zdofE(MAXEQNE,MAXbrickE),
     .          zdofV(MAXEQNV,MAXbrickV),
     .          zdofQ(MAXEQNQ,MAXbrickQ)
c
c  ...geometry
      dimension xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3),
     .          dxdt(3,2),rt(3,2),rn(3),t(2)
c
c  ...H1 shape functions
      dimension shapH(MAXbrickH),gradH(3,MAXbrickH)
c
c  ...Hdiv shape functions
      dimension shapV(3,MAXbrickV),divV(MAXbrickV)

c  ...Hcurl shape functions
      dimension shapE(3,MAXbrickE),curlE(3,MAXbrickE)
c
c  ...L2 shape functions
      dimension shapQ(MAXbrickQ)
c
c
c  ...3D quadrature data
      dimension xiloc(3,MAXNINT3ADD),waloc(MAXNINT3ADD)
c
c  ...2D quadrature data
      dimension tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
c
c  ...derivatives wrt physical coordinates, flux
      dimension dv1(3),dv2(3)
c
c  ...approximate solution
      dimension zgradHxi(3),zgradH(3),
     .          zsolVxi(3), zsolV(3)
c
c  ...directional contributions to element residual
      dimension residd(0:4),nref(3)
c
c  ...error representation function
      dimension gradpsi(3)
c
      dimension aux(10)
c  ...impressed current
      dimension zJ(3)
c
c  ...shape functions and their curl after Piola transform, rotated
c     tangential flux
      dimension E1(3),curlE1(3),E2(3),curlE2(3),rntimesE(3)
c
c  ...approximate solution
      dimension zsolExi(3,2),zsolE(3,2),zflux(3,2),
     .          zsolQ(6),
     .          zflux_exact(3,2)
c
c  ...directional contributions to element residual
      dimension zresidd(0:4)
c
c  ...error representation function
      dimension zpsi_xi(3,2),zcurl_xi_psi(3,2),
     .          zpsi(3,2),zcurl_psi(3,2)
c
c  ...metric for L2 norm
      dimension aa(3,3)

c  ...exact solution (for debugging)
      dimension  zvalH(  MAXEQNH    ),
     .          zdvalH(  MAXEQNH,3  ),
     .         zd2valH(  MAXEQNH,3,3),
     .           zvalE(3,MAXEQNE    ),
     .          zdvalE(3,MAXEQNE,3  ),
     .         zd2valE(3,MAXEQNE,3,3),
     .           zvalV(3,MAXEQNV    ),
     .          zdvalV(3,MAXEQNV,3  ),
     .         zd2valV(3,MAXEQNV,3,3),
     .           zvalQ(  MAXEQNQ    ),
     .          zdvalQ(  MAXEQNQ,3  ),
     .         zd2valQ(  MAXEQNQ,3,3)
c
      character uplo
      nk(k1,k2) = (k2-1)*k2/2+k1
c
c---------------------------------------------------------------------
c
c---------------SELECT THE CASE OF SINGLE STEP HEAT OR MAXWELL--------
c---------------------------------------------------------------------
      select case(NO_PROBLEM)
c  .....single step of transient heat equation (1)
c  .....multiple steps of transient heat equation (2)
      case(1,2)

      select case(Mdle)
      case(1)
        iprint=0
      case default
        iprint=0
      end select
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
      case('mdlb'); nordP = NODES(Mdle)%order+NORD_ADD*111
      case('mdln','mdld'); nordP = NODES(Mdle)%order+NORD_ADD
      case('mdlp'); nordP = NODES(Mdle)%order+NORD_ADD*11
      end select
c
c  ...determine edge and face orientations
      call find_orient( Mdle, norient_edge,norient_face)
c
c  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
c
c  ...determine solution dof
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
      if (iprint.eq.1) then
        write(*,7020) xnod(1,1:8),xnod(2,1:8),xnod(3,1:8)
 7020   format('elem_residual: xnod  = ',8(f8.3,2x),
     . 2(  /,'                       ',8(f8.3,2x)))
        write(*,7030) zdofH(1,1:8),zdofV(1,1:6)
 7030   format('elem_residual: zdofH = ',8(e12.5,2x),
     .     /,'               zdofV = ',6(e12.5,2x))
      endif
c
c  ...clear space for auxiliary matrices
      BLOADH = ZERO; AP_Heat = ZERO
c
c-----------------------------------------------------------------------
c
c  ...element integrals...
c
c  ...use the enriched order to set the quadrature
      INTEGRATION = NORD_ADD
      call set_3Dint(etype,norder, nint3,xiloc,waloc)
      INTEGRATION = 0
      do l=1,nint3
        xi(1:3) = xiloc(1:3,l)
        wa = waloc(l)
c
c  .....determine element H1 shape functions
        call shape3H(etype,xi,norder,norient_edge,norient_face,
     .               nrdofH,shapH,gradH)
c
c  .....determine discontinuous H1 shape functions
        call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
c
c  .....geometry
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,
     .              x,dxdxi,dxidx,rjac,iflag)
c
c  .....integration weight
        weight = rjac*wa
c
c  .....compute the approximate solution
        zsolH = 0.d0; zgradHxi(1:3) = 0.d0
        do k=1,nrdofH
          zsolH = zsolH + zdofH(1,k)*shapH(k)
          zgradHxi(1:3) = zgradHxi(1:3) + zdofH(1,k)*gradH(1:3,k)
        enddo
        zgradH(1:3) = zgradHxi(1)*dxidx(1,1:3)
     .              + zgradHxi(2)*dxidx(2,1:3)
     .              + zgradHxi(3)*dxidx(3,1:3)
c
        call exact(x,Mdle, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE,
     .                     zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
        if (iprint.eq.1) then
          write(*,8010) Mdle,x(1:3)
 8010     format('elem_residual: Mdle = ',i6,',  x = ',3f8.3)
          write(*,8020) zsolH, zgradH(1:3)
 8020     format('  approximate solution and gradient = ',
     .           e20.12,3x,3(e20.12,2x))
          write(*,8030) zvalH(1),zdvalH(1,1:3)
 8030     format('  exact solution and gradient       = ',
     .           e20.12,3x,3(e20.12,2x))
        endif
c
c  .....get the RHS
        call getf(Mdle,x, zfval,zJ)
c
c  .....loop through enriched H1 test functions
        do k1=1,nrdofHH
          v1 = shapHH(k1)
          dv1(1:3) = gradHH(1,k1)*dxidx(1,1:3)
     .             + gradHH(2,k1)*dxidx(2,1:3)
     .             + gradHH(3,k1)*dxidx(3,1:3)
c
c
c  .......single step of heat equation
            BLOADH(k1) = BLOADH(k1)
     .      + (KAPPA*DELTAT*(zgradH(1)*dv1(1)+zgradH(2)*dv1(2)+
     .                       zgradH(3)*dv1(3))
     .        +zsolH*v1-zfval*v1)*weight
c
c  .......loop through enriched H1 trial functions
          do k2=k1,nrdofHH
            v2 = shapHH(k2)
            dv2(1:3) = gradHH(1,k2)*dxidx(1,1:3)
     .               + gradHH(2,k2)*dxidx(2,1:3)
     .               + gradHH(3,k2)*dxidx(3,1:3)
c
c  .........accumulate for the test stiffness matrix
            k = nk(k1,k2)
            select case(INNER_PRODUCT)
            case(1)
              AP_Heat(k) = AP_Heat(k)
     .              + (dv1(1)*dv2(1)+dv1(2)*dv2(2)+dv1(3)*dv2(3)
     .                +v1*v2)*weight
            end select
          enddo
        enddo
      enddo
      if (iprint.eq.2) then
        do i=1,10
          do j=1,i-1
            aux(j) = AP_Heat(nk(j,i))
          enddo
          do j=i,10
            aux(j) = AP_Heat(nk(i,j))
          enddo
          write(*,7011) aux
 7011   format(10e12.5)
        enddo
        call pause
      endif
      if (iprint.ge.1) then
        write(*,7014) BLOADH(1:nrdofHH)
 7014   format('elem_residual: BLOADH AFTER VOL INT = ',
     .         20(/,10(e12.5,2x)))
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...boundary integrals
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
c
ccc P. Gatto, Jan 15 : call to new shape functions
ccc          call shape3V_1(etype,xi,norder,norient_face,
c
          call shape3V(etype,xi,norder,norient_face,
     .                 nrdofV,shapV,divV)
c
c  .......geometry
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,
     .                 x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
          weight = bjac*wtloc(l)
c
c  .......compute approximate flux at the point
          zsolVxi(1:3) = 0.d0
          do k=1,nrdofV
            zsolVxi(1:3) = zsolVxi(1:3) + zdofV(1,k)*shapV(1:3,k)
          enddo
          zsolV(1:3) = (dxdxi(1:3,1)*zsolVxi(1)
     .                 +dxdxi(1:3,2)*zsolVxi(2)
     .                 +dxdxi(1:3,3)*zsolVxi(3))/rjac
          zsolVn = ZsolV(1)*rn(1)+ZsolV(2)*rn(2)+ZsolV(3)*rn(3)
          call exact(x,Mdle,
     .              zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE,
     .              zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
          zvalVn = zvalV(1,1)*rn(1)+zvalV(2,1)*rn(2)+zvalV(3,1)*rn(3)
          if (iprint.eq.1) then
            write(*,8040) if,x(1:3)
 8040       format('elem_residual: if = ',i6,' x = ',3f8.3)
            write(*,8050) zsolVn
 8050       format('  approximate flux = ',e20.12)
            write(*,8060) zvalVn
 8060       format('  exact flux       = ',e20.12)
          endif
c
c  .......loop through enriched test functions
          do k1=1,nrdofHH
            v1 = shapHH(k1)
c
c  .........accumulate for the load vector
            BLOADH(k1) = BLOADH(k1) - zsolVn*v1*weight
ccc            BLOADH(k1) = BLOADH(k1) - zvalVn*v1*weight
          enddo
        enddo
        if (iprint.eq.1) call pause
      enddo
      if (iprint.ge.1) then
        write(*,7015) BLOADH(1:nrdofHH)
 7015   format('elem_residual: FINAL BLOADH = ',10(/,10(e12.5,2x)))
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...factorize the test stiffness matrix
      uplo = 'U'
      call ZPPTRF(uplo, nrdofHH, AP_Heat, info)
      if (info.ne.0) then
        write(*,*) 'elem_dpgH1: info = ',info
        stop1
      endif
c
c  ...save copies of the RHS to compute later the residual
      BLOADHc = BLOADH
c
c  ...compute the product of inverted test Gram matrix with RHS,
c     BLOADH is overwritten with the solution
      call ZPPTRS(uplo, nrdofHH, 1, AP_Heat, BLOADH, MAXbrickHH, info1 )
      if (info1.ne.0) then
        write(*,*) 'elem_dpgH1: info1 = ',info1
        stop1
      endif
c
c  ...compute the residual
      Resid = 0.d0
      do k=1,NRDOFHH
        Resid = Resid + BLOADHc(k)*BLOADH(k)
      enddo
c
c-----------------------------------------------------------------------
c
c  ...recompute the elemnt residual through direct integration to
c     establish anistropy flags
      residd(0:3) = 0.d0
      do l=1,nint3
        xi(1:3) = xiloc(1:3,l)
        wa = waloc(l)
c
c  .....determine element H1 shape functions (for geometry)
        call shape3H(etype,xi,norder,norient_edge,norient_face,
     .               nrdofH,shapH,gradH)
c
c  .....determine discontinuous H1 shape functions
        call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
c
c  .....geometry
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,
     .              x,dxdxi,dxidx,rjac,iflag)
c
c  .....integration weight
        weight = rjac*wa
c
c  .....compute the error representation function
        psi = 0.d0; gradpsi(1:3) = 0.d0
        do k=1,nrdofHH
          psi          = psi          + BLOADH(k)*shapHH(k)
          gradpsi(1:3) = gradpsi(1:3) + BLOADH(k)*gradHH(1:3,k)
        enddo
c
        select case(INNER_PRODUCT)
        case(1)
          do i=1,3; do j=1,3
            axi = 0.d0
            do k=1,3
              axi = axi + dxidx(i,k)*dxidx(j,k)
            enddo
            if (i.eq.j) then
              residd(i) = residd(i) + axi*gradpsi(i)**2*weight
            else
              residd(0) = residd(0) + axi*gradpsi(i)*gradpsi(j)*weight
            endif
          enddo; enddo
          residd(0) = residd(0) + psi**2*weight
        end select
c
c  ...end of loop through integration points
      enddo
      diff = residd(0)+residd(1)+residd(2)+residd(3) - Resid
      if (abs(diff).gt.1.d-8*abs(Resid)) then
        write(*,*) 'Resid = ',Resid,
     .              residd(0)+residd(1)+residd(2)+residd(3)
        write(*,*) 'residd = ',residd(0:3)
        call pause
      endif
c
c  ...determine the refinement flag
      select case(etype)
      case('mdlb')
        if (residd(0).lt..1d0*Resid) then
          nref(1:3) = 1
          do i=1,3
            if (residd(i).lt..1d0*Resid) nref(i)=0
          enddo
          Nref_flag = nref(1)*100+nref(2)*10+nref(3)
        else
          Nref_flag = 111
        endif
      case('mdln','mdld')
        Nref_flag = 1
      case('mdlp')
        if (residd(0).lt..1d0*Resid) then
          nref(1:2) = 1
          if (residd(1)+residd(2).lt..2d0*Resid) nref(1)=0
          if (residd(3).lt..1d0*Resid) nref(2)=0
          Nref_flag = nref(1)*10+nref(2)
        else
          Nref_flag = 111
        endif
      end select
c      write(*,*) 'residd = ',residd(0:3)
c      write(*,*) 'Mdle,Nref_flag = ',Mdle,Nref_flag
ccc      call pause
c
c-----------------------------------------------------------------------
c
      if (iprint.ge.1) then
        write(*,7010) Mdle, Resid
 7010   format('elem_residual: Mdle, Resid = ',i5,3x,e12.5)
        call pause
      endif
c
c  ...... time harmonic Maxwell (3)
      case(3)


      select case(Mdle)
      case(1)
        iprint=0
      case default
        iprint=0
      end select
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
      call find_orient( Mdle, norient_edge,norient_face)
c
c  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
c
c  ...determine solution dof
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
      if (iprint.eq.1) then
        write(*,9020) xnod(1,1:nrv)
 9020   format('elem_residual: xnod  = ',8(f8.3,2x))
        write(*,9025) xnod(2,1:nrv)
        write(*,9025) xnod(3,1:nrv)
 9025   format('                       ',8(f8.3,2x))
        write(*,9030) 1,zdofE(1,1:nre)
        write(*,9030) 2,zdofE(2,1:nre)
 9030   format('elem_residual: zdofE(',i1',*) = ',2(/,6(2e12.5,2x)))
        call pause
      endif
c
c  ...clear space for auxiliary matrices
      BLOADE = ZERO; AP_Maxwell = ZERO
c
c  ...shortcut: adjust frequency
c
c  ...element size:
      h = min(abs(xnod(1,2)-xnod(1,1)),abs(xnod(2,4)-xnod(2,1)),
     .        abs(xnod(3,5)-xnod(3,1)))
      omeg = min(OMEGA,6.d0/h)
c
c-----------------------------------------------------------------------
c
c  ...element integrals...
c
c  ...use the enriched order to set the quadrature
      INTEGRATION = NORD_ADD
      call set_3Dint(etype,norder, nint3,xiloc,waloc)
      INTEGRATION = 0
      do l=1,nint3
        xi(1:3) = xiloc(1:3,l)
        wa = waloc(l)
c
c  .....determine element H1 shape functions
        call shape3H(etype,xi,norder,norient_edge,norient_face,
     .               nrdofH,shapH,gradH)
c
c  .....determine element H(curl) shape functions
        call shape3E(etype,xi,norder,norient_edge,norient_face,
     .               nrdofE,shapE,curlE)
c
c  .....determine element L2 shape functions
        call shape3Q(etype,xi,norder, nrdofQ,shapQ)
c
c  .....determine discontinuous H(curl) shape functions
        call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
c
c  .....geometry
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,
     .              x,dxdxi,dxidx,rjac,iflag)
c
c  .....integration weight
        weight = rjac*wa
c
c  .....compute the approximate solution
        zsolQ = ZERO
        do k=1,nrdofQ
          zsolQ(1:6)  = zsolQ(1:6)  + zdofQ(1:6,k)*shapQ(k)
        enddo
        zsolQ = zsolQ/rjac
c
        call exact(x,Mdle, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE,
     .                     zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)

        if (iprint.eq.1) then
          write(*,3010) Mdle,x(1:3)
 3010     format('elem_residual: Mdle = ',i6,',  x = ',3f8.3)
          write(*,3020) zsolQ(1:6)
 3020     format('  approximate solution  = ', 6(2e12.5,2x))
          write(*,3030) zvalQ(1:6)
 3030     format('  exact solution        = ', 6(2e12.5,2x))
        endif
c
c  .....debugging with the exact solution...
        if (iprint.eq.10) then
          zsolQ(1:6)  = zvalQ(1:6)
        endif
c
c  .....get the RHS
        call getf(Mdle,x, zfval,zJ)
c
c  .....loop through enriched H(curl) test functions
        do k1=1,nrdofEE
          E1(1:3) = shapEE(1,k1)*dxidx(1,1:3)
     .            + shapEE(2,k1)*dxidx(2,1:3)
     .            + shapEE(3,k1)*dxidx(3,1:3)
          curlE1(1:3) = dxdxi(1:3,1)*curlEE(1,k1)
     .                + dxdxi(1:3,2)*curlEE(2,k1)
     .                + dxdxi(1:3,3)*curlEE(3,k1)
          curlE1(1:3) = curlE1(1:3)/rjac
c
c  .......accumulate for the load
          k = 2*k1
          BLOADE(k) = BLOADE(k)
     .      + (zJ(1)*E1(1)+zJ(2)*E1(2)+zJ(3)*E1(3))*weight
c
          k = 2*k1-1
          BLOADE(k) = BLOADE(k)
     .    - ((curlE1(1)*zsolQ(1)+curlE1(2)*zsolQ(2)+curlE1(3)*zsolQ(3))
     .       +ZI*OMEGA*MU*(E1(1)*zsolQ(4)+E1(2)*zsolQ(5)+E1(3)*zsolQ(6))
     .      )*weight
c
          k = 2*k1
          BLOADE(k) = BLOADE(k)
     .    - ((curlE1(1)*zsolQ(4)+curlE1(2)*zsolQ(5)+curlE1(3)*zsolQ(6))
     .       -(ZI*OMEGA*EPSILON+SIGMA)
     .        *(E1(1)*zsolQ(1)+E1(2)*zsolQ(2)+E1(3)*zsolQ(3)))*weight
c
c  .......loop through enriched H(curl) trial functions
          do k2=k1,nrdofEE
            E2(1:3) = shapEE(1,k2)*dxidx(1,1:3)
     .              + shapEE(2,k2)*dxidx(2,1:3)
     .              + shapEE(3,k2)*dxidx(3,1:3)
            curlE2(1:3) = dxdxi(1:3,1)*curlEE(1,k2)
     .                  + dxdxi(1:3,2)*curlEE(2,k2)
     .                  + dxdxi(1:3,3)*curlEE(3,k2)
            curlE2(1:3) = curlE2(1:3)/rjac
c
c  .........axiliary quantities
            CC = curlE1(1)*curlE2(1) + curlE1(2)*curlE2(2)
     .         + curlE1(3)*curlE2(3)
            EE = E1(1)*E2(1) + E1(2)*E2(2) + E1(3)*E2(3)
            CE = curlE1(1)*E2(1) + curlE1(2)*E2(2) + curlE1(3)*E2(3)
            EC = E1(1)*curlE2(1) + E1(2)*curlE2(2) + E1(3)*curlE2(3)
c
c  .........accumulate for the test stiffness matrix
            k = nk(k1,k2)
            select case(INNER_PRODUCT)
            case(1)
              k = nk(2*k1-1,2*k2-1)
              AP_Maxwell(k) = AP_Maxwell(k)
     .              + (CC + (abs(za)**2 + 1.d0)*EE)*weight
c     .              + (CC + (omeg**2*MU**2 + 1.d0)*EE)*weight
              k = nk(2*k1-1,2*k2  )
              AP_Maxwell(k) = AP_Maxwell(k)
     .              + (-ZI*omeg*MU*CE - conjg(za)*EC)*weight
cc     .              + ((ZI*omeg*EPSILON-SIGMA)*CE+ZI*omeg*MU*EC)*weight
              if (k1.ne.k2) then
                k = nk(2*k1  ,2*k2-1)
                AP_Maxwell(k) = AP_Maxwell(k)
     .                + (-za*CE + ZI*omeg*MU*EC)*weight
cc     .             - ((ZI*omeg*EPSILON+SIGMA)*EC + ZI*omeg*MU*CE)*weight
              endif
              k = nk(2*k1  ,2*k2  )
              AP_Maxwell(k) = AP_Maxwell(k)
     .         + (CC + ((omeg*MU)**2 + 1.d0)*EE)*weight
cc.         + (CC + (abs(ZI*omeg*EPSILON-SIGMA)**2 + 1.d0)*EE)*weight
            end select
          enddo
        enddo
      enddo
      !if (iprint.eq.1) call pause
      if (iprint.eq.2) then
        do i=1,10
          do j=1,i-1
            aux(j) = AP_Maxwell(nk(j,i))
          enddo
          do j=i,10
            aux(j) = AP_Maxwell(nk(i,j))
          enddo
          write(*,9011) aux
 9011   format(5(2e12.5,2x))
        enddo
        call pause
      endif
      if (iprint.ge.1) then
        write(*,9014) BLOADE(1:2*nrdofEE)
 9014   format('elem_residual: BLOADE AFTER VOL INT = ',
     .         10(/,6(2e12.5,2x)))
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...boundary integrals
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
c  .......determine discontinuous H(curl) shape functions
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
c
c  .......compute approximate fluxes at the point
          zsolExi = ZERO
          do ivar=1,2
            do k=1,nrdofE
              zsolExi(1:3,ivar)  = zsolExi(1:3,ivar)
     .                           + zdofE(ivar,k)*shapE(1:3,k)
            enddo
c
            zsolE(1:3,ivar) = zsolExi(1,ivar)*dxidx(1,1:3)
     .                      + zsolExi(2,ivar)*dxidx(2,1:3)
     .                      + zsolExi(3,ivar)*dxidx(3,1:3)
            call zcross_product(rn,zsolE(1:3,ivar), zflux(1:3,ivar))
          enddo
          call exact(x,Mdle,
     .              zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE,
     .              zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
          call zcross_product(rn,zvalQ(1:3), zflux_exact(1:3,1))
          call zcross_product(rn,zvalQ(4:6), zflux_exact(1:3,2))
          if (iprint.eq.1) then
            write(*,5040) if,x(1:3)
 5040       format('elem_residual: if = ',i6,' x = ',3f8.3)
            write(*,5050) zflux
 5050       format('  approximate flux = ',2(3(2e12.5,2x),4x))
            write(*,5060) zflux_exact
 5060       format('  exact flux       = ',2(3(2e12.5,2x),4x))
          endif
c
c  .......debugging with the exact solution...
          if (iprint.eq.10) then
            zflux = zflux_exact
          endif
c
c  .......loop through enriched test functions
          do k1=1,nrdofEE
            E1(1:3) = shapEE(1,k1)*dxidx(1,1:3)
     .              + shapEE(2,k1)*dxidx(2,1:3)
     .              + shapEE(3,k1)*dxidx(3,1:3)
c
c  .........accumulate for the load vector
            k = 2*k1-1
            BLOADE(k) = BLOADE(k)
     .      - (E1(1)*zflux(1,1)+E1(2)*zflux(2,1)+E1(3)*zflux(3,1))
     .       *weight
            k = 2*k1
            BLOADE(k) = BLOADE(k)
     .      - (E1(1)*zflux(1,2)+E1(2)*zflux(2,2)+E1(3)*zflux(3,2))
     .       *weight
          enddo
        enddo
        if (iprint.ne.0) call pause
      enddo
      if (iprint.gt.0) then
        write(*,5015) BLOADE(1:2*nrdofEE)
 5015   format('elem_residual: FINAL BLOADE = ',10(/,6(2e12.5,2x)))
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...factorize the test stiffness matrix
      uplo = 'U'
      nrtest = nrdofEE*2
      call ZPPTRF(uplo, nrtest, AP_Maxwell, info)
      if (info.ne.0) then
        write(*,*) 'elem_residual: info = ',info
        stop1
      endif
c
c  ...save copies of the RHS to compute later the residual
      BLOADEc = BLOADE
c
c  ...compute the product of inverted test Gram matrix with RHS,
c     BLOADE is overwritten with the solution
      call ZPPTRS(uplo, nrtest, 1, AP_Maxwell, BLOADE, MAXtestE, info1 )
      if (info1.ne.0) then
        write(*,*) 'elem_residual: info1 = ',info1
        stop1
      endif
c
c  ...compute the residual
      zresid = ZERO
      do k=1,nrtest
        zresid = zresid + BLOADEc(k)*conjg(BLOADE(k))
      enddo
      Resid = zresid
c
c-----------------------------------------------------------------------
c
c  ...recompute the element residual through direct integration to
c     establish anisotropy flags
      zresidd(0:3) = ZERO
      do l=1,nint3
        xi(1:3) = xiloc(1:3,l)
        wa = waloc(l)
c
c  .....determine element H1 shape functions (for geometry)
        call shape3H(etype,xi,norder,norient_edge,norient_face,
     .               nrdofH,shapH,gradH)
c
c  .....determine discontinuous H(curl) shape functions
        call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
c
c  .....geometry
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,
     .              x,dxdxi,dxidx,rjac,iflag)
c
c  .....integration weight
        weight = rjac*wa
c
c  .....metric for L2 norm
        aa = 0.d0
        do i=1,3; do j=1,3
          do k=1,3
            aa(i,j) = aa(i,j) + dxidx(i,k)*dxidx(j,k)
          enddo
        enddo; enddo
c
c  .....compute the error representation function
        zpsi_xi = ZERO; zcurl_xi_psi = ZERO
        do k=1,nrdofEE
          zpsi_xi(1:3,1) = zpsi_xi(1:3,1) + BLOADE(2*k-1)*shapEE(1:3,k)
          zpsi_xi(1:3,2) = zpsi_xi(1:3,2) + BLOADE(2*k  )*shapEE(1:3,k)
          zcurl_xi_psi(1:3,1) = zcurl_xi_psi(1:3,1)
     .                        + BLOADE(2*k-1)*curlEE(1:3,k)
          zcurl_xi_psi(1:3,2) = zcurl_xi_psi(1:3,2)
     .                        + BLOADE(2*k  )*curlEE(1:3,k)
        enddo
c
c  .....Piola transforms
        zpsi = ZERO; zcurl_psi = ZERO
        do ivar=1,2
          do j=1,3
            zpsi(1:3,ivar) = zpsi(1:3,ivar)
     .                     + zpsi_xi(j,ivar)*dxidx(j,1:3)
            zcurl_psi(1:3,ivar) = zcurl_psi(1:3,ivar)
     .                     + dxdxi(1:3,j)*zcurl_xi_psi(j,ivar)
          enddo
        enddo
        zcurl_psi = zcurl_psi/rjac
c
        select case(INNER_PRODUCT)
        case(1)
          do k=1,3
            zresidd(0) = zresidd(0)
     .      + abs(zcurl_psi(k,1)+(ZI*omeg*EPSILON-SIGMA)*zpsi(k,2))**2
     .        *weight + abs(zpsi(k,1))**2*weight
          enddo
          do k=1,3
            zresidd(0) = zresidd(0)
     .      + abs(zcurl_psi(k,2)-(ZI*omeg*MU)*zpsi(k,1))**2
     .        *weight + abs(zpsi(k,2))**2*weight
          enddo
c
c  .......the only directional contributions are coming from L2 product
          do k=1,3
            zresidd(k) = zresidd(k) + aa(k,k)*(
     .      (abs(OMEGA*MU)**2 + 1.d0)*abs(zpsi_xi(k,1))**2
     .     +(abs(ZI*omeg*EPSILON-SIGMA)**2 + 1.d0)*abs(zpsi_xi(k,2))**2)
     .       *weight
          enddo
        end select
c
c  ...end of loop through integration points
      enddo
      zresidd(0) = zresidd(0) - (zresidd(1)+zresidd(2)+zresidd(3))
      residd(0:3) = zresidd(0:3)
      diff = residd(0)+residd(1)+residd(2)+residd(3) - Resid
C       write(*,*) 'diff = ',diff
C       write(*,*) '1.d-8*abs(Resid) = ',1.d-8*abs(Resid)
C       call pause
C       if (abs(diff).gt.1.d-8*abs(Resid)) then
C         write(*,*) 'Resid = ',Resid,
C      .              residd(0)+residd(1)+residd(2)+residd(3)
C         write(*,*) 'residd = ',residd(0:3)
C         call pause
C       endif
c
c  ...determine the refinement flag
      select case(etype)
      case('mdlb')
        if (residd(0).lt..1d0*Resid) then
          nref(1:3) = 1
          do i=1,3
            if (residd(i).lt..1d0*Resid) nref(i)=0
          enddo
          Nref_flag = nref(1)*100+nref(2)*10+nref(3)
        else
          Nref_flag = 111
        endif
      case('mdln','mdld')
        Nref_flag = 1
      case('mdlp')
        if (residd(0).lt..1d0*Resid) then
          nref(1:2) = 1
          if (residd(1)+residd(2).lt..2d0*Resid) nref(1)=0
          if (residd(3).lt..1d0*Resid) nref(2)=0
          Nref_flag = nref(1)*10+nref(2)
        else
          Nref_flag = 111
        endif
      end select
ccc      write(*,*) 'residd = ',residd(0:3)
ccc      write(*,*) 'Mdle,Nref_flag = ',Mdle,Nref_flag
ccc      call pause
c
c-----------------------------------------------------------------------
c
      if (iprint.eq.1) then
        write(*,6010) Mdle, Resid
 6010   format('elem_residual: Mdle, Resid = ',i5,3x,e12.5)
        call pause
      endif
c......end select for NO_PROBLEM
      end select
c
      end subroutine elem_residual

c-----------------------------------------------------------------------
c
c     routine name      - compute_residual
c
c--------------------------------------------------------------------
c
c     latest revision:  - Mar 17
c
c     purpose:          - routine returns global residual
c                         for the Primal Poisson and UW Time Harmonic
c                         Maxwell equation
c
c---------------------------------------------------------------------
c
      subroutine compute_residual
c
      use data_structure3D
      use environment      , only : QUIET_MODE
      use problem
#include "syscom.blk"
c
c  ...visitation flag
      integer, save :: ivis = 0
c
c  ...total number of dof for the old mesh
      integer, save :: nrdof_total_old
c
c  ...residual for the old mesh
      real*8,  save :: residual_old
c
c  ...residuals and rates to display
      real*8 , dimension(2,10), save :: rwork
c
c  ...number of dof to display
      integer , dimension(10), save :: iwork

c
      iprint=0
c
c  ...compute total residual and number of dof
c  .... either Heat or Maxwell
      select case (NO_PROBLEM)
c  .... single step of Heat equation
      case(1,2)
      nrdof_total = NRDOFSH+NRDOFSV; residual = 0.d0
c  .... Maxwell case
      case(3)
      nrdof_total = 2*NRDOFSE; residual = 0.d0
      end select
      mdle = 0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_residual(mdle, resid,nref_flag)
        if (iprint.eq.1) then
          write(*,7010) iel, mdle, resid
 7010     format('compute_residual: iel, mdle = ',2i5,
     .           ' element residual = ',e12.5)
        endif
        residual = residual + resid
        call find_ndof(mdle, ndofH,ndofE,ndofV,ndofQ)
c
c  .....subtract the middle node H(div) dof from the global count
        nrdof_total = nrdof_total - ndofV
      enddo
      residual = sqrt(residual)
c
c  ...compute rate
      rate = 0.d0
      if (ivis.ne.0) then
        if (nrdof_total.gt.nrdof_total_old) then
          rate = log(residual_old/residual)
     .          /log(float(nrdof_total_old)/float(nrdof_total))
        endif
      endif
c
c  ...save current data
      ivis = ivis+1
      nrdof_total_old = nrdof_total
      residual_old = residual
c
c  ...store data to display
      iwork(ivis) = nrdof_total
      rwork(1,ivis) = residual; rwork(2,ivis) = rate
c
c  ...display the convergence history
      if (.NOT. QUIET_MODE) then
        write(*,*)''
        write(*,*)'         -- Error Report --         '
        write(*,*)'         PROBLEM:                   '
        select case(NO_PROBLEM)
        case(1,2)
        write(*,*)'         Heat equation              '
        case(3)
        write(*,*)'         Time harmonic Maxwell      '
        end select
        write(*,7100)
 7100   format(' Mesh  ','  Nrdof  ', ' Residual   ','     Rate ')
        do i=1,ivis
          write(*,7110)i,iwork(i),rwork(1:2,i)
 7110     format(i3,4x,i7,2x,e12.5,2x,f8.3)
        enddo
        write(*,*)''
      endif
c
c
      end subroutine compute_residual






