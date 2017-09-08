c--------------------------------------------------------------------
c                                                                     
c     routine name      - elem
c                                                                     
c-------------------------------------------------------------------- 
c                                                                     
c     latest revision:  - Mar 15
c                                                                     
c     purpose:          - routine returns unconstrained (ordinary) 
c                         stiffness matrix and load vector  
c                         for the primal DPG method for Maxwell
c                         eqn
c                         
c                                                                    
c     arguments:                                                     
c                                                                     
c     in:              
c             Mdle      - an element middle node number, identified
c                         with the element
c             Maxdof    - column length of Zaloc
c     out:              
c             Nrdof     - number of dof for a single component
c             Itest_e,Itest_a,Itrial_e,Itrial_a - flags indicating
c                         presence of corresponding load vectors
c                         and stiffness matrices
c             Zbloc     - load vector
c             Zaloc     - stiffness matrix
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
c  ...optimal test problem
      case(3)
c
c  .....this is user-defined routine, you are in charge of it
        Itest(1:2)=1; Itrial(1:2)=1
        call elem_dpgHcurl(Mdle,BLOC(1)%nrow,BLOC(2)%nrow,
     .                  BLOC(1)%array,ALOC(1,1)%array,ALOC(1,2)%array,
     .                  BLOC(2)%array,ALOC(2,1)%array,ALOC(2,2)%array)
c
      case default
        write(*,*) 'elem: Mdle,NODES(Mdle)%case = ',
     .              Mdle,NODES(Mdle)%case
        stop1
      end select

c
c
      end subroutine elem  

c--------------------------------------------------------------------
c                                                                     
c     routine name      - elem_dpgHcurl
c                                                                     
c-------------------------------------------------------------------- 
c                                                                     
c     latest revision:  - Mar 15
c                                                                     
c     purpose:          - routine returns unconstrained (ordinary) 
c                         stiffness matrix and load vector  
c                         for the DPGH1 formulation for Laplace
c                         equation
c                                                                    
c     arguments:                                                     
c                                                                     
c     in:              
c             Mdle      - an element middle node number, identified
c                         with the element
c             MdE       - column length of ZalocEE,ZalocEF
c             MdF       - column length of ZalocFF,ZalocFE
c     out:              
c             ZblocE,ZblocF - load vectors
c             ZalocEE,ZalocEF,ZalocFE,ZalocFF - stiffness matrices
c
c---------------------------------------------------------------------
c
      subroutine elem_dpgHcurl(Mdle,MdE,MdF,
     .                         ZblocE,ZalocEE,ZalocEF,
     .                         ZblocF,ZalocFE,ZalocFF)   
c
      use control
      use parametersDPG
      use element_data
      use data_structure3D
      use DPGHcurl
      use problem
#include "syscom.blk"
c
      dimension ZblocE(MdE),ZalocEE(MdE,MdE),ZalocEF(MdE,MdF),
     .          ZblocF(MdF),ZalocFE(MdF,MdE),ZalocFF(MdF,MdF)
c
      character(len=4) :: etype,ftype
c
c  ...element order, orientation for edges and faces
      dimension norder(19),norient_edge(12),norient_face(6)
c
c  ...face order
      dimension norderf(5)
c
c  ...geometry dof
      dimension xnod(3,MAXbrickH)
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
c  ...3D quadrature data
      dimension xiloc(3,MAXNINT3ADD),waloc(MAXNINT3ADD)
c
c  ...2D quadrature data
      dimension tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
c
c  ...BC's flags
      dimension ibc(6,NR_PHYSA)
c
c  ...impressed current
      dimension zJ(3)
c
c  ...shape functions and their curl after Piola transform, rotated
c     tangential flux
      dimension E1(3),curlE1(3),E2(3),curlE2(3),rntimesE(3)
c
c  ...for debug printing
      dimension aux(10)
c
      character uplo
      nk(k1,k2) = (k2-1)*k2/2+k1
c
c---------------------------------------------------------------------
c
      iprint=0
      write(*,*) 'elem: Mdle = ',Mdle
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
      call find_orient(Mdle, norient_edge,norient_face)
c                                                                     
c  ...determine nodes coordinates 
      call nodcor(Mdle, xnod)
c
c  ...get the element boundary conditions flags
      call find_bc(Mdle, ibc)
      if (iprint.ge.1) then
        write(*,7001) Mdle
 7001   format('elem_dpgHcurl: BC FLAGS FOR Mdle = ',i5)
        do i=1,NR_PHYSA
          write(*,7002) PHYSA(i), ibc(1:nrf,i)
 7002     format('          ATTRIBUTE = ',a6,' FLAGS = ',6i2)
        enddo
      endif
c                                                                    
c  ...clear space for stiffness matrix and rhsv:                     
      ZblocE = ZERO; ZblocF = ZERO
      ZalocEE = ZERO; ZalocEF = ZERO; ZalocFE = ZERO; ZalocFF = ZERO
c
c  ...clear space for auxiliary matrices
      BLOADE = ZERO; STIFFTE = ZERO; STIFFTF = ZERO; AP = ZERO
c
c  ...complex wave number
      zk2 = OMEGA**2*EPSILON - ZI*OMEGA*SIGMA
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
c  .....get the RHS
        call getf(Mdle,x, zJ)
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
c  .......compute the RHS
          BLOADE(k1) = BLOADE(k1) 
     .      - zi*OMEGA*(zJ(1)*E1(1)+zJ(2)*E1(2)+zJ(3)*E1(3))*weight
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
c  .........accumulate for the test stiffness matrix
            k = nk(k1,k2)
            select case(INNER_PRODUCT)
            case(1)
              AP(k) = AP(k) 
     .              + (E1(1)*E2(1)+E1(2)*E2(2)+E1(3)*E2(3)
     .                + curlE1(1)*curlE2(1)+curlE1(2)*curlE2(2)
     .                + curlE1(3)*curlE2(3))*weight
            end select
          enddo
c
c  .......loop through Hcurl trial functions
          do k2=1,nrdofE
            E2(1:3) = shapE(1,k2)*dxidx(1,1:3) 
     .              + shapE(2,k2)*dxidx(2,1:3) 
     .              + shapE(3,k2)*dxidx(3,1:3) 
            curlE2(1:3) = dxdxi(1:3,1)*curlE(1,k2)
     .                  + dxdxi(1:3,2)*curlE(2,k2)
     .                  + dxdxi(1:3,3)*curlE(3,k2)
            curlE2(1:3) = curlE2(1:3)/rjac
c
c  .........accumulate for the extended stiffness matrix
            STIFFTE(k1,k2) = STIFFTE(k1,k2)
     .      + ((curlE1(1)*curlE2(1)+curlE1(2)*curlE2(2)
     .                   +curlE1(3)*curlE2(3))/MU
     .        -zk2*(E1(1)*E2(1)+E1(2)*E2(2)+E1(3)*E2(3)))*weight
           enddo
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'elem_dpgHcurl: AP = '
        do i=1,10
          do j=1,i-1
            aux(j) = AP(nk(j,i))
          enddo
          do j=i,10
            aux(j) = AP(nk(i,j))
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
          call shape3E(etype,xi,norder,norient_edge,norient_face, 
     .                 nrdofE,shapE,curlE)
c
c  .......geometry
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,
     .                 x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
          weight = bjac*wtloc(l)
c
c  .......loop through the enriched H(curl) test functions
          do k1=1,nrdofEE
            E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) 
     .              + shapEE(2,k1)*dxidx(2,1:3) 
     .              + shapEE(3,k1)*dxidx(3,1:3) 
c
c  .........loop through H(curl) trial functions
            do k2=1,nrdofE
              E2(1:3) = shapE(1,k2)*dxidx(1,1:3) 
     .                + shapE(2,k2)*dxidx(2,1:3) 
     .                + shapE(3,k2)*dxidx(3,1:3) 
              call cross_product(rn,E2, rntimesE)
c
c  ...........accumulate for the extended stiffness matrix
              STIFFTF(k1,k2) = STIFFTF(k1,k2)
     .       + (E1(1)*rntimesE(1)+E1(2)*rntimesE(2)+E1(3)*rntimesE(3))
     .          *weight
            enddo
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c
c  ...factorize the test stiffness matrix
      uplo = 'U'
      call ZPPTRF(uplo, nrdofEE, AP, info) 
      if (info.ne.0) then
        write(*,*) 'elem_dpgHcurl: info = ',info
        stop1
      endif
c
c  ...save copies of enriched stiffness matrices
      STIFFTEc =STIFFTE; STIFFTFc =  STIFFTF
c
c  ...compute the products of inverted test matrix with RHS
c     and enriched stiffness matrices
      call ZPPTRS(uplo, nrdofEE, 1, AP, BLOADE, MAXbrickEE, info1 )
      if (info1.ne.0) then
        write(*,*) 'elem_dpgHcurl: info1 = ',info1
        stop1
      endif
      call ZPPTRS(uplo,nrdofEE,nrdofE,AP,STIFFTEc,MAXbrickEE,info2)
      if (info2.ne.0) then
        write(*,*) 'elem_dpgHcurl: info2 = ',info2
        stop1
      endif
      call ZPPTRS(uplo,nrdofEE,nrdofE,AP,STIFFTFc,MAXbrickEE,info3)
      if (info3.ne.0) then
        write(*,*) 'elem_dpgHcurl: info3 = ',info3
        stop1
      endif
c
c  ...compute the ultimate DPG load vectors and stiffness matrices
      do k1=1,nrdofE
        do k=1,nrdofEE
          ZblocE(k1) = ZblocE(k1) + BLOADE(k)*conjg(STIFFTE(k,k1))
        enddo
        do k2=1,nrdofE
          do k=1,nrdofEE
            ZalocEE(k1,k2) = ZalocEE(k1,k2) 
     .                     + STIFFTEc(k,k1)*conjg(STIFFTE(k,k2))
          enddo
        enddo
        do k2=1,nrdofE
          do k=1,nrdofEE
            ZalocEF(k1,k2) = ZalocEF(k1,k2) 
     .                     + STIFFTEc(k,k1)*conjg(STIFFTF(k,k2))
          enddo
        enddo
      enddo
      do k1=1,nrdofE
        do k=1,nrdofEE
          ZblocF(k1) = ZblocF(k1) + BLOADE(k)*conjg(STIFFTF(k,k1))
        enddo
        do k2=1,nrdofE
          do k=1,nrdofEE
            ZalocFE(k1,k2) = ZalocFE(k1,k2) 
     .                     + STIFFTFc(k,k1)*conjg(STIFFTE(k,k2))
          enddo
        enddo
        do k2=1,nrdofE
          do k=1,nrdofEE
            ZalocFF(k1,k2) = ZalocFF(k1,k2) 
     .                     + STIFFTFc(k,k1)*conjg(STIFFTF(k,k2))
          enddo
        enddo
      enddo
c
c  ...check symmetry
      diffmax = ZERO; dmax = ZERO
      do k1=1,nrdofE
        do k2=k1,nrdofE
          diffmax = max(diffmax,
     .              abs(ZalocEE(k1,k2)-conjg(ZalocEE(k2,k1))))
          dmax = max(dmax,abs(ZalocEE(k1,k2)))
        enddo
      enddo
      symmetry_tol = 1.d-9
      if (diffmax/dmax.gt.symmetry_tol) then
        write(*,7021) diffmax, dmax
 7021   format('elem_dpgHcurl: diffmax,dmax FOR ZalocEE = ',2e12.5)
        call pause
      endif
      diffmax = ZERO; dmax = ZERO
      do k1=1,nrdofE
        do k2=k1,nrdofE
          diffmax = max(diffmax,
     .                  abs(ZalocFF(k1,k2)-conjg(ZalocFF(k2,k1))))
          dmax = max(dmax,abs(ZalocFF(k1,k2)))
        enddo
      enddo
      if (diffmax/dmax.gt.symmetry_tol) then
        write(*,7022) diffmax, dmax
 7022   format('elem_dpgHcurl: diffmax,dmax FOR ZalocFF = ',2e12.5)
        call pause
      endif
      diffmax = ZERO; dmax = ZERO
      do k1=1,nrdofE
        do k2=1,nrdofv
          diffmax = max(diffmax,
     .                  abs(ZalocEF(k1,k2)-conjg(ZalocFE(k2,k1))))
          dmax = max(dmax,abs(ZalocEF(k1,k2)))
        enddo
      enddo
      if (diffmax/dmax.gt.symmetry_tol) then
        write(*,7023) diffmax, dmax
 7023   format('elem_dpgHcurl: diffmax,dmax FOR ZalocEF = ',2e12.5)
        call pause
      endif
c
c
      if (iprint.ge.1) then
        write(*,7010) 
 7010   format('elem_dpgHcurl: ZblocE,ZblocF = ')
        write(*,7011) ZblocE(1:NrdofE)
        write(*,7011) ZblocF(1:NrdofE)
 7011   format(10e12.5)
        write(*,7012)
 7012   format('elem_dpgHcurl: ZalocEE = ')
        do i=1,NrdofE
          write(*,7013) i,ZalocEE(i,1:NrdofE)
 7013     format('i = ',i3,10(/,5(2e12.5,2x)))
        enddo
        write(*,7014)
 7014   format('elem_dpgHcurl: ZalocEF = ')
        do i=1,NrdofE
          write(*,7013) i,ZalocEF(i,1:NrdofE)
        enddo
        write(*,7015)
 7015   format('elem_dpgHcurl: ZalocFF = ')
        do i=1,NrdofE
          write(*,7013) i,ZalocFF(i,1:NrdofE)
        enddo 
        call pause
      endif
c
      end subroutine        









