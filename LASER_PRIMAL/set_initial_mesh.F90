subroutine set_initial_mesh(Nelem_order)
!
      use GMP
      use parameters
      use data_structure3D
      use problem
!
      implicit none
      integer,dimension(NRELIS),intent(out) :: Nelem_order ! polynomial order for initial mesh elements
!
      integer :: iel,ndom,i,ifc,neig
      integer,dimension(6,2) :: ibc
!
!----------------------------------------------------------------------
!
!     check if have not exceeded the maximum order
      if (ORDER_APPROX.gt.MAXP) then
        write(*,*) 'set_initial_mesh: ORDER_APPROX, MAXP = ', ORDER_APPROX,MAXP
        stop
      endif
!
!     loop through initial mesh elements
      do iel=1,NRELIS
!
!  STEP 1 : set up order of approximation (elements can have different
!           orders in each direction)
        select case(ELEMS(iel)%Type)
        case('pris') ; Nelem_order(iel)=ORDER_APPROX*11
        case('bric') ; Nelem_order(iel)=ORDER_APPROX*111
        case('tetr') ; Nelem_order(iel)=ORDER_APPROX*1
        case('pyra') ; Nelem_order(iel)=ORDER_APPROX*1
        endselect
!
!  STEP 2 : set up physics
!
!       set up the number of physical attributes supported by the element
        ELEMS(iel)%nrphysics=2
        allocate(ELEMS(iel)%physics(2))
        ELEMS(iel)%physics(1)='Elfld'
        ELEMS(iel)%physics(2)='curnt'
!
!       initialize BC flags
!
!  .....single cube with Dirichlet, PEC/IBC
!!        call domain_number(iel, ndom)
        ibc = 0
        ! select case(iel)
        ! case(1)
        if (IBCFlag.eq.3) then
          ibc(1:6,1) = (/1,0,1,1,1,1/)
          ibc(1:6,2) = (/0,9,0,0,0,0/)
        else
           ibc(1:6,1) = (/1,1,1,1,1,1/)
        endif
        ! end select
!
!       allocate BC flags (1 per attribute)
        allocate(ELEMS(iel)%bcond(2))
!
!       for each attribute, encode face BC into a single BC flag
        call encodg(ibc(1:6,1),10,6, ELEMS(iel)%bcond(1))
        call encodg(ibc(1:6,2),10,6, ELEMS(iel)%bcond(2))
!
!     loop through initial mesh elements
      enddo
!
!
end subroutine set_initial_mesh
