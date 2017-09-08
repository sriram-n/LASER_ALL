
!  ...user defined routine to define the initial mesh
!
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
      integer,dimension(6,4) :: ibc
!
!----------------------------------------------------------------------
!
!  ...check if have not exceeded the maximum order
      if (ORDER_APPROX.gt.MAXP) then
        write(*,*) 'set_initial_mesh: ORDER_APPROX, MAXP = ', ORDER_APPROX,MAXP
        stop
      endif
!
!  ...loop through initial mesh elements
!      if (NRELIS.ne.1) then
!        write(*,*) 'set_initial_mesh: WRONG NRELIS = ',NRELIS
!        stop
!      endif
      !write(*,*) 'NRELIS is = ', NRELIS
      !call pause
      do iel=1,NRELIS
!
!  STEP 1 : set up order of approximation (elements may have different
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
        ELEMS(iel)%nrphysics=4
        allocate(ELEMS(iel)%physics(4))
        ELEMS(iel)%physics(1)='tempr'
        ELEMS(iel)%physics(2)='EHhat'
        ELEMS(iel)%physics(3)='hflux'
        ELEMS(iel)%physics(4)='EHfld'
!
!       initialize BC flags
		ibc = 0
      select case(GEOM_NO)
	!  .....single cube with Dirichlet, PEC BC on 5 out of 6 faces
  !  .....Impedance one one face
			case(1)
			  write(*,*) 'GEOM_NO = ', GEOM_NO
			  select case(iel)
			  case(1)
  			     ibc(1:6,1) = (/1,1,1,1,1,1/)
			  if(IBCFlag.eq.3) then
                ibc(1:6,2) = (/6,9,6,6,6,6/)
		      else
  			    ibc(1:6,2) = (/6,6,6,6,6,6/)
			  end if
			  end select

!  ..... fiber core prism -- modified input
			case (2)
			call domain_number(iel, ndom)
			  write(*,*) 'GEOM_NO = ', GEOM_NO
			  ibc(1:6,1) = (/1,1,1,0,0,0/)
		      if(IBCFlag.eq.3) then
                ibc(1:6,2) = (/6,9,6,0,0,0/)
		      else
  			    ibc(1:6,2) = (/6,6,6,0,0,0/)
			  end if
!  ..... fiber core prism -- old input
!			case (2)
!			call domain_number(iel, ndom)
!			  write(*,*) 'GEOM_NO = ', GEOM_NO
!			  ibc(1:6,1) = (/1,1,0,1,0,0/)
!		      if(IBCFlag.eq.3) then
!                ibc(1:6,2) = (/6,9,0,6,0,0/)
!		      else
!  			    ibc(1:6,2) = (/6,6,0,6,0,0/)
!			  end if

! ...... full fiber prism -- modified input
			case (3)
			call domain_number(iel, ndom)
			  write(*,*) 'GEOM_NO = ', GEOM_NO
			  select case(ndom)
			  case(1,2,3,4)
				ibc(1:6,1) = (/1,1,0,0,0,0/)
                if(IBCFlag.eq.3) then
                  ibc(1:6,2) = (/6,9,0,0,0,0/)
                else
                  ibc(1:6,2) = (/6,6,0,0,0,0/)
                end if
		      case(5,6,7,8)
		        ibc(1:6,1) = (/1,1,1,0,0,0/)
		        if(IBCFlag.eq.3) then
                  ibc(1:6,2) = (/6,9,6,0,0,0/)
		        else
  			      ibc(1:6,2) = (/6,6,6,0,0,0/)
			    end if
			endselect

! ..... full fiber prism - old input
!			case (3)
!			call domain_number(iel, ndom)
!			  write(*,*) 'GEOM_NO = ', GEOM_NO
!			  select case(ndom)
!			  case(1,2,3,4)
!				ibc(1:6,1) = (/1,1,0,0,0,0/)
!        if(IBCFlag.eq.3) then
!          ibc(1:6,2) = (/6,9,0,0,0,0/)
!        else
!        ibc(1:6,2) = (/6,6,0,0,0,0/)
!        end if
!		    case(5,6,7,8)
!		    ibc(1:6,1) = (/1,1,0,1,0,0/)
!		    if(IBCFlag.eq.3) then
!          ibc(1:6,2) = (/6,9,0,6,0,0/)
!		    else
!  			ibc(1:6,2) = (/6,6,0,6,0,0/)
!			  end if
!			  endselect

  !  ..... fiber core hexa
      case (4)
          call domain_number(iel, ndom)
          write(*,*) 'GEOM_NO = ', GEOM_NO
          select case(ndom)
          case(1)
			ibc(1:6,1) = (/1,1,0,0,0,0/)
			ibc(1:6,2) = (/6,9,0,0,0,0/)
          case(2,3,4,5)
          ibc(1:6,1) = (/1,1,1,0,0,0/)
          if(IBCFlag.eq.3) then
                ibc(1:6,2) = (/6,9,6,0,0,0/)
		      else
  			    ibc(1:6,2) = (/6,6,6,0,0,0/)
			  end if
          endselect


  !  ..... full fiber hexa
      case (5)
          call domain_number(iel, ndom)
          write(*,*) 'GEOM_NO = ', GEOM_NO
          select case(ndom)
          case(1,2,3,4,5)
			ibc(1:6,1) = (/1,1,0,0,0,0/)
			ibc(1:6,2) = (/6,9,0,0,0,0/)
          case(6,7,8,9)
          ibc(1:6,1) = (/1,1,1,0,0,0/)
          if(IBCFlag.eq.3) then
            ibc(1:6,2) = (/6,9,9,0,0,0/)
		  else
  	        ibc(1:6,2) = (/6,6,6,0,0,0/)
		  end if
          endselect

  !  .....single prism with Dirichlet, PEC BC on 5 out of 6 faces
  !  .....Impedance one one face
      case(6)
        write(*,*) 'GEOM_NO = ', GEOM_NO
        ibc(1:6,1) = (/1,1,1,1,1,0/)
        if(IBCFlag.eq.3) then
          ibc(1:6,2) = (/6,9,6,6,6,0/)
        else
          ibc(1:6,2) = (/6,6,6,6,6,0/)
        end if

!    ... end select GEOM_NO
     endselect



!
!       allocate BC flags (one per attribute)
        allocate(ELEMS(iel)%bcond(4))
!
!       for each attribute, encode face BC into a single BC flag
        do i=1,4
          call encodg(ibc(1:6,i),10,6, ELEMS(iel)%bcond(i))
        enddo
!
!  ...end of loop through initial mesh elements
      enddo
!
!
      end subroutine set_initial_mesh
