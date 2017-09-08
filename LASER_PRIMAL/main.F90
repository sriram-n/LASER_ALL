!---------------------------------------------------------
!  The code is run through a bash script
!---------------------------------------------------------
program main
!
      use environment
      use control
      use paraview
      use parametersDPG
      use GMP
      use data_structure3D
      use physics
      use uhm2
      use problem
!
      implicit none
      character(len=1024) :: argv
      real*8  :: err,rnorm,rvoid,factor,t
      integer :: mdle,i,kref,idec,nvoid,niter,ibeg,iend,nreflag,istep,nrdof_old,nrdof_new,nstop,idec_solve
      integer :: iref,numRef,itag,isolver,iso_ans,info,ref_xyz,iii
      integer, dimension(2) :: flag,iflag
!
!----------------------------------------------------------------------
!
!      L2PROJ = .TRUE.
!  ...initialize environment
      call begin_environment
!
!  ...read in HP3D input files location (if option is not present, the default value is used)
!
!                             option label      // explanation                // default value     // parameter
      call get_option_string( '-file-control'    , 'Control file'              , './files/control'  , FILE_CONTROL)
      call get_option_string( '-file-geometry'   , 'Geometry file'             , './files/cube', FILE_GEOM   )
      call get_option_string( '-file-phys'       , 'Physics file'              , './files/physics'  , FILE_PHYS   )
      call get_option_string( '-file-refinement' , 'Refinement files location' , '../../../files/ref'  , FILE_REFINE )
      call get_option_string( '-file-history'    , 'History file'              , './files/history'  , FILE_HISTORY)
      call get_option_string( '-file-err'        , 'Error file'                , './files/dump_err' , FILE_ERR    )
      call get_option_string( '-prefix'          , 'Prefix paraview file'      ,'laser1_'            , PREFIX      )
!
!     -- Parview Interface --
      ! Variables relevant to src/modules/paraview
!                        option label     // explanation                        // default value     // parameter
      call get_option_string('-file-vis-upscale','Visualization upscale file location','../../../files/vis',FILE_VIS          )
      call get_option_string('-vis-level'       ,'Visualization upscale level (0-3)'  ,'2'                 ,VLEVEL            )
      call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'./output/figures/test/' ,PARAVIEW_DIR      )
      call get_option_bool(  '-paraview-geom'   ,'Dump geom at every Paraview call'   ,.TRUE.              ,PARAVIEW_DUMP_GEOM)
      call get_option_bool(  '-paraview-attr'   ,'Dump solution to Paraview'          ,.TRUE.              ,PARAVIEW_DUMP_ATTR)
!
!  ...read in problem dependent parameters (defined in module parametersDPG,DPGH1)
!
!                              option label      // explanation                // default value     // parameter
      call get_option_int(    '-nord-add'           , 'NORD_ADD'                  , 2                  , NORD_ADD    )
      call get_option_int(    '-order-approx'       , 'ORDER_APPROX'              , 2                  , ORDER_APPROX)
      call get_option_int(    '-orderx'             , 'NPX'                       , 0                  , NPX         )
      call get_option_int(    '-ordery'             , 'NPY'                       , 0                  , NPY         )
      call get_option_int(    '-orderz'             , 'NPZ'                       , 1                  , NPZ         )
      call get_option_int(    '-isol'               , 'ISOL'                      , 1                  , ISOL        )
      call get_option_int(    '-comp'               , 'ICOMP_EXACT'               , 2                  , ICOMP_EXACT )
!
      call get_option_int(    '-inner-product'      , 'INNER_PRODUCT'             , 1                  , INNER_PRODUCT)
      call get_option_real(   '-mu'                 , 'MU'                        , 1.d0               , MU          )
      call get_option_real(   '-epsilon'            , 'EPSILON'                   , 1.d0               , EPSILON     )
      call get_option_real(   '-sigma'              , 'SIGMA'                     , 0.d0               , SIGMA       )
      call get_option_real(   '-omega'              , 'OMEGA'                     , 1.d0               , OMEGA       )
      !call get_option_real(   '-omega'              , 'OMEGA'                     , sqrt(5.d0)*PI       , OMEGA       )
       call get_option_int(    '-ibc'               , 'IBCFlag'                   , 0                  , IBCFlag     )
      GAMMA_IMP = 1.d0
      !GAMMA_IMP = dsqrt(1.d0-(PI**2/OMEGA**2))

!
!  ...finalize
      call end_environment
!
!  ...print fancy header
      write(*,*)'                      '
      write(*,*)'// --  PRIMAL DPG FOR MAXWELL EQUATION -- //'
      write(*,*)'                      '
!
!  ...initialize problem
      call initialize

!     Kyungjoo's magic...
      UHM_VERBOSE            = .FALSE.
      UHM_SOLVER_REPORT_SHOW = .FALSE.
!
      call get_command(argv)
      call uhm_initialize(argv)
!
      call uhm_option_begin
      call uhm_option_end
!
      call uhm_time_begin(UHM_TIME_LEVEL_FRONT)
      call uhm_direct_lu_piv_initialize( &
!              UHM_DOUBLE, NR_RHS_PROB, 256, UHM_SOLVER_PTR)
              UHM_DOUBLE, MY_NR_RHS, 256, UHM_SOLVER_PTR)
!
!  ...display menu in infinite loop
 10   continue
      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      write(*,*) 'Quit ...................................0'
      write(*,*) '                                         '
      write(*,*) 'Geometry graphics (X11) ................1'
      write(*,*) 'HP3D graphics (X11) ....................3'
      write(*,*) 'Paraview ...............................4'
      write(*,*) 'Print Data Structure arrays ............5'
      write(*,*) 'Dumpout Data Structure .................7'
      write(*,*) '                                         '
      write(*,*) ' --  Geometry & Refinements  --          '
      write(*,*) 'Geometry error ........................12'
      write(*,*) 'Interactive H-refinements .............31'
      write(*,*) 'Uniform     H-refinements .............32'
      write(*,*) 'Adaptive    H-refinements .............33'
      write(*,*) 'Anisotropic H-refinements .............34'
      write(*,*) 'H-refinements of Prism Core Faces......35'
      write(*,*) '                                         '
      write(*,*) 'Solve (frontal) .......................40'
      write(*,*) 'Solve (MUMPS) .........................50'
      write(*,*) 'Solve (UHM) ...........................60'
      write(*,*) '                                         '
      write(*,*) 'Compute Hcurl error ..................100'
      write(*,*) 'Compute residual .....................110'
      write(*,*) 'Adaptive DPG refinements .............120'
      write(*,*) 'Compute BC data interpolation error...130'
      write(*,*) 'Rate tests ...........................140'
      write(*,*) '                                         '
      write(*,*) 'My tests..... ........................200'
      write(*,*) 'My own tests..... ....................210'
      write(*,*) 'My Size test..... ....................220'
      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
!
      read(*,*) idec
!
!----------------------------------------------------------------------
!
      select case(idec)
!  ...quit
      case( 0) ; call finalize ; stop
!
!  ...GMP graphics
      case( 1) ; call graphg
!
!  ...hp3d graphics
      case( 3) ; call graphb

!  Paraview graphics
        case(4) ; call paraview_driver(iParAttr)
!
!  ...print data structure
      case( 5) ; call result
!
!  ...dump out data structure
      case( 7)
        write(*,*)'dumping out GMP...'
        call dumpout_GMP
        write(*,*)'dumping out physics...'
        call dumpout_physics
        write(*,*)'dumping out HP3D...'
        call dumpout_hp3d('./files/dumpc3Dhp')
!
!  ...geometry error
      case(12) ; call geom_error(err,rnorm)
!
!  ...interactive refinements
      case(31)
        write(*,*)'Active elements:'
        mdle=0
        do i=1,NRELES
          call nelcon(mdle, mdle)
!
          select case(NODES(mdle)%type)
          case('mdlp') ; write(*,7000) mdle
          case('mdlb') ; write(*,7001) mdle
          case('mdln') ; write(*,7002) mdle
          case('mdld') ; write(*,7003) mdle
          endselect
 7000     format(' mdle = ',i6,' ; PRISM')
 7001     format(' mdle = ',i6,' ; BRICK')
 7002     format(' mdle = ',i6,' ; TET')
 7003     format(' mdle = ',i6,' ; PYRAMID')
!
        enddo
        call display_ref_kind
        write(*,7010)
 7010   format(' mdle,kref =')
        read(*,*) mdle,kref
!
!       refine element
        call refine(mdle,kref)
!       recover 1-irregular mesh, update geometry and Dirichlet dof's
        call close_mesh ; call update_gdof ; call update_ddof
!
!     uniform global H-refinements
      case(32)
!       refine elements
        call global_href
!       recover 1-irregular mesh, update geometry and Dirichlet dof's
        call close_mesh ; call update_gdof ; call update_Ddof
        !  .....propagate the impedance flag
        call propagate_flag(2,9)
!
!     adaptive H-refinements
      case(33) ; call adaptivity_geom(0.3d0, nvoid,rvoid)

!  ..... anisotropic h-refinement of cylinder
        case(34)
        write(*,*)'set number of anisotropic refinements'
        write(*,*)'-1 for infinite loop'
        read(*,*), iso_ans
        info=0
        if(iso_ans.ne.-1) then
        !ref_xyz=1
        !call setAnisoRef(info,ref_xyz, info)
        ref_xyz=2
        do iii=1,iso_ans
          call setAnisoRef(info,ref_xyz, info)
          call propagate_flag(2,9)
        enddo
        !call update_gdof
        !call update_Ddof
        else
        !ref_xyz=1
        !call setAnisoRef(info,ref_xyz, info)
        ref_xyz=2
        do while(1.gt.0)
          call setAnisoRef(info,ref_xyz, info)
          call propagate_flag(2,9)
        enddo
        !call update_gdof
        !call update_Ddof
        endif

!  ..... anisotropic h-refinement of prism core face only
        case(35)
        write(*,*)'set number of uniform refinements of prism core only'
        write(*,*)'-1 for infinite loop'
        read(*,*), iso_ans
        info=0
        if(iso_ans.ne.-1) then
        ref_xyz=3
        do iii=1,iso_ans
          call setAnisoRef(info,ref_xyz, info)
        enddo
        call update_gdof
        call update_Ddof
        call propagate_flag(2,9)
        else
        ref_xyz=3
        do while(1.gt.0)
          call setAnisoRef(info,ref_xyz, info)
        enddo
        call update_gdof
        call update_Ddof
        call propagate_flag(2,9)
        endif
!
!     frontal solve
      case(40)
        !call hack_solver
        call solve1(MY_NR_RHS)
        !call unhack_solver
!
!     MUMPS solve
      case(50)
        !call hack_solver
        call mumps_sc('H')
        !call mumps_interf(MY_NR_RHS)
        !call unhack_solver
        !call mumps_solve_seq(MY_NR_RHS)
!
!     UHM solve
      case(60)
        call uhm_solve
        call uhm_solver_flush(UHM_SOLVER_PTR)
!
!     compute Hcurl error for the E-field only
      case(100)
        flag(1)=1; flag(2)=0
        call compute_error(flag,1)
!
!  ...compute residual
      case(110)
        call compute_residual
!
!  ...adaptive DPG refinements
      case(120)
  333   write(*,7011)
 7011   format('main: SET INITIAL, FINAL STEP,', &
               ' REFINEMENT FLAG (1-href,2-pref,3-hpref), FACTOR, idec_solve')
        read(*,*) ibeg,iend,nreflag,factor,idec_solve
        if (nreflag.ne.1.and.nreflag.ne.2.and.nreflag.ne.3) go to 333
        istep=ibeg-1
        do while(istep.lt.iend)
          istep=istep+1
          nrdof_old = NRDOFSE
          !!!!!!!call adapt_DPG(idec_solve,istep,nreflag,factor,nstop)
          if (nstop.eq.1) exit
          nrdof_new = NRDOFSE
          if (nrdof_new.eq.nrdof_old) then
            istep=istep-1
            factor = factor*0.25d0
            write(*,7023) factor
 7023       format('main: NEW factor = ',e12.5)
            if (factor.lt.0.000001d0) exit
          endif
        enddo
!
!     compute BC data interpolation error
      case(130)
        !!!!call compute_BC_interp_error
!


! ......... rate test for single step of Heat problem
        case(140)
        write(*,*) 'Testing h-convergence rates'
        write(*,*) 'Enter number of uniform H-refinements:'
        read(*,*) numRef
  55 continue
        write(*,*) 'MUMPS or PARDISO: MUMPS = 1, PARDISO = 2'
        read(*,*) isolver

          iref = 0
          do while(iref.lt.numRef)
            write(*,*) '1) isolver is :', isolver
!  ........ Solve the problem
!            call solve1(MY_NR_RHS)
             !call mumps_interf(MY_NR_RHS)
            call uhm_time_in
            !call hack_solver
            if(isolver.eq.1) then
              call mumps_sc('H')
            elseif (isolver.eq.2) then
              call pardiso_sc_herm
            else
               write(*,*) 'main PRIMAL_LASER: error: enter 1 or 2'
               goto 55
            endif
            !call unhack_solver
            !call mumps_solve_seq(MY_NR_RHS)
!  ........ Set error flags
            iflag(1) = 1; iflag(2) = 0; itag=1
!  ........ Compute error
            call compute_error(iflag,itag)
            call compute_residual
            call uhm_time_out(t)
!  ......... Compute residual
            !call compute_residual
!  ........ Do uniform h-refinements
            call global_href
            call close_mesh
            call update_gdof
            call update_Ddof
!  .....propagate the impedance flag
            call propagate_flag(2,9)
            iref = iref+1
          enddo





      case(200) ; call my_tests
!
      case(210) ; call my_own_tests

      case(220) ; call my_sizetest
      endselect
!
!  ...go back to menu
      goto 10
!
!
end program main

      subroutine my_tests
!
      implicit none
      integer :: idec,mdle,kref
!
      write(*,*) 'SELECT:'
      write(*,*) 'RETURN............................................0'
      write(*,*) 'PERFORM h-REFINEMENT OF A SINGLE ELEMENT..........1'
!
      read(*,*) idec
      select case(idec)
      case(0); return
      case(1)
        write(*,*) 'SET mdle,kref '
        read(*,*) mdle,kref
        write(*,*) 'mdle,kref = ',mdle,kref
        call refine(mdle,kref)
        call close_mesh ; call update_gdof ; call update_ddof
      end select
!
      end subroutine my_tests




      subroutine my_sizetest
        use data_structure3D
        implicit none
        integer :: nod
        integer*8 :: size_coord, size_h1, size_hcurl, size_hdiv, size_l2, size_tot
        size_coord = 0
        size_h1 = 0
        size_hcurl = 0
        size_hdiv = 0
        size_l2 = 0

         do nod=1,NRNODS
           !if(NODES(nod)%act.eq.1) then
             if (associated(NODES(nod)%coord).eq. .true.) then
               size_coord = size_coord + STORAGE_SIZE(NODES(nod)%coord)
             endif
             if (associated(NODES(nod)%zdofH).eq. .true.) then
               size_h1 = size_h1 + STORAGE_SIZE(NODES(nod)%zdofH)
             endif
             if (associated(NODES(nod)%zdofE).eq. .true.) then
               size_hcurl = size_hcurl + STORAGE_SIZE(NODES(nod)%zdofE)
             endif
             if (associated(NODES(nod)%zdofV).eq. .true.) then
               size_hdiv = size_hdiv + STORAGE_SIZE(NODES(nod)%zdofV)
             endif
             if (associated(NODES(nod)%zdofQ).eq. .true.) then
               size_l2 = size_l2 + STORAGE_SIZE(NODES(nod)%zdofQ)
             endif
          !endif
         size_tot = size_coord + size_h1 + size_hcurl + size_hdiv + size_l2
        enddo

         write(*,*) 'total DOF size is: ', size_tot
         call pause

      end subroutine my_sizetest

!
!  ...this is a hack to eliminate middle node H(curl) dof for the flux in celem_system by using BC flag
      subroutine hack_solver
!
      use physics, only: NR_PHYSA
      use data_structure3D
!
      implicit none
      integer :: iel,mdle
      integer, dimension(NR_PHYSA) :: nbcond
!
      nbcond(1)=0; nbcond(2)=1
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call encod(nbcond,10,NR_PHYSA, NODES(Mdle)%bcond)
        call set_index(NODES(Mdle)%case,NODES(Mdle)%bcond, NODES(Mdle)%index)
      enddo
!
      end subroutine hack_solver
!
      subroutine unhack_solver
!
      use physics, only: NR_PHYSA
      use data_structure3D
!
      implicit none
      integer :: iel,mdle
      integer, dimension(NR_PHYSA) :: nbcond
!
!---------------------------------------------------------------------------------------------------------------
!
!  ...this is a hack to eliminate middle node H(curl) dof for the flux in celem_system by using BC flag
      nbcond(1)=0; nbcond(2)=0
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call encod(nbcond,10,NR_PHYSA, NODES(Mdle)%bcond)
        call set_index(NODES(Mdle)%case,NODES(Mdle)%bcond, NODES(Mdle)%index)
      enddo
!
      end subroutine unhack_solver
!----------------------------------------------------------------------
!
!   routine name       - propagate_flag
!
!----------------------------------------------------------------------
!
!   latest revision    - Sep 17
!
!   purpose            - propagate Nflag from element faces to element
!                        edges and vertices; the flag is passed
!                        provided ALL adjacent faces share the flag
!
!   arguments:
!     in :
!              Icomp   - physics attribute number
!              Nflag   - BC flag
!
!----------------------------------------------------------------------

      subroutine propagate_flag(Icomp,Nflag)
!
      use data_structure3D
      implicit none
      integer :: Icomp,Nflag
!
      character(len=4) :: type
      integer :: iel,mdle,if,nrfn,i,j,nod,iprint
!
!  ...element nodes and orientations, face nodes
      integer :: nodesl(27),norientl(27),nface_nodes(9)
!
!  ...element face BC flags, decoded BC flag for a node
      integer :: ibc(6,NR_PHYSA), nodflag(NR_PHYSA)
!
      iprint=0
!
!  ...loop through active elements
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        type = NODES(mdle)%type
!
!  .....determine element nodes
        call elem_nodes(mdle, nodesl,norientl)
!
!  .....get the element boundary conditions flags
        call find_bc(mdle, ibc)
!
!  .....loop through element faces
        do if=1,nface(type)
!
!  .......determine face node numbers
          call face_nodes(type,if, nface_nodes,nrfn)
          if (iprint.eq.1) then
            write(*,7010) mdle,if, ibc(if,Icomp)
 7010       format('propagate_flag: mdle,if,ibc(if,Icomp) = ',i6,i3,i4)
          endif
!
!  .......loop through the face nodes
          do i=1,nrfn-1
            j = nface_nodes(i)
            nod = nodesl(j)
!
!  .........propagate the flag unless prohibited
            if (ibc(if,Icomp).eq.Nflag) then
              if (NODES(nod)%visit.ne.-Nflag) NODES(nod)%visit = Nflag
!
!  .........prohibit the flag to be passed to the node
            else
              NODES(nod)%visit = -Nflag
            endif
          enddo
        enddo
      enddo
!
!  ...change -Nflag to zero
      do nod=1,NRNODS
        if (iprint.eq.1) then
          write(*,7020) nod,NODES(nod)%visit
 7020     format('propagate_flag: nod,NODES(nod)%visit = ',i8,i3)
        endif
        if (NODES(nod)%visit.eq.0) cycle
        call decod(NODES(nod)%bcond,10,NR_PHYSA, nodflag)
        if (NODES(nod)%visit.eq.-Nflag) then
          nodflag(Icomp)=0
        elseif (NODES(nod)%visit.eq.Nflag) then
          nodflag(Icomp)=Nflag
        endif
        call encod(nodflag,10,NR_PHYSA, NODES(nod)%bcond)
!
!  .....reset the node index
        call set_index(NODES(nod)%case,NODES(nod)%bcond, NODES(nod)%index)
      enddo
      call reset_visit
!
      end subroutine propagate_flag
subroutine setAnisoRef(info,ref_xyz, info1)
! This routine anisotropically refines
! prism and hexa elements
! info = 0 if success
! ref_xyz: refinement flag to set kref
!          1 for xy plane direction
!          2 for z direction
! info = 1 otherwise
#include "implicit_none.h"
    !
    use control
    use data_structure3D
    use uhm2
    use parameters, only : ZERO,ZONE,ZIMG
      implicit none
      integer,                       intent(in)  :: info
      integer,                       intent(in)  :: ref_xyz
      integer,                       intent(out)  :: info1
      integer, allocatable, dimension(:) :: list_elem
      integer :: i,iprint,ic,mdle,iel,kref,nr_elem_to_refine
      write(*,*) 'From IANISOREF before anything'
        !call pause

        allocate (list_elem(NRELES))
        write(*,*) 'NRELES is: ', NRELES
        !call pause
        ic=0
        mdle=0
        do iel=1,NRELES
            call nelcon(mdle, mdle)
            ic=ic+1
            list_elem(ic) = mdle
        enddo
        nr_elem_to_refine = NRELES
        if (nr_elem_to_refine.gt.0) then
            !      ...refine the elements from the list
            do iel=1,nr_elem_to_refine
                mdle = list_elem(iel)
                select case(NODES(mdle)%type)
!               PRISM
                case('mdlp')
                  if(ref_xyz.eq.2) then
                    kref=01
                  elseif(ref_xyz.eq.1) then
                    kref=10
                  elseif(ref_xyz.eq.3) then
                    kref=10
                  else
                    write(*,*) 'error from IANISOREF: invalid ref_xyz in prism'
                    stop
                  endif
!
!               BRICK
                case('mdlb')
                  if(ref_xyz.eq.2) then
                    kref=001
                  elseif(ref_xyz.eq.1) then
                    kref=010
                  elseif(ref_xyz.eq.3) then
                    kref=000
                  else
                    write(*,*) 'error from IANISOREF: invalid ref_xyz in bric'
                    stop
                  endif
                end select
                call refine(mdle,kref)
            enddo
            !      ...close the mesh
            call close
            !      ...update geometry and Dirichlet flux dof after the refinements
            !call update_gdof
            !call update_Ddof
        endif
        info1 = 1
!
!
end subroutine setAnisoRef
