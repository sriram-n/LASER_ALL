!---------------------------------------------------------
!  The code is run through a bash script
!---------------------------------------------------------
program main
!
      use environment
      use paraview
      use control
      use parametersDPG
      use GMP
      use data_structure3D
      use physics
      use uhm2
      use problem
      use assembly
      use matrices
!
      implicit none
      character(len=1024) :: argv
      real*8  :: err,rnorm,rvoid,factor,t,geo_err,geo_rnorm,Ehat,gainF,xtest
      integer :: mdle,i,kref,idec,nvoid,niter,ibeg,iend,nreflag,istep,nrdof_old,nrdof_new,nstop,itag
      integer :: numRef,ians,itime,iref, numLaserIter,nLaser
      integer, dimension(5) :: iflag
      integer, allocatable, dimension(:) :: list_elem
      integer :: ic, iel, nr_elem_to_refine,iso_ans,iii,info,nonlin_steps,ref_xyz
!
!----------------------------------------------------------------------
!
!     display menu
      idec=1

!
!     initialize environment
      call begin_environment
!
!     read in HP3D input files location (if option is not present, the default value is used)
!
!                             option label      // explanation                // default value     // parameter
      call get_option_string( '-file-control'    , 'Control file'              , './files/control'  , FILE_CONTROL)
      call get_option_string( '-file-geometry'   , 'Geometry file'             , './files/cubes/cube', FILE_GEOM   )
      call get_option_string( '-file-phys'       , 'Physics file'              , './files/physics'  , FILE_PHYS   )
      call get_option_string( '-file-refinement' , 'Refinement files location' , '../../files/ref'  , FILE_REFINE )
      call get_option_string( '-file-history'    , 'History file'              , './files/history'  , FILE_HISTORY)
      call get_option_string( '-file-err'        , 'Error file'                ,'./files/dump_err'  , FILE_ERR    )
      call get_option_string( '-prefix'          , 'Prefix paraview file'      ,'laser1_'            , PREFIX      )
!
!     -- Parview Interface --
      ! Variables relevant to src/modules/paraview
!                        option label     // explanation                        // default value     // parameter
      call get_option_string('-file-vis-upscale','Visualization upscale file location','../../files/vis',FILE_VIS          )
      call get_option_string('-vis-level'       ,'Visualization upscale level (0-3)'  ,'2'                 ,VLEVEL            )
      call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'./output/figures/test_resid/' ,PARAVIEW_DIR      )
      call get_option_bool(  '-paraview-geom'   ,'Dump geom at every Paraview call'   ,.TRUE.              ,PARAVIEW_DUMP_GEOM)
      call get_option_bool(  '-paraview-attr'   ,'Dump solution to Paraview'          ,.TRUE.              ,PARAVIEW_DUMP_ATTR)
!
!
!     read in problem dependent parameters (defined in module parametersDPG)
!
!                              option label      // explanation                // default value     // parameter
      call get_option_int(    '-nord-add'           , 'NORD_ADD'                  , 2                  , NORD_ADD    )
      call get_option_int(    '-order-approx'       , 'ORDER_APPROX'              , 2                  , ORDER_APPROX)
      call get_option_int(    '-orderx'             , 'NPX'                       , 2                  , NPX         )
      call get_option_int(    '-ordery'             , 'NPY'                       , 1                  , NPY         )
      call get_option_int(    '-orderz'             , 'NPZ'                       , 0                  , NPZ         )
      call get_option_int(    '-isol'               , 'iSol'                      , 2                  , ISOL        )
      call get_option_int(    '-problem'            , 'NO_PROBLEM'                , 3                  , NO_PROBLEM  )
      call get_option_int(    '-geometry-no'        , 'Geometry file number'      , 1                  , GEOM_NO     )
      call get_option_int(    '-inner-product'      , 'INNER_PRODUCT'             , 1                  , INNER_PRODUCT)
      call get_option_real(   '-kappa'              , 'kappa'                     , 1.d0               , KAPPA       )
      call get_option_real(   '-deltat'             , 'deltaT'                    , 0.1d0              , DELTAT      )
      call get_option_real(   '-tmax'               , 'Tmax'                      , 1.0d0              , TMAX        )
      call get_option_int(    '-comp'               , 'ICOMP_EXACT'               , 2                  , ICOMP_EXACT )
      call get_option_int(    '-lasermode'          , 'LASER_MODE'                , 0                  , LASER_MODE  )
!
      call get_option_real(   '-mu'                 , 'MU'                        , 1.d0               , MU          )
      call get_option_real(   '-epsilon'            , 'EPSILON'                   , 1.d0               , EPSILON     )
      call get_option_real(   '-sigma'              , 'SIGMA'                     , 0.d0               , SIGMA       )
!
!  ...single cube problem: do not forget to reset the flag in the control file
      call get_option_real(   '-omega'              , 'OMEGA'                     , 2.5d0          , OMEGA        )
      !call get_option_real(   '-omega'              , 'OMEGA'                     , dsqrt(5.d0)*PI          , OMEGA        )
      call get_option_real(   '-gamma'              , 'GAMMA_IMP'                 , 1.0d0             , GAMMA_IMP    )
      call get_option_real(   '-waist'              , 'BEAM_WAIST'                , 0.50d0             , BEAM_WAIST   )
      call get_option_int(    '-ibc'                , 'IBCFlag'                   , 0                  , IBCFlag     )
      call get_option_int(    '-nlflag'             , 'NONLINEAR_FLAG'            , 0                  , NONLINEAR_FLAG)
      call get_option_real(   '-ntheta'             , 'NTHETA'                    , 1.d0               , NTHETA      )
      !GAMMA_IMP = dsqrt(1.d0-(PI**2/OMEGA**2))
      !GAMMA_IMP = 0.50d0
!
!aZyUta1!
!     finalize
      call end_environment

      !Number of time steps for Heat equation NSTEPS = TMAX/DELTAT
      NSTEPS = int(TMAX/DELTAT)
      NSTEPS = 1

      ! set LASER_MODE = 0
      !LASER_MODE = 0

!  ...add flags 6,7,8  to Dirichlet flag list
      call add_dirichlet_to_list(6)
      call add_dirichlet_to_list(7)
      call add_dirichlet_to_list(8)



!
!     print fancy header
      write(*,*)'                   '
      write(*,*)'                   '
      write(*,*)'                   '
      write(*,*)'// ------ THE LASER PROBLEM ------ //'
      write(*,*)'                   '
      write(*,*)'                   '
      write(*,*)'                   '
!
!     initialize problem
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
!
!     -- INTERACTIVE MODE --
!
!     display menu in infinite loop
      do while(idec /= 0)
!
        write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
        write(*,*) 'Quit ...................................0'
        write(*,*) '                                         '
        write(*,*) 'Geometry graphics (X11) ................1'
        write(*,*) 'HP3D graphics (X11) ....................3'
        write(*,*) 'Paraview ...............................4'
        write(*,*) 'Print Data Structure arrays ............5'
        write(*,*) 'DumpOut Data Structure .................7'
        write(*,*) 'DumpIn Data Structure  .................8'
        write(*,*) '                                         '
        write(*,*) ' --  Geometry & Refinements  --          '
        write(*,*) 'Geometry error ........................12'
        write(*,*) 'Interactive H-refinements .............31'
        write(*,*) 'Uniform     H-refinements .............32'
        write(*,*) 'Adaptive    H-refinements .............33'
        write(*,*) 'Anisotropic H-refinements .............34'
        write(*,*) '                                         '
        write(*,*) 'Solve (frontal) .......................40'
        write(*,*) 'Solve (MUMPS) .........................50'
        write(*,*) 'Solve (UHM) ...........................60'
        write(*,*) '                                         '
        write(*,*) 'Compute H1 error .....................100'
        write(*,*) 'Compute L2 error .....................101'
        write(*,*) 'Compute residual .....................110'
        write(*,*) 'Compute geometry error................111'
        write(*,*) 'Adaptive DPG refinements .............120'
        write(*,*) 'Compute BC data interpolation error...130'
        write(*,*) '                                         '
        write(*,*) 'Uniform refinments rate tests.........140'
        write(*,*) 'Uniform geometry error rate tests.....141'
        write(*,*) '2 uniform n anisotropic refinements...142'
        write(*,*) 'My tests..... ........................200'
        write(*,*) 'Test gain function....................201'
        write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
        read( *,*) idec
!
        select case(idec)
!  .....quit
        case(0) ; call finalize ; stop
!
!  .....GMP graphics
        case(1) ; call graphg
!
!  .....hp3d graphics
        case(3) ; call graphb

!  Paraview graphics
        case(4) ; call paraview_driver(iParAttr)
!
!  .....print data structure
        case(5) ; call result
!
!  .....dump out data structure
        case(7)
          write(*,*)'dumping out GMP...'
          call dumpout_GMP
          write(*,*)'dumping out physics...'
          call dumpout_physics
          write(*,*)'dumping out HP3D...'
          call dumpout_hp3d('./files/dumpc3Dhp')

!  .....dump in data structure
        case(8)
          write(*,*)'dumping in GMP...'
          call dumpin_GMP('./files/dumpGMP')
          write(*,*)'dumping in physics...'
          call dumpin_physics_from_file('./files/dumpPHYS')
          write(*,*)'dumping in HP3D...'
          call dumpin_hp3d('./files/dumpc3D')
!
!  .....interactive H-refinements
        case(31)
!
!         print active elements, refinement kinds
          call display_act_elem
          call display_ref_kind
!
!         select element, refinement kind
          write(*,7010)
 7010     format(' mdle,kref =')
          read(*,*) mdle,kref
!
!         refine element
          call refine(mdle,kref)
          call close_mesh
          call update_gdof
          call update_ddof
!
!  .....uniform global H-refinements
        case(32)
          call global_href
          call close_mesh
          call update_gdof
          call update_Ddof
!
!  .....adaptive H-refinements
        case(33) ; !call adaptivity_geom(0.3d0, nvoid,rvoid)

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
        enddo
        call update_gdof
        call update_Ddof
        else
        !ref_xyz=1
        !call setAnisoRef(info,ref_xyz, info)
        ref_xyz=2
        do while(1.gt.0)
          call setAnisoRef(info,ref_xyz, info)
        enddo
        call update_gdof
        call update_Ddof
        endif


!
!  .....frontal solve
        case(40)
        write(*,*) 'SELECT:'
        write(*,*) 'heat step ............. 1'
        write(*,*) 'heat stepping ......... 2'
        write(*,*) 'time Harmonic Maxwell.. 3'
        read(*,*) ians
        select case(ians)
!
!  .....single step of transient heat problem
        case(1)
!
!  .......set variables to solve for
          PHYSAm(1:4) = (/.true.,.false.,.true.,.false./)
!
!  .......set problem #
          NO_PROBLEM = 1
          call solve1(MY_NR_RHS)
!
!  .....transient heat problem, multiple time steps
        case(2)
!
!  .......set variables to solve for
          PHYSAm(1:4) = (/.true.,.false.,.true.,.false./)
!
!  .......set problem #
          NO_PROBLEM = 2
!
!  .......use the preset final time and time step
          do itime=1,NSTEPS
            call solve1(MY_NR_RHS)
          enddo

!  .....time Harmonic Maxwell
        case(3)
!  ....... set NRFL=0
          NRFL=0
!  .......set variables to solve for
          PHYSAm(1:4) = (/.false.,.true.,.false.,.true./)
!
!  .......set problem #
          NO_PROBLEM = 3
          call solve1(MY_NR_RHS)
!
        end select
!
!  .....MUMPS solve
        case(50)
!
        write(*,*) 'SELECT:'
        write(*,*) 'heat step ............. 1'
        write(*,*) 'heat stepping ......... 2'
        write(*,*) 'time Harmonic Maxwell.. 3'
        write(*,*) 'full Laser problem..... 4'
        read(*,*) ians
        select case(ians)
!
!  .....single step of transient heat problem
        case(1)
!  .......set variables to solve for
          PHYSAm(1:4) = (/.true.,.false.,.true.,.false./)
!  .......set problem #
          NO_PROBLEM = 1
          call uhm_time_in
          call mumps_interf(MY_NR_RHS)
          call uhm_time_out(t)
          write(*,*) 'time mumps = ', t
!
!  .....transient heat problem, multiple time steps
        case(2)
!  .......set variables to solve for
          PHYSAm(1:4) = (/.true.,.false.,.true.,.false./)
!  .......set problem #
          NO_PROBLEM = 2
          !call update_gdof
          !call update_ddof
!  .......use the preset final time and time step
          do itime=1,NSTEPS
          call mumps_interf(MY_NR_RHS)
          enddo
!
!  .....time Harmonic Maxwell
        case(3)
!  ....... set NRFL=0
          NRFL=0
!  .... check if solving linear or nonlinear problem
          if(NONLINEAR_FLAG.eq.0) then
!  .......set variables to solve for
            PHYSAm(1:4) = (/.false.,.true.,.false.,.true./)
!  .......set problem #
            NO_PROBLEM = 3
            !call update_gdof
            !call update_ddof
            call uhm_time_in
!  ....... allocate copy matrices
            allocate(ZFL_EE(MYE,MYE,MAXNRFL), &
            ZFL_EQ(MYE,MYQ,MAXNRFL),&
            ZFL_QQ(MYQ,MYQ,MAXNRFL))
            call mumps_sc('H')
            !call mumps_interf(MY_NR_RHS)
            deallocate(ZFL_EE,ZFL_EQ,ZFL_QQ)
            call uhm_time_out(t)
            write(*,*) 'time mumps = ', t
          else if (NONLINEAR_FLAG.eq.1) then
!  ......  check that NEXACT = 0
            if(NEXACT.ne.0) then
              write(*,*) 'NEXACT must be 0 for nonlinear Maxwell'
              stop
            endif
!  .......set variables to solve for
            PHYSAm(1:4) = (/.false.,.true.,.false.,.true./)
!  .......set problem #
            NO_PROBLEM = 3
            !call update_gdof
            !call update_ddof
            nLaser = 0
            nonlin_steps = 2
            do nLaser =1,nonlin_steps
              call uhm_time_in
              call mumps_sc('H')
              !call mumps_interf(MY_NR_RHS)
              call uhm_time_out(t)
              write(*,*) 'time mumps = ', t
            call update_gdof
            call update_ddof
            enddo
            call update_gdof
            call update_ddof
          endif

!  .....full Laser problem
        case(4)
        if(NEXACT.ne.0) then
            write(*,*) 'Must run Laser problem with NEXACT = 0 only'
            stop
        endif
        LASER_MODE = 1
        write(*,*) 'full Laser problem: set number of coupled iteration'
        read(*,*) numLaserIter
        nLaser = 0
          call uhm_time_in
          do nLaser=1,numLaserIter
  !..................... First solve maxwell .....................
  !  .......set variables to solve for:: MAXWELL
            PHYSAm(1:4) = (/.false.,.true.,.false.,.true./)
  !
  !  .......set problem #
            NO_PROBLEM = 3
            call update_gdof
            call update_ddof
            !call mumps_interf(MY_NR_RHS)
            call mumps_sc('H')
            !call mumps_solve_seq(MY_NR_RHS)


  !..................... Next solve full heat equation..............
  !  .......set variables to solve for :: HEAT
            PHYSAm(1:4) = (/.true.,.false.,.true.,.false./)
  !
  !  .......set problem #
            NO_PROBLEM = 2
  !
            call update_gdof
            call update_ddof
  !  .......use the preset final time and time step
            do itime=1,NSTEPS
              !call mumps_interf(MY_NR_RHS)
              call mumps_sc('H')
              !call mumps_solve_seq(MY_NR_RHS)
            enddo
          enddo
          call uhm_time_out(t)
          write(*,*) 'time mumps = ', t
      end select
!
!
!  .....UHM solve
        case(60)
          call uhm_solve
          call uhm_solver_flush(UHM_SOLVER_PTR)
!
!  .....compute H1 error
        case(100)
          iflag(1) = 1; iflag(2:4) = 0; itag=1
          call compute_error(iflag,itag)

!  .....compute L2 error
        case(101)
          iflag(4) = 1; iflag(1:3) = 0; itag=1
          call compute_error(iflag,itag)

!  .....compute residual
        case(110) ; !call compute_residual

!  .....compute Geometry error
        case(111)
          geo_rnorm = 0.d0
          geo_err = 0.d0
          call geometry_error(geo_err,geo_rnorm)
!
!  .....adaptive DPG refinements
        case(120)
  333     write(*,7011)
 7011     format('main: SET INITIAL, FINAL STEP,', &
                 ' REFINEMENT FLAG (1-href,2-pref,3-hpref), FACTOR')
          read(*,*) ibeg,iend,nreflag,factor
          if (nreflag.ne.1.and.nreflag.ne.2.and.nreflag.ne.3) go to 333
          istep=ibeg-1
          do while(istep.lt.iend)
            istep=istep+1
            nrdof_old = NRDOFSH
            !!!!!!!!call adapt_DPG(istep,nreflag,factor,nstop)
            if (nstop.eq.1) exit
            nrdof_new = NRDOFSH
            if (nrdof_new.eq.nrdof_old) then
              istep=istep-1
              factor = factor*0.25d0
              write(*,7023) factor
 7023         format('main: NEW factor = ',e12.5)
              if (factor.lt.0.000001d0) exit
            endif
          enddo
!
!  .....compute BC data interpolation error
        case(130) ; call compute_BC_interp_error

! ......... rate test for single step of Heat problem
        case(140)
		write(*,*) 'Testing h-convergence rates'
        write(*,*) 'SELECT:'
        write(*,*) 'heat step ..................... 1'
        write(*,*) 'time harmonic Maxwell ......... 2'
        read(*,*) ians
!  ........ Get #refinements to be done
        write(*,*) 'Enter number of uniform H-refinements:'
        read(*,*) numRef
!
        select case(ians)
!  .....single step of transient heat problem
        case(1)
        !  .......set variables to solve for
          PHYSAm(1:4) = (/.true.,.false.,.true.,.false./)
!  .......set problem #
          NO_PROBLEM = 1
          iref = 0
          do while(iref.lt.numRef)
!  ........ Solve the problem
!            call solve1(MY_NR_RHS)
            call mumps_sc('H')
            !call mumps_solve_seq(MY_NR_RHS)
!  ........ Set error flags
            iflag(1) = 1; iflag(2:4) = 0; itag=1
!  ........ Compute error
            call compute_error(iflag,itag)
!  ......... Compute residual
            !call compute_residual
!  ........ Do uniform h-refinements
            call global_href
            call close_mesh
            call update_gdof
            call update_Ddof
            iref = iref+1
          enddo

!  .....time harmonic Maxwell
        case(2)
!  ....... set NRFL=0
          NRFL=0
          if(NONLINEAR_FLAG.ne.0) then
            write(*,*) 'can test rates only for linear Maxwell'
            stop
          endif
!  .......set variables to solve for
          PHYSAm(1:4) = (/.false.,.true.,.false.,.true./)
!  .......set problem #
          NO_PROBLEM = 3
          iref = 0
          do while(iref.lt.numRef)
!  ........ Solve the problem
!            call solve1(MY_NR_RHS)
             !call mumps_interf(MY_NR_RHS)
            call uhm_time_in
            call mumps_sc('H')
            !call mumps_solve_seq(MY_NR_RHS)
!  ........ Set error flags
            iflag(4) = 1; iflag(1:3) = 0; itag=1
!  ........ Compute error
            call compute_error(iflag,itag)
            call uhm_time_out(t)
!  ......... Compute residual
            !call compute_residual
!  ........ Do uniform h-refinements
            call global_href
            call close_mesh
            call update_gdof
            call update_Ddof
            iref = iref+1
          enddo
        end select

! ......... rate test geometry error
        case(141)
    write(*,*) 'Testing h-convergence rates for geometry error'
!  ........ Get #refinements to be done
        write(*,*) 'Enter number of uniform H-refinements:'
        read(*,*) numRef
          iref = 0
          do while(iref.lt.numRef)
!  ........ Compute geometry error
          geo_rnorm = 0.d0
          geo_err = 0.d0
          call geometry_error(geo_err,geo_rnorm)
!  ........ Do uniform h-refinements
            call global_href
            call close_mesh
            call update_gdof
            call update_Ddof
            iref = iref+1
          enddo
!  .....tests
        case(200) ; call my_tests
!  .....tests
        case(201)
        xtest = 0.99d0
        do i=1,100
          Ehat = 1.d0/real(i,8)
          call get_gainFunction(Ehat, gainF)
        enddo
! ......... 1 uniform n anisotropic refinements and solve
        case(142)
        write(*,*) '2 uniform n anisotropic refinements and solve'
!  ........ Get #refinements to be done
        write(*,*) 'Enter number of anisotropic H-refinements:'
        read(*,*) numRef
!  .......set variables to solve for
            PHYSAm(1:4) = (/.false.,.true.,.false.,.true./)
!  .......set problem #
            NO_PROBLEM = 3
        info = 0
!  ........ Set error flags
        iflag(4) = 1; iflag(1:3) = 0; itag=1
        call global_href
        !call close_mesh
        call update_gdof
        call update_Ddof
        call mumps_sc('H')
        write(*,*) 'after one uniform refinement'
        call compute_error(iflag,itag)

        call global_href
        !call close_mesh
        call update_gdof
        call update_Ddof
        call mumps_sc('H')
        write(*,*) 'after 2 uniform refinement'
        call compute_error(iflag,itag)
          iref = 0
          do while(iref.lt.numRef)
!  ........ Do anisotropic h-refinements
            ref_xyz=2
            call setAnisoRef(info,ref_xyz, info)
            call update_gdof
            call update_Ddof
! ........... solve
            call mumps_sc('H')
!  ........ Compute error
            write(*,*) 'after anisotropic refinement'
            call compute_error(iflag,itag)
            iref = iref+1
          enddo


!
        end select
!
!  ...end of infinite loop
      enddo
!
!     finalize library
      call finalize
!
!
endprogram main



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
