c--------------------------------------------------------------------
c
c     routine name      - adapt_DPG
c
c--------------------------------------------------------------------
c
c     latest revision:  - May 14
c
c     purpose:          - routine performs a single h-adaptivity
c                         step using the DPG energy (residual)
c                         error indicators and the greedy strategy
c
c     arguments:
c
c     in:
c             Istep     - step number
c             Nreflag   = 1 for h-refinements
c                       = 2 for p-refinements
c                       = 3 for h- or p-refinement, depending
c                           upon the element size
c             Factor    - if element error \ge Factor*max_error
c                         then the element is refined
c     out:
c             Nstop     = 1 if no refinement has been made
c
c---------------------------------------------------------------------
c
      subroutine adapt_DPG(Istep,Nreflag,Factor, Nstop)
c
      use control
      use data_structure3D
      use uhm2
      use problem
#include "syscom.blk"
c
      parameter (max_step = 20)
      dimension nrdof_mesh(max_step),residual_mesh(max_step),
     .          rate_mesh(max_step),
     .          error_H1_mesh(max_step),rel_error_H1_mesh(max_step),
     .          ratio_H1vsEN_mesh(max_step)
c
      save nrdof_mesh,residual_mesh,rate_mesh,error_H1_mesh,
     .     rel_error_H1_mesh,ratio_L2vsEN_mesh
c
      double precision, allocatable, dimension(:) :: elem_resid
      double precision, allocatable, dimension(:) :: elem_H1_error
      double precision, allocatable, dimension(:) :: elem_H1_norm
      integer, allocatable, dimension(:) :: elem_ref_flag
      integer, allocatable, dimension(:) :: list_elem
      integer, allocatable, dimension(:) :: list_ref_flags
c
c  ...a flag indicating which variables should be integrated for the error
      integer, dimension(NR_PHYSA) :: iflag
c
c-----------------------------------------------------------------------
c
      reps = 1.d-10
      iprint=0
      if (iprint.eq.1) write(*,*) 'adapt_DPG: DEBUGGING...'
c
c  ...set up the flag for error computation
      iflag(1)=1; iflag(2)=0
c
c  ...solve the problem on the current mesh
ccc      call solve1(MY_NR_RHS)
ccc      call mumps_solve_seq(MY_NR_RHS)
      call uhm_solve
      call uhm_solver_flush(UHM_SOLVER_PTR)
      if (iprint.eq.1) write(*,*) 'adapt_DPG: HAVE SOLVED...'
c
      allocate (elem_resid(NRELES),elem_ref_flag(NRELES),
     .          list_elem(NRELES),list_ref_flags(NRELES))
      if (NEXACT.ge.1) then
        allocate (elem_H1_error(NRELES))
        allocate (elem_H1_norm(NRELES))
        elem_H1_error = 0.d0; elem_H1_norm = 0.d0
      endif
c
      ic=0
      residual = 0.d0; elem_resid_max = 0.d0
      errorH1 = 0.d0; rnormH1 = 0.d0; nrdof_total = NRDOFSH+NRDOFSV
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_residual(mdle, elem_resid(iel),elem_ref_flag(iel))
        if (NEXACT.ge.1) then
          call element_error(mdle,iflag, errorH,errorE,errorV,errorQ,
     .                                   rnormH,rnormE,rnormV,rnormQ)
          if (iprint.eq.1) then
            write(*,7005) iel,errorH,rnormH
 7005       format('adapt_DPG: iel, errorH,rnormH = ',i4,2e12.5)
          endif
          elem_H1_error(iel) = errorH
          elem_H1_norm(iel)  = rnormH
          errorH1 = errorH1 + errorH
          rnormH1 = rnormH1 + rnormH
        endif
        if (iprint.eq.1) then
          write(*,7010) iel,mdle,elem_resid(iel)
 7010     format('adapt_DPG: iel,mdle = ',i5,i6,' RESIDUAL = ',e12.5)
        endif
        elem_resid_max = max(elem_resid_max,elem_resid(iel))
        residual = residual + elem_resid(iel)
        call find_ndof(mdle, ndofH,ndofE,ndofV,ndofQ)
        nrdof_total = nrdof_total - ndofV
       enddo
c
c-----------------------------------------------------------------------
c
      nrdof_mesh(Istep) = nrdof_total
      residual_mesh(Istep) =  sqrt(residual)
      if (NEXACT.ge.1) then
         error_H1_mesh(Istep) = sqrt(errorH1)
         rel_error_H1_mesh(Istep) = sqrt(errorH1/rnormH1)
         ratio_H1vsEN_mesh(Istep)
     .    = residual_mesh(Istep)/error_H1_mesh(Istep)
      endif
c
c  ...compute decrease rate for the residual
      select case(Istep)
      case(1)
        rate_mesh(Istep) = 0.d0
      case default
        rate_mesh(Istep) =
     .         log(residual_mesh(Istep-1)/residual_mesh(Istep))/
     .         log(float(nrdof_mesh(Istep-1))/float(nrdof_mesh(Istep)))
      end select
c
c  ...check if any dof have been added
      if (Istep.gt.1) then
        if (nrdof_mesh(Istep).eq.nrdof_mesh(Istep-1)) then
          Nstop=1; return
        else
          Nstop=0
        endif
      else
        Nstop=0
      endif
c
c  ...print out the history of refinements
      write(*,7002)
 7002 format('adapt: HISTORY OF REFINEMENTS')
      do i=1,Istep
        if (NEXACT.eq.0) then
          write(*,7003) nrdof_mesh(i),residual_mesh(i),rate_mesh(i),i
 7003     format('nrdof = ',i6,' residual = ',e12.5,' rate = ',f7.2,
     .           ' iteration = ',i2)
        else
          write(*,7004) nrdof_mesh(i),residual_mesh(i),rate_mesh(i),
     .    error_H1_mesh(i),rel_error_H1_mesh(i),ratio_H1vsEN_mesh(i),i
 7004     format('nrdof = ',i6,' residual = ',e12.5,' rate = ',f7.2,
     .           ' H1 error = ',e12.5,' relative H1 error = ',e12.5,
     .           ' ratio = ',f7.2,' iteration = ',i2)
        endif
      enddo
      call pause
c
c-----------------------------------------------------------------------
c
c  ...use the greedy strategy to determine which elements to refine
      ic=0
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        if (elem_resid(iel).ge.Factor*elem_resid_max) then
          ic=ic+1
          list_elem(ic) = mdle
          list_ref_flags(ic) = elem_ref_flag(iel)
          if (iprint.eq.1) then
            write(*,7037) ic,mdle
 7037       format('adapt_DPG: ADDING TO THE LIST: ic,mdle = ',2i6)
          endif
        endif
      enddo
      nr_elem_to_refine = ic
c
c  ...refine the elements from the list
      do iel=1,nr_elem_to_refine
        mdle = list_elem(iel)
        kref = list_ref_flags(iel)
c
c  .....restrictions on h-refinements
        idec_href=1
ccc        call find_gen(mdle,ngen)
ccc        select case(NODES(mdle)%type)
ccc        case('mdlt')
ccc          if (ngen(1).eq.MAXGENT) idec_href=0
ccc        case('mdlq')
ccc         if (ngen(1).eq.MAXGENQ(1)) idec_href=0
ccc         if (ngen(2).eq.MAXGENQ(2)) idec_href=0
ccc        end select
c
c  .....p-refinements limited by the maximum order supported
c       by the code
        nord = NODES(mdle)%order
ccc        select case(NODES(mdle)%type)
ccc        case('mdlt')
ccc          nord_new = min(nord+1,MAXP-NORD_ADD-1)
ccc        case('mdlq')
ccc          call decode(nord, nordh,nordv)
ccc         nordh_new = nordh; nordv_new = nordv
ccc          if (kref(1).eq.1) then
ccc            nordh_new = min(nordh+1,MAXP-NORD_ADD-1)
ccc          endif
ccc          if (kref(2).eq.1) then
ccc            nordv_new = min(nordv+1,MAXP-NORD_ADD-1)
ccc          endif
ccc          nord_new = nordh_new*10+nordv_new
ccc        end select
ccc        if (nord_new.eq.nord) then
ccc          idec_pref=0
ccc        else
ccc          idec_pref=1
ccc        endif
c
c
        select case(Nreflag)
        case(1)
          if (idec_href.eq.1) then
ccc            call get_isoref(mdle, kref)
            call refine(mdle,kref)
            write(*,*) 'mdle,kref = ',mdle,kref
          endif
        case(2)
          write(*,*) 'p UNIFINISHED '; stop1
ccc          if (idec_pref.eq.1) then
ccc            call enrich(mdle, nord_new)
ccc            if (iprint.eq.3) then
ccc              write(*,6001) mdle,nord_new
ccc 6001         format('adapt_DPG: HAVE ENRICHED mdle = ',i6,' nordnew = ',i2)
ccc            endif
ccc          endif
        case(3)
c          write(*,*) 'hp UNIFINISHED '; stop1

c  .......call upon the user to decide whether the element should be
c         h- or p-refined
ccc          call select_refinement(Mdle, kref,iref)
ccc          select case(iref)
c
c  .......h-refinement selected
ccc          case(1)
ccc            if (idec_href.eq.1) then
ccc            if (iprint.eq.2) then
ccc                write(*,7011) iel,mdle,kref
 7011           format('adapt_DPG: REFINING     iel = ',i5,' mdle = ',i5,
     .                 ' kref = ',2i1)
ccc              endif
ccc              call refine(mdle,kref(1)*10+kref(2))
ccc             if (iprint.eq.2) then
ccc                write(*,7013) iel,mdle,kref
 7013           format('adapt_DPG: HAVE REFINED iel = ',i5,' mdle = ',i5,
     .                 ' kref = ',2i1)
ccc              endif
ccc            endif
c
c  .......p-refinement selected
ccc          case(2)
ccc            if (idec_pref.eq.1) then
ccc              if (iprint.eq.3) then
ccc                write(*,7012) iel,mdle,nord_new
 7012           format('adapt_DPG: ENRICHING iel = ',i5,' mdle = ',i5,
     .                 ' nord_new = ',i2)
ccc              endif
ccc              call enrich(mdle, nord_new)
ccc            endif
ccc         end select
        end select
      enddo
c
c  ...close the mesh
      call close
c
c  ...update geometry and Dirichlet flux dof after the refinements
      call update_gdof
      call update_Ddof
c
      if (iprint.ge.1) write(*,*) 'adapt_DPG: DEALLOCATING...'
      deallocate (elem_resid,elem_ref_flag,list_elem,list_ref_flags)
      if (NEXACT.ge.1) then
        deallocate (elem_H1_error,elem_H1_norm)
      endif


      if (iprint.ge.1) write(*,*) 'adapt_DPG: EXITING...'
c
c
      end subroutine adapt_DPG
c
