
c  ...save copies of enriched stiffness matrices
      STIFFHHc=STIFFHH ; STIFFHVc=STIFFHV
cc
cc  ...G^-1 * Load
      call ZPPTRS(uplo, nrdofHH, 1, AP_Heat, BLOADH, MAXbrickHH, info1 )
      if (info1.ne.0) then
        write(*,*) 'elem_dpgHeat: info1 = ',info1 ; stop
      endif
cc
cc  ...G^-1 * StiffHH
      call ZPPTRS(uplo,nrdofHH,nrdofH,AP_Heat,STIFFHHc,MAXbrickHH,info2)
      if (info2.ne.0) then
        write(*,*) 'elem_dpgHeat: info2 = ',info2 ; stop
      endif
c
cc  ...G^-1 * StiffHV
      call ZPPTRS(uplo,nrdofHH,nrdofV,AP_Heat,STIFFHVc,MAXbrickHH,info3)
      if (info3.ne.0) then
        write(*,*) 'elem_dpgHeat: info3 = ',info3 ; stop
      endif
cc
cc
cc  ...ultimate DPG load vector and stiffness matrices
      do k1=1,nrdofH
c
c       (G^-1 * LoadH)   * StiffHH
        do k=1,nrdofHH
          ZblocH(k1) = ZblocH(k1) + BLOADH(k)*STIFFHH(k,k1)
        enddo
c
c       (G^-1 * StiffHH) * StiffHH
        do k2=1,nrdofH
          do k=1,nrdofHH
            ZalocHH(k1,k2) = ZalocHH(k1,k2)
     .                     + STIFFHHc(k,k1)* STIFFHH(k,k2)
         enddo
        enddo
c
c       (G^-1 * StiffHH) * StiffHV
        do k2=1,nrdofV
          do k=1,nrdofHH
            ZalocHV(k1,k2) = ZalocHV(k1,k2)
     .                     + STIFFHHc(k,k1)* STIFFHV(k,k2)
          enddo
        enddo
      enddo
c
      do k1=1,nrdofV
c
c       (G^-1 * LoadH)   * StiffHV
        do k=1,nrdofHH
          ZblocV(k1) = ZblocV(k1) + BLOADH(k)*STIFFHV(k,k1)
        enddo
c
c       (G^-1 * StiffHV) * StiffHH
        do k2=1,nrdofH
          do k=1,nrdofHH
            ZalocVH(k1,k2) = ZalocVH(k1,k2)
     .                     + STIFFHVc(k,k1)* STIFFHH(k,k2)
          enddo
        enddo

c       (G^-1 * StiffHV) * StiffHV
        do k2=1,nrdofV
          do k=1,nrdofHH
            ZalocVV(k1,k2) = ZalocVV(k1,k2)
     .                     + STIFFHVc(k,k1)* STIFFHV(k,k2)
          enddo
        enddo
      enddo
cc
cc  ...check symmetry
c      diffmax = ZERO; dmax = ZERO
c      do k1=1,nrdofH
c        do k2=k1,nrdofH
c          diffmax = max(diffmax,abs(ZalocHH(k1,k2)-ZalocHH(k2,k1)))
c          dmax = max(dmax,abs(ZalocHH(k1,k2)))
c        enddo
c      enddo
c      symmetry_tol = 1.d-9
c      if (diffmax/dmax.gt.symmetry_tol) then
c        write(*,7021) diffmax, dmax
c 7021   format('elem_dpgHeat: diffmax,dmax FOR ZalocHH = ',2e12.5)
c        call pause
c      endif
c      diffmax = ZERO; dmax = ZERO
c      do k1=1,nrdofV
c        do k2=k1,nrdofV
c          diffmax = max(diffmax,abs(ZalocVV(k1,k2)-ZalocVV(k2,k1)))
c          dmax = max(dmax,abs(ZalocVV(k1,k2)))
c        enddo
c      enddo
c      if (diffmax/dmax.gt.symmetry_tol) then
c        write(*,7022) diffmax, dmax
c 7022   format('elem_dpgHeat: diffmax,dmax FOR ZalocVV = ',2e12.5)
c        call pause
c      endif
c      diffmax = ZERO; dmax = ZERO
c      do k1=1,nrdofH
c        do k2=1,nrdofv
c          diffmax = max(diffmax,abs(ZalocHV(k1,k2)-ZalocVH(k2,k1)))
c          dmax = max(dmax,abs(ZalocHV(k1,k2)))
c        enddo
c      enddo
c      if (diffmax/dmax.gt.symmetry_tol) then
c        write(*,7023) diffmax, dmax
c 7023   format('elem_dpgHeat: diffmax,dmax FOR ZalocHV = ',2e12.5)
c        call pause
c      endif
cc
cc
c      if (iprint.ge.1) then
c        write(*,7010)
c 7010   format('elem_dpgHeat: ZblocH,ZblocV = ')
c        write(*,7011) ZblocH(1:NrdofH)
c        write(*,7011) ZblocV(1:NrdofV)
#if C_MODE
 7011   format(6(2e10.3,2x))
#else
 7011   format(10e12.5)
#endif
c        write(*,7012)
c 7012   format('elem_dpgHeat: ZalocHH = ')
c        do i=1,NrdofH
c          write(*,7013) i,ZalocHH(i,1:NrdofH)
c#if C_MODE
c 7013     format('i = ',i3,10(/,6(2e10.3,2x)))
c#else
c 7013     format('i = ',i3,10(/,10e12.5))
c#endif
c        enddo
c        write(*,7014)
c 7014   format('elem_dpgHeat: ZalocHV = ')
c        do i=1,NrdofH
c          write(*,7013) i,ZalocHV(i,1:NrdofV)
c        enddo
c        write(*,7015)
c 7015   format('elem_dpgHeat: ZalocVV = ')
c        do i=1,NrdofV
c          write(*,7013) i,ZalocVV(i,1:NrdofV)
c        enddo
c      endif



























c  ...save copies of enriched stiffness matrices
      STIFFEEc =STIFFEE; STIFFEQc =  STIFFEQ
c
c  ...compute the products of inverted test matrix with RHS
c     and enriched stiffness matrices
      call ZPPTRS(uplo, nrtest, 1, AP_Maxwell, BLOADE, MAXtestE, info1 )
      if (info1.ne.0) then
        write(*,*) 'elem: info1 = ',info1
        stop1
      endif
      call ZPPTRS(uplo,nrtest,2*nrdofE,AP_Maxwell,STIFFEEc,MAXtestE,
     .           info2)
      if (info2.ne.0) then
        write(*,*) 'elem: info2 = ',info2
        stop1
      endif
      call ZPPTRS(uplo,nrtest,6*nrdofQ,AP_Maxwell,STIFFEQc,MAXtestE,
     .           info3)
      if (info3.ne.0) then
        write(*,*) 'elem: info3 = ',info3
        stop1
      endif
c
c  ...compute the ultimate DPG load vectors and stiffness matrices
      do k1=1,2*nrdofE
        do k=1,nrtest
          ZblocE(k1) = ZblocE(k1) + BLOADE(k)*conjg(STIFFEE(k,k1))
        enddo
        do k2=1,2*nrdofE
          do k=1,nrtest
            ZalocEE(k1,k2) = ZalocEE(k1,k2)
     .                     + STIFFEEc(k,k2)*conjg(STIFFEE(k,k1))
          enddo
        enddo
        do k2=1,6*nrdofQ
          do k=1,nrtest
            ZalocEQ(k1,k2) = ZalocEQ(k1,k2)
     .                     + STIFFEQc(k,k2)*conjg(STIFFEE(k,k1))
          enddo
        enddo
      enddo
      do k1=1,6*nrdofQ
        do k=1,nrtest
          ZblocQ(k1) = ZblocQ(k1) + BLOADE(k)*conjg(STIFFEQ(k,k1))
        enddo
        do k2=1,2*nrdofE
          do k=1,nrtest
            ZalocQE(k1,k2) = ZalocQE(k1,k2)
     .                     + STIFFEEc(k,k2)*conjg(STIFFEQ(k,k1))
          enddo
        enddo
        do k2=1,6*nrdofQ
          do k=1,nrtest
            ZalocQQ(k1,k2) = ZalocQQ(k1,k2)
     .                     + STIFFEQc(k,k2)*conjg(STIFFEQ(k,k1))
          enddo
        enddo
      enddo
c
c  ...check symmetry
      diffmax = ZERO; dmax = ZERO
      do k1=1,2*nrdofE
        do k2=k1,2*nrdofE
          diffmax = max(diffmax,
     .              abs(ZalocEE(k1,k2)-conjg(ZalocEE(k2,k1))))
          dmax = max(dmax,abs(ZalocEE(k1,k2)))
        enddo
      enddo
      symmetry_tol = 1.d-9
      if (diffmax/dmax.gt.symmetry_tol) then
        write(*,7021) diffmax, dmax
 7021   format('elem_dpgMaxwell: diffmax,dmax FOR ZalocEE = ',2e12.5)
        call pause
      endif
      diffmax = ZERO; dmax = ZERO
      do k1=1,6*nrdofQ
        do k2=k1,6*nrdofQ
          diffmax = max(diffmax,
     .                  abs(ZalocQQ(k1,k2)-conjg(ZalocQQ(k2,k1))))
          dmax = max(dmax,abs(ZalocQQ(k1,k2)))
        enddo
      enddo
      if (diffmax/dmax.gt.symmetry_tol) then
        write(*,7022) diffmax, dmax
 7022   format('elem_dpgMaxwell: diffmax,dmax FOR ZalocQQ = ',2e12.5)
        call pause
      endif
      diffmax = ZERO; dmax = ZERO
      do k1=1,2*nrdofE
        do k2=1,6*nrdofQ
          diffmax = max(diffmax,
     .                  abs(ZalocEQ(k1,k2)-conjg(ZalocQE(k2,k1))))
          dmax = max(dmax,abs(ZalocEQ(k1,k2)))
        enddo
      enddo
      if (diffmax/dmax.gt.symmetry_tol) then
        write(*,7023) diffmax, dmax
 7023   format('elem_dpgMaxwell: diffmax,dmax FOR ZalocEQ = ',2e12.5)
        call pause
      endif

      if (iprint.eq.10) then
        zE = ZERO; zQ = ZERO
        zE(1) = 1.d0; zE(5) = 1.d0; zE(9) = 1.d0; zE(13) = 1.d0
        zQ(1) = 1.d0
        do k1=1,2*nrdofE
          zres = ZblocE(k1)
          do k2=1,2*nrdofE
            zres = zres - ZalocEE(k1,k2)*zE(k2)
          enddo
          do k2=1,6*nrdofQ
            zres = zres - ZalocEQ(k1,k2)*zQ(k2)
          enddo
          write(*,7300) k1,zres
 7300     format('RESIDUAL FOR k1 = ',i3,2x,2e12.5)
        enddo
        do k1=1,6*nrdofQ
          zres = ZblocQ(k1)
          do k2=1,2*nrdofE
            zres = zres - ZalocQE(k1,k2)*zE(k2)
          enddo
          do k2=1,6*nrdofQ
            zres = zres - ZalocQQ(k1,k2)*zQ(k2)
          enddo
          write(*,7300) k1,zres
        enddo
        call pause
      endif
c

c
c
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
c
      end subroutine elem_dpgMaxwell











