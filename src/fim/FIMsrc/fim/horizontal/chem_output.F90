module module_chem_output

contains

  subroutine chem_output(its,nts,aod2d,exch,p10,pm25,pr3d,tk3d,tr,trfall,&
    phys2dwrf,tr1_tavg)

    use module_constants,only:rd
    use module_control,only:ArchvStep,ArchvTimeUnit,dt,filename_len,nip,ntra,&
      ntrb,nvl,nvlp1
!sms$ignore begin
    use icosio,only:icosio_out
    use module_initial_chem_namelists
!sms$ignore end
    use module_wrf_control,only:num_chem,num_moist,nvl_gocart,nvl_gocart
    use module_header,only:header

    integer,intent(in)::its,nts
!sms$distribute (dh,nip) begin
    real,intent(inout):: tr1_tavg(nvl,nip)
    real,intent(in)::aod2d(nip),exch(nvl,nip),&
      p10(nvl,nip),pm25(nvl,nip),pr3d(nvlp1,nip),&
      tk3d(nvl,nip),tr(nvl,nip,ntra+ntrb),trfall(nip,num_chem)
    real::d1st(nvl,nip),d2st(nvl,nip),d3st(nvl,nip),d4st(nvl,nip),d5st(nvl,nip),&
      dms1(nvl,nip),intaer(nip),intash(nip),intbc(nip),intdust(nip),intoc(nip),&
      intsulf(nip),qcct(nvl,nip),qict(nvl,nip),qrct(nvl,nip),qsct(nvl,nip),&
      rho_phys(nvl,nip),sea1(nvl,nip),sea2(nvl,nip),sea3(nvl,nip),sea4(nvl,nip),&
      trco(nvl,nip)
!sms$distribute end
!sms$distribute(dh,1) begin
    real,intent(in)::phys2dwrf(:,:) ! (nip,:)
!sms$distribute end
    integer::ichem_start,imoist_start,j,k
    real::dpsum
    character(len=filename_len)::filename
    integer::its2time
    integer::time

    if (mod(its,ArchvStep)==0.or.(its==nts.and.ArchvTimeUnit.eq.'ts')) then

      time=its2time(its)

      ichem_start=ntra
!     imoist_start=4-1

!SMS$IGNORE BEGIN
!TBH:  Added this IGNORE to work around a PPP core dump between here and 
!TBH:  "SMS$IGNORE END".  Mark Govett is investigating...  

      if ((.not.mp_physics==0).or.(.not.cu_physics==0)) then
        trco(:,:) = tr(:,:,5)
        trco(1,:) = phys2dwrf(:,1)
        trco(2,:) = phys2dwrf(:,2)
        trco(3,:) = phys2dwrf(:,3)
        trco(4,:) = phys2dwrf(:,4)
        trco(5,:) = phys2dwrf(:,5)
        trco(6,:) = phys2dwrf(:,6)
      endif
      if (chem_opt == 300) then
        d3st(:,:) = tr(:,:,ichem_start+p_dust_3)
        d4st(:,:) = tr(:,:,ichem_start+p_dust_4)
        d5st(:,:) = tr(:,:,ichem_start+p_dust_5)
        sea3(:,:) = tr(:,:,ichem_start+p_seas_3)
        sea4(:,:) = tr(:,:,ichem_start+p_seas_4)
      endif
      if (chem_opt == 500) then
        d1st(:,:) = tr(:,:,ichem_start+p_tr1)
        d2st(:,:) = tr(:,:,ichem_start+p_tr2)
      endif
      if (chem_opt >= 300 .and. chem_opt < 500) then
        dms1(:,:) = tr(:,:,ichem_start+p_dms)
        d1st(:,:) = tr(:,:,ichem_start+p_dust_1)
        d2st(:,:) = tr(:,:,ichem_start+p_dust_2)
        sea1(:,:) = tr(:,:,ichem_start+p_seas_1)
        sea2(:,:) = tr(:,:,ichem_start+p_seas_2)
      endif
!     if (mp_physics == 2) then
!       qcct(:,:) = tr(:,:,imoist_start+p_qc)
!       qrct(:,:) = tr(:,:,imoist_start+p_qr)
!       qict(:,:) = tr(:,:,imoist_start+p_qi)
!       qsct(:,:) = tr(:,:,imoist_start+p_qs)
!     endif
!SMS$IGNORE END
      if (chem_opt >= 300 .and. chem_opt < 500) then
!SMS$PARALLEL(dh, j) BEGIN
        do j=1,nip
          dpsum=0.
          intash(j)=0.
          intaer(j)=0.
          intbc(j)=0.
          intoc(j)=0.
          intsulf(j)=0.
          intdust(j)=0.
          do k=1,nvl
            dpsum=dpsum+(pr3d(k,j)-pr3d(k+1,j))
            rho_phys(k,j)=.5*(pr3d(k,j)+pr3d(k+1,j))&
              /(RD*tk3d(k,j)) !*(1.+.608*qv3d(k,j))
            intaer(j)=intaer(j)+pm25(k,j)*(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
            intbc(j)=intbc(j)+(tr(k,j,ichem_start+p_bc1)&
              +tr(k,j,ichem_start+p_bc2))*(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
            intoc(j)=intoc(j)+(tr(k,j,ichem_start+p_oc1)&
              +tr(k,j,ichem_start+p_oc2))*(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
            intsulf(j)=intsulf(j)+tr(k,j,ichem_start+p_sulf)*(pr3d(k,j)&
              -pr3d(k+1,j))*rho_phys(k,j)
            if (chem_opt == 300)intdust(j)=intdust(j)+(d1st(k,j)&
              +.286*d2st(k,j))*(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
            if (chem_opt==304.or.chem_opt==316.or.chem_opt==317) then
              intdust(j)=intdust(j)+d1st(k,j)*(pr3d(k,j)-pr3d(k+1,j))&
                *rho_phys(k,j)
            endif
            if (chem_opt == 316) then
              intash(j)=intash(j)+(tr(k,j,ichem_start+p_vash_1)       &
                + tr(k,j,ichem_start+p_vash_2)       &
                + tr(k,j,ichem_start+p_vash_3)       &
                + tr(k,j,ichem_start+p_vash_4)       &
                + tr(k,j,ichem_start+p_vash_5)       &
                + tr(k,j,ichem_start+p_vash_6)       &
                + tr(k,j,ichem_start+p_vash_7)       &
                + tr(k,j,ichem_start+p_vash_8)       &
                + tr(k,j,ichem_start+p_vash_9)       &
                + tr(k,j,ichem_start+p_vash_10))     &
                *(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
            endif
            if (chem_opt == 317 ) then
              intash(j)=intash(j)+(tr(k,j,ichem_start+p_vash_1)       &
                + tr(k,j,ichem_start+p_vash_2)       &
                + tr(k,j,ichem_start+p_vash_3)       &
                + tr(k,j,ichem_start+p_vash_4))      &
                *(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
            endif
          enddo
          if (chem_opt == 316 .or. chem_opt == 317 ) intash(j)=intash(j)/dpsum
          intaer(j)=intaer(j)/dpsum
          intbc(j)=intbc(j)/dpsum
          intoc(j)=intoc(j)/dpsum
          intsulf(j)=intsulf(j)/dpsum
          intdust(j)=intdust(j)/dpsum
        enddo
!SMS$PARALLEL END
      endif ! chem_opt >= 300 .and. chem_opt < 500
      if (chem_opt == 502) then
!SMS$PARALLEL(dh, j) BEGIN
        do j=1,nip
          dpsum=0.
          intash(j)=0.
            do k=1,nvl
              dpsum=dpsum+(pr3d(k,j)-pr3d(k+1,j))
              rho_phys(k,j)=.5*(pr3d(k,j)+pr3d(k+1,j))&
              /(RD*tk3d(k,j)) !*(1.+.608*qv3d(k,j))
              intash(j)=intash(j)+(tr(k,j,ichem_start+p_vash_1)       &
                + tr(k,j,ichem_start+p_vash_2)       &
                + tr(k,j,ichem_start+p_vash_3)       &
                + tr(k,j,ichem_start+p_vash_4))      &
                *(pr3d(k,j)-pr3d(k+1,j))*rho_phys(k,j)
          enddo
          intash(j)=intash(j)/dpsum
        enddo
!SMS$PARALLEL END
      endif ! chem_opt = 502

!     if (mp_physics.eq.2) then
!       call icosio_out(its,time,'qcct',qcct,nvl, filename('qcct',its), header('qcct',nvl,its))
!       call icosio_out(its,time,'qrct',qrct,nvl, filename('qrct',its), header('qrct',nvl,its))
!       call icosio_out(its,time,'qict',qict,nvl, filename('qict',its), header('qict',nvl,its))
!       call icosio_out(its,time,'qsct',qsct,nvl, filename('qsct',its), header('qsct',nvl,its))
!     endif
      if ((mp_physics.ne.0).or.(cu_physics.ne.0)) then
        call icosio_out(its,time,'trco',trco,nvl, filename('trco',its), header('trco',nvl,its))
      endif
      if (chem_opt == 500) then
!SMS$PARALLEL(dh, j) BEGIN
           if(its.gt.1)tr1_tavg(:,:) = tr1_tavg(:,:)/float(ArchvStep-1)
!SMS$PARALLEL END 
        call icosio_out(its,time,'c13D',d1st,nvl,filename('c13D',its), header('c13D',nvl,its))
        call icosio_out(its,time,'c23D',tr1_tavg,nvl,filename('c23D',its), header('c23D',nvl,its))
        tr1_tavg(:,:) = 0.
      endif 
      if (chem_opt.ge.300 .and. chem_opt.lt.500) then
        call icosio_out(its,time,'ex3D',exch,nvl, filename('ex3D',its), header('ex3D',nvl,its))
        call icosio_out(its,time,'pm25',pm25,nvl, filename('pm25',its), header('pm25',nvl,its))
        call icosio_out(its,time,'pm10',p10,nvl, filename('pm10',its), header('pm10',nvl,its))
        call icosio_out(its,time,'dms1',dms1,nvl, filename('dms1',its), header('dms1',nvl,its))
        call icosio_out(its,time,'d1st',d1st,nvl, filename('d1st',its), header('d1st',nvl,its))
        call icosio_out(its,time,'d2st',d2st,nvl, filename('d2st',its), header('d2st',nvl,its))
        call icosio_out(its,time,'s1ea',sea1,nvl, filename('s1ea',its), header('s1ea',nvl,its))
        call icosio_out(its,time,'s2ea',sea2,nvl, filename('s2ea',its), header('s2ea',nvl,its))
        if (chem_opt.eq.300) then
          call icosio_out(its,time,'s3ea',sea3,nvl, filename('s3ea',its), header('s3ea',nvl,its))
          call icosio_out(its,time,'s4ea',sea4,nvl, filename('s4ea',its), header('s4ea',nvl,its))
          call icosio_out(its,time,'d3st',d3st,nvl, filename('d3st',its), header('d3st',nvl,its))
          call icosio_out(its,time,'d4st',d4st,nvl, filename('d4st',its), header('d4st',nvl,its))
          call icosio_out(its,time,'d5st',d5st,nvl, filename('d5st',its), header('d5st',nvl,its))
        endif
!TBH:  Added these IGNOREs to work around a PPP core dump.  Mark Govett is 
!TBH:  investigating...  
!SMS$IGNORE BEGIN
        dms1(:,:) = tr(:,:,ichem_start+p_bc1)
!SMS$IGNORE END
        call icosio_out(its,time,'pbc1',dms1,nvl, filename('pbc1',its), header('pbc1',nvl,its))
!SMS$IGNORE BEGIN
        dms1(:,:) = tr(:,:,ichem_start+p_bc2)
!SMS$IGNORE END
        call icosio_out(its,time,'pbc2',dms1,nvl, filename('pbc2',its), header('pbc2',nvl,its))
!SMS$IGNORE BEGIN
        dms1(:,:) = tr(:,:,ichem_start+p_oc1)
!SMS$IGNORE END
        call icosio_out(its,time,'obc1',dms1,nvl, filename('obc1',its), header('obc1',nvl,its))
!SMS$IGNORE BEGIN
        dms1(:,:) = tr(:,:,ichem_start+p_oc2)
!SMS$IGNORE END
        call icosio_out(its,time,'obc2',dms1,nvl, filename('obc2',its), header('obc2',nvl,its))
!SMS$IGNORE BEGIN
        dms1(:,:) = tr(:,:,ichem_start+p_sulf)
!SMS$IGNORE END
        call icosio_out(its,time,'sulf',dms1,nvl, filename('sulf',its), header('sulf',nvl,its))
!SMS$IGNORE BEGIN
        dms1(:,:) = tr(:,:,ichem_start+p_so2)
!SMS$IGNORE END
        call icosio_out(its,time,'pso2',dms1,nvl, filename('pso2',its), header('pso2',nvl,its))
!SMS$IGNORE BEGIN
        dms1(:,:) = tr(:,:,ichem_start+p_msa)
!SMS$IGNORE END
        call icosio_out(its,time,'pmsa',dms1,nvl, filename('pmsa',its), header('pmsa',nvl,its))
!SMS$IGNORE BEGIN
        dms1(:,:) = tr(:,:,ichem_start+p_p25)
!SMS$IGNORE END
        call icosio_out(its,time,'pp25',dms1,nvl, filename('pp25',its), header('pp25',nvl,its))
!SMS$IGNORE BEGIN
        dms1(:,:) = tr(:,:,ichem_start+p_p10)
!SMS$IGNORE END
        call icosio_out(its,time,'pp10',dms1,nvl, filename('pp10',its), header('pp10',nvl,its))
        if (chem_opt.eq.316.or.chem_opt.eq.317) then
          print *,'p_vash_1,p_vash_4 = ',p_vash_1,p_vash_4
!SMS$IGNORE BEGIN
          dms1(:,:) = tr(:,:,ichem_start+p_vash_1)
!SMS$IGNORE END
          call icosio_out(its,time,'ash1',dms1,nvl, filename('ash1',its), header('ash1',nvl,its))
!SMS$IGNORE BEGIN
          dms1(:,:) = tr(:,:,ichem_start+p_vash_2)
!SMS$IGNORE END
          call icosio_out(its,time,'ash2',dms1,nvl, filename('ash2',its), header('ash2',nvl,its))
!SMS$IGNORE BEGIN
          dms1(:,:) = tr(:,:,ichem_start+p_vash_3)
!SMS$IGNORE END
          call icosio_out(its,time,'ash3',dms1,nvl, filename('ash3',its), header('ash3',nvl,its))
!SMS$IGNORE BEGIN
          dms1(:,:) = tr(:,:,ichem_start+p_vash_4)
!SMS$IGNORE END
          call icosio_out(its,time,'ash4',dms1,nvl, filename('ash4',its), header('ash4',nvl,its))
          if (chem_opt.eq.316) then
!SMS$IGNORE BEGIN
            dms1(:,:) = tr(:,:,ichem_start+p_vash_5)
!SMS$IGNORE END
            call icosio_out(its,time,'ash5',dms1,nvl, filename('ash5',its), header('ash5',nvl,its))
!SMS$IGNORE BEGIN
            dms1(:,:) = tr(:,:,ichem_start+p_vash_6)
!SMS$IGNORE END
            call icosio_out(its,time,'ash6',dms1,nvl, filename('ash6',its), header('ash6',nvl,its))
!SMS$IGNORE BEGIN
            dms1(:,:) = tr(:,:,ichem_start+p_vash_7)
!SMS$IGNORE END
            call icosio_out(its,time,'ash7',dms1,nvl, filename('ash7',its), header('ash7',nvl,its))
!SMS$IGNORE BEGIN
            dms1(:,:) = tr(:,:,ichem_start+p_vash_8)
!SMS$IGNORE END
            call icosio_out(its,time,'ash8',dms1,nvl, filename('ash8',its), header('ash8',nvl,its))
!SMS$IGNORE BEGIN
            dms1(:,:) = tr(:,:,ichem_start+p_vash_9)
!SMS$IGNORE END
            call icosio_out(its,time,'ash9',dms1,nvl, filename('ash9',its), header('ash9',nvl,its))
!SMS$IGNORE BEGIN
            dms1(:,:) = tr(:,:,ichem_start+p_vash_10)
!SMS$IGNORE END
            call icosio_out(its,time,'ash0',dms1,nvl, filename('ash0',its), header('ash0',nvl,its))
          endif !chem_opt=316
        endif !chem_opt=316 or chem_opt=317
      endif !chem_opt.ge.300 .and. chem_opt.lt.500
      if (chem_opt == 500) then
        call icosio_out(its,time,'fl2D',intaer,1,filename('2D__',its), header('fl2D',1,its))
      endif
! output for volcanic ash only
      if (chem_opt.eq.502) then
          call icosio_out(its,time,'iash',intash,1, filename('2D__',its), header('iash',1,its))
          print *,'p_vash_1,p_vash_4 = ',p_vash_1,p_vash_4
!SMS$IGNORE BEGIN
          dms1(:,:) = tr(:,:,ichem_start+p_vash_1)
!SMS$IGNORE END
          call icosio_out(its,time,'ash1',dms1,nvl, filename('ash1',its), header('ash1',nvl,its))
!SMS$IGNORE BEGIN
          dms1(:,:) = tr(:,:,ichem_start+p_vash_2)
!SMS$IGNORE END
          call icosio_out(its,time,'ash2',dms1,nvl, filename('ash2',its), header('ash2',nvl,its))
!SMS$IGNORE BEGIN
          dms1(:,:) = tr(:,:,ichem_start+p_vash_3)
!SMS$IGNORE END
          call icosio_out(its,time,'ash3',dms1,nvl, filename('ash3',its), header('ash3',nvl,its))
!SMS$IGNORE BEGIN
          dms1(:,:) = tr(:,:,ichem_start+p_vash_4)
!SMS$IGNORE END
          call icosio_out(its,time,'ash4',dms1,nvl, filename('ash4',its), header('ash4',nvl,its))
      endif
      if (chem_opt.ge.300 .and. chem_opt.lt.500) then
        call icosio_out(its,time,'ia2D',intaer,1, filename('2D__',its), header('ia2D',1,its))
        call icosio_out(its,time,'ib2D',intbc,1, filename('2D__',its), header('ib2D',1,its))
        call icosio_out(its,time,'io2D',intoc,1, filename('2D__',its), header('io2D',1,its))
        call icosio_out(its,time,'is2D',intsulf,1, filename('2D__',its), header('is2D',1,its))
        call icosio_out(its,time,'id2D',intdust,1, filename('2D__',its), header('id2D',1,its))
        call icosio_out(its,time,'ao2D',aod2d,1, filename('2D__',its), header('ao2D',1,its))
        if (chem_opt.eq.316.or.chem_opt.eq.317) then
          call icosio_out(its,time,'iash',intash,1, filename('2D__',its), header('iash',1,its))
        endif
      endif !chem_opt.ge.300 .and. chem_opt.lt.500
    endif

  end subroutine chem_output

end module module_chem_output
