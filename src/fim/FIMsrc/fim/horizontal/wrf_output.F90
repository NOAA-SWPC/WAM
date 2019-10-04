module module_wrf_output

contains

  subroutine wrf_output(its,nts,pr3d,tk3d,tr,phys2dwrf)

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
    real,intent(in)::&
      pr3d(nvlp1,nip),&
      tk3d(nvl,nip),tr(nvl,nip,ntra+ntrb)
    real:: qcct(nvl,nip),qict(nvl,nip),qrct(nvl,nip),qsct(nvl,nip),&
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

      ichem_start=ntra+1
      imoist_start=ntra

!SMS$IGNORE BEGIN
!TBH:  Added this IGNORE to work around a PPP core dump between here and 
!TBH:  "SMS$IGNORE END".  Mark Govett is investigating...  

      if ((.not.mp_physics==0).or.(.not.cu_physics==0)) then
        trco(:,:) = 0.
!       trco(:,:) = tr(:,:,5)
        trco(1,:) = phys2dwrf(:,1)
        trco(2,:) = phys2dwrf(:,2)
        trco(3,:) = phys2dwrf(:,3)
        trco(4,:) = phys2dwrf(:,4)
        trco(5,:) = phys2dwrf(:,5)
        trco(6,:) = phys2dwrf(:,6)
      endif
      if (mp_physics == 2) then
        qcct(:,:) = tr(:,:,imoist_start+p_qc)
        qrct(:,:) = tr(:,:,imoist_start+p_qr)
        qict(:,:) = tr(:,:,imoist_start+p_qi)
        qsct(:,:) = tr(:,:,imoist_start+p_qs)
      endif
!SMS$IGNORE END

      if (mp_physics.eq.2) then
        call icosio_out(its,time,'qcct',qcct,nvl,filename('qcct',its),header('qcct',nvl,its))
        call icosio_out(its,time,'qrct',qrct,nvl,filename('qrct',its),header('qrct',nvl,its))
        call icosio_out(its,time,'qict',qict,nvl,filename('qict',its),header('qict',nvl,its))
        call icosio_out(its,time,'qsct',qsct,nvl,filename('qsct',its),header('qsct',nvl,its))
      endif
      if ((mp_physics.ne.0).or.(cu_physics.ne.0)) then
        call icosio_out(its,time,'trco',trco,nvl,filename('trco',its),header('trco',nvl,its))
      endif
    endif

  end subroutine wrf_output

end module module_wrf_output
