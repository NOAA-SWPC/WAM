!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write_restart_phy: write physics fields to the restart file
! SMS doesn't yet properly handle Fortran derived types, so those fields
! need to go through an interface routine (writearr64).
!
! !!!!!CRITICAL!!!!! Any changes to fields read in here MUST be made in exactly the same way
! in read_restart_phy.F90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_restart_phy (unitno)
  use module_control,                 only: nip
  use module_sfc_variables
  use module_wrf_control,             only: nbands
!SMS$IGNORE BEGIN
  use gfs_physics_internal_state_mod, only: gis_phy
!SMS$IGNORE END
  use module_chem_variables,          only: ext_cof, asymp, extlw_cof

  implicit none

! Input arguments
  integer, intent(in) :: unitno   ! Unit number to write to

! Local workspace
  integer :: ipn, i, j, n, k      ! Loop indices
  integer :: ims, ime             ! memory bounds
  integer :: ips, ipe             ! owned memory bounds
! Transposed versions of arrays so writearr64 can tell SMS to do the right thing.
  real*8, allocatable :: smc_loc(:,:)
  real*8, allocatable :: stc_loc(:,:)
  real*8, allocatable :: slc_loc(:,:)
  real*8, allocatable :: hprime_loc(:,:)
  real*8, allocatable :: fluxr_loc(:,:)
!SMS$DISTRIBUTE(dh,1) BEGIN
  real*8, allocatable :: ext_cof_loc(:,:)
  real*8, allocatable :: extlw_cof_loc(:,:)
  real*8, allocatable :: asymp_loc(:,:)
  real*8, allocatable :: arr(:)
!SMS$DISTRIBUTE END

  ims = gis_phy%ims
  ime = gis_phy%ime
  ips = gis_phy%ips
  ipe = gis_phy%ipe

!SMS$SERIAL BEGIN
  write (unitno, err=90) cv2d, cvt2d, cvb2d, slmsk2d, ts2d, st3d, sheleg2d, zorl2d, snoalb2d
  write (unitno, err=90) hprm2d, hice2d, fice2d, tprcp2d, srflag2d, slc3d, sm3d, snwdph2d
  write (unitno, err=90) slope2d, shdmin2d, shdmax2d, tg32d, vfrac2d, canopy2d, vtype2d, stype2d
  write (unitno, err=90) f10m2d, ffmm2d, ffhh2d, alvsf2d, alnsf2d, alvwf2d, alnwf2d, facsf2d
  write (unitno, err=90) facwf2d, t2m2d, q2m2d
  write(6,*)'write_restart_phy: ts2d, us2d, hf2d, qf2d, sw2d, lw2d=', &
            ts2d(1), us2d(1), hf2d(1), qf2d(1), sw2d(1), lw2d(1)
  write (unitno, err=90) ts2d, us2d, hf2d, qf2d, sw2d, lw2d
!SMS$SERIAL END

! Allocate space for arrays which need to be transposed

  allocate (smc_loc(ims:ime,gis_phy%lsoil))
  allocate (stc_loc(ims:ime,gis_phy%lsoil))
  allocate (slc_loc(ims:ime,gis_phy%lsoil))
  allocate (hprime_loc(ims:ime,gis_phy%nmtvr))
  allocate (fluxr_loc(ims:ime,gis_phy%nfxr))

  allocate (ext_cof_loc(nip,nbands))
  allocate (extlw_cof_loc(nip,16))
  allocate (asymp_loc(nip,nbands))

! Write out ALL sfc_fld items. Some are definitely needed, but the list is huge
! Must call an interface routine (writearr64) until SMS can handle derived types
! 2nd arg to writearr64 is size of dimensions after ipn
! Need to transpose some fields to ipn as 1st index so writearr64 can use SMS 
! to do the right thing.

  do j=1,gis_phy%lsoil
    do ipn=ips,ipe
      smc_loc(ipn,j) = gis_phy%sfc_fld%smc(j,ipn,1)
      stc_loc(ipn,j) = gis_phy%sfc_fld%stc(j,ipn,1)
      slc_loc(ipn,j) = gis_phy%sfc_fld%slc(j,ipn,1)
    end do
  end do

  call writearr64 (smc_loc, gis_phy%lsoil, unitno)
  call writearr64 (stc_loc, gis_phy%lsoil, unitno)
  call writearr64 (slc_loc, gis_phy%lsoil, unitno)

  call writearr64 (gis_phy%sfc_fld%tsea,   1,  unitno)
  call writearr64 (gis_phy%sfc_fld%sheleg, 1,  unitno)
  call writearr64 (gis_phy%sfc_fld%sncovr, 1,  unitno)
  call writearr64 (gis_phy%sfc_fld%tg3,    1,  unitno)
  call writearr64 (gis_phy%sfc_fld%zorl,   1,  unitno)
  call writearr64 (gis_phy%sfc_fld%cv,     1,  unitno)
  call writearr64 (gis_phy%sfc_fld%cvb,    1,  unitno)
  call writearr64 (gis_phy%sfc_fld%cvt,    1,  unitno)
  write(6,*)'write_restart_phy: writing alvsf'
  call writearr64 (gis_phy%sfc_fld%alvsf,  1,  unitno)
  call writearr64 (gis_phy%sfc_fld%alvwf,  1,  unitno)
  call writearr64 (gis_phy%sfc_fld%alnsf,  1,  unitno)
  call writearr64 (gis_phy%sfc_fld%alnwf,  1,  unitno)
  call writearr64 (gis_phy%sfc_fld%slmsk,  1,  unitno)
  call writearr64 (gis_phy%sfc_fld%vfrac,  1,  unitno)
  call writearr64 (gis_phy%sfc_fld%canopy, 1,  unitno)
  call writearr64 (gis_phy%sfc_fld%f10m,   1,  unitno)
  call writearr64 (gis_phy%sfc_fld%t2m,    1,  unitno)
  call writearr64 (gis_phy%sfc_fld%q2m,    1,  unitno)
  call writearr64 (gis_phy%sfc_fld%vtype,  1,  unitno)
  call writearr64 (gis_phy%sfc_fld%stype,  1,  unitno)
  call writearr64 (gis_phy%sfc_fld%facsf,  1,  unitno)
  call writearr64 (gis_phy%sfc_fld%facwf,  1,  unitno)
  call writearr64 (gis_phy%sfc_fld%uustar, 1,  unitno)
  call writearr64 (gis_phy%sfc_fld%ffmm,   1,  unitno)
  call writearr64 (gis_phy%sfc_fld%ffhh,   1,  unitno)
  call writearr64 (gis_phy%sfc_fld%hice,   1,  unitno)
  call writearr64 (gis_phy%sfc_fld%fice,   1,  unitno)
  call writearr64 (gis_phy%sfc_fld%uustar, 1,  unitno)
  call writearr64 (gis_phy%sfc_fld%tprcp,  1,  unitno)
  call writearr64 (gis_phy%sfc_fld%srflag, 1,  unitno)
  call writearr64 (gis_phy%sfc_fld%snwdph, 1,  unitno)
  call writearr64 (gis_phy%sfc_fld%shdmin, 1,  unitno)
  call writearr64 (gis_phy%sfc_fld%shdmax, 1,  unitno)
  call writearr64 (gis_phy%sfc_fld%slope,  1,  unitno)
  call writearr64 (gis_phy%sfc_fld%snoalb, 1,  unitno)

! Write out ALL flx_fld items. Are all actually needed?
  call writearr64 (gis_phy%flx_fld%sfcdsw, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%coszen, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%tmpmin, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%tmpmax, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%dusfc,  1,  unitno)
  call writearr64 (gis_phy%flx_fld%dvsfc,  1,  unitno)
  call writearr64 (gis_phy%flx_fld%dtsfc,  1,  unitno)
  call writearr64 (gis_phy%flx_fld%dqsfc,  1,  unitno)
  call writearr64 (gis_phy%flx_fld%dlwsfc, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%ulwsfc, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%gflux,  1,  unitno)
  call writearr64 (gis_phy%flx_fld%runoff, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%ep,     1,  unitno)
  call writearr64 (gis_phy%flx_fld%cldwrk, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%dugwd,  1,  unitno)
  call writearr64 (gis_phy%flx_fld%dvgwd,  1,  unitno)
  call writearr64 (gis_phy%flx_fld%psmean, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%geshem, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%rainc,  1,  unitno)
  call writearr64 (gis_phy%flx_fld%evap,   1,  unitno)
  call writearr64 (gis_phy%flx_fld%hflx,   1,  unitno)
  call writearr64 (gis_phy%flx_fld%bengsh, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%sfcnsw, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%sfcdlw, 1,  unitno)
  call writearr64 (gis_phy%flx_fld%tsflw,  1,  unitno)
  call writearr64 (gis_phy%flx_fld%psurf,  1,  unitno)
  call writearr64 (gis_phy%flx_fld%u10m,   1,  unitno)
  call writearr64 (gis_phy%flx_fld%v10m,   1,  unitno)
  call writearr64 (gis_phy%flx_fld%hpbl,   1,  unitno)
  call writearr64 (gis_phy%flx_fld%pwat,   1,  unitno)

! These things are from cpl_dyn_to_phy. Not sure which are necessary
  call writearr64 (gis_phy%ps,   1,            unitno)
  call writearr64 (gis_phy%dp,   gis_phy%levs, unitno)
  call writearr64 (gis_phy%p,    gis_phy%levs, unitno)
  call writearr64 (gis_phy%u,    gis_phy%levs, unitno)
  call writearr64 (gis_phy%v,    gis_phy%levs, unitno)
  call writearr64 (gis_phy%dpdt, gis_phy%levs, unitno)
  call writearr64 (gis_phy%q,    gis_phy%levs, unitno)
  call writearr64 (gis_phy%oz,   gis_phy%levs, unitno)
  call writearr64 (gis_phy%cld,  gis_phy%levs, unitno)
  call writearr64 (gis_phy%t,    gis_phy%levs, unitno)

  do n=1,gis_phy%num_p3d
    do j=1,gis_phy%lats_node_r
      do i=1,gis_phy%nblck
        call writearr64 (gis_phy%phy_f3d(:,:,i,j,n), gis_phy%levs, unitno)
      end do
    end do
  end do

  do n=1,gis_phy%num_p2d
    do j=1,gis_phy%lats_node_r
      call writearr64 (gis_phy%phy_f2d(:,j,n), 1, unitno)
    end do
  end do

! These things are in gis_phy proper, not in substructures sfc_fld or flx_fld
  do j=1,gis_phy%lats_node_r
    do ipn=ips,ipe
      do n=1,gis_phy%nmtvr
        hprime_loc(ipn,n) = gis_phy%hprime(n,ipn,j)
      end do
    end do
    call writearr64 (hprime_loc, gis_phy%nmtvr, unitno)
    call writearr64 (gis_phy%coszdg(:,j), 1, unitno)
    call writearr64 (gis_phy%sfalb(:,j),  1, unitno)
    call writearr64 (gis_phy%slag(:,j),   1, unitno)
    call writearr64 (gis_phy%sdec(:,j),   1, unitno)
    call writearr64 (gis_phy%cdec(:,j),   1, unitno)
    do n=1,gis_phy%nblck
      call writearr64 (gis_phy%swh(:,:,n,j), gis_phy%levs, unitno)
      call writearr64 (gis_phy%hlw(:,:,n,j), gis_phy%levs, unitno)
    end do
    do ipn=ips,ipe
      do n=1,gis_phy%nfxr
        fluxr_loc(ipn,n) = gis_phy%fluxr(n,ipn,j)
      end do
    end do
    call writearr64 (fluxr_loc, gis_phy%nfxr, unitno)
  end do

!SMS$PARALLEL(dh, ipn) BEGIN
! These things are from chemistry, but are used in grrad. Need to transpose so
! writearr64 tells SMS to do the right thing.
  do k=1,gis_phy%levs
    do n=1,nbands
      do ipn=1,nip
        ext_cof_loc(ipn,n) = ext_cof(k,ipn,n)
        asymp_loc(ipn,n) = asymp(k,ipn,n)
      end do
    end do

    call writearr64 (ext_cof_loc, nbands, unitno)
    call writearr64 (asymp_loc,   nbands, unitno)

! The "16" is hard-wired in the allocation done in dyn_alloc
    do n=1,16
      do ipn=1,nip
        extlw_cof_loc(ipn,n) = extlw_cof(k,ipn,n)
      end do
    end do
    call writearr64 (extlw_cof_loc, 16, unitno)
  end do
  
  write (6,*) 'write_restart_phy: successfully wrote physics fields to restart file'
!SMS$PARALLEL END

  deallocate (smc_loc)
  deallocate (stc_loc)
  deallocate (slc_loc)
  deallocate (ext_cof_loc)
  deallocate (extlw_cof_loc)
  deallocate (asymp_loc)
  deallocate (hprime_loc)
  deallocate (fluxr_loc)

  return

90 write(6,*)'write_restart_phy: Error writing to unit ', unitno, '. Stopping'
  stop
end subroutine write_restart_phy
