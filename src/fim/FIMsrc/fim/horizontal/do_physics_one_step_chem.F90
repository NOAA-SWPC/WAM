module module_do_physics_one_step_chem

contains

  subroutine do_physics_one_step_chem(ipn,skip_chem,skip_cu_physics,&
    skip_mp_physics,dkt,dq3dt,dt3dt)

    use gfs_physics_internal_state_mod,only:gis_phy
    use gfs_physics_sfc_flx_mod,only:flx_var_data
    use machine,only:kind_phys,kind_rad
    use module_control,only:nip
    use module_wrf_variables,only:phys3dwrf,exch,pb2d

    integer,intent(in)::ipn
    logical,intent(in)::skip_cu_physics,skip_mp_physics,skip_chem
    real(kind=kind_phys),intent(in)::dkt(:,:)
    real(kind=kind_rad),intent(in)::dq3dt(:,:,:),dt3dt(:,:,:)
    integer::ivl

    ! chem radiation stuff, for different midbands, just change this data
    ! statement in routine colum_chem/module_ar_ra
    !
    ! data midbands/.2,.235,.27,.2875,.3025,.305,.3625,.55,1.92,1.745,6.135/

    if (skip_cu_physics.or.skip_mp_physics) then
      pb2d(ipn)=gis_phy%flx_fld%hpbl(ipn,1)
      do ivl=1,gis_phy%levs-1
        phys3dwrf(ivl,ipn,5)=dt3dt(1,ivl,3)/gis_phy%deltim                  ! bl - t
        phys3dwrf(ivl,ipn,6)=(dt3dt(1,ivl,1)+dt3dt(1,ivl,2))/gis_phy%deltim ! ra ? - t
        phys3dwrf(ivl,ipn,2)= dq3dt(1,ivl,1)/gis_phy%deltim                 ! bl - qv
      enddo
    endif
    if (skip_chem) then
      pb2d(ipn)=gis_phy%flx_fld%hpbl(ipn,1)
      do ivl=1,gis_phy%levs-1
        exch(ivl,ipn)=dkt(1,ivl)
      enddo
    endif
  end subroutine do_physics_one_step_chem

end module module_do_physics_one_step_chem
