!#include "../../../ESMFVersionDefine.h"

!
! !module: gfs_physics_getcf_mod --- configure data module of the 
!                          gridded component of the gfs physics.
!
! !description: this program uses the esmf configure class software
!               to get all parameters of the original namelists of
!               gfs physics and all trun-on/turn-off switch flags 
!               import state and the export state.
!
! !revision history:
!
!  november 2004 weiyu yang    initial code for esmf gfs model.
!  may      2005 weiyu yang    for the updated gfs version.
!  march    2006 s. moorthi    modified for the new gfs.
!  january  2007 h. juang      modified for the gfs dynamics.
!  july     2007 s. moorthi    modified for the gfs physics.
!  november 2007 h. juang      continue to finish gfs physics
!  september2011 w. yang       modified for using the ESMF 5.2.0r library.
!  may      2014 p. pegion     added stochastic physics variables
!
! !interface:
!
      module gfs_physics_getcf_mod

!
!!uses:
!

! the esmf internal state contents.
!----------------------------------
      use gfs_physics_internal_state_mod, ONLY: gfs_physics_internal_state
      USE ESMF,                           ONLY: esmf_gridcomp, esmf_vm,              &
                                                esmf_config, esmf_success,           &
                                                esmf_gridcompget, esmf_vmgetglobal,  &
                                                esmf_vmget, esmf_configgetattribute, &
!jw
                                                esmf_configFindLabel,                &
                                                esmf_vmgetcurrent
      USE gfs_physics_err_msg_mod,        ONLY: gfs_physics_err_msg_final,           &
                                                gfs_physics_err_msg,                 &
                                                gfs_physics_err_msg_var

      implicit none

      contains

!------------------------------------------------------------------------
!bop
!
! !iroutine:
!
! !interface:
!
     subroutine gfs_physics_getcf(gc_gfs_phy, int_state, rc)

      type(esmf_gridcomp),               intent(inout) :: gc_gfs_phy 
      type(gfs_physics_internal_state),  intent(inout) :: int_state 
      type(esmf_vm)                                     :: vm           
      type(esmf_vm)                                     :: vm_local    
      integer, optional,                  intent(out)   :: rc  
!
! !description: load run parameters from the resource file into internal
!               state.

!
! !revision history:
!
!  november 3, 2004 weiyu yang, for ncep gfs model.
!  may         2005 weiyu yang  for the updated gfs version.
!  january     2007 hann-ming henry juang  for the  gfs physics version.
!  july        2007 shrinivas moorthi      for the  gfs physics version.
!  november    2007 hann-ming henry juang  continue for gfs physics
!
!eop
!------------------------------------------------------------------------

      type(esmf_config) :: cf
      integer           :: ndat, i, rc1, rcfinal
      character*20      :: clab, pelab
      integer           :: pe_member,j, i1, me, tasks, totpe

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

! get the esmf config from the gfs grid component.
!-------------------------------------------------
      call esmf_gridcompget(gc_gfs_phy, config = cf, rc = rc1)

      call gfs_physics_err_msg(rc1,				&
               'gfs physics get configuration',rcfinal)

! start to get the configure data.
!---------------------------------

! the original gfs physics has namelist file, "nam_gfs_phy".
! besides to get the information from namelist file, the getcf
! process will get the "esmf_state_namelist" which contains all turn-on/turn-off
! switch flags for all possible fields in the esmf import and export states.
! the single esmf config file contains all information of the namelist files.
!-------------------------------------------------------------------------------

! for "nam_gfs_phy".
!---------------

! "esmf_configgetattribute" is an esmf config utility routine which is used
! to get information from the config file.  
! the first argument is the gfs esmf config.
! the second argument is the destination variable 
! to which the value of the namelist 
! variable will be sent. the third one is the label 
! in the config file which is used
! to identify the required namelist variable.  
! the last one is the error signal variable.
! "esmf_configgetattribute" is a generic interface name 
! which can be used to get the 
! the different data types data from the config file, 
! such as real(4), real(8), 
! integer(4), integer(8), character(...), etc..
!------------------------------------------------------------------------------
      call esmf_vmgetglobal(vm, rc = rc1)
      call esmf_vmget(vm, localpet = me, pecount = tasks, rc = rc1)

      call esmf_configgetattribute(cf, 					&
                             int_state%nam_gfs_phy%total_member,   	&
                             label = 'total_member:', rc    = rc)

      totpe = 0
      i1 = 0
      do j=1,int_state%nam_gfs_phy%total_member
        write(pelab,'("PE_MEMBER",i2.2,":")') j
        call esmf_configgetattribute(cf, 				&
                             pe_member, label = pelab, rc = rc)
        if (pe_member == 0) 						&
            pe_member = tasks / int_state%nam_gfs_phy%total_member
        totpe = totpe + pe_member
        do i = 1, pe_member
          if (me == i1) then
             int_state%nam_gfs_phy%member_id = j
          end if
          i1 = i1+1
        end do
      end do
      if (totpe /= tasks) then
       print *,' totpe=',totpe,' and tasks=',tasks, ' do not match'
       stop 9999
      endif

      call esmf_vmgetcurrent(vm_local, rc = rc1)

      call esmf_configgetattribute(cf, 					&
                              int_state%nam_gfs_phy%nlunit,         	&
                              label = 'nlunit:',  rc = rc1)
      call gfs_physics_err_msg_var(rc1,				        &
               'gfs physics getcf','nlunit',rcfinal)


      call esmf_configgetattribute(cf, 					&
                              int_state%nam_gfs_phy%gfs_phy_namelist, 	&
                              label = 'namelist:',  rc = rc1)
      call gfs_physics_err_msg_var(rc1,				        &
                 'gfs physics getcf','namelist',rcfinal)

      call esmf_configgetattribute(cf, 					&
                              int_state%nam_gfs_phy%deltim,         	&
                              label = 'deltim:', rc = rc1)
      int_state%deltim = int_state%nam_gfs_phy%deltim

      call gfs_physics_err_msg_var(rc1,				        &
               'gfs physics getcf','deltim',rcfinal)

      call esmf_configgetattribute(cf,                                 &
                              int_state%restart_run,                   &
                              label = 'restart:',  rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','restart_run',rcfinal)

      if (int_state%nam_gfs_phy%total_member <= 1) then
        int_state%nam_gfs_phy%sfc_ini  = 'sfc_ini'
        int_state%nam_gfs_phy%grd_ini  = 'grid_ini'
        int_state%nam_gfs_phy%nst_ini  = 'nst_ini'
      else
        write(int_state%nam_gfs_phy%sfc_ini, 				&
              '("sfc_ini_",i2.2)') int_state%nam_gfs_phy%member_id
        write(int_state%nam_gfs_phy%grd_ini, 				&
              '("grid_ini_",i2.2)') int_state%nam_gfs_phy%member_id
        write(Int_State%nam_gfs_phy%nst_ini,                            &
              '("nst_ini_",I2.2)') Int_State%nam_gfs_phy%member_id
      endif

!
! get output file name
!--------------------------
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%num_file,                       &
                              label = 'num_file:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','num_file',rcfinal)
!
      allocate(int_state%filename_base(int_state%num_file))
      call esmf_configFindLabel(CF,'filename_base:',rc=RC)
      Do I=1,int_state%num_file
        call esmf_configgetattribute(cf,                                &
                              int_state%filename_base(i),               &
                              rc = rc1)
      enddo
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','filename_base',rcfinal)
!
! for "esmf_state_namelist"
!--------------------------
!
! ----------- upair to couple with dynamics ----------------
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%idate1_import,    &
                              label = 'idate1_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','idate1_import',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%z_import,         &
                              label = 'z_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','z_import',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%ps_import,        &
                              label = 'ps_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','ps_import',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%u_import,         &
                              label = 'u_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','u_import',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%v_import,         &
                              label = 'v_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','v_import',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%temp_import,      &
                              label = 'temp_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','temp_import',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%tracer_import,    &
                              label = 'tracer_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','tracer_import',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%p_import, &
                              label = 'p_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','p_import',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%dp_import,        &
                              label = 'dp_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','dp_import',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%dpdt_import,      &
                              label = 'dpdt_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','dpdt_import',rcfinal)

! Added by PJP for stochastic physics

      call esmf_configgetattribute(cf,                                 &
                              int_state%esmf_sta_list%shum_wts_import, &
                              label = 'shum_wts_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','shum_wts_import',rcfinal)
                                                                               
      call esmf_configgetattribute(cf,                                 &
                              int_state%esmf_sta_list%sppt_wts_import, &
                              label = 'sppt_wts_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','sppt_wts_import',rcfinal)
                                                                               
      call esmf_configgetattribute(cf,                                 &
                              int_state%esmf_sta_list%skeb_wts_import, &
                              label = 'skeb_wts_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','skeb_wts_import',rcfinal)
                                                                               
      call esmf_configgetattribute(cf,                                 &
                              int_state%esmf_sta_list%vc_wts_import,   &
                              label = 'vc_wts_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','vc_wts_import',rcfinal)
                                                                               
                                                                               
!----------------------------------
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%idate1_export,    &
                              label = 'idate1_export:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','idate1_export',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%z_export,         &
                              label = 'z_export:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','z_export',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%ps_export,        &
                              label = 'ps_export:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','ps_export',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%u_export,         &
                              label = 'u_export:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','u_export',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%v_export,         &
                              label = 'v_export:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','v_export',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%temp_export,      &
                              label = 'temp_export:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','temp_export',rcfinal)

      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%tracer_export,    &
                              label = 'tracer_export:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','tracer_export',rcfinal)

      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%p_export, &
                              label = 'p_export:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','p_export',rcfinal)

      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%dp_export,        &
                              label = 'dp_export:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','dp_export',rcfinal)

      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%dpdt_export,      &
                              label = 'dpdt_export:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','dpdt_export',rcfinal)

!
! ----------- surface  to couple with IO and others -----
      call esmf_configgetattribute(cf,                                  &
                int_state%esmf_sta_list%orography_import,               &
                label = 'orography_import:',    rc = rc1)
      call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf','orography_import',rcfinal)
! 
       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%t_skin_import,                 &
                 label = 't_skin_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf','t_skin_import',rcfinal)
!                               
       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%soil_mois_import,              &
                 label = 'soil_mois_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','soil_mois_import',rcfinal)
                               
       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%snow_depth_import,             &
                 label = 'snow_depth_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','snow_depth_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%soil_t_import,                 &
                 label = 'soil_t_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','soil_t_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%deep_soil_t_import,            &
                 label = 'deep_soil_t_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','deep_soil_t_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%roughness_import,              &
                 label = 'roughness_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','roughness_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%conv_cloud_cover_import,       &
                 label = 'conv_cloud_cover_import:', rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','conv_cloud_cover_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%conv_cloud_base_import,        &
                 label = 'conv_cloud_base_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','conv_cloud_base_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%conv_cloud_top_import,         &
                 label = 'conv_cloud_top_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','conv_cloud_top_import',rcfinal)
       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%albedo_visible_scattered_import,&
                 label = 'albedo_visible_scattered_import:',    rc=rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'albedo_visible_scattered_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%albedo_visible_beam_import,    &
                 label = 'albedo_visible_beam_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'albedo_visible_beams_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%albedo_nearir_scattered_import,&
                 label = 'albedo_nearir_scattered_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'albedo_visible_scattered_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%albedo_nearir_beam_import,     &
                 label = 'albedo_nearir_beam_import:',    rc = rc1) 
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'albedo_nearir_beam_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%sea_level_ice_mask_import,     &
                 label = 'sea_level_ice_mask_import:',    rc = rc1) 
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'sea_level_ice_mask_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%vegetation_cover_import,       &
                 label = 'vegetation_cover_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'vegetation_cover_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%canopy_water_import,           &
                 label = 'canopy_water_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'canopy_water_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%m10_wind_fraction_import,      &
                 label = 'm10_wind_fraction_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'm10_wind_fraction_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%vegetation_type_import,        &
                 label = 'vegetation_type_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'vegetation_type_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%soil_type_import,              &
                 label = 'soil_type_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'soil_type_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%zeneith_angle_facsf_import,    &
                 label = 'zeneith_angle_facsf_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'zeneith_angle_facsf_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%zeneith_angle_facwf_import,    &
                 label = 'zeneith_angle_facwf_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'zeneith_angle_facwf_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%uustar_import,                 &
                 label = 'uustar_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'uustar_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%ffmm_import,                   &
                 label = 'ffmm_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'ffmm_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%ffhh_import,                   &
                 label = 'ffhh_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf','ffhh_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%sea_ice_thickness_import,      &
                 label = 'sea_ice_thickness_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'sea_ice_thickness_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%sea_ice_concentration_import,  &
                 label = 'sea_ice_concentration_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'sea_ice_concentration_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%tprcp_import,                  &
                 label = 'tprcp_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'tprcp_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%srflag_import,                 &
                 label = 'srflag_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'srflag_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%actual_snow_depth_import,      &
                 label = 'actual_snow_depth_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                 &
                'gfs physics getcf',                                    &
                'actual_snow_depth_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%liquid_soil_moisture_import,   &
                 label = 'liquid_soil_moisture_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'liquid_soil_moisture_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%vegetation_cover_min_import,   &
                 label = 'vegetation_cover_min_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'vegetation_cover_min_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%vegetation_cover_max_import,   &
                 label = 'vegetation_cover_max_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'vegetation_cover_max_import',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%slope_type_import,             &
                 label = 'slope_type_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'slope_type_import',rcfinal)

       call esmf_configgetattribute(cf,                                  &
                 int_state%esmf_sta_list%snow_albedo_max_import,         &
                 label = 'snow_albedo_max_import:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'snow_albedo_max_import',rcfinal)

! --------------

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%orography_export,              &
                 label = 'orography_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'orography_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%t_skin_export,                 &
                 label = 't_skin_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                't_skin_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%soil_mois_export,              &
                 label = 'soil_mois_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf','soil_mois_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%snow_depth_export,             &
                 label = 'snow_depth_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'snow_depth_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%soil_t_export,                 &
                label = 'soil_t_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf','soil_t_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%deep_soil_t_export,            &
                 label = 'deep_soil_t_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'deep_soil_t_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%roughness_export,              &
                 label = 'roughness_export:', rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf','roughness_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%conv_cloud_cover_export,       &
                label = 'conv_cloud_cover_export:', rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'conv_cloud_cover_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%conv_cloud_base_export,        &
                 label = 'conv_cloud_base_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','conv_cloud_base_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%conv_cloud_top_export,         &
                 label = 'conv_cloud_top_export:', rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
          'gfs physics getcf','conv_cloud_top_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%albedo_visible_scattered_export,&
                 label = 'albedo_visible_scattered_export:', rc= rc1)
       call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf',                                     &
               'albedo_visible_scattered_export',rcfinal)

       call esmf_configgetattribute(cf,                                  &
                 int_state%esmf_sta_list%albedo_visible_beam_export,     &
                 label = 'albedo_visible_beam_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                 &
               'gfs physics getcf',                                     &
               'albedo_visible_beam_export',rcfinal)

       call esmf_configgetattribute(cf,                                  &
                 int_state%esmf_sta_list%albedo_nearir_scattered_export, &
                 label = 'albedo_nearir_scattered_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                 &
                'gfs physics getcf',                                    &
                'albedo_nearir_scattered_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%albedo_nearir_beam_export,     &
                 label = 'albedo_nearir_beam_export:',    rc = rc1) 
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'albedo_nearir_beam_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%sea_level_ice_mask_export,     &
                 label = 'sea_level_ice_mask_export:',    rc = rc1) 
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'sea_level_ice_mask_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%vegetation_cover_export,       &
                 label = 'vegetation_cover_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'vegetation_cover_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%canopy_water_export,           &
                 label = 'canopy_water_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'canopy_water_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%m10_wind_fraction_export,      &
                 label = 'm10_wind_fraction_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'm10_wind_fraction_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%vegetation_type_export,        &
                 label = 'vegetation_type_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'vegetation_type_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%soil_type_export,              &
                 label = 'soil_type_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'soil_type_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%zeneith_angle_facsf_export,    &
                 label = 'zeneith_angle_facsf_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'zeneith_angle_facsf_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%zeneith_angle_facwf_export,    &
                 label = 'zeneith_angle_facwf_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'zeneith_angle_facwf_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%uustar_export,                 &
                 label = 'uustar_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'uustar_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%ffmm_export,                   &
                              label = 'ffmm_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf','ffmm_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%ffhh_export,                   &
                              label = 'ffhh_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','ffhh_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%sea_ice_thickness_export,      &
                 label = 'sea_ice_thickness_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'sea_ice_thickness_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%sea_ice_concentration_export,  &
                 label = 'sea_ice_concentration_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'sea_ice_concentration_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%tprcp_export,                  &
                 label = 'tprcp_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','tprcp_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%srflag_export,                 &
                 label = 'srflag_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'srflag_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%actual_snow_depth_export,      &
                 label = 'actual_snow_depth_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'actual_snow_depth_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%liquid_soil_moisture_export,   &
                 label = 'liquid_soil_moisture_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'liquid_soil_moisture_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%vegetation_cover_min_export,   &
                 label = 'vegetation_cover_min_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'vegetation_cover_min_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%vegetation_cover_max_export,   &
                 label = 'vegetation_cover_max_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
                'gfs physics getcf',                                    &
                'vegetation_cover_max_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%slope_type_export,             &
                 label = 'slope_type_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','slope_type_export',rcfinal)

       call esmf_configgetattribute(cf,                                 &
                 int_state%esmf_sta_list%snow_albedo_max_export,        &
                 label = 'snow_albedo_max_export:',    rc = rc1)
       call gfs_physics_err_msg_var(rc1,                                &
               'gfs physics getcf','snow_albedo_max_export',rcfinal)

!----------------------------------

      call gfs_physics_err_msg_final(rcfinal,'gfs_physics_getcf',rc)

      end subroutine gfs_physics_getcf

      end module gfs_physics_getcf_mod
