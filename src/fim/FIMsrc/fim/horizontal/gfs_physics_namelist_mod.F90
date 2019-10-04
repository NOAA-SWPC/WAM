!
! !module: gfs_physics_namelist_mod  ---       definition of the name list
!                                               in the esmf internal state.
!
! !description:   define the name list variables
!                 in the esmf internal state.
!---------------------------------------------------------------------------
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  february 2006      took out model namelists
!  january  2007      hann-ming henry juang for gfs dynamics only
!  july     2007      shrinivas moorthi for gfs physics only
!  november 2007      hann-ming henry juang for gfs physics
!
! !interface:
!
      module gfs_physics_namelist_mod


      implicit none

      type nam_gfs_phy_namelist
           integer                :: nlunit, total_member, member_id
           real                   :: deltim
           character(80)          :: gfs_phy_namelist
           character(20)          :: sfc_ini
      end type nam_gfs_phy_namelist
!
      type gfs_phy_state_namelist
!
!For couple between dynamics and physics
!---------------------------------------
           integer                :: idate1_import
           integer                :: z_import
           integer                :: ps_import
           integer                :: u_import
           integer                :: v_import
           integer                :: temp_import
           integer                :: q_import
           integer                :: oz_import
           integer                :: cld_import
           integer                :: p_import
           integer                :: dp_import
           integer                :: dpdt_import
                                                                                
           integer                :: idate1_export
           integer                :: z_export
           integer                :: ps_export
           integer                :: u_export
           integer                :: v_export
           integer                :: temp_export
           integer                :: q_export
           integer                :: oz_export
           integer                :: cld_export
           integer                :: p_export
           integer                :: dp_export
           integer                :: dpdt_export


!For the surface file to couple with others or io
!----------------------------------------------------
           INTEGER                :: orography_import
           INTEGER                :: t_skin_import
           INTEGER                :: soil_mois_import
           INTEGER                :: snow_depth_import
           INTEGER                :: soil_t_import
           INTEGER                :: deep_soil_t_import
           INTEGER                :: roughness_import
           INTEGER                :: conv_cloud_cover_import
           INTEGER                :: conv_cloud_base_import
           INTEGER                :: conv_cloud_top_import
           INTEGER                :: albedo_visible_scattered_import
           INTEGER                :: albedo_visible_beam_import
           INTEGER                :: albedo_nearIR_scattered_import
           INTEGER                :: albedo_nearIR_beam_import
           INTEGER                :: sea_level_ice_mask_import
           INTEGER                :: vegetation_cover_import
           INTEGER                :: canopy_water_import
           INTEGER                :: m10_wind_fraction_import
           INTEGER                :: vegetation_type_import
           INTEGER                :: soil_type_import
           INTEGER                :: zeneith_angle_facsf_import
           INTEGER                :: zeneith_angle_facwf_import
           INTEGER                :: uustar_import
           INTEGER                :: ffmm_import
           INTEGER                :: ffhh_import
           INTEGER                :: sea_ice_thickness_import
           INTEGER                :: sea_ice_concentration_import
           INTEGER                :: tprcp_import
           INTEGER                :: srflag_import
           INTEGER                :: actual_snow_depth_import
           INTEGER                :: liquid_soil_moisture_import
           INTEGER                :: vegetation_cover_min_import
           INTEGER                :: vegetation_cover_max_import
           INTEGER                :: slope_type_import
           INTEGER                :: snow_albedo_max_import

           INTEGER                :: orography_export
           INTEGER                :: t_skin_export
           INTEGER                :: soil_mois_export
           INTEGER                :: snow_depth_export
           INTEGER                :: soil_t_export
           INTEGER                :: deep_soil_t_export
           INTEGER                :: roughness_export
           INTEGER                :: conv_cloud_cover_export
           INTEGER                :: conv_cloud_base_export
           INTEGER                :: conv_cloud_top_export
           INTEGER                :: albedo_visible_scattered_export
           INTEGER                :: albedo_visible_beam_export
           INTEGER                :: albedo_nearIR_scattered_export
           INTEGER                :: albedo_nearIR_beam_export
           INTEGER                :: sea_level_ice_mask_export
           INTEGER                :: vegetation_cover_export
           INTEGER                :: canopy_water_export
           INTEGER                :: m10_wind_fraction_export
           INTEGER                :: vegetation_type_export
           INTEGER                :: soil_type_export
           INTEGER                :: zeneith_angle_facsf_export
           INTEGER                :: zeneith_angle_facwf_export
           INTEGER                :: uustar_export
           INTEGER                :: ffmm_export
           INTEGER                :: ffhh_export
           INTEGER                :: sea_ice_thickness_export
           INTEGER                :: sea_ice_concentration_export
           INTEGER                :: tprcp_export
           INTEGER                :: srflag_export
           INTEGER                :: actual_snow_depth_export
           INTEGER                :: liquid_soil_moisture_export
           INTEGER                :: vegetation_cover_min_export
           INTEGER                :: vegetation_cover_max_export
           INTEGER                :: slope_type_export
           INTEGER                :: snow_albedo_max_export
      end type gfs_phy_state_namelist

      end module gfs_physics_namelist_mod
