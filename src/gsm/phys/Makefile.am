noinst_LIBRARIES = libgfsphysf.a libgfsphys.a

libgfsphys_a_FFLAGS = $(FC_MODINC) ../libutil $(FC_MODINC) ../dyn $(FC_MODINC) ../../phys $(FC_MODINC) ../../share
libgfsphys_a_FFLAGS += $(FC_AUTODOUBLE) $(FCFLAGS_FIXED)

libgfsphys_a_SOURCES = aoicpl_prep.f compns_physics.f coordinate_def.f \
	coundummy.f d3d_def.f dfi_fixwr.f do_physics_one_step.f fix_fields.f fix_fields_idea_rst.f \
	fld_collect.f fld_collect_nst.f gbphys_adv_hyb_gc.f gbphys_adv_hyb_gc_h.f gcycle.f \
	getcon_physics.f get_lats_node_shuff_r_fix.f GFS_simple_scatter.f gg_def.f \
	glats_physics.f gloopb.f gloopr.f gribit_gsm.f hpmdummy.f info.f lats_shuff.f layout1.f \
	lon_lat_para.f mod_state.f mpi_def.f mpi_quit.f namelist_physics_def.f \
	para_fixio_w.f para_nst_w.f read_fix.f reduce_lons_grid_module.f \
	resol_def.f setlats.f setlats_r_slg.f sortrx.f vert_def.f wrt3d.f \
	wrt3d_hyb.f wrtout_physics.f function_indlsev function_indlsod

libgfsphysf_a_FFLAGS = -I $(top_srcdir)/include \
	$(FC_MODINC) ../libutil $(FC_MODINC) ../dyn $(FC_MODINC) ../../phys $(FC_MODINC) ../../share
libgfsphysf_a_FFLAGS += $(FC_AUTODOUBLE) $(FCFLAGS_FREE)

libgfsphysf_a_SOURCES = checkiauforcing.f gfs_physics_add_get_state_ESMFField_mod.f \
	gfs_physics_aoi_var_mod.f gfs_physics_err_msg_mod.f gfs_physics_finalize_mod.f \
	gfs_physics_g2d_mod.f gfs_physics_g3d_mod.f gfs_physics_getcf_mod.f gfs_physics_grid_comp_mod.f \
	gfs_physics_grid_create_mod.f gfs_physics_gridgr_mod.f gfs_physics_initialize_mod.f \
	gfs_physics_internal_state_mod.f gfs_physics_namelist_mod.f gfs_physics_nst_var_mod.f \
	gfs_physics_output.f gfs_physics_run_mod.f gfs_physics_sfc_flx_mod.f gfs_physics_sfc_flx_set_mod.f \
	gfs_physics_start_time_get_mod.f gfs_phy_states_mod.f nstio_module.f sfcio_module.f

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
#
# - gfsphys
#
libgfsphys_a-aoicpl_prep.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) \
	libgfsphys_a-resol_def.$(OBJEXT) libgfsphysf_a-gfs_physics_aoi_var_mod.$(OBJEXT) \
	libgfsphys_a-mpi_def.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-compns_physics.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) libgfsphys_a-mpi_def.$(OBJEXT)
libgfsphys_a-dfi_fixwr.$(OBJEXT) : libgfsphysf_a-gfs_physics_sfc_flx_mod.$(OBJEXT) \
	libgfsphys_a-namelist_physics_def.$(OBJEXT) libgfsphys_a-resol_def.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_nst_var_mod.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-do_physics_one_step.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) \
	libgfsphys_a-resol_def.$(OBJEXT) libgfsphysf_a-gfs_physics_g2d_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_sfc_flx_set_mod.$(OBJEXT) libgfsphys_a-mpi_def.$(OBJEXT) \
	libgfsphys_a-d3d_def.$(OBJEXT) libgfsphysf_a-gfs_physics_sfc_flx_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_gridgr_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_aoi_var_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_nst_var_mod.$(OBJEXT) libgfsphys_a-vert_def.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_g3d_mod.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-fix_fields.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) \
	libgfsphys_a-resol_def.$(OBJEXT) libgfsphysf_a-gfs_physics_sfc_flx_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_nst_var_mod.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT) libgfsphys_a-gg_def.$(OBJEXT)
libgfsphys_a-fix_fields_idea_rst.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) \
	libgfsphys_a-resol_def.$(OBJEXT) libgfsphysf_a-gfs_physics_sfc_flx_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_nst_var_mod.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT) libgfsphys_a-gg_def.$(OBJEXT)
libgfsphys_a-fld_collect.$(OBJEXT) : libgfsphys_a-resol_def.$(OBJEXT) \
	libgfsphys_a-mpi_def.$(OBJEXT) libgfsphysf_a-gfs_physics_sfc_flx_mod.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-fld_collect_nst.$(OBJEXT) : libgfsphys_a-resol_def.$(OBJEXT) \
	libgfsphys_a-mpi_def.$(OBJEXT) libgfsphysf_a-gfs_physics_nst_var_mod.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-gbphys_adv_hyb_gc.$(OBJEXT) : libgfsphys_a-coordinate_def.$(OBJEXT) libgfsphys_a-resol_def.$(OBJEXT)
libgfsphys_a-gbphys_adv_hyb_gc_h.$(OBJEXT) : libgfsphys_a-coordinate_def.$(OBJEXT) libgfsphys_a-resol_def.$(OBJEXT)
libgfsphys_a-gcycle.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) \
	libgfsphys_a-resol_def.$(OBJEXT) libgfsphysf_a-gfs_physics_sfc_flx_mod.$(OBJEXT)
libgfsphys_a-getcon_physics.$(OBJEXT) : function_indlsev function_indlsod libgfsphys_a-namelist_physics_def.$(OBJEXT) \
	libgfsphys_a-resol_def.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT) libgfsphys_a-gg_def.$(OBJEXT)
libgfsphys_a-get_lats_node_shuff_r_fix.$(OBJEXT) : libgfsphys_a-resol_def.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-GFS_simple_scatter.$(OBJEXT) : libgfsphys_a-resol_def.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-gloopb.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) libgfsphysf_a-gfs_physics_gridgr_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_nst_var_mod.$(OBJEXT) libgfsphys_a-coordinate_def.$(OBJEXT) \
	libgfsphys_a-resol_def.$(OBJEXT) libgfsphys_a-mpi_def.$(OBJEXT) libgfsphys_a-d3d_def.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_sfc_flx_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_aoi_var_mod.$(OBJEXT) \
	libgfsphys_a-vert_def.$(OBJEXT) libgfsphysf_a-gfs_physics_g3d_mod.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT) libgfsphys_a-gg_def.$(OBJEXT)
libgfsphys_a-gloopr.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_gridgr_mod.$(OBJEXT) libgfsphys_a-resol_def.$(OBJEXT) \
	libgfsphys_a-d3d_def.$(OBJEXT) libgfsphysf_a-gfs_physics_aoi_var_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_g3d_mod.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT) libgfsphys_a-gg_def.$(OBJEXT)
libgfsphys_a-gribit_gsm.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT)
libgfsphys_a-info.$(OBJEXT) : libgfsphys_a-resol_def.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-lats_shuff.$(OBJEXT) : libgfsphys_a-resol_def.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-lon_lat_para.$(OBJEXT) : libgfsphys_a-resol_def.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT) libgfsphys_a-gg_def.$(OBJEXT)
libgfsphys_a-mpi_quit.$(OBJEXT) : libgfsphys_a-mpi_def.$(OBJEXT)
libgfsphys_a-para_fixio_w.$(OBJEXT) : libgfsphys_a-resol_def.$(OBJEXT) libgfsphysf_a-gfs_physics_sfc_flx_mod.$(OBJEXT) \
	libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-para_nst_w.$(OBJEXT) : libgfsphysf_a-gfs_physics_nst_var_mod.$(OBJEXT) libgfsphys_a-resol_def.$(OBJEXT) \
	libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-read_fix.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) libgfsphys_a-mod_state.$(OBJEXT) \
	libgfsphysf_a-sfcio_module.$(OBJEXT) libgfsphysf_a-gfs_physics_nst_var_mod.$(OBJEXT) libgfsphysf_a-nstio_module.$(OBJEXT) \
	libgfsphys_a-coordinate_def.$(OBJEXT) libgfsphys_a-resol_def.$(OBJEXT) libgfsphys_a-mpi_def.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_sfc_flx_mod.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-setlats.$(OBJEXT) : libgfsphys_a-resol_def.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-setlats_r_slg.$(OBJEXT) : libgfsphys_a-resol_def.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-wrt3d.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) libgfsphys_a-resol_def.$(OBJEXT) \
	libgfsphys_a-d3d_def.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-wrt3d_hyb.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) libgfsphys_a-resol_def.$(OBJEXT) \
	libgfsphys_a-d3d_def.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphys_a-wrtout_physics.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) libgfsphys_a-mod_state.$(OBJEXT) \
	libgfsphysf_a-sfcio_module.$(OBJEXT) libgfsphysf_a-nstio_module.$(OBJEXT) libgfsphysf_a-gfs_physics_nst_var_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_output.$(OBJEXT) libgfsphys_a-resol_def.$(OBJEXT) libgfsphysf_a-gfs_physics_g2d_mod.$(OBJEXT) \
	libgfsphys_a-mpi_def.$(OBJEXT) libgfsphysf_a-gfs_physics_sfc_flx_mod.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
#
# - gfsphysf
#
libgfsphysf_a-checkiauforcing.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT)
libgfsphysf_a-gfs_physics_finalize_mod.$(OBJEXT) : libgfsphysf_a-gfs_physics_internal_state_mod.$(OBJEXT)
libgfsphysf_a-gfs_physics_getcf_mod.$(OBJEXT) : libgfsphysf_a-gfs_physics_err_msg_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_internal_state_mod.$(OBJEXT)
libgfsphysf_a-gfs_physics_grid_comp_mod.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) \
	libgfsphys_a-resol_def.$(OBJEXT) libgfsphys_a-mpi_def.$(OBJEXT) libgfsphysf_a-gfs_physics_finalize_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_initialize_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_err_msg_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_phy_states_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_internal_state_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_getcf_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_run_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_grid_create_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_output.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphysf_a-gfs_physics_grid_create_mod.$(OBJEXT) : libgfsphysf_a-gfs_physics_internal_state_mod.$(OBJEXT)
libgfsphysf_a-gfs_physics_initialize_mod.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_sfc_flx_set_mod.$(OBJEXT) libgfsphys_a-mod_state.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_gridgr_mod.$(OBJEXT) libgfsphysf_a-sfcio_module.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_nst_var_mod.$(OBJEXT) libgfsphys_a-coordinate_def.$(OBJEXT) \
	libgfsphys_a-resol_def.$(OBJEXT) libgfsphysf_a-gfs_physics_g2d_mod.$(OBJEXT) \
	libgfsphys_a-mpi_def.$(OBJEXT) libgfsphysf_a-gfs_physics_internal_state_mod.$(OBJEXT) \
	libgfsphys_a-d3d_def.$(OBJEXT) libgfsphysf_a-gfs_physics_aoi_var_mod.$(OBJEXT) \
	libgfsphys_a-vert_def.$(OBJEXT) libgfsphysf_a-gfs_physics_g3d_mod.$(OBJEXT) \
	libgfsphys_a-layout1.$(OBJEXT) libgfsphys_a-gg_def.$(OBJEXT)
libgfsphysf_a-gfs_physics_internal_state_mod.$(OBJEXT) : libgfsphysf_a-gfs_physics_gridgr_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_nst_var_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_g2d_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_sfc_flx_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_namelist_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_aoi_var_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_g3d_mod.$(OBJEXT)
libgfsphysf_a-gfs_physics_output.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) \
	libgfsphys_a-mod_state.$(OBJEXT) libgfsphys_a-coordinate_def.$(OBJEXT) libgfsphys_a-resol_def.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_err_msg_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_internal_state_mod.$(OBJEXT)
libgfsphysf_a-gfs_physics_run_mod.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) \
	libgfsphys_a-resol_def.$(OBJEXT) libgfsphysf_a-gfs_physics_internal_state_mod.$(OBJEXT) libgfsphys_a-layout1.$(OBJEXT)
libgfsphysf_a-gfs_physics_sfc_flx_set_mod.$(OBJEXT) : libgfsphysf_a-gfs_physics_sfc_flx_mod.$(OBJEXT)
libgfsphysf_a-gfs_physics_start_time_get_mod.$(OBJEXT) : libgfsphysf_a-sfcio_module.$(OBJEXT)
libgfsphysf_a-gfs_phy_states_mod.$(OBJEXT) : libgfsphys_a-namelist_physics_def.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_grid_create_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_g2d_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_err_msg_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_internal_state_mod.$(OBJEXT) \
	libgfsphysf_a-gfs_physics_namelist_mod.$(OBJEXT) libgfsphysf_a-gfs_physics_add_get_state_ESMFField_mod.$(OBJEXT) \
	libgfsphys_a-gg_def.$(OBJEXT)
