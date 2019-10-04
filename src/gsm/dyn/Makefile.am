noinst_LIBRARIES = libdynio.a libdynf.a libdyn.a libdyns.a

libdyn_a_SOURCES = bfilt_def.f common_to_model_vars.f compns_dynamics.f coundummy.f \
	damp_speed.f damp_speed_noq.f deldifs_def.f deldifs.f deldifs_noq.f \
	deldifs_tracers.f delnpe.f delnpo.f dezouv.f do_dynamics_mod.f \
	do_dynamics_one_loop.f do_dynamics_slg_loop.f do_dynamics_two_loop.f \
	dozeuv.f dyn_gocart_tracer_config_stub.f epslon.f excha.f filter1eo.f \
	filter1eo_noq.f filter2eo.f filter2eo_noq.f filtereo.f filtereo_noq.f four2fln.f \
	four2fln_gg.f four2grid_thread.f four2grid_thread.fftpack.f four_to_grid.f \
	gather_times.f gch2press.f get_am_bm_hyb.f get_am_bm_hyb_gcdp.f \
	get_am_bm_hyb_gc.f get_cd_hyb.f get_cd_hyb_gcdp.f get_cd_hyb_gc.f \
	get_cd_hyb_slg.f getcon_dynamics.f getcon_lag.f get_lats_node_shuff_a_fix.f \
	get_ls_node.f getpwatptot.f get_topo_grid_grad.f get_w_z.f \
	getysk.f getysk_gc.f gfidi_gchdp_noadv_noq.f gfidi_hyb.f \
	gfidi_hyb_gc.f gfidi_hyb_gchdp.f gfidi_hyb_gc_h.f gfidi_hyb_gc_h_w_z.f \
	gfidi_sig.f gfidi_speed.f gfidi_tracers.f \
	gfs_dficoll_dynamics.f gfs_dficoll_dynamics_slg.f \
	gfs_dyn_bfilt_def.f gfs_dyn_coordinate_def.f gfs_dyn_date_def.f \
	gfs_dyn_deldifs_def.f gfs_dyn_dfi_mod.f gfs_dyn_gg_def.f gfs_dyn_glats.f \
	gfs_dyn_io_header.f gfs_dyn_layout1.f gfs_dyn_lon_lat_para.f gfs_dyn_machine.f \
	gfs_dyn_matrix_sig_def.f gfs_dyn_mpi_def.f \
	gfs_dyn_reduce_lons_grid_module.f gfs_dyn_resol_def.f gfs_dyn_tracer_config.f \
	gfs_dyn_tracer_const.f gfs_dyn_vert_def.f gloopa_tracers.f \
	gozrineo.f grid_collect.f grid_collect_ipe.f grid_fix.f grid_to_spect.f grid_to_spect_inp.f \
	grid_to_spect_inp_slg.f grid_to_spect_rqt.f grid_to_spect_slg.f gridzz_to_spect.f hpmdummy.f hyb2press.f \
	idea_deldifs.f idea_deldifs_init.f idea_deldifs_noq.f idea_getcoef.f ifshuff.f impadj_hyb.locl.f \
	impadj_hyb.locl_gcdp.f impadj_hyb.locl_gc.f impadj_sig.locl.f impadj_slg.locl.f info.f input_fields.f \
	input_fields_rst.f input_fields_slg.f input_for_wam_ipe_rst.f lats_shuff.f layout_grid_tracers.f \
	layout_lag.f matinv.f model_to_common_vars.f mpi_GV.f mpi_quit.f namelist_dynamics_def.f \
	ndslfv_monoadvh2.f ndslfv_monoadvh.f ndslfv_monoadvv.f ndslfv_pack.f ndslfv_para.f nemsio_def.f \
	num_parthds.f omega_gcdp.f omega_gch.f omega_hyb.f omega_sig.f pln2eo.f pmgrid.f \
	redist_lats_a.f rms_spect.f set_halos_lag_sn.f setlats.f setlats_lag.f setlats_slg.f setsig.f \
	shuffle_grid.f sicdif_hyb.f sicdif_hyb_gcdp.f sicdif_hyb_gc.f sicdif_hyb_slg.f sicdif_sig.f \
	sig2press.f sig_io.f slgscan_all_redgg.f slgshr.f sltini_all.f \
	sortrx.f spectpz_to_gridxy.f spect_to_grid.f spect_to_grid_gz.f spect_to_grid_iau.f \
	spect_to_grid_inp.f spect_to_grid_rqt.f spect_to_grid_slg.f spect_to_gridxy.f sumder2.f \
	sumder2_slg_gg.f sumfln.f sumfln_slg_gg.f treadeo.io.f treadeo.io_iau.f treadeo.io_slg.f \
	treadeo_nemsio.f treadeo_nemsio_iau.f treadeo_nemsio_slg.f treadg_nemsio.f treads_nemsio.f \
	triseof.f triseori.f twriteg_rst.f twrites_hst.f twrites_rst.f twrites_rst_idea.f uveodz.f \
	uvoedz.f ver_sig.f wrtout_dynamics.f wrtout_dynamics_slg_nemsio.f ysminv.f \
	function2 function_indlsev function_indlsod

libdyn_a_FFLAGS = $(FC_MODINC) ../libutil $(FC_MODINC) ../../phys $(FC_MODINC) ../../share
libdyn_a_FFLAGS += $(FC_AUTODOUBLE) $(FCFLAGS_FIXED)


libdynf_a_SOURCES = GFS_AddParameterToStateMod.f gfs_dynamics_add_get_state_ESMFField_mod.f \
	gfs_dynamics_err_msg_mod.f gfs_dynamics_finalize_mod.f gfs_dynamics_getcf_mod.f     \
	gfs_dynamics_grid_comp_mod.f gfs_dynamics_grid_create_mod.f gfs_dynamics_initialize_mod.f \
	gfs_dynamics_initialize_slg_mod.f gfs_dynamics_internal_state_mod.f gfs_dynamics_namelist_mod.f \
	gfs_dynamics_output.f gfs_dynamics_run_mod.f gfs_dynamics_start_time_get_mod.f gfs_dynamics_states_mod.f \
	gfs_dyn_iau_mod.f gfs_dyn_physcons.f

libdynf_a_FFLAGS = $(FC_MODINC) ../libutil $(FC_MODINC) ../../phys $(FC_MODINC) ../../share
libdynf_a_FFLAGS += $(FC_AUTODOUBLE) $(FCFLAGS_FREE)


libdyns_a_SOURCES = gfs_dyn_stoch_data_mod.f gfs_get_pattern_mod.f gfs_dyn_patterngenerator.f

libdyns_a_FFLAGS = $(FC_MODINC) ../libutil $(FC_MODINC) ../../phys $(FC_MODINC) ../../share
libdyns_a_FFLAGS += $(FC_AUTODOUBLE) $(FCFLAGS_FREE)


libdynio_a_SOURCES = sigio_module.f sigio_r_module.f

libdynio_a_FFLAGS = $(FC_MODINC) ../libutil $(FC_MODINC) ../../phys $(FC_MODINC) ../../share
libdynio_a_FFLAGS += $(FCFLAGS_FREE)

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
# - dynio
libdynio_a-sigio_r_module.$(OBJEXT) : libdynio_a-sigio_module.$(OBJEXT) 

# - dyn
libdyn_a-bfilt_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT)
libdyn_a-gfs_dyn_bfilt_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-common_to_model_vars.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) 
libdyn_a-compns_dynamics.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-layout_grid_tracers.$(OBJEXT) libdyn_a-pmgrid.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) 
libdyn_a-damp_speed.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-damp_speed_noq.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-deldifs_def.$(OBJEXT) : function_indlsev function_indlsod libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_deldifs_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-deldifs.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-deldifs_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT)
libdyn_a-idea_deldifs.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_deldifs_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-deldifs_noq.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_deldifs_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-idea_deldifs_noq.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_deldifs_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-deldifs_tracers.$(OBJEXT) : function_indlsev function_indlsod libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) \
	libdyn_a-deldifs_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-delnpe.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-delnpo.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-dezouv.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-do_dynamics_mod.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) 
libdyn_a-do_dynamics_one_loop.$(OBJEXT) : function2 libdyn_a-do_dynamics_mod.$(OBJEXT) libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) \
	libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_dfi_mod.$(OBJEXT) \
	libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-do_dynamics_slg_loop.$(OBJEXT) : function_indlsev libdyn_a-do_dynamics_mod.$(OBJEXT) libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) \
	libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-layout_grid_tracers.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_bfilt_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_dfi_mod.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-do_dynamics_two_loop.$(OBJEXT) : function2 libdyn_a-do_dynamics_mod.$(OBJEXT) libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) \
	libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_bfilt_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_dfi_mod.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-dozeuv.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-epslon.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-excha.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-filter1eo.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-filter1eo_noq.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-filter2eo.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-filter2eo_noq.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-filtereo.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-filtereo_noq.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-four2fln.$(OBJEXT) : function2 libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-four2fln_gg.$(OBJEXT) : function_indlsev function_indlsod libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-four2grid_thread.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-four2grid_thread.fftpack.$(OBJEXT) : libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-four_to_grid.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gather_times.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-gch2press.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) \
	libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-get_am_bm_hyb.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-get_am_bm_hyb_gcdp.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-get_am_bm_hyb_gc.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-get_cd_hyb.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-get_cd_hyb_gcdp.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-get_cd_hyb_gc.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-get_cd_hyb_slg.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-getcon_dynamics.$(OBJEXT) : function2 libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_bfilt_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-getcon_lag.$(OBJEXT) : libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-layout_lag.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) 
libdyn_a-get_lats_node_shuff_a_fix.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-get_ls_node.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-getpwatptot.$(OBJEXT) : libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-get_topo_grid_grad.$(OBJEXT) : function_indlsev function_indlsod libdyn_a-gfs_dyn_io_header.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-get_w_z.$(OBJEXT) : libdyn_a-do_dynamics_mod.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) libdyn_a-gfs_dyn_dfi_mod.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-getysk.$(OBJEXT) : libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-getysk_gc.$(OBJEXT) : libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-gfidi_gchdp_noadv_noq.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) 
libdyn_a-gfidi_hyb.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfidi_hyb_gc.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfidi_hyb_gchdp.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) 
libdyn_a-gfidi_hyb_gc_h.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) 
libdyn_a-gfidi_hyb_gc_h_w_z.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) 
libdyn_a-gfidi_sig.$(OBJEXT) : libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfidi_speed.$(OBJEXT) : libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfidi_tracers.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dficoll_dynamics.$(OBJEXT) : libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_dfi_mod.$(OBJEXT) 
libdyn_a-gfs_dficoll_dynamics_slg.$(OBJEXT) : libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-layout_grid_tracers.$(OBJEXT) libdyn_a-layout_lag.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_dfi_mod.$(OBJEXT) 
libdyn_a-gfs_dyn_bfilt_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_date_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_deldifs_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_dfi_mod.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_gg_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_glats.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_io_header.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_lon_lat_para.$(OBJEXT) : libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_matrix_sig_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdynf_a-gfs_dyn_physcons.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_reduce_lons_grid_module.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_resol_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) 
libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) : libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gfs_dyn_vert_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-gloopa_tracers.$(OBJEXT) : function_indlsev libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-layout_grid_tracers.$(OBJEXT) \
	libdyn_a-layout_lag.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-slgshr.$(OBJEXT) libdyn_a-pmgrid.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_dfi_mod.$(OBJEXT) 
libdyn_a-gozrineo.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-grid_collect.$(OBJEXT) : libdyn_a-gfs_dyn_io_header.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-wrtout_dynamics.$(OBJEXT) libdynf_a-gfs_dynamics_output.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-grid_collect_ipe.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-wrtout_dynamics.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) 
libdyn_a-grid_fix.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-grid_to_spect.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) \
	libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-grid_to_spect_inp.$(OBJEXT) : function2 libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-grid_to_spect_inp_slg.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-layout_grid_tracers.$(OBJEXT) \
	libdyn_a-layout_lag.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-grid_to_spect_rqt.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) \
	libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-grid_to_spect_slg.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-layout_grid_tracers.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-gridzz_to_spect.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-hyb2press.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-idea_deldifs.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_deldifs_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-idea_deldifs_init.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_deldifs_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) 
libdyn_a-idea_deldifs_noq.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_deldifs_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-idea_getcoef.$(OBJEXT) : libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) 
libdyn_a-ifshuff.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-impadj_hyb.locl.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-impadj_hyb.locl_gcdp.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-impadj_hyb.locl_gc.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-impadj_sig.locl.$(OBJEXT) : libdyn_a-gfs_dyn_matrix_sig_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-impadj_slg.locl.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-info.$(OBJEXT) : libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-input_fields.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-input_fields_rst.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-input_fields_slg.$(OBJEXT) : libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-layout_grid_tracers.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-input_for_wam_ipe_rst.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-lats_shuff.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-layout_grid_tracers.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-layout_lag.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-matinv.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-model_to_common_vars.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) \
	libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-mpi_quit.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) 

libdyn_a-namelist_dynamics_def.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-ndslfv_monoadvh2.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-ndslfv_monoadvh.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-ndslfv_monoadvv.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-ndslfv_pack.$(OBJEXT) : libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) 
libdyn_a-ndslfv_para.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) 
libdyn_a-omega_gcdp.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-omega_gch.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-omega_hyb.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-omega_sig.$(OBJEXT) : libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-pln2eo.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-redist_lats_a.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-rms_spect.$(OBJEXT) : function2 libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) 
libdyn_a-set_halos_lag_sn.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-layout_lag.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-setlats.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-setlats_lag.$(OBJEXT) : libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-setlats_slg.$(OBJEXT) : libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-setsig.$(OBJEXT) : libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-shuffle_grid.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-sicdif_hyb.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-sicdif_hyb_gcdp.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-sicdif_hyb_gc.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-sicdif_hyb_slg.$(OBJEXT) : function_indlsev function_indlsod libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-sicdif_sig.$(OBJEXT) : function2 libdyn_a-gfs_dyn_matrix_sig_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-sig2press.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-sig_io.$(OBJEXT) : libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-slgscan_all_redgg.$(OBJEXT) : libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-slgshr.$(OBJEXT) libdyn_a-pmgrid.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-sltini_all.$(OBJEXT) : libdyn_a-layout_lag.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-pmgrid.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-spectpz_to_gridxy.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-spect_to_grid.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-spect_to_grid_gz.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-spect_to_grid_iau.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-spect_to_grid_inp.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-spect_to_grid_rqt.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-spect_to_grid_slg.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-layout_grid_tracers.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-spect_to_gridxy.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-sumder2.$(OBJEXT) : function2 libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-sumder2_slg_gg.$(OBJEXT) : function_indlsev function_indlsod libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-sumfln.$(OBJEXT) : function2 libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-sumfln_slg_gg.$(OBJEXT) : function_indlsev function_indlsod libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-treadeo.io.$(OBJEXT) : function2 libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_io_header.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdynio_a-sigio_module.$(OBJEXT) libdynio_a-sigio_r_module.$(OBJEXT) \
	libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-treadeo.io_iau.$(OBJEXT) : function2 libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_io_header.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdynf_a-gfs_dyn_iau_mod.$(OBJEXT) \
	libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdynio_a-sigio_module.$(OBJEXT) \
	libdynio_a-sigio_r_module.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-treadeo.io_slg.$(OBJEXT) : function2 libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_io_header.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) libdynio_a-sigio_module.$(OBJEXT) \
	libdynio_a-sigio_r_module.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-treadeo_nemsio.$(OBJEXT) : function2 libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_io_header.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-nemsio_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-treadeo_nemsio_iau.$(OBJEXT) : function2 libdyn_a-gfs_dyn_io_header.$(OBJEXT) libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) \
	libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-nemsio_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-treadeo_nemsio_slg.$(OBJEXT) : function2 libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_io_header.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-nemsio_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-treadg_nemsio.$(OBJEXT) : libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_io_header.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-treads_nemsio.$(OBJEXT) : libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_io_header.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyn_a-triseof.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-triseori.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-twriteg_rst.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_io_header.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) 
libdyn_a-twrites_hst.$(OBJEXT) : function2 libdyn_a-gfs_dyn_io_header.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdynio_a-sigio_module.$(OBJEXT) libdynio_a-sigio_r_module.$(OBJEXT) \
	libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) 
libdyn_a-twrites_rst.$(OBJEXT) : function2 libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) 
libdyn_a-twrites_rst_idea.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) 
libdyn_a-uveodz.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-uvoedz.$(OBJEXT) : function2 libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) 
libdyn_a-ver_sig.$(OBJEXT) : libdyn_a-gfs_dyn_matrix_sig_def.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdyn_a-wrtout_dynamics.$(OBJEXT) : libdyn_a-do_dynamics_mod.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) \
	libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) 
libdyn_a-wrtout_dynamics_slg_nemsio.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-layout_grid_tracers.$(OBJEXT) \
	libdyn_a-layout_lag.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) 
libdyn_a-ysminv.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 
# - dynf
libdynf_a-GFS_AddParameterToStateMod.$(OBJEXT) : libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdynf_a-gfs_dynamics_internal_state_mod.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT)
libdynf_a-gfs_dynamics_finalize_mod.$(OBJEXT) : libdynf_a-gfs_dynamics_internal_state_mod.$(OBJEXT) 
libdynf_a-gfs_dynamics_getcf_mod.$(OBJEXT) : libdynf_a-gfs_dynamics_internal_state_mod.$(OBJEXT) libdynf_a-gfs_dynamics_err_msg_mod.$(OBJEXT) 
libdynf_a-gfs_dynamics_grid_comp_mod.$(OBJEXT) : libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdynf_a-gfs_dynamics_finalize_mod.$(OBJEXT) libdynf_a-gfs_dynamics_initialize_slg_mod.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdynf_a-gfs_dynamics_initialize_mod.$(OBJEXT) libdynf_a-gfs_dynamics_states_mod.$(OBJEXT) libdynf_a-gfs_dyn_iau_mod.$(OBJEXT) \
	libdynf_a-gfs_dynamics_err_msg_mod.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) \
	libdynf_a-GFS_AddParameterToStateMod.$(OBJEXT) libdynf_a-gfs_dynamics_grid_create_mod.$(OBJEXT) libdynf_a-gfs_dynamics_output.$(OBJEXT) \
	libdynf_a-gfs_dynamics_run_mod.$(OBJEXT)
libdynf_a-gfs_dynamics_grid_create_mod.$(OBJEXT) : libdynf_a-gfs_dynamics_internal_state_mod.$(OBJEXT)
libdynf_a-gfs_dynamics_initialize_mod.$(OBJEXT) : libdynf_a-gfs_dynamics_getcf_mod.$(OBJEXT) libdyn_a-gfs_dyn_io_header.$(OBJEXT) \
	libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdynf_a-gfs_dyn_iau_mod.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-wrtout_dynamics.$(OBJEXT) \
	libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyns_a-gfs_dyn_stoch_data_mod.$(OBJEXT) 
libdynf_a-gfs_dynamics_initialize_slg_mod.$(OBJEXT) : libdynf_a-gfs_dynamics_getcf_mod.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_io_header.$(OBJEXT) libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-layout_grid_tracers.$(OBJEXT) libdyn_a-layout_lag.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdynf_a-gfs_dyn_iau_mod.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-wrtout_dynamics.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) 
libdynf_a-gfs_dynamics_internal_state_mod.$(OBJEXT) : libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dynamics_namelist_mod.$(OBJEXT) \
	libdyn_a-gfs_dyn_gg_def.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_dfi_mod.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) \
	libdyn_a-gfs_dyn_vert_def.$(OBJEXT)
libdynf_a-gfs_dynamics_output.$(OBJEXT) : libdyn_a-gfs_dyn_io_header.$(OBJEXT) libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdynf_a-gfs_dynamics_internal_state_mod.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdynf_a-gfs_dynamics_err_msg_mod.$(OBJEXT) libdyn_a-wrtout_dynamics.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) \
	libdyn_a-gfs_dyn_vert_def.$(OBJEXT)
libdynf_a-gfs_dynamics_run_mod.$(OBJEXT) : libdyns_a-gfs_get_pattern_mod.$(OBJEXT) libdyn_a-do_dynamics_mod.$(OBJEXT) \
	libdynf_a-gfs_dynamics_internal_state_mod.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT)
libdynf_a-gfs_dynamics_start_time_get_mod.$(OBJEXT) : libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) \
	libdynio_a-sigio_module.$(OBJEXT) libdynio_a-sigio_r_module.$(OBJEXT)
libdynf_a-gfs_dynamics_states_mod.$(OBJEXT) : libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) libdynf_a-gfs_dynamics_internal_state_mod.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dynamics_namelist_mod.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) \
	libdyn_a-layout_grid_tracers.$(OBJEXT) libdynf_a-gfs_dynamics_err_msg_mod.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) \
	libdynf_a-gfs_dynamics_grid_create_mod.$(OBJEXT) libdynf_a-gfs_dynamics_add_get_state_ESMFField_mod.$(OBJEXT)
libdynf_a-gfs_dyn_iau_mod.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_tracer_config.$(OBJEXT) \
	libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) libdynf_a-gfs_dynamics_namelist_mod.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) \
	libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdynf_a-gfs_dynamics_err_msg_mod.$(OBJEXT) libdynf_a-gfs_dynamics_states_mod.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) libdynf_a-gfs_dynamics_grid_create_mod.$(OBJEXT) \
	libdynf_a-gfs_dynamics_add_get_state_ESMFField_mod.$(OBJEXT) 
libdynf_a-gfs_dyn_physcons.$(OBJEXT) : libdyn_a-gfs_dyn_machine.$(OBJEXT) 

# - dyns
libdyns_a-gfs_dyn_stoch_data_mod.$(OBJEXT) : libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_machine.$(OBJEXT) libdyns_a-gfs_dyn_patterngenerator.$(OBJEXT) 
libdyns_a-gfs_get_pattern_mod.$(OBJEXT) : function_indlsev function_indlsod libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_coordinate_def.$(OBJEXT) libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdyn_a-gfs_dyn_gg_def.$(OBJEXT) \
	libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-gfs_dyn_resol_def.$(OBJEXT) libdyn_a-gfs_dyn_machine.$(OBJEXT) \
	libdyn_a-namelist_dynamics_def.$(OBJEXT) libdyn_a-gfs_dyn_date_def.$(OBJEXT) libdyn_a-gfs_dyn_tracer_const.$(OBJEXT) \
	libdyns_a-gfs_dyn_patterngenerator.$(OBJEXT) libdyns_a-gfs_dyn_stoch_data_mod.$(OBJEXT) libdyn_a-gfs_dyn_vert_def.$(OBJEXT) 
libdyns_a-gfs_dyn_patterngenerator.$(OBJEXT) : function_indlsev function_indlsod libdyn_a-gfs_dyn_mpi_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_layout1.$(OBJEXT) libdynf_a-gfs_dyn_physcons.$(OBJEXT) libdyn_a-namelist_dynamics_def.$(OBJEXT) \
	libdyn_a-gfs_dyn_machine.$(OBJEXT) 
