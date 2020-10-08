noinst_LIBRARIES = libfgfs.a libgfs.a libgfsf.a libidea.a libgfsdrv.a

libgfs_a_SOURCES = calpreciptype.f cs_conv.f funcphys.f gcm_shoc.f machine.f \
	module_nst_model.f module_nst_parameters.f module_nst_water_prop.f   \
	physcons.f physparam.f

libgfs_a_FFLAGS = $(PHYS_DEFS) $(FCFLAGS_f) $(FC_AUTODOUBLE) $(FCFLAGS_FREE)

libgfsf_a_SOURCES = cnvc90.f date_def.f dcyc2.f dcyc2.pre.rad.f gbphys.f get_prs.f \
	gfs_phy_tracer_config.f grrad.f gscond.f gscondp.f gsmddrive.f gwdc.f \
	gwdps.f iounitdef.f lrgsclr.f mfdeepcnv.f mfpbl.f mfshalcnv.f \
	module_bfmicrophysics.f moninedmf.f moninp1.f moninp.f moninq.f moninshoc.f \
	mstadb.f mstadbtn.f mstcnv.f namelist_soilveg.f ozinterp.f ozne_def.f ozphys.f \
	precpd.f precpdp.f precpd_shoc.f progt2.f progtm_module.f radiation_aerosols.f \
	radiation_astronomy.f radiation_clouds.f radiation_gases.f radiation_surface.f \
	rad_initialize.f radlw_datatb.f radlw_main.f radlw_param.f radsw_datatb.f \
	radsw_main.f radsw_param.f rascnvv2.f rayleigh_damp.f rayleigh_damp_mesopause.f \
	sascnv.f sascnvn.f set_soilveg.f sfc_cice.f sfc_diag.f sfc_diff.f sfc_drv.f \
	sfc_land.f sfc_nst.f sfc_ocean.f sfc_sice.f sfcsub.f sflx.f shalcnv.f shalcv.f \
	shalcv_opr.f tracer_const_h.f tridi2t3.f

libgfsf_a_FFLAGS = $(FC_MODINC) ../share
libgfsf_a_FFLAGS += $(FC_AUTODOUBLE) $(FCFLAGS_FIXED)

libfgfs_a_SOURCES = num_parthds.F
libfgfs_a_FFLAGS  = $(FC_AUTODOUBLE) $(FCFLAGS_FREE)

libidea_a_SOURCES = co2hc.f efield.f h2oc.f h2ohdc.f ideaca.f idea_cal_advance.f idea_co2.f \
	idea_composition.f idea_dissipation.f idea_getno_snoe.f idea_h2o.f idea_imf_input.f \
	idea_interpol_datab.f idea_ion_empirmodels.f idea_ion.f idea_ion_input.f idea_io_units.f \
	idea_merge_ipe_to_wam.f idea_mpi_def.f idea_o2_o3.f idea_phys.f idea_solar_2014.f \
	idea_solar_heating.f idea_solar_init.f idea_solar_input.f idea_tracer.f idea_tracers_input.f \
	idea_trad_geopgrav.f namelist_wamphysics_def.f w05sc_efield_merge.f

libidea_a_FFLAGS = $(FC_MODINC) ../share
libidea_a_FFLAGS += $(FC_AUTODOUBLE) $(FC_LINE_LENGTH)

libgfsdrv_a_SOURCES = nuopc_physics.F90
libgfsdrv_a_FCFLAGS = $(FC_AUTODOUBLE)

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
# gfs
libgfs_a-calpreciptype.$(OBJEXT) : libgfs_a-funcphys.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfs_a-cs_conv.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfs_a-funcphys.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfs_a-gcm_shoc.$(OBJEXT) : libgfs_a-funcphys.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfs_a-module_nst_model.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-module_nst_parameters.$(OBJEXT) libgfs_a-module_nst_water_prop.$(OBJEXT) 
libgfs_a-module_nst_parameters.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfs_a-module_nst_water_prop.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-module_nst_parameters.$(OBJEXT) 
libgfs_a-physcons.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 
libgfs_a-physparam.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 

# gfsf
libgfsf_a-cnvc90.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 
libgfsf_a-date_def.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 
libgfsf_a-dcyc2.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-dcyc2.pre.rad.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-gbphys.$(OBJEXT) : libgfs_a-module_nst_water_prop.$(OBJEXT) libgfs_a-cs_conv.$(OBJEXT) \
	libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-get_prs.$(OBJEXT) : libgfsf_a-tracer_const_h.$(OBJEXT) libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-gfs_phy_tracer_config.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 
libgfsf_a-grrad.$(OBJEXT) : libgfsf_a-radlw_main.$(OBJEXT) libgfsf_a-radsw_param.$(OBJEXT) libgfs_a-physparam.$(OBJEXT) \
	libgfsf_a-radiation_clouds.$(OBJEXT) libgfsf_a-radiation_surface.$(OBJEXT) libgfsf_a-radlw_param.$(OBJEXT) \
	libgfsf_a-radsw_main.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) \
	libgfsf_a-radiation_astronomy.$(OBJEXT) libgfsf_a-radiation_gases.$(OBJEXT) libgfsf_a-radiation_aerosols.$(OBJEXT) 
libgfsf_a-gscond.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-gscondp.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-gsmddrive.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfsf_a-module_bfmicrophysics.$(OBJEXT) 
libgfsf_a-gwdc.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 
libgfsf_a-gwdps.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 
libgfsf_a-lrgsclr.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-mfdeepcnv.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-mfpbl.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-mfshalcnv.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-module_bfmicrophysics.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-moninedmf.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-moninp.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-moninp1.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-moninq.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-moninshoc.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-mstadb.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-mstadbtn.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-mstcnv.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libidea_a-namelist_wamphysics_def.$(OBJEXT) : libidea_a-idea_mpi_def.$(OBJEXT) 
libgfsf_a-ozinterp.$(OBJEXT) : libgfsf_a-ozne_def.$(OBJEXT) libgfs_a-machine.$(OBJEXT) 
libgfsf_a-ozne_def.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 
libgfsf_a-ozphys.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-precpd.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-precpdp.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-precpd_shoc.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-progt2.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfsf_a-progtm_module.$(OBJEXT) 
libgfsf_a-progtm_module.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 
libgfsf_a-radiation_aerosols.$(OBJEXT) : libgfsf_a-radsw_param.$(OBJEXT) libgfs_a-physparam.$(OBJEXT) \
	libgfsf_a-iounitdef.$(OBJEXT) libgfsf_a-radlw_param.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) \
	libgfsf_a-gfs_phy_tracer_config.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-radiation_astronomy.$(OBJEXT) : libgfs_a-physparam.$(OBJEXT) libgfsf_a-iounitdef.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-radiation_clouds.$(OBJEXT) : libgfs_a-physparam.$(OBJEXT) libgfsf_a-iounitdef.$(OBJEXT) \
	libgfs_a-physcons.$(OBJEXT) libgfsf_a-module_bfmicrophysics.$(OBJEXT) 
libgfsf_a-radiation_gases.$(OBJEXT) : libgfsf_a-ozne_def.$(OBJEXT) libgfs_a-physparam.$(OBJEXT) \
	libgfsf_a-iounitdef.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-radiation_surface.$(OBJEXT) : libgfs_a-physparam.$(OBJEXT) libgfsf_a-iounitdef.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-rad_initialize.$(OBJEXT) : libgfs_a-physparam.$(OBJEXT) libgfsf_a-grrad.$(OBJEXT) 
libgfsf_a-radlw_datatb.$(OBJEXT) : libgfs_a-physparam.$(OBJEXT) libgfsf_a-radlw_param.$(OBJEXT) 
libgfsf_a-radlw_main.$(OBJEXT) : libgfs_a-physparam.$(OBJEXT) libgfsf_a-radlw_param.$(OBJEXT) \
	libgfs_a-physcons.$(OBJEXT) libgfsf_a-radlw_datatb.$(OBJEXT) 
libgfsf_a-radlw_param.$(OBJEXT) : libgfs_a-physparam.$(OBJEXT) 
libgfsf_a-radsw_datatb.$(OBJEXT) : libgfsf_a-radsw_param.$(OBJEXT) libgfs_a-physparam.$(OBJEXT) 
libgfsf_a-radsw_main.$(OBJEXT) : libgfsf_a-radsw_param.$(OBJEXT) libgfs_a-physparam.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) \
	libgfsf_a-radsw_datatb.$(OBJEXT) 
libgfsf_a-radsw_param.$(OBJEXT) : libgfs_a-physparam.$(OBJEXT) 
libgfsf_a-rascnvv2.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-rayleigh_damp.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 
libgfsf_a-rayleigh_damp_mesopause.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 
libgfsf_a-sascnv.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-sascnvn.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-set_soilveg.$(OBJEXT) : libgfsf_a-namelist_soilveg.$(OBJEXT) 
libgfsf_a-sfc_cice.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-sfc_diag.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-sfc_diff.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-sfc_drv.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-sfc_land.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-sfc_nst.$(OBJEXT) : libgfs_a-module_nst_model.$(OBJEXT) libgfsf_a-date_def.$(OBJEXT) libgfs_a-machine.$(OBJEXT) \
	libgfs_a-module_nst_parameters.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-module_nst_water_prop.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-sfc_ocean.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-sfc_sice.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-sfcsub.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 
libgfsf_a-sflx.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfsf_a-namelist_soilveg.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-shalcnv.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) 
libgfsf_a-shalcv.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-shalcv_opr.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-tracer_const_h.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) 
libgfsf_a-tridi2t3.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) 

# idea
libidea_a-w05sc_efield_merge.$(OBJEXT) : libidea_a-idea_io_units.$(OBJEXT)
libidea_a-efield.$(OBJEXT) : libidea_a-idea_io_units.$(OBJEXT) libidea_a-namelist_wamphysics_def.$(OBJEXT) libidea_a-w05sc_efield_merge.$(OBJEXT) 
libidea_a-idea_cal_advance.$(OBJEXT) : libidea_a-idea_io_units.$(OBJEXT) libgfsf_a-date_def.$(OBJEXT) libidea_a-idea_mpi_def.$(OBJEXT) 
libidea_a-idea_co2.$(OBJEXT) : libidea_a-idea_composition.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libidea_a-co2hc.$(OBJEXT) 
libidea_a-idea_dissipation.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libidea_a-idea_composition.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libidea_a-namelist_wamphysics_def.$(OBJEXT)
libidea_a-idea_getno_snoe.$(OBJEXT) : libidea_a-idea_io_units.$(OBJEXT) libidea_a-idea_composition.$(OBJEXT) libidea_a-idea_mpi_def.$(OBJEXT) 
libidea_a-idea_h2o.$(OBJEXT) : libidea_a-idea_composition.$(OBJEXT) 
libidea_a-idea_imf_input.$(OBJEXT) : libidea_a-idea_io_units.$(OBJEXT) libidea_a-idea_mpi_def.$(OBJEXT) libidea_a-idea_cal_advance.$(OBJEXT) 
libidea_a-idea_interpol_datab.$(OBJEXT) : 
libidea_a-idea_ion.$(OBJEXT) : libidea_a-idea_io_units.$(OBJEXT) libgfsf_a-date_def.$(OBJEXT) libidea_a-idea_composition.$(OBJEXT) libidea_a-idea_ion_input.$(OBJEXT) libidea_a-idea_imf_input.$(OBJEXT) libidea_a-idea_mpi_def.$(OBJEXT) libidea_a-efield.$(OBJEXT) 
libidea_a-idea_ion_empirmodels.$(OBJEXT) : libidea_a-idea_ion_input.$(OBJEXT) libidea_a-idea_composition.$(OBJEXT) libidea_a-idea_mpi_def.$(OBJEXT) 
libidea_a-idea_ion_input.$(OBJEXT) : libidea_a-idea_io_units.$(OBJEXT) libidea_a-idea_mpi_def.$(OBJEXT) 
libidea_a-idea_merge_ipe_to_wam.$(OBJEXT) : 
libidea_a-idea_mpi_def.$(OBJEXT) : libidea_a-idea_io_units.$(OBJEXT) 
libidea_a-idea_o2_o3.$(OBJEXT) : libidea_a-idea_composition.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libidea_a-idea_solar_init.$(OBJEXT) 
libidea_a-idea_phys.$(OBJEXT) : libidea_a-namelist_wamphysics_def.$(OBJEXT) libidea_a-idea_composition.$(OBJEXT) libidea_a-idea_imf_input.$(OBJEXT) libidea_a-idea_mpi_def.$(OBJEXT) libidea_a-idea_ion.$(OBJEXT) libidea_a-idea_cal_advance.$(OBJEXT) libidea_a-idea_solar_input.$(OBJEXT) 
libidea_a-idea_solar_2014.$(OBJEXT) : libidea_a-idea_composition.$(OBJEXT) 
libidea_a-idea_solar_heating.$(OBJEXT) : libgfsf_a-date_def.$(OBJEXT) libgfs_a-machine.$(OBJEXT) libidea_a-idea_composition.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libidea_a-idea_ion.$(OBJEXT) libidea_a-idea_cal_advance.$(OBJEXT) libidea_a-idea_solar_init.$(OBJEXT) 
libidea_a-idea_solar_init.$(OBJEXT) : libidea_a-idea_io_units.$(OBJEXT) libidea_a-idea_getno_snoe.$(OBJEXT) libidea_a-namelist_wamphysics_def.$(OBJEXT) libgfs_a-machine.$(OBJEXT) libidea_a-idea_composition.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libidea_a-idea_mpi_def.$(OBJEXT) libidea_a-idea_cal_advance.$(OBJEXT) libidea_a-idea_solar_input.$(OBJEXT) 
libidea_a-idea_solar_input.$(OBJEXT) : libidea_a-idea_io_units.$(OBJEXT) libidea_a-idea_mpi_def.$(OBJEXT) libidea_a-idea_cal_advance.$(OBJEXT) 
libidea_a-idea_tracer.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libidea_a-idea_composition.$(OBJEXT) libidea_a-idea_tracers_input.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libidea_a-idea_mpi_def.$(OBJEXT) 
libidea_a-idea_tracers_input.$(OBJEXT) : libidea_a-idea_io_units.$(OBJEXT) libidea_a-idea_composition.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libidea_a-idea_mpi_def.$(OBJEXT) 
libidea_a-idea_trad_geopgrav.$(OBJEXT) : libgfs_a-machine.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfsf_a-tracer_const_h.$(OBJEXT) libidea_a-idea_solar_init.$(OBJEXT) 

# drv
libgfsdrv_a-nuopc_physics.$(OBJEXT) : libgfsf_a-radsw_param.$(OBJEXT) libgfs_a-machine.$(OBJEXT) \
	libgfsf_a-radlw_param.$(OBJEXT) libgfs_a-physcons.$(OBJEXT) libgfs_a-funcphys.$(OBJEXT) \
	libgfsf_a-grrad.$(OBJEXT)
