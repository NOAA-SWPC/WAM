noinst_LIBRARIES = libutil.a

libutil_a_SOURCES = atmos_phy_chem_cpl_comp_mod_stub.F90 atmos_chem_phy_cpl_comp_mod_stub.F90 \
	module_gfs_machine.F90 module_include_gfs.F90 module_gfs_mpi_def.F90 module_gfs_physcons.F90 \
	module_gfs_funcphys.F90 module_gfs_tropp.F90 module_DM_PARALLEL_GFS.F90 \
	module_digital_filter_gfs.F90 module_export_import_data.F90 gfs_dyn_phy_cpl_comp_mod.F90

libutil_a_FCFLAGS = $(FC_MODINC) ../../share
libutil_a_FCFLAGS += $(FC_AUTODOUBLE) $(FCFLAGS_FREE)

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libutil_a-gfs_dyn_phy_cpl_comp_mod.$(OBJEXT) : libutil_a-module_export_import_data.$(OBJEXT)
libutil_a-module_DM_PARALLEL_GFS.$(OBJEXT) : libutil_a-module_gfs_mpi_def.$(OBJEXT) libutil_a-module_gfs_machine.$(OBJEXT)
libutil_a-module_gfs_funcphys.$(OBJEXT) : libutil_a-module_gfs_machine.$(OBJEXT) libutil_a-module_gfs_physcons.$(OBJEXT)
libutil_a-module_gfs_physcons.$(OBJEXT) : libutil_a-module_gfs_machine.$(OBJEXT)
libutil_a-module_gfs_tropp.$(OBJEXT) : libutil_a-module_gfs_machine.$(OBJEXT) libutil_a-module_gfs_physcons.$(OBJEXT)
libutil_a-module_include_gfs.$(OBJEXT) : libutil_a-module_gfs_machine.$(OBJEXT)
