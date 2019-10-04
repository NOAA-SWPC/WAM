noinst_LIBRARIES = libio.a

libio_a_SOURCES = module_GET_CONFIG_WRITE_GFS.F90 module_GFS_WRITE.F90 module_include_io.F90 \
	module_io_mpi_def.F90 module_WRITE_GRID_COMP_GFS.F90 module_WRITE_INTERNAL_STATE_GFS.F90 \
	module_WRITE_ROUTINES_GFS.F90

libio_a_FCFLAGS = $(FC_MODINC) ../share

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libio_a-module_GET_CONFIG_WRITE_GFS.$(OBJEXT) : libio_a-module_io_mpi_def.$(OBJEXT) \
	libio_a-module_WRITE_INTERNAL_STATE_GFS.$(OBJEXT)
libio_a-module_GFS_WRITE.$(OBJEXT) : libio_a-module_include_io.$(OBJEXT) \
	libio_a-module_io_mpi_def.$(OBJEXT) libio_a-module_WRITE_GRID_COMP_GFS.$(OBJEXT)
libio_a-module_WRITE_GRID_COMP_GFS.$(OBJEXT) : libio_a-module_GET_CONFIG_WRITE_GFS.$(OBJEXT) \
	libio_a-module_include_io.$(OBJEXT) libio_a-module_io_mpi_def.$(OBJEXT) \
	libio_a-module_WRITE_INTERNAL_STATE_GFS.$(OBJEXT) libio_a-module_WRITE_ROUTINES_GFS.$(OBJEXT)
libio_a-module_WRITE_ROUTINES_GFS.$(OBJEXT) : libio_a-module_include_io.$(OBJEXT) \
	libio_a-module_io_mpi_def.$(OBJEXT) \
	libio_a-module_WRITE_INTERNAL_STATE_GFS.$(OBJEXT)

