noinst_LIBRARIES = libfim.a

libfim_a_SOURCES = fim_grid_comp_stub.F90

libfim_a_FCFLAGS =

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
