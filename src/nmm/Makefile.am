noinst_LIBRARIES = libnmm.a

libnmm_a_SOURCES = module_NMM_GRID_COMP_stub.F90

libnmm_a_FCFLAGS =

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
