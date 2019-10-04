noinst_LIBRARIES = libgen.a

libgen_a_SOURCES = module_GEN_GRID_COMP_stub.F90 module_GEN_INTERNAL_STATE.F90

libgen_a_FCFLAGS = $(FC_MODINC) ../share

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libgen_a-module_GEN_GRID_COMP_stub.$(OBJEXT) : libgen_a-module_GEN_INTERNAL_STATE.$(OBJEXT)
