noinst_LIBRARIES = libpost.a

libpost_a_SOURCES = post_gfs_stub.F90  post_nems_routines_stub.F90  post_nmm_stub.F90

libpost_a_FCFLAGS =

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules
