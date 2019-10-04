noinst_LIBRARIES = libshare.a

libshare_a_SOURCES = check_esmf_pet.F90 ffsync.F90 get_variables_for_wam_ipe_coupling.F90 \
	idea_composition.f kind.inc mersenne_twister.f module_CONSTANTS.F90 module_CPLFIELDS.F90 \
	module_ERR_MSG.F90 module_INCLUDE.F90 module_IPE_to_WAM.F90

libshare_a_FCFLAGS = -I $(top_srcdir)/include
libshare_a_FFLAGS = $(FC_AUTODOUBLE)

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libshare_a-module_CONSTANTS.$(OBJEXT) : libshare_a-module_INCLUDE.$(OBJEXT)
