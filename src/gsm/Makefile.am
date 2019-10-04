SUBDIRS = libutil dyn phys .
noinst_LIBRARIES = libgsm.a

libgsm_a_SOURCES = module_GFS_CORE_SETUP.F90 module_GFS_GRID_COMP.F90 module_GFS_INTEGRATE.F90 \
	module_GFS_INTERNAL_STATE.F90 module_GOCART_ROUTINES_stub.F90

libgsm_a_FCFLAGS = -I $(top_srcdir)/include $(FC_MODINC) dyn $(FC_MODINC) libutil $(FC_MODINC) phys \
	$(FC_MODINC) ../share $(FC_MODINC) ../io
libgsm_a_FCFLAGS += $(FC_AUTODOUBLE) $(FCFLAGS_FREE)

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libgsm_a-module_GFS_GRID_COMP.$(OBJEXT) : libgsm_a-module_GFS_INTEGRATE.$(OBJEXT) \
	libgsm_a-module_GFS_INTERNAL_STATE.$(OBJEXT) libgsm_a-module_GFS_CORE_SETUP.$(OBJEXT)
libgsm_a-module_GFS_GRID_COMP.$(OBJEXT) : libgsm_a-module_GFS_CORE_SETUP.$(OBJEXT) libgsm_a-module_GFS_INTEGRATE.$(OBJEXT) \
	libgsm_a-module_GFS_INTERNAL_STATE.$(OBJEXT) libgsm_a-module_GOCART_ROUTINES_stub.$(OBJEXT)
libgsm_a-module_GFS_INTEGRATE.$(OBJEXT) : libgsm_a-module_GOCART_ROUTINES_stub.$(OBJEXT)
