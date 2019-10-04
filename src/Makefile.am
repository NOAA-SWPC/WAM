SUBDIRS = share phys io fim gen gsm nmm .

lib_LIBRARIES = libwam.a

libwam_a_SOURCES = module_ATM_GRID_COMP.F90  module_ATM_INTERNAL_STATE.F90

libwam_a_FCFLAGS = -I $(top_srcdir)/include \
	$(FC_MODINC) fim $(FC_MODINC) gen $(FC_MODINC) gsm $(FC_MODINC) nmm \
	$(FC_MODINC) share $(FC_MODINC) phys $(FC_MODINC) io

libwam_a_LIBADD = \
	fim/*.$(OBJEXT) \
	gen/*.$(OBJEXT) \
	gsm/dyn/*.$(OBJEXT) \
	gsm/libutil/*.$(OBJEXT) \
	gsm/phys/*.$(OBJEXT) \
	io/*.$(OBJEXT) \
	nmm/*.$(OBJEXT) \
	phys/*.$(OBJEXT)

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules


libwam_a-module_ATM_GRID_COMP.$(OBJEXT) : libwam_a-module_ATM_INTERNAL_STATE.$(OBJEXT)
