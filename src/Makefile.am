SUBDIRS = post share phys io gen gsm .

lib_LIBRARIES = libwam.a

libwam_a_SOURCES = module_ATM_INTERNAL_STATE.F90 wamCap.F90

libwam_a_FCFLAGS = -I $(top_srcdir)/include \
	$(FC_MODINC) gen $(FC_MODINC) gsm \
	$(FC_MODINC) share $(FC_MODINC) phys $(FC_MODINC) io

libwam_a_LIBADD = \
	gen/*.$(OBJEXT) \
	gsm/*.$(OBJEXT) \
	gsm/dyn/*.$(OBJEXT) \
	gsm/libutil/*.$(OBJEXT) \
	gsm/phys/*.$(OBJEXT) \
	io/*.$(OBJEXT) \
	phys/*.$(OBJEXT) \
	post/*.$(OBJEXT) \
	share/*.$(OBJEXT)

.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
libwam_a-wamCap.$(OBJEXT) : libwam_a-module_ATM_INTERNAL_STATE.$(OBJEXT)
