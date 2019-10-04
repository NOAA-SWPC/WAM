SUBDIRS = src
dist_doc_DATA = README

ACLOCAL_AMFLAGS = -I m4

nodist_data_DATA = wam.mk

include_HEADERS = include/wam_defs.h

wam.mk:
	@echo "# ESMF self-describing build dependency makefile fragment" > $(@)
	@echo >> $(@)
	@echo "ESMF_DEP_FRONT     = wamCap"                        >> $(@)
	@echo "ESMF_DEP_INCPATH   = $(includedir)"                 >> $(@)
	@echo "ESMF_DEP_CMPL_OBJS = "                              >> $(@)
	@echo "ESMF_DEP_LINK_OBJS = $(libdir)/libwam.a"            >> $(@)
	@echo "ESMF_DEP_SHRD_PATH = $(IPE_DEP_SHRD_PATH)"          >> $(@)
	@echo "ESMF_DEP_SHRD_LIBS = $(IPE_DEP_SHRD_LIBS)"          >> $(@)
	@echo
	@echo "Finished generating ESMF self-describing build dependency makefile fragment:" $(@)
	@echo

uninstall-local:
	-rm -f $(datarootdir)/wam.mk

CLEANFILES = wam.mk
