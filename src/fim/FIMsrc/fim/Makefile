# fim Makefile

include ../macros.make

CNTLMODDIR = ../../cntl/incmod
FIMEXE     = $(BINDIR)/fim$(P)
FIMLIBBASE = fim$(P)
FIMLIB    = $(LIBDIR)/lib$(FIMLIBBASE).a
FIMHRZL    = Horizontal$(P)
SHELL      = /bin/sh

#JR Invoke parallel capability of gmake where possible.

all:
	(cd column && ./copy.ksh && $(MAKE) $(GMAKEMINUSJ) FC=$(FC) \
		FLAGS="$(COLFLAGS)" FREEFLAG=$(FREEFLAG) FIXEDFLAG=$(FIXEDFLAG)) || \
		(echo "Make failure in column/" && exit 1)

	(cd column_chem && $(MAKE) $(GMAKEMINUSJ) FC=$(FC) FLAGS="$(FFLAGS) \
		-I../column -I$(CNTLMODDIR)" FREEFLAG=$(FREEFLAG)) || \
		(echo "Make failure in column_chem/" && exit 1)

	(cd wrfphys && $(MAKE) $(GMAKEMINUSJ) FC=$(FC) FLAGS="$(FFLAGS) -I../column \
		-I$(CNTLMODDIR)" FREEFLAG=$(FREEFLAG)) || \
		(echo "Make failure in wrfphys/" && exit 1)

	mkdir -p $(FIMHRZL)

	(cd $(FIMHRZL) && env NEMS=$(NEMS) ../horizontal/copy.ksh && $(MAKE) \
		MAKE=$(MAKE) FC=$(FC) FFLAGS="$(FFLAGS) -I../column -I../column_chem \
		-I../wrfphys -I../../w3" FIMEXE=$(FIMEXE) FREEFLAG=$(FREEFLAG) SMS=$(SMS) \
		LINKFLAGS="$(LINKFLAGS)" LINKLIBS="$(LINKLIBS)" LIBWRFP=$(LIBWRFP) \
		FIMLIB=$(FIMLIB) FIMLIBBASE=$(FIMLIBBASE)) || \
		(echo "Make failure in $(FIMHRZL)" && exit 1)

clean:
	(cd column && $(MAKE) clean)
	(cd column_chem && $(MAKE) clean)
	(cd wrfphys && $(MAKE) clean)
	($(RM) -rf Horizontal*)
