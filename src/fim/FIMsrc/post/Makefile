# post Makefile

SHELL = /bin/sh

include ../macros.make

all: 
	(cd gribio && $(MAKE) $(GMAKEMINUSJ) FC=$(FC) FFLAGS="$(FFLAGS)" LIBDIR=$(LIBDIR) CFLAGS="$(CFLAGS)")
	(cd vlint  && $(MAKE)                FC=$(FC) FFLAGS="$(FFLAGS)" LIBDIR=$(LIBDIR))
	(cd wrfio  && $(MAKE)                FC=$(FC) FFLAGS="$(FFLAGS)" LIBDIR=$(LIBDIR) LIBNETCDF="$(LIBNETCDF)" INCNETCDF=$(INCNETCDF))
	(cd pop    && $(MAKE)                FC=$(FC) FFLAGS="$(FFLAGS)" LIBDIR=$(LIBDIR) BINDIR=$(BINDIR) LIBNETCDF="$(LIBNETCDF)")

clean:
	(cd gribio && $(MAKE) clean)
	(cd vlint  && $(MAKE) clean)
	(cd wrfio  && $(MAKE) clean)
	(cd pop    && $(MAKE) clean)



