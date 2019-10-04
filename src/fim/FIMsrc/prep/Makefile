# prep Makefile

SHELL = /bin/sh

include ../macros.make

#JR Use GMAKEMINUSJ where appropriate to enable parallel capability of gmake

all:
	(cd grid && $(MAKE) $(GMAKEMINUSJ) FC=$(FC) FFLAGS="$(FFLAGS)" \
		LIBDIR=$(LIBDIR) BINDIR=$(BINDIR)) || \
		(echo "make failure in grid/" && exit 1)

	(cd sfcio && $(MAKE) FC=$(FC) BYTE_SWAP_FLAG=$(BYTE_SWAP_FLAG) \
		LIBDIR=$(LIBDIR)) || \
		(echo "make failure in sfcio/" && exit 1)

	(cd sigio && $(MAKE) $(GMAKEMINUSJ) FC=$(FC) \
		BYTE_SWAP_FLAG=$(BYTE_SWAP_FLAG) FREEFLAG=$(FREEFLAG) \
		FIXEDFLAG=$(FIXEDFLAG) LIBDIR=$(LIBDIR)) || \
		(echo "make failure in sigio/" && exit 1)

	(cd slint && $(MAKE) $(GMAKEMINUSJ) FC=$(FC) FFLAGS="$(FFLAGS)" \
		LIBDIR=$(LIBDIR) BINDIR=$(BINDIR)) || \
		(echo "make failure in slint/" && exit 1)

	(cd sp && $(MAKE) $(GMAKEMINUSJ) MAKE=$(MAKE) FC=$(FC) FFLAGS="$(SP_FFLAGS)" \
		FREEFLAG=$(FREEFLAG) FIXEDFLAG=$(FIXEDFLAG) LIBDIR=$(LIBDIR)) || \
		(echo "make failure in sp/" && exit 1)

	(cd ssfc2icos && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" \
		BYTE_SWAP_FLAG=$(BYTE_SWAP_FLAG) LIBDIR=$(LIBDIR) BINDIR=$(BINDIR)) || \
		(echo "make failure in ssfc2icos/" && exit 1)

	(cd gfsenkf && $(MAKE) FC=$(FCserial) FFLAGS="$(FFLAGS)" LIBDIR=$(LIBDIR) \
		BINDIR=$(BINDIR) FREEFLAG=$(FREEFLAG)) || \
		(echo "make failure in gfsenkf/" && exit 1)

clean:
	(cd gfsenkf   && $(MAKE) clean)
	(cd grid      && $(MAKE) clean)
	(cd sfcio     && $(MAKE) clean)
	(cd sigio     && $(MAKE) clean)
	(cd slint     && $(MAKE) clean)
	(cd sp        && $(MAKE) clean)
	(cd ssfc2icos && $(MAKE) clean)
