include ../../macros.make

all:
	(cd global_sfchdr.fd && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" LIBDIR=$(LIBDIR) \
    BINDIR=$(BINDIR) FREEFLAG=$(FREEFLAG)) || \
    (echo "make failure in global_sfchdr.fd/" && exit 1)
	(cd global_sighdr.fd && $(MAKE) FC=$(FC) FFLAGS="$(FFLAGS)" LIBDIR=$(LIBDIR) \
    BINDIR=$(BINDIR) FREEFLAG=$(FREEFLAG)) || \
    (echo "make failure in global_sighdr.fd/" && exit 1)

clean:
	(cd global_sfchdr.fd && $(MAKE) clean)
	(cd global_sighdr.fd && $(MAKE) clean)
