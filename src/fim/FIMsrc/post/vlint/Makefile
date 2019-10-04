# vlint Makefile

SHELL = /bin/sh

AR       = ar ruv
FFLAGS	 = $(FFLAGS)
LIBVLINT = $(LIBDIR)/libvlint.a

.SUFFIXES: .F90 .o .a

$(LIBVLINT): $(LIBVLINT)(vlint.o)

.F90.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $@ $*.o
	rm -f $*.o

clean:
	rm -f *.o 
