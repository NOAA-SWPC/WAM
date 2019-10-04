# wrfio Makefile

SHELL = /bin/sh

include ../../macros.make

AR       = ar ruv
LIBWRFIO = $(LIBDIR)/libwrfio.a

.SUFFIXES: .F90 .o .a

$(LIBWRFIO): $(LIBWRFIO)(wrf_io.o field_routines.o )

#JR On Mac, cpp doesn't handle concatenation (##) corectly,
#JR and "cc -E" only does the right thing when handed a .c file. Therefore
#JR need to temporarily rename the .F90 file to .c, and for Mac override
#JR the default CPP with "cc -E"

ifeq ($(MAKEFIM_TARGET),macgnu)
  CPP = cc -E
endif

.F90.a:
	cp -f $< $*.c
	$(CPP) $(CPP_FLAGS) $*.c > $*_cpp.F90
	$(FC) -c $(FFLAGS) $*_cpp.F90 $(LIBNETCDF) $(INCNETCDF)
	mv $*_cpp.o $*.o
	$(AR) $@ $*.o
	$(RM) $*.o $*.c
clean:
	$(RM) *.o *.mod *_cpp.F90 *.c
