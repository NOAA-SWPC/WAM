# sigio Makefile

SHELL = /bin/sh

include ../../macros.make

FC       = $(FC)
FLAGB    = $(BYTE_SWAP_FLAG) $(FFLAGS_NO_DEBUG)
#FFLAGS are not imported from the calling Makefile because in sigio_module.f the allocate statement 
#fails for "-g -traceback -C" and "-g -traceback -C" does not work with "-mcmodel=medium -i-dynamic"
FLAGS    = $(FFLAGS_NO_DEBUG)
LIBSIGIO = $(LIBDIR)/libsigio_4.a
OBJS     = $(addsuffix .o, $(basename $(SRCS)))
SRCS     = $(shell ls *.f *.F90)

.SUFFIXES:
.SUFFIXES: .a .o .f .F90

.F90.o:
	$(FC) -c $(FLAGB) $(FREEFLAG) $(FFLAGS) $(RCWFLAG) $<

.f.o:
	$(FC) -c $(FLAGS) $(FIXEDFLAG) $(FFLAGS) $(RCWFLAG) $<

all: $(LIBSIGIO)

$(LIBSIGIO): DEPENDENCIES $(OBJS)
	$(AR) ruv $(LIBSIGIO) $(OBJS)

DEPENDENCIES:
	$(RM) -f Filepath Srcfiles
	echo "." > Filepath
	ls -1 *.F90 *.f > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

-include DEPENDENCIES

clean:
	$(RM) -f *.o sigio*.mod $(LIBSIGIO) DEPENDENCIES
