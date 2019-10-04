# gribio Makefile

# This makefile builds library of subroutines that writes grib files.  

include ../../macros.make
SHELL = /bin/sh

AR        = ar ruv
CC        = cc
FFLAGS	  = $(FFLAGS)
LIBGRIBIO = $(LIBDIR)/libgribio.a
OBJS      = $(addsuffix .o, $(basename $(SRCS)))
SRCS      = $(shell ls *.F90 *.c)

.SUFFIXES: 
.SUFFIXES: .F90 .c .o .a

#JR Changed to enable parallel make, and to generate dependencies
#JR automatically. Parallel make enabled by building the library
#JR with a single $(AR) cmd

all: $(LIBGRIBIO)

$(LIBGRIBIO): DEPENDENCIES $(OBJS)
	$(AR) $@ $(OBJS)

#JR If OPTFLAGS not empty, hopefully compiler will override contradictory FFLAGS settings
.F90.o: 
	$(FC) -c $(FFLAGS) $(OPTFLAGS) $<

.c.o:
	$(CC) $(CFLAGS) -c $(GRIBIO_CPP_FLAGS) $<
#JR	$(CC) $(CFLAGS) -c -g -O0 -UDEBUG $(GRIBIO_CPP_FLAGS) $<

DEPENDENCIES:
	$(RM) -f Filepath Srcfiles
	echo "." > Filepath
	ls -1 *.F90 > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

-include DEPENDENCIES

clean:
	$(RM) *.o *.mod DEPENDENCIES
