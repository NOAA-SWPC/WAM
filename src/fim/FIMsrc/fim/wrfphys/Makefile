# wrfphys Makefile

SHELL = /bin/sh

include ../../macros.make
include FIM_WRFP_OBJECTS

.SUFFIXES:
.SUFFIXES: .o .f .F .a 

.f.o .F.o:
	$(FC) -c $(FLAGS) -I../../cntl -I../../utils $(FREEFLAG) $<

all: DEPENDENCIES $(OBJS)

DEPENDENCIES:
	$(RM) -f Filepath Srcfiles
	echo "." > Filepath
	ls -1 *.F > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

# The following -O0 rule ensures bitwise-exact output between serial & parallel
# runs built using ifort 9.1.

module_cu_g3.o:
	$(FC) -c $(FLAGS) $(FREEFLAG) -O0 module_cu_g3.F

-include DEPENDENCIES

clean:
	$(RM) -f *.o *.mod *.a DEPENDENCIES
