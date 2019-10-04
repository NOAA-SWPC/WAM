# column Makefile

SHELL = /bin/sh

include ../../macros.make

OBJS = $(addsuffix .o, $(basename $(SRCS)))
SRCS = $(shell ls *.f *.F90)

.SUFFIXES:
.SUFFIXES: .o .f .F90

ifeq ($(NEED_SINDCOSD),yes)
  LOCDEF = -DNEED_SINDCOSD
endif

LOCFLAGS = $(FLAGS)
ifeq ($(DEBUG),yes)
  LOCFLAGS += -g
endif

.f.o:
	$(FC) -c -I../../cntl $(GPTL_FFLAGS) $(LOCFLAGS) $(FIXEDFLAG) $(RCWFLAG) $(OPTFLAGS) $<
.F90.o:
	$(FC) -c $(GPTL_CPPFLAGS) $(GPTL_FFLAGS) $(LOCDEF) $(LOCFLAGS) $(FREEFLAG) $(RCWFLAG) $(OPTFLAGS) $<

#JR Generate FIM_COLUMN_DEPENDENCIES automatically.
#JR That is why it was removed from the repository.
#JR Also: No longer build a library from the column .o files. They are now
#JR called out directly in the "make" procedure in horizontal/.
all: FIM_COLUMN_DEPENDENCIES $(OBJS)

FIM_COLUMN_DEPENDENCIES:
	$(RM) -f Filepath Srcfiles
	echo "." > Filepath
	ls -1 *.f *.F90 > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

-include FIM_COLUMN_DEPENDENCIES

clean:
	$(RM) *.o *.mod iw3jdn.f  w3fs26.f  w3movdat.f  w3reddat.f FIM_COLUMN_DEPENDENCIES
