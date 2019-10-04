# slint Makefile

SHELL = /bin/sh

include ../../macros.make

FLAGS       = $(FFLAGS)
LIBSYSSHARE = $(LIBDIR)/libsysshare.a
SLINTEXE    = $(BINDIR)/slint
SLINTLIB    = $(LIBDIR)/libslint.a
INCMOD	 = ../incmod
FLAGS = $(FFLAGS) -I$(INCMOD)

#JR: Blank out all suffix rules. Mac does some weird .s stuff without it.

.SUFFIXES:
.SUFFIXES: .F90 .o

#JR If OPTFLAGS not empty, hopefully compiler will override contradictory FFLAGS settings
.F90.o:
	$(FC) $(FLAGS) $(OPTFLAGS) -c $<

KDOBJ = kd.o
SLINTOBJ = slint.o 

all: $(SLINTLIB) $(SLINTEXE)

$(SLINTLIB): DEPENDENCIES $(KDOBJ) $(SLINTOBJ) 
	$(AR) ruv $(SLINTLIB) $(SLINTOBJ) $(KDOBJ)
	mv -f *.mod $(INCMOD)

$(SLINTEXE): DEPENDENCIES slintest.o $(SLINTLIB)
	$(FC) $(FLAGS) -o $(SLINTEXE) slintest.o $(SLINTLIB) $(LIBSYSSHARE)

DEPENDENCIES:
	$(RM) -f Filepath Srcfiles
	echo "." > Filepath
	ls -1 *.F90 > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

-include DEPENDENCIES

clean:
	$(RM) *.o *.mod *_cpp.F90 DEPENDENCIES

# Ning Wang, Feb 2006
