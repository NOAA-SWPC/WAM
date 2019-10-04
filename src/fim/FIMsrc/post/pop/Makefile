# pop Makefile

include ../../macros.make
SHELL = /bin/sh

DEPLIBS     = $(LIBBACIO) $(LIBW3) $(LIBSLINT) $(LIBVLINT) $(LIBGRIBIO) $(LIBCNTL) $(LIBSYSSHARE) $(LIBWRFIO)
INCS        = -I ../../cntl -I../../prep/incmod -I $(UTILDIR)
LIBBACIO    = $(LIBDIR)/libbacio_4.a
LIBCNTL     = $(LIBDIR)/libcntl.a
LIBGRIBIO   = $(LIBDIR)/libgribio.a
LIBS        = -L$(LIBDIR) -lbacio_4 -lw3_4 -lslint  -lvlint -lgribio -lwrfio -lw3_4 -lcntl -lsysshare
LIBSLINT    = $(LIBDIR)/libslint.a
LIBSYSSHARE = $(LIBDIR)/libsysshare.a
LIBVLINT    = $(LIBDIR)/libvlint.a
LIBW3       = $(LIBDIR)/libw3_4.a
LIBWRFIO    = $(LIBDIR)/libwrfio.a
POP         = $(BINDIR)/pop
GET_GRIBOUT = $(BINDIR)/get_gribout
UTILDIR     = ../../utils
UTILOBJS    = $(UTILDIR)/read_queue_namelist.o \
              $(UTILDIR)/module_initial_chem_namelists.o \
              $(UTILDIR)/headers.o

.SUFFIXES:
.SUFFIXES: .F90 .o

#JR If OPTFLAGS not empty, hopefully compiler will override contradictory FFLAGS settings
.F90.o:
	$(FC) -c $(FFLAGS) $(OPTFLAGS) $(INCS) $< $(LIBNETCDF) $(INCNETCDF)

all: $(POP) $(GET_GRIBOUT) post.o 

$(POP): pop.F90 smooth.o fimnc.o postdata.o $(DEPLIBS)
	$(FC) $(FFLAGS) $(OPTFLAGS) $(INCS) -o $(POP) smooth.o fimnc.o postdata.o $(UTILOBJS) pop.F90 $(LIBS) $(LIBNETCDF) $(INCNETCDF)

$(GET_GRIBOUT): get_gribout.o postdata.o
	$(FC) $(FFLAGS) -o $@ get_gribout.o postdata.o $(LIBSYSSHARE)

clean:
	$(RM) *.o *.mod
#post.F90 is the subroutine pop interface to FIM
post.o: post.F90 fimnc.o postdata.o
get_gribout.o: get_gribout.F90 postdata.o
