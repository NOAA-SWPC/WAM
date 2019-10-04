# grid Makefile

SHELL = /bin/sh

include ../../macros.make

# Build a grid generator:

FLAGS       = $(FFLAGS) $(INCS)
GETLVL      = $(BINDIR)/getlvl
GINFO       = $(BINDIR)/ginfo
GRID        = $(BINDIR)/grid
GRID_OBJS   = icos.o bisect_triangle.o trisect_triangle.o top_gridpoints.o \
              middle.o trisect.o third.o avg.o top_triangles.o rotate.o datastru.o
GINFO_OBJS  = GridGen.o hilbert.o ijblock.o perm.o IJblockLayout.o \
              GetRegions.o SquareDecomp.o SquareLayout.o datastru.o mpi_stubs.o
INCS        = -I ../../cntl -I $(UTILDIR)
LIBCNTL     = $(LIBDIR)/libcntl.a
LIBSYSSHARE = $(LIBDIR)/libsysshare.a
UTILDIR     = ../../utils
UTILOBJS    = $(UTILDIR)/read_queue_namelist.o $(UTILDIR)/headers.o \
              $(UTILDIR)/module_initial_chem_namelists.o

#JR Only need for deleting default suffix rule is to add INCS to FFLAGS

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FLAGS) $<

all: $(GRID) $(GINFO) $(GETLVL)

$(GRID): DEPENDENCIES $(LIBCNTL)(module_control.o) $(LIBSYSSHARE)(sys_share.o) $(GRID_OBJS) $(UTILOBJS)
	 $(FCserial) $(FLAGS) -o $(GRID) $(GRID_OBJS) $(UTILOBJS) $(LIBCNTL) $(LIBSYSSHARE)

$(GINFO): DEPENDENCIES $(LIBCNTL)(module_control.o) $(LIBSYSSHARE)(sys_share.o) $(GINFO_OBJS) $(UTILOBJS)
	  $(FCginfo) $(FLAGS) -o $(GINFO) $(GINFO_OBJS) $(LIBCNTL) $(LIBSYSSHARE) $(UTILOBJS) -L$(SMS)/lib -lsms $(LDFLAGS)

$(GETLVL): DEPENDENCIES getlvl.o  $(LIBCNTL)(module_control.o) $(LIBSYSSHARE)(sys_share.o) $(UTILOBJS)
	$(FCserial) $(FLAGS) -o $(GETLVL) getlvl.o $(LIBCNTL) $(LIBSYSSHARE) $(UTILOBJS)

DEPENDENCIES:
	$(RM) -f Filepath Srcfiles
	echo "." > Filepath
	ls -1 *.F90 > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

-include DEPENDENCIES

clean:
	$(RM) *.o *.mod DEPENDENCIES
