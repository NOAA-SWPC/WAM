# ssfc2icos Makefile

include ../../macros.make

SHELL = /bin/sh

DEPLIBS    = $(LIBW3) $(LIBSP) $(LIBSIGIO) $(LIBSFCIO) $(LIBSLINT) $(LIBSYSSHARE)
FLAG1      = $(FFLAGS)
FLAGS      = $(FLAG1) $(INCS) $(DEBUG_FLAGS)
INCS       = -I ../incmod -I../../utils -I../../cntl -I../sfcio
LIBS       = -L$(LIBDIR) -lw3_4 -lsp_4 -lsigio_4 -lsfcio_4 -lslint 
LIBSFCIO   = $(LIBDIR)/libsfcio_4.a
LIBSIGIO   = $(LIBDIR)/libsigio_4.a
LIBSLINT   = $(LIBDIR)/libslint.a
LIBSP      = $(LIBDIR)/libsp_4.a
LIBSYSSHARE= $(LIBDIR)/libsysshare.a
LIBW3      = $(LIBDIR)/libw3_4.a
NEWNAME    = $(BINDIR)/newname.exe
SSFC2ICOS  = $(BINDIR)/ssfc2icos
UTILDIR    = ../../utils

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FLAGS) $<

all: $(SSFC2ICOS) $(NEWNAME)

$(NEWNAME): newname.F90 $(DEPLIBS) $(UTILDIR)/headers.o $(UTILDIR)/read_queue_namelist.o
	$(FC) $(FLAGS) -o $(NEWNAME) $(INCS) -I../../cntl/incmod newname.F90 $(LIBS) $(LIBSYSSHARE) -lcntl $(UTILDIR)/headers.o $(UTILDIR)/read_queue_namelist.o $(UTILDIR)/module_initial_chem_namelists.o

$(SSFC2ICOS): ssfc2icos.F90 $(DEPLIBS) $(UTILDIR)/headers.o $(UTILDIR)/read_queue_namelist.o read_mtnvar.o
	$(FCserial) $(FLAGS) -o $(SSFC2ICOS) $(INCS) ssfc2icos.F90 $(LIBS) $(LIBSLINT) $(LIBSYSSHARE) $(UTILDIR)/headers.o $(UTILDIR)/read_queue_namelist.o $(UTILDIR)/module_initial_chem_namelists.o read_mtnvar.o

read_mtnvar.o: read_mtnvar.F90
	$(FC) -c $(FLAGS) $(BYTE_SWAP_FLAG) read_mtnvar.F90

clean:
	$(RM) *.o *.mod

