# horizontal Makefile

include ../../macros.make
include FIM_HORIZONTAL_OBJS
include FIM_HORIZONTAL_OBJS_TOP

AR              = ar  # TODO: move into top-level makefile or macros.make
ARFLAGS         = ruv # TODO: move into top-level makefile or macros.make
COLMOBJS        = $(shell ls ../column/*.o)
COLUMN_CHEM_OBJS= $(wildcard ../column_chem/*.o)
DEPLIBS         = $(LIBDIR)/libbacio_4.a $(LIBDIR)/libsfcio_4.a                \
                  $(LIBDIR)/libsigio_4.a $(LIBDIR)/libslint.a                  \
                  $(LIBDIR)/libsp_4.a $(LIBDIR)/libsysshare.a                  \
                  $(LIBDIR)/libw3_4.a
ICOSIO          = ../../icosio/icosio.o
LIBS            = -l$(FIMLIBBASE) -L$(LIBDIR) -lbacio_4 -lgribio -lsfcio_4     \
                  -lsigio_4 -lslint -lsp_4 -lsysshare -lw3_4 -lwrfio
LOCFLAGS        = $(FFLAGS) $(FREEFLAG) -I$(SMS)/include -I$(UTILSDIR)         \
                  -I../../icosio -I../../post/pop -I../../prep/sfcio           \
                  -I../../prep/sigio -I../../prep/incmod
OBJS            = $(OBJS_BASE) $(OBJS_TOP)
POSTOBJS        = $(foreach x,fimnc post smooth postdata,../../post/pop/$(x).o)
UTILS           = $(UTILSDIR)/read_queue_namelist.o $(UTILSDIR)/headers.o
UTILSDIR        = ../../utils
WRFPHYS_OBJS    = $(shell ls ../wrfphys/*.o | grep -v module_cu_g3.o)

# Add -g for debugging.

ifeq ($(DEBUG),yes)
 LOCFLAGS += -g
endif

# For serial, use no-op 'true' (instead of 'ppp') and no extension.

ifeq ($(P),S)
  SMSEXTENSION =
  PPP = true
else
  SMSEXTENSION = _sms
  PPP = $(SMS)/bin/ppp
  # HaloSize=2000 is a guess at what will be large enough for all combinations
  # of NIP, nPEs, and decomp. GlobalSize=2621442 is for G9 and below. Higher G
  # levels will require increasing GlobalSize. GlobalSize=6000000 is for G9.5.
  # Higher G levels will require increasing GlobalSize. --Fmodule=module_decomp
  # tells SMS to insert "use module_decomp" in every file it processes. This is
  # needed to 1) not require the user to add it everywhere; 2) it is not needed
  # in serial code; 3) module decomp contains all the decomp info that SMS needs
  # in the routines it processes.
  PPP_FLAGS1 = --comment --FreeFormat --HaloSize=2000 --GlobalSize=2621442
  PPP_FLAGS = $(PPP_FLAGS1) --Fmodule=module_decomp
  PPP_HEADER_FLAGS = --header --comment
endif

# If degree-based trig intrinsics don't exist, set ifdef to define them as stmt
# functions.

ifeq ($(NEED_SINDCOSD),yes)
  CPP_FLAGS += -DNEED_SINDCOSD
endif

# Are we going to attach gdb or ddd to a running process?

ifeq ($(ATTACH_DEBUGGER),yes)
  LOCFLAGS += -DATTACH_DEBUGGER
endif

%.o: %.f90
	$(CPP) $(GPTL_CPPFLAGS) $(CPP_FLAGS) $< > $*.f
	$(FC) -c $(GPTL_FFLAGS) $(LOCFLAGS) $(OPTFLAGS) $*.f

%.o: %.F90
	$(CPP) $(GPTL_CPPFLAGS) $(CPP_FLAGS) $< > $*_cpp.f
	$(PPP) $(PPP_FLAGS) $*_cpp.f # use --debug for added PPP debugging
	$(FC) -c $(GPTL_FFLAGS) $(LOCFLAGS) $(OPTFLAGS) $*_cpp$(SMSEXTENSION).f -o $*.o

%.o: %.c
	$(CC) -c -g $(CFLAGS) $(FORTRAN_UNDERSCORING) $<

all: $(FIMEXE)

physics.o: $(DEPLIBS)

-include FIM_HORIZONTAL_DEPENDENCIES

$(FIMEXE): $(FIMLIB) $(OBJS_TOP) $(DEPLIBS)
	$(FC) -o $(FIMEXE) $(FFLAGS) -L$(SMS)/lib $(LINKFLAGS) $(OBJS_TOP) \
	$(LIBS) -lsms $(LINKLIBS) $(GPTL_LDFLAGS) $(LDFLAGS) $(LIBNETCDF)

# Build a library of all FIM object files except $(OBJS_TOP).  This is needed
# when building FIM from within the NEMS repository but would not be needed by
# other builds.  However, always building $(FIMLIB) significantly simplifies
# build automation (for example, no additional $(FTNMPI) values need to be
# supported by makefim, special logic is not needed in fim/horizontal/copy.ksh
# to avoid building $(OBJS_TOP) files, additional logic is not needed in the
# top-level Makefile, etc.).  And, since we will always build it, we should
# use it to build $(FIMEXE) to avoid using a separate but equivalent approach
# when building from within NEMS.

$(FIMLIB): $(OBJS_BASE) $(COLMOBJS) $(POSTOBJS) $(COLUMN_CHEM_OBJS) \
           $(WRFPHYS_OBJS) $(UTILS) $(ICOSIO)
	$(AR) $(ARFLAGS) $@ $?

module_control.o module_decomp.o module_wrf_control.o:
	$(CPP) $(GPTL_CPPFLAGS) $(CPP_FLAGS) $< > $*_cpp.f
	$(PPP) $(PPP_FLAGS1) $*_cpp.f
	$(FC) -c $(GPTL_FFLAGS) $(LOCFLAGS) $(OPTFLAGS) $*_cpp$(SMSEXTENSION).f -o $*.o

# The following -O0 rule ensures bitwise-exact output between serial & parallel
# runs built using ifort 9.1.

module_chem_driver.o:
	$(CPP) $(GPTL_CPPFLAGS) $(CPP_FLAGS) $< > $*_cpp.f
	$(PPP) $(PPP_FLAGS) $*_cpp.f
	$(FC) -c $(GPTL_FFLAGS) $(LOCFLAGS) $(OPTFLAGS) -O0 $*_cpp$(SMSEXTENSION).f -o $*.o
