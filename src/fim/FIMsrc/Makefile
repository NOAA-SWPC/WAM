SHELL = /bin/sh
# Makefile for the FIM system
# Please use the "makefim" script to build FIM rather than the "make" command.
# See makefim and fim_setup.ksh for details on how to build FIM.

#NOTE:  We will eventually need the following gmake subst for IBM but not yet!   
# Set up special characters
# Yes, GNU make is *required*
#null  :=
#space := $(null) $(null)
#comma := $(null),$(null)
# For IBM, use gmake "subst" to convert from cpp "-DFOO -DBAR" to xlf "-WF,-DFOO,-DBAR"
# Only one "-WF" is allowed by xlf so we must append all "-D" flags together before 
# applying this method.  Method borrowed shamelessly from CAM.
#VaporCPPDEF    := "-DFOO -DBAR"
#VaporFPPDEF    := -WF,$(subst $(space),$(comma),$(VaporCPPDEF))

# FIM_ESMF_INSTALL_LIBDIR_ABSPATH is now defined in fim_setup.ksh
ifdef FIM_ESMF_INSTALL_LIBDIR_ABSPATH
# Include makefile fragment "esmf.mk" for ESMF, used by the nems target.
# Note that all variables defined in esmf.mk begin with "ESMF_" and
# should not conflict with other make variables even when nems target
# is not being built.
include $(FIM_ESMF_INSTALL_LIBDIR_ABSPATH)/esmf.mk
endif

# TODO:  separate Intel-openmpi vs. Intel-mvapich vs. ifort here...
FCserialintel = ifort
FC      = mpif90

# TODO:  separate Intel-openmpi vs. Intel-mvapich vs. ifort here...
EJET_VEC_FFLAGS = -O3 -xW -vec-report3
EJET_VEC_FFLAGS = -O3 -xP
#JET_VEC_FFLAGS = -O3 -xT -fp-model precise
#JET_VEC_FFLAGS = -O3 -xT
#JET_VEC_FFLAGS = -O3 -xW
#JET_VEC_FFLAGS = -O3 -xT -fp-model precise -vec-report3
JET_VEC_FFLAGS = -O3
#JET_VEC_FFLAGS = -O0 -g -traceback -C -fpe0
#JET_VEC_FFLAGS = -O0 -g -traceback -C
#VEC_FFLAGS = $(EJET_VEC_FFLAGS)
VEC_FFLAGS = $(JET_VEC_FFLAGS)
#INTEL_FFLAGS = -O2 -DVERBOSE
#INTEL_FFLAGS = $(VEC_FFLAGS)
INTEL_FFLAGS = $(VEC_FFLAGS)
#INTEL_FFLAGS = $(VEC_FFLAGS) -r8

# generic settings
# netcdf:  rely on NETCDF environment variable by default, if present
ifdef NETCDF
LIBNETCDF = "-L$(NETCDF)/lib -lnetcdf"
INCNETCDF = "-I$(NETCDF)/include"
endif
# ESMF
LINK_FLAGS_NOESMF = ""
LINK_LIBS_NOESMF = ""
LINK_FLAGS_ESMF = "$(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS)"
LINK_LIBS_ESMF = "$(ESMF_F90ESMFLINKLIBS)"
#JR changed this: LINK_LIBS puts at end of mpif90 stmt which works
#JR LINK_FLAGS puts at beginning which fails

JAGUAR_LINK_LIBS_NOESMF = ""
FROST_LINK_LIBS_NOESMF = ""

# Intel settings
INTEL_FMVAPICH = $(VEC_FFLAGS)
INTEL_COL_FFLAGS = $(VEC_FFLAGS) -r8
INTEL_ESMF_FFLAGS = "$(INTEL_FFLAGS) $(ESMF_F90COMPILEPATHS)"
PREP_FFLAGS = "-O2"
INTEL_SP_FFLAGS = "-O2"
#PREP_FFLAGS= "-O0 -g -traceback -C -fpe0"
POST_FFLAGS = "-O2"
CCintel = "cc"
CFLAGS= "-O2"
INTEL_FDEBUG = "-O0 -g -traceback -C -fpe0 -ftz"
INTEL_COL_FDEBUG = "-O0 -g -traceback -C -fpe0 -ftz -r8"
INTEL_PREP_DEBUG = "-O0 -g -traceback -C -fpe0 -ftz"
INTEL_POST_DEBUG = "-O0 -g -traceback -C -fpe0 -ftz"
INTEL_GRIBIO_CPP_FLAGS = "-DSUN"
INTEL_FREEFLAG = "-free"
INTEL_FIXEDFLAG = "-fixed"
CDEBUG = "-g"
INTEL_STATIC_FLAG = "-static"
# byte-swapping handled by env flags
INTEL_BYTE_SWAP_FLAG = ""
INTEL_FLUSH_DEF = ""

# Lahey settings
# The sigio library and other libaries are not TKR consistent so they cannot
# use --chkglobal. However --chkglobal must be used in all routines if it is
# used at all so it cannot be used anywhere.
# The sigio and other libraries should be made TKR consistent.
FCseriallahey = lf95
#LAHEY_FDEBUG = "--chk -g --trace --info"
#LAHEY_FDEBUG = "--chk a,e,s,u,x -O0 -g --trace --info"
#LAHEY_FDEBUG = "--chkglobal -g --info"
#LAHEY_COL_FDEBUG= "--chkglobal -g --info --dbl"
LAHEY_FDEBUG = "--chk -g --trap --trace --info"
LAHEY_COL_FDEBUG = "--chk -g --trap --trace --info --dbl"
LAHEY_PREP_DEBUG = "--chk -g --trap --trace --info"
# sp lib has underflow...
LAHEY_SP_FFLAGS = ""
LAHEY_FFLAGS_NO_DEBUG = "-g"
LAHEY_POST_DEBUG = "--chk -g --trap --trace --info"
CClahey = "cc"
CFLAGSlahey = "-O2"
LAHEY_GRIBIO_CPP_FLAGS = "-DSUN"
LAHEY_FREEFLAG = "-Free"
LAHEY_FIXEDFLAG = "-Fixed"
LAHEY_STATIC_FLAG = "-static"
# byte-swapping handled by env flags
LAHEY_BYTE_SWAP_FLAG = ""
LAHEY_FLUSH_DEF = ""

# PGI settings
# TODO:  integrate and then test PGI settings
PGI_FCserial = pgf90
PGI_OPT_FFLAGS = -O2
PGI_FREEFLAG = "-Mfree"
PGI_FIXEDFLAG = "-Mfixed"
PGI_BYTE_SWAP_FLAG = "-byteswapio"
PGI_STATIC_FLAG = ""

# IBM settings
FCserialibm = xlf90_r
FCibm = mpxlf90_r
CCibm = mpcc_r
FCserialvapor = $(FCserialibm)
FCvapor = $(FCibm)
FCserialbluefire = $(FCserialibm)
FCbluefire = $(FCibm)
# Use -qnostrict for similar behavior to other compilers (loosen IEEE restrictions)
IBM_archflags = -q64 -qnostrict -qarch=auto -qspill=32767 -qmaxmem=32767
IBM_FFLAGS = $(IBM_archflags)
IBM_COL_FFLAGS = $(IBM_FFLAGS) -qrealsize=8 -qintsize=4 -O3
IBM_ESMF_FFLAGS ="$(IBM_FFLAGS) $(ESMF_F90COMPILEPATHS)"
IBM_FIM_FFLAGS = $(IBM_FFLAGS) -O3
IBM_ICOSIO_FFLAGS = $(IBM_FFLAGS) -O3
IBM_PREP_FFLAGS = $(IBM_FFLAGS) -O2
IBM_SP_FFLAGS = $(IBM_FFLAGS) -O3
IBM_POST_FFLAGS = $(IBM_FFLAGS) -O2
IBM_GRIBIO_CPP_FLAGS = "-DIBM"
IBM_CFLAGS = $(IBM_archflags) -O2
IBM_FREEFLAG = "-qfree=f90"
IBM_FIXEDFLAG = "-qfixed"
IBM_CDEBUG = "$(IBM_archflags) -g"
# no need to byte swap on IBM
IBM_BYTE_SWAP_FLAG = ""
# try IBM without static linking
IBM_STATIC_FLAG = ""
IBM_FLUSH_DEF = "-WF,-DNO_FLUSH"
# NCEP netcdf paths are special
IBM_LIBNETCDFvapor = "-L/nwprod/lib -lnetcdf_64"
IBM_INCNETCDFvapor = "-I/nwprod/lib/sorc/netcdf/netcdf-3.5.0/include"

#JR gfortran flags (used by linuxpcgnu and macgnu targets)
GF_SP_FFLAGS = "-O2"
GF_FFLAGS = -O2 -ffree-line-length-none
GF_POST_FFLAGS = "-O2"
GF_CC = "cc"
GF_CFLAGS = "-O2"
#JR Most users won't have a gfortran named "gfortran44".
#JR This is special for me on nix, where I needed to have
#JR the most recent compiler installed.
GF_FCserial = gfortran44
GF_COL_FFLAGS = -fdefault-real-8 -O2
GF_FREEFLAG = "-ffree-form"
GF_FIXEDFLAG = "-ffixed-form"
#JR Putting these flags in 2 places is BAD (it's also in macros.make.linuxpcgnu)
#JR For now it'll have to wait till we do some cleanup of the targets
GF_PREP_FFLAGS = "-O2 -frecord-marker=4 -fno-range-check"

MAC_FCserial = gfortran

# generic paths
CWD = $(shell pwd)
LIBDIR = $(CWD)/lib
BINDIR = $(CWD)/bin

# path to mkDepends
MKDEPENDS = $(CWD)/tools/mkDepends

#TODO:  Fix ginfo to eliminate need to link it with SMS at all.  This
#TODO:  will require *all* input and output files to *always* use the
#TODO:  "ij" curve, ALWAYS.  Which will make Rainer and Mike Fiorino
#TODO:  happier.  And inv_perm will have to be used on input if
#TODO:  curve != "ij".  Non-"ij" curves will only be used for performance
#TODO:  testing unless there is a compelling reason to do otherwise
#TODO:  (like a large performance boost due to other curves).
#TODO:  Then eliminate FCginfo here and replace FCginfo with FCserial
#TODO:  in prep/grid/Makefile .

.PHONY: openmpi
openmpi:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			PREP_FFLAGS=$(PREP_FFLAGS)                                \
			FFLAGS_NO_DEBUG="$(PREP_FFLAGS)"                          \
			SP_FFLAGS="$(INTEL_SP_FFLAGS)"                            \
			BYTE_SWAP_FLAG="$(INTEL_BYTE_SWAP_FLAG)"                  \
			STATIC_FLAG="$(INTEL_STATIC_FLAG)"                        \
			FIM_FFLAGS="$(INTEL_FFLAGS)"                              \
			POST_FFLAGS=$(POST_FFLAGS)                                \
			CC=$(CCintel) CFLAGS=$(CFLAGS)                            \
			GRIBIO_CPP_FLAGS=$(INTEL_GRIBIO_CPP_FLAGS)                \
			FC=$(FC) LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)               \
			FCginfo=$(FC)                                             \
			FCserial=$(FCserialintel)                                 \
			LIBNETCDF=$(LIBNETCDF) INCNETCDF=$(INCNETCDF)             \
			COL_FFLAGS="$(INTEL_COL_FFLAGS)"                          \
			FREEFLAG=$(INTEL_FREEFLAG)                                \
			FIXEDFLAG=$(INTEL_FIXEDFLAG)                              \
			LINKFLAGS=$(LINK_FLAGS_NOESMF)                            \
			LINKLIBS=$(LINK_LIBS_NOESMF) NEMS=""                      \
			FLUSH_DEF=$(INTEL_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: mvapich
mvapich:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			PREP_FFLAGS=$(PREP_FFLAGS)                                \
			FFLAGS_NO_DEBUG="$(PREP_FFLAGS)"                          \
			SP_FFLAGS="$(INTEL_SP_FFLAGS)"                            \
			BYTE_SWAP_FLAG="$(INTEL_BYTE_SWAP_FLAG)"                  \
			STATIC_FLAG="$(INTEL_STATIC_FLAG)"                        \
			FIM_FFLAGS="$(INTEL_FMVAPICH)"                            \
			POST_FFLAGS=$(POST_FFLAGS)                                \
			CC=$(CCintel) CFLAGS=$(CFLAGS)                            \
			GRIBIO_CPP_FLAGS=$(INTEL_GRIBIO_CPP_FLAGS)                \
			FC=$(FC) LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)               \
			FCginfo=$(FC)                                             \
			FCserial=$(FCserialintel)                                 \
			LIBNETCDF=$(LIBNETCDF) INCNETCDF=$(INCNETCDF)             \
			COL_FFLAGS="$(INTEL_COL_FFLAGS)"                          \
			FREEFLAG=$(INTEL_FREEFLAG)                                \
			FIXEDFLAG=$(INTEL_FIXEDFLAG)                              \
			LINKFLAGS=$(LINK_FLAGS_NOESMF)                            \
			LINKLIBS=$(LINK_LIBS_NOESMF) NEMS=""                      \
			FLUSH_DEF=$(INTEL_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: debug
debug:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			PREP_FFLAGS=$(INTEL_PREP_DEBUG)                           \
			FFLAGS_NO_DEBUG="$(PREP_FFLAGS)"                          \
			SP_FFLAGS="$(INTEL_SP_FFLAGS)"                            \
			BYTE_SWAP_FLAG="$(INTEL_BYTE_SWAP_FLAG)"                  \
			STATIC_FLAG="$(INTEL_STATIC_FLAG)"                        \
			FIM_FFLAGS=$(INTEL_FDEBUG)                                \
			POST_FFLAGS=$(INTEL_POST_DEBUG)                           \
			CC=$(CCintel) CFLAGS=$(CDEBUG)                            \
			GRIBIO_CPP_FLAGS=$(INTEL_GRIBIO_CPP_FLAGS)                \
			FC=$(FC) LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)               \
			FCginfo=$(FC)                                             \
			FCserial=$(FCserialintel)                                 \
			LIBNETCDF=$(LIBNETCDF) INCNETCDF=$(INCNETCDF)             \
			COL_FFLAGS=$(INTEL_COL_FDEBUG)                            \
			FREEFLAG=$(INTEL_FREEFLAG)                                \
			FIXEDFLAG=$(INTEL_FIXEDFLAG)                              \
			LINKFLAGS=$(LINK_FLAGS_NOESMF)                            \
			LINKLIBS=$(LINK_LIBS_NOESMF) NEMS=""                      \
			FLUSH_DEF=$(INTEL_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: nems
nems:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			PREP_FFLAGS=$(PREP_FFLAGS)                                \
			FFLAGS_NO_DEBUG="$(PREP_FFLAGS)"                          \
			SP_FFLAGS="$(INTEL_SP_FFLAGS)"                            \
			BYTE_SWAP_FLAG="$(INTEL_BYTE_SWAP_FLAG)"                  \
			STATIC_FLAG="$(INTEL_STATIC_FLAG)"                        \
			FIM_FFLAGS=$(INTEL_ESMF_FFLAGS)                           \
			POST_FFLAGS=$(POST_FFLAGS)                                \
			CC=$(CCintel) CFLAGS=$(CFLAGS)                            \
			GRIBIO_CPP_FLAGS=$(INTEL_GRIBIO_CPP_FLAGS)                \
			FC=$(FC) LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)               \
			FCginfo=$(FC)                                             \
			FCserial=$(FCserialintel)                                 \
			LIBNETCDF=$(LIBNETCDF) INCNETCDF=$(INCNETCDF)             \
			COL_FFLAGS="$(INTEL_COL_FFLAGS)"                          \
			FREEFLAG=$(INTEL_FREEFLAG)                                \
			FIXEDFLAG=$(INTEL_FIXEDFLAG)                              \
			LINKFLAGS=$(LINK_FLAGS_ESMF)                              \
			LINKLIBS=$(LINK_LIBS_ESMF) NEMS="nems"                    \
			FLUSH_DEF=$(INTEL_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: lahey
lahey:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			PREP_FFLAGS=$(LAHEY_PREP_DEBUG)                           \
			FFLAGS_NO_DEBUG="$(LAHEY_FFLAGS_NO_DEBUG)"                \
			SP_FFLAGS="$(LAHEY_SP_FFLAGS)"                            \
			BYTE_SWAP_FLAG="$(LAHEY_BYTE_SWAP_FLAG)"                  \
			STATIC_FLAG="$(LAHEY_STATIC_FLAG)"                        \
			FIM_FFLAGS=$(LAHEY_FDEBUG)                                \
			POST_FFLAGS=$(LAHEY_POST_DEBUG)                           \
			CC=$(CClahey) CFLAGS=$(CDEBUG)                            \
			GRIBIO_CPP_FLAGS=$(LAHEY_GRIBIO_CPP_FLAGS)                \
			FC=$(FC) LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)               \
			FCginfo=$(FC)                                             \
			FCserial=$(FCseriallahey)                                 \
			LIBNETCDF=$(LIBNETCDF) INCNETCDF=$(INCNETCDF)             \
			COL_FFLAGS=$(LAHEY_COL_FDEBUG)                            \
			FREEFLAG=$(LAHEY_FREEFLAG)                                \
			FIXEDFLAG=$(LAHEY_FIXEDFLAG)                              \
			LINKFLAGS=$(LINK_FLAGS_NOESMF)                            \
			LINKLIBS=$(LINK_LIBS_NOESMF) NEMS=""                      \
			FLUSH_DEF=$(LAHEY_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: ranger
ranger:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			PREP_FFLAGS=$(PREP_FFLAGS)                                \
			FFLAGS_NO_DEBUG="$(PREP_FFLAGS)"                          \
			SP_FFLAGS="$(INTEL_SP_FFLAGS)"                            \
			BYTE_SWAP_FLAG="$(INTEL_BYTE_SWAP_FLAG)"                  \
			STATIC_FLAG="$(INTEL_STATIC_FLAG)"                        \
			FIM_FFLAGS="$(INTEL_FFLAGS)"                              \
			POST_FFLAGS=$(POST_FFLAGS)                                \
			CC=$(CCintel) CFLAGS=$(CFLAGS)                            \
			GRIBIO_CPP_FLAGS=$(INTEL_GRIBIO_CPP_FLAGS)                \
			FC=$(FC) LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)               \
			FCginfo=$(FC)                                             \
			FCserial=$(FCserialintel)                                 \
			LIBNETCDF="-L$(TACC_NETCDF_DIR)/lib -lnetcdf"             \
			INCNETCDF="-I$(TACC_NETCDF_DIR)/include"                  \
			COL_FFLAGS="$(INTEL_COL_FFLAGS)"                          \
			FREEFLAG=$(INTEL_FREEFLAG)                                \
			FIXEDFLAG=$(INTEL_FIXEDFLAG)                              \
			LINKFLAGS=$(LINK_FLAGS_NOESMF)                            \
			LINKLIBS=$(LINK_LIBS_NOESMF) NEMS=""                      \
			FLUSH_DEF=$(INTEL_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: devccs
devccs: vapor

.PHONY: vapor
vapor:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			ICOSIO_FFLAGS="$(IBM_ICOSIO_FFLAGS)"                      \
			PREP_FFLAGS="$(IBM_PREP_FFLAGS)"                          \
			FFLAGS_NO_DEBUG="$(IBM_PREP_FFLAGS)"                      \
			SP_FFLAGS="$(IBM_SP_FFLAGS)"                              \
			BYTE_SWAP_FLAG="$(IBM_BYTE_SWAP_FLAG)"                    \
			STATIC_FLAG="$(IBM_STATIC_FLAG)"                          \
			FIM_FFLAGS="$(IBM_FIM_FFLAGS)"                            \
			POST_FFLAGS="$(IBM_POST_FFLAGS)"                          \
			GRIBIO_CPP_FLAGS=$(IBM_GRIBIO_CPP_FLAGS)                  \
			CC=$(CCibm) CFLAGS="$(IBM_CFLAGS)"                        \
			FC=$(FCvapor) LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)          \
			FCginfo=$(FCserialvapor)                                  \
			FCserial=$(FCserialvapor)                                 \
			LIBNETCDF=$(IBM_LIBNETCDFvapor)                           \
			INCNETCDF=$(IBM_INCNETCDFvapor)                           \
			COL_FFLAGS="$(IBM_COL_FFLAGS)"                            \
			FREEFLAG=$(IBM_FREEFLAG)                                  \
			FIXEDFLAG=$(IBM_FIXEDFLAG)                                \
			LINKFLAGS=$(LINK_FLAGS_NOESMF)                            \
			LINKLIBS=$(LINK_LIBS_NOESMF) NEMS=""                      \
			FLUSH_DEF=$(IBM_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: bluefire
bluefire:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			ICOSIO_FFLAGS="$(IBM_ICOSIO_FFLAGS)"                      \
			PREP_FFLAGS="$(IBM_PREP_FFLAGS)"                          \
			FFLAGS_NO_DEBUG="$(IBM_PREP_FFLAGS)"                      \
			SP_FFLAGS="$(IBM_SP_FFLAGS)"                              \
			BYTE_SWAP_FLAG="$(IBM_BYTE_SWAP_FLAG)"                    \
			STATIC_FLAG="$(IBM_STATIC_FLAG)"                          \
			FIM_FFLAGS="$(IBM_FIM_FFLAGS)"                            \
			POST_FFLAGS="$(IBM_POST_FFLAGS)"                          \
			GRIBIO_CPP_FLAGS=$(IBM_GRIBIO_CPP_FLAGS)                  \
			CC=$(CCibm) CFLAGS="$(IBM_CFLAGS)"                        \
			FC=$(FCbluefire) LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)       \
			FCginfo=$(FCserialbluefire)                               \
			FCserial=$(FCserialbluefire)                              \
			LIBNETCDF=$(LIBNETCDF) INCNETCDF=$(INCNETCDF)             \
			COL_FFLAGS="$(IBM_COL_FFLAGS)"                            \
			FREEFLAG=$(IBM_FREEFLAG)                                  \
			FIXEDFLAG=$(IBM_FIXEDFLAG)                                \
			LINKFLAGS=$(LINK_FLAGS_NOESMF)                            \
			LINKLIBS=$(LINK_LIBS_NOESMF) NEMS=""                      \
			FLUSH_DEF=$(IBM_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: linuxpcgnu
linuxpcgnu:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			PREP_FFLAGS=$(GF_PREP_FFLAGS)                             \
			FFLAGS_NO_DEBUG="$(PREP_FFLAGS)"                          \
			SP_FFLAGS="$(GF_SP_FFLAGS)"                               \
			BYTE_SWAP_FLAG=""                                         \
			STATIC_FLAG=""                                            \
			FIM_FFLAGS="$(GF_FFLAGS)"                                 \
			POST_FFLAGS=$(GF_POST_FFLAGS)                             \
			CC=$(GF_CC) CFLAGS=$(GF_CFLAGS)                           \
			GRIBIO_CPP_FLAGS=$(INTEL_GRIBIO_CPP_FLAGS)                \
			FC=$(FC) LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)               \
			FCginfo=$(FC)                                             \
			FCserial=$(GF_FCserial) LIBNETCDF=$(LIBNETCDF)            \
			LIBNETCDF=$(LIBNETCDF) INCNETCDF=$(INCNETCDF)             \
			COL_FFLAGS="$(GF_COL_FFLAGS)"                             \
			FREEFLAG=$(GF_FREEFLAG)                                   \
			FIXEDFLAG=$(GF_FIXEDFLAG)                                 \
			LINKFLAGS=$(LINK_FLAGS_NOESMF)                            \
			LINKLIBS=$(LINK_LIBS_NOESMF) NEMS=""                      \
			FLUSH_DEF=$(INTEL_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: macgnu
macgnu:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			PREP_FFLAGS=$(GF_PREP_FFLAGS)                             \
			FFLAGS_NO_DEBUG="$(PREP_FFLAGS)"                          \
			SP_FFLAGS="$(GF_SP_FFLAGS)"                               \
			BYTE_SWAP_FLAG=""                                         \
			STATIC_FLAG=""                                            \
			FIM_FFLAGS="$(GF_FFLAGS)"                                 \
			POST_FFLAGS=$(GF_POST_FFLAGS)                             \
			CC=$(GF_CC) CFLAGS=$(GF_CFLAGS)                           \
			GRIBIO_CPP_FLAGS=$(INTEL_GRIBIO_CPP_FLAGS)                \
			FC=$(FC) LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)               \
			FCginfo=$(FC)                                             \
			FCserial=$(MAC_FCserial) LIBNETCDF=$(LIBNETCDF)           \
			LIBNETCDF=$(LIBNETCDF) INCNETCDF=$(INCNETCDF)             \
			COL_FFLAGS="$(GF_COL_FFLAGS)"                             \
			FREEFLAG=$(GF_FREEFLAG)                                   \
			FIXEDFLAG=$(GF_FIXEDFLAG)                                 \
			LINKFLAGS=$(LINK_FLAGS_NOESMF)                            \
			LINKLIBS=$(LINK_LIBS_NOESMF) NEMS=""                      \
			FLUSH_DEF=$(INTEL_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: jaguarintel
jaguarintel:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			PREP_FFLAGS=$(PREP_FFLAGS)                                \
			FFLAGS_NO_DEBUG="$(PREP_FFLAGS)"                          \
			SP_FFLAGS="$(INTEL_SP_FFLAGS)"                            \
			BYTE_SWAP_FLAG="$(INTEL_BYTE_SWAP_FLAG)"                  \
			STATIC_FLAG="$(INTEL_STATIC_FLAG)"                        \
			FIM_FFLAGS="$(INTEL_FFLAGS)"                              \
			POST_FFLAGS=$(POST_FFLAGS)                                \
			CC=cc CFLAGS=$(CFLAGS)                                    \
			GRIBIO_CPP_FLAGS=$(INTEL_GRIBIO_CPP_FLAGS)                \
			FC=ftn LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)                 \
			FCginfo=ftn                                               \
			FCserial=ftn                                              \
			LIBNETCDF=""                                              \
			INCNETCDF=""                                              \
			COL_FFLAGS="$(INTEL_COL_FFLAGS)"                          \
			FREEFLAG=$(INTEL_FREEFLAG)                                \
			FIXEDFLAG=$(INTEL_FIXEDFLAG)                              \
			LINKFLAGS=$(JAGUAR_LINK_FLAGS_NOESMF)                     \
			LINKLIBS=$(JAGUAR_LINK_LIBS_NOESMF) NEMS=""               \
			FLUSH_DEF=$(INTEL_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: jaguargnu
jaguargnu:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			PREP_FFLAGS=$(GF_PREP_FFLAGS)                             \
			FFLAGS_NO_DEBUG="$(PREP_FFLAGS)"                          \
			SP_FFLAGS="$(GF_SP_FFLAGS)"                               \
			BYTE_SWAP_FLAG=""                                         \
			STATIC_FLAG=""                                            \
			FIM_FFLAGS="$(GF_FFLAGS)"                                 \
			POST_FFLAGS=$(GF_POST_FFLAGS)                             \
			CC=cc CFLAGS=$(GF_CFLAGS)                                 \
			GRIBIO_CPP_FLAGS=$(INTEL_GRIBIO_CPP_FLAGS)                \
			FC=ftn LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)                 \
			FCginfo=ftn                                               \
			FCserial=ftn                                              \
			LIBNETCDF=""                                              \
			INCNETCDF=""                                              \
			COL_FFLAGS="$(GF_COL_FFLAGS)"                             \
			FREEFLAG=$(GF_FREEFLAG)                                   \
			FIXEDFLAG=$(GF_FIXEDFLAG)                                 \
			LINKFLAGS=$(JAGUAR_LINK_FLAGS_NOESMF)                     \
			LINKLIBS=$(JAGUAR_LINK_LIBS_NOESMF) NEMS=""               \
			FLUSH_DEF=$(INTEL_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: frostintel
frostintel:
	$(MAKE) -f Makesub MAKE=$(MAKE)                               \
			PREP_FFLAGS=$(PREP_FFLAGS)                                \
			FFLAGS_NO_DEBUG="$(PREP_FFLAGS)"                          \
			SP_FFLAGS="$(INTEL_SP_FFLAGS)"                            \
			BYTE_SWAP_FLAG="$(INTEL_BYTE_SWAP_FLAG)"                  \
			STATIC_FLAG="$(INTEL_STATIC_FLAG)"                        \
			FIM_FFLAGS="$(INTEL_FFLAGS)"                              \
			POST_FFLAGS=$(POST_FFLAGS)                                \
			CC=cc CFLAGS=$(CFLAGS)                                    \
			GRIBIO_CPP_FLAGS=$(INTEL_GRIBIO_CPP_FLAGS)                \
			FC=ifort LIBDIR=$(LIBDIR)  BINDIR=$(BINDIR)               \
			FCginfo=ifort                                             \
			FCserial=ifort                                            \
			LIBNETCDF="-L/ccs/home/rosinski/frost/intel/lib -lnetcdf" \
			INCNETCDF="-I/ccs/home/rosinski/frost/intel/include"      \
			COL_FFLAGS="$(INTEL_COL_FFLAGS)"                          \
			FREEFLAG=$(INTEL_FREEFLAG)                                \
			FIXEDFLAG=$(INTEL_FIXEDFLAG)                              \
			LINKFLAGS=$(FROST_LINK_FLAGS_NOESMF)                      \
			LINKLIBS=$(FROST_LINK_LIBS_NOESMF) NEMS=""                \
			FLUSH_DEF=$(INTEL_FLUSH_DEF) MKDEPENDS=$(MKDEPENDS)

.PHONY: cleanall
cleanall: clean
.PHONY: clean
clean:
	$(RM) -rf ../FIMsrc_* macros.make

# Do not use this unless you know what you are doing.
.PHONY: cleanbelow
cleanbelow:
	(cd bin    && rm -f *)
	(cd lib    && rm -f *)
	(cd icosio && $(MAKE) clean)
	(cd cntl   && $(MAKE) clean)
	(cd bacio  && $(MAKE) clean)
	(cd w3     && $(MAKE) clean)
	(cd utils  && $(MAKE) clean)
	(cd prep   && $(MAKE) clean)
	(cd fim    && $(MAKE) clean)
	(cd post   && $(MAKE) clean)
