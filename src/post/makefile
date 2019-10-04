include ../../conf/configure.nems

MAKEFILE = makefile

UTILINCSPOST = -I../share -I../nmm/ -I../io $(POST_INC)
UTILINCSPOSTGFS = -I../share -I../io $(POST_INC)
UTILINCSPOSTNMM = -I../share -I../nmm $(POST_INC)
UTILINCSPOSTSTUB = -I../share

LIBRARY  = libpost.a

MODULES      = post_gfs.o post_nmm.o post_nems_routines.o

MODULES_STUB  = post_gfs_stub.o post_nmm_stub.o post_nems_routines_stub.o
MODULES_GFS_STUB  = post_nmm.o post_gfs_stub.o post_nems_routines.o
MODULES_NMM_STUB  = post_nmm_stub.o post_gfs.o post_nems_routines.o

DEPEND_FILES = ${MODULES:.o=.F90}

all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(MODULES)
	$(AR) $(ARFLAGS) $@ $(MODULES)

stub: $(MODULES_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_STUB)
	
gfs_stub: $(MODULES_GFS_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_GFS_STUB)
	
nmm_stub: $(MODULES_NMM_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_NMM_STUB)
	
$(MODULES): %.o: %.f90
	$(FC) $(FFLAGS_GFS) $(TRAPS) $(UTILINCSPOST) -c $*.f90

$(MODULES_GFS_STUB): %.o: %.f90
	$(FC) $(FFLAGS) $(TRAPS) $(UTILINCSPOSTNMM) -c $*.f90

$(MODULES_NMM_STUB): %.o: %.f90
	$(FC) $(FFLAGS_GFS) $(TRAPS) $(UTILINCSPOSTGFS) -c $*.f90

$(MODULES_STUB): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCSPOSTSTUB) -c $*.f90

clean:
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst lm map depend

MKDEPENDS = ../../../exe/mkDepends.pl

include ../../conf/make.rules

include depend
