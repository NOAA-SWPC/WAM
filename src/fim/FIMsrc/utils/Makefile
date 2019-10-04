# utils Makefile

include ../macros.make

SHELL = /bin/sh

COMMONOBJS       = read_queue_namelist.o module_initial_chem_namelists.o 
GETCHEM          = $(BINDIR)/GetChemEnabled
GETCOMPUTETASKS  = $(BINDIR)/GetComputeTasks
GETDATADIR       = $(BINDIR)/GetDATADIR
GETDATADR2       = $(BINDIR)/GetDATADR2
GETFIMDIR        = $(BINDIR)/GetFIMDIR
GETGLVL          = $(BINDIR)/GetGLVL
GETNIP           = $(BINDIR)/GetNIP 
GETNVL           = $(BINDIR)/GetNVL
GETPARALLELISM   = $(BINDIR)/GetParallelism
GETPREPDIR       = $(BINDIR)/GetPREPDIR
GETQUEUETIME     = $(BINDIR)/GetQueueTime
GETSRCDIR        = $(BINDIR)/GetSRCDIR
GETWRFCU         = $(BINDIR)/GetWRFcuEnabled
GETWRFMP         = $(BINDIR)/GetWRFmpEnabled
GETWRITETASKINFO = $(BINDIR)/GetWriteTaskInfo
GET_NUM_CORES    = $(BINDIR)/get_num_cores
LIBSYSSHARE      = $(LIBDIR)/libsysshare.a
REDUCE           = $(BINDIR)/reduce
EXTRACT_ATCF     = $(BINDIR)/extract_atcf
TARGETS = $(GETCOMPUTETASKS) $(GETDATADIR) $(GETDATADR2) $(GETFIMDIR) $(GETGLVL) $(GETNIP) $(GETNVL) \
          $(GETPARALLELISM) $(GETPREPDIR) $(GETQUEUETIME) $(GETSRCDIR) $(GETWRITETASKINFO) $(REDUCE) \
          $(GETCHEM) $(GETWRFCU) $(GETWRFMP) $(GET_NUM_CORES) $(EXTRACT_ATCF) headers.o

.SUFFIXES:
.SUFFIXES: .o .f .F90

.f.o:
	$(FC) -c $(DEBUG_FLAGS) $(FFLAGS) $(FREEFLAG) $<
.F90.o:
	$(FC) -c $(DEBUG_FLAGS) $(FFLAGS) $(FREEFLAG) $<

all: $(TARGETS)

$(GETWRITETASKINFO): DEPENDENCIES GetWriteTaskInfo.o wtinfo.o
	$(FC) $(FFLAGS) -o $@ GetWriteTaskInfo.o wtinfo.o $(LIBSYSSHARE)
$(GETCHEM): DEPENDENCIES GetChemEnabled.o read_queue_namelist.o
	$(FC) $(FFLAGS) -o $@ GetChemEnabled.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETDATADIR): DEPENDENCIES GetDATADIR.o read_queue_namelist.o
	$(FC) $(FFLAGS) -o $@ GetDATADIR.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETDATADR2): DEPENDENCIES GetDATADR2.o read_queue_namelist.o
	$(FC) $(FFLAGS) -o $@ GetDATADR2.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETFIMDIR): DEPENDENCIES GetFIMDIR.o read_queue_namelist.o
	$(FC) $(FFLAGS) -o $@ GetFIMDIR.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETPREPDIR): DEPENDENCIES GetPREPDIR.o $(COMMONOBJS)
	$(FC) $(FFLAGS) -o $@ GetPREPDIR.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETSRCDIR): DEPENDENCIES GetSRCDIR.o read_queue_namelist.o
	$(FC) $(FFLAGS) -o $@ GetSRCDIR.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETCOMPUTETASKS): DEPENDENCIES GetComputeTasks.o read_queue_namelist.o
	$(FC) $(FFLAGS) $(STATIC_FLAG) -o $@ GetComputeTasks.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETQUEUETIME): DEPENDENCIES GetQueueTime.o read_queue_namelist.o
	$(FC) $(FFLAGS) $(STATIC_FLAG) -o $@ GetQueueTime.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETPARALLELISM): DEPENDENCIES GetParallelism.o read_queue_namelist.o
	$(FC) $(FFLAGS) $(STATIC_FLAG) -o $@ GetParallelism.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETGLVL): DEPENDENCIES GetGLVL.o read_queue_namelist.o
	$(FC) $(FFLAGS) -o $@ GetGLVL.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETNIP): DEPENDENCIES GetNIP.o read_queue_namelist.o
	$(FC) $(FFLAGS) -o $@ GetNIP.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETNVL): DEPENDENCIES GetNVL.o read_queue_namelist.o
	$(FC) $(FFLAGS) -o $@ GetNVL.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETWRFCU): DEPENDENCIES GetWRFcuEnabled.o read_queue_namelist.o
	$(FC) $(FFLAGS) -o $@ GetWRFcuEnabled.o $(COMMONOBJS) $(LIBSYSSHARE)
$(GETWRFMP): DEPENDENCIES GetWRFmpEnabled.o read_queue_namelist.o
	$(FC) $(FFLAGS) -o $@ GetWRFmpEnabled.o $(COMMONOBJS) $(LIBSYSSHARE)
$(REDUCE): DEPENDENCIES reduce.o
	$(FC) $(FFLAGS) $(STATIC_FLAG) -o $@ reduce.o $(LIBSYSSHARE)
$(GET_NUM_CORES): DEPENDENCIES get_num_cores.o read_queue_namelist.o wtinfo.o
	$(FC) $(FFLAGS) -o $@ get_num_cores.o $(COMMONOBJS) wtinfo.o $(LIBSYSSHARE)
$(EXTRACT_ATCF): DEPENDENCIES extract_atcf.o
	$(FC) $(FFLAGS) $(STATIC_FLAG) -o $@ extract_atcf.o

DEPENDENCIES:
	$(RM) -f Filepath Srcfiles
	echo "." > Filepath
	ls -1 *.F90 > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

-include DEPENDENCIES

clean:
	$(RM) *.o *.mod DEPENDENCIES
