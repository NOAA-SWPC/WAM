# cntl Makefile

SHELL = /bin/sh

include ../macros.make

FLAGS   = $(FFLAGS) $(INCS)
INCS    = -I $(UTILDIR)
LIBCNTL = $(LIBDIR)/libcntl.a
OBJS    = module_control.o units.o module_constants.o module_wrf_control.o
UTILDIR = ../utils

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FLAGS) $<

all: $(LIBCNTL)

$(LIBCNTL): $(OBJS)
	$(AR) ruv $@ $(OBJS)

module_control.o: module_wrf_control.o units.o

clean:
	$(RM) -r *.o *.mod
