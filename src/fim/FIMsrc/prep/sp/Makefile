#!/bin/sh
###############################################################
#
#   AUTHOR:    Vuong - W/NP11
#
#   DATE:      12/04/2000
#
#   PURPOSE:   This script uses the make utility to update the libsp 
#              archive libraries.
#              It first reads a list of source files in the library and
#              then generates a makefile used to update the archive
#              libraries.  The make command is then executed for each
#              archive library, where the archive library name and 
#              compilation flags are passed to the makefile through 
#              environment variables.
#
#   REMARKS:   Only source files that have been modified since the last
#              library update are recompiled and replaced in the object
#              archive libraries.  The make utility determines this
#              from the file modification times.
#
#              New source files are also compiled and added to the object 
#              archive libraries.
#
###############################################################

#
#     Generate a list of object files that corresponds to the
#     list of Fortran ( .f ) files in the current directory
#

include ../../macros.make

SRCS = $(shell ls *.f)
OBJS = $(addsuffix .o, $(basename $(SRCS)))
LIB = $(LIBDIR)/libsp_4.a

.SUFFIXES:
.SUFFIXES: .o .f

.f.o:
	$(FC) -c $(FFLAGS) $(FIXEDFLAG) $<

all: $(LIB)

$(LIB): $(OBJS)
	ar ruv $(AFLAGS)  $@ $(OBJS)

clean:
	$(RM) *.o *.mod
