###############################################################
#
#   AUTHOR:    Gilbert - W/NP11
#
#   DATE:      01/11/1999
#
#   DATE:      03/25/2010
#   Modified by Rosinski (NOAA/ESRL) to use only a Makefile,
#   eliminate the need for a script
#
###############################################################

#
#     Generate a list of object files that corresponds to the
#     list of Fortran ( .f, .F90 ) files in the current directory
#

include ../macros.make

LIB  = $(LIBDIR)/libw3_4.a
OBJS = $(addsuffix .o, $(basename $(SRCS)))
SRCS = $(shell ls *.f *.F90 | grep -v jdate.F90)

.SUFFIXES:
.SUFFIXES: .o .f .F90

.F90.o:
	$(FC) $(FFLAGS_NO_DEBUG) $(FREEFLAG) -c $<

.f.o:
	$(FC) -c $(FFLAGS) $(FIXEDFLAG) $<

all: $(LIB) $(BINDIR)/jdate

$(LIB): $(OBJS)
	ar ruv $(AFLAGS) $@ $(OBJS)

$(BINDIR)/jdate: $(LIB)
	$(FCserial) $(FFLAGS) $(FREEFLAG) -o $@ jdate.F90 $(LIB)

clean:
	$(RM) *.o *.mod
