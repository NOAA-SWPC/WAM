# sys_share Makefile

SHELL = /bin/sh

include ../macros.make

FLAGS       = $(FFLAGS) $(FLUSH_DEF)
LIBSYSSHARE = $(LIBDIR)/libsysshare.a
OBJ         = sys_share.o

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FLAGS) $<

all: $(LIBSYSSHARE)

$(LIBSYSSHARE): $(LIBSYSSHARE)($(OBJ))
	$(RM) *.o

clean:
	$(RM) *.o *.mod
