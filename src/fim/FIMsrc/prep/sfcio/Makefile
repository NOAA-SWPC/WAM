# sfcio Makefile

SHELL = /bin/sh

include ../../macros.make

#"-O0 -g -traceback -C" does not work with -mcmodel=medium -i-dynamic - JFM 12/12/07
FLAGS    = $(FFLAGS_NO_DEBUG) $(BYTE_SWAP_FLAG) $(DEBUG_FLAGS)
LIBSFCIO = $(LIBDIR)/libsfcio_4.a

.SUFFIXES:
.SUFFIXES: .a .o .F90

all: $(LIBSFCIO)

$(LIBSFCIO): $(LIBSFCIO)(sfcio_module.o)

.F90.o: 
	$(FC) -c $(FLAGS) $<

clean:
	rm -f *.o sfcio*.mod
