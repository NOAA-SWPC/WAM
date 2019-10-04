include ../../../macros.make

CMD=$(BINDIR)/global_sighdr
INCMOD=../../sigio
LIBS=-L$(LIBDIR) -lsigio_4 -lw3_4
SHELL=/bin/sh
SRCM=sighdr.f

$(CMD): $(SRCM)
	$(FC) $(FREEFLAG) -I $(INCMOD) $(FFLAGS) $(LDFLAGS) $(SRCM) $(LIBS) -o $(CMD)

clean:
	$(RM) $(CMD)
