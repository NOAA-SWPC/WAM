include ../../../macros.make

CMD=$(BINDIR)/global_sfchdr
INCMOD=../../sfcio
LIBS=-L$(LIBDIR) -lsfcio_4 -lw3_4
SHELL=/bin/sh
SRCM=sfchdr.f

$(CMD): $(SRCM)
	$(FC) $(FREEFLAG) -I $(INCMOD) $(FFLAGS) $(LDFLAGS) $(SRCM) $(LIBS) -o $(CMD)

clean:
	$(RM) $(CMD)
