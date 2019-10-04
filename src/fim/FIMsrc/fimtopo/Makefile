# fimtopo Makefile

include ../macros.make

#FC      = g95
#LD      = g95
CMD      = fimtopo.x
CPPFLAGS =
FC       = /opt/lahey/lf6480/bin/f95
FC       = ifort
FFLAGS   =
LD       = /opt/lahey/lf6480/bin/f95
LD       = ifort
LDFLAGS  = -g  
LIBS     =	
OBJS     = fimtopo.o
SHELL    = /bin/sh

# Lines from here on down should not need to be changed. They are the actual
# rules which make uses to build the executable.

$(CMD): $(OBJS)
	$(LD) -o $(CMD) $(OBJS) $(LIBS) $(LDFLAGS)

clean:
	rm -f *.o $(CMD) *.mod

%.o: %.f90
	@echo Compiling $<; $(FC) -o $@ -c $< 







