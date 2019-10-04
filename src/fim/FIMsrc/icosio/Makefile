# icosio Makefile

# DEFS may include: -DSERIAL for no-MPI build
#                   -DNOGRIB for no-GRIB build

all: icosio.o

icosio.o: icosio_cpp.F90
	$(FC) -c $(FFLAGS) $(GPTL_FFLAGS) icosio_cpp.F90 -o icosio.o

icosio_cpp.F90: icosio.F90
	$(CPP) $(CPP_FLAGS) $(DEFS) icosio.F90 > icosio_cpp.F90

clean:
	$(RM) icosio_cpp.F90 icosio.o *.mod
