# Makefile for Fortran programs
# It is assumed that the HDF5 library is available, and includes h5fc
# as a wrapper for the Fortran compiler

FC = h5fc
FCFLAGS = -O2 -Wall -fdefault-real-8
LDFLAGS =

PROGRAMS = md_springs md_springs_mts md_constraints

all: $(PROGRAMS)

hdf5_module.o: hdf5_module.f90
hdf5_module.mod: hdf5_module.f90
maths_module.o: maths_module.f90
maths_module.mod: maths_module.f90
config_io_module.o: config_io_module.f90
config_io_module.mod: config_io_module.f90
md_module.o: md_module.f90 maths_module.mod
md_module.mod: md_module.f90 maths_module.mod
force_module.o: force_module.f90
force_module.mod: force_module.f90
md_springs.o: md_springs.f90 maths_module.mod md_module.mod force_module.mod config_io_module.mod hdf5_module.mod
md_springs_mts.o: md_springs_mts.f90 maths_module.mod md_module.mod force_module.mod config_io_module.mod hdf5_module.mod
md_constraints.o: md_constraints.f90 maths_module.mod md_module.mod force_module.mod config_io_module.mod hdf5_module.mod
md_springs: md_module.o force_module.o config_io_module.o maths_module.o hdf5_module.o
md_springs_mts: md_module.o force_module.o config_io_module.o maths_module.o hdf5_module.o
md_constraints: md_module.o force_module.o config_io_module.o maths_module.o hdf5_module.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o %.mod: %.f90
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o core *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
