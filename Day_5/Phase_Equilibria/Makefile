# Makefile for Fortran programs
# It is assumed that the HDF5 library is available, and includes h5fc
# as a wrapper for the Fortran compiler

FC = h5fc
FCFLAGS = -O2 -Wall -fdefault-real-8
LDFLAGS =

PROGRAMS = mc_muca mc_gibbs

all: $(PROGRAMS)

maths_module.o: maths_module.f90
maths_module.mod: maths_module.f90
hdf5_module.o: hdf5_module.f90
hdf5_module.mod: hdf5_module.f90
potential_module.o: potential_module.f90 maths_module.mod
potential_module.mod: potential_module.f90 maths_module.mod
config_io_module.o: config_io_module.f90
config_io_module.mod: config_io_module.f90
mc_module.o: mc_module.f90 potential_module.mod maths_module.mod
mc_module.mod: mc_module.f90 potential_module.mod maths_module.mod
mc_gibbs_module.o: mc_gibbs_module.f90 potential_module.mod maths_module.mod
mc_gibbs_module.mod: mc_gibbs_module.f90 potential_module.mod maths_module.mod
mc_muca_module.o: mc_muca_module.f90 potential_module.mod maths_module.mod
mc_muca_module.mod: mc_muca_module.f90 potential_module.mod maths_module.mod
mc_gibbs.o: mc_gibbs.f90 mc_module.mod mc_gibbs_module.mod potential_module.mod config_io_module.mod hdf5_module.mod
mc_muca.o: mc_muca.f90 mc_module.mod mc_muca_module.mod potential_module.mod config_io_module.mod hdf5_module.mod
mc_gibbs: mc_module.o mc_gibbs_module.o potential_module.o config_io_module.o maths_module.o hdf5_module.o
mc_muca: mc_module.o mc_muca_module.o potential_module.o config_io_module.o maths_module.o hdf5_module.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o %.mod: %.f90
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o core *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)





