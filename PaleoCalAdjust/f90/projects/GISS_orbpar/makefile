# GISS_orbpar

# (relative) paths to source files -- needs localization
prog = ../../main_programs/
mod = ../../modules/

# compiler and flags
FC = gfortran
FFLAGS = -O3 -Wall

# target, object and .mod files
target = GISS_orbpar
objects = GISS_orbpar_subs.o GISS_orbpar_driver.o
modules = GISS_orbpar_subs.mod

$(target): $(objects) 
	$(FC) $(FFLAGS) -o $(target) $(objects)

GISS_orbpar_subs.o: $(mod)GISS_orbpar_subs.f90
	$(FC) $(FFLAGS) -c $(mod)GISS_orbpar_subs.f90

GISS_orbpar_driver.o: $(prog)GISS_orbpar_driver.f90
	$(FC) $(FFLAGS) -c $(prog)GISS_orbpar_driver.f90

clean: 
	rm -f $(objects)
	rm -f $(modules)

clobber: clean
	rm -f $(target)