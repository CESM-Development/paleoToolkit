# month_length

# (relative) paths to source files -- needs localization
prog = ../../main_programs/
mod = ../../modules/

# compiler and flags
FC = gfortran
FFLAGS = -O3 -Wall

# target, object and .mod files
target = month_length
objects = GISS_orbpar_subs.o GISS_srevents_subs.o month_length_subs.o month_length.o
modules = GISS_orbpar_subs.mod GISS_srevents_subs.mod month_length_subs.mod

$(target): $(objects) 
	$(FC) $(FFLAGS) -o $(target) $(objects)

GISS_orbpar_subs.o: $(mod)GISS_orbpar_subs.f90
	$(FC) $(FFLAGS) -c $(mod)GISS_orbpar_subs.f90

GISS_srevents_subs.o: $(mod)GISS_srevents_subs.f90
	$(FC) $(FFLAGS) -c $(mod)GISS_srevents_subs.f90

month_length_subs.o: $(mod)month_length_subs.f90
	$(FC) $(FFLAGS) -c $(mod)month_length_subs.f90

month_length.o: $(prog)month_length.f90
	$(FC) $(FFLAGS) -c $(prog)month_length.f90

clean: 
	rm -f $(objects)
	rm -f $(modules)

clobber: clean
	rm -f $(target)