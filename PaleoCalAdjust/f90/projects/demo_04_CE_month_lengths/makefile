# demo_04_CE_month_lengths

# (relative) paths to source files -- needs localization
prog = ../../main_programs/
mod = ../../modules/

# compiler and flags
FC = gfortran
FFLAGS = -O3 -Wall

# target, object and .mod files
target = demo_04_CE_month_lengths
objects = GISS_orbpar_subs.o GISS_srevents_subs.o month_length_subs.o demo_04_CE_month_lengths.o
modules = GISS_orbpar_subs.mod GISS_srevents_subs.mod month_length_subs.mod

$(target): $(objects) 
	$(FC) $(FFLAGS) -o $(target) $(objects)

GISS_orbpar_subs.o: $(mod)GISS_orbpar_subs.f90
	$(FC) $(FFLAGS) -c $(mod)GISS_orbpar_subs.f90

GISS_srevents_subs.o: $(mod)GISS_srevents_subs.f90
	$(FC) $(FFLAGS) -c $(mod)GISS_srevents_subs.f90

month_lengths_subs.o: $(mod)month_length_subs.f90
	$(FC) $(FFLAGS) -c $(mod)month_length_subs.f90

demo_04_CE_month_lengths.o: $(prog)demo_04_CE_month_lengths.f90
	$(FC) $(FFLAGS) -c $(prog)demo_04_CE_month_lengths.f90

clean: 
	rm -f $(objects)
	rm -f $(modules)

clobber: clean
	rm -f $(target)