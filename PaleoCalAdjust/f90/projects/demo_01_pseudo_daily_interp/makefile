# demo_01_pseudo_daily_interp

# (relative) paths to source files
prog = ../../main_programs/
mod = ../../modules/

# compiler and flags
FC = gfortran
FFLAGS = -O3 -Wall

# target, object and .mod files
target = demo_01_pseudo_daily_interp
objects = pseudo_daily_interp_subs.o demo_01_pseudo_daily_interp.o
modules = pseudo_daily_interp_subs.mod

$(target): $(objects) 
	$(FC) $(FFLAGS) -o $(target) $(objects)

pseudo_daily_interp_subs.o: $(mod)pseudo_daily_interp_subs.f90
	$(FC) $(FFLAGS) -c $(mod)pseudo_daily_interp_subs.f90

demo_01_pseudo_daily_interp.o: $(prog)demo_01_pseudo_daily_interp.f90
	$(FC) $(FFLAGS) -c $(prog)demo_01_pseudo_daily_interp.f90

clean: 
	rm -f $(objects)
	rm -f $(modules)

clobber: clean
	rm -f $(target)