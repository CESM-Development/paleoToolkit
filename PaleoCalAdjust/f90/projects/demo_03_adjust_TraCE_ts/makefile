# demo_03_adjust_TraCE_ts

# (relative) paths to source files -- needs localization
prog = ../../main_programs/
mod = ../../modules/

# compiler and flags
FC = gfortran
FFLAGS = -O3 -Wall

# target, object and .mod files
target = demo_03_adjust_TraCE_ts
objects = pseudo_daily_interp_subs.o calendar_effects_subs.o demo_03_adjust_TraCE_ts.o
modules = pseudo_daily_interp_subs.mod calendar_effects_subs.mod 

$(target): $(objects) 
	$(FC) $(FFLAGS) -o $(target) $(objects)

pseudo_daily_interp_subs.o: $(mod)pseudo_daily_interp_subs.f90
	$(FC) $(FFLAGS) -c $(mod)pseudo_daily_interp_subs.f90

calendar_effects_subs.o: $(mod)calendar_effects_subs.f90
	$(FC) $(FFLAGS) -c $(mod)calendar_effects_subs.f90

demo_03_adjust_TraCE_ts.o: $(prog)demo_03_adjust_TraCE_ts.f90
	$(FC) $(FFLAGS) -c $(prog)demo_03_adjust_TraCE_ts.f90

clean: 
	rm -f $(objects)
	rm -f $(modules)

clobber: clean
	rm -f $(target)	
