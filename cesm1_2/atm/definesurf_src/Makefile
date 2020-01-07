# Makefile to build definesurf on various platforms
# Note: If netcdf library is not built in the standard location, you must set the environment
# variables INC_NETCDF and LIB_NETCDF

EXEDIR = .
EXENAME = definesurf
RM = rm

.SUFFIXES:
.SUFFIXES: .f90 .F90 .o

FC         := ifort 
LIB_NETCDF := -L/glade/apps/opt/netcdf/4.2/intel/default/lib -lnetcdf -lnetcdff
INC_NETCDF := -I/glade/apps/opt/netcdf/4.2/intel/default/include


.f90.o:
	$(FC) $(LIB_NETCDF) $(INC_NETCDF) -c $< -o $@ 

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := ao.o ao_i.o area_ave.o binf2c.o cell_area.o  \
         chkdims.o endrun.o fmain.o handle_error.o inimland.o \
         lininterp.o map_i.o max_ovr.o shr_kind_mod.o sghphis.o sm121.o \
         terrain_filter.o varf2c.o wrap_nf.o interplandm.o map2f.o

definesurf: $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIB_NETCDF) 

debug: 
	echo "OBS = " $(OBJS)

clean:
	$(RM) -f $(OBJS)  *.mod $(EXEDIR)/$(EXENAME)

ao.o : ao.f90 shr_kind_mod.o
ao_i.o : ao_i.f90 shr_kind_mod.o
area_ave.o : area_ave.f90 shr_kind_mod.o
binf2c.o : binf2c.f90 shr_kind_mod.o
cell_area.o : cell_area.f90 shr_kind_mod.o
chkdims.o : chkdims.f90
endrun.o : endrun.f90
fmain.o : fmain.f90 shr_kind_mod.o
handle_error.o : handle_error.f90
inimland.o : inimland.f90 shr_kind_mod.o
lininterp.o : lininterp.f90 shr_kind_mod.o
map_i.o : map_i.f90 shr_kind_mod.o
max_ovr.o : max_ovr.f90 shr_kind_mod.o
shr_kind_mod.o : shr_kind_mod.f90
sghphis.o : sghphis.f90 shr_kind_mod.o
sm121.o : sm121.f90 shr_kind_mod.o
terrain_filter.o : terrain_filter.f90
map2f.o : map2f.f90 
varf2c.o : varf2c.f90 shr_kind_mod.o
wrap_nf.o : wrap_nf.f90
interplandm.o : interplandm.f90

