#-----------------------------------------------------------------------
# This Makefile is for building any of the executables in this directory.
# USAGE: eg. gmake EXENAME=paleo_mkraw_cesm1
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# NWSC
#------------------------------------------------------------------------

FC         = ifort -assume byterecl
LIB_NETCDF = -L/glade/apps/opt/netcdf/4.2/intel/default/lib -lnetcdf
INC_NETCDF = -I/glade/apps/opt/netcdf/4.2/intel/default/include

TARGETS  = topo2rdirc check_inf_loop

all: $(TARGETS)

topo2rdirc : topo2rdirc.F90
	$(FC) -o topo2rdirc topo2rdirc.F90 $(INC_NETCDF) $(LIB_NETCDF)

check_inf_loop : check_inf_loop.F90
	$(FC) -o check_inf_loop check_inf_loop.F90 $(INC_NETCDF) $(LIB_NETCDF)

clean:
	-rm $(TARGETS)
	










