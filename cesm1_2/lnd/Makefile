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

TARGETS  = paleo_mkraw_cesm1

all: $(TARGETS)

paleo_mkraw_cesm1 : paleo_mkraw_cesm1.F90
	$(FC) -o paleo_mkraw_cesm1 paleo_mkraw_cesm1.F90 $(INC_NETCDF) $(LIB_NETCDF)

clean:
	-rm $(TARGETS)
	










