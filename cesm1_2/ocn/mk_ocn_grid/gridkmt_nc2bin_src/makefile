#
# Makefile for gridkmt_nc2bin
#
#===============================================================================


COMP = ifort -assume byterecl -convert big_endian
# COMP = ifort -assume byterecl -convert little_endian
FLAGS = -I/usr/local/include
# FLAGS = -g -qrealsize=8 -qfullpath -O2 -qstrict -qmaxmem=-1 -qarch=auto -Q -qsigtrap=xl__trcedump -I/usr/local/include
LIB  = -L /glade/apps/opt/netcdf/4.2/intel/default/lib -lnetcdf
INCL = /glade/u/apps/ch/opt/netcdf/4.7.1/intel/19.0.5/include

SRCDIR  = .
# SRCDIR  = /glade/p/cgd/ccr/people/nanr/permian/mapping/makegrids
EXEDIR  = .

gridkmt_nc2bin: gridkmt_nc2bin.o
	$(COMP) $(FLAGS) gridkmt_nc2bin.o $(LIB) -o $(EXEDIR)/gridkmt_nc2bin

gridkmt_nc2bin.o: $(SRCDIR)/gridkmt_nc2bin.f90 $(INCL)
	$(COMP) $(FLAGS) -c $(SRCDIR)/gridkmt_nc2bin.f90
clean: 
	/bin/rm *.o *.mod

