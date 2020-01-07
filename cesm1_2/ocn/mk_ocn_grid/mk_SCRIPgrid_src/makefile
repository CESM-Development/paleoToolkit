#
# Makefile for interpolation code
# CVS:$Id: makefile,v 1.3 2000/04/19 22:05:59 pwjones Exp $
#
#===============================================================================
# SVN $Id$
# SVN $URL$
#===============================================================================

COMP = ifort -assume byterecl -convert big_endian
# COMP = ifort -assume byterecl -convert little_endian
FLAGS = -I/usr/local/include
# FLAGS = -g -qrealsize=8 -qfullpath -O2 -qstrict -qmaxmem=-1 -qarch=auto -Q -qsigtrap=xl__trcedump -I/usr/local/include
LIB  = -L /glade/apps/opt/netcdf/4.2/intel/default/lib -lnetcdf
INCL = /glade/u/apps/ch/opt/netcdf/4.7.1/intel/19.0.5/include

SRCDIR  = ./
# SRCDIR  = /glade/p/cgd/ccr/people/nanr/permian/mapping/makegrids
EXEDIR  = .
OBJ  = \
	kinds_mod.o \
	constants.o \
	iounits.o \
	netcdf.o

all: myconvertPOPTkmt myconvertPOPT convertPOPT convertgauss create_latlon

convertold: $(OBJ) convert_old.o
	$(COMP) $(FLAGS) $(OBJ) convert_old.o $(LIB) -o $(EXEDIR)/convertold

myconvertPOPT: $(OBJ) myconvertPOPT.o
	$(COMP) $(FLAGS) $(OBJ) myconvertPOPT.o $(LIB) -o $(EXEDIR)/myconvertPOPT

myconvertPOPTkmt: $(OBJ) myconvertPOPTkmt.o
	$(COMP) $(FLAGS) $(OBJ) myconvertPOPTkmt.o $(LIB) -o $(EXEDIR)/myconvertPOPTkmt

convertPOPT: $(OBJ) convertPOPT.o
	$(COMP) $(FLAGS) $(OBJ) convertPOPT.o $(LIB) -o $(EXEDIR)/convertPOPT

convertgauss: $(OBJ) convertgauss.o
	$(COMP) $(FLAGS) $(OBJ) convertgauss.o $(LIB) -o $(EXEDIR)/convertgauss

create_latlon: $(OBJ) create_latlon.o
	$(COMP) $(FLAGS) $(OBJ) create_latlon.o $(LIB) -o $(EXEDIR)/create_latlon

kinds_mod.o: $(SRCDIR)/kinds_mod.f $(INCL)
	$(COMP) $(FLAGS) -c $(SRCDIR)/kinds_mod.f

constants.o: $(SRCDIR)/constants.f kinds_mod.o $(INCL)
	$(COMP) $(FLAGS) -c $(SRCDIR)/constants.f

iounits.o: $(SRCDIR)/iounits.f kinds_mod.o constants.o $(INCL)
	$(COMP) $(FLAGS) -c $(SRCDIR)/iounits.f

netcdf.o: $(SRCDIR)/netcdf.f kinds_mod.o constants.o $(INCL)
	$(COMP) $(FLAGS) -c $(SRCDIR)/netcdf.f

convert_old.o: convert_old.f kinds_mod.o constants.o iounits.o netcdf.o $(INCL)
	$(COMP) $(FLAGS) -c convert_old.f

myconvertPOPTkmt.o: myconvertPOPTkmt.f kinds_mod.o constants.o iounits.o netcdf.o $(INCL)
	$(COMP) $(FLAGS) -c myconvertPOPTkmt.f

myconvertPOPT.o: myconvertPOPT.f kinds_mod.o constants.o iounits.o netcdf.o $(INCL)
	$(COMP) $(FLAGS) -c myconvertPOPT.f

convertPOPT.o: convertPOPT.f kinds_mod.o constants.o iounits.o netcdf.o $(INCL)
	$(COMP) $(FLAGS) -c convertPOPT.f

convertgauss.o: convertgauss.f kinds_mod.o constants.o iounits.o netcdf.o $(INCL)
	$(COMP) $(FLAGS) -c convertgauss.f

create_latlon.o: create_latlon.f kinds_mod.o constants.o iounits.o netcdf.o $(INCL)
	$(COMP) $(FLAGS) -c create_latlon.f

clean: 
	/bin/rm *.o *.mod

