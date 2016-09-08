#!/bin/csh -f
set echo
set verbose
#----------------------- rdirc.csh --------------------------------------
# Revised - July 2013 - nanr
# -------------------------------------------------------------------
#
# this cshell script runs the program topo2rdirc.F90 
#
# -------------------------------------------------------------------
# Notes
# -------------------------------------------------------------------
# Start with topography and end with river directions map.
# Simpler procedure than Graham et al. (1999), because
# we have no idea where the rivers were in the geologic past.
#
# The output consists of two ascii files: fort.10 and fort.11.
# The former is the river directions map in the format required by CLM.
# The latter lists all grid cells involved in infinite river loops.
# The user must redirect rivers in the vicinity of infinite loops.
# Alternatively, the user may start with topography which always slopes to
# the ocean, doesn't have internal basins, and doesn't contain completely
# flat plateaus.
#
# Currently written for Cretaceous topography.
#
# slevis, feb 2003
#
#-----------------------------------------------------------------
# YOU NEED  : Makefile,topo2rdirc_sed.F90,rdirc.csh
# TO RUN    : rdirc.csh
# NOTE      : make sure "limits" on your machine are set appropriately
#           : if not set large enough, the program may not have enough
#             space/memory to run
#-----------------------------------------------------------------
#
#----------------------------------------------------------------------------
# specify input/output pathnames with enviromental variables set in this script 
#----------------------------------------------------------------------------

set INFILE  = /myPath/myTopo.nc
set CASE    = myCase

#----------------------------------------------------------------------------
# use "sed" to insert pathnames into fortran code 
#----------------------------------------------------------------------------

cat >! pathname.sed << EOF
s#input_topo_data#$INFILE#
EOF
cp  topo2rdirc_sed.F90 topo2rdirc.F90 
sed -f pathname.sed topo2rdirc_sed.F90 >! topo2rdirc.F90

#----------------------------------------------------------------------------
#  compile and run 
#----------------------------------------------------------------------------
cp  $INFILE . 
make EXENAME=topo2rdirc
echo '' > fort.10
echo '' > fort.11
echo '' > fort.12
echo '' > fort.13
echo '' > fort.14
echo '' > fort.15
topo2rdirc

#---------------------------------------------------------------
# output fort.10 and fort.11 will in working directory 
# cp fort.10 to original file to save 
mv fort.10 fort.10_$CASE # original direction file
mv fort.11 fort.11_$CASE # original loops file
mv fort.12 fort.12_$CASE # original internal basin file
mv fort.14 fort.14_$CASE # modified direction file
mv fort.13 fort.13_$CASE # modified direction file
mv fort.15 fort.15_$CASE # modified internal basin file

#----------------------------------------------------------------
#----------------------------------------------------------------------------


exit



