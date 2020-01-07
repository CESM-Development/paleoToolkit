#!/bin/csh

# Notes:
# - will build the CCSM runoff correcting/smoothing code in ./obj
# - must specify location of src code, Makefile, Macros file, dependancy generator 
#===============================================================================

setenv SRCDIR `pwd`/src
setenv TOOLDIR `pwd`/tools

echo source dir: $SRCDIR

if !(-d obj) mkdir obj
cd obj

cc -o makdep $TOOLDIR/makdep.c

echo $SRCDIR  >! Filepath

gmake VPFILE=Filepath THREAD=TRUE -f $SRCDIR/Makefile MACFILE=$SRCDIR/Macros.cheyenne || exit -1

cd ..
rm              runoff_map
ln -s obj/a.out runoff_map

