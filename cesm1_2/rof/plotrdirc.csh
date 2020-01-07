#! /bin/csh -f

set echo
set verbose

#---------------------------------------------------------------------------
# Make a plot of runoff vectors from rdirc file: 
#---------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Begin User Modify : specify input/output pathnames 
#----------------------------------------------------------------------------

set IFILE  = <topography-bathymetry_file>
setenv NLAT      180
setenv NLON      360
set RESOLN   =   1x1

# set IFILE  = topo.1x1deg.nc
# setenv NLAT      180
# setenv NLON      360
# set RESOLN   =   1x1

# set IFILE  = topo.2x2deg.nc
# setenv NLAT      90
# setenv NLON      180
# set RESOLN  =    2x2

setenv CASE   <casename>
setenv PLOTNAME rdirc_${IFILE}

#----------------------------------------------------------------------------
# End User Modify
#----------------------------------------------------------------------------

# RFILE is the fort.10 output by rdirc.csh + topo2rdirc.F90
setenv RFILE1    fort.10_$CASE
setenv EFILE1    fort.11_$CASE

setenv RFILE2    fort.13_$CASE
setenv EFILE2    fort.14_$CASE

#---------------------------------------------------------------------------
#  End of user input
#---------------------------------------------------------------------------

if (-e rdirc.$RESOLN.$CASE) then
	echo 'removing old file: ' rdirc.$RESOLN.$CASE
	rm rdirc.$RESOLN.$CASE
endif

cp fort.13_$CASE rdirc.$RESOLN.$CASE

ncl < plot_rdirc.ncl


#---------------------------------------------------------------------------
exit(0)

