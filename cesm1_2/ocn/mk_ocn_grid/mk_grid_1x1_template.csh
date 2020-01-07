#!/bin/csh -fv 

set echo
set verbose
###################################################################
# Required input file:  
# topography/bathymetry file at 0.5,1 or 2deg (others not tested)
# variable in topo/bath file is assumed to be named "topo"
###################################################################
#Define global variables:

set CASE  = <casename>
set ITER   = 1   			# iteration 
set topopath = <path_to_topo_file>
set topofile  = <topo-bath_file>
set vrtgridpath = <path_to_vrtgrid_file>
set vrtgrid = <vertical_grid_file>
set DATE  = "`date +%y%m%d`"

if ( (! -e $topopath/$topofile) ) then
	echo ''
	echo 'FATAL ERROR: topo file does not exist',$topopath/$topofile
	exit
   endif
endif
if ( (! -e $vrtgridpath/$vrtgrid) ) then
	echo ''
	echo 'FATAL ERROR: vertical grid file does not exist',$vrtgridpath/$vrtgrid
	exit
   endif
endif

###################################################################
# User changes:  For each new horizontal grid you must edit the following parameters:
###################################################################

# Change for your continent distribution 
# ncview your topo_0.5x0.5 map to see where to place the poles over land.
# NOTES:  latitude and  longitude must be btwn -180 : 180
set lonnp   =   50.  		# longitude of North Pole (must be btwn -180 : 180)
set latnp   =   75.   		# latitude  of North Pole
set lonsp   =   50.   		# longitude of South Pole (must be btwn -180 : 180)
set latsp   =  -69.   		# latitude  of South Pole

# change iteratively after you look at your grid.
# Use these to adjust gridcell distribution so you don't have big
# cells in the N. Hem and small cells in the S. Hem, for example.
set nlatn   = 205   		# number of j grid lines in NH
set nlats   = 179   		# number of j grid lines in SH

set nx      = 320  		# number of i grid lines
set nz      = 60		# number of vertical grid levels (gx1v6 = 60)
set dyeq    = .26   		# dy in degrees at Equator 	        (roughly equal to modern TLAT)
set dsig    = 23.  		# Gaussian efold-scale at equator
set jcon    = 11    		# jcon = # rows of constant dy at poles


# Don't change
set popgrid = grid.$ITER.pop.da	# binary grid file
set pltgrid = grid.$ITER.plot.da	# binary grid file
set cdfgrid = gridkmt.$ITER.nc		# netcdf grid file
set depgrid = h.$ITER.da		# binary depth array
set kmtgrid = kmt.$ITER.da		# binary kmt file
set minz     = 5.		# min z in meters (Don't lower!!!)
set mink     = 3		# minimum allowed km value (Don't lower!!!)
@ ny      = $nlatn + $nlats  	# number of j grid lines

# not used
set iocean   = 50		# ocean (i,j) point used to define
set jocean   = 60		#	contiguous ocean.

# !!!SET PLOTTING VIEWPOINTS BELOW!!!

#################################
# SET HERE FILE INPUT
#################################

cat >! input.$ITER.ns_dipole << EOF
$nx $nlatn $nlats
$lonnp $latnp
$lonsp $latsp
$dyeq $dsig $jcon
1
$pltgrid
1
$popgrid
EOF

cat >! input.$ITER.paleotopo << EOF
$nx $ny $nz
$pltgrid
$vrtgrid
$topofile
1
$depgrid
$minz
1
0
$mink
1
$kmtgrid
EOF

cat >! input.$ITER.grid_bin2nc << EOF
$nx
$ny
$cdfgrid
$vrtgrid
$kmtgrid
$pltgrid
$popgrid
EOF

#################################
#get input data:
ln -s $topopath/$topofile .
ln -s $vrtgridpath/$vrtgrid .

#################################
# RUN GRID ROUTINES
# This must have been compiled. 
#################################
./ns_dipole < input.$ITER.ns_dipole  >! output.$ITER.ns_dipole
./paleotopo < input.$ITER.paleotopo  >! output.$ITER.paleotopo
./grid_bin2nc < input.$ITER.grid_bin2nc >! output.$ITER.bin2nc

#################################
# CLEAN UP
#################################
#/usr/bin/rm -f input.ns_dipole input.paleotopo input.display_hgrid_land

exit(0)
