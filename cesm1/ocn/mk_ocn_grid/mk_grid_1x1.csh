#!/bin/csh -fv 

set echo
set verbose
###################################################################
# Required input file:  
# topography/bathymetry file at .5x.5 degree resolution.
# This can be obtained from the 2x2 degree resolution data
# using create.05deg.ncl in this ocn/src/gridgen.
# usage:  % ncl create.05degree.ncl
# be sure to compile code in src dir and place executable in bin dir 
# makefile for ocean code found in ocn/bin

###################################################################
#Define global variables:

set MYDIR = USERPATH/paleo_setup/
set CASE  = permian
set ITER  = 4   			# iteration 
set TOPO  = topo.1deg.nc
set DATE  = "`date +%y%m%d`"

setenv DATADIR  USERPATH/paleo_setup/ocn/griddata
setenv SRCDIR   $MYDIR/ocn/mk_ocn_grid/
setenv WKDIR    $MYDIR/ocn/mk_ocn_grid/$CASE
setenv MYDATA   $MYDIR/topobath/   # cp topo data here

if !(-d ${WKDIR}) mkdir -p ${WKDIR}

if ( (! -e $MYDATA/$TOPO) ) then
	echo ''
	echo 'FATAL ERROR: topo file does not exist locally or in ' $MYDATA '. Create new file using create.05degree.ncl'
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
set latsp   =  -69.   		# latitude  of South Pole (kmt.3.da)

# change iteratively after you look at your grid.
# Use these to adjust gridcell distribution so you don't have big
# cells in the N. Hem and small cells in the S. Hem, for example.
set nlatn   = 214   		# number of j grid lines in NH	
set nlats   = 170   		# number of j grid lines in SH
# set nlatn   = 63   		# number of j grid lines in NH
# set nlats   = 53   		# number of j grid lines in SH

set nx      = 320  		# number of i grid lines
set nz      = 60		# number of vertical grid levels (gx1v6 = 60)
set dyeq    = .26   		# dy in degrees at Equator 	        (roughly equal to modern TLAT)
set dsig    = 23.  		# Gaussian efold-scale at equator	(Dan Lunt cret)
set jcon    = 11    		# jcon = # rows of constant dy at poles


# Don't change
set popgrid = grid.$ITER.pop.da	# binary grid file
set pltgrid = grid.$ITER.plot.da	# binary grid file
set cdfgrid = gridkmt.$ITER.nc		# netcdf grid file
set depgrid = h.$ITER.da		# binary depth array
set kmtgrid = kmt.$ITER.da		# binary kmt file
set vrtgrid = gx1v6_vert_grid	# ascii vertical grid input file
set topo     = $TOPO  	        # netcdf topography input dataset
set minz     = 5.		# min z in meters (Don't lower!!!)
set mink     = 3		# minimum allowed km value (Don't lower!!!)
@ ny      = $nlatn + $nlats  	# number of j grid lines

# not used
set iocean   = 50		# ocean (i,j) point used to define
set jocean   = 60		#	contiguous ocean.

# !!!SET PLOTTING VIEWPOINTS BELOW!!!

#################################
# change to work directory
#################################
chdir $WKDIR
#################################
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
$topo
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
chdir $WKDIR
ln -s $MYDATA/$TOPO .
cp $DATADIR/$vrtgrid . 
cp $DATADIR/$vrtgrid . 
cp $SRCDIR/ns_dipole .
cp $SRCDIR/paleotopo .
cp $SRCDIR/grid_bin2nc .

#################################
# RUN GRID ROUTINES
# This must have been compiled. 
# Executable is copied into $WKDIR.
#################################
ns_dipole < input.$ITER.ns_dipole  >! output.$ITER.ns_dipole
paleotopo < input.$ITER.paleotopo  >! output.$ITER.paleotopo
grid_bin2nc < input.$ITER.grid_bin2nc >! output.$ITER.bin2nc

#################################
# CLEAN UP
#################################
#/usr/bin/rm -f input.ns_dipole input.paleotopo input.display_hgrid_land

exit(0)
