#!/bin/csh -fv 

set echo off
set verbose

set DATE    = "`date +%y%m%d`"

set popgriddir = <workspace>

set scripdir   = <workspace>


set nx      = 320  			# number of i ocn gridpoints (longitudes)
set ny      = 384   	 		# number of j ocn gridpoints (latitudes)
# -----------------------------------------------------------------
# Input files -----------------------------------------------------
# -----------------------------------------------------------------
set popgrid     = grid.<iter>.pop.da			# binary grid file from mk_grid.csh
set kmtgrid     = kmt.<iter>.da			# binary kmt  file from mk_grid.csh

set popgridfile = $popgriddir/$popgrid		
set kmtgridfile = $popgriddir/$kmtgrid

# -----------------------------------------------------------------
# Output labels -----------------------------------------------------
# -----------------------------------------------------------------
set ocngridname = <gridname>

set ocnres = gx1<case>

# -----------------------------------------------------------------
# Output filenames ------------------------------------------------
# -----------------------------------------------------------------
set ocngridfile = ${ocnres}_$DATE.nc

#################################
# SET HERE FILE INPUT
#################################

cat >! input.scrip_grid << EOF
$nx $ny
${ocngridname}
${popgridfile}
${kmtgridfile}
${ocngridfile}
1
EOF


#################################
# RUN SCRIP ROUTINES
#################################
 ${scripdir}/myconvertPOPT < input.scrip_grid


#################################
# ADD DOCUMENTATION
#################################
foreach i ($ocngridfile )
  ncatted -O -h -a Created_by,global,a,c,"`whoami` on `date`; using pop input gridfile: $popgrid" $i
  ncatted -O -h -a 1D_grid_indexing,global,a,c,"if n is 1D index, i runs fast, j runs slow: n=(j-1)*fast_grid_dim+i" $i
end

#################################
# CLEAN UP
#################################


exit(0)
