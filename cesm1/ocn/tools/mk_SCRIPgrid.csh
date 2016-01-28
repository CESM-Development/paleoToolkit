#!/bin/csh -fv 

set echo off
set verbose

set DATE    = "`date +%y%m%d`"

set CASE    = USER_DEFINED

set popgriddir = USERPATH/paleo_setup/ocn/mk_ocn_grid/$CASE

set scripdir   = USERPATH/paleo_setup/cpl_mapping/scrip1.4/


set nx      = 320  			# number of i ocn gridpoints (longitudes)
set ny      = 384   	 		# number of j ocn gridpoints (latitudes)
# -----------------------------------------------------------------
# Input files -----------------------------------------------------
# -----------------------------------------------------------------
set popgrid     = grid.pop.da			# binary grid file from mk_grid.csh
set kmtgrid     = kmt.da			# binary kmt  file from mk_grid.csh

set popgridfile = $popgriddir/$popgrid		
set kmtgridfile = $popgriddir/$kmtgrid

set atmdomain   = fv1.9x2.5_090205.nc		# cesm1 FV2 atm domain
set atmgridfile = /glade/p/cesm/cseg/mapping/grids/$atmdomain
# -----------------------------------------------------------------
# Output labels -----------------------------------------------------
# -----------------------------------------------------------------
set map1_name   = 'gx1v6-perm to fv_0.9x1.25 Mapping('${kmtgrid}' to '${atmdomain}')'
set map2_name   = 'fv_0.9x1.25 to gx1v6-perm Mapping('${atmdomain}' to '${kmtgrid}')'
set ocngridname = 'permian'

set ocnres = 'gx1PT'
set atmres = 'fv19_25'

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
 ${scripdir}/grids/myconvertPOPT < input.scrip_grid


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
