#!/bin/csh -fv 

set echo off
set verbose

set DATE    = "`date +%y%m%d`"

set CASE    = LGM

set popgriddir = /glade/p/cesm/cseg/inputdata/ocn/pop/gx1v6/grid/

set scripdir   = /glade/u/home/nanr/setup_tools_cesm1/cpl_mapping/scrip1.4/


set nx      = 320  			# number of i ocn gridpoints (longitudes)
set ny      = 384   	 		# number of j ocn gridpoints (latitudes)
# -----------------------------------------------------------------
# Input files -----------------------------------------------------
# -----------------------------------------------------------------
set popgrid     = horiz_grid_20010402.ieeer8	        # binary grid file from mk_grid.csh
set kmtgrid     = kmt_gx1v6_lgm21ka.110315.ieeei4	# binary kmt  file from mk_grid.csh

set popgridfile = $popgriddir/$popgrid		
set kmtgridfile = $popgriddir/$kmtgrid

set atmdomain   = fv0.9x1.25_070727.nc		# cesm1 atm domain
set atmgridfile = /glade/p/cesm/cseg/mapping/grids/$atmdomain
# -----------------------------------------------------------------
# Output labels -----------------------------------------------------
# -----------------------------------------------------------------
set map1_name   = 'gx1v6LGM to fv_0.9x1.25 Mapping('${kmtgrid}' to '${atmdomain}')'
set map2_name   = 'fv_0.9x1.25 to gx1v6LGM Mapping('${atmdomain}' to '${kmtgrid}')'
set ocngridname = 'LGM'

set ocnres = 'gx1LGM'
set atmres = 'fv09_1.25'

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
