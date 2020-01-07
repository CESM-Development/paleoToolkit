#!/bin/csh -f
set echo
set verbose
#----------------------- paleo_mkraw_cesm1.csh --------------------------------------
#
# all work is done in directory where script and code resides 
#
# this cshell script runs the program paleo_mkraw_cesm1.F90 
#
# see code for details (comments)
# input:  
# input data: (2x2) LSM vegetation types from desired period
#           : (mksrf_soitex.10level.nc) IGBP soil texture for 0k
#           : zonally averaged model organic soil density distribution **
#           : 2x2 topography/bathymetry file (topo_depth) file 
#
# output data: raw data files necessary to create surface-data file
#              at model runtime: mksrf_glacier_paleo.nc
#                              : mksrf_urban_paleo.nc
#                              : mksrf_lanwat_paleo.nc
#                              : mksrf_soicol_paleo.nc
#                              : mksrf_lai_paleo.nc
#                              : mksrf_soitex_paleo.nc
#CAS
# updates for cesm1/ccsm4 includes 
#			       : mksrf_organic_paleo.nc
#			       : mksrf_fmax_paleo.nc
#                              : mksrf_landuse_paleo.nc (replaces mksrf_pft_paleo.nc)
#                              : mksrf_topo_paleo.nc
#			       : mksrf_vocef_paleo.nc
#
#
#  mkrsf_* files created here can be used for offline cesm1/ccsm4 tool, mksurfdata,
#  to create the surface_data set for cesm1/ccsm4
#
#  original:  cesm1/ccsm4 code base:  models/lnd/clm/tools/mksurfdata
#
#-----------------------------------------------------------------
# YOU NEED  : Makefile,paleo_mkraw_cesm1_sed.F90,paleo_mkraw_cesm1.csh,input data
# TO RUN    : paleo_mkraw_cesm1.csh
# NOTE      : make sure "limits" on your machine are set appropriately
#           : if not set large enough, the program may not have enough
#             space/memory to run
#-----------------------------------------------------------------
#
#----------------------------------------------------------------------------
# specify input/output pathnames with enviromental variables set in this script 
#----------------------------------------------------------------------------

setenv CASE  <casename>
set DATE  = "`date +%y%m%d`"
set INFILE_FORMAT = netcdf              # [netcdf, ascii]

# ----------------
#  -- INPUT FILES
# ----------------
# Paleo veg and topo (2 degree)
set INPUT_LSM_DATA = <lsm_file>
set INPUT_TOP_DATA = <topo-bath_file>

# Soil texture (use present day)
set INPUT_SOI_DATA = /glade/p/cesmdata/cseg/inputdata/lnd/clm2/rawdata/mksrf_soitex.10level.c010119.nc 
set INPUT_ORG_DATA = mksrf_zon_organic.10level.nc 

# ----------------
# -- OUTPUT files
# ----------------
set OUTPUT_GLACIER = mksrf_glacier_$CASE.c$DATE.nc
set OUTPUT_URBAN   = mksrf_urban_$CASE.c$DATE.nc
set OUTPUT_LANWAT  = mksrf_lanwat_$CASE.c$DATE.nc
set OUTPUT_SOICOL  = mksrf_soicol_$CASE.c$DATE.nc
set OUTPUT_PFT     = mksrf_landuse_$CASE.c$DATE.nc
set OUTPUT_LAI     = mksrf_lai_$CASE.c$DATE.nc
set OUTPUT_SOITEX  = mksrf_soitex_$CASE.c$DATE.nc
set OUTPUT_ORGANIC = mksrf_organic_$CASE.c$DATE.nc
set OUTPUT_FMAX    = mksrf_fmax_$CASE.c$DATE.nc
set OUTPUT_TOPO    = mksrf_topo_$CASE.c$DATE.nc
set OUTPUT_VOCEF   = mksrf_vocef_$CASE.c$DATE.nc

#----------------------------------------------------------------------------
# use "sed" to insert pathnames into fortran code 
#----------------------------------------------------------------------------

cat >! pathnames.sed << EOF
s#input_lsm_data#$INPUT_LSM_DATA#
s#input_top_data#$INPUT_TOP_DATA#
s#input_soi_data#$INPUT_SOI_DATA#
s#input_org_data#$INPUT_ORG_DATA#
s#output_glacier#$OUTPUT_GLACIER#
s#output_urban#$OUTPUT_URBAN#
s#output_lanwat#$OUTPUT_LANWAT#
s#output_soicol#$OUTPUT_SOICOL#
s#output_pft#$OUTPUT_PFT#
s#output_lai#$OUTPUT_LAI#
s#output_soitex#$OUTPUT_SOITEX#
s#output_organic#$OUTPUT_ORGANIC#
s#output_fmax#$OUTPUT_FMAX#
s#output_topo#$OUTPUT_TOPO#
s#output_vocef#$OUTPUT_VOCEF#
EOF
sed -f pathnames.sed paleo_mkraw_cesm1_sed.F90 >! paleo_mkraw_cesm1.F90

#----------------------------------------------------------------------------
#  compile and run 
#----------------------------------------------------------------------------
make EXENAME=paleo_mkraw_cesm1
./paleo_mkraw_cesm1


exit



