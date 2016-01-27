#! /bin/csh -f

set echo 
set verbose

#---------------------------------------------------------------------------
# Make the input_template directory and input_template files for ocn and ice
#---------------------------------------------------------------------------

# ----------------------------------------
# - User defined -------------------------
# ----------------------------------------
setenv CASE MAA_mod
setenv GRID grid.1.pop.da
#setenv KMT kmt.m.140502.ieeei4
setenv KMT kmt.ieeei4

setenv CSMROOT  /glade/p/work/crtabor/cesm1_2_0/models/ocn/pop2/input_templates
setenv CSMCODE  $CSMROOT/models

#Set New Ocn/Ice Grid parameters:
#set    NX = 100
#set    NY = 116
set    NX = 320
set    NY = 384

# Set OCNGRID to new pop grid name 
setenv OCNGRID  MAA_mod_gx1v6

# ----------------------------------------
# - End User defined -------------------------
# ----------------------------------------
#Set up for making regionmask file:

setenv GRIDFILE   /glade/p/work/crtabor/paleo_setup/ocn/mk_ocn_grid/MAA_mod/$CASE/$GRID
setenv KMTFILE    /glade/p/work/crtabor/paleo_setup/ocn/mk_ocn_grid/MAA_mod/$CASE/$KMT

setenv CSMOCN  /glade/p/work/crtabor/cesm1_2_0/models/ocn/pop2/input_templates
setenv CSMICE  /glade/p/work/crtabor/cesm1_2_0/models/ice/cice/src/source/ice_domain_size.F90

# OLDGRID is the pop grid name for input_template files to be copied
# setenv OLDGRID  gx3v5 
setenv OLDGRID  gx1v6 

  if ($OLDGRID == gx3v5) then
    set NXOLD = 100 
    set NYOLD = 116 
    set  NCAT = 5
  endif
  if ($OLDGRID == gx1v6) then
    set NXOLD = 320 
    set NYOLD = 384 
    set  NCAT = 10		# ??
  endif

# Number of Ice categories:
# Don't change! (unless you really want to)

setenv ICEGRID  ${NX}x${NY}x${NCAT}

##########################################################################
#1) Make regionmask file region.ieeei4:

#get files for input:
 cp $GRIDFILE grid.pop.da
 cp $KMTFILE  kmt.da

#get region code, edit, compile and run:

cat >! commands.sed << EOF
s#NX#${NX}#
s#NY#${NY}#

EOF
 sed -f commands.sed modregmsk_edit.f >! modregmsk.f 

 ifort -convert big_endian -o modregmsk.be modregmsk.f
#ifort -convert big_endian -o modregmsk.le modregmsk.f

 ./modregmsk.be
 	mv region.ieeei4 region.$CASE.be.ieeei4
 # ./modregmsk.le					!! to prove that kmt.da and grid.da are BigEndian
 # 	mv region.ieeei4 region.$CASE.le.ieeei4		!! to prove that kmt.da and grid.da are BigEndian


##########################################################################
#copy over unchanged input_template files ice and ocn: 

foreach FILE (depth_accel history_contents \
              movie_contents tavg_contents \
              vert_grid )

#cp ${OLDGRID}_$FILE ${OCNGRID}_$FILE.$CASE
cp $CSMOCN/${OLDGRID}_$FILE ${OCNGRID}_$FILE
#cp $CSMOCN/input_templates/${OLDGRID}_$FILE $OCNINPUT/${OCNGRID}_$FILE

end

#copy over the ones to change locally:
cp ice_model_size.F.${NXOLD}x${NYOLD}x5 ice_model_size.F 
# cp $CSMICE/input_templates/ice_model_size.F.${NXOLD}x${NYOLD}x5 ice_model_size.F 

# cp ${OLDGRID}_model_size.F ocn_model_size.F 
cp $CSMOCN/${OLDGRID}_model_size.F ocn_model_size.F 


##########################################################################

#Create new input files:
#Create region_ids
# Note: for negative valued region_ids (i.e. marginal seas), 
# you need to set a lat, lon location and area for the
# redistribution of net freshwater from the marginal sea
# (when marginal sea balance is turned on.)
# The example below has these set to zero.
# See user's guide for more info on how to set this.
   
cat >! region_ids << EOF
   1  'Southern Hemisphere '        0.0   0.0      0.0
   2  'Northern Hemisphere '        0.0   0.0      0.0
EOF
! $OCNGRID 

#Create transport_contents:
# This file is for diagnostic calculation of section transports.
# First line tells how many sections to integrate.
# Next lines are i,j,k limits of sections, whether
# section is meridional or zonal {(merid,zonal)} and a section name.
# Transport output appears in dt files.
# 

#imin imax jmin jmax kmin kmax sect   Name - NOTE:  Spacing is EXACT!
cat >! transport_contents << EOF
2
  30   36  215  215    1   60 zonal  Equatorial bay
 195  195  100  125    1   60 merid  Mid-ocean ridge
EOF

#2345678901234567890123456789012345678901234567890123456789012345678901234567890
# cat >! transport_contents << EOF
# 3
# 33   33    3   11    1   25 merid  Drake Passage
# 95   95   20   41    1   25 merid  Indonesian Throughflow
# 75   92   57   57    1   25 zonal  Tethys Seaway

# 10   10   55   65    1   25 merid  TESTM
# 45   55   10   10    1   25 zonal  TESTZ

# 64   64    1  128    1   20 zonal  sample zonal section
#  1  192   64   64    1   10 merid  sample meridional section

#EOF

#END OF USER EDITING...


#using sed, new NX,NY, create model_size.F for ice and ocn
cat >! commands.sed << EOF
s#${NXOLD}#${NX}#
s#${NYOLD}#${NY}#

EOF

sed -f commands.sed ocn_model_size.F >! model_size.F 
sed -f commands.sed ice_model_size.F >! ice_size.F 

#s/cp new files over to new input_templates
foreach FILE (region_ids model_size.F transport_contents)
 cp $FILE ${OCNGRID}_$FILE.$CASE
 # cp $FILE $OCNINPUT/${OCNGRID}_$FILE
end
 cp $CSMICE ice_model_size.F.${NX}x${NY}x${NCAT}

exit(0)
