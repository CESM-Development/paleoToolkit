#! /bin/csh -f

set echo 
set verbose

#---------------------------------------------------------------------------
# Make the input_template directory and input_template files for ocn and ice
#---------------------------------------------------------------------------

# ----------------------------------------
# - User defined -------------------------
# ----------------------------------------
setenv CASE <casename>
setenv GRID grid.<iter>.pop.da
setenv KMT kmt.<iter>.da

set    NX = 320
set    NY = 384

# ----------------------------------------
# - End User defined -------------------------
# ----------------------------------------
#Set up for making regionmask file:
setenv GRIDFILE   ./$GRID
setenv KMTFILE    ./$KMT

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

#Create new input files:
#Create region_ids
# Note: for negative valued region_ids (i.e. marginal seas), 
# you need to set a lat, lon location and area for the
# redistribution of net freshwater from the marginal sea
# (when marginal sea balance is turned on.)
# The example below has these set to zero.
# See user's guide for more info on how to set this.
   
cat >! region_ids_$CASE << EOF
   1  'Southern Hemisphere '        0.0   0.0      0.0
   2  'Northern Hemisphere '        0.0   0.0      0.0
EOF

#cat >! region_ids_$CASE <<EOF
# EOF
#   1  'Arctic Ocean        '        0.0   0.0      0.0
#   2  'North Tethys Ocean  '        0.0   0.0      0.0
#   3  'North Atlantic Ocean'        0.0   0.0      0.0
#   4  'Southern Ocean      '        0.0   0.0      0.0
#   5  'South Atlantic Ocean'        0.0   0.0      0.0
#   6  'South Tethys Ocean  '        0.0   0.0      0.0
#   7  'Indian Ocean        '        0.0   0.0      0.0
#   8  'Pacific Ocean       '        0.0   0.0      0.0
#EOF

#Create transport_contents:
# This file is for diagnostic calculation of section transports.
# First line tells how many sections to integrate.
# Next lines are i,j,k limits of sections, whether
# section is meridional or zonal {(merid,zonal)} and a section name.
# Transport output appears in dt files.
# 

#imin imax jmin jmax kmin kmax sect   Name - NOTE:  Spacing is EXACT!
cat >! transport_contents_$CASE << EOF
2
  30   36  215  215    1   60 zonal  Equatorial bay
 195  195  100  125    1   60 merid  Mid-ocean ridge
EOF

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


exit(0)
