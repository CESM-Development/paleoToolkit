#!/bin/bash
#
#
# Batch script to submit to create mapping files for all standard
# resolutions.  If you provide a single resolution via "$RES", only
# that resolution will be used. In that case: If it is a regional or
# single point resolution, you should set 'BSUB -n' to 1, and be sure
# that '-t regional' is specified in cmdargs.
#
# Currently only setup to run on yellowstone/caldera/geyser. Note that
# geyser is needed for very high resolution files (e.g., 1 km) because
# of its large memory per node, so that is set as the default.
# However, for coarser resolutions, you may get better performance on
# caldera or yellowstone.
# 
# yellowstone specific batch commands:
#BSUB -P P93300624
#BSUB -n 8
#BSUB -R "span[ptile=8]" 
#BSUB -o regrid.%J.out   # ouput filename
#BSUB -e regrid.%J.err   # error filename
#BSUB -J regrid          # job name
#BSUB -W 24:00
#BSUB -q geyser          # queue

#----------------------------------------------------------------------
# Set parameters
#----------------------------------------------------------------------
# point to user modified input script
# e.g., cp mkmapdata_paleo.sh mkmapdata-mycase.sh
# modify mkmapdata-mycase.sh to point to your SCRIPgrid_myRun_RESOLTION (look for "Change me")
mycase="mkmapdata_paleo.sh"

# Which version of CLM to generate mapping files for
# Can be clm4_0 or clm4_5
phys="clm4_0"
#phys="clm4_5"

# hardcode paleo model resolution for output grid
# e.g. "1.9x2.5", "0.9x1.25"
resols="1.9x2.5"
echo "Run for paleo  $resols resolution"

#----------------------------------------------------------------------
# Begin main script
#----------------------------------------------------------------------

# if [ -z "$RES" ]; then
   # echo "Run for all valid resolutions"
   # resols=`../../../bld/queryDefaultNamelist.pl -res list -silent -phys $phys`
# else
   # resols="$RES"
# fi
# echo "Create mapping files for this list of resolutions: $resols"


#----------------------------------------------------------------------

for res in $resols; do
   echo "Create mapping files for: $res"
#----------------------------------------------------------------------
   cmdargs="--phys $phys -r $res -b"

   # For single-point and regional resolutions, tell mkmapdata that
   # output type is regional
   if [[ `echo "$res" | grep -c "1x1_"` -gt 0 || `echo "$res" | grep -c "5x5_"` -gt 0 ]]; then
       res_type="regional"
   else
       res_type="global"
   fi

   cmdargs="$cmdargs -t $res_type"

   if [ $res_type = "regional" ]; then
       # For regional and (especially) single-point grids, we can get
       # errors when trying to use multiple processors - so just use 1.
       # We also do NOT set batch mode in this case, because some
       # machines (e.g., yellowstone) do not listen to REGRID_PROC, so to
       # get a single processor, we need to run mkmapdata.sh in
       # interactive mode.
       regrid_num_proc=1
   else
       regrid_num_proc=8
       if [ ! -z $LSF_PJL_TYPE ]; then
	   cmdargs="$cmdargs -b"
       fi
   fi

   time env REGRID_PROC=$regrid_num_proc $mycase $cmdargs
   # time env REGRID_PROC=$regrid_num_proc ./mkmapdata_paleo.sh $cmdargs
done
