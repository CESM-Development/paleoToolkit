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

# cheyenne specific batch commands:
#PBS -A #########
#PBS -N regridbatch_paleo
#PBS -q regular
#PBS -l select=8:ncpus=36:mpiprocs=36:ompthreads=1
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -S /bin/csh -V


#----------------------------------------------------------------------
# Set parameters
#----------------------------------------------------------------------
# Which version of CLM to generate mapping files for
# Can be clm4_0 or clm4_5
phys="clm4_0"
#phys="clm4_5"

# hardcode paleo model resolution for output grid
# e.g. "1.9x2.5", "0.9x1.25"
resols="<output_resolution>"
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

   csmsrc="<directory_with_model_code>"   
   time env CSMSRC=$csmsrc REGRID_PROC=$regrid_num_proc ./mkmapdata_paleo.sh $cmdargs
done
