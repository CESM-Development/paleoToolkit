;=========================================================; 
; reads in data merged topo kmt file
; reads in landfrac from cesm1 surface_data file 
; (run mksurfdata offline first)
; read in template cam bnd_topo (Netcdf) file 
; makes modifications, writes out new file
; 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;====================user defined parameters

dd= systemfunc("date -u +%y%m%d")

scale_Factor = 0.16
write_landm_coslat = True		; LANDM_COSLAT is not used in CAM5; set to True for CAM4; 


 case = "permian"
 print("Running Case = " + case)

 ifile1 = "bnd_topo_permian_1.9x2.5_remap.150331.nc"
 ofile  = "bnd_topo_permian_1.9x2.5_remap_sgh30.150331.nc"

system("rm "+ofile)
fout = addfile(ofile,"c")

; copy input variabls to new output file

in1 = addfile(ifile1,"r")
in1_global= getvaratts(in1)
in1_vars = getfilevarnames(in1)
; copy global attributes
if (.not.all(ismissing(in1_global))) then
do i = 0, dimsizes(in1_global) - 1
	print("copy_fileatts: global attributes->" + in1_global(i) )
	fout@$in1_global(i) $ = in1@$in1_global(i)$
end do
end if
; copy variables
if (.not.all(ismissing(in1_vars))) then
do i = 0, dimsizes(in1_vars) - 1
	if (in1_vars(i) .eq."LANDM_COSLAT" .and. write_landm_coslat .eq. False) then
		print("EXCLUDING var ->" + in1_vars(i) )
	else
		fout->$in1_vars(i)$ = in1->$in1_vars(i)$
		print("writing var ->" + in1_vars(i) )
	end if
end do
end if

; create and write new sgh30 variable
; For paleo, we approximate sgh30 as hi_res_topo*0.1
; But we don't have the hi_res_topo, and it caused problems
; to try to use the 10min topo and interpolate, now we
; approximate using SGH30=SGH*0.16, based on a global mean PD sgh30/sgh ratio.
sgh = in1->SGH
sgh30 = sgh
sgh30 = sgh*scale_Factor

fout->SGH30 = sgh30

; create file attributes
date_written = systemfunc("date -u +%y%m%d")
user  = systemfunc("whoami")
scriptFile = get_script_prefix_name()
fout@history = "modified by " + user + " with "+scriptFile+ " on " + date_written 
fout@sgh30_scaleFactor = "sgh30 = sgh*"+scale_Factor+" based on PD sgh30/sgh ratio."


end
    


