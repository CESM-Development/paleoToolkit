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

sgh30_2_sgh_SF = 0.16
write_landm_coslat = True	


cases  = (/"<casename>"/)
ncases = dimsizes(cases)

do c = 0,ncases-1

 case = cases(c)
 print("Running Case = " + case)

ifile1 = "<bnd_topo_file_in>"
ofile  = "<bnd_topo_file_out>"

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

;======= create landm_coslat field
;cas
;creating landm_coslat for appropriate topography

if(write_landm_coslat.eq.True) then

 landfrac = in1->LANDFRAC 
 lat = in1->lat
 lon = in1->lon
 nlat = dimsizes(lat)
 nlon = dimsizes(lon)
 landm_sm = smth9(landfrac,.5,.25,True)
 landm_coslat = new((/nlat,nlon/),double)
 degrad = 0.0174533  ; conversion of degrees to radians
 latr = lat*degrad
 latrC = conform(landfrac,latr,0)
 landm_coslat = landfrac
 do nlt = 0, nlat-1
 do nln = 0, nlon-1
  if(landfrac(nlt,nln).lt.0.999)then
   landm_coslat(nlt,nln) = landm_sm(nlt,nln)*cos(latrC(nlt,nln))
  end if
 end do
 end do
 landm_coslat!0 = "lat"
 landm_coslat!1 = "lon"
 landm_coslat&lat = lat
 landm_coslat&lon = lon
 fout->LANDM_COSLAT = landm_coslat
 print("LANDM_COSLAT fields created from interpolated LANDFRAC")

end if

;---

; create and write new sgh30 variable
; For paleo, we approximate sgh30 as hi_res_topo*0.1
; But we don't have the hi_res_topo, and it caused problems
; to try to use the 10min topo and interpolate, so I will
; approximate from PHIS, and from SGH, which are already smoothed.
sgh = in1->SGH
sgh30 = sgh
sgh30 = sgh*sgh30_2_sgh_SF

fout->SGH30 = sgh30

; create file attributes
date_written = systemfunc("date -u +%y%m%d")
user  = systemfunc("whoami")
scriptFile = get_script_prefix_name()
fout@history = "modified by " + user + " with "+scriptFile+ " on " + date_written 
fout@sgh30_scaleFactor = "sgh30 = sgh*"+sgh30_2_sgh_SF+" based on PD sgh30/sgh ratio."
if(write_landm_coslat.eq.True)then
 fout@landm_coslat = "changed to appropriate topography"
end if

end do

end
    


