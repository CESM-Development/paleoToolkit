;=========================================================; 
; reads in data merged topo kmt file
; reads in landfrac from cesm1 surface_data file 
; (run mksurfdata offline first)
; read in template cam bnd_topo (Netcdf) file 
; makes modifications, writes out new file
  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;====================user defined parameters

dd= systemfunc("date -u +%y%m%d")

cases = (/"<casename>"/)

do c = 0,dimsizes(cases)-1

usecase = cases(c)

cesmfiles = "/glade/p/cesm/cseg/inputdata/atm/cam/topo/"
path = "<directory_with_topo-bath_file>"

;=============== The code reads landfrac from paleo surface_data
; so you should have already created this offline using mksurfdata tool.
; =====================================================================

topoinput = "<topo-bath_file>"

i10min = "USGS-gtopo30_10min_c050419.nc"

filename_topo    = path+topoinput
new_cdf1         = usecase+"_10min_topo_4input2definesurf."+dd+".nc"		

;============== read atm topo file

msg = 1.e+36
topof = addfile(filename_topo,"r")
; topf = topof->topo_depth
topf = flt2dble(topof->topo)
ilat = topof->lat
ilon = topof->lon
topf@_FillValue = 0.
topf = mask(topf,topf.le.0,False)
delete(topf@_FillValue)
topf@_FillValue = msg
printVarSummary(topf)
print("Atm data read")

;=============== read in atm template netcdf file
 
; USING linint2: interpolate from one grid to another grid
; using bilinear interpolation. Missing values are allowed.
; produces similar interpolation used by jim's "definesurf"

usetopo = i10min
in3   = addfile(cesmfiles+i10min,"r")
olat  = in3->lat
olon  = in3->lon
nlat = dimsizes(olat) 
nlon = dimsizes(olon) 
print("nlon = " + nlon)
print("nlat = " + nlat)

; interpolate to 10min, b/c we need a 10min grid for definesurf:
htopo  = linint2_Wrap(ilon,ilat,topf,True, olon,olat, 0)      ; resample to 10min
print("Atm topography interpolated to 10min for definesurf")
htopo@long_name = "Topography"

;======= fix any interpolation problems at poles

doFix = True
if (doFix .eq. True) then

if (any(ismissing(htopo(2,:)))) then
	htopo(0,:) = (/htopo(3,:)/)
	htopo(1,:) = (/htopo(3,:)/)
	htopo(2,:) = (/htopo(3,:)/)
end if
if (any(ismissing(htopo(dimsizes(olat)-3,:)))) then
	htopo(dimsizes(olat)-1,:)=(/htopo(dimsizes(olat)-4,:)/)
	htopo(dimsizes(olat)-2,:)=(/htopo(dimsizes(olat)-4,:)/)
	htopo(dimsizes(olat)-3,:)=(/htopo(dimsizes(olat)-4,:)/)
end if
if (any(ismissing(htopo(1,:)))) then
	htopo(0,:) = (/htopo(2,:)/)
	htopo(1,:) = (/htopo(2,:)/)
end if
if (any(ismissing(htopo(dimsizes(olat)-2,:)))) then
	htopo(dimsizes(olat)-1,:)=(/htopo(dimsizes(olat)-3,:)/)
	htopo(dimsizes(olat)-2,:)=(/htopo(dimsizes(olat)-3,:)/)
end if
if (any(ismissing(htopo(0,:)))) then
	htopo(0,:) = (/htopo(1,:)/)
end if
if (any(ismissing(htopo(dimsizes(olat)-1,:)))) then
	htopo(dimsizes(olat)-1,:)=(/htopo(dimsizes(olat)-2,:)/)
end if

end if ; doFix loop
if (any(ismissing(htopo))) then
	print("FATAL:  Missing data in htopo; check the poles")
	exit
end if


;======= create landfrac field 
htopo@_FillValue = msg
ftopo = htopo		; create landfrac (ftopo) variable
ftopo = where(htopo .gt. 0, 1, 0)
ftopo@long_name = "landfrac"

landm_coslat = htopo

;========create new file surface file

;----> write new .nc file
dimNames = (/"lat","lon"/)
dimSizes = (/ nlat, nlon/)
dimUnlim = (/ False,  False/)

; define dimensions
system("rm "+new_cdf1)
cdf1 = addfile(new_cdf1, "c")
print("Creating new landm_coslat file....")
filedimdef(cdf1,dimNames,dimSizes,dimUnlim)

;define variable names
filevardef(cdf1,"lat","double",(/"lat"/))
filevardef(cdf1,"lon","double",(/"lon"/))
filevardef(cdf1,"htopo","double",(/"lat","lon"/))
filevardef(cdf1,"ftopo","double",(/"lat","lon"/))

delete(htopo@_FillValue)
delete(ftopo@_FillValue)
htopo@_FillValue = msg
ftopo@_FillValue = msg

; write variables
cdf1->lat    = olat
cdf1->lon    = olon
cdf1->htopo  = htopo
cdf1->ftopo  = ftopo
print(new_cdf1 + " has been written")

; create file attributes
date_written = systemfunc("date -u +%y%m%d")
user  = systemfunc("whoami")
scriptFile = get_script_prefix_name()
cdf1@history = "created by " + user + " with "+scriptFile+ " on " + date_written 
cdf1@template = "template = " + usetopo 
cdf1@paleotopo = "paleo topo = " +  filename_topo 

end do		; case loop

end
    


