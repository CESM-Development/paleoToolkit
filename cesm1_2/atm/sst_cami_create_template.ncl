;====================================================================
; reads a CAM annual cycle climo file, modifies a selected variable,
; and writes out a new CCSM sst climo file with mods included 
; see string "MODIFICATIONS", to find where to insert your modfications
; cshields, sep 03
; mods for interpolated sst file to t170, sep 07
; updated for fv0.47x0.63 petm, sep 2010

;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl";


;==========================================================

begin

new_cdf = "./sst_<out_res>_<casename>_c<date>.nc"  ; name for new .nc file
cam_cdf = "./jan-dec.nc"
;fv_cdf = "/glade/p/cesm/cseg/inputdata/ocn/docn7/SSTDATA/sst_HadOIBl_bc_0.23x0.31_clim_c061106.nc"
;fv_cdf = "/glade/p/cesm/cseg/inputdata/ocn/docn7/SSTDATA/sst_HadOIBl_bc_0.47x0.63_clim_c061106.nc
;fv_cdf = "/glade/p/cesm/cseg/inputdata/ocn/docn7/SSTDATA/sst_HadOIBl_bc_0.9x1.25_clim_c040926.nc"
;fv_cdf = "/glade/p/cesm/cseg/inputdata/ocn/docn7/SSTDATA/sst_HadOIBl_bc_1.9x2.5_clim_c061031.nc"
modtext =  "SST/ICEFRAC interplolate up to fv<out_res> from fv<in_res> yrs <years_of_model_run> of <fully_coupled_run> "

; =====================> read netcdf files

localfile = cam_cdf
fvfile  = fv_cdf
print("Reading "+ cam_cdf)
print("Reading "+ fv_cdf)

;==================> get dimension sizes from files and read data

fv  = addfile(fvfile,"r")
lat = fv->lat
lon = fv->lon
date = fv->date
datesec = fv->datesec
time = fv->time
ice_cov = fv->ice_cov
ice_cov_prediddle = fv->ice_cov_prediddle
SST_cpl = fv->SST_cpl
SST_cpl_prediddle = fv->SST_cpl_prediddle

nlon = dimsizes(lon)
nlat = dimsizes(lat)
ntimes= dimsizes(time)

cam  = addfile(localfile,"r")
latt31 = cam->lat
lont31 = cam->lon
ts      = cam->TS
icefrac = cam->ICEFRAC

nlont31  = dimsizes(lont31)
nlatt31  = dimsizes(latt31)

print("nlon = " + nlont31)
print("nlat = " + nlatt31)
print("ntimes = " + ntimes)

;============= make changes to var_to_change =================
; MAKE YOUR MODIFICATIONS IN THIS SECTION

; MODIFICATIONS GO HERE
newsst = linint2 (lont31,latt31,ts, True, lon,lat, 0)
newice = linint2 (lont31,latt31,icefrac, True, lon,lat, 0)

;--HARDWIRE, for t31->T170: fix rows 0-3 and nlon-4 from missing to row 4  and nlon-5
; otherwise missing data and model will croak
;print("Fixing inherent interpolation error from t31 to t170...")
;do n = 0,3
 ;newsst(:,n,:) = newsst(:,4,:)
 ;newice(:,n,:) = newice(:,4,:)
;end do
;do n = nlat-4,nlat-1
 ;newsst(:,n,:) = newsst(:,nlat-5,:)
 ;newice(:,n,:) = newice(:,nlat-5,:)
;end do

;--HARDWIRE, for t31->f05 fix rows 0-6 and nlat-6 from missing to row 7  and nlat-7
; otherwise missing data and model will croak
;print("Fixing inherent interpolation error from t31 to f05...")
;do n = 0,6
; newsst(:,n,:) = newsst(:,7,:)
; newice(:,n,:) = newice(:,7,:)
;end do
;do n = nlat-7,nlat-1
; newsst(:,n,:) = newsst(:,nlat-8,:)
; newice(:,n,:) = newice(:,nlat-8,:)
;end do

;assign coords
newsst!0 = "time"
newsst!1 = "lat"
newsst!2 = "lon"
newsst&lat = lat
newsst&lon = lon
newsst&time = time
newice!0 = "time"
newice!1 = "lat"
newice!2 = "lon"
newice&lat = lat
newice&lon = lon
newice&time = time


;---convert to degC
newsst = newsst - 273.15

; make sure original attributes preserved 
copy_VarAtts(SST_cpl,newsst)
copy_VarAtts(ice_cov,newice)

; just copy to prediddle..not sure it matters for PETM...talk to john truesdale?
newsst_pred = newsst
newice_pred = newice

print("Changes complete")
printVarSummary(newsst)
printVarSummary(newice)


;========create new file surface file, preserve old file
;        attributes, write out variables with no changes

;----> write new .nc file

system("rm "+new_cdf)
cdf = addfile(new_cdf, "c")

; define dimensions
dimNames = (/"lat","lon","time"/)
dimSizes = (/ nlat, nlon, -1/)
dimUnlim = (/ False,  False, True/) 
filedimdef(cdf,dimNames,dimSizes,dimUnlim)

;define variable names
filevardef(cdf,"lat","double",(/"lat"/))
filevardef(cdf,"lon","double",(/"lon"/))
filevardef(cdf,"time","double",(/"time"/))
filevardef(cdf,"date","integer",(/"time"/))
filevardef(cdf,"datesec","integer",(/"time"/))
filevardef(cdf,"SST_cpl","float",(/"time","lat","lon"/))
filevardef(cdf,"ice_cov","float",(/"time","lat","lon"/))
filevardef(cdf,"SST_cpl_prediddle","float",(/"time","lat","lon"/))
filevardef(cdf,"ice_cov_prediddle","float",(/"time","lat","lon"/))

; write variables
cdf->lat    = lat
cdf->lon    = lon
cdf->time   = time 

cdf->date   = date 
cdf->datesec= datesec 
cdf->SST_cpl = newsst 
cdf->ice_cov = newice 
cdf->SST_cpl_prediddle = newsst_pred 
cdf->ice_cov_prediddle = newice_pred 


; preserve file attributes
cdf@history_t31_sst = cam@history
cdf@history_f05_template = fv@history

; add new text
cdf@Modifications = modtext
cdf@nclcode = "sst_cami_create_f05.ncl" 

print(new_cdf + " has been written")


end 
