;=========================================================; 
; 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; =================
; User defined
; =================

case = "<casename>"
fili = "<input_topo-bath_filename>"
filo = "topo.1x1deg."+case+".nc"

; =================
; =================

fin  = addfile(fili,"r")
system("/bin/rm -f " + filo)
fout = addfile(filo,"c")

ilat = fin->lat
ilon = fin->lon
top = fin->topo
fmsg = default_fillvalue("float")


 nlat = 180
 nlon = 360
 lon1x1deg = new((/nlon/),"float",fmsg)
 lat1x1deg = new((/nlat/),"float",fmsg)

; =====================================
 ; create lon
 sl = 0.5
 do i = 0,nlon-1
        lon1x1deg(i) = sl 
	sl = sl + 1.0
  end do

  ; create lat
  sl = -89.5
  do i = 0,nlat-1
        lat1x1deg(i) = sl
        sl = sl + 1.0
  end do

lon1x1deg!0 = "lon"
lat1x1deg!0 = "lat"

lon1x1deg@units = "degrees_east"
lat1x1deg@units = "degrees_north"

lon1x1deg@long_name = "longitude"
lat1x1deg@long_name = "latitude"

; =====================================

 top1x1 = new((/dimsizes(lat1x1deg), dimsizes(lon1x1deg)/),"float",fmsg)

 top1x1!0   = "lat"
 top1x1&lat = lat1x1deg
 top1x1!1   = "lon"
 top1x1&lon = lon1x1deg

; =====================================

  top1x1 = linint2_Wrap (top&lon,top&lat,top, True, lon1x1deg,lat1x1deg, 0)
  top1x1@long_name = case + " topography 1x1 resolution"


  fileAtt               = True
  fileAtt@creation_date = systemfunc("date")
  fileAtt@input         = fili
  fileAtt@src           = "create-topo.1x1.ncl"
  fileAtt@contents      = case + " topography 1x1 degree"
  fileattdef(fout, fileAtt)

  fout->topo = top1x1

  delete(filo)
  delete(fin)
end
    


