;=========================================================; 
; 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

case = "<casename>"
fili = "<input_topo-bath_filename>"
filo = "topo.2x2deg."+case+".nc"

fin  = addfile(fili,"r")
system("/bin/rm -f " + filo)
fout = addfile(filo,"c")

ilat = fin->lat
ilon = fin->lon
top = fin->topo
fmsg = default_fillvalue("float")


 nlat = 90
 nlon = 180
 lon2x2deg = new((/nlon/),"float",fmsg)
 lat2x2deg = new((/nlat/),"float",fmsg)

; =====================================
 ; create lon
 sl = 1.0
 do i = 0,nlon-1
        lon2x2deg(i) = sl 
	sl = sl + 2.0
  end do

  ; create lat
  sl = -89
  do i = 0,nlat-1
        lat2x2deg(i) = sl
        sl = sl + 2.0
  end do

lon2x2deg!0 = "lon"
lat2x2deg!0 = "lat"

lon2x2deg@units = "degrees_east"
lat2x2deg@units = "degrees_north"

lon2x2deg@long_name = "longitude"
lat2x2deg@long_name = "latitude"

; =====================================

 top2x2 = new((/dimsizes(lat2x2deg), dimsizes(lon2x2deg)/),"float",fmsg)

 top2x2!0   = "lat"
 top2x2&lat = lat2x2deg
 top2x2!1   = "lon"
 top2x2&lon = lon2x2deg

; =====================================

  top2x2 = linint2_Wrap (top&lon,top&lat,top, True, lon2x2deg,lat2x2deg, 0)
  top2x2@long_name = case + " topography 2x2 resolution"


  fileAtt               = True
  fileAtt@creation_date = systemfunc("date")
  fileAtt@input         = fili
  fileAtt@src           = "create-topo.2x2.ncl"
  fileAtt@contents      = case + " topography 2x2 degree"
  fileattdef(fout, fileAtt)

  fout->topo = top2x2

  delete(filo)
  delete(fin)
end
    


