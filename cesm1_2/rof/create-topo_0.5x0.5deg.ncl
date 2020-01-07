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
filo = "topo.0.5x0.5deg."+case+".nc"

fin  = addfile(fili,"r")
system("/bin/rm -f " + filo)
fout = addfile(filo,"c")

fmsg = default_fillvalue("float")

ilat = fin->lat
ilon = fin->lon
topf = fin->topo
topf@_FillValue = fmsg

 nlat = 360
 nlon = 720
 lon05deg = new((/nlon/),"float",fmsg)
 lat05deg = new((/nlat/),"float",fmsg)

; =====================================
 ; create lon
 sl = 0.25
 do i = 0,nlon-1
        lon05deg(i) = sl 
	sl = sl + 0.5
  end do

  ; create lat
  sl = -89.75
  do i = 0,nlat-1
        lat05deg(i) = sl
        sl = sl + 0.5
  end do

lon05deg!0 = "lon"
lat05deg!0 = "lat"

lon05deg@units = "degrees_east"
lat05deg@units = "degrees_north"

lon05deg@long_name = "longitude"
lat05deg@long_name = "latitude"


; =====================================

 top05 = new((/dimsizes(lat05deg), dimsizes(lon05deg)/),"float",fmsg)

 top05!0   = "lat"
 top05&lat = lat05deg
 top05!1   = "lon"
 top05&lon = lon05deg

; =====================================

  top05 = linint2_Wrap (topf&lon,topf&lat,topf, True, lon05deg,lat05deg, 0)
  top05@long_name = "permian topo/bathy 0.5x0.5 resolution"
  top05@_FillValue = fmsg

  ; SH
  ll = 4
  if(any(ismissing(top05(ll,:)))) then

		latfill = ll+1
		top05(ll,:) = top05(latfill,:)
		ll = ll -1
        else if(any(ismissing(top05(ll,:)))) then
		latfill = ll+1
		top05(ll,:) = top05(latfill,:)
		ll = ll - 2
        else if(any(ismissing(top05(ll,:)))) then
		latfill = ll+1
		top05(ll,:) = top05(latfill,:)
		ll = ll - 3
        else if(any(ismissing(top05(ll,:)))) then
		latfill = ll+1
		top05(ll,:) = top05(latfill,:)
		ll = ll - 4
        else if(any(ismissing(top05(ll,:)))) then
		latfill = ll+1
		top05(ll,:) = top05(latfill,:)
  end if
  end if
  end if
  end if
  end if

  ; NH
  NH1 = dimsizes(lat05deg)-5
  NH2 = dimsizes(lat05deg)-1
  do ll = NH1,NH2
      if(any(ismissing(top05(ll,:)))) then
	latfill = ll-1
	top05(ll,:) = top05(latfill,:)
      end if
  end do


  fileAtt               = True
  fileAtt@creation_date = systemfunc("date")
  fileAtt@input         = "bathtopo-fin.grd.nc"
  fileAtt@src           = "create0.5degree.ncl"
  fileAtt@contents      = "permian topo + bathymetry 0.5 x 0.5 degree"
  fileattdef(fout, fileAtt)

  fout->topo = top05

  delete(filo)
  delete(fin)
end
    


