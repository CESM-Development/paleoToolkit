; ***********************************************
; 10min2_1x1deg.ncl
; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
;************************************************
; read in data
;************************************************

 ncol = 3

 in1   = addfile("myfile_10min.nc","r")
 ice   = in1->ice
 topo  = in1->htopo
 lmask = in1->landfract
 lfrac = in1->landfract
 lon   = in1->lon
 lat   = in1->lat

 ilon = ispan(0,359,1)
 ilat = fspan(-89.5,89.5,180)
 nlon = dimsizes(ilon)
 nlat = dimsizes(ilat)

 ; create 1 degree dataset for lat/lon

 ; interpolating to 1 degree space
 ice_1   = linint2_Wrap(  ice&lon,  ice&lat,  ice,True,ilon,ilat,0) 
 lmask_1 = linint2_Wrap(lmask&lon,lmask&lat,lmask,True,ilon,ilat,0) 
 lfrac_1 = linint2_Wrap(lfrac&lon,lfrac&lat,lfrac,True,ilon,ilat,0) 
 topo_1  = linint2_Wrap( topo&lon, topo&lat, topo,True,ilon,ilat,0) 

 ; -----------------------------------------------------------------------------
 ; linint places missing values near poles (b/c the grids are interpolated and the poles have no values).
 ; This messes up the definesurf-tool routine and the PS, PHIS, and SGH, adding spurious structure.
 ; So I fill the missing values with ocn in NH; fill with land in SH.
 ; -----------------------------------------------------------------------------
 ; -----------------------------------------------------------------------------
 ; NH:
 ; -----------------------------------------------------------------------------
 ; replace Arctic missing values with ocn

 maxlat = dimsizes(ice&lat)-1
 ice(    maxlat-2:maxlat,:) = 0       ; no land

 ; -----------------------------------------------------------------------------
 ; SH
 ; -----------------------------------------------------------------------------

 ice(    0:2,:) = 100       ; all ice

 omask = floattoint(lmask_1)
 omask@longname  = "landmask"
 topo_1@longname = "topography"
 ice_1@longname = "land_ice"
 lfrac_1@longname = "land_fraction"
   omask!0        = "lat"
  topo_1!0        = "lat"
   ice_1!0        = "lat"
 lfrac_1!0        = "lat"
   omask&lat      = ilat
  topo_1&lat      = ilat
   ice_1&lat      = ilat
 lfrac_1&lat      = ilat
   omask!1        = "lon"
  topo_1!1        = "lon"
   ice_1!1        = "lon"
 lfrac_1!1        = "lon"
   omask&lon      = ilon
  topo_1&lon      = ilon
   ice_1&lon      = ilon
 lfrac_1&lon      = ilon

 dd = systemfunc("date -u +%y%m%d")
 ofile1 = "topo_1deg."+dd+".nc"
 system("/bin/rm -f " + ofile1)
 fout1 = addfile(ofile1,"c")
 fileAtt             = True
 fileAtt@create_date = systemfunc("date")
 fileAtt@conventions = "CCSM"
 fileAtt@author      = "nanr@ucar.edu"
 fileAtt@srcCode     = "10min2_1x1deg.ncl"
 fileAtt@project     = "CESM1"
 fileAtt@in2         = "biome_veg_v1.2.nc (land ice)"
 fileAtt@in1         = "topo_v1.4.nc"
 fileAtt@history     = "topo_v1.4.nc diffs added to USGS-gtopo30_10min_c050419.nc"
 fileattdef(fout1, fileAtt)
   if(any(ismissing(lfrac_1))) then
       print("Your data contains some missing values. Beware.")
   end if


 fout1->TOP       = topo_1
 fout1->ICE       = ice_1
 fout1->LANDMASK  = omask
 fout1->LANDFRAC  = lfrac_1

end
