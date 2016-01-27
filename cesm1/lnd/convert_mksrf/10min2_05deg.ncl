; ***********************************************
; mkIceSheet_aixue.ncl
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
 ipath = "./"


 ifile1 = "ice5g_v1.2_15.0k_1deg.nc"            ; Peltier ice5g 15ka topo
 in1   = addfile(ipath + ifile1,"r")
 topo15ka = in1->orog
 ice  = in1->sftgif
 dlat = in1->lat
 dlon = in1->long

 ; in1   = addfile(ipath + "ice5g_v1.2_00.0k_1deg.nc","r")
 ; ice   = in1->sftgif				; glacier fraction
 ; ilat  = in1->lat
 ; ilon  = in1->long

 in2   = addfile(ipath + "topo_15ka_10min.140226.nc","r")
 topo  = in2->htopo				; topo
 lmask = in2->landmask				; land fraction (0-1)
 tlat  = in2->lat
 tlon  = in2->lon


 ; Read 05 degree dataset for lat/lon
 gpath = "/glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/"
 in4 = addfile(gpath + "mksrf_glacier.060929.nc","r")
 icemksrf = in4->PCT_GLACIER
 lat05    = in4->LAT
 lon05    = in4->LON
 icemksrf!0   = "lat"
 icemksrf&lat = lat05
 icemksrf!1   = "lon"
 icemksrf&lon = lon05
 icemksrf = lonFlip(icemksrf)

 ; interpolating to 05 degree space
 ice_05   = linint2_Wrap(ice&long,ice&lat,ice,True,icemksrf&lon,icemksrf&lat,0) 
 lmask_05 = linint2_Wrap(lmask&lon,lmask&lat,lmask,True,icemksrf&lon,icemksrf&lat,0) 
 topo_05  = linint2_Wrap(topo&lon,topo&lat,topo,True,icemksrf&lon,icemksrf&lat,0) 

 dd = systemfunc("date -u +%y%m%d")
 ofile1 = "topo-ice.05.15ka."+dd+".nc"
 system("/bin/rm -f " + ofile1)
 fout1 = addfile(ofile1,"c")
 fileAtt             = True
 fileAtt@create_date = systemfunc("date")
 fileAtt@conventions = "CCSM"
 fileAtt@author      = "nanr@ucar.edu"
 fileAtt@srcCode     = "10min2_05deg.ncl"
 fileAtt@project     = "Aixue Hu 15ka"
 fileattdef(fout1, fileAtt)


 ; -----------------------------------------------------------------------------
 ; linint places missing values near poles b/c the grids are offset by 0.25 deg
 ; fill with ocn in NH; fill with land in SH.
 ; -----------------------------------------------------------------------------
 ; first create 2D lat array to use in where statement.
 uselat = lmask_05&lat
 nlat = dimsizes(lmask_05&lat)
 nlon = dimsizes(lmask_05&lon)
 USELAT = conform_dims( (/nlat,nlon/), uselat, 0)

 ; -----------------------------------------------------------------------------
 ; NH:
 ; -----------------------------------------------------------------------------
 ; replace Arctic missing values with ocn
 lmask_05 = where((ismissing(lmask_05) .and. USELAT .gt. 0), 0.,lmask_05)	; no land
 ice_05   = where((ismissing(ice_05)   .and. USELAT .gt. 0), 0.,ice_05)		; no land-ice
 topo_05  = where((ismissing(topo_05)  .and. USELAT .gt. 0), 0.,topo_05)	; no topography

 ; -----------------------------------------------------------------------------
 ; SH
 ; -----------------------------------------------------------------------------
 lmask_05 = where((ismissing(lmask_05) .and. USELAT .lt. 0), 1.,lmask_05)	; landmask = 1
 ice_05 = where((ismissing(ice_05) .and. USELAT .lt. 0), 100.,ice_05)	        ; land-ice = 1
 topo_05(0,:) = topo_05(1,:) 							; nearest neighbor

 ; -----------------------------------------------------------------------------
 ; linint places missing values near poles (b/c the grids are interpolated and the poles have no values).
 ; This messes up the definesurf-tool routine and the PS, PHIS, and SGH, adding spurious structure.
 ; So I fill the missing values with ocn in NH; fill with land in SH.
 ; -----------------------------------------------------------------------------
 ; -----------------------------------------------------------------------------
 ; NH:
 ; -----------------------------------------------------------------------------
 ; replace Arctic missing values with ocn

 maxlat = dimsizes(lmask_05&lat)-1
 lmask_05(maxlat-2:maxlat,:) = 0       ; no land
 ice_05(  maxlat-2:maxlat,:) = 0       ; no land
 topo_05( maxlat-2:maxlat,:) = 0       ; no land

 ; -----------------------------------------------------------------------------
 ; SH
 ; -----------------------------------------------------------------------------

 lmask_05(0:2,:) = 1       ; all land
 ice_05(  0:2,:) = 100     ; all ice
 topo_05(0,:) = topo_05(3,:)                                                    ; nearest neighbor
 topo_05(1,:) = topo_05(3,:)                                                    ; nearest neighbor
 topo_05(2,:) = topo_05(3,:)                                                    ; nearest neighbor

 ; -----------------------------------------------------------------------------
 ; Now check for missing points in all fields.
 ; -----------------------------------------------------------------------------
  if(any(ismissing(lmask_05))) then
    print("lmask contains some missing values. Beware.")
  end if
  if(any(ismissing(topo_05))) then
    print("topo_05 contains some missing values. Beware.")
  end if
  if(any(ismissing(ice_05))) then
    print("ice_05 contains some missing values. Beware.")
  end if

  newlmask_05 = lmask_05
  newlmask_05 = 0
  newlmask_05 = where(lmask_05 .gt. 0, 1, 0)

  lfrac_05 = lmask_05
  lfrac_05 = lfrac_05 

  printVarSummary(lfrac_05)
  printVarSummary(newlmask_05)

  fout1->TOP       = topo_05
  fout1->ICE       = ice_05
  fout1->LANDFRAC  = lfrac_05
  fout1->LANDMASK  = newlmask_05

end
