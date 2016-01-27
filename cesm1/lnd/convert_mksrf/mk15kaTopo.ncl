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


 ifile1 = "ice5g_v1.2_15.0k_1deg.nc"		; Peltier ice5g 15ka topo
 in1   = addfile(ipath + ifile1+".nc","r")
 topo15ka = in1->orog
 ice15ka  = in1->sftgif
 dlat = in1->lat
 dlon = in1->long

 ifile2 = "ice5g_v1.2_00.0k_1deg.nc" 		; Peltier ice5g 0ka topo
 in2   = addfile(ipath + ifile2+".nc","r")
 topo0ka = in2->orog
 ice0ka  = in2->sftgif
 dlat = in2->lat
 dlon = in2->long

 ; create delT from the Peltier ICE5g 15ka-0ka, then add the delta
 ; to the USGS control orography.

 delT = topo15ka 
 delT = 0
 delT = topo15ka - topo0ka

 deli = ice15ka
 deli = 0
 deli = ice15ka - ice0ka

 lfrac = topo15ka
 lfrac = 0
 lfrac = where(topo15ka .gt. 0, 1, 0)

 ifile3 = "/fs/cgd/csm/inputdata/atm/cam/topo/USGS-gtopo30_10min_c050419.nc"
 in3   = addfile(ifile3+".nc","r")
 htopo = in3->htopo
 hlat  = in3->lat
 hlon  = in3->lon
 ivar= in3->variance
 nhlat = dimsizes(hlat)
 nhlon = dimsizes(hlon)

 hmask = htopo
 hmask = 0
 hmask = where(htopo .gt. 0, 1, 0)

 ; interpolate the PMIP3 deltas to 10min, b/c we need a 10min grid for definesurf:
 
 ilfrac = htopo
 ilfrac = 0

 idelT  = linint2_Wrap(delT&long,delT&lat,delT,True, hlon,hlat, 0) 	; resample to 10min
 ideli  = linint2_Wrap(deli&long,deli&lat,deli,True, hlon,hlat, 0) 	; resample to 10min
 ilfrac = linint2_Wrap(lfrac&long,lfrac&lat,lfrac,True, hlon,hlat, 0)    ; resample to 10min
 iice15ka = linint2_Wrap(ice15ka&long,ice15ka&lat,ice15ka,True, hlon,hlat, 0)   	; resample to 10min

 ilfrac@_FillVal  = -1e34
 iice15ka@_FillValue  = -1e34
 idelT@_FillValue = -1e34
 ideli@_FillValue = -1e34

 ; -----------------------------------------------------------------------------
 ; linint places missing values near poles (b/c the grids are interpolated and the poles have no values).
 ; This messes up the definesurf-tool routine and the PS, PHIS, and SGH, adding spurious structure.
 ; So I fill the missing values with ocn in NH; fill with land in SH.
 ; -----------------------------------------------------------------------------
 ; -----------------------------------------------------------------------------
 ; NH:
 ; -----------------------------------------------------------------------------
 ; replace Arctic missing values with ocn
 
 maxlat = dimsizes(iice15ka&lat)-1 
 ilfrac(maxlat-2:maxlat,:) = 0       ; no land
 iice15ka(  maxlat-2:maxlat,:) = 0       ; no land
 idelT( maxlat-2:maxlat,:) = 0       ; no land
 ideli( maxlat-2:maxlat,:) = 0       ; no land

 ; -----------------------------------------------------------------------------
 ; SH
 ; -----------------------------------------------------------------------------

 ilfrac(0:2,:) = 1         ; all land
 iice15ka(0:2,:) = 100       ; all land
 idelT(0,:) = idelT(3,:)                                                    ; nearest neighbor
 idelT(1,:) = idelT(3,:)                                                    ; nearest neighbor
 idelT(2,:) = idelT(3,:)                                                    ; nearest neighbor

 ideli(0,:) = ideli(3,:)                                                    ; nearest neighbor
 ideli(1,:) = ideli(3,:)                                                    ; nearest neighbor
 ideli(2,:) = ideli(3,:)                                                    ; nearest neighbor

 ; -----------------------------------------------------------------------------
 ; check for missing values
 ; -----------------------------------------------------------------------------
 if(any(ismissing(ilfrac))) then
    print("lmask contains some missing values. Beware.")
 end if
 if(any(ismissing(iice15ka))) then
    print("ice contains some missing values. Beware.")
 end if
 if(any(ismissing(idelT))) then
    print("idelT contains some missing values. Beware.")
 end if
 if(any(ismissing(ideli))) then
    print("ideli contains some missing values. Beware.")
 end if

 ; add topography delta to NEW ice areas (areas with 15ka delta in ice thickness)
 ; This does not include the himalaya.
 topo15ka10min = htopo
 topo15ka10min = 0.
 topo15ka10min = where(iice15ka .eq. 100, htopo + idelT, htopo)

 iceout = iice15ka 

 domask = hmask
 ; fill in ocn where new glaciers extend beyond continent
; domask(540:1043,:) = where(iice15ka(540:1043,:) .eq. 100, 1, hmask(540:1043,:))

 domask = where (topo15ka10min .gt. 0, 1, domask)


 domask@longname = "dlandmask"
 domask!0	= "lat"
 domask&lat	= hlat
 domask!1	= "lon"
 domask&lon	= hlon
 domask@_FillValue = -9999
 domask@valid_max = 1
 domask@valid_min = 0

 dd = systemfunc("date -u +%y%m%d")

 ofile1 = "topo_15ka_10min."+dd+".nc"
 system("/bin/rm -f " + ofile1)

 fout1 = addfile(ofile1,"c")
 fileAtt             = True
 fileAtt@create_date = systemfunc("date")
 fileAtt@author      = "nanr@ucar.edu"
 fileAtt@srcCode     = "mk15kaTopo.ncl"
 fileAtt@infile1     = ifile1
 fileAtt@infile2     = ifile2
 fileAtt@infile3     = ifile3
 fileattdef(fout1, fileAtt)

 fout1->htopo      = topo15ka10min
 fout1->ice        = iceout
 fout1->landfract  = domask
 fout1->landmask   = domask
 fout1->variance   = ivar

;*************************************************
; create plot
;*************************************************
  ofile2 = "topo_15ka_10min."+dd
  ; ofile2 = "topo_10min.lmask."+tol+"."+dd
  wks  = gsn_open_wks("ps",ofile2)     ; open a ps file
  gsn_define_colormap(wks,"BlAqGrYeOrRe")     ; choose colormap

  res                      = True     ; plot mods desired

  res@cnFillOn             = True     ; turn on color fill
  res@cnFillMode           = "RasterFill"     
  res@cnLineLabelsOn       = False    ; turn off line labels
  res@cnLinesOn            = False    ; turn of contour lines
  res@cnMissingValPerimOn  = True     ; outline missing values
  res@cnInfoLabelOn        = False    ; turn off info label
  res@tmLabelAutoStride    = True     ; nice stride on labels
  res@lbLabelAutoStride    = True     ; nice label bar labels
  ;res@gsnSpreadColors      = True     ; use full colormap
  res@gsnMaximize          = True
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnFillColors                = (/"white","blue"/)
  ; res@cnFillColors                = (/"white","white","blue"/)
  res@cnFillDrawOrder             = "PreDraw"
  res@cnExplicitLegendLabelsOn    = True
  res@cnLevels = (/50/)
  ; res@cnLevels = (/0,50/)
  res@lbLabelStrings        = " " 
  ; res@lbLabelStrings        = (/"ocean","land"/)
  res@lbLabelStride         = 1
  res@cnInfoLabelOn = False 
  res@tiMainString   = "10min interpolated Peltier ice5g deltas (15ka-0ka)"
  res@tiXAxisString   = ""
  res@tiYAxisString   = ""
  res@gsnRightString  = " "
  res@gsnLeftString = " "
 
  plot = gsn_csm_contour(wks,domask,res)    ; contour the variable

end
