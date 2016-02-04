;*************************************************
; NCL Graphics: polar_5.ncl
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
;************************************************
; read in netCDF file and access data
;************************************************

  ; open output file.
  dd = systemfunc("date -u +%y%m%d")

  ifile  = "topo_myrun_remap_09x25.sm9.DATE.nc"
  ofile  = "topo_myrun_remap_09x25.mod.sm9.DATE.nc"
  in = addfile(ifile,"r")
  out = addfile(ofile,"c")

  landm = in->LANDM_COSLAT
  landm_coslat = in->LANDM_COSLAT
  lat = in->lat
  lon = in->lon
  landmask = in->LANDFRAC
  ; landmask = landmask * 0.01		; convert to fraction for mkgriddata - nanr 8oct10

  ; Bering strait 
  ; nloni = noi
  ; xloni = xoi
  ; nlatj = naj
  ; xlatj = xaj
  noi = 134
  xoi = 165
  naj = 152
  xaj = 176
  landm(naj:xaj,noi:xoi) = where(landmask(naj:xaj,noi:xoi) .gt. 0, 1, landm(naj:xaj,noi:xoi))

  ; Hudson Bay
  noi = 208
  xoi = 231
  naj = 148
  xaj = 170
  landm(naj:xaj,noi:xoi) = where(landmask(naj:xaj,noi:xoi) .gt. 0, 1, landm(naj:xaj,noi:xoi))

  ; FennoScandia
  noi = 0
  xoi = 134
  naj = 158
  xaj = 185
  landm(naj:xaj,noi:xoi) = where(landmask(naj:xaj,noi:xoi) .gt. 0, 1, landm(naj:xaj,noi:xoi))

  ; Australia
  noi = 100
  xoi = 120
  naj = 78
  xaj = 90
  landm(naj:xaj,noi:xoi) = where(landmask(naj:xaj,noi:xoi) .gt. 0, 1, landm(naj:xaj,noi:xoi))

  ; Indonesia
  noi = 75
  xoi = 95
  naj = 83
  xaj = 115
  landm(naj:xaj,noi:xoi) = where(landmask(naj:xaj,noi:xoi) .gt. 0, 1, landm(naj:xaj,noi:xoi))

; ==============================================
; == Read in temporary file and copy all variables to 
;    file.  Substitute in new landm.
; ==============================================
 system ("rm "+ ofile)             ; remove any pre-existing file
 ncout = addfile(ofile,"c")

 att_names = getvaratts(in)           ; get CCM file's global attributes
 if(.not.all(ismissing(att_names))) then
    do i = 0,dimsizes(att_names)-1
      print("copy_fileatts: global attributes->" + att_names(i))
      ncout@$att_names(i)$ = in@$att_names(i)$     ; copy file's global attributes
    end do
  end if

  fileAtt             = True
  fileAtt@title       = "Landm_coslat modified for LGM21ka"
  fileAtt@create_date = systemfunc("date")
  fileAtt@src_code = "mod_landm_coslat.ncl"
  fileAtt@author = systemfunc("echo $USER")
  fileAtt@input1 = ifile
  fileattdef(ncout,fileAtt)
 
  names = getfilevarnames(in)       ; get CCM file's variable names and print them
  ; print (names)
;
; loop over variables and copy to netCDF file

  ; degrad = 0.0174533  ; conversion of degrees to radians
  ; latr   = lat*degrad 
  ; latrC  = conform(landm,latr,0)
  ; landm  = where(landm_coslat.lt.0.999, landm_coslat/cos(latrC),landm_coslat)

  do i=0, dimsizes(names)-1
    print("writing ncdf: i,name="+ i+ "   "+ names(i))
       if (names(i) .eq. "LANDM_COSLAT") then
            ncout->$names(i)$ = landm                          ; new LGM landm_coslat 
            ncout->LANDM_COSLAT_ORIG = in->$names(i)$ 
       else
            if (names(i) .eq. "LANDFRAC") then
                 ncout->$names(i)$ = landmask                  ; output in 0-1 for mkgriddata
	    else
                 ncout->$names(i)$ = in->$names(i)$                  ; output in 0-1 for mkgriddata
            end if
       end if
  end do

  print (ncout)      ; print overview of contents of netCDF file

end

