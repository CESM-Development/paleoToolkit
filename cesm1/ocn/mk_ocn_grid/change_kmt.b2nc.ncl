load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"   

; ----------------------------------------------------
; usage:  ncl change_kmt.b2nc.ncl
; Description: Read binary direct access file and rewrite to netcdf
; also remaps ocn grid to lat/lon using user defined binary mapping file
; required: mapping file and ocn domain output from scrip (mk_SCRIPgrid.csh).
; ----------------------------------------------------

begin
 print ("=========================================")
 print ("Starting:  change_kmt.b2nc.ncl")
 print ("Start Time: "+systemfunc("date") )
 print ("=========================================")

 dd= systemfunc("date -u +%y%m%d")

 case  = "PT"

; ==============================================
; ==============================================
; REMAP pop grid  to lat/lon topo file 
; Remapping requires a coupler mapping file:: 
;	e.g., map_gx1PT_TO_fv19_25_blin.{DATE}.nc
; The map file name must follow the exact syntax:  map_XXXX_to_XXXXX_bilin_da_DATE.nc
; However, this syntax differs from that produced by the cesm1_2 mapping tools.
; Therefore, you need to rename, or create a softlink to, 
; a mapping filename that matches the syntax shown above.
; e.g., 
; > ln -s /MYPATH/cpl_mapping/maps/map_gx1PT_TO_fv19_25_blin.{DATE}.nc ./map_gx1PT_to_fv19_25_bilin_da_{DATE}.nc
; ==============================================
; ==============================================

 REMAP = False

 OUTPATH  = "./permian/"
 TOPOPATH = "/MYPATH/topobath/"

 ifile1    = OUTPATH + "/kmt.4.da"	     		     ; kmt in bin format
 ifile2    = OUTPATH + "/grid.4.pop.da"	     ; region_mask in bin format
 ifile3    = TOPOPATH + "/myTopo.1deg.DATE.nc"
 ifile4    = OUTPATH + "/gx1PT_DATE.nc"
 ofile     = "kmt."+dd+".nc"
 ofilebin1 = "kmt."+dd+".ieeei4"

; =========================================
; Warning: ================================
; getting corners from netcdf file;  Need to be sure that grid.1.pop.da 
; is consistent with netcdf file.
; =========================================

 vtype = "integer"
 dtype = "double"
 npts  = 0

 PI = 4.0 * atan(1.0)
 RadToDeg = 180.0/PI

; NOTE lat and lon are hardcoded to regular POP Grid
 nlon = 320
 nlat = 384 + npts
 dims = (/nlat,nlon/)

; ==============================================
; == Read in kmt from topography file
; ==============================================
; NOTE:  fortran is column major .... NCL [C] is row major
; fortran(nlon,nlat)  <===> NCL(nlat,nlon)
; ==============================================

 ; ===================================================================
 ; NOTE: kmt and grid files must be bigEndian
 ; ===================================================================

  setfileoption ("bin", "ReadByteOrder", "BigEndian")
  setfileoption ("bin", "WriteByteOrder","BigEndian")
  print("Big Endian")

 ; ===================================================================
 ; ===================================================================
 
 kmt    = fbindirread(ifile1,0,dims,vtype)
 ulati  = fbindirread(ifile2,0,dims,dtype)
 uloni  = fbindirread(ifile2,1,dims,dtype)
 ulati  = ulati*RadToDeg
 uloni  = uloni*RadToDeg
 ulat   = ulati
 ulon   = uloni
 ulat0  = ulati(:,0)
 ulon0  = uloni(0,:)
 spolelat = ulat0(0)
 spolelon = ulon0(0)
 npolelat = ulat0(dimsizes(ulat0) - 1)
 npolelon = ulon0(dimsizes(ulon0) - 1)
print("NH lat = " + npolelat)
print("NH lon = " + npolelon)
print("SH lat = " + spolelat)
print("SH lon = " + spolelon)

if (REMAP .eq. True)	; ----------------->


  dkmt = kmt * 1.d
  rkmt=PopLatLon(dkmt,"gx1PT","fv19_25","bilin","da","130909")

end if 		; REMAP loop -------------->

; ==============================================
; Add coordinate vars
; ==============================================

  kmt@units = "unitless"
  ulat@units = "degrees_north"
  ulon@units = "degrees_east"

  kmt@lat2d = ulati
  kmt@lon2d = uloni 

; ==============================================
; create a new output KMT map:  opaleo
; ==============================================

  opaleo = kmt

; ==============================================
; Example of how to change points, e.g., convert land to ocean,
; deepen channels by increasing KMT, convert ocean to
; land, etc.
; ==============================================
  ; minloni = 150
  ; maxloni = 160
  ; minlatj = 242
  ; maxlatj = 250
  ; opaleo(minlatj:maxlatj,minloni:maxloni) = 5	; turn these points to KMT=5

  ; loni = 74
  ; latj = 74
  ; opaleo(latj,loni) = 0		; turning KMT in land

; ==============================================
;  write out new kmt and region_mask in direct access binary 
; ==============================================

  system("/bin/rm -f " + ofilebin1)
  fbindirwrite(ofilebin1,opaleo)

; ==============================================
;  write out new kmt and region_mask in netcdf to check 
; ==============================================

 system("/bin/rm -f " + ofile)
 fout = addfile(ofile,"c")

 fileAtt             = True
 fileAtt@title       = case
 fileAtt@create_date = systemfunc("date")
 fileAtt@src_code = "change_kmt.b2nc.ncl"
 fileAtt@author = systemfunc("echo $USER")
 fileAtt@input3 = ifile3
 fileAtt@input2 = ifile2
 fileAtt@input1 = ifile1
 fileattdef(fout,fileAtt)

 dimNames = (/"lat","lon"/)
 dimSizes = (/nlat,nlon/)
 dimUnlim = (/False,False/)

 otype = (/"lat","lon"/)

 filedimdef(fout,dimNames,dimSizes,dimUnlim)
 filevardef(fout,"kmt",vtype,otype)

 fout->Modifiedkmt = opaleo
 fout->inputKMT = kmt

 if (REMAP .eq. True) then
 	ofile2    = "kmt.remap."+dd+".nc"
 	system("/bin/rm -f " + ofile2)
 	fout2 = addfile(ofile2,"c")
 	fout2->OrigKmtRemap = rkmt
 end if
; ==============================================
; create plot
; ==============================================
  wks1  = gsn_open_wks("ps","kmt_paleoSphere."+case+"."+dd)     ; open a ps file
  gsn_define_colormap(wks1,"BlAqGrYeOrRe")     ; choose colormap

  plot1 = new(2,graphic)

  res                      = True     ; plot mods desired
  res2                     = True     ; plot mods desired
  res@gsnAddCyclic = True
  res2@gsnAddCyclic = True


  res@mpProjection         = "Orthographic"
  res@mpOutlineOn            = False
  res@mpFillOn               = False
; North America
  ; res@mpCenterLatF      =  30.
  ; res@mpCenterLonF      =  0.
; South America
  ; res@mpCenterLatF      =  -30.
  ; res@mpCenterLonF      =  0.
  res@mpMinLatF             =  31.         ; min lat
  res@mpMaxLatF             =  50.         ; max lat
  res@mpMinLonF             =  -78.        ; min lon
  res@mpMaxLonF             =  -70.        ; max lon

  res@mpGridLineColor = "transparent"
  res@mpGridAndLimbOn = True
  res@mpPerimOn       = False

  res2@mpGridLineColor = "transparent"
  res2@mpGridAndLimbOn = True
  res2@mpPerimOn              = False
  res2@mpOutlineOn            = False
  res2@mpFillOn               = False

  res@vpXF      = 0.05
  res@vpYF      = 0.9
  res@vpWidthF  = 0.9
  res@vpHeightF = 0.8

  res@gsnDraw  = False       ; don't draw the plots now
  res@gsnFrame = False       ; or advance the frame
  res2@gsnDraw              = False
  res2@gsnFrame             = False


  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillMode           = "RasterFill"         ; turn on raster mode
  res@cnLineLabelsOn       = False    ; turn off line labels
  res@cnLinesOn            = False    ; turn of contour lines
  res@cnMissingValPerimOn  = True     ; outline missing values
  res@tmLabelAutoStride    = True     ; nice stride on labels
  res@lbLabelAutoStride    = True     ; nice label bar labels
  res@gsnSpreadColors      = True     ; use full colormap

  res2@cnFillOn             = True     ; turn on color fill
  res2@cnLineLabelsOn       = False    ; turn off line labels
  res2@cnLinesOn            = False    ; turn of contour lines
  res2@cnMissingValPerimOn  = True     ; outline missing values
  res2@tmLabelAutoStride    = True     ; nice stride on labels
  res2@lbLabelAutoStride    = True     ; nice label bar labels
  res2@gsnSpreadColors      = True     ; use full colormap

  res@tiMainString         = case + " KMT"
  res@gsnRightString       = ""
  res2@gsnRightString      = ""

  plot1(0) = gsn_csm_contour_map(wks1,opaleo,res)    ; contour the variable
  plot1(1) = gsn_csm_contour_map(wks1,opaleo,res2)    ; contour the variable

  pres     		= True
  pres@gsnPanelRowSpec	= True
  pres@gsnMaximize	= True

  gsn_panel(wks1,plot1,(/1,2/),pres)

end
