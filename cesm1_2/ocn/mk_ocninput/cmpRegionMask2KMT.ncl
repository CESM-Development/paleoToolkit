load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
 print ("=========================================")
 print ("Starting:  chk_grid.ncl")
 print ("Start Time: "+systemfunc("date") )
 print ("=========================================")

 dd= systemfunc("date -u +%y%m%d")

 IPATH  = "./"
 ifile1 = IPATH + "/kmt.<iter>.da"           	      ; kmt in bin format
 ifile2 =  "region.<case>.be.ieeei4"          ; region_mask in bin format
 ifile11 = "/glade/p/cesmdata/cseg/inputdata/ocn/pop/gx1v6/grid/topography_20090204.ieeei4"
 ifile22 = "/glade/p/cesmdata/cseg/inputdata/ocn/pop/gx1v6/grid/region_mask_20090205.ieeei4"
 ofile  = "./region_mask-kmt.nc"

;NOTE lat and lon are hardcoded to regular POP Grid
 nlon = 320
 nlat = 384 
 dims = (/nlat,nlon/)

 vtype = "float"
; ==============================================
; == Read in kmt from topography file
; == Read binary direct access file and rewrite to netcdf
; ----------------------------------------------------
; ==============================================
; NOTE:  fortran is column major .... NCL [C] is row major
; fortran(nlon,nlat)  <===> NCL(nlat,nlon)
; ==============================================

  bigEnd = True
  ; if (bigEnd .eq. True) then
   ; if (.not. isbigendian()) then
          setfileoption ("bin", "ReadByteOrder", "BigEndian")
          setfileoption ("bin", "WriteByteOrder","BigEndian")
   ; end if
   ikmt1be = fbindirread(ifile1,0,dims,vtype)
   ireg2be = fbindirread(ifile2,0,dims,vtype)
   ikmt11be = fbindirread(ifile11,0,dims,vtype)
   ireg22be = fbindirread(ifile22,0,dims,vtype)
  ; else
          setfileoption ("bin", "ReadByteOrder", "LittleEndian")
          setfileoption ("bin", "WriteByteOrder","LittleEndian")
   ikmt1le = fbindirread(ifile1,0,dims,vtype)
   ireg2le = fbindirread(ifile2,0,dims,vtype)
   ikmt11le = fbindirread(ifile11,0,dims,vtype)
   ireg22le = fbindirread(ifile22,0,dims,vtype)
  ; end if
 

; ==============================================
; Mask kmt and region_masks to compare records.
; ==============================================

print("RegionMask + bigEnd (paleo first; default second)" )
printMinMax(ireg2be,0) 
printMinMax(ireg22be,0)
print("RegionMask + litEnd (paleo first; default second)" )
printMinMax(ireg2le,0) 
printMinMax(ireg22le,0)
print("kmt + bigEnd (paleo first; default second)" )
 printMinMax(ikmt1be,0) 
printMinMax(ikmt11be,0)
print("kmt + ligEnd (paleo first; default second)" )
 printMinMax(ikmt1le,0) 
 printMinMax(ikmt11le,0)
 dff = ikmt1be;
 ireg = ireg2be
 ikmt = ikmt1be
 ; ikmt = where(ikmt.eq.0,0,1)
 ; ireg = where(ireg.eq.1,1,0)
 ; ikmt = where(ikmt.eq.0,0,1)
 ; ireg = where(ireg.eq.0,0,1)
 dff = ireg - ikmt;
 ; dff = ireg 
 if((any(dff.ne.0)) )

 	nlon  = 320
 	nlat  = 384 

	lat = ispan(1,384,1)
	lon = ispan(1,320,1)

	print("Found some errors.  Printing output file " + ofile)
  	system("/bin/rm -f " + ofile)
 	fout = addfile(ofile,"c")
 	fileAtt             = True
 	fileAtt@title       = ireg + " - " + ikmt
 	fileAtt@conventions = "CCSM"
 	fileAtt@create_date = systemfunc("date")

 	dimNames = (/"lat","lon"/)
 	dimSizes = (/nlat,nlon/)
 	dimUnlim = (/False,False/)

 	otype = (/"lat","lon"/)

 	filedimdef(fout,dimNames,dimSizes,dimUnlim)
 	filevardef(fout,"kmt",vtype,otype)
 	filevardef(fout,"reg",vtype,otype)
	
 	fout->dff = dff
 	fout->kmt = ikmt
 	fout->reg = ireg

	; plot results

 	; wks = gsn_open_wks("X11","diffMap")
 	wks = gsn_open_wks("ps","diffMap")
 	gsn_define_colormap(wks,"BlueDarkRed18")

 	res = True
 	res@gsnMaximize = True
 	res@gsnPaperOrientation = "portrait"
 	res@cnFillOn = True
	res@mpOutlineOn            = False         ; Turn on map outlines
	res@mpFillOn               = False         ; Turn on map fill
 	; res@lbLabelBarOn = False
 	res@gsnCenterString = "("+ ireg + ") - (" + ikmt + ")"
 	res@gsnLeftString =  " "
 	res@gsnRightString = " "
 	plot = gsn_csm_contour_map_ce(wks,dff,res)
 else
	print("Region mask and KMT match exactly")
 end if

       ; plot results

        ; wks = gsn_open_wks("X11","diffMap2")
        wks = gsn_open_wks("ps","diffMap2")
        gsn_define_colormap(wks,"BlueDarkRed18")

        res = True
        res@gsnMaximize = True
        res@gsnPaperOrientation = "portrait"
        res@cnFillOn = True
        res@mpOutlineOn            = False         ; Turn on map outlines
        res@mpFillOn               = False         ; Turn on map fill
        res@lbLabelBarOn = False
        res@gsnCenterString = "("+ ikmt + ") - (" + ireg + ")"
        res@gsnLeftString =  " "
        res@gsnRightString = " "
        ; plot = gsn_csm_contour_map_ce(wks,dff,res)

end
