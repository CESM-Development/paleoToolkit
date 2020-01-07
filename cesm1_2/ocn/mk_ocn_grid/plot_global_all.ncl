load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"

; modified from gridplot.ncl - author:  Mike Levy
; revised by:  nanr, dan lunt


begin

; *********************
;  refine_type = "q3"
; *********************

print(plot_type)
print(expt)
print(type)
print(date)
print(space)

out_type = "ps"
; out_type = "pdf"
;  out_type = "X11"


  gridfile = expt+"_"+date+".nc"
  ; gridfile = "gx1PT_"+date+".plot.nc"
  ; gridfile = "gx1PT_"+date+"_"+expt+"_"+type+".nc"
  if (.not.isfilepresent(gridfile)) then
    print((/"Can not find "+gridfile+"!"/))
    exit
  end if

  print((/"Plotting mesh from "+gridfile/))
  ; system("ln -s "+gridfile+" "+gridfile+".nc")

  f = addfile(gridfile+".nc","r")

  lon = f->grid_corner_lon*90.d/acos(0.d)
  lat = f->grid_corner_lat*90.d/acos(0.d)

  clon = f->grid_center_lon*90.d/acos(0.d)
  clat = f->grid_center_lat*90.d/acos(0.d)

  print("number of verticies = "+dimsizes(lon))
  print("lat min/max = "+min(lat)+" "+max(lat))

  nelem = dimsizes(lon(:,0))

  print("number of elements = "+nelem)


  ; Orthographic Projection
  gridname = "Grid_"+expt+"_"+date+"_"+type
  wks = gsn_open_wks(out_type,gridname+"_"+plot_type)

  res = True

  res@tiMainString = gridname
  res@mpOutlineOn            = False
  res@mpFillOn               = False

  res@mpProjection      = "Orthographic"
  
  res@mpLimitMode = "LatLon"

  if (plot_type.eq."sp") then
    ; South Pole
    res@mpCenterLatF      =  -90.
    res@mpCenterLonF      =  0.
    res@mpMinLonF         =  -180.
    res@mpMaxLonF         =  180.
    res@mpMaxLatF         =  -50.
    res@mpMinLatF         =  -90.
  end if

  if (plot_type.eq."np") then
    ; South Pole
    res@mpCenterLatF      =  90.
    res@mpCenterLonF      =  0.
    res@mpMinLonF         =  -180.
    res@mpMaxLonF         =  180.
    res@mpMaxLatF         =  90.
    res@mpMinLatF         =  50.
  end if

  if (plot_type.eq."q1") then
    ; 1st quadrant
    res@mpCenterLatF      =  0.
    res@mpCenterLonF      =  0.
  end if

  if (plot_type.eq."q2") then
    ; 2nd quadrant
    res@mpCenterLatF      =  0.
    res@mpCenterLonF      =  90.
  end if

  if (plot_type.eq."q3") then
    ; 3rd quadrant
    res@mpCenterLatF      =  0.
    res@mpCenterLonF      =  180.
  end if

  if (plot_type.eq."q4") then
    ; 4th quadrant
    res@mpCenterLatF      =  0.
    res@mpCenterLonF      =  270.
  end if

  res@vpXF      = 0.05
  res@vpYF      = 0.9
  res@vpWidthF  = 0.9
  res@vpHeightF = 0.8

  res@gsnDraw  = False       ; don't draw the plots now
  res@gsnFrame = False       ; or advance the frame

  res@cnFillOn = True
  res@cnFillMode = "RasterFill"
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels   = (/0.5/)
  res@cnFillColors = (/"lightblue","grey"/)
  res@cnLineLabelsOn = False			; turn off contour labels
  res@mpGridLineColor = "transparent"
  res@mpGridAndLimbOn = True
  res@mpPerimOn              = False

  ; res@mpInlandWaterFillColor = "Blue"
  ; res@mpGreatCircleLinesOn = True

  gsn_define_colormap(wks,(/"White","Black","Tan","LightBlue","Blue"/))

  res@sfXArray = f->grid_center_lon*90./acos(0.)
  res@sfYArray = f->grid_center_lat*90./acos(0.)
  ; ikmt = (/f->grid_ikmt/)
  ikmt = (/f->grid_imask/)
  ikmt = where(ikmt .eq. 0, 1, 0)

  plot = gsn_csm_contour_map(wks,ikmt, res)
  ; plot = gsn_csm_contour_map(wks,(/f->grid_ikmt/), res)
  ; plot = gsn_csm_contour_map(wks,kmt, res)
  draw(plot)

  ; polygon resources
  res_p             = True

  res_p@gsLineThicknessF = 0.2
  res_p@gsLineColor   = "black"
  res@mpPerimOn              = True

  j1 = ispan(0,3,1)
  j2 = (j1+1)%4

  xlon = new ( (/5/), "double")
  xlat = new ( (/5/), "double")

  xlon = new( (/5/), double)
  xlat = new( (/5/), double)
 
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.005                     ; font smaller. default big

  space2 = 5000
  space = 20
  do i=0,nelem-1,space
    if ( mod(i,space2).eq.0) then
      print ("i = "+i+"/"+(nelem-1))
    end if
      xlon(0:3) = lon(i,0:3)
      xlat(0:3) = lat(i,0:3)

    do j=0,2
      if ( abs(xlon(j+1)-xlon(j)) .gt. 180.0) then
        if (xlon(j+1) .gt. xlon(j) ) then
          xlon(j)=xlon(j)+360.
        else
          xlon(j+1)=xlon(j+1)+360.
        end if
      end if
    end do

    xlon(4)=xlon(0)
    xlat(4)=xlat(0)

    gsn_polyline(wks, plot, xlon,xlat,res_p)
    res_p@gsFillColor = ikmt(i)*10
    gsn_polygon(wks, plot, xlon,xlat,res_p)
;    txt = gsn_add_text(wks,plot,sprintf("%5.2f",clon(i)),clon(i),clat(i),txres) 
  end do

;  draw(plot)

delete(wks)





end


