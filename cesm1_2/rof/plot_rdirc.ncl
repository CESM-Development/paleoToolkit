;======================================================================;
; reads in fort.## ascii files
; compares original and corrected files and plots
; 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;====================user defined parameters

filename1    = getenv("RFILE1")
filename_e1  = getenv("EFILE1")
filename2    = getenv("RFILE2")
filename_e2  = getenv("EFILE2")
case         = getenv("CASE")
plotname     = getenv("PLOTNAME")

print("processing case = " + case)

ilat         = getenv("NLAT")
ilon         = getenv("NLON")

flat         = stringtofloat(ilat)
flon         = stringtofloat(ilon)
nlat         = floattoint(flat)
nlon         = floattoint(flon)


;============== read atm ascii file
;2x2 resolution:
; nlat = 90
; nlon = 180

;1x1 resolution:
; nlat = 180
; nlon = 360

;=============== read in lat, lon, idir from rdirc file
print("Rdirc data being read from file:" + filename1 + " and " + filename2)

 lat   = new((/nlat/),float)
 lon   = new((/nlon/),float)
 lat = latGlobeFo(nlat,"lat","latitude","degrees_north")
 lon = lonGlobeFo(nlon,"lon","longitude","degrees_east")-180. 

lons=new((/nlat,nlon/),float)
lats=new((/nlat,nlon/),float)




 UR1    = new((/nlat,nlon/),float)
 VR1    = new((/nlat,nlon/),float)
 MI1    = new((/nlat,nlon/),float)
 UR2    = new((/nlat,nlon/),float)
 VR2    = new((/nlat,nlon/),float)
 MI2    = new((/nlat,nlon/),float)
idirc = new((/nlat,nlon/),integer)

lat!0= "lat"
lat&lat = lat 
lon!0= "lon"
lon&lon = lon 

do i=0,nlon-1
lats(:,i)=lat(:)
end do

do j=0,nlat-1
lons(j,0:nlon/2-1)=lon(0:nlon/2-1)+360
lons(j,nlon/2:nlon-1)=lon(nlon/2:nlon-1)
end do


printVarSummary(lon)

;idirc:    0        1    2        3     4          5     6         7     8 
 U0    = (/0.,      0.,  1., sqrt(2.),  1.,        0.,  -1., -sqrt(2.), -1./)
 V0    = (/0., sqrt(2.), 1.,       0., -1., -sqrt(2.),  -1.,       0.,   1./)
;normalize:
 U    = U0/sqrt(2.) 
 V    = V0/sqrt(2.)

ksoi1=asciiread(filename1,(/nlon*nlat,3/),"float")
ksoi2=asciiread(filename2,(/nlon*nlat,3/),"float")

f11=asciiread(filename_e1,-1,"float")
f12=asciiread(filename_e2,-1,"float")

if (dimsizes(f11).ne.1) then
  print("number of infinite loops: "+dimsizes(f11)/4.0)
  xx1=dimsizes(f11)
else
  print("no infinite loops!")
  xx1=0
end if
if (dimsizes(f12).ne.1) then
  print("number of infinite loops: "+dimsizes(f12)/4.0)
  xx2=dimsizes(f12)
else
  print("no infinite loops!")
  xx2=0
end if

MI1(:,:)=0
MI2(:,:)=0
do zz1=0,(xx1/4.0)-1
MI1=MI1+where(lats.eq.(f11(zz1*4+2)+0.25).and.lons.eq.(f11(zz1*4+3)+0.25),1,0)
end do
do zz2=0,(xx2/4.0)-1
MI2=MI2+where(lats.eq.(f12(zz2*4+2)+0.25).and.lons.eq.(f12(zz2*4+3)+0.25),1,0)
end do

print("first max/min/total:")
print(min(MI1))
print(max(MI1))
print(sum(MI1))

print("second max/min/total:")
print(min(MI2))
print(max(MI2))
print(sum(MI2))


if (max(MI1).ne.1 .and. min(MI1).ne.0) then
print("max value of first MI is not 1 and min value of first MI is not 0")
end if
if (max(MI2).ne.1 .and. min(MI2).ne.0) then
print("max value of second MI is not 1 and min value of second MI is not 0")
end if

print(" rdirc data read")

  id0      = floattointeger(ksoi1(:,2))
  id1      = onedtond(id0,(/nlat,nlon/))

  id02     = floattointeger(ksoi2(:,2))
  id2      = onedtond(id02,(/nlat,nlon/))


do i=0, nlon-1
do j=0, nlat-1

  coslat = cos(lat(j)*2*3.14159265/360.)
  UA1 = (/U0(id1(j,i))/)*coslat
  VA1 = (/V0(id1(j,i))/)

  UA2 = (/U0(id2(j,i))/)*coslat
  VA2 = (/V0(id2(j,i))/)
  
  AMP1 = sqrt(UA1*UA1 + VA1*VA1)/(MI1(j,i)+1.1)
  AMP2 = sqrt(UA2*UA2 + VA2*VA2)/(MI2(j,i)+1.1)

  if( AMP1 .ne. 0.) then
   UR1(j,i) = UA1/AMP1
   VR1(j,i) = VA1/AMP1
  else
   UR1(j,i) = 0.
   VR1(j,i) = 0.
  end if
  if( AMP2 .ne. 0.) then
   UR2(j,i) = UA2/AMP2
   VR2(j,i) = VA2/AMP2
  else
   UR2(j,i) = 0.
   VR2(j,i) = 0.
  end if

end do
end do

lon@units = "degrees_east"
UR1!0 = "lat"
UR1!1 = "lon"
UR1&lat = lat
UR1&lon = lon  

VR1!0 = "lat"
VR1!1 = "lon"
VR1&lat = lat
VR1&lon = lon  


UR2!0 = "lat"
UR2!1 = "lon"
UR2&lat = lat
UR2&lon = lon  

VR2!0 = "lat"
VR2!1 = "lon"
VR2&lat = lat
VR2&lon = lon  

;================= plot
res = True
PSName= "ps"
wks = gsn_open_wks(PSName,plotname)        ; open a ncgm file

res@vcRefMagnitudeF = 3 
res@vcPositionMode  = "ArrowTail" 
res@vcRefLengthF    = 0.035
res@vcLineArrowThicknessF = 0.4 
res@vcLineArrowHeadMinSizeF = 0.001
res@vcLevelSpacingF = 0.5
res@vcMinLevelValF = 0
res@vcMaxLevelValF = 2
res@gsnDraw         = False               ; Do not draw plot
res@gsnFrame        = False               ; Do not advance frome
res@mpFillOn        = False
res@mpOutlineOn        = False
res@vcMonoLineArrowColor = False
res@vcLevelSelectionMode = "ManualLevels" 

; Zoom in on data
;res@mpMinLatF = -80
;res@mpMaxLatF = -60
;res@mpMinLonF = 335
;res@mpMaxLonF = 355

gsn_define_colormap(wks,"temp1")

resP = True
resP@gsnMaximize = True
res@gsnCenterString = case+" Runoff Vectors with interior basins"

plt= new(2,graphic)
plt(0) = gsn_csm_vector_map(wks,UR1,VR1,res)
res@gsnCenterString = case+" redirected Runoff Vectors"
plt(1) = gsn_csm_vector_map(wks,UR2,VR2,res)
;gsn_panel(wks,plt(1),(/1,1/),resP)
gsn_panel(wks,plt,(/2,1/),resP)
print("Data plotted")

end
    
