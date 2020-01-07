;this procedure converts the rtm text file to netcdf
pro rtm,pnum,op=post,cp=closepost

close,/all
if !d.name eq 'X' then device, decomp=0

;-- User Modify infile

rtmfile1='fort.13_<casename>'
outfile='rdirc.1x1.<casename>.nc'
;resnum=1		; half-degree
resnum=2		; one-degree 

;--  original direction file files
nm=file_lines(rtmfile1)

openr, lun, rtmfile1, /get_lun
str=fltarr(3,nm)
readf,lun,str
free_lun,lun
rtm_lat=reform(str[0,*])
rtm_lon=reform(str[1,*])
rdir=reform(str[2,*])

del=abs(rtm_lon[0]-rtm_lon[1])

case resnum of
   1:begin ;half-degree
      im0=long(720)
      jm0=long(360)
   end
   2:begin ;one-degree
      im0=long(360)
      jm0=long(180)
   end 
   3:begin ;two-degree
      im0=long(180)
      jm0=long(90)
   end
endcase

;--  netcdf file  --------------------------------------------
id=ncdf_create(outfile,/clob)

 xid = NCDF_DIMDEF(id, 'ni', im0) ; Define the X dimension.  
 yid = NCDF_DIMDEF(id, 'nj', jm0) ; Define the Y dimension.  
 wid1 = NCDF_VARDEF(id, 'RTM_FLOW_DIRECTION', [xid,yid], /FLOAT)  
 lnid = NCDF_VARDEF(id, 'xc', [xid,yid], /FLOAT)  
 ltid = NCDF_VARDEF(id, 'yc', [xid,yid], /FLOAT)  

;--  Setup variable attributes  --------------------------
 NCDF_ATTPUT, id, /GLOBAL, 'title', $
   'River Transport Model (RTM) flow directions'  
 NCDF_ATTPUT, id, /GLOBAL, 'conventions', 'CF-1.0'
 NCDF_ATTPUT, id, /GLOBAL, 'SVD_ID', 'none'
 NCDF_ATTPUT, id, /GLOBAL, 'SVN_URL', 'none'
 NCDF_ATTPUT, id, /GLOBAL, 'history', 'created by <user> <date>'
 NCDF_ATTPUT, id, /GLOBAL, 'source', 'fort.13_<casename>'

 NCDF_ATTPUT, id, lnid, 'long_name', 'longitude of grid cell center'
 NCDF_ATTPUT, id, lnid, 'units', 'degrees_east'
 NCDF_ATTPUT, id, lnid, 'mode', 'time-invariant'
    
 NCDF_ATTPUT, id, ltid, 'long_name', 'latitude of grid cell center'
 NCDF_ATTPUT, id, ltid, 'units', 'degrees_north'
 NCDF_ATTPUT, id, ltid, 'mode', 'time-invariant'

 NCDF_ATTPUT, id, wid1, 'long_name', 'RTM flow direction'
 NCDF_ATTPUT, id, wid1, 'units', 'unitless'
 NCDF_ATTPUT, id, wid1, 'mode', 'time-invariant'
 NCDF_ATTPUT, id, wid1, 'comment', 'N,NE,E,SE,S,SW,W,NW = 1,2,3,4,5,6,7,8'

NCDF_CONTROL, id, /ENDEF        ; Put the file into data mode.  

n=long(0)
rdir2=fltarr(im0,jm0)
rslope2=fltarr(im0,jm0)
rlon2=fltarr(im0,jm0)
rlat2=fltarr(im0,jm0)
;--  inner loop over longitude
for j=0,jm0-1 do begin
    for i=0,im0-1 do begin
        rdir2[i,j]=rdir[n]
        rlon2[i,j]=rtm_lon[n]
        rlat2[i,j]=rtm_lat[n]
        n+=1
    endfor
endfor
ncdf_varput,id,wid1,rdir2
ncdf_varput,id,lnid,rlon2
ncdf_varput,id,ltid,rlat2
ncdf_close,id

!p.multi=0
!p.region=0
; print,crashnow
end
