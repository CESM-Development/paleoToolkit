load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; ----------------------------------------------------
; Read binary direct access file and rewrite to netcdf
; ----------------------------------------------------

begin
 print ("=========================================")
 print ("Starting:  mk_TSic.ncl")
 print ("Start Time: "+systemfunc("date") )
 print ("=========================================")

 setfileoption ("bin", "ReadByteOrder", "BigEndian")
 setfileoption ("bin", "WriteByteOrder","BigEndian")

 PROC  = "/myDir/"
 HIST  = "/myHistPath/"

;  in this example I am using the monthly ocean history values from year December of year 451 to initialize.
 ocnTSfile = "myCase.pop.h.0451-12.nc"
 fkmti  = PROC  + "/topography_20090204.nc"			; kmt on old tgrid 
 fkmtn  = PROC  + "/topography_plioWAIS_gx1v6_141022.nc"	; kmt on new tgrid 
 fTSi   = HIST  + "/"+ocnTSfile		; TS history file 
 casename  = "b40.plio.FV1.003"

 dd= systemfunc("date -u +%y%m%d")
 ofile    = "TSic_myCase.pop.h.0451-12_ITF_"+dd+".nc"
 ofilebin = "TSic_myCase.pop.h.0451-12_ITF_"+dd+".ieeer8"
 vtype = "double"
 slat = 50		; need to fill in all of lower hemisphere 

km = 60
;new grid dimensions:
nlat = 384
nlng = 320
;old grid dimensions:
nlato = 384
nlngo = 320

 dims3n = (/km,nlat,nlng/)
 dims2n = (/nlat,nlng/)

; ==============================================
; == Read in kmt from topography files
; ==============================================
 ifile = addfile(fkmti,"r")		; old KMT
 ikmt = ifile->kmt                
 delete(ifile)

 ifile = addfile(fkmtn,"r")		; new KMT
 nkmt = ifile->kmt_wais                
 delete(ifile)
; ==============================================
; == Read in T and S from old history file
; ==============================================
 ifile = addfile(fTSi,"r")
 TEMPi = ifile->TEMP(0,:,:,:)
 SALTi = ifile->SALT(0,:,:,:)
 if(isatt(SALTi,"scale_factor")) then
   SALTi = SALTi*SALTi@scale_factor  
   delete(SALTi@scale_factor)
 end if
 delete(ifile)

; ==========================================================
; == Replace new T and S with T and S from old history file
; ==========================================================
  TEMPn = new(dims3n,vtype)
  SALTn = new(dims3n,vtype)

  TEMPn=TEMPi(:,:,:)
  SALTn=SALTi(:,:,:)

; CHANGE YOUR GRID POINTS HERE. using this logic
; old ocean point =  new ocean point--> gets same T and S:
;    TEMPn(:,j,i)=TEMPo(:,j,i)
;    SALTn(:,j,i)=SALTo(:,j,i)
;new ocean point, old land point, gets zonal ave T and S:
;     TEMPn(k,j,i)=(/zTEMPo(k,j)/)
;     SALTn(k,j,i)=(/zSALTo(k,j)/)


  system("/bin/rm -f " + ofilebin)

; Note:  NCL switches dimension order:  (/lev,lat,lon/)  => binary(lon,lat,lev) which is how the ocn files are stored.

   fbindirwrite(ofilebin,TEMPn)
   fbindirwrite(ofilebin,SALTn)

; ==============================================
;  write out new  T and S in netcdf to check it
;  NOTE: TEMP-> TEMPERATURE; SALT->SALINITY for netcdf input to pop2
; ==============================================

 system("/bin/rm -f " + ofile)
 fout = addfile(ofile,"c")

 fileAtt             = True
 fileAtt@title       = "TEMP and SALT Initial Condition file from "+casename
 fileAtt@conventions = "none"
 fileAtt@create_date = systemfunc("date")
 fileAtt@src_code = "mk_TSic.ncl"
 fileAtt@author = "Esther Brady and Nan Rosenbloom"
 fileAtt@input = ocnTSfile 
 fileattdef(fout,fileAtt)

 dimNames = (/"depth","lat","lon"/)
 dimSizes = (/km,nlat,nlng/)
 dimUnlim = (/False,False,False/)

 otype = (/"depth","lat","lon"/)

 filedimdef(fout,dimNames,dimSizes,dimUnlim)
 filevardef(fout,"TEMPERATURE","double",otype)
 filevardef(fout,"SALINITY","double",otype)

 fout->TEMPERATURE = (/TEMPn/) 
 fout->SALINITY = (/SALTn/) 

end
