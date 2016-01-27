;*************************************************
; 
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

  timeperiod = "15ka"
  ifile1 = "/glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/mksrf_urban_3den_0.5x0.5_simyr2000.c090223_v1.nc"
  ifile2 = "mksrf_pft_15ka.140501.noCrop.nn.nc"


  ofile = "mksrf_urban_"+timeperiod+"." + dd + ".nc"
  system("/bin/rm -f " + ofile)
  fout = addfile(ofile,"c")

  in1 = addfile(ifile1,"r")
  in2 = addfile(ifile2,"r")

  landmask = in2->LANDMASK

  ; copy variables from standard file
  in1 = addfile(ifile1,"r")
  in1_global= getvaratts(in1)
  in1_vars = getfilevarnames(in1)

  ; copy global attributes
  print(" creating mksrf_urban_lgm21ka")
  if (.not.all(ismissing(in1_global))) then
  	do i = 0, dimsizes(in1_global) - 1
  		print("copy_fileatts: global attributes->" + in1_global )
  		fout@$in1_global$ = in1@$in1_global$
  	end do
  end if

  ; copy variables
  if (.not.all(ismissing(in1_vars))) then
  	do i = 0, dimsizes(in1_vars) - 1
	    print("i = " + i + " in1vars = " + in1_vars(i))
	    if (in1_vars(i) .eq."PCT_URBAN") then
  			pct_urban = in1->$in1_vars(i)$ 
			pct_urban = 0
  		        fout->$in1_vars(i)$ = pct_urban
	    else
	    	if (in1_vars(i) .eq."LANDMASK") then
  			fout->$in1_vars(i)$ = landmask
		else
  			fout->$in1_vars(i)$ = in1->$in1_vars(i)$ 
                end if
            end if
  	end do
  end if

end

