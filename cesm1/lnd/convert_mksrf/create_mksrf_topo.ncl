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

  timePeriod = "15ka"
  ifile1 = "topo_15ka_10min.140226.nc"
  ifile2 = "/glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/mksrf_vocef_0.5x0.5_simyr2000.c110531.nc"	; template file
  in2 = addfile(ifile2,"r")
  LAT  = in2->lat
  LON  = in2->lon
  LATIXY  = in2->LATIXY
  LONGXY  = in2->LONGXY

  ofile = "mksrf_topo_15ka.0.5x0.5."  + dd + ".nc"
  system("/bin/rm -f " + ofile)
  fout = addfile(ofile,"c")

  in1 = addfile(ifile1,"r")
  htopo180 = lonFlip(in1->htopo)
  ilat  = htopo180&lat;
  ilon  = htopo180&lon;

  htopo05 = linint2_Wrap(ilon,ilat,htopo180,False,LON,LAT,0)

  fout->TOPO = htopo05
  fout->TOPO_ICE = htopo05
  fout->LONGXY = LONGXY
  fout->LATIXY = LATIXY

  fout@source = "create_mksrf_topo.ncl"
  fout@date = dd
  fout@input1 = ifile1

end

