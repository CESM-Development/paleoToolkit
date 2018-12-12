;
; Regrid Glacier and PFT datasets to half degree
;
begin
   glacierfile = "$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_glacier_3x3min_simyr2000.c120926.nc"
   pftfile     = "$DIN_LOC_ROOT/lnd/clm2/rawdata/pftcftdynharv.0.25x0.25.LUH2.simyrPtVg.c181106/mksrf_landuse_potvegclm50_LUH2.c181106.nc"
   gridfile    = "$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_lanwat.050425.nc"


   print( "Open: "+glacierfile );
   ncglc = addfile( glacierfile, "r" );
   print( "Open: "+pftfile );
   ncpft = addfile( pftfile, "r" );
   print( "Open: "+gridfile );
   ncgrd = addfile( gridfile, "r" );

   sdate    = systemfunc( "date +%y%m%d" );

   outglacierfile = "mksrf_glacier_0.5x0.5_simyr2000.c"+sdate+".nc";
   outpftfile     = "mksrf_landuse_0.5x0.5_potvegclm50_LUH2.c"+sdate+".nc";

end
