;
; Regrid Glacier and PFT datasets to half degree
;
procedure duplicate_file4newgrid ( outfilename:string, ncin:file, ncout:file, outlat, outlon, outmsk, outlatixy, outlongxy )
;
; Duplicate the input file on the output file with a new grid
;
begin
   print( "Work on "+outfilename );
   ;
   ; Define Dimensions
   ;
   print( "Define dims" );
   nlat = dimsizes(outlat);
   nlon = dimsizes(outlon)
   dimnames = getfiledims( ncin );
   dsizes   = getfiledimsizes( ncin );
   is_unlim = new( (/ dimsizes(dimnames) /), "logical" );
   is_unlim = False;
   do i = 0, dimsizes(dimnames)-1
      if ( dimnames(i) .eq. "time" )then
         is_unlim(i) = True;
      else
         if ( dimnames(i) .eq. "lat" )then
            dsizes(i) = nlat;
         else
            if ( dimnames(i) .eq. "lon" )then
               dsizes(i) = nlon;
            end if
         end if
      end if
   end if
   filedimdef( ncoout, dimnames, dsizes, is_unlim );
end

begin
   ;
   ; Set Input files and open them
   ;
   glacierfile = "$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_glacier_3x3min_simyr2000.c120926.nc"
   pftfile     = "$DIN_LOC_ROOT/lnd/clm2/rawdata/pftcftdynharv.0.25x0.25.LUH2.simyrPtVg.c181106/mksrf_landuse_potvegclm50_LUH2.c181106.nc"
   gridfile    = "$DIN_LOC_ROOT/lnd/clm2/rawdata/mksrf_lanwat.050425.nc"

   print( "Open: "+glacierfile );
   ncglc = addfile( glacierfile, "r" );
   print( "Open: "+pftfile );
   ncpft = addfile( pftfile, "r" );
   print( "Open: "+gridfile );
   ncgrd = addfile( gridfile, "r" );
   ;
   ; Open output files for half degree versions
   ;
   sdate    = systemfunc( "date +%y%m%d" );

   outglacierfile = "mksrf_glacier_0.5x0.5_simyr2000.c"+sdate+".nc";
   outpftfile     = "mksrf_landuse_0.5x0.5_potvegclm50_LUH2.c"+sdate+".nc";

   print( "Open : "+outglacierfile );
   if ( systemfunc("test -f "+outglacierfile+"; echo $?" ) .eq. 0 )then
      system( "/bin/rm "+outglacierfile );
   end if
   nco_glc = addfile( outglacierfile, "c" );
   print( "Open : "+outpftfile );
   if ( systemfunc("test -f "+outpftfile+"; echo $?" ) .eq. 0 )then
      system( "/bin/rm "+outpftfile );
   end if
   nco_pft = addfile( outpftfile, "c" );
   ;
   ; Get grid information on the half degree to regrid to
   ;
   outlat = ncgrd->LAT;
   outlon = ncgrd->LON;
   outmsk = ncgrd->LANDMASK;
   outlatixy = ncgrd->LATIXY;
   outlongxy = ncgrd->LONGXY;

   duplicate_file4newgrid ( outglacierfile, ncglc, nco_glc, outlat, outlon, outmsk, outlatixy, outlongxy )
   duplicate_file4newgrid ( outpftfile,     ncpft, nco_pft, outlat, outlon, outmsk, outlatixy, outlongxy )

   print( "" );
   print( "Successfully regrid files" );
end
