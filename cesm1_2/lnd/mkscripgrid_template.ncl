;
;  mkscripgrid.ncl
;
;  Create SCRIP grid and mapping file for a land-only point or region.
;  Requires NCL 6.1.0 or later for the ESMF regridding functions
;
;  Erik Kluzek
;  Dec/07/2011
;  $Id: mkscripgrid.ncl 47629 2013-05-31 08:59:50Z erik $
;  $HeadURL;
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
begin
  ; ===========================================================================================================
  ; Set a few constants needed later
  cdate  = systemfunc( "date +%y%m%d" );
  ldate  = systemfunc( "date" );
  ;
  ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE or use ENV VARIABLE SETTINGS
  ; Edit the following as needed to interpolate to a new resolution.
  ;
  ; Input resolution and position
  ;

  name  = "<casename>"                    ; name time period
  ipath = "./"
  fn1   = "<lanwat_file>"
  f1    = addfile(ipath + fn1,"r")
  ilat   = f1->lat
  ilon   = f1->lon
  latS  = f1->EDGES
  latN  = f1->EDGEN
  lonE  = f1->EDGEE 
  lonW  = f1->EDGEW 

  imask= f1->LANDMASK
  printVarSummary(imask)
  printVarSummary(lonW)

  nx    = dimsizes(ilon)
  ny    = dimsizes(ilat)

  outfilename = "SCRIPgrid_"+name+"_"+ny+"x"+nx+"."+cdate+".nc"

     printn = False;            ; set to True for debugging
     printn = True;             ; set to True for debugging


  system( "/bin/rm -f "+outfilename );
  if ( printn )then
     print( "output file: "+outfilename );
  end if

  ;
  ; Compute derived quantities
  ;

  ; ns = nx * ny;

  delX       = (lonE - lonW) / int2dble(nx);
  delY       = (latN - latS) / int2dble(ny);

  lonCenters = fspan( (lonW + delX/2.d0), (lonE - delX/2.d0), nx)
  latCenters = fspan( (latS + delY/2.d0), (latN - delY/2.d0), ny)
  lon   = new( (/ny, nx/), "double" );
  lat   = new( (/ny, nx/), "double" );
  lonCorners = new( (/ny, nx, 4/), "double" );
  latCorners = new( (/ny, nx, 4/), "double" );
  do i = 0, nx-1
    lat(:,i)          =  latCenters;
    latCorners(:,i,0) =  latCenters - delY/2.d0;
    latCorners(:,i,1) =  latCenters + delY/2.d0;
    latCorners(:,i,2) =  latCenters + delY/2.d0;
    latCorners(:,i,3) =  latCenters - delY/2.d0;
  end do
  do j = 0, ny-1
    lon(j,:)          =  lonCenters;
    lonCorners(j,:,0) =  lonCenters - delX/2.d0;
    lonCorners(j,:,1) =  lonCenters - delX/2.d0;
    lonCorners(j,:,2) =  lonCenters + delX/2.d0;
    lonCorners(j,:,3) =  lonCenters + delX/2.d0;
  end do

  ; for some reason, "No_FillValue" isn't working in the case where imask=1
  Mask2D = new( (/ny,nx/), "integer", "No_FillValue" )
  Mask2D(:,:) = floattointeger(imask)
  gridSize = delX+"x"+delY

  ;
  ; Create SCRIP grid file
  ;
  
  Opt = True
  Opt@Mask2D   = Mask2D
  Opt@GridCornerLat = latCorners
  Opt@GridCornerLon = lonCorners
  Opt@Title = "SCRIP grid file for "+name
  if (printn) then
    Opt@Debug = True
  end if
  curvilinear_to_SCRIP(outfilename, lat, lon, Opt)

  ;
  ; Add global attributes to file
  ;

  nc = addfile( outfilename, "w" );
  nc@history = ldate+": create using mkscripgrid_paleo.ncl";
  nc@comment = "Ocean is assumed to non-existant at this point";
  nc@Version = "$HeadURL: https://svn-ccsm-models.cgd.ucar.edu/clm2/branch_tags/cesm1_2_rel_tags/cesm1_2_0_n02_clm4_5_07/models/lnd/clm/tools/shared/mkmapgrids/mkscripgrid.ncl $";
  nc@Revision = "$Id: mkscripgrid.ncl 47629 2013-05-31 08:59:50Z erik $";
  if ( printn )then
    print( "================================================================================================" );
    print( "Successfully created SCRIP grid file: "+outfilename);
  end if

  ; ===========================================================================================================

end
