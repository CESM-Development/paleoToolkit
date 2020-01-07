program paleo_mkraw_ccsm4

  implicit none
  include 'netcdf.inc'

!----------------------------------------------------------------- 
!
! make raw datasets for Paleo CCSM simulation
! input data: (2x2 nc) LSM vegetation types from desired period
!	    : make sure lat/lon in netcdf is compatible with 2dlatlon
!	      definition below (see code)  
!           : (mksrf_soitex.10level.nc) IGBP soil texture for 0k
!	    : zonally averaged model organic soil density distribution
!	    : 2x2 topography/bathymetry file (topo_depth) file 
!
! output data: raw data files necessary to create surface-data file
!              at model runtime: mksrf_glacier.nc
!			       : mksrf_urban.nc
!			       : mksrf_lanwat.nc
!			       : mksrf_soicol_clm2.nc 
!			       : mksrf_lai.nc
!			       : mksrf_soitex.10level_paleo.nc
! updates for ccsm4
!                              : mksrf_organic_paleo.nc
!                              : mksrf_fmax_paleo.nc
!                              : mksrf_landuse_paleo.nc (replaces mksrf_pft_paleo.nc)
!                              : mksrf_topo_paleo.nc
!                              : mksrf_vocef_paleo.nc
!
!            : output file names specifed below
!
! from 2x2 LSM veg type data, this program will create pfts, glacier,
! urban, lanwat, soilcol, lai  fields. the soil texture fields will
! place soil profiles on
! the paleo grid. color will be assigned to a value of 15. (similar to ccsm3 soilcolor=4).
! glacier,urban,and lakes will all be set to zero. 
!
! ccsm4 updates:
! --------------
! ccsm4 updates include an updated pft file which includes pft harvest and grazing
! variables, organic soil density, max saturated areas for hydrology, topography
! placeholders for eventual land-ice model, and tree emmissions factors for atmos
! chemistry.
!
! currently we set harvest, grazing, and tree emissions to zero. 
! if using waccm (atmo chem), and want to specify emissions factors,
! you will need to develop an equivalence algorithm to modern values (see default
! mksrf_vocef file in inputdata) and create your own code.  also assuming
! assuming harvest and grazing to be zero for deep time. if not, develop personal
! code.  placeholders exist in below subroutines.
!
! organic soil density is assumed to approximate zonally averaged modern values
! max saturation areas (proxy for wetlands and used in the hydrology) is
! currently given a global number (like soil color) obtained by global avg
! modern value. 
!
! this program is designed to provide a framework to create the
! raw data files. any modifications to the aligorthims in this code
! can be applied by the user appropriate for the time period. 
!
! also note that ccsm4 has 17 (0 to 16) pfts as opposed to 16 (0 to 15) 
! for ccsm3. however, the new pft is crop2 and will NOT be used for 
! paleo. therefore, no modifcations to this code are necessary.  
!
! after paleo_mkraw datasets have been created, run mksurfdata clm tool to
! create surface_dataset offline.
! 
! cshields, slevis,  feb 2003, may 2010
!-----------------------------------------------------------------


  integer, parameter :: r8 = selected_real_kind(12)

  integer, parameter :: nlon = 180        !input grid : longitude points
  integer, parameter :: nlat =  90        !input grid : latitude  points
  integer, parameter :: numpft = 16       !number of plant types
  integer, parameter :: time   = 12       !number of months in year 
  integer, parameter :: nlay = 10         !input grid : number of soil layers
 
  integer, parameter :: nlon_st = 4320      !soil input grid : longitude points
  integer, parameter :: nlat_st = 2160      !soil input grid : latitude  points
  integer, parameter :: nmapunits = 4931 !soil input grid : # of igbp soil 'mapunits'
  integer, parameter :: mapunitmax = 6998!soil input grid : max value of 'mapunits'

  integer, parameter :: nlat_org =  360      !organics input grid : latitude points

  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)  
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid
  real(r8) :: dx,dy                       !grid increments


  real(r8) :: landmask(nlon,nlat)         !fraction of land

  real(r8) :: dzsoi(10), zsoi(10)        !soil layer thickness and depth
  real(r8) :: pct_sand0(nlay,mapunitmax)    !original percent sand 
  real(r8) :: pct_clay0(nlay,mapunitmax)    !original percent clay 
  real(r8) :: mapunits0(nlon_st,nlat_st)     !mapunits 

  real(r8) :: pct_pft(nlon,nlat,0:numpft) !percent pft 
  real(r8) :: pct_glacier(nlon,nlat)      !percent glacier
  real(r8) :: pct_urban(nlon,nlat)        !percent urban 
  real(r8) :: pct_wetland(nlon,nlat)      !percent wetland 
  real(r8) :: pct_lake(nlon,nlat)         !percent lake 
  real(r8) :: sumpctland(nlon,nlat)       !percent lake 

!cas additions for ccsm4 stuff

  real(r8) :: organic(nlon,nlat,nlay)     !organic soil density
  real(r8) :: fmax(nlon,nlat)             !maximum fractional saturated area
  real(r8) :: harvest_vh1(nlon,nlat)      !harvest from primary forest
  real(r8) :: harvest_vh2(nlon,nlat)      !harvest from primary non-forest 
  real(r8) :: harvest_sh1(nlon,nlat)      !harvest from secondary mature forest 
  real(r8) :: harvest_sh2(nlon,nlat)      !harvest from secondary young forest 
  real(r8) :: harvest_sh3(nlon,nlat)      !harvest from secondary non-forest 
  real(r8) :: grazing(nlon,nlat)          !grazing of herbacous pfts 
  real(r8) :: topo_bedrock(nlon,nlat)     !topography height (m) 
  real(r8) :: topo_ice(nlon,nlat)         !topography height  +  land ice (m)) 
  real(r8) :: ef_btr(nlon,nlat)           !micrograms isoprene (m-2h-1); broadleaf tree emission factor
  real(r8) :: ef_fet(nlon,nlat)           !""; fineleaf evergreen tree emission factor
  real(r8) :: ef_fdt(nlon,nlat)           !""; fineleaf deciduous tree emission factor
  real(r8) :: ef_shr(nlon,nlat)           !""; shrub emission factor
  real(r8) :: ef_grs(nlon,nlat)           !""; grass, non-vascular plants and other ground cover em 
  real(r8) :: ef_crp(nlon,nlat)           !""; crop emission factor 

  real(r8) :: topoin(nlon,nlat)	          !input topo veg array  
  real(r8) :: z_organic(nlat_org,nlay)	  !input zonal organic array  
  real(r8) :: zorg_lat(nlat_org)          !lats for organic soil density input


  integer :: i,j,m,ret,cnt                !indices
  integer :: ncid_soitex_i                !netCDF file id for soitex input
  integer :: ncid_sur_i                   !netCDF file id for surf type input"
  integer :: ncid_top_i                   !netCDF file id for topo input "
  integer :: ncid_zorg_i                  !netCDF file id for zonal organics input "
  integer :: ncid_glacier                 !netCDF file id for glacier output
  integer :: ncid_urban                   !netCDF file id for urban " 
  integer :: ncid_lanwat                  !netCDF file id for lanwat "
  integer :: ncid_lai                     !netCDF file id for lai "
  integer :: ncid_pft                     !netCDF file id for pft  "
  integer :: ncid_soicol                  !netCDF file id for soicol " 
  integer :: ncid_soitex                  !netCDF file id for soitex "
  integer :: ncid_organic                 !netCDF file id for organic "
  integer :: ncid_fmax                    !netCDF file id for fmax "
  integer :: ncid_vocef                   !netCDF file id for vocef "
  integer :: ncid_topo                    !netCDF file id for topo "

  integer :: dzsoi_id		          !dzsoi id
  integer :: zsoi_id		          !zsoi id
  integer :: pct_clay0_id	          !pct_clay0 id
  integer :: pct_sand0_id                 !pct_sand0 id
  integer :: mapunits0_id                 !mapunits id
  integer :: sur_id                       !surface type id
  integer :: top_id                       !input topo id
  integer :: z_organic_id                 !zonal organic id
  integer :: zorg_lat_id                  !zonal lat id

  integer :: veg(nlon,nlat)	          !input lsm veg array  


  character(len=80) :: filei_lsmveg,filei_soitex,filei_org,filei_topo  ! input filenames
  character(len=80) :: fileo_glacier              !output filenames
  character(len=80) :: fileo_urban                !output filenames
  character(len=80) :: fileo_lanwat               !output filenames
  character(len=80) :: fileo_lai                  !output filenames
  character(len=80) :: fileo_pft                  !output filenames
  character(len=80) :: fileo_soicol               !output filenames
  character(len=80) :: fileo_soitex               !output filenames
  character(len=80) :: fileo_organic              !output filenames
  character(len=80) :: fileo_fmax                 !output filenames
  character(len=80) :: fileo_vocef                !output filenames
  character(len=80) :: fileo_topo                 !output filenames
 
  character(len=78) :: char_string        !character string for ascii file


!-----------------------------------------------------------------

! Determine input and output file names

  filei_lsmveg =  'input_lsm_data'
  filei_topo   =  'input_top_data'
  filei_org    =  'input_org_data'
  filei_soitex =  'input_soi_data'

  fileo_glacier= 'output_glacier'
  fileo_urban  = 'output_urban' 
  fileo_lanwat = 'output_lanwat'
  fileo_lai    = 'output_lai'
  fileo_pft    = 'output_pft'
  fileo_soicol = 'output_soicol'
  fileo_soitex = 'output_soitex'
  fileo_organic= 'output_organic'
  fileo_fmax   = 'output_fmax'
  fileo_topo   = 'output_topo'
  fileo_vocef  = 'output_vocef'

  print *,'Output file names defined '


! -------------------------------------------------------------------
! Read in netcdf 2x2 land surface file
! note: data starts from north, need to flip lats...see below code
! -------------------------------------------------------------------
!
  ret = nf_open (filei_lsmveg, nf_nowrite, ncid_sur_i)
  if (ret == nf_noerr) then
    write(6,*)'Successfully opened netcdf surf ',trim(filei_lsmveg)
    call wrap_inq_varid (ncid_sur_i, 'SUR', sur_id   )
    call wrap_get_var_int (ncid_sur_i, sur_id, veg)
  else
    write(6,*)'cannot open surface type file successfully'
    call endrun
  endif
  ret = nf_close (ncid_sur_i)
  print *, 'Netcdf 2x2 Surface Type file read '

! -------------------------------------------------------------------
! Read in netcdf 2x2 topo file
! -------------------------------------------------------------------
!
  ret = nf_open (filei_topo, nf_nowrite, ncid_top_i)
  if (ret == nf_noerr) then
    write(6,*)'Successfully opened netcdf topo ',trim(filei_topo)
    call wrap_inq_varid (ncid_top_i, 'topo', top_id   )
    call wrap_get_var_int (ncid_top_i, top_id, topoin)
  else
    write(6,*)'cannot open topo type file successfully'
    call endrun
  endif
  ret = nf_close (ncid_top_i)
  print *, 'Netcdf 2x2 topo_depth file read '


! --------------------------------------------------------
! Define landmask
! --------------------------------------------------------

  do j = 1,nlat
  do i = 1,nlon
   if (veg(i,j)/=0) landmask(i,j) = 1.
  end do
  end do

! -------------------------------------------------------
! Define lat/lon 2d arrays
! -------------------------------------------------------
! Define North, West, South, East edges of grid
! nan + danlunt fall 2013

 edge(1) =   90.
 edge(2) =  360.
 edge(3) =  -90.
 edge(4) =    0.

! CCSM3/CLM3.5
! edge(1) =   90.
! edge(2) =  180.
! edge(3) =  -90.
! edge(4) = -180.

! Make latitudes and longitudes at center of grid cell

  dx = (edge(2)-edge(4)) / nlon
  dy = (edge(1)-edge(3)) / nlat

  do j = 1, nlat
     do i = 1, nlon
        latixy(i,j) = (edge(3)+dy/2.) + (j-1)*dy
        longxy(i,j) = (edge(4)+dx/2.) + (i-1)*dx
       end do
  end do

  lat(:) = latixy(1,:)
  lon(:) = longxy(:,1)

  print *, 'LAT= ',lat
  print *, 'LON= ',lon

  print *, 'Lat and Lon arrays defined '

! -------------------------------------------------------------------
! Read in netcdf soil texture file
! -------------------------------------------------------------------
!
  ret = nf_open (filei_soitex, nf_nowrite, ncid_soitex_i)
  if (ret == nf_noerr) then
    write(6,*)'Successfully opened netcdf soitex ',trim(filei_soitex)
    call wrap_inq_varid (ncid_soitex_i, 'DZSOI', dzsoi_id   )
    call wrap_inq_varid (ncid_soitex_i, 'ZSOI', zsoi_id   )
    call wrap_inq_varid (ncid_soitex_i, 'PCT_CLAY', pct_clay0_id   )
    call wrap_inq_varid (ncid_soitex_i, 'PCT_SAND', pct_sand0_id   )
    call wrap_inq_varid (ncid_soitex_i, 'MAPUNITS', mapunits0_id   )
    call wrap_get_var8 (ncid_soitex_i, dzsoi_id, dzsoi) 
    call wrap_get_var8 (ncid_soitex_i, zsoi_id, zsoi)
    call wrap_get_var8 (ncid_soitex_i, pct_clay0_id, pct_clay0)
    call wrap_get_var8 (ncid_soitex_i, pct_sand0_id, pct_sand0)
    call wrap_get_var8 (ncid_soitex_i, mapunits0_id, mapunits0)
  else
    write(6,*)'cannot open soil texture file successfully'
    call endrun
  endif
  ret = nf_close (ncid_soitex_i)
  print *, 'Netcdf IGBC Soil Texture file read '

! -------------------------------------------------------------------
! Read in netcdf zonal organic file
! -------------------------------------------------------------------
!
  ret = nf_open (filei_org, nf_nowrite, ncid_organic)
  if (ret == nf_noerr) then
    write(6,*)'Successfully opened netcdf zon organic ',trim(filei_org)
    call wrap_inq_varid (ncid_organic, 'z_ORGANIC', z_organic_id   )
    call wrap_inq_varid (ncid_organic, 'LAT', zorg_lat_id   )
    call wrap_get_var8 (ncid_organic, z_organic_id, z_organic)
    call wrap_get_var8 (ncid_organic, zorg_lat_id, zorg_lat)
  else
    write(6,*)'cannot open organic soil density file successfully'
    call endrun
  endif
  ret = nf_close (ncid_organic)
  print *, 'Netcdf zonal organic soil density file read '


! -------------------------------------------------------------------
! Call subroutines to create each raw data file
! --------------------------------------------------------------------

  call create_mksrf_glacier(veg,landmask,                         &
      		            nlon,nlat,lon,longxy,lat,latixy,edge, &
		            fileo_glacier,ncid_glacier,	          &
		            pct_glacier)
  print *, 'Netcdf File mksfr_glacier_paleo created'

  call create_mksrf_urban(veg,landmask,			          &
	 		  nlon,nlat,lon,longxy,lat,latixy,edge,   &
	                  fileo_urban,ncid_urban,		  &
		          pct_urban)
  print *, 'Netcdf File mksfr_urban_paleo created'

  call create_mksrf_lanwat(veg,landmask,		          &
	             	   nlon,nlat,lon,longxy,lat,latixy,edge,  &
  		           fileo_lanwat,ncid_lanwat,		  &
			   pct_wetland,pct_lake)
  print *, 'Netcdf File mksfr_lanwat_paleo created'

  call create_mksrf_pft(veg,landmask,numpft,			  &
	               nlon,nlat,lon,longxy,lat,latixy,edge,      &
                       pct_pft,harvest_vh1,harvest_vh2,harvest_sh1, &
                       harvest_sh2,harvest_sh3,grazing,          &
                       fileo_pft,ncid_pft) 
  print *, 'Netcdf File mksfr_pft_paleo created'

  call create_mksrf_lai(veg,landmask,numpft,		 	  &
                        nlon,nlat,lon,longxy,lat,latixy,edge,     &
 			fileo_lai,ncid_lai,                       &
                        pct_pft) 
  print *, 'Netcdf File mksfr_lai_paleo created'

  call create_mksrf_soicol(veg,landmask,			  &
			   nlon,nlat,lon,longxy,lat,latixy,edge , &
		           fileo_soicol,ncid_soicol)
  print *, 'Netcdf File mksfr_soicol_paleo created'

  call create_mksrf_soitex(veg,landmask,            	             &
		           dzsoi,zsoi,mapunits0,pct_clay0,pct_sand0, &
			   nlay,nlon_st,nlat_st,nmapunits,mapunitmax,&
		           nlon,nlat,lon,longxy,lat,latixy,edge,     &
                           fileo_soitex,ncid_soitex)
  print *, 'Netcdf File mksfr_soitex_paleo created'

  call create_mksrf_organic(veg,landmask,nlay,			  &
			   dzsoi,zsoi,z_organic,zorg_lat,nlat_org,&
			   nlon,nlat,lon,longxy,lat,latixy,edge , &
		           fileo_organic,ncid_organic)
  print *, 'Netcdf File mksfr_organic_paleo created'

  call create_mksrf_fmax(veg,landmask,			  &
			   nlon,nlat,lon,longxy,lat,latixy,edge , &
		           fileo_fmax,ncid_fmax)
  print *, 'Netcdf File mksfr_fmax_paleo created'

  call create_mksrf_topo(landmask,topoin,                               &
                               nlon,nlat,lon,longxy,lat,latixy,edge,     &
                               fileo_topo,ncid_topo)
  print *, 'Netcdf File mksfr_topo_paleo created'

  call create_mksrf_vocef(veg,landmask,		  &
			   nlon,nlat,lon,longxy,lat,latixy,edge , &
		           fileo_vocef,ncid_vocef)
  print *, 'Netcdf File mksfr_vocef_paleo created'

! error check 
! pct_glacier+pct_lake+pct_wetland+pct_urban and checks .le. 100.
  sumpctland = pct_glacier+pct_lake+pct_wetland+pct_urban 
  do j = 1, nlat
    do i = 1, nlon
      if (landmask(i,j)==1 .and. sumpctland(i,j)>100.) then
        write(*,*) 'ERROR: sumpctland (glacier+urban+wetland+lake)= ',&
	            sumpctland(i,j),i,j
      end if
    end do
  end do

  print *, 'End program' 

end program paleo_mkraw_ccsm4




!===============================================================================
subroutine create_mksrf_glacier(veg,landmask,			      &
		                nlon,nlat,lon,longxy,lat,latixy,edge, &
			        fileo,ncid)

  implicit none
  include 'netcdf.inc'


! ---------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer  :: veg(nlon,nlat)	          !lsm veg array  (input)
  real(r8) :: landmask(nlon,nlat)	  !landmask array  (input)
  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)  
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid

  integer :: nlon,nlat,ncid              !number lats/lons, nc fileo id 
  character(len=80) :: fileo             !output filenames

! ------------------------------------------------------------------
! Define local variables
! ------------------------------------------------------------------

  real(r8) :: pct_glacier(nlon,nlat)      !percent glacier

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id

  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: pct_glacier_id               !pct_glacier id
  integer :: landmask_id                  !landmask id

  integer :: i,j                          !indices
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: status                       !status

  character(len=256) :: name,unit            !netCDF attributes


! ------------------------------------------------------------
! Create skeleton netcdf
! --------------------------------------------------------------

! Define global attributes

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'pct_glacier_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat, dimlat_id)

! Define grid variables

  name = 'lon'
  unit = 'degrees east'
  dim1_id(1) = dimlon_id
  call wrap_def_var (ncid,'lon', nf_float, 1, dim1_id, lon_id)
  call wrap_put_att_text (ncid, lon_id, 'long_name', name)
  call wrap_put_att_text (ncid, lon_id, 'units'    , unit)

  name = 'lat'
  unit = 'degrees north'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'lat', nf_float, 1, dim1_id, lat_id)
  call wrap_put_att_text (ncid, lat_id, 'long_name', name)
  call wrap_put_att_text (ncid, lat_id, 'units'    , unit)

  name = 'longitude-2d'
  unit = 'degrees east'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)

  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)

  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)

  name = 'western edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)

! Define input file specific variables

  name = 'percent glacier'
  unit = 'unitless'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'PCT_GLACIER' ,nf_float, 2, dim2_id, pct_glacier_id)
  call wrap_put_att_text (ncid, pct_glacier_id, 'long_name', name)
  call wrap_put_att_text (ncid, pct_glacier_id, 'units'    , unit)

  name = 'land mask'
  unit = 'unitless'
  call wrap_def_var (ncid ,'LANDMASK' ,nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

! End of definition

  status = nf_enddef(ncid)

! -----------------------------------------------------------------------
! Create pct_glacier from LSM vegtypes 
! -----------------------------------------------------------------------

  do j = 1, nlat
    do i = 1, nlon
      if (veg(i,j).eq.1) then
        pct_glacier(i,j) = 100._r8
      else
        pct_glacier(i,j) = 0._r8
      end if
    end do
  end do

! --------------------------------------------------------------------------
! Write variables
! --------------------------------------------------------------------------

  call wrap_put_var_realx (ncid, lon_id        , lon)
  call wrap_put_var_realx (ncid, lat_id        , lat)
  call wrap_put_var_realx (ncid, longxy_id     , longxy)
  call wrap_put_var_realx (ncid, latixy_id     , latixy)
  call wrap_put_var_realx (ncid, edgen_id      , edge(1))
  call wrap_put_var_realx (ncid, edgee_id      , edge(2))
  call wrap_put_var_realx (ncid, edges_id      , edge(3))
  call wrap_put_var_realx (ncid, edgew_id      , edge(4))
  call wrap_put_var_realx (ncid, pct_glacier_id, pct_glacier)
  call wrap_put_var_realx (ncid, landmask_id   , landmask)

  call wrap_close(ncid)

end subroutine create_mksrf_glacier



!========================================================================
subroutine create_mksrf_urban(veg,landmask,      	            &
			      nlon,nlat,lon,longxy,lat,latixy,edge, &
			      fileo,ncid) 

  implicit none
  include 'netcdf.inc'

! ---------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer  :: veg(nlon,nlat)	          !lsm veg array  (input)
  real(r8) :: landmask(nlon,nlat)	  !landmask array  (input)
  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)  
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid

  integer :: nlon,nlat,ncid              !number lats/lons, nc fileo id 
  character(len=80) :: fileo             !output filenames

! ------------------------------------------------------------------
! Define local variables
! ------------------------------------------------------------------

  real(r8) :: pct_urban(nlon,nlat)        !pct urban

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id

  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: pct_urban_id                 !percent urban id
  integer :: landmask_id                  !landmask id

  integer :: i,j                          !indices
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: status                       !status

  character(len=256) :: name,unit         !netCDF attributes

! ------------------------------------------------------------
! Create skeleton netcdf
! --------------------------------------------------------------

! Define global attributes

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'pct_urban_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat, dimlat_id)

! Define input file independent variables 

  name = 'lon'
  unit = 'degrees east'
  dim1_id(1) = dimlon_id
  call wrap_def_var (ncid,'lon', nf_float, 1, dim1_id, lon_id)
  call wrap_put_att_text (ncid, lon_id, 'long_name', name)
  call wrap_put_att_text (ncid, lon_id, 'units'    , unit)

  name = 'lat'
  unit = 'degrees north'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'lat', nf_float, 1, dim1_id, lat_id)
  call wrap_put_att_text (ncid, lat_id, 'long_name', name)
  call wrap_put_att_text (ncid, lat_id, 'units'    , unit)

  name = 'longitude-2d'
  unit = 'degrees east'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)

  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)

  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)

  name = 'western edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)

! Define input file specific variables

  name = 'percent urban'
  unit = 'unitless'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'PCT_URBAN' ,nf_float, 2, dim2_id, pct_urban_id)
  call wrap_put_att_text (ncid, pct_urban_id, 'long_name', name)
  call wrap_put_att_text (ncid, pct_urban_id, 'units'    , unit)

  name = 'land mask'
  unit = 'unitless'
  call wrap_def_var (ncid ,'LANDMASK' ,nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

! End of definition

  status = nf_enddef(ncid)


! -----------------------------------------------------------------------
! Create pct_urban from LSM vegtypes
! -----------------------------------------------------------------------

  pct_urban = 0.

! --------------------------------------------------------------------------
! Write variables
! ----------------------------------------------------------------------

  call wrap_put_var_realx (ncid, lon_id      , lon)
  call wrap_put_var_realx (ncid, lat_id      , lat)
  call wrap_put_var_realx (ncid, longxy_id   , longxy)
  call wrap_put_var_realx (ncid, latixy_id   , latixy)
  call wrap_put_var_realx (ncid, edgen_id    , edge(1))
  call wrap_put_var_realx (ncid, edgee_id    , edge(2))
  call wrap_put_var_realx (ncid, edges_id    , edge(3))
  call wrap_put_var_realx (ncid, edgew_id    , edge(4))
  call wrap_put_var_realx (ncid, pct_urban_id, pct_urban)
  call wrap_put_var_realx (ncid, landmask_id, landmask)

  call wrap_close(ncid)

end subroutine create_mksrf_urban



!===========================================================================
subroutine create_mksrf_lanwat(veg,landmask,			     &
			       nlon,nlat,lon,longxy,lat,latixy,edge, &
                               fileo,ncid)
 

  implicit none
  include 'netcdf.inc'

! ---------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer  :: veg(nlon,nlat)	          !lsm veg array  (input)
  real(r8) :: landmask(nlon,nlat)	  !landmask array  (input)
  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)  
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid

  integer :: nlon,nlat,ncid              !number lats/lons, nc fileo id 
  character(len=80) :: fileo             !output filenames

! ------------------------------------------------------------------
! Define local variables
! ------------------------------------------------------------------


  real(r8) :: pct_lake(nlon,nlat)         !pct lake
  real(r8) :: pct_wetland(nlon,nlat)      !pct wetland

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id

  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: pct_lake_id                  !pct_lake id
  integer :: pct_wetland_id               !pct_wetland id
  integer :: landmask_id                  !landmask id

  integer :: i,j                          !indices
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: status                       !status

  character(len=256) :: name,unit         !netCDF attributes


! ------------------------------------------------------------
! Create skeleton netcdf
! --------------------------------------------------------------

! Define global attributes

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'lanwat_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat, dimlat_id)

! Define input file independent variables 

  name = 'lon'
  unit = 'degrees east'
  dim1_id(1) = dimlon_id
  call wrap_def_var (ncid,'lon', nf_float, 1, dim1_id, lon_id)
  call wrap_put_att_text (ncid, lon_id, 'long_name', name)
  call wrap_put_att_text (ncid, lon_id, 'units'    , unit)

  name = 'lat'
  unit = 'degrees north'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'lat', nf_float, 1, dim1_id, lat_id)
  call wrap_put_att_text (ncid, lat_id, 'long_name', name)
  call wrap_put_att_text (ncid, lat_id, 'units'    , unit)

  name = 'longitude-2d'
  unit = 'degrees east'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)

  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)

  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)

  name = 'western edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)

! Define percent lake and wetland variables

  name = 'percent lake'
  unit = 'unitless'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'PCT_LAKE' ,nf_float, 2, dim2_id, pct_lake_id)
  call wrap_put_att_text (ncid, pct_lake_id, 'long_name', name)
  call wrap_put_att_text (ncid, pct_lake_id, 'units'    , unit)

  name = 'percent wetland'
  unit = 'unitless'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'PCT_WETLAND' ,nf_float, 2, dim2_id, pct_wetland_id)
  call wrap_put_att_text (ncid, pct_wetland_id, 'long_name', name)
  call wrap_put_att_text (ncid, pct_wetland_id, 'units'    , unit)

  name = 'land mask'
  unit = 'unitless'
  call wrap_def_var (ncid ,'LANDMASK' ,nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

! End of definition

  status = nf_enddef(ncid)


! -----------------------------------------------------------------------
! Create pct_lake and pct_wetland from LSM vegtypes 
! -----------------------------------------------------------------------

 pct_lake = 0.

 do j = 1, nlat
   do i = 1, nlon
     if (veg(i,j).eq.27) then
       pct_wetland(i,j) = 20._r8
     elseif (veg(i,j).eq.28) then
       pct_wetland(i,j) = 100._r8
     else
       pct_wetland(i,j) = 0._r8
     end if
   end do
 end do

! --------------------------------------------------------------------------
! Write variables
! --------------------------------------------------------------------------


  call wrap_put_var_realx (ncid, lon_id        , lon)
  call wrap_put_var_realx (ncid, lat_id        , lat)
  call wrap_put_var_realx (ncid, longxy_id     , longxy)
  call wrap_put_var_realx (ncid, latixy_id     , latixy)
  call wrap_put_var_realx (ncid, edgen_id      , edge(1))
  call wrap_put_var_realx (ncid, edgee_id      , edge(2))
  call wrap_put_var_realx (ncid, edges_id      , edge(3))
  call wrap_put_var_realx (ncid, edgew_id      , edge(4))
  call wrap_put_var_realx (ncid, pct_lake_id   , pct_lake)
  call wrap_put_var_realx (ncid, pct_wetland_id, pct_wetland)
  call wrap_put_var_realx (ncid, landmask_id   , landmask)

  call wrap_close(ncid)

end subroutine create_mksrf_lanwat



!=============================================================================
subroutine create_mksrf_pft(veg,landmask,numpft,		  &
			    nlon,nlat,lon,longxy,lat,latixy,edge, &
                      pct_pft,harvest_vh1,harvest_vh2,harvest_sh1, &
                       harvest_sh2,harvest_sh3,grazing,            &
		       fileo,ncid)

  implicit none
  include 'netcdf.inc'

! ---------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer  :: veg(nlon,nlat)	          !lsm veg array  (input)
  real(r8) :: landmask(nlon,nlat)	  !landmask array  (input)
  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)  
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid

  integer :: nlon,nlat,numpft,ncid       !number lats/lons/pfts, nc fileo id 
  character(len=80) :: fileo             !output filenames

  real(r8) :: pct_pft(nlon,nlat,0:numpft) !percent pft 
  real(r8) :: harvest_vh1(nlon,nlat) !harvest from primary forest 
  real(r8) :: harvest_vh2(nlon,nlat) !harvest from primary non-forest 
  real(r8) :: harvest_sh1(nlon,nlat) !harvest from secondary mature forest 
  real(r8) :: harvest_sh2(nlon,nlat) !harvest from secondary young forest 
  real(r8) :: harvest_sh3(nlon,nlat) !harvest from secondary non-forest 
  real(r8) :: grazing(nlon,nlat)  !grazing of herbacous pfts 

! ------------------------------------------------------------------
! Define local variables
! ------------------------------------------------------------------

  real(r8) :: sumpctpft(nlon,nlat)        
  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id
  integer :: dimpft_id                    !netCDF dimension id

  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: landmask_id                  !landmask id
  integer :: pct_pft_id                   !pct_pft id
  integer :: harvest_vh1_id              !harvest id 
  integer :: harvest_vh2_id              !harvest id 
  integer :: harvest_sh1_id              !harvest id 
  integer :: harvest_sh2_id              !harvest id 
  integer :: harvest_sh3_id              !harvest id 
  integer :: grazing_id                   !grazing id 

  integer :: i,j                          !indices
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: dim3_id(3)                   !netCDF dimension id for 3-d variables
  integer :: status                       !status

  character(len=256) :: name,unit         !netCDF attributes
  
! ------------------------------------------------------------
! Create skeleton netcdf
! --------------------------------------------------------------

! Define global attributes

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'pft_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon', nlon    , dimlon_id)
  call wrap_def_dim (ncid, 'lat', nlat    , dimlat_id)
  call wrap_def_dim (ncid, 'pft', numpft+1, dimpft_id)

! Define grid variables

  name = 'lon'
  unit = 'degrees east'
  dim1_id(1) = dimlon_id
  call wrap_def_var (ncid,'lon', nf_float, 1, dim1_id, lon_id)
  call wrap_put_att_text (ncid, lon_id, 'long_name', name)
  call wrap_put_att_text (ncid, lon_id, 'units'    , unit)

  name = 'lat'
  unit = 'degrees north'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'lat', nf_float, 1, dim1_id, lat_id)
  call wrap_put_att_text (ncid, lat_id, 'long_name', name)
  call wrap_put_att_text (ncid, lat_id, 'units'    , unit)

  name = 'longitude-2d'
  unit = 'degrees east'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)

  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)

  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)

  
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)

! Define pft variables

  name = 'land mask'
  unit = 'unitless'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid ,'LANDMASK' ,nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

  name = 'percent pft'
  unit = 'unitless'
  dim3_id(1) = dimlon_id
  dim3_id(2) = dimlat_id
  dim3_id(3) = dimpft_id
  call wrap_def_var (ncid ,'PCT_PFT' ,nf_float, 3, dim3_id, pct_pft_id)
  call wrap_put_att_text (ncid, pct_pft_id, 'long_name', name)
  call wrap_put_att_text (ncid, pct_pft_id, 'units'    , unit)

  name = 'harvest from primary forest'
  unit = 'unitless'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid ,'HARVEST_VH1' ,nf_float, 2, dim2_id, harvest_vh1_id)
  call wrap_put_att_text (ncid, harvest_vh1_id, 'long_name', name)
  call wrap_put_att_text (ncid, harvest_vh1_id, 'units'    , unit)

  name = 'harvest from primary non-forest'
  unit = 'unitless'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid ,'HARVEST_VH2' ,nf_float, 2, dim2_id, harvest_vh2_id)
  call wrap_put_att_text (ncid, harvest_vh2_id, 'long_name', name)
  call wrap_put_att_text (ncid, harvest_vh2_id, 'units'    , unit)

  name = 'harvest from secondary marture-forest'
  unit = 'unitless'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid ,'HARVEST_SH1' ,nf_float, 2, dim2_id, harvest_sh1_id)
  call wrap_put_att_text (ncid, harvest_sh1_id, 'long_name', name)
  call wrap_put_att_text (ncid, harvest_sh1_id, 'units'    , unit)

  name = 'harvest from secondary young-forest'
  unit = 'unitless'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid ,'HARVEST_SH2' ,nf_float, 2, dim2_id, harvest_sh2_id)
  call wrap_put_att_text (ncid, harvest_sh2_id, 'long_name', name)
  call wrap_put_att_text (ncid, harvest_sh2_id, 'units'    , unit)

  name = 'harvest from secondary non-forest'
  unit = 'unitless'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid ,'HARVEST_SH3' ,nf_float, 2, dim2_id, harvest_sh3_id)
  call wrap_put_att_text (ncid, harvest_sh3_id, 'long_name', name)
  call wrap_put_att_text (ncid, harvest_sh3_id, 'units'    , unit)

  name = 'grazing of herbacous pfts'
  unit = 'unitless'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid ,'GRAZING' ,nf_float, 2, dim2_id, grazing_id)
  call wrap_put_att_text (ncid, grazing_id, 'long_name', name)
  call wrap_put_att_text (ncid, grazing_id, 'units'    , unit)

! End of definitions

  status = nf_enddef(ncid)

! -----------------------------------------------------------------------
! Create pfts from LSM vegtypes 
! -----------------------------------------------------------------------

  print *, 'inside pft'
  pct_pft = 0.                        !initialize
  do j = 1, nlat
    do i = 1, nlon
      if (veg(i,j).eq.1) then         !land ice
        pct_pft(i,j,0) = 100._r8      !is bare
        pct_pft(i,j,1:numpft) = 0._r8
      elseif (veg(i,j).eq.2) then     !desert
        pct_pft(i,j,0) = 100._r8      !is bare
      elseif (veg(i,j).eq.3) then     !cool needleleaf evergreen forest
        pct_pft(i,j,0) = 25._r8       !includes bare ground and
        pct_pft(i,j,2) = 75._r8       !needleleaf evergreen boreal tree
      elseif (veg(i,j).eq.4) then     !cool needleleaf deciduous forest
        pct_pft(i,j,0) = 50._r8       !includes bare ground and
        pct_pft(i,j,3) = 50._r8       !needleleaf deciduous boreal tree
      elseif (veg(i,j).eq.5) then     !cool broadleaf deciduous forest
        pct_pft(i,j,0) = 25._r8       !includes bare ground and
        pct_pft(i,j,8) = 75._r8       !broadleaf deciduous boreal tree
      elseif (veg(i,j).eq.6) then     !cool mixed forest
        pct_pft(i,j,0) = 26._r8       !includes bare ground and
        pct_pft(i,j,2) = 37._r8       !needleleaf evergreen boreal tree
        pct_pft(i,j,8) = 37._r8       !broadleaf deciduous boreal tree
      elseif (veg(i,j).eq.7) then     !warm needleleaf evergreen forest
        pct_pft(i,j,0) = 25._r8       !includes bare ground and
        pct_pft(i,j,1) = 75._r8       !needleleaf evergreen temperate tree
      elseif (veg(i,j).eq.8) then     !warm broadleaf deciduous forest
        pct_pft(i,j,0) = 25._r8       !includes bare ground and
        pct_pft(i,j,7) = 75._r8       !broadleaf deciduous temperate tree
      elseif (veg(i,j).eq.9) then     !warm mixed forest
        pct_pft(i,j,0) = 26._r8       !includes bare ground and
        pct_pft(i,j,1) = 37._r8       !needleleaf evergreen temperate tree
        pct_pft(i,j,7) = 37._r8       !broadleaf deciduous temperate tree
      elseif (veg(i,j).eq.10) then    !tropical broadleaf evergreen forest
        pct_pft(i,j,0) =  5._r8       !includes bare ground and
        pct_pft(i,j,4) = 95._r8       !broadleaf evergreen tropical tree
      elseif (veg(i,j).eq.11) then    !tropical seasonal deciduous forest
        pct_pft(i,j,0) = 25._r8       !includes bare ground and
        pct_pft(i,j,6) = 75._r8       !broadleaf deciduous tropical tree
      elseif (veg(i,j).eq.12) then    !savanna
        pct_pft(i,j,14) = 70._r8      !includes warm c4 grass and
        pct_pft(i,j,6) = 30._r8       !broadleaf deciduous tropical tree
      elseif (veg(i,j).eq.13) then    !evergreen forest tundra
        pct_pft(i,j,0) = 50._r8       !includes bare ground and
        pct_pft(i,j,2) = 25._r8       !needleleaf evergreen boreal tree
        pct_pft(i,j,12) = 25._r8      !arctic c3 grass
      elseif (veg(i,j).eq.14) then    !deciduous forest tundra
        pct_pft(i,j,0) = 50._r8       !includes bare ground and
        pct_pft(i,j,3) = 25._r8       !needleleaf deciduous boreal tree
        pct_pft(i,j,12) = 25._r8      !arctic c3 grass
      elseif (veg(i,j).eq.15) then    !cool forest crop
        pct_pft(i,j,15) = 40._r8      !crop
        pct_pft(i,j,2) = 30._r8       !needleleaf evergreen boreal tree
        pct_pft(i,j,8) = 30._r8       !broadleaf deciduous boreal tree
      elseif (veg(i,j).eq.16) then    !warm forest crop
        pct_pft(i,j,15) = 40._r8      !crop
        pct_pft(i,j,1) = 30._r8       !needleleaf evergreen temperate tree
        pct_pft(i,j,7) = 30._r8       !broadleaf deciduous temperate tree
      elseif (veg(i,j).eq.17) then    !cool grassland
        pct_pft(i,j,0) = 20._r8       !includes bare ground and
        pct_pft(i,j,13) = 60._r8      !cool c3 grass
        pct_pft(i,j,14) = 20._r8      !warm c4 grass
      elseif (veg(i,j).eq.18) then    !warm grassland
        pct_pft(i,j,0) = 20._r8       !includes bare ground and
        pct_pft(i,j,13) = 20._r8      !cool c3 grass
        pct_pft(i,j,14) = 60._r8      !warm c4 grass
      elseif (veg(i,j).eq.19) then    !tundra
        pct_pft(i,j,0) = 40._r8       !includes bare ground and
        pct_pft(i,j,11) = 30._r8      !broadleaf deciduous boreal shrub
        pct_pft(i,j,12) = 30._r8      !arctic c3 grass
      elseif (veg(i,j).eq.20) then    !evergreen shrubland
        pct_pft(i,j,0) = 20._r8       !includes bare ground and
        pct_pft(i,j,9) = 80._r8       !broadleaf evergreen temperate shrub
      elseif (veg(i,j).eq.21) then    !deciduous shrubland
        pct_pft(i,j,0) = 20._r8       !includes bare ground and
        pct_pft(i,j,10) = 80._r8      !broadleaf deciduous temperate shrub
      elseif (veg(i,j).eq.22) then    !semi-desert
        pct_pft(i,j,0) = 90._r8       !includes bare ground and
        pct_pft(i,j,10) = 10._r8      !broadleaf deciduous temperate shrub
      elseif (veg(i,j).eq.23) then    !cool irrigated crop
        pct_pft(i,j,0) = 15._r8       !includes bare ground and
        pct_pft(i,j,15) = 85._r8      !crop
      elseif (veg(i,j).eq.24) then    !cool crop
        pct_pft(i,j,0) = 15._r8       !includes bare ground and
        pct_pft(i,j,15) = 85._r8      !crop
      elseif (veg(i,j).eq.25) then    !warm irrigated crop
        pct_pft(i,j,0) = 15._r8       !includes bare ground and
        pct_pft(i,j,15) = 85._r8      !crop
      elseif (veg(i,j).eq.26) then    !warm crop
        pct_pft(i,j,0) = 15._r8       !includes bare ground and
        pct_pft(i,j,15) = 85._r8      !crop
      elseif (veg(i,j).eq.27) then    !forest wetland
        pct_pft(i,j,0) = 20._r8       !includes bare ground and
        pct_pft(i,j,5) = 80._r8       !broadleaf evergreen temperate tree
      elseif (veg(i,j).eq.28) then    !non-forest wetland
        pct_pft(i,j,0) = 100._r8      !is bare
      elseif (veg(i,j) < 0 .or. veg(i,j) > 28) then
        write(*,*) 'ERROR mapping veg to pct_pft: veg < 0 OR veg > 28',veg(i,j),i,j
      end if
    end do
  end do

! error check
  sumpctpft = sum(pct_pft, dim=3) !sum of %pft for every grid cell
  do j = 1, nlat
    do i = 1, nlon
      if (landmask(i,j)==1 .and. sumpctpft(i,j)/=100.) then
        write(*,*) 'ERROR: sumpctpft =',sumpctpft(i,j),i,j
      end if
    end do
  end do

!---------------------------------------------------------------------------
! for deep time assign harvest and grazing to zero
! change code here if this is not desired
!---------------------------------------------------------------------------

  do j = 1,nlon
  do i = 1,nlat
   if(landmask(j,i) == 1.) harvest_vh1(j,i) = 0._r8
   if(landmask(j,i) == 1.) harvest_vh2(j,i) = 0._r8
   if(landmask(j,i) == 1.) harvest_sh1(j,i) = 0._r8
   if(landmask(j,i) == 1.) harvest_sh2(j,i) = 0._r8
   if(landmask(j,i) == 1.) harvest_sh3(j,i) = 0._r8
   if(landmask(j,i) == 1.) grazing(j,i) = 0._r8
  enddo
  enddo
  print *, 'Harvest sample at point (1,1)',harvest_vh1(1,1)


! -------------------------------------------------------------------------- 
! Write variables
! --------------------------------------------------------------------------

  call wrap_put_var_realx (ncid, lon_id     , lon)
  call wrap_put_var_realx (ncid, lat_id     , lat)
  call wrap_put_var_realx (ncid, longxy_id  , longxy)
  call wrap_put_var_realx (ncid, latixy_id  , latixy)
  call wrap_put_var_realx (ncid, edgen_id   , edge(1))
  call wrap_put_var_realx (ncid, edgee_id   , edge(2))
  call wrap_put_var_realx (ncid, edges_id   , edge(3))
  call wrap_put_var_realx (ncid, edgew_id   , edge(4))
  call wrap_put_var_realx (ncid, landmask_id, landmask)
  call wrap_put_var_realx (ncid, pct_pft_id , pct_pft)
  call wrap_put_var_realx (ncid, harvest_vh1_id , harvest_vh1)
  call wrap_put_var_realx (ncid, harvest_vh2_id , harvest_vh2)
  call wrap_put_var_realx (ncid, harvest_sh1_id , harvest_sh1)
  call wrap_put_var_realx (ncid, harvest_sh2_id , harvest_sh2)
  call wrap_put_var_realx (ncid, harvest_sh3_id , harvest_sh3)
  call wrap_put_var_realx (ncid, grazing_id , grazing)

  call wrap_close(ncid)

end subroutine create_mksrf_pft



!============================================================================
subroutine create_mksrf_lai(veg,landmask,numpft,		  &
			    nlon,nlat,lon,longxy,lat,latixy,edge, &
                            fileo,ncid,                           &
                            pct_pft) 


  implicit none
  include 'netcdf.inc'

! ---------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer  :: veg(nlon,nlat)	          !lsm veg array  (input)
  real(r8) :: landmask(nlon,nlat)	  !landmask array  (input)
  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)  
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid


  integer :: nlon,nlat,numpft,ncid       !number lats/lons/pfts, nc fileo id 
  character(len=80) :: fileo             !output filenames

  real(r8) :: pct_pft(nlon,nlat,0:numpft) !percent pft 

! ------------------------------------------------------------------
! Define local variables
! ------------------------------------------------------------------

  real(r8) mlai (nlon,nlat,0:numpft) !monthly lai in
  real(r8) msai (nlon,nlat,0:numpft) !monthly sai in
  real(r8) mhgtt(nlon,nlat,0:numpft) !monthly height (top) in
  real(r8) mhgtb(nlon,nlat,0:numpft) !monthly height (bottom) in

  integer :: dimlon_id               !netCDF dimension id
  integer :: dimlat_id               !netCDF dimension id
  integer :: dimpft_id               !netCDF dimension id
  integer :: dimtim_id               !netCDF dimension id

  integer :: lon_id                  !1d longitude array id
  integer :: lat_id                  !1d latitude array id
  integer :: longxy_id               !2d longitude array id
  integer :: latixy_id               !2d latitude array id
  integer :: edgen_id                !northern edge of grid (edge(1)) id
  integer :: edgee_id                !eastern  edge of grid (edge(2)) id
  integer :: edges_id                !southern edge of grid (edge(3)) id
  integer :: edgew_id                !western  edge of grid (edge(4)) id
  integer :: landmask_id             !landmask id
  integer :: mlai_id                 !monthly mlai id
  integer :: msai_id                 !monthly msai id  
  integer :: mhgtt_id                !monthly mghtt id
  integer :: mhgtb_id                !monthly mhgtb id 
 
  integer :: ntim                    !month time index
  integer :: i,j,l                   !indices
  integer :: beg4d(4),len4d(4)       !netCDF edge
  integer :: dim1_id(1)              !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)              !netCDF dimension id for 2-d variables
  integer :: dim4_id(4)              !netCDF dimension id for 4-d variables
  integer :: status                  !status

  character(len=80) :: name,unit     !netCDF attributes

  integer :: month
  real(r8) gai(14,12)      !leaf area index (one sided)      |taken from LSM's
  real(r8) tai(14,12)      !leaf+stem area index (one sided) |subroutine
  real(r8) hvt(14)         !top of canopy (m)                |vegconi.F
  real(r8) hvb(14)         !bottom of canopy (m)             

  data (tai(1,i),i=1,12)  /4.5,4.7,5.0,5.1,5.3,5.5,5.3,5.3,5.2,4.9,4.6,4.5/
  data (tai(2,i),i=1,12)  /0.3,0.3,0.3,1.0,1.6,2.4,4.3,2.9,2.0,1.3,0.8,0.5/
  data (tai(3,i),i=1,12)  /5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0/
  data (tai(4,i),i=1,12)  /0.4,0.4,0.7,1.6,3.5,5.1,5.4,4.8,3.8,1.7,0.6,0.4/
  data (tai(5,i),i=1,12)  /1.2,1.0,0.9,0.8,0.8,1.0,2.0,3.7,3.2,2.7,1.9,1.2/
  data (tai(6,i),i=1,12)  /0.7,0.8,0.9,1.0,1.5,3.4,4.3,3.8,1.8,1.0,0.9,0.8/
  data (tai(7,i),i=1,12)  /1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3/
  data (tai(8,i),i=1,12)  /1.0,1.0,0.8,0.3,0.6,0.0,0.1,0.3,0.5,0.6,0.7,0.9/
  data (tai(9,i),i=1,12)  /0.1,0.1,0.1,0.1,0.1,0.3,1.5,1.7,1.4,0.1,0.1,0.1/
  data (tai(10,i),i=1,12) /0.7,0.8,0.9,1.0,1.5,3.4,4.3,3.8,1.8,1.0,0.9,0.8/
  data (tai(11,i),i=1,12) /0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,1.5,0.0,0.0,0.0/
  data (tai(12,i),i=1,12) /0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,1.5,0.0,0.0,0.0/
  data (tai(13,i),i=1,12) /0.7,0.8,0.9,1.0,1.5,3.4,4.3,3.8,1.8,1.0,0.9,0.8/
  data (tai(14,i),i=1,12) /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/

  data (gai(1,i),i=1,12)  /4.1,4.2,4.6,4.8,4.9,5.0,4.8,4.7,4.6,4.2,4.0,4.0/
  data (gai(2,i),i=1,12)  /0.0,0.0,0.0,0.6,1.2,2.0,2.6,1.7,1.0,0.5,0.2,0.0/
  data (gai(3,i),i=1,12)  /4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5/
  data (gai(4,i),i=1,12)  /0.0,0.0,0.3,1.2,3.0,4.7,4.5,3.4,1.2,0.3,0.0,0.0/
  data (gai(5,i),i=1,12)  /0.8,0.7,0.4,0.5,0.5,0.7,1.7,3.0,2.5,1.6,1.0,1.0/
  data (gai(6,i),i=1,12)  /0.4,0.5,0.6,0.7,1.2,3.0,3.5,1.5,0.7,0.6,0.5,0.4/
  data (gai(7,i),i=1,12)  /1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/
  data (gai(8,i),i=1,12)  /0.9,0.8,0.2,0.2,0.0,0.0,0.0,0.2,0.4,0.5,0.6,0.8/
  data (gai(9,i),i=1,12)  /0.0,0.0,0.0,0.0,0.0,0.2,1.4,1.2,0.0,0.0,0.0,0.0/
  data (gai(10,i),i=1,12) /0.4,0.5,0.6,0.7,1.2,3.0,3.5,1.5,0.7,0.6,0.5,0.4/
  data (gai(11,i),i=1,12) /0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,1.5,0.0,0.0,0.0/
  data (gai(12,i),i=1,12) /0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,1.5,0.0,0.0,0.0/
  data (gai(13,i),i=1,12) /0.4,0.5,0.6,0.7,1.2,3.0,3.5,1.5,0.7,0.6,0.5,0.4/
  data (gai(14,i),i=1,12) /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/

  data hvt /17.0,14.0,35.0,20.0,18.0, 0.5, 0.5, &
             0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0/

  data hvb / 8.50, 7.00, 1.00,11.50,10.00, 0.01, 0.10, &
             0.10, 0.10, 0.01, 0.01, 0.01, 0.01, 0.00/

! ------------------------------------------------------------
! Create skeleton netcdf
! --------------------------------------------------------------

! Define global attributes
 
  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'lai_sai_data')
  
! Define dimensions
  
  call wrap_def_dim (ncid, 'lon' , nlon        , dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat        , dimlat_id)
  call wrap_def_dim (ncid, 'pft' , numpft+1    , dimpft_id)
  call wrap_def_dim (ncid, 'time', nf_unlimited, dimtim_id)
  
! Define input file independent variables 
  
  name = 'lon'
  unit = 'degrees east'
  dim1_id(1) = dimlon_id
  call wrap_def_var (ncid,'lon', nf_float, 1, dim1_id, lon_id)
  call wrap_put_att_text (ncid, lon_id, 'long_name', name)
  call wrap_put_att_text (ncid, lon_id, 'units'    , unit)
  
  name = 'lat'
  unit = 'degrees north'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'lat', nf_float, 1, dim1_id, lat_id)
  call wrap_put_att_text (ncid, lat_id, 'long_name', name)
  call wrap_put_att_text (ncid, lat_id, 'units'    , unit)
  
  name = 'longitude-2d'
  unit = 'degrees east'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)
  
  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)
  
  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)
  
  name = 'western edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)
     
! Define input file specific variables
  
  name = 'land mask'
  unit = 'unitless'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LANDMASK', nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

  dim4_id(1) = dimlon_id
  dim4_id(2) = dimlat_id
  dim4_id(3) = dimpft_id
  dim4_id(4) = dimtim_id

  name = 'monthly leaf area index'
  unit = 'unitless'
  call wrap_def_var (ncid ,'MONTHLY_LAI', nf_float, 4, dim4_id, mlai_id)
  call wrap_put_att_text (ncid, mlai_id, 'long_name', name)
  call wrap_put_att_text (ncid, mlai_id, 'units'    , unit)
     
  name = 'monthly stem area index'
  unit = 'unitless'
  call wrap_def_var (ncid ,'MONTHLY_SAI', nf_float, 4, dim4_id, msai_id)
  call wrap_put_att_text (ncid, msai_id, 'long_name', name)
  call wrap_put_att_text (ncid, msai_id, 'units'    , unit)
     
  name = 'monthly height top'
  unit = 'meters'
  call wrap_def_var (ncid ,'MONTHLY_HEIGHT_TOP', nf_float, 4, dim4_id, mhgtt_id)
  call wrap_put_att_text (ncid, mhgtt_id, 'long_name', name)
  call wrap_put_att_text (ncid, mhgtt_id, 'units'    , unit)
  
  name = 'monthly height bottom'
  unit = 'meters'
  call wrap_def_var (ncid ,'MONTHLY_HEIGHT_BOT', nf_float, 4, dim4_id, mhgtb_id)
  call wrap_put_att_text (ncid, mhgtb_id, 'long_name', name)
  call wrap_put_att_text (ncid, mhgtb_id, 'units'    , unit)
     
  status = nf_enddef(ncid)

! -----------------------------------------------------------------------
! Create monthly_lai, monthly_sai, monthly_height_top, monthly_height_bot
! from LSM vegtypes, and LSM annual cycle lookup table. LSM subroutine
! vegconi.F contains lookup tables for the variables gai, tai, hvt, hvb:
! 
! dimensioned 14 LSM pfts by 12 months
! ------------------------------------
! gai is LAI
! tai is LAI+SAI, so SAI=tai-gai
! 
! dimensioned 14 LSM pfts (constant in time)
! ------------------------------------------
! hvt: height at top of canopy
! hvb: height at bottom of canopy 
! -----------------------------------------------------------------------

  mlai = 0. !initialize four variables globally
  msai = 0. !these values will be used for ocean and for bare ground
  mhgtt = 0.
  mhgtb = 0.

! --------------------------------------------------------------------------
! Write variables
! --------------------------------------------------------------------------

  call wrap_put_var_realx (ncid, lon_id     , lon)
  call wrap_put_var_realx (ncid, lat_id     , lat)
  call wrap_put_var_realx (ncid, longxy_id  , longxy)
  call wrap_put_var_realx (ncid, latixy_id  , latixy)
  call wrap_put_var_realx (ncid, edgen_id   , edge(1))
  call wrap_put_var_realx (ncid, edgee_id   , edge(2))
  call wrap_put_var_realx (ncid, edges_id   , edge(3))
  call wrap_put_var_realx (ncid, edgew_id   , edge(4))


! write out landmask

  call wrap_put_var_realx (ncid, landmask_id, landmask)

! now enter time loop

  do ntim = 1,12

  do j = 1,nlat !this loop written by slevis
    if (j <= nlat/2 .and. ntim <= 6) then
      month = ntim + 6
    else if (j <= nlat/2 .and. ntim > 6) then
      month = ntim - 6
    else
      month = ntim
    end if
    do i = 1, nlon
      if (pct_pft(i,j,1) > 0) then
        mlai(i,j,1) = gai(1,month)
        msai(i,j,1) = tai(1,month)-gai(1,month)
        mhgtt(i,j,1) = hvt(1)
        mhgtb(i,j,1) = hvb(1)
      end if
      if (pct_pft(i,j,2) > 0) then
        mlai(i,j,2) = gai(1,month)
        msai(i,j,2) = tai(1,month)-gai(1,month)
        mhgtt(i,j,2) = hvt(1)
        mhgtb(i,j,2) = hvb(1)
      end if
      if (pct_pft(i,j,3) > 0) then
        mlai(i,j,3) = gai(2,month)
        msai(i,j,3) = tai(2,month)-gai(2,month)
        mhgtt(i,j,3) = hvt(2)
        mhgtb(i,j,3) = hvb(2)
      end if
      if (pct_pft(i,j,4) > 0) then
        mlai(i,j,4) = gai(3,month)
        msai(i,j,4) = tai(3,month)-gai(3,month)
        mhgtt(i,j,4) = hvt(3)
        mhgtb(i,j,4) = hvb(3)
      end if
      if (pct_pft(i,j,5) > 0) then
        mlai(i,j,5) = gai(3,month)
        msai(i,j,5) = tai(3,month)-gai(3,month)
        mhgtt(i,j,5) = hvt(3)
        mhgtb(i,j,5) = hvb(3)
      end if
      if (pct_pft(i,j,6) > 0) then
        mlai(i,j,6) = gai(5,month)
        msai(i,j,6) = tai(5,month)-gai(5,month)
        mhgtt(i,j,6) = hvt(5)
        mhgtb(i,j,6) = hvb(5)
      end if
      if (pct_pft(i,j,7) > 0) then
        mlai(i,j,7) = gai(4,month)
        msai(i,j,7) = tai(4,month)-gai(4,month)
        mhgtt(i,j,7) = hvt(4)
        mhgtb(i,j,7) = hvb(4)
      end if
      if (pct_pft(i,j,8) > 0) then
        mlai(i,j,8) = gai(4,month)
        msai(i,j,8) = tai(4,month)-gai(4,month)
        mhgtt(i,j,8) = hvt(4)
        mhgtb(i,j,8) = hvb(4)
      end if
      if (pct_pft(i,j,9) > 0) then
        mlai(i,j,9) = gai(7,month)
        msai(i,j,9) = tai(7,month)-gai(7,month)
        mhgtt(i,j,9) = hvt(7)
        mhgtb(i,j,9) = hvb(7)
      end if
      if (pct_pft(i,j,10) > 0) then
        mlai(i,j,10) = gai(8,month)
        msai(i,j,10) = tai(8,month)-gai(8,month)
        mhgtt(i,j,10) = hvt(8)
        mhgtb(i,j,10) = hvb(8)
      end if
      if (pct_pft(i,j,11) > 0) then
        mlai(i,j,11) = gai(9,month)
        msai(i,j,11) = tai(9,month)-gai(9,month)
        mhgtt(i,j,11) = hvt(9)
        mhgtb(i,j,11) = hvb(9)
      end if
      if (pct_pft(i,j,12) > 0) then
        mlai(i,j,12) = gai(10,month)
        msai(i,j,12) = tai(10,month)-gai(10,month)
        mhgtt(i,j,12) = hvt(10)
        mhgtb(i,j,12) = hvb(10)
      end if
      if (pct_pft(i,j,13) > 0) then
        mlai(i,j,13) = gai(6,month)
        msai(i,j,13) = tai(6,month)-gai(6,month)
        mhgtt(i,j,13) = hvt(6)
        mhgtb(i,j,13) = hvb(6)
      end if
      if (pct_pft(i,j,14) > 0) then
        mlai(i,j,14) = gai(13,month)
        msai(i,j,14) = tai(13,month)-gai(13,month)
        mhgtt(i,j,14) = hvt(13)
        mhgtb(i,j,14) = hvb(13)
      end if
      if (pct_pft(i,j,15) > 0) then
        mlai(i,j,15) = gai(11,month)
        msai(i,j,15) = tai(11,month)-gai(11,month)
        mhgtt(i,j,15) = hvt(11)
        mhgtb(i,j,15) = hvb(11)
      end if
    end do
  end do !end slevis loops

  beg4d(1) = 1     ; len4d(1) = nlon
  beg4d(2) = 1     ; len4d(2) = nlat
  beg4d(3) = 1     ; len4d(3) = numpft+1
  beg4d(4) = ntim  ; len4d(4) = 1
     
  call wrap_put_vara_realx (ncid, mlai_id , beg4d, len4d, mlai )
  call wrap_put_vara_realx (ncid, msai_id , beg4d, len4d, msai )
  call wrap_put_vara_realx (ncid, mhgtt_id, beg4d, len4d, mhgtt)
  call wrap_put_vara_realx (ncid, mhgtb_id, beg4d, len4d, mhgtb)

  end do   

  ! Close output file

  call wrap_close(ncid)

end subroutine create_mksrf_lai


!=================================================================================
subroutine create_mksrf_soicol(veg,landmask,			     &
			       nlon,nlat,lon,longxy,lat,latixy,edge, &
                               fileo,ncid)


  implicit none
  include 'netcdf.inc'

! ---------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer  :: veg(nlon,nlat)	          !lsm veg array  (input)
  real(r8) :: landmask(nlon,nlat)	  !landmask array  (input)
  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)  
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid

  integer :: nlon,nlat,ncid              !number lats/lons, nc fileo id 
  character(len=80) :: fileo             !output filenames


! ------------------------------------------------------------------
! Define local variables
! ------------------------------------------------------------------

  real(r8) :: soil_color(nlon,nlat)       !lsm soil color

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id

  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: soil_color_id                !soil color id
  integer :: landmask_id                  !landmask id

  integer :: i,j                          !indices
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: status                       !status

  character(len=256) :: name,unit         !netCDF attributes


! ------------------------------------------------------------
! Create skeleton netcdf
! --------------------------------------------------------------

! Define global attributes

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'soil_color_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat, dimlat_id)

! Define input file independent variables 

  name = 'lon'
  unit = 'degrees east'
  dim1_id(1) = dimlon_id
  call wrap_def_var (ncid,'lon', nf_float, 1, dim1_id, lon_id)
  call wrap_put_att_text (ncid, lon_id, 'long_name', name)
  call wrap_put_att_text (ncid, lon_id, 'units'    , unit)

  name = 'lat'
  unit = 'degrees north'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'lat', nf_float, 1, dim1_id, lat_id)
  call wrap_put_att_text (ncid, lat_id, 'long_name', name)
  call wrap_put_att_text (ncid, lat_id, 'units'    , unit)

  name = 'longitude-2d'
  unit = 'degrees east'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)

  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)

  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)

  name = 'western edge of surface grid'
  
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)

! Define input file specific variables

  name = 'soil color'
  unit = 'unitless'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'SOIL_COLOR' ,nf_float, 2, dim2_id, soil_color_id)
  call wrap_put_att_text (ncid, soil_color_id, 'long_name', name)
  call wrap_put_att_text (ncid, soil_color_id, 'units'    , unit)

  name = 'land mask'
  unit = 'unitless'
  call wrap_def_var (ncid ,'LANDMASK' ,nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

! End of definition

  status = nf_enddef(ncid)

! -----------------------------------------------------------------------
! Create soil_color from LSM vegtypes (if user so desires)
! -----------------------------------------------------------------------

  do j = 1,nlon
  do i = 1,nlat
   if(landmask(j,i) == 1.) soil_color(j,i) = 10._r8   ! see Table 3.3. in CLM4 doc to choose your colr
						      ! http://www.cesm.ucar.edu/models/cesm1.0/clm/CLM4_Tech_Note.pdf
						      ! or add code here to modify soil color by geography 
   !! if(landmask(j,i) == 1.) soil_color(j,i) = 15._r8   ! used to be 4 for ccsm3
  enddo
  enddo
  print *, 'Soil Color sample at point (1,1)',soil_color(1,1)

! --------------------------------------------------------------------------
! Write variables
! --------------------------------------------------------------------------

  call wrap_put_var_realx (ncid, lon_id       , lon)
  call wrap_put_var_realx (ncid, lat_id       , lat)
  call wrap_put_var_realx (ncid, longxy_id    , longxy)
  call wrap_put_var_realx (ncid, latixy_id    , latixy)
  call wrap_put_var_realx (ncid, edgen_id     , edge(1))
  call wrap_put_var_realx (ncid, edgee_id     , edge(2))
  call wrap_put_var_realx (ncid, edges_id     , edge(3))
  call wrap_put_var_realx (ncid, edgew_id     , edge(4))
  call wrap_put_var_realx (ncid, landmask_id  , landmask)
  call wrap_put_var_realx (ncid, soil_color_id, soil_color)

! Close output file

  call wrap_close(ncid)

end subroutine create_mksrf_soicol



!===============================================================================
subroutine create_mksrf_soitex(veg,landmask,     		         &
			       dzsoi,zsoi,mapunits0,pct_clay0,pct_sand0, &
		               nlay,nlon_st,nlat_st,nmapunits,mapunitmax,&
		               nlon,nlat,lon,longxy,lat,latixy,edge,     &
                               fileo,ncid)


  implicit none
  include 'netcdf.inc'

! ---------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer  :: veg(nlon,nlat)	          !lsm veg array  (input)
  real(r8) :: landmask(nlon,nlat)	  !landmask array  (input)
  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)  
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid

  real(r8) :: dzsoi(10), zsoi(10)        !soil layer thickness and depth
  real(r8) :: pct_sand0(nlay,mapunitmax)    !original percent sand 
  real(r8) :: pct_clay0(nlay,mapunitmax)    !original percent clay 
  real(r8) :: mapunits0(nlon_st,nlat_st)    !mapunits 
  

  integer :: nlon_st,nlat_st,nlay	 !original soil tex dims
  integer :: nmapunits,mapunitmax        !# igbp mapunits, max value mu 
  integer :: nlon,nlat,ncid 		 !number lats/lons, nc fileo id 
  character(len=80) :: fileo             !output filenames


! ------------------------------------------------------------------
! Define local variables
! ------------------------------------------------------------------


  integer  :: countloam                   !layers of loam per mapunit
  integer  :: mu                          !current mapunit
  real(r8) :: pct_sand(mapunitmax,nlay)   !pct sand 
  real(r8) :: pct_clay(mapunitmax,nlay)   !pct clay
  real(r8) :: mapunits(nlon,nlat)         !new mapunits

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id
  integer :: dimlay_id                    !netCDF dimension id
  integer :: dimmapunits_id               !netCDF dimension id
  integer :: dimmapunitmax_id             !netCDF dimension id

  integer :: dzsoi_id                     !soil thickness by layer
  integer :: zsoi_id                      !soil depth by layer
  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: lay_id                       !1d layer array id
  integer :: mapunit_id                   !2d mapunits array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: pct_sand_id                  !sand id
  integer :: pct_clay_id                  !clay id
  integer :: landmask_id                  !landmask id

  integer :: i,j,k                        !indices
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: status                       !status

  character(len=256) :: name,unit         !netCDF attributes


! ------------------------------------------------------------
! Create skeleton netcdf
! --------------------------------------------------------------

! Define global attributes

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'igbp_soil_texture_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat, dimlat_id)
  call wrap_def_dim (ncid, 'number_of_layers'   , nlay      , dimlay_id)
  call wrap_def_dim (ncid, 'number_of_mapunits' , nmapunits , dimmapunits_id)
  call wrap_def_dim (ncid, 'max_value_mapunit'  , mapunitmax, dimmapunitmax_id)

! Define input file independent variables 

  name = 'lon'
  unit = 'degrees east'
  dim1_id(1) = dimlon_id
  call wrap_def_var (ncid,'lon', nf_float, 1, dim1_id, lon_id)
  call wrap_put_att_text (ncid, lon_id, 'long_name', name)
  call wrap_put_att_text (ncid, lon_id, 'units'    , unit)

  name = 'lat'
  unit = 'degrees north'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'lat', nf_float, 1, dim1_id, lat_id)
  call wrap_put_att_text (ncid, lat_id, 'long_name', name)
  call wrap_put_att_text (ncid, lat_id, 'units'    , unit)

  name = 'longitude-2d'
  unit = 'degrees east'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'land mask'
  unit = 'unitless'
  call wrap_def_var (ncid ,'LANDMASK' ,nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

! to possibly replace the next two variables
! find out about dimensioned variables
! (eg, see how pressure levels are treated)

  name = 'soil layer thickness'
  unit = 'm'
  dim1_id(1) = dimlay_id
  call wrap_def_var (ncid,'DZSOI', nf_float, 1, dim1_id, dzsoi_id)
  call wrap_put_att_text (ncid, dzsoi_id, 'long_name', name)
  call wrap_put_att_text (ncid, dzsoi_id, 'units'    , unit)

  name = 'soil layer depth'
  unit = 'm'
  dim1_id(1) = dimlay_id
  call wrap_def_var (ncid,'ZSOI', nf_float, 1, dim1_id, zsoi_id)
  call wrap_put_att_text (ncid, zsoi_id, 'long_name', name)
  call wrap_put_att_text (ncid, zsoi_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)

  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)

  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)

  name = 'western edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)

! Define soil type and soil texture variables

  name = 'igbp soil mapunit'
  unit = 'unitless'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid ,'MAPUNITS' ,nf_float, 2, dim2_id, mapunit_id)
  call wrap_put_att_text (ncid, mapunit_id, 'long_name', name)
  call wrap_put_att_text (ncid, mapunit_id, 'units'    , unit)

  name = 'percent sand'
  unit = 'unitless'
  dim2_id(1) = dimmapunitmax_id
  dim2_id(2) = dimlay_id
  call wrap_def_var (ncid ,'PCT_SAND' ,nf_float, 2, dim2_id, pct_sand_id)
  call wrap_put_att_text (ncid, pct_sand_id, 'long_name', name)
  call wrap_put_att_text (ncid, pct_sand_id, 'units'    , unit)

  name = 'percent clay'
  unit = 'unitless'
  call wrap_def_var (ncid ,'PCT_CLAY' ,nf_float, 2, dim2_id, pct_clay_id)
  call wrap_put_att_text (ncid, pct_clay_id, 'long_name', name)
  call wrap_put_att_text (ncid, pct_clay_id, 'units'    , unit)

! End of definition

  status = nf_enddef(ncid)


! -----------------------------------------------------------------------
! pct_clay and pct_sand data 
! Think of mapunits as soil profiles. Each profile corresponds to a unique
! set of %sand and %clay values from the surface layer to the deepest one.
!
! Here we look for the mapunit that corresponds to a loam (an average soil)
! for all layers. We assign that mapunit to all land points.
! Users with more extensive knowledge of the soil texture for the period
! of their interest, should expand this code accordingly.
! -----------------------------------------------------------------------

  pct_sand = 0.  !initialize
  pct_clay = 0.  !three
  mapunits = 0.  !output variables

  do mu = 1,mapunitmax
    countloam = 0
    do k = 1,nlay
      if (pct_sand0(k,mu) > 39 .and. pct_sand0(k,mu) < 47 .and. &
          pct_clay0(k,mu) > 14 .and. pct_clay0(k,mu) < 22) then
        countloam = countloam + 1
      end if
      if (countloam == 10) then  !all layers are loam
        write(*,*) 'Found loam to be mapunit:',mu
        do j = 1,nlat
          do i = 1,nlon
            if (landmask(i,j) == 1) then
              mapunits(i,j) = mu !set mapunits globally to current mu (loam)
            end if
          end do !lon loop
        end do   !lat loop
      end if
      pct_sand(mu,k) = pct_sand0(k,mu) !the mapping from mapunit to %sand
      pct_clay(mu,k) = pct_clay0(k,mu) !and %clay remains the same
    end do       !loop through soil layers
  end do         !loop through mapunits

! --------------------------------------------------------------------------
! Write variables
! --------------------------------------------------------------------------

  call wrap_put_var_realx (ncid, lon_id     , lon)
  call wrap_put_var_realx (ncid, lat_id     , lat)
  call wrap_put_var_realx (ncid, longxy_id  , longxy)
  call wrap_put_var_realx (ncid, latixy_id  , latixy)
  call wrap_put_var_realx (ncid, landmask_id, landmask)
  call wrap_put_var_realx (ncid, edgen_id   , edge(1))
  call wrap_put_var_realx (ncid, edgee_id   , edge(2))
  call wrap_put_var_realx (ncid, edges_id   , edge(3))
  call wrap_put_var_realx (ncid, edgew_id   , edge(4))
  call wrap_put_var_realx (ncid, dzsoi_id   , dzsoi)
  call wrap_put_var_realx (ncid, zsoi_id    , zsoi)
  call wrap_put_var_realx (ncid, mapunit_id , mapunits)
  call wrap_put_var_realx (ncid, pct_sand_id, pct_sand)
  call wrap_put_var_realx (ncid, pct_clay_id, pct_clay)

  call wrap_close(ncid)

end subroutine create_mksrf_soitex 


!===============================================================================
subroutine create_mksrf_organic(veg,landmask,nlay,                       &
                               dzsoi,zsoi,z_organic,zorg_lat,nlat_org,	 &
                               nlon,nlat,lon,longxy,lat,latixy,edge,     &
                               fileo,ncid)


  implicit none
  include 'netcdf.inc'

! ---------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer  :: veg(nlon,nlat)	          !lsm veg array  (input)
  real(r8) :: landmask(nlon,nlat)	  !landmask array  (input)
  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)  
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid

  real(r8) :: dzsoi(10), zsoi(10)        !soil layer thickness and depth
  real(r8) :: z_organic(nlat_org,10)     !zonal organic soil density values 
  real(r8) :: zorg_lat(nlat_org)         !lat for zonal file

  integer :: nlon,nlat,nlat_org,nlay,ncid  !number lats/lons, nc fileo id 
  character(len=80) :: fileo             !output filenames


! ------------------------------------------------------------------
! Define local variables
! ------------------------------------------------------------------


  real(r8) :: organic(nlon,nlat,10) 	  !organic soil density
  real(r8) :: sum80,sum70,sum60,sum50,sum40,sum30,sum20,sum10,sum0
  real(r8) :: summ0,summ10,summ20,summ30,summ40,summ50,summ60,summ70,summ80,summ90
  real(r8) :: cnt80,cnt70,cnt60,cnt50,cnt40,cnt30,cnt20,cnt10,cnt0,cntm0,cntm10
  real(r8) :: cntm20,cntm30,cntm40,cntm50,cntm60,cntm70,cntm80,cntm90

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id
  integer :: dimlay_id                    !netCDF dimension id

  integer :: dzsoi_id                     !soil thickness by layer
  integer :: zsoi_id                      !soil depth by layer
  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: lay_id                       !1d layer array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: organic_id                   !organic id
  integer :: landmask_id                  !landmask id

  integer :: i,j,k,jj                     !indices
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: dim3_id(3)                   !netCDF dimension id for 2-d variables
  integer :: status                       !status

  character(len=256) :: name,unit         !netCDF attributes


! ------------------------------------------------------------
! Create skeleton netcdf
! --------------------------------------------------------------

! Define global attributes

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'organic_soil_density')

! Define dimensions

  call wrap_def_dim (ncid, 'lon' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat, dimlat_id)
  call wrap_def_dim (ncid, 'number_of_layers'   , nlay      , dimlay_id)

! Define input file independent variables 

  name = 'lon'
  unit = 'degrees east'
  dim1_id(1) = dimlon_id
  call wrap_def_var (ncid,'lon', nf_float, 1, dim1_id, lon_id)
  call wrap_put_att_text (ncid, lon_id, 'long_name', name)
  call wrap_put_att_text (ncid, lon_id, 'units'    , unit)

  name = 'lat'
  unit = 'degrees north'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'lat', nf_float, 1, dim1_id, lat_id)
  call wrap_put_att_text (ncid, lat_id, 'long_name', name)
  call wrap_put_att_text (ncid, lat_id, 'units'    , unit)

  name = 'longitude-2d'
  unit = 'degrees east'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'land mask'
  unit = 'unitless'
  call wrap_def_var (ncid ,'LANDMASK' ,nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

  name = 'soil layer thickness'
  unit = 'm'
  dim1_id(1) = dimlay_id
  call wrap_def_var (ncid,'DZSOI', nf_float, 1, dim1_id, dzsoi_id)
  call wrap_put_att_text (ncid, dzsoi_id, 'long_name', name)
  call wrap_put_att_text (ncid, dzsoi_id, 'units'    , unit)

  name = 'soil layer depth'
  unit = 'm'
  dim1_id(1) = dimlay_id
  call wrap_def_var (ncid,'ZSOI', nf_float, 1, dim1_id, zsoi_id)
  call wrap_put_att_text (ncid, zsoi_id, 'long_name', name)
  call wrap_put_att_text (ncid, zsoi_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)

  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)

  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)

  name = 'western edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)

  name = 'Organic soil density at soil levels'
  unit = 'kg/m3 (assumed carbon content 0.58 gC per gOM)'
  dim3_id(1) = dimlon_id
  dim3_id(2) = dimlat_id
  dim3_id(3) = dimlay_id
  call wrap_def_var (ncid ,'ORGANIC' ,nf_float, 3, dim3_id, organic_id)
  call wrap_put_att_text (ncid, organic_id, 'long_name', name)
  call wrap_put_att_text (ncid, organic_id, 'units'    , unit)

! End of definition

  status = nf_enddef(ncid)


! -----------------------------------------------------------------------
! As a first guess for any deep time case, we assume the organic soil
! density to approximate the zonal average modern values. Essentially,
! this assumes the latitudinal structure of organic soils, i.e. pole to
! equator stays the same.  If this is not desirable, change code here.
! -----------------------------------------------------------------------

 organic = 0.  ! initialize all point to zero

! clunky...feel free to change
    ! do k = 1,nlay
      sum80 = 0.
      sum70 = 0.
      sum60 = 0.
      sum50 = 0.
      sum40 = 0.
      sum30 = 0.
      sum20 = 0.
      sum10 = 0.
      sum0 = 0.
      summ10 = 0.
      summ20 = 0.
      summ30 = 0.
      summ40 = 0.
      summ50 = 0.
      summ60 = 0.
      summ70 = 0.
      summ80 = 0.
      summ90 = 0.
	cnt80 = 0.
	cnt70 = 0.
	cnt60 = 0.
	cnt50 = 0.
	cnt40 = 0.
	cnt30 = 0.
	cnt20 = 0.
	cnt10 = 0.
	cnt0 = 0.
	cntm10 = 0.
	cntm20 = 0.
	cntm30 = 0.
	cntm40 = 0.
	cntm50 = 0.
	cntm60 = 0.
	cntm70 = 0.
	cntm80 = 0.
	cntm90 = 0.
        do j = 1,nlat
          do i = 1,nlon
            if (landmask(i,j) == 1) then
		do jj = 1,nlat_org
 		 if(lat(j).ge.80.and.zorg_lat(jj).ge.80)then
                   sum80 = z_organic(jj,k) + sum80 
                   cnt80 = cnt80 + 1.
                 elseif (lat(j).ge.70.and.lat(j).lt.80.and.zorg_lat(jj).ge.70.and.zorg_lat(jj).lt.80)then
                   sum70 = z_organic(jj,k) + sum70 
                   cnt70 = cnt70 + 1.
                 elseif (lat(j).ge.60.and.lat(j).lt.70.and.zorg_lat(jj).ge.60.and.zorg_lat(jj).lt.70)then
                   sum60 = z_organic(jj,k) + sum60 
                   cnt60 = cnt60 + 1.
                 elseif (lat(j).ge.50.and.lat(j).lt.60.and.zorg_lat(jj).ge.50.and.zorg_lat(jj).lt.60)then
                   sum50 = z_organic(jj,k) + sum50 
                   cnt50 = cnt50 + 1.
                 elseif (lat(j).ge.40.and.lat(j).lt.50.and.zorg_lat(jj).ge.40.and.zorg_lat(jj).lt.50)then
                   sum40 = z_organic(jj,k) + sum40 
                   cnt40 = cnt40 + 1.
                 elseif (lat(j).ge.30.and.lat(j).lt.40.and.zorg_lat(jj).ge.30.and.zorg_lat(jj).lt.40)then
                   sum30 = z_organic(jj,k) + sum30 
                   cnt30 = cnt30 + 1.
                 elseif (lat(j).ge.20.and.lat(j).lt.30.and.zorg_lat(jj).ge.20.and.zorg_lat(jj).lt.30)then
                   sum20 = z_organic(jj,k) + sum20 
                   cnt20 = cnt20 + 1.
                 elseif (lat(j).ge.10.and.lat(j).lt.20.and.zorg_lat(jj).ge.10.and.zorg_lat(jj).lt.20)then
                   sum10 = z_organic(jj,k) + sum10 
                   cnt10 = cnt10 + 1.
                 elseif (lat(j).ge.0.and.lat(j).lt.10.and.zorg_lat(jj).ge.0.and.zorg_lat(jj).lt.10)then
                   sum0 = z_organic(jj,k) + sum0 
                   cnt0 = cnt0 + 1.
                 elseif (lat(j).ge.-10.and.lat(j).lt.0.and.zorg_lat(jj).ge.-10.and.zorg_lat(jj).lt.0)then
                   summ10 = z_organic(jj,k) + summ10 
                   cntm10 = cntm10 + 1.
                 elseif (lat(j).ge.-20.and.lat(j).lt.-10.and.zorg_lat(jj).ge.-20.and.zorg_lat(jj).lt.-10)then
                   summ20 = z_organic(jj,k) + summ20 
                   cntm20 = cntm20 + 1.
                 elseif (lat(j).ge.-30.and.lat(j).lt.-20.and.zorg_lat(jj).ge.-30.and.zorg_lat(jj).lt.-20)then
                   summ30 = z_organic(jj,k) + summ30 
                   cntm30 = cntm30 + 1.
                 elseif (lat(j).ge.-40.and.lat(j).lt.-30.and.zorg_lat(jj).ge.-40.and.zorg_lat(jj).lt.-30)then
                   summ40 = z_organic(jj,k) + summ40 
                   cntm40 = cntm40 + 1.
                 elseif (lat(j).ge.-50.and.lat(j).lt.-40.and.zorg_lat(jj).ge.-50.and.zorg_lat(jj).lt.-40)then
                   summ50 = z_organic(jj,k) + summ50 
                   cntm50 = cntm50 + 1.
                 elseif (lat(j).ge.-60.and.lat(j).lt.-50.and.zorg_lat(jj).ge.-60.and.zorg_lat(jj).lt.-50)then
                   summ60 = z_organic(jj,k) + summ60 
                   cntm60 = cntm60 + 1.
                 elseif (lat(j).ge.-70.and.lat(j).lt.-60.and.zorg_lat(jj).ge.-70.and.zorg_lat(jj).lt.-60)then
                   summ70 = z_organic(jj,k) + summ70 
                   cntm70 = cntm70 + 1.
                 elseif (lat(j).ge.-80.and.lat(j).lt.-70.and.zorg_lat(jj).ge.-80.and.zorg_lat(jj).lt.-70)then
                   summ80 = z_organic(jj,k) + summ80 
                   cntm80 = cntm80 + 1.
                 elseif (lat(j).ge.-90.and.lat(j).lt.-80.and.zorg_lat(jj).ge.-90.and.zorg_lat(jj).lt.-80)then
                   summ90 = z_organic(jj,k) + summ90 
                   cntm90 = cntm90 + 1.
	         end if
	       end do
	    end if 
	 end do
	end do
        sum80 = sum80/cnt80     
        sum70 = sum70/cnt70     
        sum60 = sum60/cnt60     
        sum50 = sum50/cnt50     
        sum40 = sum40/cnt40     
        sum30 = sum30/cnt30     
        sum20 = sum20/cnt20     
        sum10 = sum10/cnt10     
        sum0 = sum0/cnt0     
        summ10 = summ10/cntm10     
        summ20 = summ20/cntm20     
        summ30 = summ30/cntm30     
        summ40 = summ40/cntm40     
        summ50 = summ50/cntm50     
        summ60 = summ60/cntm60     
        summ70 = summ70/cntm70     
        summ80 = summ80/cntm80     
        summ90 = summ90/cntm90     
        do j = 1,nlat
          do i = 1,nlon
            if (landmask(i,j) == 1) then
 		 if(lat(j).ge.80)then
		  organic(i,j,k) = sum80
                 elseif (lat(j).ge.70.and.lat(j).lt.80)then
		  organic(i,j,k) = sum70
                 elseif (lat(j).ge.60.and.lat(j).lt.70)then
		  organic(i,j,k) = sum60
                 elseif (lat(j).ge.50.and.lat(j).lt.60)then
		  organic(i,j,k) = sum50
                 elseif (lat(j).ge.40.and.lat(j).lt.50)then
		  organic(i,j,k) = sum40
                 elseif (lat(j).ge.30.and.lat(j).lt.40)then
		  organic(i,j,k) = sum30
                 elseif (lat(j).ge.20.and.lat(j).lt.30)then
		  organic(i,j,k) = sum20
                 elseif (lat(j).ge.10.and.lat(j).lt.20)then
		  organic(i,j,k) = sum10
                 elseif (lat(j).ge.0.and.lat(j).lt.10)then
		  organic(i,j,k) = sum0
                 elseif (lat(j).ge.-10.and.lat(j).lt.0)then
		  organic(i,j,k) = summ10
                 elseif (lat(j).ge.-20.and.lat(j).lt.-10)then
		  organic(i,j,k) = summ20
                 elseif (lat(j).ge.-30.and.lat(j).lt.-20)then
		  organic(i,j,k) = summ30
                 elseif (lat(j).ge.-40.and.lat(j).lt.-30)then
		  organic(i,j,k) = summ40
                 elseif (lat(j).ge.-50.and.lat(j).lt.-40)then
		  organic(i,j,k) = summ50
                 elseif (lat(j).ge.-60.and.lat(j).lt.-50)then
		  organic(i,j,k) = summ60
                 elseif (lat(j).ge.-70.and.lat(j).lt.-60)then
		  organic(i,j,k) = summ70
                 elseif (lat(j).ge.-80.and.lat(j).lt.-70)then
		  organic(i,j,k) = summ80
                 elseif (lat(j).ge.-90.and.lat(j).lt.-80)then
		  organic(i,j,k) = summ90
	        endif
               endif
	     enddo
	  enddo
        print *,' k =  ', k
   	print *, 'sum80 = ', sum80
   	print *, 'sum70 = ', sum70
   	print *, 'sum60 = ', sum60
   	print *, 'sum50 = ', sum50
   	print *, 'sum40 = ', sum40
   	print *, 'sum30 = ', sum30
   	print *, 'sum20 = ', sum20
   	print *, 'sum10 = ', sum10
   	print *, 'sum0 = ', sum0
   	print *, 'summ10 = ', summ10
   	print *, 'summ20 = ', summ20
   	print *, 'summ30 = ', summ30
   	print *, 'summ40 = ', summ40
   	print *, 'summ50 = ', summ50
   	print *, 'summ60 = ', summ60
   	print *, 'summ70 = ', summ70
   	print *, 'summ80 = ', summ80
   	print *, 'summ90 = ', summ90


! --------------------------------------------------------------------------
! Write variables
! --------------------------------------------------------------------------

  call wrap_put_var_realx (ncid, lon_id     , lon)
  call wrap_put_var_realx (ncid, lat_id     , lat)
  call wrap_put_var_realx (ncid, longxy_id  , longxy)
  call wrap_put_var_realx (ncid, latixy_id  , latixy)
  call wrap_put_var_realx (ncid, landmask_id, landmask)
  call wrap_put_var_realx (ncid, edgen_id   , edge(1))
  call wrap_put_var_realx (ncid, edgee_id   , edge(2))
  call wrap_put_var_realx (ncid, edges_id   , edge(3))
  call wrap_put_var_realx (ncid, edgew_id   , edge(4))
  call wrap_put_var_realx (ncid, dzsoi_id   , dzsoi)
  call wrap_put_var_realx (ncid, zsoi_id    , zsoi)
  call wrap_put_var_realx (ncid, organic_id, organic)

  call wrap_close(ncid)


! ---------------------------------------------------------------------

end subroutine create_mksrf_organic



!===============================================================================
subroutine create_mksrf_fmax(veg,landmask,                               &
                               nlon,nlat,lon,longxy,lat,latixy,edge,     &
                               fileo,ncid)


  implicit none
  include 'netcdf.inc'

! ---------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer  :: veg(nlon,nlat)	          !lsm veg array  (input)
  real(r8) :: landmask(nlon,nlat)	  !landmask array  (input)
  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)  
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid

  integer :: nlon,nlat,ncid              !number lats/lons, nc fileo id 
  character(len=80) :: fileo             !output filenames


! ------------------------------------------------------------------
! Define local variables
! ------------------------------------------------------------------

  real(r8) :: fmax(nlon,nlat)       !maximum fractiona sat area
  real(r8) :: global_avg            !global value of fmax

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id

  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: gfmax_id                     ! array id
  integer :: fmax_id                      !fmax id
  integer :: landmask_id                  !landmask id

  integer :: i,j                          !indices
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: status                       !status

  character(len=256) :: name,unit         !netCDF attributes


! ------------------------------------------------------------
! Create skeleton netcdf
! --------------------------------------------------------------

! Define global attributes

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'fmax_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat, dimlat_id)

! Define input file independent variables 

  name = 'lon'
  unit = 'degrees east'
  dim1_id(1) = dimlon_id
  call wrap_def_var (ncid,'lon', nf_float, 1, dim1_id, lon_id)
  call wrap_put_att_text (ncid, lon_id, 'long_name', name)
  call wrap_put_att_text (ncid, lon_id, 'units'    , unit)

  name = 'lat'
  unit = 'degrees north'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'lat', nf_float, 1, dim1_id, lat_id)
  call wrap_put_att_text (ncid, lat_id, 'long_name', name)
  call wrap_put_att_text (ncid, lat_id, 'units'    , unit)

  name = 'longitude-2d'
  unit = 'degrees east'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)

  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)

  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)

  name = 'western edge of surface grid'
  
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)

! Define input file specific variables

  name = 'global average FMAX'
  unit = 'fraction [0-1]'
  call wrap_def_var (ncid,'global_avg', nf_float, 0, 0, gfmax_id)
  call wrap_put_att_text (ncid, gfmax_id, 'long_name', name)
  call wrap_put_att_text (ncid, gfmax_id, 'units'    , unit)

  name = 'Maximum Fractional Saturated Area'
  unit = 'fraction [0-1]'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'FMAX' ,nf_float, 2, dim2_id, fmax_id)
  call wrap_put_att_text (ncid, fmax_id, 'long_name', name)
  call wrap_put_att_text (ncid, fmax_id, 'units'    , unit)

  name = 'land mask'
  unit = 'unitless'
  call wrap_def_var (ncid ,'LANDMASK' ,nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

! End of definition

  status = nf_enddef(ncid)

! -----------------------------------------------------------------------
! Create fmax from LSM vegtypes (if user so desires)
! currently setting number to global value 
! global value computed as global weighted average value over
! land from default file: mksrf_fmax.070406.nc
! 0.370548  mk_deeptime_frommod.ncl 
! -----------------------------------------------------------------------


  fmax = -999.99  
  do j = 1,nlon
  do i = 1,nlat
   if(landmask(j,i) == 1.) fmax(j,i) =  0.370548 
  enddo
  enddo
  global_avg =  0.370548
  print *, 'FMAX sample at point (20,40)',fmax(20,40)


! --------------------------------------------------------------------------
! Write variables
! --------------------------------------------------------------------------

  call wrap_put_var_realx (ncid, lon_id       , lon)
  call wrap_put_var_realx (ncid, lat_id       , lat)
  call wrap_put_var_realx (ncid, longxy_id    , longxy)
  call wrap_put_var_realx (ncid, latixy_id    , latixy)
  call wrap_put_var_realx (ncid, edgen_id     , edge(1))
  call wrap_put_var_realx (ncid, edgee_id     , edge(2))
  call wrap_put_var_realx (ncid, edges_id     , edge(3))
  call wrap_put_var_realx (ncid, edgew_id     , edge(4))
  call wrap_put_var_realx (ncid, landmask_id  , landmask)
  call wrap_put_var_realx (ncid, fmax_id     , fmax)
  call wrap_put_var_realx (ncid, gfmax_id     ,global_avg)

! Close output file

  call wrap_close(ncid)

end subroutine create_mksrf_fmax


!===============================================================================
subroutine create_mksrf_topo(landmask,topoin,                               &
                             nlon,nlat,lon,longxy,lat,latixy,edge,          &
                             fileo,ncid)


  implicit none
  include 'netcdf.inc'

! ---------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer  :: topoin(nlon,nlat)	          !topo array  (input)
  real(r8) :: landmask(nlon,nlat)	  !landmask array  (input)
  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)  
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid

  integer :: nlon,nlat,ncid              !number lats/lons, nc fileo id 
  character(len=80) :: fileo             !output filenames


! ------------------------------------------------------------------
! Define local variables
! ------------------------------------------------------------------

  real(r8) :: topo_bedrock(nlon,nlat)       !topographical heights (m)
  real(r8) :: topo_ice(nlon,nlat)           !topographical heights + ice sheet hgt (m)
  integer  :: numlon(nlat)                  !number of lons per lat

  integer :: dimlsmlon_id                    !netCDF dimension id
  integer :: dimlsmlat_id                    !netCDF dimension id

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id

  integer :: numlon_id                    !1d numlon array id
  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: topobed_id                   ! array id
  integer :: topoice_id                   !fmax id
  integer :: landmask_id                  !landmask id

  integer :: i,j                          !indices
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: status                       !status

  character(len=256) :: name,unit         !netCDF attributes


! ------------------------------------------------------------
! Create skeleton netcdf
! --------------------------------------------------------------

! Define global attributes

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'topoice_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat, dimlat_id)

! Define input file independent variables 

  name = 'lon'
  unit = 'degrees east'
  dim1_id(1) = dimlon_id
  call wrap_def_var (ncid,'lon', nf_float, 1, dim1_id, lon_id)
  call wrap_put_att_text (ncid, lon_id, 'long_name', name)
  call wrap_put_att_text (ncid, lon_id, 'units'    , unit)

  name = 'lat'
  unit = 'degrees north'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'lat', nf_float, 1, dim1_id, lat_id)
  call wrap_put_att_text (ncid, lat_id, 'long_name', name)
  call wrap_put_att_text (ncid, lat_id, 'units'    , unit)

  name = 'longitude-2d'
  unit = 'degrees east'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)

  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)

  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)

  name = 'western edge of surface grid'
  
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)

! Define input file specific variables

  name = 'number of grid cells at each latitude'
  unit = 'unitless'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'NUMLON', nf_float, 1, dim1_id, numlon_id)
  call wrap_put_att_text (ncid, numlon_id, 'long_name', name)
  call wrap_put_att_text (ncid, numlon_id, 'units'    , unit)

  name = 'topographical height'
  unit = 'm'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid,'TOPO', nf_float, 2, dim2_id, topobed_id)
  !! call wrap_def_var (ncid,'TOPO_BEDROCK', nf_float, 2, dim2_id, topobed_id)
  call wrap_put_att_text (ncid, topobed_id, 'long_name', name)
  call wrap_put_att_text (ncid, topobed_id, 'units'    , unit)

  name = 'topographical height + ice'
  unit = 'm'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'TOPO_ICE' ,nf_float, 2, dim2_id, topoice_id)
  call wrap_put_att_text (ncid, topoice_id, 'long_name', name)
  call wrap_put_att_text (ncid, topoice_id, 'units'    , unit)

! End of definition

  status = nf_enddef(ncid)

! -----------------------------------------------------------------------
! read topo from topo_depth file
! add code to specify ice sheet if so desired
! if you wish to read in a seperaete ice sheet topofile, 
! you will need to modify the code
! -----------------------------------------------------------------------

  do i = 1,nlat
   numlon(i) =  nlon
  end do 
 print *, numlon

 topo_ice = topoin
 topo_bedrock = topoin
 
! elevation values are for land only, replace bathymetry values with 0.
  do j = 1,nlon
  do i = 1,nlat
   if(landmask(j,i) == 0.) topo_ice(j,i) =  0. 
   if(landmask(j,i) == 0.) topo_bedrock(j,i) =  0. 
  end do
  end do

! ---- add ice mods here
!

! --------------------------------------------------------------------------
! Write variables
! --------------------------------------------------------------------------

  call wrap_put_var_int (ncid, numlon_id    , numlon)
  call wrap_put_var_realx (ncid, lon_id        , lon)
  call wrap_put_var_realx (ncid, lat_id        , lat)
  call wrap_put_var_realx (ncid, longxy_id    , longxy)
  call wrap_put_var_realx (ncid, latixy_id    , latixy)
  call wrap_put_var_realx (ncid, edgen_id     , edge(1))
  call wrap_put_var_realx (ncid, edgee_id     , edge(2))
  call wrap_put_var_realx (ncid, edges_id     , edge(3))
  call wrap_put_var_realx (ncid, edgew_id     , edge(4))
  call wrap_put_var_realx (ncid, topobed_id, topo_bedrock)
  call wrap_put_var_realx (ncid, topoice_id, topo_ice)

! Close output file

  call wrap_close(ncid)

! ---------------------------------------------------------------------

end subroutine create_mksrf_topo

!===============================================================================
subroutine create_mksrf_vocef(veg,landmask,                              &
                               nlon,nlat,lon,longxy,lat,latixy,edge,     &
                               fileo,ncid)

  implicit none
  include 'netcdf.inc'

! ---------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer  :: veg(nlon,nlat)	          !lsm veg array  (input)
  real(r8) :: landmask(nlon,nlat)	  !landmask array  (input)
  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)  
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid

  integer :: nlon,nlat,ncid              !number lats/lons, nc fileo id 
  character(len=80) :: fileo             !output filenames


! ------------------------------------------------------------------
! Define local variables
! ------------------------------------------------------------------
					 
  real(r8) :: ef_btr(nlon,nlat)          !broadleaf tree emission factor
  real(r8) :: ef_fet(nlon,nlat)          !Fineleaf evergreen tree tree emission factor
  real(r8) :: ef_fdt(nlon,nlat)          !Fineleaf deciduous tree emission factor
  real(r8) :: ef_shr(nlon,nlat)          !shrub emission factor
  real(r8) :: ef_grs(nlon,nlat)          !grass, non-vascular plants and other ground cover emission factor
  real(r8) :: ef_crp(nlon,nlat)          !crop emission factor

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id

  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: ef_btr_id                    !ef_btr id
  integer :: ef_fet_id                    !ef_fet id
  integer :: ef_fdt_id                    !ef_fdt id
  integer :: ef_shr_id                    !ef_shr id
  integer :: ef_grs_id                    !ef_grs id
  integer :: ef_crp_id                    !ef_crp id
  integer :: landmask_id                  !landmask id

  integer :: i,j                          !indices
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: status                       !status

  character(len=256) :: name,unit         !netCDF attributes


! ------------------------------------------------------------
! Create skeleton netcdf
! --------------------------------------------------------------

! Define global attributes

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'vocef_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat, dimlat_id)

! Define input file independent variables 

  name = 'lon'
  unit = 'degrees east'
  dim1_id(1) = dimlon_id
  call wrap_def_var (ncid,'lon', nf_float, 1, dim1_id, lon_id)
  call wrap_put_att_text (ncid, lon_id, 'long_name', name)
  call wrap_put_att_text (ncid, lon_id, 'units'    , unit)

  name = 'lat'
  unit = 'degrees north'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'lat', nf_float, 1, dim1_id, lat_id)
  call wrap_put_att_text (ncid, lat_id, 'long_name', name)
  call wrap_put_att_text (ncid, lat_id, 'units'    , unit)

  name = 'longitude-2d'
  unit = 'degrees east'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)

  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)

  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)

  name = 'western edge of surface grid'
  
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)

! Define input file specific variables

  name = 'broadleaf tree emission factor'
  unit = 'micrograms isoprene m-2 h-1'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'ef_btr' ,nf_float, 2, dim2_id, ef_btr_id)
  call wrap_put_att_text (ncid, ef_btr_id, 'long_name', name)
  call wrap_put_att_text (ncid, ef_btr_id, 'units'    , unit)

  name = 'Fineleaf evergreen tree tree emission factor'
  unit = 'micrograms isoprene m-2 h-1'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'ef_fet' ,nf_float, 2, dim2_id, ef_fet_id)
  call wrap_put_att_text (ncid, ef_fet_id, 'long_name', name)
  call wrap_put_att_text (ncid, ef_fet_id, 'units'    , unit)

  name = 'Fineleaf deciduous tree emission factor'
  unit = 'micrograms isoprene m-2 h-1'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'ef_fdt' ,nf_float, 2, dim2_id, ef_fdt_id)
  call wrap_put_att_text (ncid, ef_fdt_id, 'long_name', name)
  call wrap_put_att_text (ncid, ef_fdt_id, 'units'    , unit)

  name = 'shrub tree emission factor'
  unit = 'micrograms isoprene m-2 h-1'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'ef_shr' ,nf_float, 2, dim2_id, ef_shr_id)
  call wrap_put_att_text (ncid, ef_shr_id, 'long_name', name)
  call wrap_put_att_text (ncid, ef_shr_id, 'units'    , unit)

  name = 'grass, non-vascular plants and other ground cover emission factor'
  unit = 'micrograms isoprene m-2 h-1'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'ef_grs' ,nf_float, 2, dim2_id, ef_grs_id)
  call wrap_put_att_text (ncid, ef_grs_id, 'long_name', name)
  call wrap_put_att_text (ncid, ef_grs_id, 'units'    , unit)

  name = 'crop emission factor'
  unit = 'micrograms isoprene m-2 h-1'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'ef_crp' ,nf_float, 2, dim2_id, ef_crp_id)
  call wrap_put_att_text (ncid, ef_crp_id, 'long_name', name)
  call wrap_put_att_text (ncid, ef_crp_id, 'units'    , unit)

  name = 'land mask'
  unit = 'unitless'
  call wrap_def_var (ncid ,'LANDMASK' ,nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

! End of definition

  status = nf_enddef(ncid)

! -----------------------------------------------------------------------
! set all emission factors to zero, i.e. this will set isoprene emissions 
! due to vegetation in model to zero. 
!
! don't add factors unless you are using atm chemistry and know what you
! are doing. thought needs to be put into how to define emissions factors for
! paleo. for present day reference, see default modern file... a starting
! place would be to develop an equivalence between modern definitions 
! and what you expect from your paleo world. 
! -----------------------------------------------------------------------

  ef_btr = 0.
  ef_fet = 0.
  ef_fdt = 0.
  ef_shr = 0.
  ef_grs = 0.
  ef_crp = 0.

! --------------------------------------------------------------------------
! Write variables
! --------------------------------------------------------------------------

  call wrap_put_var_realx (ncid, lon_id       , lon)
  call wrap_put_var_realx (ncid, lat_id       , lat)
  call wrap_put_var_realx (ncid, longxy_id    , longxy)
  call wrap_put_var_realx (ncid, latixy_id    , latixy)
  call wrap_put_var_realx (ncid, edgen_id     , edge(1))
  call wrap_put_var_realx (ncid, edgee_id     , edge(2))
  call wrap_put_var_realx (ncid, edges_id     , edge(3))
  call wrap_put_var_realx (ncid, edgew_id     , edge(4))
  call wrap_put_var_realx (ncid, landmask_id  , landmask)
  call wrap_put_var_realx (ncid, ef_btr_id, ef_btr)
  call wrap_put_var_realx (ncid, ef_fet_id, ef_fet)
  call wrap_put_var_realx (ncid, ef_fdt_id, ef_fdt)
  call wrap_put_var_realx (ncid, ef_shr_id, ef_shr)
  call wrap_put_var_realx (ncid, ef_grs_id, ef_grs)
  call wrap_put_var_realx (ncid, ef_crp_id, ef_crp)

! Close output file

  call wrap_close(ncid)


end subroutine create_mksrf_vocef



!================================================================================
!================================================================================

subroutine wrap_create (path, cmode, ncid)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  character(len=*) path
  integer cmode, ncid, ret
  ret = nf_create (path, cmode, ncid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_create

!===============================================================================

subroutine wrap_def_dim (nfid, dimname, len, dimid)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: nfid, len, dimid
  character(len=*) :: dimname
  integer ret
  ret = nf_def_dim (nfid, dimname, len, dimid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_def_dim

!===============================================================================

subroutine wrap_def_var (nfid, name, xtype, nvdims, vdims, varid)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: nfid, xtype, nvdims, varid
  integer :: vdims(nvdims)
  character(len=*) :: name
  integer ret
  ret = nf_def_var (nfid, name, xtype, nvdims, vdims, varid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_def_var

!===============================================================================

subroutine wrap_put_att_text (nfid, varid, attname, atttext)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: nfid, varid
  character(len=*) :: attname, atttext
  integer :: ret, siz
  siz = len_trim(atttext)
  ret = nf_put_att_text (nfid, varid, attname, siz, atttext)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_put_att_text

!===============================================================================

subroutine wrap_put_var_realx (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: nfid, varid
  real(r8) :: arr(*)
  integer :: ret
#ifdef CRAY
  ret = nf_put_var_real (nfid, varid, arr)
#else
  ret = nf_put_var_double (nfid, varid, arr)
#endif
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_put_var_realx

!===============================================================================

subroutine wrap_put_var_int (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: nfid, varid
  integer :: arr(*)
  integer :: ret
  ret = nf_put_var_int (nfid, varid, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_put_var_int

!===============================================================================

subroutine wrap_put_vara_realx (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: nfid, varid
  integer :: start(*), count(*)
  real(r8) arr(*)
  integer ret
#ifdef CRAY
  ret = nf_put_vara_real (nfid, varid, start, count, arr)
#else
  ret = nf_put_vara_double (nfid, varid, start, count, arr)
#endif
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_put_vara_realx

  
!===============================================================================

subroutine wrap_close (ncid)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: ncid
  integer :: ret
  ret = nf_close (ncid)
  if (ret.ne.NF_NOERR) then
     write(6,*)'WRAP_CLOSE: nf_close failed for id ',ncid
     call handle_error (ret)
  end if
end subroutine wrap_close

!===============================================================================

subroutine handle_error(ret)
  implicit none
  include 'netcdf.inc'
  integer :: ret
  if (ret .ne. nf_noerr) then
     write(6,*) 'NCDERR: ERROR: ',nf_strerror(ret)
     call abort
  endif
end subroutine handle_error

!==============================================================================

subroutine wrap_inq_varid (nfid, varname, varid)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) varname

  integer ret

  ret = nf_inq_varid (nfid, varname, varid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_inq_varid

!=============================================================================

subroutine wrap_get_var8 (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  real*8 arr(*)

  integer ret

  ret = nf_get_var_double (nfid, varid, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_get_var8

!==============================================================================

   subroutine wrap_get_var_int (nfid, varid, arr)
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(out):: arr(*)

   integer ret      ! NetCDF return code

   ret = nf_get_var_int (nfid, varid, arr)
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_GET_VAR_INT: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_var_int

!=============================================================================

subroutine endrun
  implicit none
  include 'netcdf.inc'

  call abort
  stop 999
end subroutine endrun


!===============================================================================


