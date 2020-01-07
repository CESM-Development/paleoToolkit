program cal_adjust_PMIP
! Calculates calendar (month-length) adjustments of data in a CMIP/PMIP-formatted netCDF file.
! Creates a new netCDF file by copying dimension variables and global attributes from the input file.
! This version supports 3-D (e.g. longitude, latitude, time) and 4-D (e.g. longitude, latitude, level, time)
! long-term mean (AClim, Oclim), monthly (Amon, Omon, OImon) or daily input files.

! The program requires the modules: calendar_effects_subs.f90, CMIP_netCDF_subs.f90, GISS_orbpar_subs.f90,
! GISS_srevents_subs.f90, month_length_subs.f90 and pseudo_daily_interp_subs.f90
! The program must be compiled with local netCDF support, and it will use OpenMP if available

! An info .csv file (with an appropriate header) containing the following information is read:
! activity      :: (string) PMIP3 or PMIP4
! variable      :: (string) CMIP/PMIP variable name (e.g. "tas", "pr")
! time_freq     :: (string) CMIP/PMIP output time-frequency type (e.g. "Amon", "Aclim",...)
! model         :: (string) CMIP/PMIP model name
! experiment    :: (string) CMIP/PMIP experiment name (e.g. "midHolocene")
! ensemble      :: (string) CMIP/PMIP ensemble member code (e.g. "r1i1p1")
! grid_label    :: (string) CMIP6/PMIP4 grid label (blank if CMIP5/PMIP3)
! begdate       :: (string) beginning date of simulation as YYYYMM or YYYYMMDD string
! enddate       :: (string) ending date of simulation as YYYYMM or YYYYMMDD string
! suffix        :: (string) input filename suffix (usually blank, or "clim" for Aclim-type files)
! adj_name      :: (string) output filename suffix for adjusted netCDF file (e.g. "_cal_adj")
! calendar_type :: (string) calendar type (e.g. "noleap", "proleptic_gregorian", etc.)
! begageBP      :: (integer) beginning simulation age (year BP) (e.g. 21000 (= 21 ka)
! endageBP      :: (integer) ending simulation age (year BP) 
! agestep       :: (integer) interval between age calculations
! begyrCE       :: (integer) beginning calendar year of simulation for multi-year simulations at each age
! nsimyrs       :: (integer) number of calendar years
! source_path   :: (string) path to (input) source files (enclosed in quotation marks)
! adjusted_path :: (string) path to (output) adjusted files (enclosed in quotation marks)

! (The above information could also be gotten by parsing the netCDF file name, and reading the calendar attribute.)

! Author: Patrick J. Bartlein, Univ. of Oregon (bartlein@uoregon.edu), with contributions by S.L. Shafer (sshafer@usgs.gov)
!
! Version: 1.0d  (Generalized to handle 4-D files)
! Last update: 2019-08-05

use calendar_effects_subs
use pseudo_daily_interp_subs
use month_length_subs
use CMIP_netCDF_subs
use omp_lib

implicit none

! past ages are negative, e.g. 21 ka = 21,000 cal yr BP = -21000, and 1950 CE = 0 cal yr BP = 0
! simulation age-related variables (controls of orbital parameters)
integer(4)              :: begageBP             ! beginning year (BP) (negative, e.g. 10 ka = -10000 BP)
integer(4)              :: endageBP             ! ending year (BP)
integer(4)              :: agestep              ! age step size
integer(4)              :: nages                ! number of simulation ages
integer(4), allocatable :: iageBP(:)            ! (ny) year BP 1950 (negative, e.g. 1900 CE = -50 years BP 1950)

! simulation year-related variables (controls of vernal equinox day and leap-year status)
integer(4)              :: begyrCE              ! beginning (pseudo-) year of individual model simulation
integer(4)              :: nsimyrs              ! number of years of simulation
integer(4), allocatable :: iyearCE(:)           ! (ny) yearCE simulation year (e.g. 1850CE, 850CE, etc.)

! month-length variables
integer(4), allocatable :: imonlen_0ka(:,:)     ! (ny,nm) integer-value month lengths -- 0ka
integer(4), allocatable :: imonmid_0ka(:,:)     ! (ny,nm) integer-value mid days -- 0ka
integer(4), allocatable :: imonbeg_0ka(:,:)     ! (ny,nm) integer-value month beginning days -- 0ka
integer(4), allocatable :: imonend_0ka(:,:)     ! (ny,nm) integer-value month ending days -- 0ka
integer(4), allocatable :: imonlen(:,:)         ! (ny,nm) integer-value month lengths (paleo)
integer(4), allocatable :: imonmid(:,:)         ! (ny,nm) integer-value mid days (paleo)
integer(4), allocatable :: imonbeg(:,:)         ! (ny,nm) integer-value month beginning (paleo)
integer(4), allocatable :: imonend(:,:)         ! (ny,nm) integer-value month ending days (paleo)
real(8), allocatable    :: rmonlen(:,:)         ! (ny,nm) real-value month lengths (paleo)
real(8), allocatable    :: rmonmid(:,:)         ! (ny,nm) real-value mid days (paleo)
real(8), allocatable    :: rmonbeg(:,:)         ! (ny,nm) real-value month beginning (paleo)
real(8), allocatable    :: rmonend(:,:)         ! (ny,nm) real-value month ending days (paleo)
real(8), allocatable    :: VE_day(:)            ! (ny) vernal equinox day in simulation year
real(8), allocatable    :: SS_day(:)            ! (ny) (northern) summer solstice day in simulation year
integer(4), allocatable :: ndays(:)             ! (ny) number of days in year

integer(4), allocatable :: imonlen_0ka_ts(:)    ! (nt) integer-value month lengths at present as time series
real(8), allocatable    :: rmonmid_ts(:)        ! (nt) real-value paleo month mid days as time series
real(8), allocatable    :: rmonbeg_ts(:)        ! (nt) real-value paleo month beginning days as time series
real(8), allocatable    :: rmonend_ts(:)        ! (nt) real-value paleo month ending as time series
integer(4), allocatable :: imonmid_ts(:)        ! (nt) integer-value paleo month mid days as time series
integer(4), allocatable :: imonbeg_ts(:)        ! (nt) integer-value paleo month beginning days as time series
integer(4), allocatable :: imonend_ts(:)        ! (nt) integer-value paleo month ending as time series
integer(4), allocatable :: ndays_ts(:)          ! (nt) integer-value times series of paleo year lengths
character(256)          :: time_comment         ! source of new monthly time values
real(8), allocatable    :: mon_time(:)          ! (nt) new monthly time values for daily-input files
real(8), allocatable    :: mon_time_bnds(:,:)   ! (2,nt) new monthly time-bounds values for daily input files

! other names
character(32)           :: calendar_type        ! calendar type

! components of file names
character(8)            :: activity             ! PMIP version (PMIP3 or PMIP4)
character(64)           :: variable             ! variable name
character(8)            :: time_freq            ! type of CMIP/PMIP time frequency (e.g. Aclim, Amon, day, etc.)
character(8)            :: time_freq_output     ! time_freq output label
character(64)           :: model                ! model name
character(64)           :: experiment           ! experiment name
character(16)           :: ensemble             ! ensemble designator
character(64)           :: grid_label           ! grid type (PMIP4 only)
character(8)            :: begdate, enddate     ! string beginning and ending dates of simulation
character(32)           :: suffix               ! file name suffix (e.g. "-clim")
character(32)           :: adj_name             ! adjustment name (e.g. "_adj")

! data
integer(4)              :: ny, nt               ! number of years and total nmber of months (nt = ny*nm)
integer(4)              :: ndtot,ndtot_0ka,ndyr ! total number of days
real(4), allocatable    :: var3d_in(:,:,:)      ! 3-D (e.g. nlon,nlat,ndtot) input data 
real(4), allocatable    :: var3d_out(:,:,:)     ! 3-D (e.g. nlon,nlat,nt) output adjusted data 
real(8), allocatable    :: xdh(:,:)             ! 3-D (e.g. ivar_dimlen(1),ndtot) pseudo- or actual daily data
real(8), allocatable    :: var_adj(:,:)       ! 3-D (e.g. ivar_dimlen(1),nt) adjusted data 
real(8), allocatable    :: var4d_in(:,:,:,:)    ! 4-D (e.g. nlon,nlat,nlev,nt) input data
real(4), allocatable    :: var4d_out(:,:,:,:)   ! 4-D (e.g. nlon,nlat,nlev,nt) output adjusted data
real(8)                 :: vfill                ! fill value

! smoothing parameters for multi-year pseudo-daily interpolation
integer(4)              :: nw_tmp=21, nsw_tmp=20    ! smoothing parameters
logical                 :: smooth=.true., restore=.true.
logical                 :: no_negatives = .false.   ! restrict pseudo-daily interpolated values to positive values?

integer(4)              :: n,m,j,k,l,i          ! indices
integer(4)              :: max_threads          ! if OpenMP enabled
integer(4)              :: iostatus             ! IOSTAT value

! file paths and names
character(2048)         :: source_path, ncfile_in, adjusted_path, ncfile_out, nc_fname, infopath
character(64)           :: infofile
character(1)            :: csvheader            ! info .csv file header

! timers
real(4)                 :: overall_secs, monlen_secs, define_secs, read_secs, allocate_secs
real(4)                 :: openmp_secs, write_secs, total_secs

! progress bar
integer(4)              :: ntotalpts, numberdone, nprogress, nextra
real(4)                 :: onepercent, nextonepercent

! if OpenMP enabled
max_threads = omp_get_max_threads()
write (*,'("OMP max_threads: ",i4)') max_threads
max_threads = max_threads !- 2 ! subtract 2 to free up threads for other computational processes
call omp_set_num_threads(max_threads)

! info files
infopath = "/project/cseg/mvr/angular/PaleoCalAdjust/data/info_files/" ! Windows path
!infopath = "/Users/bartlein/Projects/Calendar/PaleoCalAdjust/data/info_files/"  ! Mac path
infofile = "cal_adj_info.csv"

! open the info file, and loop over specified calendar tables
open (3,file=trim(infopath)//trim(infofile))
write (*,'(a)') trim(infopath)//trim(infofile)
read (3,'(a)') csvheader

! main loop (over individual model output files)
iostatus = 1
overall_secs = secnds(0.0)
do
    total_secs = secnds(0.0); monlen_secs = secnds(0.0)
    ! Step 1:  Read info file, construct variable and file names, allocate month-length arrays, etc.
    
    ! read a line from the info file, and construct netCDF file names
    write (*,'(125("*"))')
    suffix = ""
    read (3,*,iostat=iostatus) activity, variable, time_freq, model, experiment, ensemble, grid_label, &
        begdate, enddate, suffix, adj_name, calendar_type, begageBP, endageBP, agestep, begyrCE, nsimyrs, &
        source_path, adjusted_path
    if (iostatus.lt.0) then
        write (*,'(a)') "*** Done ***"
        exit
    end if
    write (*,'(a)') trim(activity)
    write (*,'(a)') trim(source_path)

    varinname = trim(variable); varoutname = varinname
    time_freq_output = trim(time_freq)
    if (trim(time_freq_output) .eq. 'day') time_freq_output = "Amon2"
    if (suffix .eq. 'clim') suffix = "-"//trim(suffix)

    select case (trim(activity))
    case ('PMIP3', 'pmip3')
        ncfile_in = trim(variable)//"_"//trim(time_freq)//"_"//trim(model)//"_"//trim(experiment)//"_"// &
            trim(ensemble)//"_"//trim(begdate)//"-"//trim(enddate)//trim(suffix)//".nc"
        ncfile_out = trim(variable)//"_"//trim(time_freq_output)//"_"//trim(model)//"_"//trim(experiment)//"_"// &
            trim(ensemble)//"_"//trim(begdate)//"-"//trim(enddate)//trim(suffix)//trim(adj_name)//".nc"
    case ('PMIP4', 'pmip4')
        ncfile_in = trim(variable)//"_"//trim(time_freq)//"_"//trim(model)//"_"//trim(experiment)//"_"// &
            trim(ensemble)//"_"//trim(grid_label)//"_"//trim(begdate)//"-"//trim(enddate)//trim(suffix)//".nc"
        ncfile_out = trim(variable)//"_"//trim(time_freq_output)//"_"//trim(model)//"_"//trim(experiment)//"_"// &
            trim(ensemble)//"_"//trim(grid_label)//"_"//trim(begdate)//"-"//trim(enddate)//trim(suffix)//trim(adj_name)//".nc"

     case ('CAM', 'cam')
        ncfile_in = trim(experiment)//".cam.h0."//trim(variable)//"."//trim(begdate)//"-"//trim(enddate)//".nc"
        ncfile_out = trim(experiment)//".cam.h0."//trim(variable)//"."//trim(begdate)//"-"//trim(enddate)//trim(adj_name)//".nc"

     case ('CAM2', 'cam2')
        ncfile_in = trim(experiment)//".cam2.h0."//trim(variable)//"."//trim(begdate)//"-"//trim(enddate)//".nc"
        ncfile_out = trim(experiment)//".cam2.h0."//trim(variable)//"."//trim(begdate)//"-"//trim(enddate)//trim(adj_name)//".nc"

     case ('CICE', 'cice')
        ncfile_in = trim(experiment)//".cice.h."//trim(variable)//"."//trim(begdate)//"-"//trim(enddate)//".nc"
        ncfile_out = trim(experiment)//".cice.h."//trim(variable)//"."//trim(begdate)//"-"//trim(enddate)//trim(adj_name)//".nc"

     case ('CSIM', 'csim')
        ncfile_in = trim(experiment)//".csim.h."//trim(variable)//"."//trim(begdate)//"-"//trim(enddate)//".nc"
        ncfile_out = trim(experiment)//".csim.h."//trim(variable)//"."//trim(begdate)//"-"//trim(enddate)//trim(adj_name)//".nc"

     case ('POP', 'pop')
        ncfile_in = trim(experiment)//".pop.h."//trim(variable)//"."//trim(begdate)//"-"//trim(enddate)//".nc"
        ncfile_out = trim(experiment)//".pop.h."//trim(variable)//"."//trim(begdate)//"-"//trim(enddate)//trim(adj_name)//".nc"

     case default
        stop "activity type"
    end select

    write (*,'(" ncfile_in: ",a)') trim(ncfile_in)
    write (*,'("ncfile_out: ",a)') trim(ncfile_out)
    write (*,'(a, 1x, a, 1x, 5i7)') trim(variable),trim(calendar_type), begageBP, endageBP, agestep, begyrCE, nsimyrs

    ! set no_negatives flag
    no_negatives = .false.
    select case (trim(variable))
    case ('pr','clt','sic')
        no_negatives = .true.
    case default
        continue
    end select
    write (*,'("no_negatives: ", l1)') no_negatives

    ! allocate month-length arrays
    nages = (endageBP - begageBP)/agestep + 1
    ny = nages * nsimyrs
    nt = ny * nm
    write (*,'("nages, nsimyrs, nt: ",3i8)') nages, nsimyrs, nt
    allocate (iageBP(ny), iyearCE(ny))
    allocate (imonlen_0ka(ny,nm),imonmid_0ka(ny,nm),imonbeg_0ka(ny,nm),imonend_0ka(ny,nm))
    allocate (imonlen(ny,nm), imonmid(ny,nm), imonbeg(ny,nm), imonend(ny,nm))
    allocate (rmonlen(ny,nm), rmonmid(ny,nm), rmonbeg(ny,nm), rmonend(ny,nm))
    allocate (VE_day(ny), SS_day(ny), ndays(ny), imonlen_0ka_ts(nt),ndays_ts(nt))
    allocate (imonmid_ts(nt), imonbeg_ts(nt),imonend_ts(nt), rmonmid_ts(nt), rmonbeg_ts(nt), rmonend_ts(nt))
    allocate (mon_time(nt), mon_time_bnds(2,nt))

    ! Step 2:  get month lengths

    ! 0 ka month lengths (used in pseudo-daily interpolation)
    if (trim(time_freq) .ne. 'day') then
        write (*,'(a)') "0 ka month lengths for pseudo-daily interpolation..."
        call get_month_lengths(calendar_type, 0, 0, agestep, nages, begyrCE, nsimyrs, &
            iageBP, iyearCE, imonlen_0ka, imonmid_0ka, imonbeg_0ka, imonend_0ka, rmonlen, rmonmid, rmonbeg, rmonend, &
                VE_day, SS_day, ndays)
    end if

    ! paleo month lengths
    write (*,'(a)') "Paleo month-lengths for aggregation of daily data..."
    call get_month_lengths(calendar_type, begageBP, endageBP, agestep, nages, begyrCE, nsimyrs, &
        iageBP, iyearCE, imonlen, imonmid, imonbeg, imonend, rmonlen, rmonmid, rmonbeg, rmonend, VE_day, SS_day, ndays)

    ! reshape month lengths into time series, and get cumulative number of days in years
    ndtot_0ka = 0; ndyr = 0
    do n=1,ny
        ndtot_0ka=ndtot_0ka + ndays(n)
        ndyr = ndyr + ndays(n)
        do m=1,nm
            i = (n-1)*nm + m
            imonlen_0ka_ts(i) = imonlen_0ka(n,m)
            rmonmid_ts(i) = rmonmid(n,m); rmonbeg_ts(i) = rmonbeg(n,m); rmonend_ts(i) = rmonend(n,m)
            imonmid_ts(i) = imonmid(n,m); imonbeg_ts(i) = imonbeg(n,m); imonend_ts(i) = imonend(n,m)
            ndays_ts(i) = ndyr
        end do
    end do

    ! total number of days in simulation
    ndtot = 0
    do n=1,ny
        ndtot = ndtot + ndays(n)
        !write (10,'(2i6,f12.6,2i8)') iageBP(n),iyearCE(n),VE_day(n),ndays(n),ndtot
    end do
    write (*,'("ny, nt, ndtot: ",4i8)') ny, nt, ndtot, ndtot_0ka
    write (*,'(a,f7.2)') "Month-length time: ", secnds(monlen_secs)
    
    ! Step 3:  Open input and create output netCDF files 

    ! input netCDF file
    define_secs = secnds(0.0)
    write (*,'(a)') trim(source_path)
    nc_fname = trim(source_path)//trim(ncfile_in)
    print '(" nc_fname (in) = ",a)', trim(nc_fname)
    call check( nf90_open(nc_fname, nf90_nowrite, ncid_in) )
    if (nc_print) print '("  ncid_in = ",i8)', ncid_in

    ! output netCDF file
    nc_fname = trim(adjusted_path)//trim(ncfile_out)
    print '(" nc_fname (out) = ",a)', trim(nc_fname)
    call check( nf90_create(nc_fname, 0, ncid_out) )
    if (nc_print) print '("  ncid_out = ",i8)', ncid_out
    
    ! Step 4:  Redefine the time-coordinate variable

    ! redefine the time coordinate
    if (trim(time_freq) .eq. 'day') then
        time_comment = "paleo monthly mid, beginning and ending dates from original daily values"
        call new_time_day(ncid_in, ny, nm, nt, ndtot, &
            imonmid_ts, imonbeg_ts, imonend_ts, ndays_ts, mon_time, mon_time_bnds)
    else
        time_comment = "paleo monthly values from get_month_lengths()"
        call new_time_month(calendar_type, ncid_in, ny, nm, nt, &
            rmonmid_ts, rmonbeg_ts, rmonend_ts, ndays_ts, mon_time, mon_time_bnds)
    end if
    
    ! Step 5:  Define the new netCDF file
    
    ! get number of input variable dimensions, dimension names and lengths
    call get_var_diminfo(ncid_in, varinname, invar_ndim, invar_dimid, invar_dimlen, invar_dimname)

    ! define the new netCDF file, and copy dimension variables and global attributes
    call current_time(current)
    addglattname = "paleo_calendar_adjustment"
    addglatt = trim(current)//" paleo calendar adjustment by cal_adjust_PMIP.f90"
    call copy_dims_and_glatts(ncid_in, ncid_out, addglattname, addglatt, nt, &
        mon_time, mon_time_bnds, time_comment, varid_out)
    
    ! define the output (adjusted) variable, and copy attributes
    write (*,'(a)') "Defining (new) adjusted variable..."
    addvarattname = "paleo_calendar_adjustment"
    select case (trim(time_freq))
    case ('Aclim', 'Oclim', 'OIclim', 'LIclim')
        addvaratt = "long-term mean values adjusted for the appropriate paleo calendar"
    case ('Amon', 'Omon', 'OImon', 'LImon')
        addvaratt = "monthly values adjusted for the appropriate paleo calendar"
    case ('day')
        addvaratt = "daily values aggregated to monthly values using the appropriate paleo month lengths"
    case default 
        stop "defining new variable"
    end select

    ! define the output variable
    call define_outvar(ncid_in, ncid_out, varinname, varid_out, varoutname, addvarattname, addvaratt, varid_in) 
    write (*,'(a,f7.2)') "Define time: ", secnds(define_secs)

    ! Step 6:  Get the input variable to be adjusted
    
    ! allocate variables
    allocate_secs = secnds(0.0)
    write (*,'(a)') "Allocating arrays"
    select case(invar_ndim)
    case (3)
        write (*,*) invar_dimlen(1),invar_dimlen(2),nt,invar_dimlen(1)*invar_dimlen(2)*nt
        if (trim(time_freq) .eq. 'day') then
            allocate(var3d_in(invar_dimlen(1),invar_dimlen(2),ndtot))
        else
            allocate(var3d_in(invar_dimlen(1),invar_dimlen(2),nt))
        end if
        allocate(var3d_out(invar_dimlen(1),invar_dimlen(2),nt))
    case (4)
        write (*,*) invar_dimlen(1),invar_dimlen(2),invar_dimlen(3),nt,invar_dimlen(1)*invar_dimlen(2)*invar_dimlen(3)*nt
        if (trim(time_freq) .eq. 'day') then
            allocate(var4d_in(invar_dimlen(1),invar_dimlen(2),invar_dimlen(3),ndtot))
        else
            allocate(var4d_in(invar_dimlen(1),invar_dimlen(2),invar_dimlen(3),nt))
        end if
        allocate(var4d_out(invar_dimlen(1),invar_dimlen(2),invar_dimlen(3),nt))
    case default
        stop "allocating variables"
    end select
    allocate(xdh(invar_dimlen(1),ndtot), var_adj(invar_dimlen(1),nt))
   write (*,'(a,f7.2)') "Allocate time: ", secnds(allocate_secs)
    
    ! get input data
    write (*,'(a)') "Reading input data..."
    read_secs = secnds(0.0)
    select case (invar_ndim)
    case (3)
        call check( nf90_get_var(ncid_in, varid_in, var3d_in) )
    case (4)
        call check( nf90_get_var(ncid_in, varid_in, var4d_in) )
    case default
        stop "input data"
    end select

    ! get _FillValue
    if ( nf90_get_att(ncid_in, varid_in, '_FillValue', vfill).ne.nf90_noerr) then
       vfill = -900.
    endif
    write (*,'("_FillValue:", g14.6)') vfill
    write (*,'(a,f7.2)') "Read time: ", secnds(read_secs)
    
    ! Step 7:  Get calendar-adjusted values
    
    openmp_secs = 0.0
    if (trim(time_freq) .eq.'day') then
        write (*,'(a)') "Aggregating daily data ..."
    else
        write (*,'(a)') "Interpolating and aggregating ..."
    end if
    
    openmp_secs = secnds(0.0)
    select case(invar_ndim)
    case (3) 
        ntotalpts = invar_dimlen(1) * invar_dimlen(2)
        write (*,'("Number of grid points = ",i6)') ntotalpts
    case (4)
        ntotalpts = invar_dimlen(1) * invar_dimlen(2) * invar_dimlen(3)
        write (*,'("Number of grid points x number of levels = ",i6)') ntotalpts
    end select
        
    numberdone = 0; onepercent = (ntotalpts/100); nextonepercent = onepercent; nprogress = 0
    write (*,'("|",9("---------+"),"---------|")'); write (*,'(a1,$)') " "
    
    ! Note: OpenMP seems to work best if only innermost loop is parallelized
    select case (invar_ndim)
    case (3)
        do k=1,invar_dimlen(2)      
            ! if daily data, copy into xdh
            if (trim(time_freq) .eq.'day') xdh(:,:) = var3d_in(:,k,:)
        
            !$omp parallel do
            do j=1,invar_dimlen(1)
                !write (*,'(/2i5)') j,k
                if (trim(time_freq) .eq.'day') xdh(j,:) = var3d_in(j,k,:)
                ! unless the input data are daily, do pseudo-daily interpolation of the monthly input data
                if (trim(time_freq) .ne. 'day') then
                    ! interpolate
                    call mon_to_day_ts(nt, imonlen_0ka_ts, dble(var3d_in(j,k,:)), dble(vfill), &
                        no_negatives, smooth, restore, ndtot, nw_tmp, nsw_tmp, xdh(j,:))
                    ! reaggregate daily data using correct calendar
                    call day_to_mon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xdh(j,:),dble(vfill),var_adj(j,:))
                else
                    ! input data are already daily, so just reaggregate using correct calendar
                    call day_to_mon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xdh(j,:),dble(vfill),var_adj(j,:))
                end if
            
                var3d_out(j,k,:)=sngl(var_adj(j,:))

            end do
            !$omp end parallel do
          
            ! update progress bar
            numberdone = numberdone + invar_dimlen(1)
            if (numberdone .ge. nextonepercent) then
                if (nprogress .lt. 99) then
                    write (*,'(a1,$)') "="
                    nextonepercent = nextonepercent + onepercent
                    nprogress = nprogress + 1
                    !write (*,*) numberdone, nextonepercent
                end if
            end if

        end do
        
    case(4)
        do l=1,invar_dimlen(3)
            do k=1,invar_dimlen(2)      
                ! if daily data, copy into xdh
                if (trim(time_freq) .eq.'day') xdh(:,:) = var4d_in(:,k,l,:)
        
                !$omp parallel do
                do j=1,invar_dimlen(1)
                    !write (*,'(/3i5)') j,k,l
                    if (trim(time_freq) .eq.'day') xdh(j,:) = var4d_in(j,k,l,:)
                    ! unless the input data are daily, do pseudo-daily interpolation of the monthly input data
                    if (trim(time_freq) .ne. 'day') then
                        ! interpolate
                        call mon_to_day_ts(nt, imonlen_0ka_ts, dble(var4d_in(j,k,l,:)), dble(vfill), &
                            no_negatives, smooth, restore, ndtot, nw_tmp, nsw_tmp, xdh(j,:))
                        ! reaggregate daily data using correct calendar
                        call day_to_mon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xdh(j,:),dble(vfill),var_adj(j,:))
                    else
                        ! input data are already daily, so just reaggregate using correct calendar
                        call day_to_mon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xdh(j,:),dble(vfill),var_adj(j,:))
                    end if
            
                    var4d_out(j,k,l,:)=sngl(var_adj(j,:))

                end do
                !$omp end parallel do
          
                ! update progress bar
                numberdone = numberdone + invar_dimlen(1)
                if (numberdone .ge. nextonepercent) then
                    if (nprogress .lt. 99) then
                        write (*,'(a1,$)') "="
                        nextonepercent = nextonepercent + onepercent
                        nprogress = nprogress + 1
                        !write (*,*) numberdone, nextonepercent
                    end if
                end if

            end do
        end do
    case default
        stop "adjusting"
    end select
    
    ! finish progress bar if short
    if (nprogress .lt. 99) then
        do nextra=1,(99-nprogress)
            write (*,'(a1,$)') "="
        end do
    end if
    write (*,'(a)') " "
    
    write (*,'(a,f7.2)') "OpenMP time: ", secnds(openmp_secs)

    ! Step 8:  Write out the adjusted data, and close the output netCDF file.
    
    ! write out adjusted data
    write_secs = secnds(0.0)
    write (*,'(/a)') "Writing adjusted data..."
    select case(invar_ndim)
    case (3)
        call check( nf90_put_var(ncid_out, varid_out, var3d_out) )
    case (4)
        call check( nf90_put_var(ncid_out, varid_out, var4d_out) )
    case default
        stop "adjusted output"
    end select
    ! close the output file
    call check( nf90_close(ncid_out) )
    write (*,'(a,f7.2)') "write time: ", secnds(write_secs)

    deallocate (iageBP, iyearCE)
    deallocate (imonlen_0ka,imonmid_0ka,imonbeg_0ka,imonend_0ka)
    deallocate (imonlen, imonmid, imonbeg, imonend)
    deallocate (rmonlen, rmonmid, rmonbeg, rmonend)
    deallocate (VE_day, SS_day, ndays, imonlen_0ka_ts,ndays_ts)
    deallocate (imonmid_ts, imonbeg_ts,imonend_ts, rmonmid_ts, rmonbeg_ts, rmonend_ts)
    deallocate (mon_time, mon_time_bnds)
    select case (invar_ndim)
    case (3)
        deallocate (var3d_in, var3d_out)
    case (4)
        deallocate (var4d_in, var4d_out)
    case default
        stop "deallocating"
    end select
    deallocate (xdh, var_adj)

    !write (*,'(a)') " "
    write (*,'(a,f7.2)') "Total time: ", secnds(Total_secs)
    
end do

    write (*,'(a,f7.2)') "Overall time: ", secnds(overall_secs)

end program

