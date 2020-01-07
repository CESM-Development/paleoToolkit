module month_length_subs

    implicit none

    integer(4), parameter   :: nm = 12              ! number of months in the year
    integer(4), parameter   :: nd_360 = 360         ! number of days in a 360-day year
    integer(4), parameter   :: daysinmonth360 = 30  ! number of days in a month of a 360-day year
    integer(4), parameter   :: nd_365 = 365         ! number of days in a 365-day "noleaps" year
    integer(4), parameter   :: nd_366 = 366         ! number of days in a 366-day leap year

    ! other calendar-related variables
    real(8)     :: veqday_360 = 80.0d0              ! fixed vernal equinox day, 360-day year
    real(8)     :: veqday_365 = 80.5d0              ! fixed vernal equinox day, 365-day year
    real(8)     :: veqday_366 = 81.5d0              ! fixed vernal equinox day, 366-day year
    real(8)     :: ssday_360 = 170.5d0              ! fixed (northern) summer solstice day, 360-day year
    real(8)     :: ssday_365 = 173.0d0              ! fixed (northern) summer solstice day, 365-day year
    real(8)     :: ssday_366 = 173.5d0              ! fixed (northern) summer solstice day, 366-day year
    real(8)     :: tropical_year = 365.24219876d0   ! length of a tropical year (days)
    real(8)     :: progreg_year = 365.2425d0        ! length of a Gregorian year (days)

    integer(4)  :: nd_progreg                       ! number of days in a 365 or 366-day year proleptic_gregorian calendar
    real(8)     :: veqday_progreg                   ! vernal equinox day in a 365 or 366-day year proleptic_gregorian calendar
    real(8)     :: midMarch                         ! mid-March day
    real(8)     :: perihelion                       ! perihelion day (in a proleptic Gregorian calendar)
    real(8)     :: aphelion                         ! aphelion day (in a proleptic Gregorian calendar)

    ! month-length definitions
    real(8)     :: present_mon_360(nm) = 30.0d0     ! present-day month lengths in 360-day year
    real(8)     :: present_mon_noleap(nm) = &       ! present-day month lengths in 365-day (noleap) year
        (/ 31.0d0, 28.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    real(8)     :: present_mon_leap(nm) = &         ! present-day month lengths in 366-day (leap) year
        (/ 31.0d0, 29.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    real(8)     :: present_mon_365_trop(nm) = &     ! present-day month lengths in a tropical year (note Feb.)
        (/ 31.0d0, 28.24219876d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    real(8)     :: present_mon_365_progreg(nm) = &  ! present-day month lengths in a Gregorian year (note Feb.)
        (/ 31.0d0, 28.2425d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    
    real(8)     :: present_beg_360(nm) = &          ! present-day month beginning day in 360-day year
        (/  0.0d0, 30.0d0, 60.0d0, 90.0d0, 120.0d0, 150.0d0, 180.0d0, 210.0d0, 240.0d0, 270.0d0, 300.0d0, 330.0d0 /)
    real(8)     :: present_mid_360(nm) = &          ! present-day month middle day in 360-day year
        (/ 15.0d0, 45.0d0, 75.0d0, 105.0d0, 135.0d0, 165.0d0, 195.0d0, 225.0d0, 255.0d0, 285.0d0, 315.0d0, 345.0d0 /)
    real(8)     :: present_end_360(nm) = &          ! present-day month ending day in 360-day year
        (/ 30.0d0, 60.0d0, 90.0d0, 120.0d0, 150.0d0, 180.0d0, 210.0d0, 240.0d0, 270.0d0, 300.0d0, 330.0d0, 360.0d0 /)
    real(8)     :: present_beg_365(nm) = &          ! present-day month beginning day in 365-day (noleap) year
        (/  0.0d0, 31.0d0, 59.0d0, 90.0d0, 120.0d0, 151.0d0, 181.0d0, 212.0d0, 243.0d0, 273.0d0, 304.0d0, 334.0d0 /)
    real(8)     :: present_mid_365(nm) = &          ! present-day month middle day in 365-day (noleap) year
        (/ 15.5d0, 45.0d0, 74.5d0, 105.0d0, 135.5d0, 166.0d0, 196.5d0, 227.5d0, 258.0d0, 288.5d0, 319.0d0, 349.5d0 /)
    real(8)     :: present_end_365(nm) = &          ! present-day month ending day in 365-day (noleap) year
        (/ 31.0d0, 59.0d0, 90.0d0, 120.0d0, 151.0d0, 181.0d0, 212.0d0, 243.0d0, 273.0d0, 304.0d0, 334.0d0, 365.0d0 /)
    real(8)     :: present_beg_366(nm) = &          ! present-day month beginning day in 366-day (leap) year
        (/  0.0d0, 31.0d0, 60.0d0, 91.0d0, 121.0d0, 152.0d0, 182.0d0, 213.0d0, 244.0d0, 274.0d0, 305.0d0, 335.0d0 /)
    real(8)     :: present_mid_366(nm) = &          ! present-day month beginning day in 366-day (leap) year
        (/ 15.5d0, 45.5d0, 75.5d0, 106.0d0, 136.5d0, 167.0d0, 197.5d0, 228.5d0, 259.0d0, 289.5d0, 320.0d0, 350.5d0 /)
    real(8)     :: present_end_366(nm) = &          ! present-day month beginning day in 366-day (leap) year
        (/ 31.0d0, 60.0d0, 91.0d0, 121.0d0, 152.0d0, 182.0d0, 213.0d0, 244.0d0, 274.0d0, 305.0d0, 335.0d0, 366.0d0 /)
    real(8)     :: present_beg_365_trop(nm) = &     ! present-day month beginning day in a tropical year
        (/  0.0000d0, 31.0000d0, 59.2422d0, 90.2422d0, 120.2422d0, 151.2422d0, 181.2422d0, 212.2422d0, 243.2422d0, &
            273.2422d0, 304.2422d0, 334.2422d0 /)
    real(8)     :: present_mid_365_trop(nm) = &     ! present-day month middle day in a tropical year
        (/ 15.5000d0, 45.1211d0, 74.7422d0, 105.2422d0, 135.7422d0, 166.2422d0, 196.7422d0, 227.7422d0, 258.2422d0, &
            288.7422d0, 319.2422d0, 349.7422d0 /)
    real(8)     :: present_end_365_trop(nm) = &     ! present-day month ending day in a tropical year
        (/ 31.0000d0, 59.2422d0, 90.2422d0, 120.2422d0, 151.2422d0, 181.2422d0, 212.2422d0, 243.2422d0, 273.2422d0, &
            304.2422d0, 334.2422d0, 365.2422d0 /)
    real(8)     :: present_beg_365_progreg(nm) = &  ! present-day month beginning day in a Gregorian year
        (/  0.0000d0, 31.0000d0, 59.2425d0, 90.2425d0, 120.2425d0, 151.2425d0, 181.242d0, 212.2425d0, 243.2425d0, &
            273.2425d0, 304.2425d0, 334.2425d0 /)
    real(8)     :: present_mid_365_progreg(nm) = &  ! present-day month beginning day in a Gregorian year
        (/ 15.5000d0, 45.1213d0, 74.7425d0, 105.2425d0, 135.7425d0, 166.2425d0, 196.7425d0, 227.7425d0, 258.2425d0, &
        288.7425d0, 319.2425d0, 349.7425d0 /)
    real(8)     :: present_end_365_progreg(nm) = &  ! present-day month beginning day in a Gregorian year
        (/ 31.0000d0, 59.2425d0, 90.2425d0, 120.2425d0, 151.2425d0, 181.2425d0, 212.2425d0, 243.2425d0, 273.2425d0, &
        304.2425d0, 334.2425d0, 365.2425d0 /)

    logical                 :: debug_write = .false.  ! write additional diagnostic/debugging info
    logical                 :: debug_monlen = .false.
    
contains

subroutine get_month_lengths(calendar_type, begageBP, endageBP, agestep, nages, begyrCE, nsimyrs, &
    iageBP, iyearCE, imonlen, imonmid, imonbeg, imonend, rmonlen, rmonmid, rmonbeg, rmonend, VE_day, SS_day, ndays)

    use GISS_orbpar_subs
    use GISS_srevents_subs

    implicit none

    ! calendar type
    character(23), intent(in)   :: calendar_type

    ! simulation age-related variables (controls orbital elements)
    integer(4), intent(in)  :: begageBP                     ! beginning year (BP) (negative, e.g. 10 ka = -10000 BP)
    integer(4), intent(in)  :: endageBP                     ! ending year (BP)
    integer(4), intent(in)  :: agestep                      ! age step size
    integer(4), intent(in)  :: nages                        ! number of simulation ages

    ! individual model simulation year-related variables (controls equinox and solstice days and leap-year status)
    integer(4), intent(in)  :: begyrCE                      ! beginning (pseudo-) year of individual model simulation
    integer(4), intent(in)  :: nsimyrs                      ! number of years of simulation

    ! (output) month-length variables
    integer(4), intent(out) :: iageBP(nages*nsimyrs)        ! year BP 1950 (negative, e.g. 1900 CE = -50.0d0 BP)
    integer(4), intent(out) :: iyearCE(nages*nsimyrs)       ! yearCE simulation year
    integer(4), intent(out) :: imonlen(nages*nsimyrs,nm)    ! integer-value month lengths
    integer(4), intent(out) :: imonmid(nages*nsimyrs,nm)    ! integer-value mid-month days
    integer(4), intent(out) :: imonbeg(nages*nsimyrs,nm)    ! integer-value beginning days
    integer(4), intent(out) :: imonend(nages*nsimyrs,nm)    ! integer-value ending days
    real(8), intent(out)    :: rmonlen(nages*nsimyrs,nm)    ! real-value month lengths
    real(8), intent(out)    :: rmonmid(nages*nsimyrs,nm)    ! real-value mid-month days
    real(8), intent(out)    :: rmonbeg(nages*nsimyrs,nm)    ! real-value month beginning day
    real(8), intent(out)    :: rmonend(nages*nsimyrs,nm)    ! real-value month ending days
    real(8), intent(out)    :: VE_day(nages*nsimyrs)        ! real-value vernal equinox day in simulation year
    real(8), intent(out)    :: SS_day(nages*nsimyrs)        ! real-value (northern) summer solstice in simulation year
    integer(4), intent(out) :: ndays(nages*nsimyrs)         ! integer number of days in simulation year

    ! subroutine GISS_orbpars() and GISS_srevents() input and output arguments
    character(2)            :: year_type = 'BP'             ! AD (AD/BC), CE (CE/BCE), BP (before 1950)
    real(8)                 :: AgeBP                        ! age (BP 1950) (input)
    real(8)                 :: eccen                        ! eccentricity of orbital ellipse
    real(8)                 :: obliq_deg                    ! obliquity (degrees)
    real(8)                 :: perih_deg                    ! longitude of perihelion (degrees)
    real(8)                 :: precc                        ! climatological precession parameter = eccen * sin(omegvp)
    real(8)                 :: veqday                       ! (real) day of vernal equinox

    ! monlen() subroutine arguments
    real(8)                 :: yrlen                        ! real number of days in year (year length)
    integer(4)              :: ndyr                         ! integer number of days in year

    ! present-day month-length variables
    real(8)                 :: rmonlen_0ka(nm)              ! real calculated month lengths at 0 ka
    real(8)                 :: rmonlen_0ka_leap(nm)         ! real calculated month lengths at 0 ka in a leap year
    real(8)                 :: present_monlen(nm)           ! "present day" month lengths
    real(8)                 :: rmonbeg_0ka(nm)              ! real calculated month beginning day at 0ka
    real(8)                 :: rmonmid_0ka(nm)              ! real calculated month mid-month day at 0ka
    real(8)                 :: rmonend_0ka(nm)              ! real calculated month ending day at 0ka
    real(8)                 :: rmonbeg_0ka_leap(nm)         ! real calculated month beginning day at 0ka in a leap year
    real(8)                 :: rmonmid_0ka_leap(nm)         ! real calculated month mid-month day at 0ka in a leap year
    real(8)                 :: rmonend_0ka_leap(nm)         ! real calculated month ending day at 0ka in a leap year

    ! arrays for calculating various month-length statistics
    real(8)                 :: rmonlen_rel(nm)              ! difference between real month lengths and present
    real(8)                 :: rmonlen_ratio(nm)            ! ratio of real month lengths and present
    real(8)                 :: ryeartot(nages*nsimyrs)      ! real total number of days in year
    integer(4)              :: iyeartot(nages*nsimyrs)      ! integer total number of days in year
    
    ! indices
    integer(4)              :: n        ! simulation-age index
    integer(4)              :: i        ! simulation-year index
    integer(4)              :: ii       ! age and year index
    integer(4)              :: m        ! month index

    logical                 :: adjust_to_0ka = .true.   ! adjust month-lengths and start days to common values at 0 ka/1950 CE
    
    character(2048)         :: debugpath

    ! debugging output
    debugpath = "/project/cseg/mvr/angular/"
!    debugpath = "/Projects/Calendar/debug/"  ! Windows path  ! Windows path
    !debugpath = "/Users/bartlein/Projects/Calendar/PaleoCalAdjust/data/debug_files/"  ! Mac path
    if (debug_write) then
        open (22,file=trim(debugpath)//trim(calendar_type)//"_monlen_debug.dat")
        open (11,file=trim(debugpath)//trim(calendar_type)//"_cal_rmonlen_raw.dat")
        open (12,file=trim(debugpath)//trim(calendar_type)//"_cal_rmonlen_rel.dat")
        open (13,file=trim(debugpath)//trim(calendar_type)//"_cal_rmonlen_ratio.dat")
        open (23, file=trim(debugpath)//"kepler_test.dat")
        open (24, file=trim(debugpath)//"month_test.dat")
    end if

    ! check for supported calendar types
    select case (trim(calendar_type))
    case ('360_day','noleap','365_day','365.2425','proleptic_gregorian','progreg','gregorian','standard')
        continue
    case default
        stop "Calendar type not supported"
    end select
    
    ! ===================================================================================================================

    ! Step 1:  generate target years -- experiment ages (in yrs BP) and simulation years (in yrs CE)
    write (*,'(a)') "Generating target years..."
    if (debug_write) write (*,'("begageBP, endageBP, agestep, nages, begyrCE, nsimyrs: ",6i7)') &
        begageBP, endageBP, agestep, nages, begyrCE, nsimyrs
    if (debug_write) write (22,'("begageBP, endageBP, agestep, nages, begyrCE, nsimyrs: ",6i7)') &
            begageBP, endageBP, agestep, nages, begyrCE, nsimyrs

    ii=0
    do n = 1, nages
        do i = 1, nsimyrs
            ii = ii+1
            iageBP(ii) = begageBP + (n - 1) * agestep
            iyearCE(ii) = begyrCE + (i - 1)
            if (debug_write) write (22,'("n,i,ii,iageBP,iyearCE ", 5i8)') n,i,ii,iageBP(ii),iyearCE(ii)
        end do
    end do

    ! initialize arrays
    imonlen = 0; imonmid = 0; imonbeg = 0; imonend = 0
    rmonlen = 0.0d0; rmonmid = 0.0d0; rmonbeg = 0.0d0; rmonend = 0.0d0
    VE_day = 0.0d0; ndays = 0
    ryeartot = 0.0d0; iyeartot = 0

    ! ===================================================================================================================

    ! Step 2:  get 0 ka (1950CE) calculated month lengths and beginning, middle and end days, which can be used to adjust 
    ! all other calculated month lengths to nominal "present-day" values

    ! orbital elements for 0 ka
    write (*,'(a)') " 0 ka orbital elements..."
    ageBP = 0.0d0
    call GISS_orbpars('BP', ageBP, eccen, obliq_deg, perih_deg, precc)
    if (debug_write) write (22,'("ageBP, eccen, obliq_deg, perih_deg, precc: ",f10.1,4f17.12)') &
        ageBP, eccen, obliq_deg, perih_deg, precc

    ! Step 3:  0 ka calculated month lengths, also set subroutine arguments
    write (*,'(a)') " 0 ka month lengths..."
    select case (trim(calendar_type))
    case ('360_day')
        yrlen = dble(nd_360); ndyr = nd_360; veqday = veqday_360; present_monlen = present_mon_360
        call monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, rmonbeg_0ka, rmonmid_0ka, rmonend_0ka)
    case ('noleap', '365_day')
        yrlen = dble(nd_365); ndyr = nd_365; veqday = veqday_365; present_monlen = present_mon_noleap
        call monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, rmonbeg_0ka, rmonmid_0ka, rmonend_0ka)
    case ('366_day')
        yrlen = dble(nd_366); ndyr = nd_366; veqday = veqday_366; present_monlen = present_mon_leap
        call monlen(yrlen, ndyr, veqday,  int(present_monlen), eccen, perih_deg, rmonlen_0ka, rmonbeg_0ka, rmonmid_0ka, rmonend_0ka)
    case ('365.2425')
        yrlen = dble(progreg_year); ndyr = nd_365; veqday = veqday_365; present_monlen = present_mon_365_progreg
        call monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, rmonbeg_0ka, rmonmid_0ka, rmonend_0ka)
    case ('proleptic_gregorian','progreg','gregorian','standard')
        ! 365-day year
        yrlen = dble(nd_365); ndyr = nd_365; veqday = veqday_365; present_monlen = present_mon_noleap
        call monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, rmonbeg_0ka, rmonmid_0ka, rmonend_0ka)
        ! 366-day year
        yrlen = dble(nd_366); ndyr = nd_366; veqday = veqday_366; present_monlen = present_mon_leap
        call monlen(yrlen, ndyr, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka_leap, & 
            rmonbeg_0ka_leap, rmonmid_0ka_leap, rmonend_0ka_leap)
    case default
        stop "calendar type not supported"
    end select

    if (debug_write) write (22,'("rmonlen_0ka: ",21x,12f12.6)') rmonlen_0ka(1:nm)
    
    ! loop over simulation ages and years

    write (*,'(a)') "Looping over ages and years..."
    if (debug_write) write (*,'("nages, nsimyrs, nages*nsimyrs:  ",3i8)') nages, nsimyrs, nages*nsimyrs
    ii = 0
    do n = 1, nages
        do i = 1, nsimyrs
            ii = ii + 1
            if (mod(n,1000) .eq. 0) write (*,'(i8,$)') n; if (mod(n,15000).eq.0) write (*,'(" ")')
            if (debug_write) write (24,'(2i8)') iageBP(ii),iyearCE(ii)

            ! Step 4:  orbital elements for simulation age (e.g. 6 ka)
            call GISS_orbpars('BP', dble(iageBP(n)), eccen, obliq_deg, perih_deg, precc)
            
            if (debug_write) write (23,'("ageBP, eccen, obliq_deg, perih_deg, precc: ",f10.1,4f17.12)') &
                dble(iageBP(n)), eccen, obliq_deg, perih_deg, precc
            
            ! Steps 5&6:  get real-valued month lengths for different calendars (step 5) and adjust values to 0 ka (step 6)

            select case (trim(calendar_type))
                
            ! non time-varying calendars
            case ('360_day')
                yrlen = dble(nd_360); ndyr = nd_360; veqday = veqday_360
                call monlen(yrlen, ndyr, veqday, int(present_mon_360), eccen, perih_deg, rmonlen(ii,:), &
                    rmonbeg(ii,:), rmonmid(ii,:), rmonend(ii,:))
                ! adjust values so that 0 ka (1950 CE) will have nominal present-day month lengths
                if (adjust_to_0ka) call adjust_to_ref_length(rmonlen(ii,:), rmonlen_0ka, present_mon_360)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonbeg(ii,:), rmonbeg_0ka, present_beg_360)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonmid(ii,:), rmonmid_0ka, present_mid_360)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonend(ii,:), rmonend_0ka, present_end_360)
            case ('noleap', '365_day')
                yrlen = dble(nd_365); ndyr = nd_365; veqday = veqday_365
                call monlen(yrlen, ndyr, veqday, int(present_mon_noleap), eccen, perih_deg, rmonlen(ii,:), &
                    rmonbeg(ii,:), rmonmid(ii,:), rmonend(ii,:))
                ! adjust values so that 0 ka (1950 CE) will have nominal present-day month lengths
                if (adjust_to_0ka) call adjust_to_ref_length(rmonlen(ii,:), rmonlen_0ka, present_mon_noleap)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonbeg(ii,:), rmonbeg_0ka, present_beg_365)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonmid(ii,:), rmonmid_0ka, present_mid_365)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonend(ii,:), rmonend_0ka, present_end_365)
            case ('366_day')
                yrlen = dble(nd_366); ndyr = nd_366; veqday = veqday_366
                call monlen(yrlen, ndyr, veqday, int(present_mon_leap), eccen, perih_deg, rmonlen(ii,:), &
                    rmonbeg(ii,:), rmonmid(ii,:), rmonend(ii,:))
                ! adjust values so that 0 ka (1950 CE) will have nominal present-day month lengths
                if (adjust_to_0ka) call adjust_to_ref_length(rmonlen(ii,:), rmonlen_0ka, present_mon_leap)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonbeg(ii,:), rmonbeg_0ka, present_beg_366)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonmid(ii,:), rmonmid_0ka, present_mid_366)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonend(ii,:), rmonend_0ka, present_end_366)                
            case ('365.2425')
                yrlen = dble(progreg_year); ndyr = nd_365; veqday = veqday_365
                call monlen(yrlen, ndyr, veqday, int(present_mon_365_progreg), eccen, perih_deg, rmonlen(ii,:), &
                    rmonbeg(ii,:), rmonmid(ii,:), rmonend(ii,:))
                ! adjust values so that 0 ka (1950 CE) will have nominal present-day month lengths
                if (adjust_to_0ka) call adjust_to_ref_length(rmonlen(ii,:), rmonlen_0ka, present_mon_365_progreg)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonbeg(ii,:), rmonbeg_0ka, present_beg_365_progreg)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonmid(ii,:), rmonmid_0ka, present_mid_365_progreg)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonend(ii,:), rmonend_0ka, present_end_365_progreg)                
                
            ! proleptic_gregorian-like calendars
            case ('proleptic_gregorian','progreg','gregorian','standard')
                ! check for leap year
                nd_progreg = yearlen_CE(iyearCE(ii))
                if (debug_write) write (22,'("ii, iyearCE, nd_progreg: ",3i6,i4)') ii, iageBP(ii), iyearCE(ii), nd_progreg
                if (nd_progreg .eq. 365) then ! (non leap year)
                    yrlen = dble(nd_365); ndyr = nd_365; veqday = veqday_365
                    ! get real-valued month lengths
                    call monlen(yrlen, ndyr, veqday, int(present_mon_noleap), eccen, perih_deg, rmonlen(ii,:), &
                        rmonbeg(ii,:), rmonmid(ii,:), rmonend(ii,:))
                    ! adjust values so that 0 ka (1950 CE) will have nominal present-day month lengths
                    if (adjust_to_0ka) call adjust_to_ref_length(rmonlen(ii,:), rmonlen_0ka, present_mon_noleap)
                    if (adjust_to_0ka) call adjust_to_ref_day(rmonbeg(ii,:), rmonbeg_0ka, present_beg_365)
                    if (adjust_to_0ka) call adjust_to_ref_day(rmonmid(ii,:), rmonmid_0ka, present_mid_365)
                    if (adjust_to_0ka) call adjust_to_ref_day(rmonend(ii,:), rmonend_0ka, present_end_365)
                else ! ndprogreg = 366 (leap year)
                    yrlen = dble(nd_366); ndyr = nd_366; veqday = veqday_366
                    call monlen(yrlen, ndyr, veqday, int(present_mon_leap), eccen, perih_deg, rmonlen(ii,:), &
                        rmonbeg(ii,:), rmonmid(ii,:), rmonend(ii,:))
                    ! adjust values so that 0 ka (1950 CE) will have nominal present-day month lengths
                    if (adjust_to_0ka) call adjust_to_ref_length(rmonlen(ii,:), rmonlen_0ka_leap, present_mon_leap)
                    if (adjust_to_0ka) call adjust_to_ref_day(rmonbeg(ii,:), rmonbeg_0ka_leap, present_beg_366)
                    if (adjust_to_0ka) call adjust_to_ref_day(rmonmid(ii,:), rmonmid_0ka_leap, present_mid_366)
                    if (adjust_to_0ka) call adjust_to_ref_day(rmonend(ii,:), rmonend_0ka_leap, present_end_366)
                end if
            case default ! 
                stop "calendar_type"
            end select

            if (debug_write) write (22,'("ageBP, rmonlength:        ",i8,13f12.6)') iageBP(ii),rmonlen(ii,1:nm)
            if (debug_write) write (11,'(2i8,13f12.6)') iageBP(ii),iyearCE(ii),rmonlen(ii,1:nm) ! raw month lengths

            ! Step 7:  require the sum of month lengths each year to exactly equal the year length
            if (adjust_to_0ka) call adjust_to_yeartot(rmonlen(ii,:), yrlen, ryeartot(ii))

            ! Step 8:  integer month lengths
            call integer_monlen(rmonlen(ii,:), ndyr, imonlen(ii,:), iyeartot(ii))
            
            ! Step 9: integer beginning, middle and ending days
            call imon_begmidend(imonlen(ii,:), rmonbeg(ii,:), imonbeg(ii,:), imonmid(ii,:), imonend(ii,:))       
            
            ! Step 10: get vernal equinox and summer solstice days for ii-th simulation year
            select case (trim(calendar_type))
            case ('360_day')
                VE_day(ii) = veqday_360; SS_day(ii) = ssday_360; ndays(ii) = nd_360
            case ('noleap', '365_day', '365.2425')
                VE_day(ii) = veqday_365; SS_day(ii) = ssday_365; ndays(ii) = nd_365
            case ('366_day')
                VE_day(ii) = veqday_366; SS_day(ii) = ssday_366; ndays(ii) = nd_366
            case ('proleptic_gregorian','progreg','gregorian','standard')
                ! get vernal equinox and summer solstice days for model simulation year (not age)
                year_type = 'CE'
                call GISS_srevents(year_type, iyearCE(ii), progreg_year, VE_day(ii), SS_day(ii), perihelion, aphelion, ndays(ii))
                if (debug_write) write (24,'(i8,1x,5f12.6,1x,i4)') iyearCE(ii), progreg_year, VE_day(ii), SS_day(ii), &
                    perihelion, aphelion, ndays(ii)
            case default
                stop "paleo calendar type not defined"
            end select
            
            ! various month-length statistics
            do m=1,nm
                rmonlen_rel(m) = rmonlen(ii,m)-present_monlen(m)
                rmonlen_ratio(m) = rmonlen(ii,m)/present_monlen(m)
            end do

            if (debug_write) then
                write (12,'(2i8,13f12.8)') iageBP(ii),iyearCE(ii),rmonlen_rel
                write (13,'(2i8,13f12.8)') iageBP(ii),iyearCE(ii),rmonlen_ratio
            end if
            
        end do
    end do
    write (*,'(a)') " "

    if (debug_write) close(22); close(11); close(12); close(13)
    if (debug_monlen) close(23); close(24)

end subroutine get_month_lengths

subroutine monlen(yrlen, ndyr, veqday, imonlen, eccen, perih, rmonlen, rmonbeg, rmonmid, rmonend)
! calculate month lengths, and beginning, middle and ending day times for a particular orbital configuration, 
! calendar and vernal equinox day using a "traverse-time/time-of-flight" expression based on Kepler's equation:
! (Curtis, H.D. (2014) Orbital Mechanics for Engineering Students, Elsevier, Ch. 3)

    implicit none

    real(8), intent(in)     :: yrlen                ! total length of year (days)
    integer(4), intent(in)  :: ndyr                 ! number of days in year
    real(8), intent(in)     :: veqday               ! real vernal equinox day
    integer(4), intent(in)  :: imonlen(nm)          ! number of days in each month at present (calendar dependent)
    real(8), intent(in)     :: eccen, perih         ! orbital parameters

    real(8), intent(out)    :: rmonlen(nm)          ! real month lengths
    real(8), intent(out)    :: rmonbeg(nm)          ! real month beginning day
    real(8), intent(out)    :: rmonmid(nm)          ! real month middle day
    real(8), intent(out)    :: rmonend(nm)          ! real month ending day
    
    real(8)                 :: angle_veq_to_Jan1    ! angle between vernal equinox and Jan 1
    real(8)                 :: tt_perih_to_veq      ! traverse time, perihelion to vernal equinox
    real(8)                 :: angle_perih_to_Jan1  ! angle, perihelion to Jan 1
    real(8)                 :: angle_perih_to_veq   ! angle, perihelion to vernal equinox
    
    real(8)                 :: mon_angle(nm)        ! month angle (degrees)
    real(8)                 :: beg_angle(nm), mid_angle(nm), end_angle(nm) ! beginning, middle and end of month, relative to Jan 1
    
    ! angles and traverse times for beginning, middle and ending days of each month 
    real(8)                 :: perih_angle_month_beg(nm), tt_month_beg(nm), t_month_beg(nm)
    real(8)                 :: perih_angle_month_mid(nm), tt_month_mid(nm), t_month_mid(nm)
    real(8)                 :: perih_angle_month_end(nm), tt_month_end(nm), t_month_end(nm)

    real(8)     :: pi                               ! pi

    integer(4)  :: m                                ! indices

    pi=4.0d0*datan(1.0d0)

    rmonlen = 0.0d0; rmonbeg = 0.0d0; rmonmid = 0.0d0; rmonend = 0.0d0 
    
    ! number of days in year
    if (debug_write) write (*,'("number of days in year: ",i4)') ndyr
    
    ! perihelion (perihelion longitude = longitude of vernal equinox + 180 degrees, by definition)
    if (debug_write) write (*,'("perihelion: ",f14.8)') perih

    ! angle and traverse time, perihelion to vernal equinox 
    angle_perih_to_veq = 360.00d0 - perih
    if (debug_write) write (*,'("angle_perih_to_veq: ",f14.8)') angle_perih_to_veq

    ! traverse time (along elliptical orbit) perihelion to vernal equinox
    call kepler_time(eccen, yrlen, angle_perih_to_veq, tt_perih_to_veq)
    if (debug_write) write (*,'("tt_perih_to_veq: ",f14.8)') tt_perih_to_veq

    ! angle, perihelion to Jan1
    angle_perih_to_Jan1 = angle_perih_to_veq - veqday * (360.0d0/yrlen)
    if (debug_write) write (*,'("angle_perih_to_Jan1: ",f14.8)') angle_perih_to_Jan1

    ! angle, vernal equinox to Jan1 (Jan 1 begins at 0.0 (or 360.0 degrees))
    ! for checking, angle_perih_to_Jan1 * (yrlen/360.0d0) should equal 0.0
    angle_veq_to_Jan1 = 360.0d0 - veqday * (360.0d0/yrlen)
    if (debug_write) write (*,'("angle_veq_to_Jan1, check: ",2f14.8)') angle_veq_to_Jan1, angle_perih_to_Jan1 * (yrlen/360.0d0) 

    ! angles, traverse times, month lengths, etc. for individual months

    ! month angle (and beginning, middle and end of month, relative to Jan1)
    if (debug_write) write (*,'(a)') " "
    mon_angle(1) = dble(imonlen(1) * (360.0d0/yrlen))
    beg_angle(1) = 0.0d0
    mid_angle(1) = beg_angle(1) + mon_angle(1)/2.0d0
    end_angle(1) = mon_angle(1)
    do m=2,nm
        mon_angle(m) = dble(imonlen(m)) * (360.0d0/yrlen)
        beg_angle(m) = beg_angle(m-1) + mon_angle(m-1)
        mid_angle(m) = beg_angle(m) + mon_angle(m) / 2.0d0
        end_angle(m) = beg_angle(m) + mon_angle(m)
    end do
    if (debug_write) then
        do m = 1,nm
            write (*,'(i4,5f9.4)') m,dble(imonlen(m)), mon_angle(m), beg_angle(m), mid_angle(m), end_angle(m)
        end do
    end if

    ! angles from perihelion etc. for month beginning, mid, and ending days
    do m = 1, nm
    
        ! angle from perihelion for each month's beginning, middle and ending days (on circular orbit)
        perih_angle_month_beg(m) = angle_perih_to_Jan1 + beg_angle(m) 
        perih_angle_month_mid(m) = angle_perih_to_Jan1 + mid_angle(m) 
        perih_angle_month_end(m) = angle_perih_to_Jan1 + end_angle(m)  
    
        ! traverse times from perihelion for each month's beginning, middle and ending days (on elliptical orbit)
        call kepler_time(eccen, yrlen, perih_angle_month_beg(m), tt_month_beg(m))
        call kepler_time(eccen, yrlen, perih_angle_month_mid(m), tt_month_mid(m))
        call kepler_time(eccen, yrlen, perih_angle_month_end(m), tt_month_end(m))
    
        ! traverse time from vernal equinox for each month's beginning, middle and ending days (on elliptical orbit)
        t_month_beg(m) = tt_month_beg(m) - tt_perih_to_veq
        t_month_mid(m) = tt_month_mid(m) - tt_perih_to_veq
        t_month_end(m) = tt_month_end(m) - tt_perih_to_veq
    
        ! beginning, middle and ending days of each month (relative to Jan1)
        rmonbeg(m) = dmod(t_month_beg(m) + veqday, yrlen)
        if (perih_angle_month_beg(m) .gt. 360.0d0) rmonbeg(m) = rmonbeg(m) + yrlen
        rmonmid(m) = dmod(t_month_mid(m) + veqday, yrlen)
        if (perih_angle_month_mid(m) .gt. 360.0d0) rmonmid(m) = rmonmid(m) + yrlen
        rmonend(m) = dmod(t_month_end(m) + veqday, yrlen) 
        if (perih_angle_month_end(m) .gt. 360.0d0) rmonend(m) = rmonend(m) + yrlen
    
    end do
    
    ! fixup for ending day of last month (when end of last month appers in beginning of year)
    if (rmonend(nm) .lt. 30.0d0) rmonend(nm) = rmonend(nm) + yrlen  

    if (debug_write) then
        do m=1,nm
            write (*,'("m, perih_angle_month_beg, tt_month_beg, t_month_beg, beg_day: ", i3, 4f14.8)') &
                m, perih_angle_month_beg(m), tt_month_beg(m), t_month_beg(m), rmonbeg(m)
        end do
        do m=1,nm
            write (*,'("m, perih_angle_month_mid, tt_month_mid, t_month_mid, mid_day: ", i3, 4f14.8)') &
                m, perih_angle_month_mid(m), tt_month_mid(m), t_month_mid(m), rmonmid(m)
        end do
        do m=1,nm
            write (*,'("m, perih_angle_month_end, tt_month_end, t_month_end, end_day: ", i3, 4f14.8)') &
                m, perih_angle_month_end(m), tt_month_end(m), t_month_end(m), rmonend(m)
        end do
    end if

    ! month length (month beginning to next month beginning)
    do m = 1,nm-1
        rmonlen(m) = t_month_beg(m+1) - t_month_beg(m)
        if (rmonlen(m) .le. 0.0d0) rmonlen(m) = rmonlen(m) + yrlen
        if (debug_write) write (*,'("m, rmonlen: ", i4,4f14.8)') m, rmonlen(m)
    end do
    rmonlen(nm) = t_month_beg(1) - t_month_beg(nm)
    if (rmonlen(nm) .le. 0.0d0) rmonlen(nm) = rmonlen(nm) + yrlen
    if (debug_write) write (*,'("m, rmonlen: ", i4,4f14.8)') nm, rmonlen(nm)

end subroutine monlen
    
subroutine kepler_time(eccen, T, theta_deg, time)
! travel time along orbit relative to periapsis/perihelion (i.e. 0.0 at perihelion)
! input:  true anomaly (theta, degrees); output:  traverse time since perihelion (time, same units as T)
! (Curtis, H.D. (2014) Orbital Mechanics for Engineering Students, Elsevier, Ch. 3)

    implicit none
    
    real(8), intent(in)     :: eccen
    real(8), intent(in)     :: T                ! orbital period (yrlen)
    real(8), intent(in)     :: theta_deg        ! "true" anomaly (angle from perihelion (degrees)
    real(8), intent(out)    :: time             ! traverse time from periapsis (e.g. perihelion)
    
    !real(8)                 :: radians
    real(8)                 :: theta_rad        ! theta_rad (radians)
    real(8)                 :: M                ! mean anomaly
    real(8)                 :: E                ! eccentric anomaly
    real(8)                 :: pi
    !real(8)                 :: a, r             ! semi-major axis length, polar-coordinate r (for plotting)
    
    pi=4.0d0*datan(1.0d0)  
    !a = 1.0d0
    
    theta_rad = radians(theta_deg)
    E = 2.0d0 * datan( dsqrt((1.0d0 - eccen) / (1.0d0 + eccen)) * dtan(theta_rad/ 2.0d0) )  ! eq 3.13b
    M = E - eccen*dsin(E)                                                                   ! eq 3.14
    time = (M / (2.0d0 * pi)) * T                                                           ! eq 3.15
    if (time .lt. 0.0d0) time = time + T

    ! for plotting orbit using polar coordinates
    !r = a * (1.0d0 - eccen**2) / (1.0d0 + eccen * dcos(theta_rad))    
    
end subroutine kepler_time    
    
subroutine kepler_theta(eccen, M_angle, theta_deg)
! angular position along orbit
! input:  Mean anomaly (M_angle, degrees); output:  true anomaly (theta, degrees)
! (Curtis, H.D. (2014) Orbital Mechanics for Engineering Students, Elsevier, Ch. 3)

    implicit none
    
    real(8), intent(in)     :: eccen
    real(8), intent(in)     :: M_angle          ! Mean anomaly (degrees)
    real(8), intent(out)    :: theta_deg        ! true anomaly angular position along elliptical orbit (degrees)
    
    real(8)                 :: M                ! Mean anomaly (radians)
    real(8)                 :: E                ! Eccentric anomaly (radians)
    
    real(8)                 :: theta            ! angular position along orbit (radians)
    real(8)                 :: pi
    real(8)                 :: a, r             ! semi-major axis length, polar-coordinate r (for plotting)
    real(8)                 :: tol = 1.0d-08
    real(8)                 :: E_ratio, fE1, fE2

    pi=4.0d0*datan(1.0d0)
    a = 1.0d0
    
    M = radians(M_angle)
    
    ! Newton's method
    ! initial value of E
    if (M .lt. pi) then 
        E = M + eccen/2.0d0
    else
        E = M - eccen/2.0d0
    endif 
    
    ! iterate
    E_ratio = 1.0d0
    do while (dabs(E_ratio) .gt. tol)
        fE1 = (E - eccen * dsin(E) - M)
        fE2 = 1.0d0 - eccen * dcos(E)
        E_ratio = fE1 / fE2
        E = E - E_ratio
        ! write (*,*) fE1, fE2, E_ratio, E
    end do 
    
    if (debug_write) then
        ! check
        M = E - eccen * dsin(E)
        write (*,*) radians(M_angle), M
    end if
    
    theta= 2.0d0 * datan(dsqrt((1.0d0 + eccen)/(1 - eccen)) * dtan(E/2.0d0))
    r = a * (1.0d0 - eccen**2) / (1.0d0 + eccen * dcos(theta))

    theta_deg = degrees(theta)
    if (theta_deg .lt. 0.0d0) then
        theta_deg = theta_deg + 360.0d0
    else
        theta_deg = theta_deg
    end if
    
end subroutine kepler_theta


subroutine kepler_daylen(ndyr, yrlen, eccen, veqday, angle_perih_to_Jan1, tt_perih_to_veq, & 
    perih_angle_day, tt_day, t_day, daylen, start_day)

    implicit none
    
    integer(4), intent(in)  :: ndyr                 ! integer number of days in year
    real(8), intent(in)     :: yrlen                ! length of year
    real(8), intent(in)     :: eccen                ! eccentricity
    real(8), intent(in)     :: veqday               ! vernal equinox day
    real(8), intent(in)     :: angle_perih_to_Jan1  ! angle (on circle) from perihelion to Jan1
    real(8), intent(in)     :: tt_perih_to_veq      ! traverse time, from perihelion to veq (same units as yrlen)
    
    real(8), intent(out)    :: perih_angle_day(0:ndyr)  ! angle from Jan 1
    real(8), intent(out)    :: tt_day(0:ndyr)           ! traverse time, from perihelion to each day (same units as yrlen)
    real(8), intent(out)    :: t_day(0:ndyr)            ! traverse time, from vernal equinox (same units as yrlen)
    real(8), intent(out)    :: daylen(0:ndyr)           ! length of day (same units as yrlen)
    real(8), intent(out)    :: start_day(0:ndyr)        ! start time of each day (same units as yrlen)
    
    integer(4)              :: i
    
    ! angle from perihelion ("mean anomaly") for each day, beginning Jan 1
    do i=1,ndyr
        perih_angle_day(i) = angle_perih_to_Jan1 + (dble(i-1)*(360.0d0/yrlen)) ! i - 1 puts Jan 1 at 0.0 degrees
        perih_angle_day(i) = dmod(perih_angle_day(i), 360.0d0)
    end do
    perih_angle_day(0) = perih_angle_day(ndyr)

    ! traverse time of each day from perihelion (tt_day) and traverse time from vernal equinox (t)
    do i = 1,ndyr
        call kepler_time(eccen, yrlen, perih_angle_day(i), tt_day(i))
        t_day(i) = tt_day(i) - tt_perih_to_veq
        !if (i .gt. veqday) t_day(i) = t_day(i) + yrlen
    end do
    tt_day(0) = tt_day(ndyr); t_day(0) = t_day(ndyr)

    ! day lengths, and start of each day
    daylen(0) = t_day(0) - t_day(ndyr-1)
    start_day(0) = t_day(0) + veqday
    if (start_day(0) .gt. yrlen) start_day(0) = start_day(0) - yrlen
    do i = 1,ndyr
        daylen(i) = t_day(i) - t_day(i-1)
        if (daylen(i) .le. 0.0d0) daylen(i) = daylen(i) + yrlen
        start_day(i) = start_day(i-1) + daylen(i)
    end do
    
end subroutine kepler_daylen

subroutine adjust_to_ref_length(rmonlen, rmonlenref, rmonlentarg)

    implicit none

    real(8), intent(inout)  :: rmonlen(nm)      ! real month lengths
    real(8), intent(in)     :: rmonlenref(nm)   ! reference month lengths (usually calculated 0 ka)
    real(8), intent(in)     :: rmonlentarg(nm)  ! target month lengths (usually conventional 0 ka)
    
    real(8)                 :: rmonlentot_in, rmonlentot_out

    integer(4)              :: m

    rmonlentot_in = 0.0d0; rmonlentot_out = 0.0d0
    ! adjust rmonlen to reference year
    do m=1,nm
        rmonlentot_in = rmonlentot_in + rmonlen(m)
        rmonlen(m) = rmonlen(m) * (rmonlentarg(m)/rmonlenref(m))
        rmonlentot_out = rmonlentot_out + rmonlen(m)
    end do
    
end subroutine adjust_to_ref_length

subroutine adjust_to_ref_day(rmonday, rdayref, rdaytarg)

    implicit none

    real(8), intent(inout)  :: rmonday(nm)      ! real beginning, middle or ending days
    real(8), intent(in)     :: rdayref(nm)      ! reference day (usually calculated 0 ka)
    real(8), intent(in)     :: rdaytarg(nm)     ! target day (usually conventional 0 ka)

    integer(4)              :: m

    ! adjust rmonday to reference year 
    do m=1,nm
        rmonday(m) = rmonday(m) - (rdayref(m) - rdaytarg(m))
    end do

end subroutine adjust_to_ref_day

subroutine adjust_to_yeartot(rmonlen, ryeartottarg, ryeartot)

    implicit none

    real(8), intent(inout)  :: rmonlen(nm)      ! real month lengths
    real(8), intent(in)     :: ryeartottarg     ! total annual month lengths target
    real(8), intent(out)    :: ryeartot         ! total annual month lengths

    integer(4)              :: m

    ! get ryeartot
    ryeartot=0.0d0
    do m=1,nm
        ryeartot = ryeartot + rmonlen(m)
    end do

    ! adjust rmonlen to ryeartot target
    do m=1,nm
        rmonlen(m) = rmonlen(m) * (ryeartottarg/ryeartot)
    end do

    ! recalc ryeartot
    ryeartot=0.0d0
    do m=1,nm
        ryeartot = ryeartot + rmonlen(m)
    end do

end subroutine adjust_to_yeartot
                
subroutine compare_monthdefs(nyrs, imondef, rmondef)
! used for debugging

    implicit none

    integer(4), intent(in)  :: nyrs
    integer(4), intent(in)  :: imondef(nyrs,nm)    ! integer month definition
    real(8), intent(in)     :: rmondef(nyrs,nm)    ! real month definition

    real(8)                 :: diff(nm),ssq(nm),mse(nm)
    integer(4)              :: n,m

    diff=0.0d0; ssq=0.0d0; mse=0.0d0

    do m=1,nm
        do n=1,nyrs
            diff(m) = diff(m)+dble(imondef(n,m))-rmondef(n,m)
            ssq(m) = ssq(m)+(dble(imondef(n,m))-rmondef(n,m))*(dble(imondef(n,m))-rmondef(n,m))
        end do
    end do
    mse(:) = dsqrt(ssq(:)/dble(nyrs))

    write (24,'("diff: ",12f12.6)') diff(1:nm)
    write (24,'(" ssq: ",12f12.6)') ssq(1:nm)
    write (24,'(" mse: ",12f12.6)') mse(1:nm)

end subroutine compare_monthdefs

integer(4) function yearlen_BP(ageBP)
! gets number of days in a BP year -- no year zero or Gregorian-Julian adjustment

    implicit none

    integer(4), intent(in)  :: ageBP   ! Year BP (negative, e.g. 1900 CE = -50 BP)
    integer(4)              :: yearCE

    yearCE = ageBP + 1950 ! no Year 0 BCE/CE adjustment

    yearlen_BP = 365
    if (mod(yearCE, 4) .eq. 0) yearlen_BP = 366
    if (mod(yearCE, 100) .eq. 0) yearlen_BP = 365
    if (mod(yearCE, 400) .eq. 0) yearlen_BP = 366

end function yearlen_BP

integer(4) function yearlen_CE(yearCE)
! gets number of days in a CE year -- no year zero or Gregorian-Julian adjustment

    implicit none

    integer(4), intent(in)  :: yearCE   ! YearCE

    yearlen_CE = 365
    if (mod(yearCE, 4) .eq. 0) yearlen_CE = 366
    if (mod(yearCE, 100) .eq. 0) yearlen_CE = 365
    if (mod(yearCE, 400) .eq. 0) yearlen_CE = 366

end function yearlen_CE

subroutine integer_monlen(rmonlen, ndtarg, imonlen, iyeartot)

    implicit none

    real(8), intent(in)     :: rmonlen(nm)      ! real month lengths
    integer(4), intent(in)  :: ndtarg           ! total annual integer number of days target
    integer(4), intent(out) :: imonlen(nm)      ! integer month lengths
    integer(4), intent(out) :: iyeartot

    integer(4)              :: diff_sign
    real(8)                 :: ryeartot
    real(8)                 :: inc
    integer(4)              :: i,m

    ! integer month lengths
    iyeartot=0; ryeartot=0.0d0
    do  m=1,nm
        imonlen(m)=idint(dnint(rmonlen(m)))
        iyeartot=iyeartot+imonlen(m)
        ryeartot=ryeartot+rmonlen(m)
    end do
    !write (23,'(13(i4,7x))') imonlen(1:nm),iyeartot

    ! force integer month lengths to sum to ndtarg
    if (iyeartot.ne.ndtarg) then

        if (debug_write) write (22,'(13f12.6)') rmonlen(1:nm),ryeartot
        if (debug_write) write (22,'(13(i4,7x))') imonlen(1:nm),iyeartot

        ! add (diff_sign = 1) or subtract (diff_sign = -1) a small increment to each rmonlen value
        ! iterate until iyeartot = ndtarg, incrementing small increment via i (i*inc)
        diff_sign=1
        if ((iyeartot-ndtarg) .gt. 0) diff_sign=-1
        !write (20,*) diff_sign
        inc=0.000001; i=0
        do
            if (iyeartot.eq.ndtarg) exit
            i=i+1
            iyeartot=0
            do m=1,nm
                imonlen(m)=idint(dnint(rmonlen(m)+i*inc*diff_sign))
                iyeartot=iyeartot+imonlen(m)
            end do

            diff_sign=1
            if ((iyeartot-ndtarg) .gt. 0) diff_sign=-1

        end do

        if (debug_write) write (22,'(13(i4,7x))') imonlen(1:nm),iyeartot
        if (debug_write) write (22,'(a)')
    end if

end subroutine integer_monlen

subroutine imon_begmidend(imonlen, rmonbeg, imonbeg, imonmid, imonend)

    implicit none

    integer(4), intent(in)  :: imonlen(nm)  ! integer month length
    real(8), intent(in)     :: rmonbeg(nm)  ! real month beginning day
    integer(4), intent(out) :: imonbeg(nm)  ! integer beginning day of month
    integer(4), intent(out) :: imonmid(nm)  ! integer mid-month day
    integer(4), intent(out) :: imonend(nm)  ! integer ending day of month

    integer(4)              :: m

    imonbeg(1) = idint(rmonbeg(1))
    imonend(1) = imonbeg(1) + imonlen(1) - 1
    do m=2,nm
        imonbeg(m) = imonend(m-1) + 1
        imonend(m) = imonend(m-1) + imonlen(m)
    end do
    do m=1,nm
        imonmid(m) = imonbeg(m) + imonlen(m)/2
    end do

end subroutine imon_begmidend

real(8) function radians(d)
! decimal degrees to radians

    implicit none

    real(8) d, pi
    pi=4.0d0*datan(1.0d0)
    radians=d*((2.0d0*pi)/360.0d0)

end function radians

real(8) function degrees(r)
! radians to decimal degrees

    implicit none

    real(8) r, pi
    pi=4.0d0*datan(1.0d0)
    degrees=r/((2.0d0*pi)/360.0d0)

end function degrees

end module month_length_subs
    
