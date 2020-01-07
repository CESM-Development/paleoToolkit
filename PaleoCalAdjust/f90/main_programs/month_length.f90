program month_length
! Calculates the length of months under different paleo calendars, using 
! - orbital-element parameters calculated by GISS ORBPAR.FOR (https://data.giss.nasa.gov/ar5/SOLAR/ORBPAR.FOR) using 
!   Berger (1978, JAS 35:2362-2367) algorithm and tables;
! - "solar events" (i.e. vernal equinox day) calculated by GISS SREVENTS.FOR (https://data.giss.nasa.gov/ar5/SOLAR/SREVENTS.FOR)
!   (Orbital-element parameter calculations are thought to be valid over the past and future 1.0 Myr)
! - generalization of the 360-day month-length algorithm from Kutzbach and Gallimore (1988, JGR 93(D1):803-821), and
! - application of Kepler's "travel time" or "time-of-flight" equation, after
!   Curtis, H.D. (2014, Orbital Mechanics for Engineering Students, Elsevier, Ch. 3).
   
! The program requires the modules month_length_subs.f90, GISS_orbpar_subs.f90 and GISS_srevents_subs.f90
! GISS_orbpar_subs.f90 and GISS_srevents_subs.f90 are based on GISS ORBPAR.FOR and SREVENTS.FOR, which are retrievable from:
! https://web.archive.org/web/20150920211936/http://data.giss.nasa.gov/ar5/solar.html

! An info .csv file (with an appropriate header) containing the following information is read:
! prefix        :: (string) a short name for the output month-length table
! calendar_type :: (string) calendar type (e.g. "noleap", "proleptic_gregorian", etc.)
! begageBP      :: (integer) beginning age of the table (e.g. 21000 (= 21 ka)
! endageBP      :: (integer) ending age of the table (e.g. 0 (= 0ka)
! agestep       :: (integer) interval between age calculations
! begyrCE       :: (integer) beginning calendar year of simulation for multi-year simulations at each age
! nsimyrs       :: (integer) number of calendar years
    
! Author: Patrick J. Bartlein, Univ. of Oregon (bartlein@uoregon.edu), with contributions by S.L. Shafer (sshafer@usgs.gov)
!
! Version: 1.0d
! Last update: 2019-08-04

use month_length_subs
    
implicit none

! simulation age-related variables (controls orbital elements)
integer(4)              :: begageBP             ! beginning year (BP) (negative, e.g. 10 ka = -10000 BP)
integer(4)              :: endageBP             ! ending year (BP) 
integer(4)              :: agestep              ! age step size
integer(4)              :: nages                ! number of simulation ages 
integer(4), allocatable :: iageBP(:)            ! (nages*nsimyrs) year BP 1950 (negative, e.g. 1900 CE = -50 years BP 1950)

! individual model simulation year-related variables (controls vernal equinox day and leap-year status)
integer(4)              :: begyrCE              ! beginning (pseudo-) year of individual model simulation
integer(4)              :: nsimyrs              ! number of years of simulation
integer(4), allocatable :: iyearCE(:)           ! (nages*nsimyrs) yearCE simulation year (e.g. 1850CE, 850CE, etc.)

! month-length variables
integer(4), allocatable :: imonlen(:,:),imonmid(:,:)    ! (nages*nsimyrs,nm) integer-value month lengths and mid days
integer(4), allocatable :: imonbeg(:,:),imonend(:,:)    ! (nages*nsimyrs,nm) integer-value month beginning and ending days 
real(8), allocatable    :: rmonlen(:,:),rmonmid(:,:)    ! (nages*nsimyrs,nm) real-value month lengths and mid days 
real(8), allocatable    :: rmonbeg(:,:),rmonend(:,:)    ! (nages*nsimyrs,nm) real-value month beginning and ending days 
real(8), allocatable    :: VE_day(:)                    ! (nages*nsimyrs) vernal equinox day in simulation year
real(8), allocatable    :: SS_day(:)                    ! (nages*nsimyrs) (northern) summer solstice day in simulation year
integer(4), allocatable :: ndays(:)                     ! (nages*nsimyrs) number of days in year

! calendar type
character(32)           :: calendar_type               

! other variables
character(3)    :: monname(nm)                  ! names of months
character(256)  :: header(8)                    ! headers
character(64)   :: prefix                       ! file name prefix
character(1)    :: csvheader                    ! info .csv file header
character(2048) :: infopath                     ! info file path
character(64)   :: infofile                     ! info file name
character(2048) :: outpath                      ! output file path

! indices
integer(4)              :: n        ! simulation-age index
integer(4)              :: i        ! simulation-year index
integer(4)              :: ii       ! age and year index
integer(4)              :: m        ! month index
integer(4)              :: iostatus ! IOSTAT value

data monname /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

outpath = "/Projects/Calendar/PaleoCalAdjust/data/month_lengths/"  ! Windows path
infopath = "/Projects/Calendar/PaleoCalAdjust/data/info_files/"    ! Windows path
!outpath = "/Users/bartlein/Projects/Calendar/PaleoCalAdjust/data/month_lengths/"  ! Mac path
!infopath = "/Users/bartlein/Projects/Calendar/PaleoCalAdjust/data/info_files/"    ! Mac path
infofile = 'month_length_info.csv' ! 

! open the info file, and loop over specified calendar tables
! past ages are negative, e.g. 21 ka = 21,000 cal yr BP = -21000, and 1950 CE = 0 cal yr BP = 0 here

open (3,file=trim(infopath)//trim(infofile))
read (3,'(a)') csvheader

iostatus = 1
do
    read (3,*,iostat=iostatus) prefix, calendar_type, begageBP, endageBP, agestep, begyrCE, nsimyrs
    !write (*,'("iostatus = ",i2)') iostatus
    if (iostatus.lt.0) then
        write (*,'(a)') "Done (month_length)"
        exit
    end if
    write (*,'(125("="))')   
    write (*,'("prefix: ", a)') trim(prefix)
    write (*,'("begageBP, endageBP, agestep, begyrCE, nsimyrs, calendar_type: ", 5i7, 2x, a)') &
        begageBP, endageBP, agestep, begyrCE, nsimyrs, trim(adjustl(calendar_type))

    ! output files
    ! generate file headers
    header = ""; header(1) = "AgeBP,   YearCE,"; header(5) = "AgeBP,   YearCE,";
    do m=1,nm
        header(1)=trim(header(1))//" "//trim(monname(m))//"Days     ,"  ! real month-length .csv file header
        header(2)=trim(header(2))//" "//trim(monname(m))//"Mid      ,"
        header(3)=trim(header(3))//" "//trim(monname(m))//"Beg      ,"
        header(4)=trim(header(4))//" "//trim(monname(m))//"End      ,"
        header(5)=trim(header(5))//" "//trim(monname(m))//"Days,"       ! integer month-length .csv file header
        header(6)=trim(header(6))//" "//trim(monname(m))//"Mid ,"
        header(7)=trim(header(7))//" "//trim(monname(m))//"Beg ,"
        header(8)=trim(header(8))//" "//trim(monname(m))//"End ,"
    end do

    ! open files
    open (1,file=trim(outpath)//trim(prefix)//"_cal_"//trim(calendar_type)//"_rmonlen.csv")
    open (2,file=trim(outpath)//trim(prefix)//"_cal_"//trim(calendar_type)//"_imonlen.csv")
    write (1,'(a)') "   "//trim(header(1))//trim(header(2))//trim(header(3))//trim(header(4))//" VE_day , SS_day , ndays"
    write (2,'(a)') "   "//trim(header(5))//trim(header(6))//trim(header(7))//trim(header(8))//" VE_day , SS_day , ndays"

    ! allocate arrays
    nages = (endageBP - begageBP)/agestep + 1
    allocate (iageBP(nages*nsimyrs), iyearCE(nages*nsimyrs))
    allocate (imonlen(nages*nsimyrs,nm),imonmid(nages*nsimyrs,nm),imonbeg(nages*nsimyrs,nm),imonend(nages*nsimyrs,nm))
    allocate (rmonlen(nages*nsimyrs,nm),rmonmid(nages*nsimyrs,nm),rmonbeg(nages*nsimyrs,nm),rmonend(nages*nsimyrs,nm))
    allocate (VE_day(nages*nsimyrs),SS_day(nages*nsimyrs),ndays(nages*nsimyrs))

    ! initialize arrays
    imonlen = 0; imonmid = 0; imonbeg = 0; imonend = 0
    rmonlen = 0.0d0; rmonmid = 0.0d0; rmonbeg = 0.0d0; rmonend = 0.0d0
    VE_day = 0.0d0; SS_day = 0.0d0; ndays = 0

    ! get month lengths
    write (*,'(a)') "Getting month lengths..."
    call get_month_lengths(calendar_type, begageBP, endageBP, agestep, nages, begyrCE, nsimyrs, & 
        iageBP, iyearCE, imonlen, imonmid, imonbeg, imonend, rmonlen, rmonmid, rmonbeg, rmonend, VE_day, SS_day, ndays)

    ! write out the data
    write (*,'(a)') "Writing month lengths..."
    ii = 0
    do n = 1, nages 
        do i = 1, nsimyrs
            ii = ii + 1
            select case(prefix)
                case ('tr21', 'PMIP4', 'GENMOM', 'PMIP360')
                    iyearCE(ii) = iageBP(ii) + 1950      ! put actual YearCE into table
                case default
                    continue
            end select
            write (1,'(i8,", ",i8,48(", ",f12.8),2(", ",f7.3),", ",i4)') &
                iageBP(ii),iyearCE(ii),rmonlen(ii,1:nm),rmonmid(ii,1:nm), rmonbeg(ii,1:nm),rmonend(ii,1:nm),VE_day(ii), &
                    SS_day(ii),ndays(ii)
            write (2,'(i8,", ",i8,48(", ",3x,i4),2(", ",f7.3),", ",i4)') &
                iageBP(ii),iyearCE(ii),imonlen(ii,1:nm),imonmid(ii,1:nm), imonbeg(ii,1:nm),imonend(ii,1:nm),VE_day(ii), &
                    SS_day(ii),ndays(ii)
        end do
    end do

    close(1); close(2)
    
    deallocate(iageBP, iyearCE, imonlen, imonmid, imonbeg, imonend, rmonlen, rmonmid, rmonbeg, rmonend, VE_day, SS_day, ndays)
    
    write (*,*) " "

end do

end program month_length
