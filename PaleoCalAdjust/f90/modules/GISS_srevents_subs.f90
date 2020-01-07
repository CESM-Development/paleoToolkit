module GISS_srevents_subs
! subroutines based on
! SREVENTS.FOR    Solar EVENTS each year    2012/05/29
! https://data.giss.nasa.gov/ar5/SOLAR/SREVENTS.FOR downloaded 2017-09-12
! also retrievable from https://web.archive.org/web/20150920211936/http://data.giss.nasa.gov/ar5/solar.html   
! conversion to .f90 by P.J. Bartlein

    implicit none

    ! components of event dates
    integer(4)      :: JVEYR,JVEMON,JVEDAT,JVEHR,JVEMIN, JSSYR,JSSMON,JSSDAT,JSSHR,JSSMIN
    integer(4)      :: JAEYR,JAEMON,JAEDAT,JAEHR,JAEMIN, JWSYR,JWSMON,JWSDAT,JWSHR,JWSMIN
    integer(4)      :: JPRYR,JPRMON,JPRDAT,JPRHR,JPRMIN, JAPYR,JAPMON,JAPDAT,JAPHR,JAPMIN
    integer(4)      :: JEXYR,JEXMON,JEXDAT,JEXHR,JEXMIN
    integer(4)      :: KPERIH, KAPHEL

    real(8)         :: vereqx

contains

subroutine GISS_srevents(year_type, iyear, EDAYzY, veqday, ssday, perihelion, aphelion, ndays_in_year)
! subroutines based on
! SREVENTS.FOR    Solar EVENTS each year    2012/05/29
! https://data.giss.nasa.gov/ar5/SOLAR/SREVENTS.FOR downloaded 2017-09-12

    use giss_orbpar_subs

    implicit none

    integer(4), parameter       :: nm=12
    character(2), intent(in)    :: year_type
    integer(4), intent(in)      :: iyear
    real(8), intent(in)         :: EDAYzY
    real(8), intent(out)        :: veqday, ssday, perihelion, aphelion
    integer(4), intent(out)     :: ndays_in_year

    integer(4)                  :: YearCE, YearBP

!  The unit (days) means days measured since 2000 January 1, hour 0

    real(8)         :: year, eccen, obliq, omegvp
    real(8)         :: bsemi
    real(8)         :: TAofVE, EAofVE, MAofVE
    real(8)         :: PERIH1, PERIH2, APHEL1, APHEL2
    real(8)         :: TAofSS, EAofSS, MAofSS
    real(8)         :: SUMSOL
    real(8)         :: TAofAE, EAofAE, MAofAE
    real(8)         :: AUTEQX
    real(8)         :: TAofWS, EAofWS, MAofWS
    real(8)         :: WINSOL

    real(8)         :: pi, twopi

    real(8)     :: nd_noleap(nm) = &       ! present-day month lengths in 365-day (noleap) year
        (/ 31.0d0, 28.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    real(8)     :: nd_leap(nm) = &         ! present-day month lengths in 366-day (leap) year
        (/ 31.0d0, 29.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    real(8)     :: accumday_noleap(nm) = &       ! present-day accumulated day lengths in 365-day (noleap) year
        (/  0.0d0, 31.0d0, 59.0d0, 90.0d0,120.0d0,151.0d0,181.0d0,212.0d0,243.0d0,273.0d0,304.0d0,334.0d0 /)
    real(8)     :: accumday_leap(nm) = &         ! present-day accumulated day lengths in 366-day (leap) year
        (/  0.0d0, 31.0d0, 60.0d0, 91.0d0,121.0d0,152.0d0,182.0d0,213.0d0,244.0d0,274.0d0,305.0d0,335.0d0 /)


    pi=4.0d0*datan(1.0d0)
    twopi = pi * 2.0d0

    ! NOTE:  Year CE/AD = 0 is assumed to exist, and is equivalent to 1950 BP (-1950)
    ! subroutine orbpar() expects real-valued YearCE, but converts to YearBP for calculations
    select case (year_type)
    case ('CE', 'AD', 'ce', 'ad')
        YearCE = iyear
        YearBP = iyear - 1950
    case ('BP', 'bp')
        YearCE = iyear + 1950
        YearBP = iyear
    case default
        stop 'year_type'
    end select

!  Determine orbital parameters
    YEAR = dble(YearCE)
    CALL ORBPAR (YEAR, ECCEN,OBLIQ,OMEGVP) ! orbpar() expects YearCE input
    !write (*,'("YEAR, ECCEN,OBLIQ,OMEGVP:" 4f14.7)') YEAR, ECCEN,OBLIQ,OMEGVP 
    BSEMI  = dSQRT (1.0d0-ECCEN*ECCEN)
!  Vernal Equinox
    VEREQX = VERNAL (YearCE, EDAYzY)! NOTE:  made EDAYzY an argument
    !VEREQX = dmod(VERNAL (YearCE, EDAYzY), EDAYzy) ! NOTE:  made EDAYzY an argument
    CALL DtoYMDHM (VEREQX, JVEYR,JVEMON,JVEDAT,JVEHR,JVEMIN)
    TAofVE = - OMEGVP
    EAofVE = dATAN2 (dSIN(TAofVE)*BSEMI, dCOS(TAofVE)+ECCEN)
    MAofVE = EAofVE - ECCEN*dSIN(EAofVE)
    IF(MAofVE.lt.0.0d0)  MAofVE = MAofVE + TWOPI
    !write (*,'("VEREQX, TAofVE, EAofVE, MAofVE: ",4f14.7)') VEREQX, TAofVE, EAofVE, MAofVE
!  Perihelion
    KPERIH = 0
    PERIH1 = VEREQX - MAofVE*EDAYzY/TWOPI
    PERIH2 = PERIH1 + EDAYzY
    !write (*,'("EDAYzY, PERIH1, PERIH2: ",3f14.7)') EDAYzY, PERIH1, PERIH2
    Call DtoYMDHM (PERIH1, JPRYR,JPRMON,JPRDAT,JPRHR,JPRMIN)
    If (JPRYR /= IYEAR)  GoTo 210
    KPERIH = 1
    Call DtoYMDHM (PERIH2, JEXYR,JEXMON,JEXDAT,JEXHR,JEXMIN)
    If (JEXYR == IYEAR)  KPERIH = 2
    GoTo 220
210 Call DtoYMDHM (PERIH2, JPRYR,JPRMON,JPRDAT,JPRHR,JPRMIN)
    If (JPRYR == IYEAR)  KPERIH = 1
!  Aphelion
220 KAPHEL = 0
    APHEL1 = PERIH2 - 0.5d0*EDAYzY
    APHEL2 = PERIH2 + 0.5d0*EDAYzY
    Call DtoYMDHM (APHEL1, JAPYR,JAPMON,JAPDAT,JAPHR,JAPMIN)
    If (JAPYR /= IYEAR)  GoTo 230
    KAPHEL = 1
    If (KPERIH == 2)  GoTo 240
    Call DtoYMDHM (APHEL2, JEXYR,JEXMON,JEXDAT,JEXHR,JEXMIN)
    If (JEXYR == IYEAR)  KAPHEL = 2
    GoTo 240
230 Call DtoYMDHM (APHEL2, JAPYR,JAPMON,JAPDAT,JAPHR,JAPMIN)
    If (JAPYR == IYEAR)  KAPHEL = 1
!  Summer Solstice
240 TAofSS = TAofVE + 0.25d0*TWOPI
    EAofSS = dATan2 (dSin(TAofSS)*BSEMI, dCos(TAofSS)+ECCEN)
    MAofSS = EAofSS - ECCEN*dSin(EAofSS)
    If (MAofSS-MAofVE < 0)  MAofSS = MAofSS + TWOPI
    SUMSOL = VEREQX + (MAofSS-MAofVE)*EDAYzY/TWOPI
    Call DtoYMDHM (SUMSOL, JSSYR,JSSMON,JSSDAT,JSSHR,JSSMIN)
!  Autumnal Equinox
    TAofAE = TAofVE + 0.5d0*TWOPI
    EAofAE = dATan2 (dSin(TAofAE)*BSEMI, dCos(TAofAE)+ECCEN)
    MAofAE = EAofAE - ECCEN*dSin(EAofAE)
    If (MAofAE-MAofVE < 0)  MAofAE = MAofAE + TWOPI
    AUTEQX = VEREQX + (MAofAE-MAofVE)*EDAYzY/TWOPI
    Call DtoYMDHM (AUTEQX, JAEYR,JAEMON,JAEDAT,JAEHR,JAEMIN)
!  Winter Solstice
    TAofWS = TAofVE + 0.75d0*TWOPI
    EAofWS =dATan2 (dSin(TAofWS)*BSEMI, dCos(TAofWS)+ECCEN)
    MAofWS = EAofWS - ECCEN*dSin(EAofWS)
    If (MAofWS-MAofVE < 0)  MAofWS = MAofWS + TWOPI
    WINSOL = VEREQX + (MAofWS-MAofVE)*EDAYzY/TWOPI
    Call DtoYMDHM (WINSOL, JWSYR,JWSMON,JWSDAT,JWSHR,JWSMIN)

! vernal equinox and northern summer solstice days

    veqday = 0.0d0; ssday = 0.0d0
    if(isleap(YearCE)) then
        ndays_in_year = 366
        veqday = veqday + nd_leap(1) + nd_leap(2) + dble(JVEDAT) + dble(JVEHR)/24.0d0 + dble(JVEMIN)/1440.0d0
        ssday = ssday + nd_leap(1) + nd_leap(2) + nd_leap(3) + nd_leap(4) + nd_leap(5) &
            + dble(JSSDAT) + dble(JSSHR)/24.0d0 + dble(JSSMIN)/1440.0d0
        perihelion = dble(accumday_leap(JPRMON)) + dble(JPRDAT) + dble(JPRHR)/24.0d0 + dble(JPRMIN)/1440.0d0
        aphelion   = dble(accumday_leap(JAPMON)) + dble(JAPDAT) + dble(JAPHR)/24.0d0 + dble(JAPMIN)/1440.0d0
    else
        ndays_in_year = 365
        veqday = veqday + nd_noleap(1) + nd_noleap(2) + dble(JVEDAT) + dble(JVEHR)/24.0d0 + dble(JVEMIN)/1440.0d0
        ssday = ssday + nd_noleap(1) + nd_noleap(2) + nd_noleap(3) + nd_noleap(4) + nd_noleap(5) &
            + dble(JSSDAT) + dble(JSSHR)/24.0d0 + dble(JSSMIN)/1440.0d0
        perihelion = dble(accumday_noleap(JPRMON)) + dble(JPRDAT) + dble(JPRHR)/24.0d0 + dble(JPRMIN)/1440.0d0
        aphelion   = dble(accumday_noleap(JAPMON)) + dble(JAPDAT) + dble(JAPHR)/24.0d0 + dble(JAPMIN)/1440.0d0
    end if

end subroutine GISS_srevents

logical(4) function isleap(yearCE)
! is yearCE a leap year? -- no year zero or Gregorian-Julian adjustment

    ! NOTE:  Year CE/AD = 0 is assumed to exist, and is equivalent to 1950 BP (-1950)

    implicit none

    integer(4), intent(in)  :: yearCE

    isleap = .false.
    if (mod(yearCE, 4) .eq. 0) isleap = .true.
    if (mod(yearCE, 100) .eq. 0) isleap = .false.
    if (mod(yearCE, 400) .eq. 0) isleap = .true.

end function isleap

integer function before_leap(yearCE)
! number of years before next leap year

    ! NOTE:  Year CE/AD = 0 is assumbed to exist, and is equivalent to 1950 BP (-1950)

    implicit none

    integer(4), intent(in)  :: yearCE

    before_leap = 0
    do
        if(isleap(YearCE + before_leap)) exit
        before_leap = before_leap + 1
    end do

end function before_leap

real(8) function vernal (iyear, edayzy)
! subroutines based on
! SREVENTS.FOR    Solar EVENTS each year    2012/05/29
! https://data.giss.nasa.gov/ar5/SOLAR/SREVENTS.FOR downloaded 2017-09-12

!  For a given year, vernal calculates an approximate time of vernal
!  equinox in days measured from 2000 January 1, hour 0.
!
!  Vernal assumes that vernal equinoxes from one year to the next
!  are separated by exactly 365.2425 days, a tropical year
!  [Explanatory Supplement to the Astronomical Ephemeris].  If the
!  tropical year is 365.2422 days, as indicated by other references,
!  then the time of the vernal equinox will be off by 2.88 hours in
!  400 years.
!
!  Time of vernal equinox for year 2000 A.D. is March 20, 7:36 GMT
!  [NASA Reference Publication 1349, Oct. 1994].  Vernal assumes
!  that vernal equinox for year 2000 will be on March 20, 7:30, or
!  79.3125 days from 2000 January 1, hour 0.  Vernal equinoxes for
!  other years returned by vernal are also measured in days from
!  2000 January 1, hour 0.  79.3125 = 31 + 29 + 19 + 7.5/24.

    implicit none

    real(8), parameter      :: ve2000=79.3125d0
    real(8), intent(in)     :: edayzy
    integer(4), intent(in)  :: iyear

    vernal = ve2000 + dble((iyear-2000))*edayzy
    ! write (*,*) iyear, edayzy, ve2000, vernal, dble((iyear-2000))

end function vernal

SUBROUTINE DtoYMDHM (DAY, IYEAR,IMONTH,IDATE,IHOUR,IMINUT)
! subroutines based on
! SREVENTS.FOR    Solar EVENTS each year    2012/05/29
! https://data.giss.nasa.gov/ar5/SOLAR/SREVENTS.FOR downloaded 2017-09-12
!
!  DtoYMDHM receives DAY measured since 2000 January 1, hour 0 and
!  returns YEAR, MONTH, DATE, HOUR and MINUTE based on the Gregorian
!  calendar.

    implicit none

    real(8), intent(in)     :: day
    integer(4), intent(out) :: IYEAR,IMONTH,IDATE,IHOUR,IMINUT

    real(8)                 :: date

    call DtoYMD (DAY, IYEAR,IMONTH,DATE)
    IDATE  = int(DATE+1.)
    IMINUT = Nint ((DATE-IDATE+1)*24*60)
    IHOUR  = IMINUT / 60
    IMINUT = IMINUT - IHOUR*60
    Return

end subroutine DtoYMDHM

SUBROUTINE DtoYMD (DAY, IYEAR,IMONTH,DATE)
! subroutines based on
! SREVENTS.FOR    Solar EVENTS each year    2012/05/29
! https://data.giss.nasa.gov/ar5/SOLAR/SREVENTS.FOR downloaded 2017-09-12
!
!  For a given DAY measured from 2000 January 1, hour 0, determine
!  the IYEAR (A.D.), IMONTH and DATE (between 0. and 31.).

    implicit none

    real(8), intent(in)     :: day
    integer(4), intent(out) :: iyear, imonth
    real(8), intent(out)    :: date

    real(8), PARAMETER  :: JDAY4C = 365*400 + 97, &     !  number of days in 4 centuries
                            JDAY1C = 365*100 + 24, &    !  number of days in 1 century
                            JDAY4Y = 365*  4 +  1, &    !  number of days in 4 years
                            JDAY1Y = 365               !  number of days in 1 year

    integer(4)  :: m
    INTEGER(4)  :: JDSUMN(12),JDSUML(12), n4year, n1year
    INTEGER(8)  :: N4CENT, n1cent
    real(8)     :: day4c, day1c, day4y, day1y

    DATA JDSUMN /0,31,59, 90,120,151, 181,212,243, 273,304,334/
    DATA JDSUML /0,31,60, 91,121,152, 182,213,244, 274,305,335/

    N4CENT = FLOOR(DAY/JDAY4C)
    DAY4C  = DAY - N4CENT*JDAY4C
    N1CENT = int8((DAY4C-1)/JDAY1C)
    If (N1CENT > 0)  GoTo 10
!  First of every fourth century: 16??, 20??, 24??, etc.
    DAY1C  = DAY4C
    N4YEAR = int(DAY1C/JDAY4Y)
    DAY4Y  = DAY1C - N4YEAR*JDAY4Y
    N1YEAR = int((DAY4Y-1)/JDAY1Y)
    If (N1YEAR > 0)  GoTo 200
    GoTo 100
!  Second to fourth of every fourth century: 21??, 22??, 23??, etc.
10 DAY1C  = DAY4C - N1CENT*JDAY1C - 1
    N4YEAR = int((DAY1C+1)/JDAY4Y)
    If (N4YEAR > 0)  GoTo 20
!  First four years of every second to fourth century when there is
!  no leap year: 2100-2103, 2200-2203, 2300-2303, etc.
    DAY4Y  = DAY1C
    N1YEAR = int(DAY4Y/JDAY1Y)
    DAY1Y  = DAY4Y - N1YEAR*JDAY1Y
    GoTo 210
!  Subsequent four years of every second to fourth century when
!  there is a leap year: 2104-2107, 2108-2111 ... 2204-2207, etc.
20 DAY4Y  = DAY1C - N4YEAR*JDAY4Y + 1
    N1YEAR = int((DAY4Y-1)/JDAY1Y)
    If (N1YEAR > 0)  GoTo 200
!
!  Current year is a leap frog year
!
100 DAY1Y = DAY4Y
    Do 120 M=1,11
120 If (DAY1Y < JDSUML(M+1))  GoTo 130
!     M=12
130 IYEAR  = 2000 + int(N4CENT*400) + int(N1CENT*100) + N4YEAR*4 + N1YEAR
    IMONTH = M
    DATE   = DAY1Y - JDSUML(M)
    Return
!
!  Current year is not a leap frog year
!
200 DAY1Y  = DAY4Y - N1YEAR*JDAY1Y - 1
210 Do 220 M=1,11
220 If (DAY1Y < JDSUMN(M+1))  GoTo 230
!     M=12
230 IYEAR  = 2000 + int(N4CENT)*400 + int(N1CENT)*100 + N4YEAR*4 + N1YEAR
    IMONTH = M
    DATE   = DAY1Y - JDSUMN(M)
    Return

end subroutine DtoYMD

end module GISS_srevents_subs
