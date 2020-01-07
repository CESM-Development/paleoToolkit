program srevents
! f90 version of SREVENTS.FOR driver program
! subroutines based on GISS program SREVENTS.FOR -- Solar EVENTS each year    2003/10/22
! downloaded 2016-10-07 16:20
! also retrieveable from https://web.archive.org/web/20150920211936/http://data.giss.nasa.gov/ar5/solar.html
    
! NOTE:  Year CE/AD = 0 is assumed to exist, and is equivalent to 1950 BP (-1950)
    
use GISS_orbpar_subs
use GISS_srevents_subs
      
implicit none

real(8), parameter      :: EDAYzY=365.2425d0 ! 365.24219876d0 ) ! 
integer(4), parameter   :: nm=12

integer(4)      :: iyminBP, iymaxBP, iyinc, iyminCE, iymaxCE, iyearCE, ndays, i4yrs
real(8)         :: veqday, veqday2, ssday, ssday2, perihelion, aphelion
character(2)    :: year_type
character(2056) :: outpath
character(64)   :: sreventsfile, VPAdayfile

outpath = "/Projects/Calendar/PaleoCalAdjust/data/GISS_orbital/" ! Windows path
!outpath = "/Users/bartlein/Projects/Calendar/PaleoCalAdjust/data/GISS_orbital/" ! Mac Path
sreventsfile="srevents_150ka_1kyr.txt" 
VPAdayfile="VPAday_150ka_1kyr.csv" 

open (1, file=trim(outpath)//trim(sreventsfile))
open (2, file=trim(outpath)//trim(VPAdayfile))

iyminBP = -150000 ! 
iymaxBP =  0 
iyinc = 1000 

year_type = 'CE'
iyminCE = iyminBP + 1950
iymaxCE = iymaxBP + 1950

!  Write header information
WRITE (1,920) EDAYzY
write (2,'(a)') "YearBP, YearCE, ndays, before_leap, veq_day, veq_day2, ss_day, ss_day2, perihelion_day, aphelion_day"

! loop over years

i4yrs = before_leap(iyminCE + iyinc)
do iyearCE = iyminCE,iymaxCE,iyinc
    
    call GISS_srevents(year_type, iyearCE, edayzy, veqday, ssday, perihelion, aphelion, ndays)

    if(KPERIH==1 .and. KAPHEL==1)  then
        write (1,927) iyearCE, &
            jvemon,jvedat,jvehr,jvemin, jssmon,jssdat,jsshr,jssmin, &
            jaemon,jaedat,jaehr,jaemin, jwsmon,jwsdat,jwshr,jwsmin, &
            jprmon,jprdat,jprhr,jprmin, japmon,japdat,japhr,japmin 
    end if

    if(KPERIH==0) then 
        write (1,928) iyearCE, &
            jvemon,jvedat,jvehr,jvemin, jssmon,jssdat,jsshr,jssmin, &
            jaemon,jaedat,jaehr,jaemin, jwsmon,jwsdat,jwshr,jwsmin, &
                                        japmon,japdat,japhr,japmin
    end if

    if(KPERIH==2 .or. KAPHEL==2) then 
        write (1,927) iyearCE, &
            jvemon,jvedat,jvehr,jvemin, jssmon,jssdat,jsshr,jssmin, &
            jaemon,jaedat,jaehr,jaemin, jwsmon,jwsdat,jwshr,jwsmin, &
            jprmon,jprdat,jprhr,jprmin
    end if
    
    veqday2 =  80.0d0 + (dble(i4yrs) * 0.25d0)
    ssday2 = 170.5d0 + (dble(i4yrs) * 0.25d0)
    i4yrs = i4yrs + 1
    if (i4yrs.gt.3) i4yrs = 0

    write (2,'(i8,", ",i8,2(", ",i4),4(", ",f10.6), 2(", ",f12.6))') & 
        iyearCE-1950, iyearCE, ndays, i4yrs, veqday, veqday2, ssday, ssday2, perihelion, aphelion

end do

GO TO 999

920 format ('solar events' / &
    '  Tropical Year = ', f13.9,' (days)' // &
    '         Vernal      Summer      Autumnal     Winter      '  / &
    '   YrCE     Equinox    Solstice     Equinox     Solstice   Perihelion    Aphelion ' / &
    '   ----  ----------  ----------  ----------   ----------   ----------   ----------')
927 format (i7,3(i3,'/',i2.2,i3,':',i2.2),3(i4,'/',i2.2,i3,':',i2.2))
928 format (i7,3(i3,'/',i2.2,i3,':',i2.2),  i4,'/',i2.2,i3,':',i2.2, &
            13x, i4,'/',i2.2,i3,':',i2.2)
! 981 format (a / a)
999 continue
    
write (*,'(a)') "Done (GISS_srevents_driver)"
    
end program srevents 


      