program giss_orbpar
! calls GISS subroutine ORBPAR(YEAR, ECCEN,OBLIQ,OMEGVP) to compute orbital parameters
! using Berger (1978, JAS 35:2362-2367) algorithm and tables
! https://data.giss.nasa.gov/ar5/SOLAR/ORBPAR.FOR downloaded 2017-09-04 17:17
! also retrievable from https://web.archive.org/web/20150920211936/http://data.giss.nasa.gov/ar5/solar.html
    
use GISS_orbpar_subs

implicit none

! ages -- input argument is YearBP.  YearCE is calculated for reference
character(2)    :: year_type        ! year type (AD, CE, BP)
real(8)         :: yearBP           ! Year BP 1950
real(8)         :: yearCE           ! Year CE (input argument to GISS_orbpars())

! subroutine get_GISS_orbpars output arguments
real(8)         :: eccen            ! eccentricity of orbital ellipse
real(8)         :: obliq_deg        ! obliquity (degrees)
real(8)         :: perih_deg        ! longitude of perihelion (degrees)
real(8)         :: precc            ! climatological precession parameter = eccen * sin(omegvp)

! indices
integer(4)      :: begyr            ! beginning year
integer(4)      :: endyr            ! ending year
integer(4)      :: yrstep           ! year step size
integer(4)      :: n                ! year index

character(2056) :: outpath, outfile ! output file name

outpath = "/Projects/Calendar/PaleoCalAdjust/data/GISS_orbital/"   ! Windows path
!outpath = "/Users/bartlein/Projects/Calendar/PaleoCalAdjust/data/GISS_orbital/"    ! Mac path
outfile="orb_elt_150ka_1kyr.csv"
write (*,'(a)') trim(outpath)//trim(outfile)
open (1, file=trim(outpath)//trim(outfile))
write (1,'(a)') "YearCE, YearBP, Eccen, Obliq_deg, Perih_deg, ClimPrecc"

! past ages are negative, e.g. 21 ka = 21,000 cal yr BP = -21000, and 1950 CE = 0 cal yr BP = 0 here

year_type = "BP"
begyr=-150000; endyr=0; yrstep=1000

do n=begyr,endyr,yrstep
    ! NOTE:  Year CE/AD = 0 is assumed to exist, and is equivalent to 1950 BP (-1950)
    ! subroutine orbpar() expects real-valued YearCE, but converts to YearBP for calculations
    select case (year_type)
    case ('CE', 'AD', 'ce', 'ad')
        YearCE = dble(n)
        YearBP = dble(n) - 1950.0d0
    case ('BP', 'bp')
        YearCE = dble(n) + 1950.0d0
        YearBP = dble(n)
    case default
        stop 'year_type'
    end select

    call GISS_orbpars(year_type, yearBP, eccen, obliq_deg, perih_deg, precc)

    write (1,'(f10.1,", ",f10.1,4(", ",f17.12))') yearCE, yearBP, eccen, obliq_deg, perih_deg, precc

end do

write (*,'(a)') "Done (GISS_orbpar_driver)"

end program giss_orbpar


