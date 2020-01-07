program demo_02_adjust_1yr
! demonstrates (paleo) calendar adjustment for one year (of monthly values)
! 6 ka example using a 365-day "noleaps" calendar
    
use pseudo_daily_interp_subs    ! pseudo-daily interpolation module
use calendar_effects_subs       ! calendar effects subroutines

implicit none

integer, parameter  :: nm=12, nd=365, ny=1    ! number of months in year, number of days in year, number of years

! month length, beginning and ending days (from month_length.f90), present (1950 CE) and 6 ka
real(8)         :: rmonlen_00(nm)
real(8)         :: rmonlen_06(nm), rmonbeg_06(nm), rmonend_06(nm)

! number of days in year
integer(4)      :: ndays(ny)

! data, "raw" (i.e. on 0ka calendar) and adjusted to 6 ka calendar
real(8)         :: xm_0ka_cal(nm), xm_adj_cal(nm)

! weighted annual averages
real(8)         :: xann_0ka_cal, xann_adj_cal

! pseudo-daily interpolated data
real(8)         :: xdh(nd),vfill

! smoothing parameters for multi-year pseudo-daily interpolation
integer(4)              :: nw_tmp=21, nsw_tmp=20    ! smoothing parameters
logical                 :: smooth=.false., restore=.true.
logical                 :: no_negatives = .false. ! restrict pseudo-daily interpolated values to positive values?      

! month-length data, 0 ka
data rmonlen_00 / 31.0d0,  28.0d0,  31.0d0,  30.0d0,  31.0d0,  30.0d0,  31.0d0,  31.0d0,  30.0d0,  31.0d0,  30.0d0,  31.0d0 /

! paleo month-length data (e.g. 6 ka)

data rmonlen_06 / 32.48411855d0,  29.53995521d0,  32.47169095d0,  30.80519959d0,  30.97070708d0,  29.16304027d0,  & 
    29.54741870d0,  29.34677894d0,  28.62575603d0,  30.18400184d0,  30.01008258d0,  31.85125027d0 / 
data rmonbeg_06 / -4.05417526d0,  28.40375842d0,  57.92559617d0,  90.40319079d0, 121.23431371d0, 152.22542876d0, & 
    181.38183671d0, 210.90212219d0, 240.22885383d0, 268.86102212d0, 299.07135244d0, 329.10071578d0 / 
data rmonend_06 / 28.40375842d0,  57.92559617d0,  90.40319079d0, 121.23431371d0, 152.22542876d0, 181.38183671d0, & 
    210.90212219d0, 240.22885383d0, 268.86102212d0, 299.07135244d0, 329.10071578d0, 360.94582474d0 /

! data to adjust -- example tas values c65r36 from TraCE-21ka data
data xm_0ka_cal /-7.68829d0, -2.08282d0, -2.47379d0, 0.732910d0, 3.59814d0, 10.4207d0, &
    15.6256d0, 13.9701d0, 10.3976d0, 7.65222d0, -1.86514d0, -4.61826d0 / 
vfill = 1.0e32
ndays(1) = nd

write (out_unit,'(a)') "Paleo calendar adjustment of monthly means:"

write (out_unit,'(/a)') "Input data (monthly mean tas on 0 ka 365-day noleaps calendar):"
write (out_unit,'("   xm_0ka_cal: ",12f9.3)') xm_0ka_cal

! weighted annual mean of monthly input
call ann_wmean(nm, xm_0ka_cal, rmonlen_00, xann_0ka_cal)
write (out_unit,'("Weighted annual mean of input data: ", f12.6)') xann_0ka_cal

! interpolate monthly data to pseudo-daily values
call mon_to_day_ts(nm, int(rmonlen_00), xm_0ka_cal, vfill, no_negatives, smooth, restore, nd, nw_tmp, nsw_tmp, xdh)
write (*,*) nd,ndays(1)

! reaggregate pseudo-daily values to new calendar
!call adjusted_monthly_means(nm, nd, xdh, rmonbeg_06, rmonend_06, xm_adj_cal)
call day_to_mon_ts(ny, ndays, rmonbeg_06, rmonend_06, nd, xdh, vfill, xm_adj_cal)
write (out_unit,'(/a)') "Paleo month-length adjusted values on 6ka calendar:"
write (out_unit,'("   xm_adj_cal: ",12f9.3)') xm_adj_cal

! weighted annual mean of adjusted data
call ann_wmean(nm, xm_adj_cal, rmonlen_06, xann_adj_cal)
write (out_unit,'("Weighted annual mean of adjusted data: ", f12.6)') xann_adj_cal

end program