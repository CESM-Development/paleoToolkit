program demo01_pseudo_daily_interp
! compares linear and harmonic (Epstein, 1991) pseudo-daily interpolation

use pseudo_daily_interp_subs

implicit none

integer, parameter  :: nm=12, nd=365
real(8)             :: xm(nm)
real(8)             :: xd_linear(nd), xd_harmonic(nd)
real(8)             :: xm_linear(nm), xm_harmonic(nm)
real(8)             :: xann, xann_linear, xann_harmonic
integer(4)          :: imonlen(nm)

logical             :: no_negatives
logical             :: debug_write=.true.

debug_unit=10; out_unit=1
if (debug_write) open (10, file="/Projects/Calendar/data/work01/debug_pseudo_daily.csv")
open (out_unit, file="/Projects/Calendar/PaleoCalAdjust/data/figure_data/pseudo_daily_plots/Fig15.dat")

write (out_unit,'(a)') "Comparison of linear and harmonic pseudo-daily interpolation"

! month length 
data imonlen /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

! Example 1: linear and harmonic (Epstein, 1991) interpolation
! example input data (e.g. tas)
xm = (/ -6.77, -4.33, 1.52, 8.35, 14.74, 20.29, 23.22, 22.17, 16.94, 10.35, 3.08, -4.29/)

! write out the monthly "control" data
write (out_unit,'(/a)') "Example 1, Madison CFSR tas (Fig. 15a):"
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)

! weighted annual mean of monthly input
call ann_wmean(nm, xm, dble(imonlen), xann)
write (out_unit,'("Weighted annual mean of input data: ", f12.6)') xann

! linear interpolation (not mean preserving)
write (out_unit,'(/a)') "Linear interpolation (not mean preserving):"
call dayinterp(nm, nd, imonlen, xm, xd_linear)
call monmean(nm, nd, imonlen, xd_linear, xm_linear)
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)
write (out_unit,'("    xm_linear: ",12f9.3)') xm_linear(1:nm)
write (out_unit,'("   difference: ",12f9.3)') xm_linear-xm
call ann_mean(nd, xd_linear, xann_linear)
write (out_unit,'("Annual mean of interpolated data: ", f12.6, "    Difference:", f12.6)') xann_linear, xann-xann_linear

! mean-preserving (harmonic) interpolation (Epstein, 1991) 
write (out_unit,'(/a)') "Mean-preserving (harmonic) interpolation (Epstein, 1991):"
no_negatives = .false.
call  hdaily(nm, nd, xm, imonlen, no_negatives, xd_harmonic)
call monmean(nm, nd, imonlen, xd_harmonic, xm_harmonic)
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)
write (out_unit,'("  xm_harmonic: ",12f9.3)') xm_harmonic(1:nm)
write (out_unit,'("   difference: ",12f9.3)') xm_harmonic-xm
call ann_mean(nd, xd_harmonic, xann_harmonic)
write (out_unit,'("Annual mean of interpolated data: ", f12.6, "    Difference:", f12.6)') xann_harmonic, xann-xann_harmonic

! Example 2: linear and harmonic (Epstein, 1991) interpolation
! example input data (e.g. tas)
xm = (/ 25.320, 24.930, 22.630, 19.730, 16.170, 13.000, 12.380, 14.040, 17.360, 19.780, 22.520, 24.040/)

! write out the monthly "control" data
write (out_unit,'(/a)') "Example 2, Australia CFSR tas (Fig. 15b):"
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)

! weighted annual mean of monthly input
call ann_wmean(nm, xm, dble(imonlen), xann)
write (out_unit,'("Weighted annual mean of input data: ", f12.6)') xann

! linear interpolation (not mean preserving)
write (out_unit,'(/a)') "Linear interpolation (not mean preserving):"
call dayinterp(nm, nd, imonlen, xm, xd_linear)
call monmean(nm, nd, imonlen, xd_linear, xm_linear)
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)
write (out_unit,'("    xm_linear: ",12f9.3)') xm_linear(1:nm)
write (out_unit,'("   difference: ",12f9.3)') xm_linear-xm
call ann_mean(nd, xd_linear, xann_linear)
write (out_unit,'("Annual mean of interpolated data: ", f12.6, "    Difference:", f12.6)') xann_linear, xann-xann_linear

! mean-preserving (harmonic) interpolation (Epstein, 1991) 
write (out_unit,'(/a)') "Mean-preserving (harmonic) interpolation (Epstein, 1991):"
no_negatives = .false.
call  hdaily(nm, nd, xm, imonlen, no_negatives, xd_harmonic)
call monmean(nm, nd, imonlen, xd_harmonic, xm_harmonic)
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)
write (out_unit,'("  xm_harmonic: ",12f9.3)') xm_harmonic(1:nm)
write (out_unit,'("   difference: ",12f9.3)') xm_harmonic-xm
call ann_mean(nd, xd_harmonic, xann_harmonic)
write (out_unit,'("Annual mean of interpolated data: ", f12.6, "    Difference:", f12.6)') xann_harmonic, xann-xann_harmonic

! Example 3: linear and harmonic (Epstein, 1991) interpolation
! example input data (e.g. precip (rate))
xm = (/ 0.900, 0.980, 1.570, 2.490, 2.710, 3.350, 3.070, 3.160, 2.700, 2.070, 1.960, 1.170 /)

! write out the monthly "control" data
write (out_unit,'(/a)') "Example 3, Madison CMAP precip (Fig. 15c):"
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)

! weighted annual mean of monthly input
call ann_wmean(nm, xm, dble(imonlen), xann)
write (out_unit,'("Weighted annual mean of input data: ", f12.6)') xann

! linear interpolation (not mean preserving)
write (out_unit,'(/a)') "Linear interpolation (not mean preserving):"
call dayinterp(nm, nd, imonlen, xm, xd_linear)
call monmean(nm, nd, imonlen, xd_linear, xm_linear)
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)
write (out_unit,'("    xm_linear: ",12f9.3)') xm_linear(1:nm)
write (out_unit,'("   difference: ",12f9.3)') xm_linear-xm
call ann_mean(nd, xd_linear, xann_linear)
write (out_unit,'("Annual mean of interpolated data: ", f12.6, "    Difference:", f12.6)') xann_linear, xann-xann_linear

! mean-preserving (harmonic) interpolation (Epstein, 1991) 
write (out_unit,'(/a)') "Mean-preserving (harmonic) interpolation (Epstein, 1991):"
no_negatives = .false.
call  hdaily(nm, nd, xm, imonlen, no_negatives, xd_harmonic)
call monmean(nm, nd, imonlen, xd_harmonic, xm_harmonic)
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)
write (out_unit,'("  xm_harmonic: ",12f9.3)') xm_harmonic(1:nm)
write (out_unit,'("   difference: ",12f9.3)') xm_harmonic-xm
call ann_mean(nd, xd_harmonic, xann_harmonic)
write (out_unit,'("Annual mean of interpolated data: ", f12.6, "    Difference:", f12.6)') xann_harmonic, xann-xann_harmonic

! Example 4: linear and harmonic (Epstein, 1991) interpolation
! pathological precipitation case
xm = (/ 0.21, 0.07, 0.03, 0.09, 2.53, 8.23, 2.0, 1.09, 1.73, 2.23, 2.48, 1.02/)

! write out the monthly "control" data
write (out_unit,'(/a)') "Example 4, Indian Ocean CMAP precip (pathological precipitation case) (Fig. 15d):"
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)

! weighted annual mean of monthly input
call ann_wmean(nm, xm, dble(imonlen), xann)
write (out_unit,'("Weighted annual mean of input data: ", f12.6)') xann

! linear interpolation (not mean preserving)
write (out_unit,'(/a)') "Linear interpolation (not mean preserving):"
call dayinterp(nm, nd, imonlen, xm, xd_linear)
call monmean(nm, nd, imonlen, xd_linear, xm_linear)
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)
write (out_unit,'("    xm_linear: ",12f9.3)') xm_linear(1:nm)
write (out_unit,'("   difference: ",12f9.3)') xm_linear-xm
call ann_mean(nd, xd_linear, xann_linear)
write (out_unit,'("Annual mean of interpolated data: ", f12.6, "    Difference:", f12.6)') xann_linear, xann-xann_linear

! harmonic (Epstein, 1991) interpolation mean preserving)
write (out_unit,'(/a)') "Harmonic (Epstein, 1991) interpolation (mean preserving):"
no_negatives = .true.
call  hdaily(nm, nd, xm, imonlen, no_negatives, xd_harmonic)
call monmean(nm, nd, imonlen, xd_harmonic, xm_harmonic)
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)
write (out_unit,'("  xm_harmonic: ",12f9.3)') xm_harmonic(1:nm)
write (out_unit,'("   difference: ",12f9.3)') xm_harmonic-xm
call ann_mean(nd, xd_harmonic, xann_harmonic)
write (out_unit,'("Annual mean of interpolated data: ", f12.6, "    Difference:", f12.6)') xann_harmonic, xann-xann_harmonic


! Example 5: linear and harmonic (Epstein, 1991) interpolation
! example input data with zeros (e.g. pre)
xm = (/ 44., 60.,  182.,  103.,  7.,  0.,  0.,  0.,  6., 36., 78., 66./)

! write out the monthly "control" data
write (out_unit,'(/a)') "Example 5,  pre with zeros (Not in Fig. 15):"
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)

! weighted annual mean of monthly input
call ann_wmean(nm, xm, dble(imonlen), xann)
write (out_unit,'("Weighted annual mean of input data: ", f12.6)') xann

! linear interpolation (not mean preserving)
write (out_unit,'(/a)') "Linear interpolation (not mean preserving):"
call dayinterp(nm, nd, imonlen, xm, xd_linear)
call monmean(nm, nd, imonlen, xd_linear, xm_linear)
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)
write (out_unit,'("    xm_linear: ",12f9.3)') xm_linear(1:nm)
write (out_unit,'("   difference: ",12f9.3)') xm_linear-xm
call ann_mean(nd, xd_linear, xann_linear)
write (out_unit,'("Annual mean of interpolated data: ", f12.6, "    Difference:", f12.6)') xann_linear, xann-xann_linear

! harmonic (Epstein, 1991) interpolation mean preserving)
write (out_unit,'(/a)') "Harmonic (Epstein, 1991) interpolation (mean preserving):"
no_negatives = .true.
call  hdaily(nm, nd, xm, imonlen, no_negatives, xd_harmonic)
call monmean(nm, nd, imonlen, xd_harmonic, xm_harmonic)
write (out_unit,'("           xm: ",12f9.3)') xm(1:nm)
write (out_unit,'("  xm_harmonic: ",12f9.3)') xm_harmonic(1:nm)
write (out_unit,'("   difference: ",12f9.3)') xm_harmonic-xm
call ann_mean(nd, xd_harmonic, xann_harmonic)
write (out_unit,'("Annual mean of interpolated data: ", f12.6, "    Difference:", f12.6)') xann_harmonic, xann-xann_harmonic

end program demo01_pseudo_daily_interp
    
