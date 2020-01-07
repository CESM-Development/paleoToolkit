;======================================================
load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"
;======================================================
; script from clay tabor
; to achieve desired affect for cam5 solar spectrum file

;=== examples of values for different periods
;cret values in w/m2
;MAA = 1357.18
;TUR = 1354.49
;CEN = 1353.90
;BER = 1348.65
;PT  = 1338.
;PETM = 1355.00
;1850 CAM4 default:  1361.27
;other examples: if want to increase scon by 10%, so ~ 1497.4
;sol_con_per = 1.1  ; amount to scale ; took this line out of code: it is NOT USED


begin


;=====================================================
;Load
;=====================================================
; copy default solar to infil directory and chmod u+w (to give write permissions)

 infil =  "/glade/p/cgd/ccr/people/shields/misc/SOLAR_SPECTRAL_Lean_1610-2008_annual_c090324.nc"
 outfil = "./solar_scon_<casename>.<date>.nc"


 system("/bin/rm -f " + outfil)
 system("cp " + infil + " " + outfil) ; copy file to output

 sol_in = addfile(infil,"r")
 solar_out = addfile(outfil,"w")
 band_in = sol_in->band_width
 wvl_val = sol_in->wvl
 date_in = sol_in->date
 wvl_in = wvl_val / 1.e9

 ssi_in = sol_in->ssi
 ssi_out = ssi_in 


;=====================================================
;Variables
;=====================================================

 c = 300000000d                   ;speed of light
 h = 6.625e-34                    ;planck's constant
 k = 1.37e-23                     ;boltzmann's constant
 sig = 5.67e-8                    ;stefan-Boltzmann constant
 d = 149597870700d                ;mean earth sun distance
 pi = acos(-1.)                   ;3.14
 rp = 696342000d                  ;radius of photosphere
 So_adj = <value>                 ;desired solar constant 

;=====================================================
;Calculate new solar constants
;=====================================================

 v = c / wvl_in                  ;frequency
 int_w_in = 0. * ssi_in

 do i = 0,dimsizes(ssi_in(:,0))-1
   int_w_in(i,:) = ssi_in(i,:) * band_in(:)       ;band width weighted irradiance
 end do

 So_mod = dim_sum_n(int_w_in,1) / 1000. ;modern solar constants
 So_new = So_mod * 0. + So_adj          ;new solar constant
 So_pct_chng = So_new / So_mod          ;percent difference in solar constant
 Lo_mod = So_mod * 4. * pi * d^2        ;modern luminosity
 Lo_new = So_new * 4. * pi * d^2        ;new luminosity
 Fd_mod = Lo_mod / (4. * pi * rp^2)     ;modern flux density
 Fd_new = Lo_new / (4. * pi * rp^2)     ;new flux density
 T_mod = (Fd_mod / sig) ^ 0.25          ;modern bb T
 T_new = (Fd_new / sig) ^ 0.25          ;new bb T

 Bv_mod = 0. * ssi_in 
 Bv_new = 0. * ssi_in

 do i = 0,dimsizes(ssi_in(:,0))-1
   do j = 0,dimsizes(ssi_in(0,:))-1
     Bv_mod(i,j) = 1.e12 * ((2. * h * v(j)^3)) / (c^2 * (exp((h * v(j)) / (k * T_mod(i)))) - 1.) ;modern bb intensity
     Bv_new(i,j) = 1.e12 * ((2. * h * v(j)^3)) / (c^2 * (exp((h * v(j)) / (k * T_new(i)))) - 1.) ;new bb intensity
   end do
 end do

 pct_chng = Bv_new / Bv_mod        ;% difference between bb spectra
 ssi_iout = pct_chng * ssi_in      ;scaling new irradiance by % bb difference
 ;print("pct_chng = " + pct_chng)

 int_w_iout = 0. * ssi_in

 do i = 0,dimsizes(ssi_in(:,0))-1
   int_w_iout(i,:) = ssi_iout(i,:) * band_in(:) ;new weighted irradiance
 end do

 So_iout = dim_sum_n(int_w_iout,1) / 1000.     ;new solar constants
 scl_corr = So_adj / So_iout                   ;solar constant difference between conversion and desired
 ;print("scl_corr = " + scl_corr)

 ssi_fout = 0. * ssi_in
 ssi_sout = 0. * ssi_in

 do i = 0,dimsizes(ssi_in(:,0))-1
   ssi_fout(i,:) = ssi_iout(i,:) * scl_corr(i)       ;scale remaining solar constant difference
   ssi_sout(i,:) = ssi_in(i,:) * So_pct_chng(i)      ;uniformally scaled irradiance (for comparison)
 end do

 int_w_fout = 0. * ssi_in

 do i = 0,dimsizes(ssi_in(:,0))-1
   int_w_fout(i,:) = ssi_fout(i,:) * band_in(:)       ;scale corrected new weighted irradiance
 end do

 So_fout = dim_sum_n(int_w_fout,1) / 1000.    ;scaled corrected new solar constants

 ssi_out = ssi_fout

;=====================================================
;Check / compare
;=====================================================

 print(So_iout(0))
 print(So_fout(0))

;=====================================================
;Create netCDF file
;=====================================================

; solar_out->ssi_sout = ssi_sout
; solar_out->ssi_iout = ssi_iout
; solar_out->ssi_fout = ssi_fout
; solar_out->Bv_mod = Bv_mod
; solar_out->Bv_new = Bv_new
  solar_out->ssi = ssi_out
  solar_out->band_width = band_in
  solar_out->date = date_in
  solar_out@paleo_comment = outfil +  " data at year 1850 has been overwritten for paleo period"
  solar_out@scon = " data is appropriate for " + So_adj + " W/m2"

end
