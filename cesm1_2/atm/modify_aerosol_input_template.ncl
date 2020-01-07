
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"



function zero_var(Fil, varName)
local X
begin
 X = Fil->$varName$
 X = (/0.0 * X/)
 return(X)
end



function get_vars_from_list(fvars, listvars)
local flags, i, varlist
begin
 flags = new( dimsizes(fvars), integer )
 do i = 0,dimsizes(fvars)-1 ;; fvars probably shorter than listvars
   if any( fvars(i).eq.listvars ) then
     flags(i) = 1
   end if
 end do
 if .not.all(ismissing(flags)) then
   varlist = fvars( ind( flags.eq.1 ) )
   print("CHANGING: " + varlist)
   varlist2 = fvars( ind( ismissing(flags) ) )
   print("NOT CHANGING: " + varlist2)
 else
   varlist = new( 1, typeof(fvars) )
   print("WARNING!! NO VARIABLES FOUND WHICH MATCH!")
   varlist2 = fvars( ind( ismissing(flags) ) )
   print("NOT CHANGING: " + varlist2)
 end if
 return(varlist)
end



function my_wrap(Xtmp,my_rank,my_z)
local tz
begin
 if (my_rank.eq.4) then
 tz = dim_avg_n_Wrap( Xtmp(time|:, $my_z$|:, lat|:, lon|:), (/3/) )
 end if
 if (my_rank.eq.3) then
 tz = dim_avg_n_Wrap( Xtmp(time|:, lat|:, lon|:), (/2/) )
 end if
 return(tz)
end



function my_sym(tz,my_rank,my_z)
local tz_sym
begin

 if (my_rank.eq.4) then

do i = 0,dimsizes(tz(:,0,0))-1
  do j = 0,dimsizes(tz(0,:,0))-1
    do k = 0,dimsizes(tz(0,0,:))-1
      if (ismissing(tz(i,j,k))) then
        if ((i+6.le.11)) then
        month = i+6
        else
        month = i-6
        end if
        tmp = tz(month,j,95-k)
        tz(i,j,k) = tmp
      end if
    end do
  end do
end do

 end if

 if (my_rank.eq.3) then

  do i = 0,dimsizes(tz(:,0))-1
    do j = 0,dimsizes(tz(0,:))-1
      if (ismissing(tz(i,j))) then
        if ((i+6.le.11)) then
        month = i+6
        else
        month = i-6
        end if
        tmp = tz(month,95-j)
        tz(i,j) = tmp
      end if
    end do
  end do

 end if

 tz_sym = tz

 if (my_rank.eq.4) then
do i = 0,5
do j = 0,47
 tz_sym(time|i,$my_z$|:,lat|j) = (/ 0.5*( tz(time|i,$my_z$|:,lat|j)+tz(time|i+6,$my_z$|:,lat|95-j) ) /)
 tz_sym(time|i,$my_z$|:,lat|j+48) = (/ 0.5*( tz(time|i,$my_z$|:,lat|48+j)+tz(time|i+6,$my_z$|:,lat|47-j) ) /)
end do
end do
do i = 6,11
do j = 0,47
 tz_sym(time|i,$my_z$|:,lat|j) = (/ 0.5*( tz(time|i,$my_z$|:,lat|j)+tz(time|i-6,$my_z$|:,lat|95-j) ) /)
 tz_sym(time|i,$my_z$|:,lat|j+48) = (/ 0.5*( tz(time|i,$my_z$|:,lat|48+j)+tz(time|i-6,$my_z$|:,lat|47-j) ) /)
end do
end do
 end if

 if (my_rank.eq.3) then
do i = 0,5
do j = 0,47
 tz_sym(time|i,lat|j) = (/ 0.5*( tz(time|i,lat|j)+tz(time|i+6,lat|95-j) ) /)
 tz_sym(time|i,lat|j+48) = (/ 0.5*( tz(time|i,lat|48+j)+tz(time|i+6,lat|47-j) ) /)
end do
end do
do i = 6,11
do j = 0,47
 tz_sym(time|i,lat|j) = (/ 0.5*( tz(time|i,lat|j)+tz(time|i-6,lat|95-j) ) /)
 tz_sym(time|i,lat|j+48) = (/ 0.5*( tz(time|i,lat|48+j)+tz(time|i-6,lat|47-j) ) /)
end do
end do
 end if

 return(tz_sym)
end


function paleotize(X,my_z,my_rank,my_type)
local Xtmp, lat, coslat, tz, tz_sym, out
begin

 Lfil = addfile("/glade/p/cesm/cseg/inputdata/share/domains/domain.lnd.fv1.9x2.5_gx1v6.090206.nc", "r")
 LAND = Lfil->frac

 Pfil = addfile("<path/domain.lnd_file>", "r")
 PAND = Pfil->frac

;; printVarSummary(LAND)
;; printVarSummary(X)

 dimNames = getvardims(X)
;; print(dimNames)
 lat = X&lat
 coslat = lat
 coslat = (/ cos( lat*4.*atan(1.)/180.) /)

 if (my_rank.eq.4) then
 my_dims1=(/0,1,2/)
 my_dims2=(/2,3/)
 end if
 if (my_rank.eq.3) then
 my_dims1=(/0,1/)
 my_dims2=(/1,2/)
 end if

 ;; ALL:

 X@_FillValue = getVarFillValue(X)
;; print(X@_FillValue)

 land = conform(X, LAND, my_dims2)
 pand = conform(X, PAND, my_dims2)
 my_mask = (/ where(pand.gt.0, 1, 0) /)

 ;; OCN:
 Xtmp = X
 Xtmp = (/ where(land.gt.0, X@_FillValue, Xtmp) /)
 tz=my_wrap(Xtmp,my_rank,my_z)
 delete(Xtmp)
 tz_sym = my_sym(tz,my_rank,my_z)
 v_ocn = (/conform( X, tz_sym, my_dims1 ) /) * (1-my_mask)

 ;; LND:
 Xtmp = X
 Xtmp = (/ where(land.le.0, X@_FillValue, Xtmp) /)
 tz=my_wrap(Xtmp,my_rank,my_z)
 delete(Xtmp)
 tz_sym = my_sym(tz,my_rank,my_z)

if (my_rank.eq.4) then
 tz_sym(:,:,15) = 0.25 * (3.*tz_sym(:,:,14)+1.*tz_sym(:,:,18))
 tz_sym(:,:,16) = 0.25 * (2.*tz_sym(:,:,14)+2.*tz_sym(:,:,18))
 tz_sym(:,:,17) = 0.25 * (1.*tz_sym(:,:,14)+3.*tz_sym(:,:,18))
 tz_sym(:,:,78) = 0.25 * (3.*tz_sym(:,:,77)+1.*tz_sym(:,:,81))
 tz_sym(:,:,79) = 0.25 * (2.*tz_sym(:,:,77)+2.*tz_sym(:,:,81))
 tz_sym(:,:,80) = 0.25 * (1.*tz_sym(:,:,77)+3.*tz_sym(:,:,81)) 
end if
if (my_rank.eq.3) then
 tz_sym(:,15) = 0.25 * (3.*tz_sym(:,14)+1.*tz_sym(:,18))
 tz_sym(:,16) = 0.25 * (2.*tz_sym(:,14)+2.*tz_sym(:,18))
 tz_sym(:,17) = 0.25 * (1.*tz_sym(:,14)+3.*tz_sym(:,18))
 tz_sym(:,78) = 0.25 * (3.*tz_sym(:,77)+1.*tz_sym(:,81))
 tz_sym(:,79) = 0.25 * (2.*tz_sym(:,77)+2.*tz_sym(:,81))
 tz_sym(:,80) = 0.25 * (1.*tz_sym(:,77)+3.*tz_sym(:,81))
end if

 v_lnd = (/conform( X, tz_sym, my_dims1 ) /) * (my_mask)

 out=v_lnd+v_ocn
 out@info = "time-mean zonal-mean symmetrized about equator, separate values applied to land and ocean"

 return(out)
end


function set_PS(X)
local out, Fp
begin

  Fp = addfile("<path/cam_initial_file>","r")

  ptmp = Fp->$"PS"$
  ptmp2=doubletofloat(rm_single_dims(ptmp))
;;  printVarSummary(X)
;;  printVarSummary(ptmp)
;;  printVarSummary(ptmp2)
  out=X
  out = conform(X,ptmp2,(/1,2/)) 
 return(out)
end


function set_dust(X)
local out, erode_value
begin
 erode_value=0.19
 Pfil = addfile("<path/domain.lnd_file>", "r")
 PAND = Pfil->frac
 my_mask = (/ where(PAND.gt.0, 1, 0) /)
  out=my_mask*erode_value
 return(out)
end

function set_mask(X)
local out
begin
 top_in = addfile("<path/topo-bath_file>", "r")
 ptopo = top_in->topo
 my_mask = (/ where(ptopo.gt.0, 1, 0) /)
  out=my_mask
 return(out)
end

function set_trop(X)
local out
begin
 trop_in = addfile("/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/trop_mozart/ub/clim_p_trop.nc", "r")
 trop = trop_in->trop_p 
 trop_zon = dim_avg_n_Wrap( trop(time|:, lat|:, lon|:), (/2/) )
;trop_zon = dim_avg_n_wrap(trop,2)
    do i = 0,11
      do j = 0,72
        trop(i,j,:) = trop_zon(i,j)
      end do
    end do

out = (/trop/)

;do i = 0,5
;do j = 0,35
; out(time|i,lat|j) = (/ 0.5*( trop(time|i,lat|j)+trop(time|i+6,lat|72-j) ) /)
; out(time|i,lat|j+37) = (/ 0.5*( trop(time|i,lat|37+j)+trop(time|i+6,lat|37-j) ) /)
;end do
;end do
;do i = 6,11
;do j = 0,35
; out(time|i,lat|j) = (/ 0.5*( trop(time|i,lat|j)+trop(time|i-6,lat|72-j) ) /)
; out(time|i,lat|j+37) = (/ 0.5*( trop(time|i,lat|37+j)+trop(time|i-6,lat|37-j) ) /)
;end do
;end do
;
; out = (/trop/)
 return(out)
end

begin

;;; INPUT 
istem=new(20,"string") 
istem(0:16) = "/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/trop_mozart_aero/"
istem(17) = "/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/trop_mozart/"
istem(18) = "/glade/p/cesmdata/cseg/inputdata/atm/cam/"
istem(19) = "/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/trop_mozart/"

 ifil_list = (/ \
 "oxid_1.9x2.5_L26_1850clim_c091123.nc",   \
 "ar5_mam3_so2_elev_1850_c090726.nc",      \  
 "ar5_mam3_bc_elev_1850_c090726.nc",       \
 "ar5_mam3_num_a1_elev_1850_c090726.nc" ,  \
 "ar5_mam3_num_a2_elev_1850_c090726.nc" ,  \
 "ar5_mam3_oc_elev_1850_c090726.nc",       \
 "ar5_mam3_so4_a1_elev_1850_c090726.nc" ,  \
 "ar5_mam3_so4_a2_elev_1850_c090726.nc",   \
 "aerocom_mam3_dms_surf_2000_c090129.nc",  \
 "ar5_mam3_so2_surf_1850_c090726.nc",      \
 "ar5_mam3_soag_1.5_surf_1850_c100217.nc", \
 "ar5_mam3_bc_surf_1850_c090726.nc",       \
 "ar5_mam3_num_a1_surf_1850_c090726.nc" ,  \
 "ar5_mam3_num_a2_surf_1850_c090726.nc" ,  \
 "ar5_mam3_oc_surf_1850_c090726.nc",       \
 "ar5_mam3_so4_a1_surf_1850_c090726.nc" ,  \
 "ar5_mam3_so4_a2_surf_1850_c090726.nc",   \
 "clim_p_trop.nc",   \
 "dst_1.9x2.5_c090203.nc", \ 
 "regrid_vegetation.nc" \
  /)

ifil_type = (/0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4/)
ifil_name = (/"oxid/","emis/","ub/","dst/","dvel/"/)

;; This is a list of ALL the variables we want to change in all the files above.
 varNames = (/ "emiss_ene" , "emiss_ind" , "forestfire" , "grassfire" , "contvolc" , \ ;; so2_elev
             "BC_forestfire" , "BC_grassfire" , "OC_forestfire" , "OC_grassfire" , "SO4_emiss_ene" , "SO4_emiss_ind" , "SO4_forestfire" , "SO4_grassfire" , "SO4_contvolc" , \ ;; num_a1_elev ;;
             "emiss_awb" , "emiss_dom" , "emiss_tra" , "emiss_wst" , "emiss_shp" , \  ;; so2_surf
             "SOAG_BIGALK", "SOAG_BIGENE", "SOAG_ISOPRENE", "SOAG_TERPENE", "SOAG_TOLUENE" , \ soag_1.5_surf
             "BC_emiss_awb" , "BC_emiss_dom" , "BC_emiss_ene" , "BC_emiss_ind" , "BC_emiss_tra" , "BC_emiss_wst" , "BC_emiss_shp" , \ ;; num_a1_surf
             "OC_emiss_awb" , "OC_emiss_dom" , "OC_emiss_ene" , "OC_emiss_ind" , "OC_emiss_tra" , "OC_emiss_wst" , "OC_emiss_shp" , \ ;; num_a1_surf
             "SO4_emiss_awb" , "SO4_emiss_wst" , "SO4_emiss_shp" , \ ;; num_a1_surf
             "SO4_emiss_dom" , "SO4_emiss_tra" , \ ;; num_a1_surf
             "OH" , "HO2" , "O3" , "NO3" , "H2O2" , "PS" ,\ ;; oxid
             "DMS" ,\ ;; dms_surf           
             "trop_p", \ ;; clim_p_trop
             "mbl_bsn_fct_geo", \ ;; dust
             "LANDMASK" , "PCT_LAKE" , "PCT_WETLAND" , "PCT_URBAN" , "PCT_PFT" , "MONTHLY_LAI" \ ;; regrid_vegetation
             /)

 anthro = (/ 1,1,0,0,0, \ ;; so2_elev
             0,0,0,0,1,1,0,0,0, \ ;; num_a1_elev ;;
             1,1,1,1,1, \  ;; so2_surf
             1,1,1,1,1, \ soag_1.5_surf
             1,1,1,1,1,1,1, \ ;; num_a1_surf
             1,1,1,1,1,1,1, \ ;; num_a1_surf
             1,1,1, \ ;; num_a1_surf
             1,1, \ ;; num_a1_surf
             0,0,0,0,0,0, \ ;; oxid
             0, \ ;; dms_surf           
             0, \ ;; clim_p_trop
             0, \ ;; dust
             0,1,1,1,1,1 \ ;; regrid_vegetation
             /)


;; OUTPUT

 output_stem = "./"
 system( "mkdir " + output_stem )

 do i = 0,dimsizes(ifil_list)-1

   ifil_list_new=str_insert(ifil_list(i),"_for_<casename>",-4) ; Name of case/run?
   outfile=output_stem+ifil_list_new

   my_filename = istem(i) + ifil_name(ifil_type(i)) + ifil_list(i)
   system( "cp " + my_filename + " " + outfile ) ;; move original files into place
   system( "chmod a+w " + outfile )
 ;; identify the variables in the file

   print("MODIFYING: " + outfile )
   F = addfile( outfile ,"w")
   fil_vars = getfilevarnames(F)
   vars = get_vars_from_list( fil_vars, varNames) ;; return the variable names that are in the file
   if .not.all(ismissing(vars)) then
     do j = 0,dimsizes(vars)-1 ;; Loop over all then named variables in this file.
       xtmp = F->$vars(j)$
       rank = dimsizes(dimsizes(xtmp))

       xtmp1 = (/0.*xtmp/)

       if (anthro(ind(varNames.eq.vars(j))) .eq. 0) then
       if rank.eq.4 then
         if ifil_type(i).eq.0 then
           xtmp1 = paleotize(xtmp,"lev",4,"hyb")
         end if
         if ifil_type(i).eq.1.or.ifil_type(i).eq.2.or.ifil_type(i).eq.3 then
           xtmp1 = paleotize(xtmp,"altitude",4,"hyb")
         end if
         if ifil_type(i).eq.4 then
           xtmp1 = paleotize(xtmp,"pft",4,"hyb")
         end if
         if (vars(j).eq."SO4_contvolc") then
           xtmp1 = xtmp1*0.
         end if
         if (vars(j).eq."contvolc") then
           xtmp1 = xtmp1*0.
         end if
         print("paleotized: "+vars(j))
       end if

       if rank.eq.3 then 
         if (vars(j).ne."PS") then
           if ifil_type(i).eq.2 then
             xtmp1 = set_trop(xtmp)
           else
             xtmp1 = paleotize(xtmp,"dummy",3,"hyb")
           end if
         else
           xtmp1 = set_PS(xtmp)
           print("PS set to initial condition")
         end if
       end if

       if rank.eq.2 then 
         if (vars(j).eq."mbl_bsn_fct_geo") then
           xtmp1 = set_dust(xtmp)
           print("dust set")
         end if
         if (vars(j).eq."LANDMASK") then
           xtmp1 = set_mask(xtmp)
           print("mask set")
         end if
       end if

       else
         print("set to zero: "+vars(j))
       end if

       xtmp = (/xtmp1/)

if ((vars(j).ne."mbl_bsn_fct_geo").and.(vars(j).ne."LANDMASK")) then
  p = .5
  q = .0    ;  -.25 = light local smoothing... +.25 =heavy
  cyclic = True
  xtmp = smth9_Wrap(xtmp,p,q,cyclic)
end if

       F->$vars(j)$ = xtmp
       delete(xtmp)
       delete(xtmp1)
     end do
     F@brianpm_comments = (/ "Modified values with "+get_script_name()+" on " + systemfunc("date") /)
   else
     print("DID NOT MODIFY FILE: " +  outfile )
   end if    
   delete( vars)
   delete( fil_vars )
   delete( F )
     print("WRITTEN: " + outfile )
 end do

end
