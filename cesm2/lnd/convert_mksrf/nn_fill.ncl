;*************************************************
; 
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
;************************************************
; read in netCDF file and access data
;************************************************

  ; open output file.
  dd = systemfunc("date -u +%y%m%d")

  sourceCode = "nn_fill.ncl"

  case   = "21ka"
  path   = "./"
  input  = "$DIN_LOC_ROOT"

  finglc = "mksrf_glacier_lgm21ka.example"
  finlan = "mksrf_lanwatnDepth_lgm21ka.example.nc"
  finpft = "mksrf_pft_lgm21ka.example"

  ifile1 = path+finglc+".nc"
  ifile2 = path+finlan+".nc"
  ifile3 = path+finpft+".nc"
  ifile4 = input+"/lnd/clm2/rawdata/mksrf_soitex.10level.c010119.nc"
  ifile5 = input+"/lnd/clm2/rawdata/mksrf_organic.10level.0.5deg.081112.nc"   ;;; Use old dataset so it's at half degree
  ;ifile5 = input+"/lnd/clm2/rawdata/mksrf_organic_10level_5x5min_ISRIC-WISE-NCSCD_nlev7_c120830.nc";
  ifile6 = input+"/lnd/clm2/rawdata/mksrf_fmax.070406.nc"   ;;; Use old dataset so it's at half degree
  ;ifile6 = input+"/lnd/clm2/rawdata/mksrf_fmax_3x3min_USGS_c120911.nc"
  ;ifile7 = input+"/lnd/clm2/rawdata/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/mksrf_soilcol_global_c090324.nc"
  ifile7 = input+"/lnd/clm2/rawdata/pftcftlandusedynharv.0.25x0.25.MODIS.simyr1850-2015.c170412/mksrf_soilcolor_CMIP6_simyr2005.c170623.nc"
  ;ifile8 = input+"/lnd/clm2/rawdata/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/mksrf_lai_global_c090506.nc"
  ifile8 = input+"/lnd/clm2/rawdata/pftcftlandusedynharv.0.25x0.25.MODIS.simyr1850-2015.c170412/mksrf_lai_78pfts_simyr2005.c170413.nc"
  ifile9 = input+"/lnd/clm2/rawdata/mksrf_vocef_0.5x0.5_simyr2000.c110531.nc"

  ofile1 = path+finglc+".nn.nc"
  ofile2 = path+finlan+".nn.nc"
  ofile3 = path+finpft+".nn.nc"
  ofile4 = path+"mksrf_soitex.10level_"+case+"."+dd+".nn.nc"
  ofile5 = path+"mksrf_organic.10level.0.5deg_"+case+"."+dd+".nn.nc"
  ofile6 = path+"mksrf_fmax_"+case+"."+dd+".nn.nc"
  ofile7 = path+"mksrf_soilcol_"+case+"."+dd+".nn.nc"
  ofile8 = path+"mksrf_lai_global_"+case+"."+dd+".nn.nc"
  ofile9 = path+"mksrf_vocef_0.5x0.5_simyr2000_"+case+"."+dd+".nn.nc"


  ; glaciers
  print( "Open: "+ifile1 );
  in1 = addfile(ifile1,"r")
  pctgla  = in1->PCT_GLACIER
  newmask = in1->LANDMASK
  lgmlat  = in1->LAT
  lgmlon  = in1->LON
  in1_global = getvaratts(in1)
  in1_vars   = getfilevarnames(in1)

  ; lakes and wetlands
  print( "Open: "+ifile2 );
  in2 = addfile(ifile2,"r")
  pctlak = in2->PCT_LAKE
  pctwet = in2->PCT_WETLAND
  in2_global = getvaratts(in2)
  in2_vars   = getfilevarnames(in2)
  if ( any(newmask .ne. in2->LANDMASK) )then
     print( "Lake/wetlands landmask does NOT match the glacier mask" )
     exit;
  end if

  ; pft
  print( "Open: "+ifile3 );
  in3 = addfile(ifile3,"r")
  pctnatvg = in3->PCT_NATVEG
  pctcrop  = in3->PCT_CROP
  pctpft   = in3->PCT_NAT_PFT
  pctcft   = in3->PCT_CFT
  in3_global = getvaratts(in3)
  in3_vars   = getfilevarnames(in3)
  if ( any(newmask .ne. in3->LANDMASK) )then
     print( "PFT landmask does NOT match the glacier mask" )
     exit;
  end if

  npctpft = pctpft
  opctpft = pctpft



  ; soitex
  domapunits = 0
  if (domapunits .eq. 1) then
    print( "Open: "+ifile4 );
    in4 = addfile(ifile4,"r")
    mapunits   = in4->MAPUNITS
    maskmap    = in4->LANDMASK
    latmap     = in4->LAT
    lonmap     = in4->LON
    in4_global = getvaratts(in4)
    in4_vars   = getfilevarnames(in4)
  end if

  ; organic
  print( "Open: "+ifile5 );
  in5 = addfile(ifile5,"r")
  orig_organic3d = in5->ORGANIC
  organic3d = in5->ORGANIC
  ZSOI      = in5->ZSOI
  in5_global = getvaratts(in5)
  in5_vars   = getfilevarnames(in5)
  nlayers = dimsizes(ZSOI)


  ; fmax
  print( "Open: "+ifile6 );
  in6 = addfile(ifile6,"r")
  fmax   = in6->FMAX
  pdmask = in6->LANDMASK
  in6_global = getvaratts(in6)
  in6_vars   = getfilevarnames(in6)

  ; soilcolor
  print( "Open: "+ifile7 );
  in7 = addfile(ifile7,"r")
  soic2d     = area_conserve_remap ( in7->LON, in7->LAT, in7->SOIL_COLOR, lgmlon, lgmlat, False );
  soicmask   = area_conserve_remap ( in7->LON, in7->LAT, in7->LANDMASK,   lgmlon, lgmlat, False );
  in7_global = getvaratts(in7)
  in7_vars   = getfilevarnames(in7)
  newsoicmask = where(newmask .eq. 0, 1, soicmask)

  ; lai
  print( "Open: "+ifile8 );
  in8 = addfile(ifile8,"r")
  hbot = new( (/dimsizes(in8->time), dimsizes(in8->pft), dimsizes(lgmlat), dimsizes(lgmlon)/), "double" );
  htop = hbot;
  lai  = hbot;
  sai  = hbot;
  do itime = 0, dimsizes(in8->time)-1
     print( "time: "+in8->time(itime) );
     do ipft = 0, dimsizes(in8->pft)-1
        ; uniform for given PFT
        hbot(itime,ipft,:,:)       = area_conserve_remap ( in8->LON, in8->LAT, in8->MONTHLY_HEIGHT_BOT(itime,ipft,:,:), lgmlon, lgmlat, False );
        ; uniform for given PFT
        htop(itime,ipft,:,:)       = area_conserve_remap ( in8->LON, in8->LAT, in8->MONTHLY_HEIGHT_TOP(itime,ipft,:,:), lgmlon, lgmlat, False );
        ; nn interpolate
        lai(itime,ipft,:,:)        = area_conserve_remap ( in8->LON, in8->LAT, in8->MONTHLY_LAI(itime,ipft,:,:),        lgmlon, lgmlat, False );
        ; nn interpolate
        sai(itime,ipft,:,:)        = area_conserve_remap ( in8->LON, in8->LAT, in8->MONTHLY_SAI(itime,ipft,:,:),        lgmlon, lgmlat, False );
     end do
  end do
  in8_global = getvaratts(in8)
  in8_vars   = getfilevarnames(in8)

  ; vocef
  print( "Open: "+ifile9 );
  in9 = addfile(ifile9,"r")
  in9_vars   = getfilevarnames(in9)

  numpft = 15
  pdmaskfill = pdmask

  minlat = 0
  minlon = 0
  maxlat = 360-1
  maxlon = 720-1
  ; fill in 0 and 359

  ; New land points
doThis = 1
if (doThis .eq. 1) then
  do i=minlat,maxlat  				; lat
  do j=minlon,maxlon                            ; lon
    if(pdmask(i,j) .eq. 0 .and. newmask(i,j) .eq. 1) then

  	; fill in  fmax on longitude seams
  	fmax(i,j) = where(pdmask(i,j) .eq. 0 .and. newmask(i,j) .eq.1, 0.365782, fmax(i,j))
  	;fmax(i,j) = where(pdmask(i,j) .eq. 0 .and. newmask(i,j) .eq.1, 0.365782, fmax(:,minlon))

	assigned = 0

	; assign bareground under new glaciers
       	if (pctgla(i,j) .eq. 100) then
        	soic2d(i,j)  	 = 15
		fmax(i,j)    	 = 0.365241	; constant over glaciers
                pctwet(i,j)      = 0.
                pctlak(i,j)      = 0.
                organic3d(:,i,j) = 0.
                npctpft(0,i,j)    = 100.		; bareground under new glaciers
                npctpft(1:numpft-1,i,j) = 0.
                hbot(:,:,i,j)    = 0.
                htop(:,:,i,j)    = 0.
                lai( :,:,i,j)    = 0.
                sai( :,:,i,j)    = 0.
		assigned = 1
 	end if
    end if
  end do
  end do
end if		; doThis loop
	



  ; nearest neighbor interpolation ---------------------------------------

  do_nn = 1

  if (do_nn .eq. 1)
  print( "Do the nearest neighbor interpolation" );
  maxlatm2 = 360-2
  maxlonm2 = 720-2
  minlatp1 = 1
  minlonp1 = 1
  ; fill in  fmax on longitude seams
  fmax(:,maxlon) = where(pdmask(:,maxlon) .eq. 0 .and. newmask(:,maxlon) .eq.1, 0.365241, fmax(:,maxlon))
  fmax(:,minlon) = where(pdmask(:,minlon) .eq. 0 .and. newmask(:,minlon) .eq.1, 0.365241, fmax(:,minlon))

  maxiter = 4
  do niter = 0,maxiter
  do i=minlatp1,maxlatm2				; lat (ignore 1st and last cells)
  do j=minlonp1,maxlonm2                                ; lon (ignore 1st and last cells)

    ; create new land over land ice points (e.g., Hudson bay)
    if(pdmaskfill(i,j) .eq. 0 .and. newmask(i,j) .eq. 1) then

	fmax(i,j) = 0.365241
       	if (pctgla(i,j) .eq. 100) then
                print( "Fill: i="+i+" j="+j );
                soic2d(i,j)      = 15
                pctwet(i,j)      = 0.
                pctlak(i,j)      = 0.
                organic3d(:,i,j) = 0.
                npctpft(0,i,j)    = 100.
                npctpft(1:numpft-1,i,j) = 0.
                hbot(:,:,i,j)    = 0.
                htop(:,:,i,j)    = 0.
                lai( :,:,i,j)    = 0.
                sai( :,:,i,j)    = 0.
		found = 1
	else 
  	     noluck = 0
  	     found = 0
             ; look west

             if (pdmaskfill(i,j-1) .eq. 1 .and. npctpft(0,i,j-1) .ne. 100.) then
                usei = i
                usej = j-1
                found = 1
		print("found cell to west i,j, usei, usej = " + i + " " + j + " " + usei + " " + usej)
             else 
		; look east
		if (pdmaskfill(i,j+1) .eq. 1 .and. npctpft(0,i,j+1) .ne. 100.) then
                   usei = i
                   usej = j+1
                   found = 1
		   print("found cell to east i,j,usei,usej = " + i + " " + j + " " + usei + " " + usej)
                else 
		   ; look north
		   if (pdmaskfill(i+1,j) .eq. 1 .and. npctpft(0,i+1,j) .ne. 100.) then
                	usei = i+1
                	usej = j
                	found = 1
		   	print("found cell to North i,j,usei,usej = " + i + " " + j + " " + usei + " " + usej)
             	   else 
			; look south
			if (pdmaskfill(i-1,j) .eq. 1 .and. npctpft(0,i-1,j) .ne. 100.) then
                		usei = i-1
                		usej = j
                		found = 1
		   		print("found cell to south i,j,usei,usej = " + i + " " + j + " " + usei + " " + usej)
             		else
                		if (niter .eq. maxiter) then
					noluck = 1
				end if
			end if
		   end if
		 end if
             end if
             if (noluck .eq. 1) then
		print("Cell will remain bareground, i,j " + i + " " + j)
                soic2d(i,j)      = 15
                pctwet(i,j)      = 0.
                pctlak(i,j)      = 0.
                organic3d(:,i,j) = 0.
                hbot(:,:,i,j)    = 0.
                htop(:,:,i,j)    = 0.
                lai( :,:,i,j)    = 0.
                sai( :,:,i,j)    = 0.
		pdmaskfill(i,j)      = 1.
		soicmask(i,j)    = 1.
             else 
		if (found  .eq. 1) then
                	; assign new cell
			pdmaskfill(i,j) = 1
			if (pctgla(usei,usej) .eq. 100.) then
                		npctpft(0,i,j) = 100.
                		npctpft(1:numpft-1,i,j) = 0.
                		pctgla(i,j)   = 100.
                		hbot(:,:,i,j) = 0.
                		htop(:,:,i,j) = 0.
                		lai( :,:,i,j) = 0.
                		sai( :,:,i,j) = 0.
			end if
                	soic2d(i,j)   = soic2d(usei,usej)
			fmax(i,j)     = fmax(usei,usej)
                	organic3d(:,i,j) = organic3d(:,usei,usej)
                	hbot(:,:,i,j) = hbot(:,:,usei,usej)
                	htop(:,:,i,j) = htop(:,:,usei,usej)
                	lai( :,:,i,j) = lai( :,:,usei,usej)
                	sai( :,:,i,j) = sai( :,:,usei,usej)
			
			; fill seams at 0 and 360
			if (j .eq. minlonp1 .and. pdmaskfill(i,minlonp1) .eq. 0 .and. newmask(i,minlonp1) .eq. 1) then
                		lai( :,:,i,minlonp1) = lai( :,:,usei,usej)
                		sai( :,:,i,minlonp1) = sai( :,:,usei,usej)
                		hbot(:,:,i,minlonp1) = hbot(:,:,usei,usej)
                		htop(:,:,i,minlonp1) = htop(:,:,usei,usej)
				fmax(    i,minlonp1) = 0.365783
                		organic3d(:,i,minlonp1) = organic3d(:,usei,usej)
			end if
			if (j .eq. maxlonm2 .and. pdmaskfill(i,maxlonm2) .eq. 0 .and. newmask(i,maxlonm2) .eq. 1) then
                		lai( :,:,i,maxlonm2) = lai( :,:,usei,usej)
                		sai( :,:,i,maxlonm2) = sai( :,:,usei,usej)
                		hbot(:,:,i,maxlonm2) = hbot(:,:,usei,usej)
                		htop(:,:,i,maxlonm2) = htop(:,:,usei,usej)
				fmax(    i,maxlonm2) = 0.365783
                		organic3d(:,i,maxlonm2) = organic3d(:,usei,usej)
			end if
                end if
              end if

	      delete(sumpft)

	      sumpft = 0.
              do m = 0, numpft-1
                          sumpft = sumpft + npctpft(m,i,j)
              end do
	      
              ddiff = 100 - sumpft
              if (sumpft .ne. 100) then
			
			print("sumpft = "+sumpft + " for cell " + i + " " +j)
			; zero cell is set to bareground
			if (sumpft .eq. 0) then
				pctlak(i,j) = 0.
				pctwet(i,j) = 0.
				pctgla(i,j) = 0.
				npctpft(0,i,j) = 100.
				npctpft(1:numpft-1,i,j) = 0.
			else 
			     ; partially glaciated cell is set to part bareground, part glacier
			     ; add any missing ground to bareground
			     if (pctgla(i,j) .eq. 0. .and. sumpft .lt. 100.) then
				oldbare = npctpft(0,i,j)
			   	npctpft(0,i,j)  = oldbare + 100.-sumpft
			     else
			        ; bareground is adjusted
			     	if (ddiff .gt. 0.) then
			   		npctpft(0,i,j)  = npctpft(0,i,j) + 100.-sumpft
				else
			   		npctpft(0,i,j)  = npctpft(0,i,j) - 100.-sumpft
	      		     	end if
	      		     end if
                        end if
			; check again
              		sumpft = 0.
              		do m = 0, numpft-1
                             sumpft = sumpft + npctpft(m,i,j)
              		end do
			print("new sumpft = " + sumpft + " for cell " + i + " " +j)
               end if
		delete(sumpft)
        end if 		; 100% glacier if loop
    end if 		; matching landmask / nn loop
        ; check all pfts
        sumpft = sum(npctpft(:,i,j))
        if(newmask(i,j) .eq. 1. .and. pctgla(i,j) .gt. 0. .and.  sumpft .ne. 100. ) then    ;this solves problem for land points with 0<pctgla<100
              npctpft(0,i,j) = 100.
              npctpft(1:numpft-1,i,j) = 0.
              sumpft = sum(npctpft(:,i,j))
        else if(newmask(i,j) .eq. 1. .and. pctgla(i,j) .eq. 0. .and.  sumpft .ne. 100. ) then    ;this solves problem for land points with pctgla=0
              npctpft(0,i,j) = 100.
              npctpft(1:numpft-1,i,j) = 0.
              sumpft = sum(npctpft(:,i,j))
        end if
        end if
        if (sumpft .ne. 100.0 .and. newmask(i,j) .eq. 1) then
		sumpft = 0.
		osumpft = 0.
                m = numpft-1
                sumpft  = dim_sum_n(npctpft(m,i,j),0)
                osumpft = dim_sum_n(opctpft(m,i,j),0)
		if (sumpft .ne. 100.) then
                	print("Exit WARNING sumpft was ne 100 : " + i + " " + j + " corrected sum=  " + sumpft )
                	print("Exit WARNING sumpft was ne 100 : " + i + " " + j + " old sum=  " + osumpft )
			print(" ------------ ")
			exit
		end if
        end if
    end do			; lon loop
    end do			; lat loop
  end do      			; iteration loop

  end if			; do_nn loop

; -----------------------------------------------
; Write output files
; -----------------------------------------------
  system("/bin/rm -f " + ofile1)
  system("/bin/rm -f " + ofile2)
  system("/bin/rm -f " + ofile3)
  system("/bin/rm -f " + ofile4)
  system("/bin/rm -f " + ofile5)
  system("/bin/rm -f " + ofile6)
  system("/bin/rm -f " + ofile7)
  system("/bin/rm -f " + ofile8)
  system("/bin/rm -f " + ofile9)
  o1 = addfile(ofile1, "c")
  o2 = addfile(ofile2, "c")
  o3 = addfile(ofile3, "c")
  o4 = addfile(ofile4, "c")
  o5 = addfile(ofile5, "c")
  o6 = addfile(ofile6, "c")
  o7 = addfile(ofile7, "c")
  o8 = addfile(ofile8, "c")
  o9 = addfile(ofile9, "c")

; ========================================
  print(" glacier")
; ========================================
 ; copy variables from standard file
  ; copy global attributes
  if (.not.all(ismissing(in1_global))) then
        do i = 0, dimsizes(in1_global) - 1
                print("copy_fileatts: global attributes->" + in1_global )
                o1@$in1_global$ = in1@$in1_global$
        end do
        o1@revisionCode = sourceCode
        o1@revisionDate = dd
    	o1@input  = ifile1
  end if

  ; copy variables
  if (.not.all(ismissing(in1_vars))) then
        do i = 0, dimsizes(in1_vars) - 1
                usevar = in1_vars(i)
                print("i = " + i + " in1vars = " + in1_vars(i))
                if (in1_vars(i) .ne."PCT_GLACIER") then
                        o1->$in1_vars(i)$ = in1->$in1_vars(i)$
                        print("writing var ->" + in1_vars(i) )
                else
                        o1->$in1_vars(i)$ = pctgla
                        print("writing var ->" + in1_vars(i) )
                end if
        end do
  end if

; ========================================
  print(" LANWAT")
; ========================================
 ; copy variables from standard file
  ; copy global attributes
  if (.not.all(ismissing(in2_global))) then
        do i = 0, dimsizes(in2_global) - 1
                print("copy_fileatts: global attributes->" + in2_global )
                o2@$in2_global$ = in2@$in2_global$
        end do
        o2@revisionCode = sourceCode
        o2@revisionDate = dd
    	o2@input  = ifile2
  end if

  ; copy variables
  if (.not.all(ismissing(in2_vars))) then
        do i = 0, dimsizes(in2_vars) - 1
                usevar = in2_vars(i)
                print("i = " + i + " in2vars = " + in2_vars(i))
                if (in2_vars(i) .ne."PCT_LAKE".and.in2_vars(i) .ne."PCT_WETLAND" .and. in2_vars(i) .ne. "LANDMASK") then
                        o2->$in2_vars(i)$ = in2->$in2_vars(i)$
                        print("writing var ->" + in2_vars(i) )
                else 
                   if (in2_vars(i) .eq. "PCT_LAKE") then
                        o2->$in2_vars(i)$ = pctlak
                        print("writing var ->" + in2_vars(i) )
                   else 
                        if (in2_vars(i) .eq. "PCT_WETLAND") then
                                o2->$in2_vars(i)$ = pctwet
                                print("writing var ->" + in2_vars(i) )
                        else
                        	if (in2_vars(i) .eq. "LANDMASK") then
					tmp1 = in2->$in2_vars(i)$
					tmp2 = tmp1
					tmp2 = (/newmask(:,:)/)
                        		o2->$in2_vars(i)$ = tmp2
                        		print("writing var ->" + in2_vars(i) )
					delete([/tmp1,tmp2/])
                        	else
                                	print("Done Writing")
                        	end if
                        end if
                   end if
                end if
        end do
  end if

; ========================================
  print(" PFT")
; ========================================
 ; copy global attributes
  if (.not.all(ismissing(in3_global))) then
        do i = 0, dimsizes(in3_global) - 1
                print("copy_fileatts: global attributes->" + in3_global )
                o3@$in3_global$ = in3@$in3_global$
        end do
        o3@revisionCode = sourceCode
        o3@revisionDate = dd
    	o3@input  = ifile3
  end if

  ; copy variables
  if (.not.all(ismissing(in3_vars))) then
        do i = 0, dimsizes(in3_vars) - 1
                usevar = in3_vars(i)
                print("i = " + i + " in3vars = " + in3_vars(i))
                if (in3_vars(i) .ne."PCT_PFT" .and. in3_vars(i) .ne. "LANDMASK") then
                        o3->$in3_vars(i)$ = in3->$in3_vars(i)$
                        print("writing var ->" + in3_vars(i) )
                else 
                	if (in3_vars(i) .eq."PCT_PFT") then
                        	o3->$in3_vars(i)$ = npctpft
			else
                		if (in3_vars(i) .eq."LANDMASK") then
					tmp1 = in3->$in3_vars(i)$
					tmp2 = tmp1
					tmp2 = (/newmask(:,:)/)
                        		o3->$in3_vars(i)$ = tmp2
                        		print("writing var ->" + in3_vars(i) )
					delete([/tmp1,tmp2/])
				else
                        		print("writing var ->" + in3_vars(i) )
				end if
                        end if
                end if
        end do
  end if

; ========================================
  print(" soitex ")
; ========================================

if (domapunits .eq. 1) then

 ; copy global attributes
  if (.not.all(ismissing(in4_global))) then
        do i = 0, dimsizes(in4_global) - 1
                print("copy_fileatts: global attributes->" + in4_global )
                o4@$in4_global$ = in4@$in4_global$
        end do
        o4@revisionCode = sourceCode
        o4@revisionDate = dd
    	o4@input  = ifile4
  end if


  ; copy variables
  if (.not.all(ismissing(in4_vars))) then
        do i = 0, dimsizes(in4_vars) - 1
                usevar = in4_vars(i)
                print("i = " + i + " in4_vars = " + in4_vars(i))
                if (in4_vars(i) .ne."MAPUNITS" .and.in4_vars(i) .ne."LANDMASK") then
                        o4->$in4_vars(i)$ = in4->$in4_vars(i)$
                        print("writing var ->" + in4_vars(i) )
                else 
                   if (in4_vars(i) .eq. "MAPUNITS") then
                        o4->$in4_vars(i)$ = mapunits
                        print("writing var ->" + in4_vars(i) )
	 	   else
                        if (in4_vars(i) .eq. "LANDMASK") then
                        	o4->$in4_vars(i)$ = maskmap
                        	print("writing var ->" + in4_vars(i) )
                   	end if
                   end if
                end if
        end do
  end if

end if		; do map units look


; ========================================
  print(" organic ")
; ========================================


 ; copy global attributes
  if (.not.all(ismissing(in5_global))) then
        do i = 0, dimsizes(in5_global) - 1
                print("copy_fileatts: global attributes->" + in5_global )
                o5@$in5_global$ = in5@$in5_global$
        end do
        o5@revisionCode = sourceCode
        o5@revisionDate = dd
    	o5@input  = ifile5
  end if


  ; copy variables
  if (.not.all(ismissing(in5_vars))) then
        do i = 0, dimsizes(in5_vars) - 1
                usevar = in5_vars(i)
                print("i = " + i + " in5vars = " + in5_vars(i))
                if (in5_vars(i) .ne."ORGANIC" .and. in5_vars(i) .ne. "LANDMASK") then
                        o5->$in5_vars(i)$ = in5->$in5_vars(i)$
                        print("writing var ->" + in5_vars(i) )
                else
                	if (in5_vars(i) .eq."ORGANIC") then
                        	o5->$in5_vars(i)$ = organic3d
                        	print("writing var ->" + in5_vars(i) )
			else
                		if (in5_vars(i) .eq."LANDMASK") then
					tmp1 = in5->$in5_vars(i)$
					tmp2 = tmp1
					tmp2 = (/newmask(:,:)/)
                        		o5->$in5_vars(i)$ = tmp2
                        		print("writing var ->" + in5_vars(i) )
					delete([/tmp1,tmp2/])
				end if
                	end if
                end if
        end do
  end if

; ========================================
  print(" fmax ")
; ========================================


 ; copy global attributes
  if (.not.all(ismissing(in6_global))) then
        do i = 0, dimsizes(in6_global) - 1
                print("copy_fileatts: global attributes->" + in6_global )
                o6@$in6_global$ = in6@$in6_global$
        end do
        o6@revisionCode = sourceCode
        o6@revisionDate = dd
    	o6@input  = ifile6
  end if


  ; copy variables
  if (.not.all(ismissing(in6_vars))) then
        do i = 0, dimsizes(in6_vars) - 1
                usevar = in6_vars(i)
                print("i = " + i + " in6vars = " + in6_vars(i))
                if (in6_vars(i) .ne."FMAX" .and. in6_vars(i) .ne."LANDMASK") then
                        o6->$in6_vars(i)$ = in6->$in6_vars(i)$
                        ; o6->$in6_vars(i)$ = (/in6->$in6_vars(i)$/)
                        print("writing var ->" + in6_vars(i) )
                else
                	if (in6_vars(i) .eq."FMAX") then
                        	o6->$in6_vars(i)$ = fmax
                        	print("writing var ->" + in6_vars(i) )
			else
                		if (in6_vars(i) .eq."LANDMASK") then
					tmp1 = in6->$in6_vars(i)$
					tmp2 = tmp1
					tmp2 = (/newmask(:,:)/)
                        		o6->$in6_vars(i)$ = tmp2
                        		print("writing var ->" + in6_vars(i) )
					delete([/tmp1,tmp2/])
				end if
                	end if
                end if
        end do
  end if


; ========================================
  print(" Soilcolor ")
; ========================================


 ; copy global attributes
  if (.not.all(ismissing(in7_global))) then
        do i = 0, dimsizes(in7_global) - 1
                print("copy_fileatts: global attributes->" + in7_global )
                o7@$in7_global$ = in7@$in7_global$
        end do
        o7@revisionCode = sourceCode
        o7@revisionDate = dd
    	o7@input  = ifile7
  end if


  ; copy variables
  if (.not.all(ismissing(in7_vars))) then
        o7->LATIXY = in3->LATIXY;
        o7->LONGXY = in3->LONGXY;
        do i = 0, dimsizes(in7_vars) - 1
                usevar = in7_vars(i)
                print("i = " + i + " in7vars = " + in7_vars(i))
               	if (in7_vars(i) .eq."SOIL_COLOR") then
                       	o7->$in7_vars(i)$ = soic2d
                       	print("writing var ->" + in7_vars(i) )
		else
               		if (in7_vars(i) .eq."LANDMASK") then
                       		o7->$in7_vars(i)$ = newmask
                       		print("writing var ->" + in7_vars(i) )
                        else
               		    if (in7_vars(i) .eq."LAT") then
                       	    	o7->$in7_vars(i)$ = lgmlat
                       		print("writing var ->" + in7_vars(i) )
                            else
               		        if (in7_vars(i) .eq."LON") then
                       	            o7->$in7_vars(i)$ = lgmlon
                       	            print("writing var ->" + in7_vars(i) )
                                 else
                                           print( "Don't know how to handle this variable:"+in7_vars(i) );
                                 end if
                            end if
               		end if
                end if
        end do
  end if

; ========================================
  print(" Lai ")
; ========================================


 ; copy global attributes
  if (.not.all(ismissing(in8_global))) then
        do i = 0, dimsizes(in8_global) - 1
                print("copy_fileatts: global attributes->" + in8_global )
                o8@$in8_global$ = in8@$in8_global$
        end do
        o8@revisionCode = sourceCode
        o8@revisionDate = dd
    	o8@input  = ifile8
  end if


  ; copy variables
  if (.not.all(ismissing(in8_vars))) then
        o8->LATIXY = in3->LATIXY;
        o8->LONGXY = in3->LONGXY;
        o8->pft    = in8->pft;
        o8->time   = in8->time;
        do i = 0, dimsizes(in8_vars) - 1
                usevar = in8_vars(i)
                print("i = " + i + " in8vars = " + in8_vars(i))
                if (in8_vars(i) .ne."MONTHLY_HEIGHT_BOT" .and. in8_vars(i) .ne. "LANDMASK" .and. \
		    in8_vars(i) .ne."MONTHLY_HEIGHT_TOP" .and. in8_vars(i) .ne. "MONTHLY_LAI" .and. \
		    in8_vars(i) .ne."MONTHLY_SAI") then
                        if ( in8_vars(i) .eq. "LAT" )then
                           o8->$in8_vars(i)$ = lgmlat
                        else
                           if ( in8_vars(i) .eq. "LON" )then
                              o8->$in8_vars(i)$ = lgmlon
                           else
                              print( "Can't handle this variable" );
                           end if
                        end if
                        print("writing var ->" + in8_vars(i) )
                else
                  if (in8_vars(i) .eq."MONTHLY_HEIGHT_BOT") then
                        o8->$in8_vars(i)$ = hbot
                        print("writing var ->" + in8_vars(i) )
                  else
                     if (in8_vars(i) .eq."LANDMASK") then
                        o8->$in8_vars(i)$ = newmask
                        print("writing var ->" + in8_vars(i) )
                     else
                        if (in8_vars(i) .eq."MONTHLY_HEIGHT_TOP") then
                           o8->$in8_vars(i)$ = htop
                           print("writing var ->" + in8_vars(i) )
                        else
                           if (in8_vars(i) .eq."MONTHLY_LAI") then
                              o8->$in8_vars(i)$ = lai
                              print("writing var ->" + in8_vars(i) )
                           else
                              if (in8_vars(i) .eq."MONTHLY_SAI") then
                                 o8->$in8_vars(i)$ = sai
                                 print("writing var ->" + in8_vars(i) )
                              end if
                           end if
                        end if
                     end if
                  end if
                end if
        end do
  end if

; ========================================
  print(" vocef ")
; ========================================


 ; global attributes
    o9@revisionCode = sourceCode
    o9@revisionDate = dd
    o9@input  = ifile9


  ; copy variables
  if (.not.all(ismissing(in9_vars))) then
        do i = 0, dimsizes(in9_vars) - 1
                usevar = in9_vars(i)
                print("i = " + i + " in9vars = " + in9_vars(i))
                if (in9_vars(i) .ne."LANDMASK") then
                        o9->$in9_vars(i)$ = in9->$in9_vars(i)$
                        print("writing var ->" + in9_vars(i) )
                else
			tmp1 = in9->$in9_vars(i)$
			tmp2 = tmp1
			tmp2 = (/newmask(:,:)/)
                        o9->$in9_vars(i)$ = tmp2
                        print("writing var ->" + in9_vars(i) )
			delete([/tmp1,tmp2/])
                end if
        end do
  end if


  print( "Successfully ran nn_fill" );  

end

