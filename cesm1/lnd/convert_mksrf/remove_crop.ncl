begin
; author:  Nan Rosenbloom
; may2014
; Description:
; 1)  remove crop (pft=15).  redistribute crop to other pfts w/in cell.  If
; the cell is 100% crop, use nearest neighbors to fill in pft distribution.
; 2)  add harvest variables and set them to 0.0

 dd = systemfunc("date -u +%y%m%d")

ifilename = "mksrf_pft_15ka.140501.nc"
ofilename = "mksrf_pft_15ka.140501.noCrop."+dd+".nc"
ifile = addfile(ifilename , "r")
system("/bin/rm -f " + ofilename)
ofile = addfile(ofilename , "c")

ilat = ifile->LAT
ilon = ifile->LON


; copy global variables
in1_global= getvaratts(ifile)
in1_vars = getfilevarnames(ifile)

; copy global attributes
if (.not.all(ismissing(in1_global))) then
do i = 0, dimsizes(in1_global) - 1
	print("copy_fileatts: global attributes->" + in1_global(i) )
	ofile@$in1_global(i)$ = ifile@$in1_global(i)$
end do
end if

; copy variables
if (.not.all(ismissing(in1_vars))) then
do i = 0, dimsizes(in1_vars) - 1
if (in1_vars(i) .ne."PCT_PFT") then
	ofile->$in1_vars(i)$ = ifile->$in1_vars(i)$
	print("writing var ->" + in1_vars(i) )
else 
  
  crops  = (ifile->PCT_PFT(15,:,:))
  
  ; read in pctpfts and set crops to 0.
  oldpft = (ifile->PCT_PFT)
  newpft = oldpft
  newpft(15,:,:) = 0
  newpft(16,:,:) = 0


  ; sum over non-crop, vegetated area (not bareground)
  sumpft = dim_sum_n(newpft(1:14,:,:),0)
  nlat   = dimsizes(sumpft(:,0))
  nlon   = dimsizes(sumpft(0,:))
  ; first loop:  replace crops in cells where pft != 100% crop
  do lat = 0,nlat-1
  do lon = 0,nlon-1
     ; where possible, divide crops into pfts 1:14.
     if (sumpft(lat,lon) .ne. 0) then
	dlat = 77
	dlon = 215
  	do pft = 1,14
     		; newpft(pft,lat,lon) = oldpft(pft,lat,lon) + (oldpft(pft,lat,lon)/sumpft(lat,lon))*oldpft(15,lat,lon)
     		t = oldpft(pft,lat,lon) + (oldpft(pft,lat,lon)/sumpft(lat,lon))*oldpft(15,lat,lon)
		tt = (round(t,1) * 10.)/10.
		newpft(pft,lat,lon) = tt
		delete([/t,tt/])
  	end do
     ; otherwise, put the crop area into bareground
     else 
     	barecrop = oldpft(0,lat,lon) + oldpft(15,lat,lon)
     	if (barecrop .eq. 100 )
		newpft(0,lat,lon)  = 100.
     	end if	
     end if
  end do
  end do

  ; second loop:  replace crops in cells where pft == 100% crop with nearest neighbor
  do lat = 0,nlat-1
  do lon = 0,nlon-1
     if (crops(lat,lon) .eq. 100) then
        npts = 4
        do p = 1,npts
           lonp1 = lon+p
           lonm1 = lon-p
	   if(crops(lat,lonp1) .ne. 100 .and. sumpft(lat,lonp1) .ne. 0) then
	        print("looking east ---")
		; crops(lat,lon) = crops(lat,lonp1)
		newpft(:,lat,lon) = newpft(:,lat,lonp1)
	   else
		if(crops(lat,lonm1) .ne. 100 .and. sumpft(lat,lonm1) .ne. 0) then
	        	print("looking west ---")
			; crops(lat,lon) = crops(lat,lonm1)
			newpft(:,lat,lon) = newpft(:,lat,lonm1)
		else
			latm1 = lat - 1
			latp1 = lat + 1
			if(crops(latm1,lon) .ne. 100. .and. sumpft(latm1,lon) .ne. 0) then
	        		print("looking south ---")
				; crops(lat,lon) = crops(latm1,lon)
				newpft(:,lat,lon) = newpft(:,latm1,lon)
			else
				if(crops(latp1,lon) .ne. 100. .and. sumpft(latp1,lon) .ne. 0) then
	        			print("looking north ---")
					; crops(lat,lon) = crops(latp1,lon)
					newpft(:,lat,lon) = newpft(:,latp1,lon)
				else
					print("surrounded by crops E and W. ;  last resort: bareground for lat= "+lat + " lon= " + lon)
					pctpft(0,lat,lon) = crops(lat,lon)
					crops(lat,lon) = 0
					exit
                		end if	
                	end if	
                end if	
           end if	
        end do
     end if	
  end do
  end do
  do lat = 0,nlat-1
  do lon = 0,nlon-1
     newsum = 0.
     oldsum = 0.
     newsum = dim_sum_n(newpft(:,lat,lon),0)
     oldsum = dim_sum_n(oldpft(:,lat,lon),0)
     if (newsum .ne. 100  .and. newsum .ne. 0 ) then
	if (lat .eq. 77 .and. lon .eq. 215) then
		print("pfts don't sum to 100%.: " + newsum +  " ... adjusting them ..")
	end if
	t = (round(newpft(:,lat,lon),1) * 10.)/10.
  	nsumpft = dim_sum_n(t,0)
	resid = 100.-nsumpft
	if (resid .ne. 0.) then
		newmax = 0.
		oldmax = 0.
		maxp   = 0
		do p  = 0,16
			newmax = t(p)
			if (newmax .gt. oldmax) then
				oldmax = newmax
				maxp   = p	
			end if
		end do
		
		tt = t
		tt(maxp) = t(maxp) + resid
		newpft(:,lat,lon) = (/tt/)
	else
		tt = t
		newpft(:,lat,lon) = tt
	end if
  	fsumpft = dim_sum_n(tt,0)
	if (fsumpft .ne. 100.) then
		print("pfts don't sum to 100 " + fsumpft + " " + resid )
		exit
	end if
	delete([/t,nsumpft/])
     end if
  end do
  end do

end if
end do
end if			; end pft loop

  outpft = oldpft	
  outpft(:,:,:) = (/newpft(:,:,:)/)	; trick to add pft coordinate variable
  ofile->PCT_PFT = outpft
  ofile@source   = "remove_crop.ncl"
  ofile@in1   = ifilename
  ofile@date  = dd
  

; add harvest variables from standard CESM input file, and zero them out.
ifile1 = "/glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/mksrf_landuse_rc2000_c090630.nc"
ivars = (/"HARVEST_VH1","HARVEST_VH2", "HARVEST_SH1","HARVEST_SH2","HARVEST_SH3","GRAZING"/)
in1 = addfile(ifile1,"r")
do i = 0, dimsizes(ivars) - 1
                print("writing var ->" + ivars(i) )
                usevar = in1->$ivars(i)$
                usevar = usevar * 0.0
                usevar!0 = "lat"
                usevar!1 = "lon"
                usevar&lat = ilat
                usevar&lon = ilon
                ofile->$ivars(i)$ = usevar
                delete(usevar)
end do


end
