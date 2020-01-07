load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Take a BIOME4 output file and map each biome to it's equivelant LSM biome 
; for input into the paleo setup of CLM. Regrids the data to 2x2 degrees.
;
; BIOME4: http://pmip2.lsce.ipsl.fr/synth/biome4.shtml
; CLM: http://www.cesm.ucar.edu/models/cesm1.0/clm/
; Good reference for biome real world descriptions: Walter, H., Vegetation of the Earth, and ecological systems of the geo-biosphere.
;
; nherold, 2014

begin
; ***********************************************
; BEGIN USER INPUT
; ***********************************************

        infile = "herold_etal_eocene_biome_1x1.nc"
        outfile = "LSM_types_1x1_from_biome4.nc"
        desired_lat = fspan(-89.5,89.5,180)        ; array of latitudes for output grid. Usually this is 2 degrees, but 1deg here.
        desired_lon = fspan(0.5,359.5,360)         ; array of longitudes for output grid. Usually this is 2 degrees, but 1deg here.

; ***********************************************
; END USER INPUT
; ***********************************************

; ***********************************************
; Read in files and variables
        add0 = addfile(infile,"r")

        biome = add0->eocene_biome

        biomedim = dimsizes(biome)
        lat = add0->lat
        lon = add0->lon
        lon = where(lon.lt.0,lon+360,lon)
        lon = lonFlip(lon)

        pft = biome
        pft = lonFlip(pft)
        pfttmp = pft

; ***********************************************
; Convert BIOME4 biomes to LSM biomes.
;
; See biome4.f for a list of biome names used in BIOME4.
; See paleo_mkraw_ccsm4_public_sed.F90_eo in the paleo_mkraw folder for a list of LSM biome names used in the CLM paleo setup process. 
        pft = where(pfttmp.eq.1,10,pft)
        pft = where(pfttmp.eq.2,11,pft)
        pft = where(pfttmp.eq.3,12,pft)
        pft = where(pfttmp.eq.4,8,pft)
        pft = where(pfttmp.eq.5,7,pft)
        pft = where(pfttmp.eq.6,9,pft)
        pft = where(pfttmp.eq.7,6,pft)
        pft = where(pfttmp.eq.8,3,pft)  ; ? 3/4
        pft = where(pfttmp.eq.9,6,pft)  ; ?
        pft = where(pfttmp.eq.10,3,pft) ; ? 3/6
        pft = where(pfttmp.eq.11,4,pft) ; ? 4/6
        pft = where(pfttmp.eq.12,12,pft)
        pft = where(pfttmp.eq.13,20,pft); ?
        pft = where(pfttmp.eq.14,21,pft); ?
        pft = where(pfttmp.eq.15,9,pft) ; ?
        pft = where(pfttmp.eq.16,12,pft)
        pft = where(pfttmp.eq.17,13,pft); ? 13/14
        pft = where(pfttmp.eq.18,17,pft); ?
        pft = where(pfttmp.eq.19,18,pft)
        pft = where(pfttmp.eq.20,17,pft)
        pft = where(pfttmp.eq.21,2,pft)
        pft = where(pfttmp.eq.22,19,pft)
        pft = where(pfttmp.eq.23,19,pft)
        pft = where(pfttmp.eq.24,19,pft)
        pft = where(pfttmp.eq.25,20,pft); ? 20/21
        pft = where(pfttmp.eq.26,17,pft); ?
        pft = where(pfttmp.eq.27,2,pft) ; ? 2/22
        pft = where(pfttmp.eq.28,1,pft)

; ***********************************************
; Create output grid
        lsmbiome = new((/dimsizes(desired_lat),dimsizes(desired_lon)/),"float")
        lsmbiome!0 = "lat"
        lsmbiome!1 = "lon"
        lsmbiome&lat = desired_lat
        lsmbiome&lon = desired_lon
        lsmbiome&lat@units = "degrees_north"
        lsmbiome&lon@units = "degrees_east"
        lsmbiome&lat@long_name = "latitude"
        lsmbiome&lon@long_name = "longitude"
        lsmbiome@long_name = "LSM biomes for CLM paleo install. Converted from BIOME4 output."

        do a = 0,dimsizes(desired_lat)-1
                do b = 0,dimsizes(desired_lon)-1
                        latind = ind_nearest_coord(desired_lat(a),lat,0)
                        lonind = ind_nearest_coord(desired_lon(b),lon,0)
                        lsmbiome(a,b) = (/pft(latind,lonind)/)
                end do
        end do

; make ocean = 0
        lsmbiome = where(ismissing(lsmbiome),0.,lsmbiome)

; make a land mask
        landmask = lsmbiome
        landmask = where(lsmbiome.gt.0,1.,0.)

; ***********************************************
; Write out data
        system("rm -f "+outfile)
        addout = addfile(outfile,"c")
        addout->pft = lsmbiome
        addout->landmask = landmask

        print("LSM BIOME FILE CREATED: "+outfile)
end
