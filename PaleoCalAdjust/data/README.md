## Data folders and files, and program input and output ##

Data folders:

	/figure_data	! data files (or links) to data used to make the figures
	/GISS_orbital	! orbital-parameter output
	/info_files	! info files
	/month_lengths	! month-length tables 
	/nc_files	! download instructions for netCDF files
	/TraCE_example	! time series of TraCE-21k data for two example grid cells
	/debug_files	! used for debugging -- not included here

The folders are used as follows:

- `GISS_orbar_driver.f90` and `GISS_srevents.f90` write orbital-parameter output to the folder `/GISS_orbital`;
- `month_length.f90` reads the info file `month_length_info.csv` in the folder `/info_files` and writes month-length tables to the folder `/month_lengths`;
- `cal_adjust_PMIP.f90` reads the info file `cal_adjust_info.csv` in the folder `/info_files`, and source netCDF files from the folders `/nc_files/PMIP3_source` and `/nc_files/PMIP4_source` and writes paleo calendar-adjusted netCDF files to the folder `/nc_files/PMIP3_adjusted` and `/nc_files/PMIP4_adjusted`;
- `demo_03_adjust_TraCE_ts.f90` reads, for example, the file `TraCE_c30r40_tas_land_monlen0ka_Jan-Dec.csv` in the folder `TraCE_example` and writes a paleo calendar-adjusted output file into the same folder.

Note:  The example source and adjusted netCDF files, as well as those used to create figures, are not included in this repository.  The folder `/nc_files` instead includes instructions on how they may be downloaded from Globus or Dropbox sources.
