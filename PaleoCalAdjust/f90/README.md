## Programs (PaleoCalAdjust v1.0) ##

Main programs, including `month_lenth.f90` and `cal_adjust_PMIP.f90`, plus additional driver and demonstration programs are in the folder `/main_programs`:

	month_length.f90			! month-length tables
	cal_adjust_PMIP.f90			! paleo calendar adjustment
	GISS_orbpar_driver.f90			! orbital parameters (eccentricity, climatic precession)
	GISS_srevents_driver.f90		! equinox, solstice, perihelion dates
	demo_01_pseudo_daily_interp.f90		! demo of pseudo-daily interpolation
	demo_02_adjust_1yr.f90			! demo of paleo calendar adjustment
	demo_03_adjust_TraCE_ts.f90		! demo of transient simulation adjustment
	demo_04_CE_month_lengths.f90	! calculation of month-length variations over the Common Era (on a proleptic Gregorian calendar)

The `/modules` folder contains the following:

	calendar_effects_subs.f90
	CMIP_netCDF_subs.f90
	GISS_orbpar_subs.f90
	GISS_srevents_subs.f90
	month_length_subs.f90
	pseudo_daily_interp_subs.f90

The `/projects` folder contains a set of subfolders, one for each main program, containing GNU Make makefiles for the individual main programs.

The programs are used as follows:

- `GISS_orpar_driver.f90` and `GISS_srevents_driver.f90` write orbital-parameter output to the folder `/GISS_orbital`, using specific parameter values set in the programs;
- `month_length.f90` reads the info file `month_length_info.csv` in the folder `/info_files` and writes month-length tables to the folder `/month_lengths`;
- `cal_adjust_PMIP.f90` reads the info file `cal_adj_info.csv` in the folder `/info_files`, and source netCDF files from the folders `/nc_files/PMIP3_source` and `/nc_files/PMIP4_source` and writes paleo calendar-adjusted netCDF files to the folder `/nc_files/PMIP3_adjusted` and `/nc_files/PMIP4_adjusted`;
- `demo_01_pseudo_daily_interp.f90`, `demo_02_adjust_1yr.f90`, and `demo_04_CE_month_lengths.f90` are stand-alone programs, writing only to the console;
- `demo_03_adjust_TraCE_ts.f90` reads, for example, the file `TraCE_c30r40_tas_land_monlen0ka_Jan-Dec.csv` in the folder `TraCE_example` and writes a paleo calendar-adjusted output file into the same folder.

Most of the programs write progress information to the console.

With the exception of directory paths, the same code compiles and runs on the following systems:

- Windows 10: Intel Parallel Studio XE 2019 Update 4 Composer Edition for Fortran Windows, with netCDF version 4.1.3, using the Visual Studio 2019 IDE; 
- MacOS: gfortran version 9.1.0, with netCDF version 4.6.3 (from Homebrew), using the Eclipse IDE for Scientific Computing, Version 4.12.0 (2019-06), with Eclipse for Parallel Application Developers, Version 9.2.1.