PaleoCalAdjust
===================

This is the repository that accompanies the paper:

Bartlein, P. J. and Shafer, S. L.: Paleo calendar-effect adjustments in time-slice and transient climate-model simulations (PaleoCalAdjust v1.0): impact and strategies for data analysis, *Geosci. Model Dev.*,  [https://doi.org/10.5194/gmd-2018-283](https://doi.org/10.5194/gmd-2018-283), 2019.

### Abstract ###

The “paleo calendar effect” is a common expression for the impact that the changes in the length of months or seasons over time, related to changes in the eccentricity of Earth’s orbit and precession, have on the analysis or summarization of climate-model output. This effect can have significant implications for paleoclimate analyses. In particular, using a “fixed-length” definition of months (i.e. defined by a fixed number of days), as opposed to a “fixed-angular” definition (i.e. defined by a fixed number of degrees of the Earth’s orbit), leads to comparisons of data from different positions along the Earth’s orbit when comparing paleo with modern simulations. This effect can impart characteristic spatial patterns or signals in comparisons of time-slice simulations that otherwise might be interpreted in terms of specific paleoclimatic mechanisms, and we provide examples for 6, 97, 116, and 127 ka. The calendar effect is exacerbated in transient climate simulations, where, in addition to spatial or map-pattern effects, it can influence the apparent timing of extrema in individual time series and the characterization of phase relationships among series. We outline an approach for adjusting paleo simulations that have been summarized using a modern fixed-length definition of months and that can also be used for summarizing and comparing data archived as daily data. We describe the implementation of this approach in a set of Fortran 90 programs and modules (PaleoCalAdjust v1.0).

## Overview ##

The PaleoCalAdjust programs implement an approach for adjusting paleoclimate model output for the "paleo calendar effect" -- the impact that the changes in the length of months or seasons over time (related to the changes in the eccentricity of Earth's orbit and to precession), have on the summarization of the model output.  The current version of the key program is `cal_adjust_PMIP.f90` (in the folder `/f90/main_programs`), which applies the paleo calendar effect adjustment to CMIP5/PMIP3- or CMIP6/PMIP4-formatted files.  There is a related program, `month_length.f90`, that can be used to produce tables of the changing length of months over time. (The original main program `cal_adjust_PMIP3.f90` was modified to accommodate CMIP6/PMIP4-formatted files and renamed.)  Figures illustrating the paleo calendar effect are in the folder `/figures`, and relevant data sets for exercising the programs are in the folder `/data`.  

Several minor modifications to the main program `cal_adjust_PMIP.f90` and its modules were made since the original *Geoscientific Model Development Discussions (GMDD)* manuscript submission to accommodate the adjustment of CMIP6-PMIP4 files, the filenames of which contain an additional "grid_label" field not present in CMIP5-PMIP3 filenames.  Additionally, following a referee's suggestion, we replaced the approach for calculating month lengths using the approximation of Kutzbach and Gallimore (1988, *J. Geophys. Res.* 93(D1):803-821), with a direct approach based on Kepler's equation. This substitution of approaches had no practical significance.  Several other code modifications were made in the interests of transparency.

This release (v1.0d) is consistent with the revised and accepted version of the paper in the GMDD discussion, but has been modified to also handle 4-D files (e.g. *lon* x *lat* x *level* x *time*), and rotated-pole files (e.g. *Oclim* or *OIclim* files).

The main changes from the original submission therefore include:

- the main program `cal_adjust_PMIP3.f90` and module `CMIP5_netCDF_subs.f90` were renamed as `cal_adjust_PMIP.f90` and `CMIP_netCDF_subs.f90` respectivly, because they are now generic.  No changes were made to `CMIP_netCDF_subs.f90`, and only minor changes were made to `cal_adjust_PMIP.f90` to allow reading and writing of both CMIP5-PMIP3 and CMIP6-PMIP4 formatted files;
- the info file read by `cal_adjust_PMIP.f90` (e.g. `cal_adj_info.csv`) was modified to include an "activity" field (either "PMIP3" or "PMIP4") and a "grid_label" field (blank for PMIP3 files);
- the paths to the "source" and "adjusted" netCDF files are now explicitly given in the info file, as opposed to being set in the main program;
- the subroutine `monlen(...)` in the module file `month_length_subs.f90` now computes month length, beginning, middle and ending days using an approach based on Kepler's equation.  To view the code changes and their impact on the results, the code and figures (plus the data used to plot the figures) in this release can be compared with those in previous releases (e.g. v1.0b);
- The main program `cal_adjust_PMIP.f90` was modified to allow the adjustment of 4-D data sets and generalized to work with "rotated-pole" data sets (by eliminating explicit references to longitude and latitude).  This change required the addition of a subroutine `get_var_diminfo(...)` to `CMIP_netCDF_subs.f90` for getting dimension-variable information to `CMIP_netCDF_subs.f90`     

## Version history ##

### v1.0 ###

Original release.

### v1.0a ###

Minor modifications for consistency with the Bartlein and Shafer (2018, *Geoscientific Model Development Discussions*) submission.

### v1.0b ###

Several minor modifications to the main program `cal_adjust_PMIP.f90` and its modules were made to accommodate CMIP6/PMIP4 files, the filenames of which contain an additional "grid_label" field not present in CMIP5/PMIP3 filenames.  

- The main program `cal_adjust_PMIP3.f90` and the module `CMIP5_netCDF_subs.f90` were renamed as `cal_adjust_PMIP.f90` and `CMIP_netCDF_subs.f90` respectively, because they now accommodate CMIP5/PMIP3- and CMIP6/PMIP4-formatted files.  No changes were made to `CMIP_netCDF_subs.f90`, and only minor changes were made to `cal_adjust_PMIP.f90` to allow reading and writing of both CMIP5/PMIP3- and CMIP6/PMIP4-formatted files;

- The info file read by `cal_adjust_PMIP.f90` (e.g. `cal_adj_info.csv`) was modified to include an "activity" field (either "PMIP3" or "PMIP4") and a "grid_label" field (blank for PMIP3 files), and the paths to the "source" and "adjusted" netCDF files are now explicitly given in the info file, as opposed to being set in the main program.

### v1.0c ###

Following a referee's suggestion, we replaced the approach for calculating month lengths using the approximation of Kutzbach and Gallimore (1988, *J. Geophys. Res.* 93(D1):803-821), with a direct approach based on Kepler's equation.  This substitution of approaches had no practical significance.  Several other code modifications were made in the interests of transparency.

- The subroutine `monlen(...)` in the module file `month_length_subs.f90` now computes month length, beginning, middle and ending days using an approach based on Kepler's equation.  To view the code changes and their negligible impact on the results, the code and figures (plus the data used to plot the figures) in this release can be compared with those in previous releases (e.g. v1.0b).    

### v1.0d ###

- The main program `cal_adjust_PMIP.f90`, was modified to handle rotated-pole and 4-dimensional (e.g. time x level x latitude x longitude) files.  This change required development of an additional subroutine, `get_var_diminfo(...)` in the module `CMIP_netCDF_subs.f90`. 

