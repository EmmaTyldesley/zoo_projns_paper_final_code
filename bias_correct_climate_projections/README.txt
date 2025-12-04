zoo_projns_paper_final_code/bias_correct_climate_projections

EJT Dec 2025

This MATLAB code reads monthly projection data for now, near future and far future decadal slices and calculates a monthly climatology. These are used to calculate laterally gridded monthly relative change Delta values. For physical variables, Deltas are calculated additively. For biogeochemical variables, Deltas are calculated multiplicatively since variables must be non-negative. The Deltas are applied to the AMM7 reanalysis to get bias-corrected projections, i.e. bias corrected future = AMM7 now +(or *) Delta. The bias corrected predictors are saved in MATLAB struct format and converted to timestep tables to be read into R as data frames for the random forest models to predict on. Also includes MATLAB scripts to make manuscript plots and tables.

-- Scripts:

Calculating Delta values
- bias_correction_bgc_<projn_model_name>.m - derives Deltas for bgc variables
- bias_correction_phys_<projn_model_name>.m - for phys variable
- bloom_phenology_projns.m - for RECICLE phytoplankton bloom dynamics
Note: calculated separately because of grid differences between the projection models but could be generalised.
- combine_monthly_Deltas.m - combines Deltas by variable and model into one structure

Applying Delta values
apply_Delta_monthly_all_projns.m - loads monthly AMM7 reanalysis and applies Delta change values for each projection model

Make plots/tables for the paper:
- figures_and_tables\plot_maps_bias_corrected_predictors.m
- figures_and_tables\plot_maps_Deltas.m
- figures_and_tables\table_predictor_Deltas_by_region.m

-- Inputs:
- Monthly mean climate projections & mesh files. Projections data can be obtained as outlined below. 
	HADGEM - Zenodo at https://zenodo.org/records/3953801
	RECICLE IPSL - BODC at https://doi.org/10.5285/0786d770-ba54-0e57-e063-6c86abc09fdf
	RECICLE GFDL - BODC at https://doi.org/10.5285/07877700-0e22-5fb5-e063-6c86abc03058
Note that the RECICLE data are only hosted as *daily* outputs. Mean monthly data were derived using NCO's ncra command (https://nco.sourceforge.net/nco.html). Example bash script included here.

- Daily RECICLE chl files. See links above.

-- Outputs:
Delta change values
- Delta_CHL_PFT_NO3_<projn_model_name> - bgc Deltas
- Delta_SST_NBT_SSS_MLD_<projn_model_name> - phys Deltas
- bloom_phenology_<projn_model_name> - phytoplankton bloom start and duration
- Delta_bloom_updated.mat - bloom Deltas

Bias corrected predictors
MATLAB format
- bias_corrected_predictors_bloom.mat - bloom dynamics obtained by applying Deltas to reanalysis data
- bias_corrected_<projn_model_name>_predictors_incl_no3.mat
Textfile format for reading into R as dataframe
- projn_timestep_<projn_model_name>_mo_NEAtlantic_<YYYYMM>.txt - table of predictors, one file per month, grid-cells in rows rather than gridded






