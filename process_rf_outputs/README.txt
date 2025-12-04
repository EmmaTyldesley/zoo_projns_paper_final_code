zoo_projns_paper_final_code/process_rf_outputs

EJT Dec 2025

This MATLAB code takes the random forest predictions generated in R (R\zoo_projections_train_RF\scripts\main\predict_rf_projns.R) and stored temporarily (temp_rf_output_files) and collates into time-series of gridded monthly projections.

Scripts:
- process_rf_all_steps_monthly.m - high level file which sets model parameters (domain, projection, time periods)
- process_rf_predictions_monthly - function called by the above to load RF predictions

Outputs:
predicted_<RF_model_name>_projn_<projn_model_name>_NEAtlantic_monthly_<decade_start>_<decade_end>.mat - table of date, lon, lat, depth, AMM7 region and predicted zooplankton concentration by taxa (in units of log10(ZE+offset) when offset = 0.5*min(ZE>0)
