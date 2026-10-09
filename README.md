# ifs_aifs_windens_comparison

Repository for the paper "A global evaluation of ECMWF's operational IFS and AIFS wind speed ensemble forecasts"

## Overview
This repository provides R codes used for fitting EMOS and quantile regression post-processing models to IFS and AIFS wind speed ensemble forecasts.

## Scripts
- `DualModelEmos_functions.R`: Functions used for fitting EMOS models
- `DualModelEmos_fit.R`: Main script to fit local EMOS models
- `DualModelQR_functions.R`: Functions used for fitting quantile regression models
- `DualModelQR_fit.R`: Main script to fit semi-local quantile regression models

## Other files
- `XY.csv`: Information of SYNOP stations used in the study

## Data availability
The wind speed ensemble and observation data used in the study is not provided due to data sharing policies.

## Input data structure of post-processing modeling codes

The post-processing modeling codes require 4-dimensional arrays with the following dimensions:
1. **Date:** Forecast initialization date
2. **Station:** SYNOP station
3. **Variable:** Wind speed ensemble summary statistics (mean, variance, 0.1-, 0.5- and 0.9-quantiles) and the corresponding observation:
   - Ensemble mean (`EnsMean`)
   - Ensemble variance (`EnsVar`)
   - Quantiles: 10th, 50th and 90th percentiles (`EnsQ10`, `EnsQ50`, `EnsQ90`)
   - Observation (`obs`)
4. **Lead time:** Forecast lead time in hours
