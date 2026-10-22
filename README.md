# climate-netcdf-aggregation

Utilities and scripts for aggregating gridded NetCDF meteorological forcing data to **gauge/catchment polygons** (spatial aggregation), and optionally **resampling in time** (temporal aggregation). The typical output is one NetCDF (and optional ASCII metadata) per gauge, suitable for hydrologic modelling workflows.

## What’s included

- `aggregate_driver.py`  
  
  Example/driver script to run the end-to-end workflow across many gauges.
  
- `spatial_aggregation.py`  
  
  Core spatial aggregation routines (grid to polygon weighting and aggregation).
  
- `temporal_aggregation.py`  
  
  Temporal aggregation utilities (e.g., hourly to daily), including accumulation vs mean handling.
  
- `forcing_utils.py`  
  
  Shared helper functions (naming, metadata export, zipping outputs, convenience wrappers).
  
- `derive_grid_weights_modified_202312.py`  
  
  Helper script to derive grid-cell weights for a target polygon geometry (modified from Juliane Mai's gridweights generator; added more outputs for aggregation use).

## Requirements

This repo is written in Python 3.11 (recommended: 3.11/3.12/3.13).

Typical dependencies:
- `numpy`
- `pandas`
- `xarray`
- `netCDF4`
- `tqdm` (optional but recommended)

- `geopandas`
- `shapely`
- `pyproj`
