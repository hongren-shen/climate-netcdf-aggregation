# Workflow configuration

Edit the defaults at the top of `workflow/aggregate_driver.py`, or override them
with command-line options for one run. Every gauge worker uses the same settings.
Run `python workflow/aggregate_driver.py --help` from the repository root to see
all options and their current defaults.

On/off options come in pairs, such as `--temporal` / `--no-temporal`.
The tables below list the repository defaults, before any local edits.

## Paths and prepared inputs

Default paths use the repository folder, so they also work from an IDE or
another working directory. Relative paths given on the command line start from
the folder where you run the command. Paths below are relative to the repository
root unless stated otherwise.

The four bundled forcing files for `02KA015` cover the complete 2020 calendar
year in UTC, with hourly precipitation and three-hourly temperature. For other
gauges or regions, supply prepared inputs for that domain. Put forcing files in
`<input-root>/<region>/<filename>.nc` and routing files in
`<routing-root>/<region>/<gauge>/`. See [the input guide](inputs.md) for source
and preparation details.

Before clearing outputs, the driver checks the gauge CSV's required columns and
selected gauges. It also checks the forcing paths, required routing files and
shapefile companions, and the grid-weight helper. It does not check NetCDF
contents or routing geometry at startup. If processing fails, the run stops with
a nonzero exit status before packaging that region's outputs and leaves the files
for inspection. ZIP failures also return a nonzero exit status. See
[inputs.md](inputs.md) for input requirements.

## Inputs and gauge selection

| Option | Top-level setting | Default | Meaning |
| --- | --- | --- | --- |
| `--regions REGION ...` | `region_for_process` | `["SE"]` | Region subfolders to process. |
| `--gauges GAUGE ...` | `extract_gauge_list` | `["02KA015"]` | Gauge IDs to select; `--gauges` without IDs selects all matching inventory rows. |
| `--input-root PATH` | `nc_data_folder` | `data/casr_3.1` | Folder containing prepared forcing NetCDFs, organized by region. |
| `--variables FILE ...` | `nc_files_process` | Four files listed below | NetCDF filenames to process in each selected region. |
| `--gauge-table CSV` | `gauge_info_path` | `data/gauge_inventory/SE/24H/gauge_info.csv` | Gauge inventory with `Region` and `Obs_NM` columns. |
| `--routing-root PATH` | `gauge_routing_product_folder` | `data/Routing_product_simplified_extraction` | Folder containing the region/gauge routing subfolders. |

The default forcing filenames are:

```text
CaSR_v3_1_A_PR0_SFC.nc
CaSR_v3_1_A_RAIN_SFC.nc
CaSR_v3_1_A_SNOW_SFC.nc
CaSR_v3_1_A_TT_1_5m.nc
```

The bundled table is the OWL-Hydro SE **24H hydrometric inventory**. Selecting it
does not enable daily forcing aggregation. The driver selects rows by `Region`
and `Obs_NM`; it does not filter by calibration role. It replaces spaces in gauge
names with underscores before matching the requested IDs.

The inventory contains all 120 rows. Bundled routing and forcing cover only
`02KA015`. To select all gauges, supply forcing and routing for every selected
gauge, including both HRU layouts where needed for output weights.

## Spatial aggregation and OWL-Hydro versions

| Option | Top-level setting | Default | Meaning |
| --- | --- | --- | --- |
| `--targets TARGET ...` | `want_aggregated_targets` | `["subbasin"]` | One or more of `basin`, `subbasin`, and `hru`. |
| `--olrrp-version VERSION` | `olrrp_version` | `v2-0-1` | Base routing-product release used in shapefile names. |
| `--hru-versions VERSION ...` | `hru_versions` | `["v2-0-1a", "v2-0-1b"]` | HRU layouts to process when the `hru` target is selected. |
| `--gridded` | `want_gridded_forcing` | `False` | Also export cropped gridded forcing and its grid weights. |

The `basin`, `subbasin`, and `hru` targets correspond to weight Versions 1, 2,
and 3. See [grid-weight versions](../README.md#grid-weight-versions) for how each
version assigns forcing to HRUs and how gridded and aggregated weights differ.

OWL-Hydro uses `v2-0-1a` for the semi-distributed layout and `v2-0-1b` for the
routing-only layout. Their input files are `final_hru_v2-0-1a.shp` and
`final_hru_v2-0-1b.shp`.

The base version in each full HRU code must match `olrrp_version`. The earlier
names `semi-distributed` and `routing-only` are also accepted. Changing a version selects filenames; it
does not convert the routing data. Basin and subbasin outputs include weight
mappings for both HRU layouts, regardless of `--hru-versions`.

## Temporal aggregation

| Option | Top-level setting | Default | Meaning |
| --- | --- | --- | --- |
| `--temporal` | `want_temporal_aggregation` | `False` | Enable temporal aggregation. |
| `--frequency FREQ` | `temporal_agg_freq` | `D` | One output frequency per run: daily (`D`) or a positive whole-hour interval, such as `3H`, `6H`, or `12H`. |
| `--timezone ZONE` | `temporal_tz_name` | `UTC-5` | Clock used for time grouping: UTC, a fixed offset, or a named time zone such as `America/Toronto`. |
| `--anchor-hour HOUR` | `temporal_anchor_hour` | `0` | Local daily boundary hour, from 0 to 23; 0 means midnight. |
| `--period-ending` | `period_ending2starting` | `True` | Shift period-ending timestamps back by one native time step before grouping. |
| `--daily-temp-minmax` | `derive_daily_temp_minmax` | `False` | Add daily temperature MIN and MAX from input time steps. |
| `--label-avg-when-enabled-only` | `label_avg_when_enabled_only` | `True` | Add AVG to daily mean temperature filenames only when MIN/MAX are enabled. Turning this off always adds AVG for daily temperature. |

With temporal aggregation disabled, the driver keeps source timestamps and native
time steps. It applies the timezone, period-ending correction, and daily
temperature options only when temporal aggregation is enabled. UTC-5 is a fixed
offset without daylight saving time; a named time zone can include seasonal clock changes.

Daily aggregation sums the configured precipitation amounts and averages
temperature. With `--daily-temp-minmax`, it also calculates daily MIN and MAX
from the same input time steps. MIN/MAX are produced only for a daily frequency.
The driver accepts the daily aliases `1D` and `24H`, and hourly frequency letters
are case insensitive.

Selecting `1H` does not interpolate native 3-hourly data into a complete hourly
series. The driver still applies the configured timezone and period-ending
correction when temporal aggregation is enabled. The command-line options do not
accept monthly frequencies. Use separate runs and output folders for different frequencies.

## Outputs and execution

| Option | Top-level setting | Default | Meaning |
| --- | --- | --- | --- |
| `--output-root PATH` | `output_folder` | `results` | Folder for outputs; each processed regional folder is cleared first. |
| `--workers N` | `num_cores` | `1` | Number of gauge workers; use 1 to process gauges one at a time. Must be positive. |
| `--debug` | `debug_mode` | `False` | Process the first two selected gauges per region; does not reduce the worker count. |
| `--compress` | `aggregate_compress` | `True` | Compress aggregated NetCDF files. |
| `--compression-level LEVEL` | `aggregate_compress_level` | `1` | NetCDF compression level, from 0 to 9. |
| `--overwrite` | `aggregate_overwrite` | `True` | Allow overwriting individual aggregation files; does not control regional cleanup. |
| `--time-block-records N` | `time_block_hours` | `8760` | Time records per spatial processing block. Must be positive. |
| `--output-time-chunk-records N` | `output_time_chunk_hours` | `8760` | Maximum time records per output chunk. Must be positive. |
| `--zip` | `if_zip_gauge_folder` | `True` | Package each gauge folder as a ZIP archive. |
| `--delete-after-zip` | `if_delete_original_gauge_folder` | `True` | Remove the uncompressed gauge folder after zipping; has no effect with `--no-zip`. |

The two settings whose names end in `_hours` count time records. A block of 8760
records is not always one year of data. Processing blocks and output chunks
control memory use and file layout. They do not set the output frequency.

After startup checks pass, the driver deletes and recreates
`<output-root>/<region>` for each processed region, even with `--no-overwrite`.
Save any outputs you need before rerunning, and keep the output folder separate
from all inputs. Write test runs to a separate folder, such as `results_test/native`.

## Advanced shapefile and grid settings

| Option | Top-level setting | Default | Meaning |
| --- | --- | --- | --- |
| `--subbasin-prefix PREFIX` | `subbasin_name_prefix` | `finalcat_info_` | Prefix before the base version in subbasin shapefile names. |
| `--hru-prefix PREFIX` | `hru_name_prefix` | `final_hru_` | Prefix before the full HRU version in HRU shapefile names. |
| `--outline-name FILE` | `outline_shp_name` | `outline.shp` | Gauge outline filename. |
| `--gridweights-prefix PREFIX` | `grid_weights_file_prefix` | `gridweights_` | Prefix for gridded-product weight files. |
| `--gridweights-generator SCRIPT` | `grid_weights_generator` | `src/derive_grid_weights_modified_202312.py` | Grid-weight helper script. |
| `--crop-buffer-cells N` | `crop_buffer_cell` | `5` | Nonnegative cell padding for the optional gauge crop when `--gridded` is enabled. |

Keep the supplied subbasin prefix and outline name. The writer for aggregated
weight files still uses `outline.shp` and `finalcat_info_<version>.shp` internally.
Those two options do not change the filenames in every part of the workflow.

## Script-only setting

`prefer_precomputed_regional = False` remains at the top of the driver, with no
command-line option. It looks for existing regional aggregations. The driver
clears regional outputs before looking for those files, so enabling this setting
alone does not preserve or reuse outputs from a previous run.

## Examples from the repository root

After installing the dependencies, run the driver on the bundled inputs. The
examples below use the single-gauge data and write to separate test folders.
For another prepared dataset, use `--input-root PATH` to select its input folder.

```bash
# Preserve native timing and retain uncompressed files for inspection.
python workflow/aggregate_driver.py --regions SE --gauges 02KA015 --output-root results_test/native --workers 1 --no-zip

# Process both OWL-Hydro HRU layouts.
python workflow/aggregate_driver.py --regions SE --gauges 02KA015 --targets hru --hru-versions v2-0-1a v2-0-1b --output-root results_test/hru --workers 1 --no-zip

# Create daily subbasin forcing and daily temperature AVG, MIN, and MAX.
python workflow/aggregate_driver.py --regions SE --gauges 02KA015 --targets subbasin --temporal --frequency D --daily-temp-minmax --output-root results_test/daily --workers 1 --no-zip
```
