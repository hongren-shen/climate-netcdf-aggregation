# climate-netcdf-aggregation

Python tools for area-weighted aggregation of prepared gridded NetCDF
meteorological forcing to gauge basins, subbasins, and hydrologic response
units (HRUs). The workflow writes aggregated NetCDF files, Raven grid weights,
and text metadata, with optional temporal aggregation and ZIP packaging.

The example uses CaSR v3.1 forcing and OWL-Hydro routing products. The forcing
was downloaded and merged into one NetCDF file per variable, precipitation was
partitioned into rain and snow, and NetCDF attributes were harmonized before
the example was extracted. The aggregation script starts from these prepared files.

Code for CaSR download, rain/snow partitioning, and NetCDF harmonization is
planned for future updates.

## Setup

Use Python 3.11. The tests and example spatial and temporal aggregations were
run with Python 3.11.14. Package versions are listed in
[requirements.txt](requirements.txt) and [environment.yml](environment.yml).

From the repository root, create and activate the environment:

```bash
conda env create -f environment.yml
conda activate climate-netcdf-aggregation
python workflow/aggregate_driver.py --help
```

GDAL and its `osgeo` Python bindings are required by the grid-weight helper.
The environment file installs them with the scientific Python packages from
conda-forge. The helper runs through `python`, so keep this environment activated
when running the driver, including from an IDE. See the
[GDAL installation instructions](https://gdal.org/en/stable/api/python/python_bindings.html#installation)
for other installation methods. A pip installation from `requirements.txt`
requires a compatible native GDAL installation.

SciPy speeds up aggregation using sparse matrices and is included in the
environment. The code uses dense matrices if SciPy is unavailable. The driver
does not require ArcGIS or NCO.

## Included SE example

The repository includes an example for gauge **02KA015** in southeastern
Ontario (`SE`): routing shapefiles, four prepared gridded forcing files, and
reference results in a ZIP file. The forcing covers the full **2020 calendar
year in UTC**, with hourly total precipitation, rainfall, and snowfall, and
three-hourly temperature. It includes all grid cells intersecting the gauge
basin, plus a buffer of five grid cells on each side. The resulting grid is
13 by 13 cells, and the four input files total approximately 7 MB.
See [the input guide](docs/inputs.md#bundled-2020-forcing-example)
for the original grid indices and data preparation steps.

The default gauge table is an unchanged copy of the 120-row OWL-Hydro SE
**24H** inventory. The default settings select only `02KA015`. The driver
selects gauges using `Region` and `Obs_NM`; it does not filter by calibration
role. The included forcing and routing files cover only the example gauge.

```text
climate-netcdf-aggregation/
|-- CITATION.cff
|-- LICENSE
|-- README.md
|-- requirements.txt
|-- environment.yml
|-- docs/
|   |-- configuration.md
|   `-- inputs.md
|-- data/
|   |-- gauge_inventory/SE/24H/gauge_info.csv
|   |-- Routing_product_simplified_extraction/SE/02KA015/
|   `-- casr_3.1/SE/                 # four prepared 2020 forcing files
|-- results/SE/02KA015.zip
|-- src/
|-- tests/
`-- workflow/aggregate_driver.py
```

The four forcing files included under `data/casr_3.1/SE/` are:

```text
CaSR_v3_1_A_PR0_SFC.nc
CaSR_v3_1_A_RAIN_SFC.nc
CaSR_v3_1_A_SNOW_SFC.nc
CaSR_v3_1_A_TT_1_5m.nc
```

After installing the dependencies, run the example from the repository root:

```bash
python workflow/aggregate_driver.py
```

This recreates `results/SE/02KA015.zip`. To keep the individual output files
in a separate folder:

```bash
python workflow/aggregate_driver.py --regions SE --gauges 02KA015 --targets subbasin --workers 1 --no-zip --output-root results_test
```

This writes subbasin forcing to `results_test/SE/02KA015/`, keeping the source
timestamps and time steps. The example ZIP remains in `results/`.

Default paths are relative to the repository, so the script can be run from
another directory. Relative paths supplied on the command line are interpreted
from the directory where you run the command.

## Inputs for another region

To process another region, prepare its forcing and routing data, then set
`--regions`, `--gauges`, `--gauge-table`, `--input-root`, and `--routing-root`
to match those files. The included forcing covers only the area around `02KA015`.
The driver expects this folder structure:

```text
<input-root>/<REGION>/<VARIABLE>.nc
<routing-root>/<REGION>/<GAUGE>/outline.shp
<routing-root>/<REGION>/<GAUGE>/finalcat_info_v2-0-1.shp
<routing-root>/<REGION>/<GAUGE>/final_hru_v2-0-1a.shp
<routing-root>/<REGION>/<GAUGE>/final_hru_v2-0-1b.shp
```

Each shapefile needs its companion files and a defined coordinate reference
system (CRS). Gauge IDs must match `Obs_NM` and the gauge folder names. Region
values must match the input folder names. NetCDF files must use the expected
grid coordinates, time attributes, and units. See [input requirements](docs/inputs.md)
for details and the corresponding folders in the OWL-Hydro archive.

## Configuration

You can edit the settings at the top of `workflow/aggregate_driver.py` or pass
options when running it in a terminal. Command-line options change the settings
for that run, including all parallel processes. They do not change the defaults
saved in the script. Use `--help` to see the available options, or read the
[parameter reference](docs/configuration.md) for definitions and examples.

| Setting | Default behavior |
| --- | --- |
| Gauge inventory | OWL-Hydro SE 24H table; only `02KA015` selected |
| Spatial target | Subbasin aggregation |
| HRU layouts | `v2-0-1a` and `v2-0-1b`, when the `hru` target is selected |
| Temporal aggregation | Disabled; preserve source timestamps and native time steps |
| Daily temperature MIN/MAX | Disabled |
| Parallel workers | 1; increase for a larger gauge selection when resources permit |
| Output packaging | One ZIP per gauge; remove the uncompressed gauge folder after zipping |

The HRU suffixes follow the **OWL-Hydro convention**: `v2-0-1a` is the
semi-distributed layout and `v2-0-1b` is the routing-only layout. The base
version must match `--olrrp-version`. The earlier names `semi-distributed` and
`routing-only` are also accepted.

The **24H gauge inventory does not set the forcing output frequency**. To
create daily forcing, enable temporal aggregation explicitly:

```bash
python workflow/aggregate_driver.py --gauges 02KA015 --targets subbasin --temporal --frequency D --daily-temp-minmax --workers 1 --no-zip --output-root results_daily
```

This command processes all four forcing files using the default fixed UTC-5
time zone. It shifts period-ending timestamps back by one source time step
and groups values into days starting at midnight. Because the inputs cover a
calendar year in UTC, the first and last days in UTC-5 have fewer than 24 hours
of data. Check these partial days before using the results in a model.

## Processing rules

- Spatial aggregation uses grid-cell and polygon overlap weights. The targets
  are a complete gauge basin (`basin`), individual subbasins (`subbasin`), or
  individual HRUs (`hru`). Source units are preserved.
- With temporal aggregation disabled, each variable keeps its original time
  step: hourly precipitation and three-hourly temperature in the CaSR example.
- With temporal aggregation enabled, the driver uses the filename to choose
  whether to sum or average values. It sums variables recognized as
  precipitation and averages temperature and other variables. Inputs must
  already contain the required interval amounts or temperatures. The driver
  does not convert precipitation rates to amounts or Kelvin to Celsius.
- With the default temporal settings, period-ending timestamps are shifted
  back by one source time step, inferred from the time coordinate. Times are
  converted to fixed UTC-5 without daylight saving time. Each day starts at
  the configured local hour, and its output timestamp marks that start.
- Daily temperature MIN and MAX are calculated from the spatially aggregated
  subdaily series. They are not spatial averages of grid-cell daily extrema.
- Each run selects one output frequency. A finer requested frequency does not
  interpolate missing time steps or create hourly temperature from a
  three-hourly source.

Use inputs on matching grids that cover all selected basins, subbasins, or HRUs
and have regular timestamps. The code does not fill data gaps, require a minimum number
of records per output period, or supply weights for polygons with no grid
coverage. See [the input guide](docs/inputs.md) before using another dataset.

## Outputs

Without ZIP packaging, generated files are written to
`<output-root>/<REGION>/<GAUGE>/`. They include:

- Aggregated NetCDFs with dimensions `time` and `stations`. The `station_id`
  variable contains basin, subbasin, or HRU IDs, not hydrometric gauge IDs.
- Raven grid-weight files linking model HRUs to the aggregated series.
- `netcdf_metadata_*.txt` files describing dimensions, variables, and attributes.
- Gridded NetCDFs cropped to the gauge area, with their grid weights, when
  `--gridded` is enabled.

When writing grid-weight headers, the code reads the software version from
`CITATION.cff`. Update its `version` field for a new release. If the citation file
is missing or its version cannot be read, the header shows only the repository
name and link. Existing output files retain their original headers.
The forcing label comes from the referenced NetCDF's global `title` attribute,
or `source` when no title is provided. It is omitted if both are unavailable.
The header also records the referenced NetCDF filename. Weights shared across
variables refer to the first processed file used to generate them.

Filenames identify the variable, gauge, spatial target, routing version, date
range, time step, and time zone. Daily temperature statistics add AVG, MIN,
or MAX labels when requested. With the default ZIP settings, each gauge folder
is packaged as `<output-root>/<REGION>/<GAUGE>.zip` and then removed.

The reference file `results/SE/02KA015.zip` was generated from the included
2020 inputs using the default settings: subbasin aggregation, original time
steps, and one worker process. It contains ten files: four aggregated NetCDFs,
four text files describing their contents, and two Raven grid-weight files
for the `v2-0-1a` and `v2-0-1b` HRU layouts. The aggregated series keep the
source UTC timestamps. The original gridded inputs are in `data/casr_3.1/SE/`.

The driver checks the gauge table, matching gauge selection, forcing files,
required routing files and shapefile companions, and grid-weight helper before
changing outputs. After these checks pass, **each selected regional output
directory is cleared**, including with `--no-overwrite`. Use a separate folder
for outputs and save any previous results you need. If processing or ZIP
creation fails, the run stops and returns a nonzero exit status to indicate
an error. If processing fails, the current region is not packaged into ZIP
files. Partial results are kept so you can inspect them.

### Grid-weight versions

The `version1`, `version2`, and `version3` prefixes in weight filenames show
which area is averaged to assign forcing to model HRUs:

| Weight version | Aggregated target | How forcing is assigned to HRUs |
| --- | --- | --- |
| Version 1 | `basin` | Average over the complete gauge basin. All HRUs in that basin share the same weights and forcing series. |
| Version 2 | `subbasin` | Average over each subbasin. HRUs with the same `SubId` share its weights and forcing series. |
| Version 3 | `hru` | Average over each HRU polygon. Each HRU uses weights derived from its own geometry. |

For aggregated NetCDFs, `--targets basin`, `--targets subbasin`, and `--targets hru`
select Versions 1, 2, and 3, respectively. These `*_agg_*` weight files use
weight `1.0` to select the appropriate already averaged NetCDF series for each
HRU. The default example uses Version 2.

With `--gridded`, the driver also creates all three weight versions for each
available HRU layout (`a` and `b`). These files contain area-overlap weights
referencing grid cells. Use each weight file with the NetCDF named in its header.

These weight version numbers are separate from the software and routing
versions. For example, `version2_gridweights_agg_subbasin_v2-0-1a.txt` uses
subbasin averaging (Version 2), routing release `v2-0-1`, and the semi-distributed
HRU layout `a`. The `Created by` header separately records the software version,
such as `v1.0.0`.

## Structure and checks

```text
workflow/   main script and command-line options
src/        spatial and temporal aggregation, metadata, and grid weights
tests/      checks for settings, errors, metadata, and example data
```

From the repository root, with the environment activated:

```bash
python -m unittest discover -s tests -v
python src/derive_grid_weights_modified_202312.py --help
```

The tests check command-line options, HRU names, the gauge table, input checks
before output deletion, settings passed to parallel processes, and errors
during processing or ZIP creation. They also check software and forcing
metadata, preservation of packed values and missing-value codes during time
cropping, and the example files' grids, time steps, and full 2020 coverage.
Run the example to check the complete workflow. Expected files are described
in [Outputs](#outputs).

## Citation and data sources

### Software

Software author: Hongren Shen. Citation metadata for version 1.0.0 is provided
in [CITATION.cff](CITATION.cff). When using GitHub code, record the exact commit.
For an archived software release, cite the DOI of the version used. The
repository citation is:

Shen, H. (2026). *Climate NetCDF aggregation tools* (Version 1.0.0)
[Computer software]. [GitHub](https://github.com/hongren-shen/climate-netcdf-aggregation).

### Grid-weight generator

The workflow uses a modified copy of Juliane Mai's
[Grid Weights Generator for Raven](https://github.com/julemai/GridWeightsGenerator).
Please acknowledge the original tool when using this workflow:

Mai, J. *Grid Weights Generator for Raven* [Computer software].
[GitHub](https://github.com/julemai/GridWeightsGenerator).

The included `src/derive_grid_weights_modified_202312.py` calculates weights
from the area of overlap between grid cells and polygons. Our modification
adds a table of polygon IDs, grid-cell indices, row and column positions, and
weights. `spatial_aggregation.py` uses this table to aggregate the forcing
time series. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for the
original source, license notice, and changes to the script.

### OWL-Hydro dataset

For the OWL-Hydro gauge inventory, routing configurations, or prepared data
products used with these tools, please cite:

Shen, H., Tolson, B., Han, M., Craig, J., Metcalfe, R. (2026).
*OWL-Hydro v1.0: An Integrated, Model-Ready Dataset for Semi-Distributed
Hydrologic Modelling in Ontario, Canada*. Federated Research Data Repository.
[https://doi.org/10.20383/103.01878](https://doi.org/10.20383/103.01878)

### Routing and meteorological sources

The routing inputs are based on OLRRP v2. Please also cite:

Shen, H., Han, M., Tolson, B., and Metcalfe, R. A. (2024).
*The Ontario Lake and River Routing Product Version 2 (OLRRP V2)* [Data set].
[https://doi.org/10.5281/zenodo.14225103](https://doi.org/10.5281/zenodo.14225103)

Acknowledge Environment and Climate Change Canada for CaSR. Record the actual
source version, variables, units, time convention, and data preparation steps.
The [ECCC CaSR documentation](https://hpfx.collab.science.gc.ca/~scar700/rcas-casr/dataset_specifics.html)
describes the dataset. This repository's example uses CaSR v3.1.

## License and acknowledgements

Code is provided under the [GNU Lesser General Public License, version 3 or
later](LICENSE). The modified grid-weight generator retains Juliane Mai's
copyright and license notices. See [third-party notices](THIRD_PARTY_NOTICES.md).

The included OWL-Hydro gauge inventory is provided under CC BY 4.0, separately
from the code license. Routing data and the CaSR-derived forcing and outputs
retain their source attribution and applicable data terms. See the
[input guide](docs/inputs.md) for the gauge table and forcing data sources, and cite
the dataset records listed above.
