# Input requirements

The driver uses prepared regional forcing and routing layers extracted for
individual gauges. The bundled example is `02KA015` in southeastern Ontario
(`SE`), with CaSR v3.1 forcing and OLRRP `v2-0-1` routing products. Paths below are
relative to the repository unless an input root is given.

## Gauge inventory

The default table, `data/gauge_inventory/SE/24H/gauge_info.csv`, contains all
120 station records from the OWL-Hydro SE 24H inventory. The default gauge list
selects only `02KA015`.

The source file within OWL-Hydro v1.0 is
`data/regions/SE/04_gauge_inventory/24H/gauge_info.csv`. The table is copied
without changes and is licensed under CC BY 4.0. This data license is separate
from the software license. Cite the
[OWL-Hydro dataset](../README.md#owl-hydro-dataset) when using it.

| Field | Driver use |
| --- | --- |
| `Region` | Matches `--regions` and the regional input folder. |
| `Obs_NM` | Gauge identifier, matched against `--gauges` and the routing folder name. |

The other columns retain the source station information. The driver does not
filter on `cal_target`, `data_type`, or fields describing data availability.
Keep gauge identifiers as text, including leading zeros. Matching is case
sensitive; spaces in `Obs_NM` are replaced by underscores.

For another regional inventory, pass its `gauge_info.csv` to `--gauge-table`.
The 1H or 24H inventory supplies station information. Set the forcing output
time step separately with `--temporal` and `--frequency`.

Use `--gauges` without IDs to select all matching rows. This repository includes
the routing folder and forcing subset for `02KA015` only. Some inventory stations
do not have a separate gauge package. Provide routing files and forcing data
that cover every selected basin.

The public OWL-Hydro archive divides SW into SWN and SWS, while source
inventories can retain `Region = SW`. The driver does not translate these
labels. For SWN and SWS, edit a working copy of the table so its region values
match the forcing and routing folder names.

## Routing products

Place routing layers under `<routing-root>/<REGION>/<GAUGE>/`. The default
shapefile names are:

| File | Required identifier fields | Use |
| --- | --- | --- |
| `outline.shp` | `ID` | Basin aggregation and optional gridded crop. |
| `finalcat_info_v2-0-1.shp` | `SubId` | Subbasin aggregation. |
| `final_hru_v2-0-1a.shp` | `HRU_ID`, `SubId` | OWL-Hydro semi-distributed HRUs and their forcing mappings. |
| `final_hru_v2-0-1b.shp` | `HRU_ID`, `SubId` | OWL-Hydro routing-only HRUs and their forcing mappings. |

Keep `.shx`, `.dbf`, `.prj`, and any supplied `.cpg` files alongside each `.shp`.
Geometries must be valid and have a defined CRS. Identifiers must match the
model routing products. Basin and subbasin outputs also need HRU layers to
write their weight mappings, even when `--targets hru` is not selected.

OWL-Hydro stores routing layers inside each gauge package. Extract the needed
layers into the layout above before setting `--routing-root`. The driver reads
this folder structure directly and does not rearrange the archive. The bundled
routing folder also contains example shapefiles that the driver does not use.

## Prepared forcing NetCDFs

The driver expects `<input-root>/<REGION>/<FILENAME>`. The default root is
`data/casr_3.1/`. All four files below are bundled under its `SE/` directory:

| Filename | Prepared quantity | Example units | Input time step |
| --- | --- | --- | --- |
| `CaSR_v3_1_A_PR0_SFC.nc` | Total precipitation amount | `kg m-2` | 1 hour |
| `CaSR_v3_1_A_RAIN_SFC.nc` | Rainfall amount from precipitation partitioning | `kg m-2` | 1 hour |
| `CaSR_v3_1_A_SNOW_SFC.nc` | Snowfall water-equivalent amount from precipitation partitioning | `kg m-2` | 1 hour |
| `CaSR_v3_1_A_TT_1_5m.nc` | Air temperature at 1.5 m | `degC` | 3 hours |

### Bundled 2020 forcing example

The four files were extracted from the prepared SE regional files in `SE.zip`.
They cover the full 2020 calendar year in UTC, with hourly precipitation,
rainfall, and snowfall, and three-hourly temperature. The source timestamps
still mark the end of each interval. The intervals span 1 January 2020 00:00 UTC
through 1 January 2021 00:00 UTC. The precipitation files contain 8,784 records,
labelled from 1 January 2020 01:00 UTC to 1 January 2021 00:00 UTC. Temperature
contains 2,928 records, labelled from 1 January 2020 03:00 UTC to the same final
timestamp. The January 2021 label closes the last interval of 2020.

The source regional grid has 70 by 70 cells. Seven cells overlap the `02KA015`
basin, within zero-based source rows 51 through 53 and columns 11 through 13.
Adding five grid cells on each side of that bounding rectangle gives the
Python slices `rows[46:59]` and `columns[6:19]`, producing a 13 by 13 grid.
The subset preserves source coordinates and values, with no spatial resampling
or temporal interpolation. Its coverage is limited to the example basin and
buffer. Check spatial coverage before selecting other inventory gauges; the
subset does not cover all of SE.

Before extraction, the source data had been downloaded and merged into one
NetCDF per variable. Precipitation had been partitioned into rain and snow,
and the NetCDF attributes had been standardized. Acknowledge Environment and
Climate Change Canada for CaSR v3.1, along with the dataset citations in the
[main README](../README.md#citation-and-data-sources).

CaSR download, rain/snow partitioning, and NetCDF harmonization workflows are
planned additions to the repository. The current driver does not perform
these steps or convert units. Prepare other input products with the intended
units, names, and time conventions before running the driver. Use `--variables`
to select which supplied files to process.

### Grid and time coordinates

The grid-weight wrapper uses `rlon,rlat` as the spatial dimensions and `lon,lat`
as the geographic coordinate variables. In the CaSR example:

- Meteorological variables use `(time, rlat, rlon)`.
- `rlat` and `rlon` locate the rotated grid; two-dimensional `lat` and `lon`
  give the cells' geographic coordinates.
- `time` has CF-compatible numeric units and a Gregorian-compatible calendar,
  with ordered, regular timestamps. Temporal processing converts dates to Python
  datetimes, so non-Gregorian dates such as 30 February are not supported.
  The configured temporal workflow expects source UTC times.
- Each file covers the complete spatial extent needed by the selected gauges.

The core spatial routines can read other dimensions. The driver and grid-weight
wrapper use the CaSR layout above, so check and update coordinate handling
when using another grid.

## Time and data preparation

Spatial aggregation preserves source units and timestamps. For temporal
aggregation, precipitation-related filename patterns select a sum; temperature
and other variables use a mean. Check which calculation the filename selects
before processing another variable, especially a rate or accumulated field.

The default period-ending correction applies only when temporal aggregation
is enabled. Set `--no-period-ending` for source timestamps that already mark
the start of each interval. The default `UTC-5` is a fixed offset without
daylight saving time. Daily output timestamps mark the start of each day.

Use regular time series with appropriate missing-value encoding. Spatial
aggregation does not renormalize weights when cells have missing values.
Temporal aggregation does not require a minimum number of records per output
period and does not fill gaps.
The first and last output periods can be incomplete. Check missing values,
gaps, and these edge periods before using the results in a model.

Grid weights must cover the selected polygons. The workflow does not assign
nearest-neighbour weights to uncovered polygons. Use valid input geometries
and check that all expected spatial units appear in the output.

## Startup checks and common errors

Before clearing regional outputs, the driver checks the gauge CSV and required
columns, finds matching stations, and checks the forcing files, routing files,
shapefile companions, and grid-weight helper. These startup checks do not
validate all NetCDF values or routing geometries. If processing fails, the run
returns a nonzero exit status and stops before creating ZIP files for the
current region. Partial outputs are kept for inspection. ZIP failures also
return a nonzero exit status.

| Message or symptom | Action |
| --- | --- |
| Gauge table not found | Set `--gauge-table` to an existing regional inventory. The bundled default is the SE 24H table. |
| No gauges match | Check the table's `Region` and `Obs_NM` values against the selected region and gauge IDs. |
| Prepared forcing NetCDF files not found | Restore the four bundled files for the example, or set `--input-root` to a directory containing the region folder and selected files. |
| Required routing files not found | Set `--routing-root` to matching region/gauge folders and retain the required `.shp`, `.shx`, `.dbf`, and `.prj` files. |
| Grid-weight helper not found | Restore the bundled script or set `--gridweights-generator` to an existing copy. |
| `No module named osgeo` or grid-weight helper fails | Activate the environment containing GDAL and confirm that its `python` is used by the driver and helper. |
| No output for a target | Check the required shapefile, its companion files, IDs, geometry, grid coverage, and the processing messages. |

See [configuration](configuration.md) for the full parameter list and the
[main README](../README.md#citation-and-data-sources) for dataset citations.
