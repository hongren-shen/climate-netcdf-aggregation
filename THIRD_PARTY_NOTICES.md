# Third-party notices

## Grid-weight generator

`src/derive_grid_weights_modified_202312.py` is a modified copy of Juliane
Mai's [Grid Weights Generator for Raven](https://github.com/julemai/GridWeightsGenerator).
The original script is
[`derive_grid_weights.py`](https://github.com/julemai/GridWeightsGenerator/blob/main/derive_grid_weights.py).
The included copy keeps the original copyright notice:

```text
Copyright 2016-2020 Juliane Mai - juliane.mai(at)uwaterloo.ca
```

The bundled script is used under GNU Lesser General Public License, version 3
or any later version (`LGPL-3.0-or-later`). The local [LICENSE](LICENSE)
includes the LGPL v3 and incorporated GPL v3 terms.

### Role and local adaptation

The script calculates weights from the area of overlap between grid cells
and polygons. Our modification adds a second output through `-o2` /
`--output_file2`, containing `HRU_ID`, `Cell_index`, `Row`, `Col`, and `w_kl`.
In this table, `HRU_ID` stores the basin, subbasin, or HRU identifier selected
for aggregation.

`spatial_aggregation.py` reads this table and uses the weights to aggregate
the NetCDF time series. Juliane Mai developed the original grid-weight
algorithm. Our changes add the output needed by this workflow. Keep the
copyright and license notices when redistributing the script.

Suggested acknowledgement: Mai, J. *Grid Weights Generator for Raven*
[Computer software]. [GitHub](https://github.com/julemai/GridWeightsGenerator).

The other aggregation source files identify Hongren Shen, 2019-2026, under
`LGPL-3.0-or-later` in their SPDX headers.

## Data

Source data and derived example outputs are separate from the software license.
See [data sources and preparation](docs/inputs.md), [dataset citations](README.md#citation-and-data-sources),
and [reference outputs](README.md#outputs) for attribution and applicable terms.

The OWL-Hydro gauge inventory is CC BY 4.0. The routing and meteorological
products retain their provider attribution and source notices.
