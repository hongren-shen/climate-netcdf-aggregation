#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2019-2026 Hongren Shen <hongren.shen@uwaterloo.ca>
# SPDX-License-Identifier: LGPL-3.0-or-later

"""
workflow/aggregate_driver.py

"""

from __future__ import annotations

import os
import sys
import argparse
import copy
import re
import shutil
import time
import gc
from pathlib import Path
from typing import Dict, List, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

import pandas as pd
from tqdm import tqdm

dir_path   = os.path.dirname(os.path.realpath(__file__))
REPO_ROOT = Path(dir_path).parent
python_lib = "../src"
sys.path.insert(0, os.path.join(dir_path, python_lib))

from forcing_utils import (
    format_elapsed_time, 
    zip_and_cleanup_gauge_folder,
    discover_precomputed_region_agg,
    cleanup_intermediates_in_gauge_dir,
    export_metadata,
    is_daily_freq,
    is_temperature_var,
    temporal_aggregate,
    temporal_aggregate_specific,
    rename_with_label,
)

from spatial_aggregation import (
    aggregate_netcdf_forcing,
    normalize_hru_version,
    unit_ids_for_gauge,
    crop_gridded_nc_for_gauge,
    derive_gridweights_for_gridded_products,
    extract_gauge_from_regional_agg_nc,
    _write_version_gridweights_for_aggregated,
)

##############################################################################
debug_mode                          = False
# Maximum gauge workers for the current input file; regions and files run in sequence.
# Each task processes one gauge and one file. Use 1 to process gauges one at a time.
num_cores                           = 1

# Processing scope: inputs and outputs are organized by region.
region_for_process                  = ["SE"]

# Use [] to process all gauges
extract_gauge_list                  = ["02KA015"]

# Defaults use REPO_ROOT so they also work when launched from an IDE or another directory.
# Prepare forcing NetCDFs and convert source units before running, when needed.
nc_data_folder                      = [str(REPO_ROOT / 'data/casr_3.1')]
nc_files_process                    = [
                                        "CaSR_v3_1_A_PR0_SFC.nc",
                                        "CaSR_v3_1_A_RAIN_SFC.nc",
                                        "CaSR_v3_1_A_SNOW_SFC.nc",
                                        "CaSR_v3_1_A_TT_1_5m.nc",
                                    ]

# OWL-Hydro SE 24H gauge inventory, matching the bundled SE example.
# Inventory selection is independent of the forcing output frequency below.
gauge_info_path                     = [str(REPO_ROOT / "data/gauge_inventory/SE/24H/gauge_info.csv")]
gauge_routing_product_folder        = [str(REPO_ROOT / "data/Routing_product_simplified_extraction")]

# OLRRP base release used by the OWL-Hydro routing products
olrrp_version                       = "v2-0-1"
subbasin_name_prefix                = "finalcat_info_"
hru_name_prefix                     = "final_hru_"
outline_shp_name                    = "outline.shp"

# Outputs
output_folder                       = [str(REPO_ROOT / "results")]
grid_weights_file_prefix            = "gridweights_"

# Spatial aggregation controls
want_gridded_forcing                = False           # True: produce gauge-cropped gridded *.nc
want_aggregated_targets             = ["subbasin"]     # any of: "basin", "subbasin", "hru"

# OWL-Hydro convention: a = semi-distributed; b = routing-only.
# Full HRU version codes must match olrrp_version above.
# This list applies only to the "hru" target. Basin/subbasin outputs retain both HRU weight mappings.
hru_versions                        = ["v2-0-1a", "v2-0-1b"]
prefer_precomputed_regional         = False

# Temporal aggregation: one output frequency per run.
want_temporal_aggregation           = False    # False: preserve source timestamps and native time steps.
temporal_agg_freq                   = "D"      # Used only when enabled; examples: "3H", "6H", "12H", "D".
temporal_anchor_hour                = 0        # Local daily boundary hour; 0 = midnight.

# Daily temperature products: active only with temporal aggregation and a daily frequency.
derive_daily_temp_minmax            = False    # True: add daily MIN and MAX alongside the mean.
label_avg_when_enabled_only         = True     # Add AVG to daily temperature names when MIN/MAX are enabled.

# Time conventions: applied only when temporal aggregation is enabled.
temporal_tz_name                    = "UTC-5"  # Fixed UTC offset; no daylight saving time.
# Bundled CaSR inputs are period-ending; disable this shift for period-starting inputs.
period_ending2starting              = True     # With temporal aggregation, shift back by one source time step.

# IO / performance
aggregate_compress                  = True
aggregate_compress_level            = 1
aggregate_overwrite                 = True
time_block_hours                    = 24 * 365
output_time_chunk_hours             = 24 * 365

# Output packaging
if_zip_gauge_folder                 = True
if_delete_original_gauge_folder     = True

# Grid-weight helper and crop padding
grid_weights_generator              = str(REPO_ROOT / "src/derive_grid_weights_modified_202312.py")
crop_buffer_cell                    = 5
##############################################################################

# Top-level settings remain the defaults; command-line options override them for one run.
_CONFIG_DEFAULTS = {
    name: copy.deepcopy(globals()[name])
    for name in (
        "debug_mode", "num_cores", "region_for_process", "extract_gauge_list",
        "nc_data_folder", "nc_files_process", "gauge_info_path", "gauge_routing_product_folder",
        "olrrp_version", "subbasin_name_prefix", "hru_name_prefix", "outline_shp_name",
        "output_folder", "grid_weights_file_prefix", "want_gridded_forcing",
        "want_aggregated_targets", "hru_versions", "prefer_precomputed_regional",
        "want_temporal_aggregation", "temporal_agg_freq", "temporal_anchor_hour",
        "derive_daily_temp_minmax", "label_avg_when_enabled_only", "temporal_tz_name",
        "period_ending2starting", "aggregate_compress", "aggregate_compress_level",
        "aggregate_overwrite", "time_block_hours", "output_time_chunk_hours",
        "if_zip_gauge_folder", "if_delete_original_gauge_folder",
        "grid_weights_generator", "crop_buffer_cell",
    )
}


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Show defaults and preserve example line breaks in command-line help."""


def build_argument_parser() -> argparse.ArgumentParser:
    """Describe public workflow options while retaining editable top-level defaults."""
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate prepared gridded forcing NetCDFs to gauge polygons, "
            "with optional temporal aggregation. Default paths are anchored to the "
            "repository; relative command-line paths use the current working directory."
        ),
        epilog=(
            "Examples (from workflow/):\n"
            "  python aggregate_driver.py --help\n"
            "  python aggregate_driver.py --regions SE --gauges 02KA015 --workers 1\n"
            "  python aggregate_driver.py --targets hru --hru-versions v2-0-1a v2-0-1b\n"
            "  python aggregate_driver.py --temporal --frequency D --daily-temp-minmax\n\n"
            "Use a dedicated output directory: each selected output/<region> directory "
            "is cleared before processing, including when --no-overwrite is set."
        ),
        formatter_class=_HelpFormatter,
    )
    parser.set_defaults(**copy.deepcopy(_CONFIG_DEFAULTS))
    inputs = parser.add_argument_group("Inputs and gauge selection")
    inputs.add_argument("--regions", dest="region_for_process", nargs="+", metavar="REGION",
                        help="Region subfolders to process.")
    inputs.add_argument("--gauges", dest="extract_gauge_list", nargs="*", metavar="GAUGE",
                        help="Gauge IDs to select; supply --gauges with no IDs to select all gauges in the regions.")
    inputs.add_argument("--input-root", dest="nc_data_folder", nargs=1, metavar="PATH",
                        help="Prepared forcing root; inputs are PATH/REGION/VARIABLE.nc.")
    inputs.add_argument("--variables", dest="nc_files_process", nargs="+", metavar="FILE",
                        help="Prepared NetCDF filenames. Convert source units before running when needed.")
    inputs.add_argument("--gauge-table", dest="gauge_info_path", nargs=1, metavar="CSV",
                        help="Gauge CSV containing Region and Obs_NM columns. Defaults to the bundled OWL-Hydro SE 24H inventory; this does not set the forcing output frequency.")
    inputs.add_argument("--routing-root", dest="gauge_routing_product_folder", nargs=1, metavar="PATH",
                        help="Gauge routing root; shapefiles are under PATH/REGION/GAUGE/.")

    spatial = parser.add_argument_group("Spatial aggregation and OWL-Hydro versions")
    spatial.add_argument("--targets", dest="want_aggregated_targets", nargs="+",
                         choices=("basin", "subbasin", "hru"), help="Polygon levels to aggregate.")
    spatial.add_argument("--olrrp-version", dest="olrrp_version", metavar="VERSION",
                         help="Base routing-product release used in shapefile names.")
    spatial.add_argument("--hru-versions", dest="hru_versions", nargs="+", metavar="VERSION",
                         help="HRU targets only. OWL-Hydro convention: v2-0-1a = semi-distributed; v2-0-1b = routing-only. Must match --olrrp-version.")
    spatial.add_argument("--gridded", dest="want_gridded_forcing", action=argparse.BooleanOptionalAction,
                         help="Also export gauge-cropped gridded forcing and its grid weights.")

    temporal = parser.add_argument_group("Temporal aggregation")
    temporal.add_argument("--temporal", dest="want_temporal_aggregation", action=argparse.BooleanOptionalAction,
                          help="Enable time aggregation; when disabled, preserve source timestamps and native time steps.")
    temporal.add_argument("--frequency", dest="temporal_agg_freq", metavar="FREQ",
                          help="One output frequency per run: D or a positive whole-hour interval such as 3H, 6H, or 12H. Does not interpolate finer time steps.")
    temporal.add_argument("--timezone", dest="temporal_tz_name", metavar="ZONE",
                          help="Target clock for temporal aggregation: UTC, a fixed offset such as UTC-5, or an IANA location such as America/Toronto. Fixed offsets exclude daylight saving time.")
    temporal.add_argument("--anchor-hour", dest="temporal_anchor_hour", type=int, metavar="HOUR",
                          help="Local boundary hour for daily aggregation, from 0 to 23.")
    temporal.add_argument("--period-ending", dest="period_ending2starting", action=argparse.BooleanOptionalAction,
                          help="With --temporal, treat input labels as period ends and shift them back by one source time step. Disable for period-starting inputs.")
    temporal.add_argument("--daily-temp-minmax", dest="derive_daily_temp_minmax", action=argparse.BooleanOptionalAction,
                          help="Add temperature MIN and MAX from input time steps when temporal aggregation and a daily frequency are enabled.")
    temporal.add_argument("--label-avg-when-enabled-only", dest="label_avg_when_enabled_only", action=argparse.BooleanOptionalAction,
                          help="Label daily mean temperature AVG only when MIN/MAX are enabled; the negative option always adds AVG for daily temperature.")

    outputs = parser.add_argument_group("Outputs and execution")
    outputs.add_argument("--output-root", dest="output_folder", nargs=1, metavar="PATH",
                         help="Dedicated output root. Existing PATH/REGION directories are cleared before processing.")
    outputs.add_argument("--workers", dest="num_cores", type=int, metavar="N",
                         help="Maximum gauge worker processes for the current input file. Regions and files run in sequence; use 1 to process gauges one at a time.")
    outputs.add_argument("--debug", dest="debug_mode", action=argparse.BooleanOptionalAction,
                         help="Process only the first two selected gauges in each region; worker count is unchanged.")
    outputs.add_argument("--compress", dest="aggregate_compress", action=argparse.BooleanOptionalAction,
                         help="Compress aggregated NetCDF outputs.")
    outputs.add_argument("--compression-level", dest="aggregate_compress_level", type=int, metavar="LEVEL",
                         help="NetCDF compression level from 0 to 9.")
    outputs.add_argument("--overwrite", dest="aggregate_overwrite", action=argparse.BooleanOptionalAction,
                         help="Allow overwriting individual aggregation files; this does not disable initial regional output cleanup.")
    outputs.add_argument("--time-block-records", dest="time_block_hours", type=int, metavar="N",
                         help="Time records per spatial processing block; the legacy setting name is time_block_hours.")
    outputs.add_argument("--output-time-chunk-records", dest="output_time_chunk_hours", type=int, metavar="N",
                         help="Maximum time records per output chunk; the legacy setting name is output_time_chunk_hours.")
    outputs.add_argument("--zip", dest="if_zip_gauge_folder", action=argparse.BooleanOptionalAction,
                         help="Create one ZIP archive per gauge after processing.")
    outputs.add_argument("--delete-after-zip", dest="if_delete_original_gauge_folder", action=argparse.BooleanOptionalAction,
                         help="Remove the uncompressed gauge folder after zipping; inactive with --no-zip.")

    advanced = parser.add_argument_group("Advanced shapefile and grid settings")
    advanced.add_argument("--subbasin-prefix", dest="subbasin_name_prefix", metavar="PREFIX",
                          help="Prefix before the base version in subbasin shapefile names.")
    advanced.add_argument("--hru-prefix", dest="hru_name_prefix", metavar="PREFIX",
                          help="Prefix before the full version in HRU shapefile names.")
    advanced.add_argument("--outline-name", dest="outline_shp_name", metavar="FILE",
                          help="Gauge outline shapefile name.")
    advanced.add_argument("--gridweights-prefix", dest="grid_weights_file_prefix", metavar="PREFIX",
                          help="Output prefix for gridded-product weight files.")
    advanced.add_argument("--gridweights-generator", dest="grid_weights_generator", metavar="SCRIPT",
                          help="Path to the grid-weight helper script.")
    advanced.add_argument("--crop-buffer-cells", dest="crop_buffer_cell", type=int, metavar="N",
                          help="Nonnegative grid-cell padding around the gauge crop.")
    return parser


def parse_runtime_config(argv: Optional[List[str]] = None) -> dict:
    """Parse and validate settings before any processing or output cleanup."""
    parser = build_argument_parser()
    config = vars(parser.parse_args(argv))
    for key, option in (
        ("num_cores", "--workers"), ("time_block_hours", "--time-block-records"),
        ("output_time_chunk_hours", "--output-time-chunk-records"),
    ):
        if config[key] < 1:
            parser.error(f"{option} must be a positive integer.")
    if not 0 <= config["temporal_anchor_hour"] <= 23:
        parser.error("--anchor-hour must be between 0 and 23.")
    if not 0 <= config["aggregate_compress_level"] <= 9:
        parser.error("--compression-level must be between 0 and 9.")
    if config["crop_buffer_cell"] < 0:
        parser.error("--crop-buffer-cells must be nonnegative.")

    frequency = config["temporal_agg_freq"].strip().upper()
    if frequency not in {"D", "1D"} and not re.fullmatch(r"(?:[1-9]\d*)?H", frequency):
        parser.error("--frequency must be D or a positive whole-hour interval such as 3H, 6H, or 24H.")
    config["temporal_agg_freq"] = frequency
    timezone = config["temporal_tz_name"]
    try:
        if "/" in timezone:
            ZoneInfo(timezone)
        elif timezone.upper() != "UTC":
            if not re.fullmatch(r"UTC[+-]\d+(?:\.\d+)?", timezone.upper()):
                raise ValueError("invalid fixed UTC offset")
            if not -24 < float(timezone[3:]) < 24:
                raise ValueError("UTC offset outside range")
    except (ValueError, ZoneInfoNotFoundError):
        parser.error("--timezone must be UTC, a fixed offset such as UTC-5, or a valid IANA location such as America/Toronto.")

    if "hru" in config["want_aggregated_targets"]:
        try:
            for version in config["hru_versions"]:
                normalize_hru_version(version, olrrp_version=config["olrrp_version"])
        except ValueError as exc:
            parser.error(str(exc))

    for region in config["region_for_process"]:
        if not region.strip() or region in {".", ".."} or Path(region).name != region or "/" in region or "\\" in region:
            parser.error("--regions values must be individual subfolder names, such as SE.")
        output = (Path(config["output_folder"][0]) / region).resolve()
        inputs = [
            (Path(config["nc_data_folder"][0]) / region).resolve(),
            (Path(config["gauge_routing_product_folder"][0]) / region).resolve(),
            Path(config["gauge_info_path"][0]).resolve(),
            Path(config["grid_weights_generator"]).resolve(),
        ]
        if any(path == output or output in path.parents or path in output.parents for path in inputs):
            parser.error("--output-root must keep regional outputs separate from forcing inputs, routing products, the gauge table, and the grid-weight helper.")
    return config


def load_gauge_table_and_check_forcing(config: dict) -> pd.DataFrame:
    """Check the inventory, forcing, routing, and helper before output cleanup."""
    table_path = Path(config["gauge_info_path"][0]).resolve()
    if not table_path.is_file():
        raise FileNotFoundError(
            f"Gauge table not found: {table_path}\n"
            "Set gauge_info_path or --gauge-table to an OWL-Hydro regional gauge_info.csv "
            "containing Region and Obs_NM columns."
        )
    gauge_df = pd.read_csv(table_path, dtype={"Region": "string", "Obs_NM": "string"})
    missing_columns = {"Region", "Obs_NM"} - set(gauge_df.columns)
    if missing_columns:
        raise ValueError(f"Gauge table {table_path} is missing columns: {', '.join(sorted(missing_columns))}.")
    gauge_df["Obs_NM"] = gauge_df["Obs_NM"].str.replace(" ", "_", regex=False)

    selected = gauge_df.loc[gauge_df["Region"].isin(config["region_for_process"])]
    if config["extract_gauge_list"]:
        selected = selected.loc[selected["Obs_NM"].isin(config["extract_gauge_list"])]
    if selected.empty:
        raise ValueError("No gauges match --regions and --gauges in the selected gauge table.")
    if config["debug_mode"]:
        selected = selected.groupby("Region", sort=False).head(2)
    missing_files = []
    for region in selected["Region"].drop_duplicates():
        for filename in config["nc_files_process"]:
            path = (Path(config["nc_data_folder"][0]) / region / filename).resolve()
            if not path.is_file():
                missing_files.append(str(path))
    if missing_files:
        raise FileNotFoundError(
            "Prepared forcing NetCDF files not found:\n  " + "\n  ".join(missing_files) +
            "\nSet nc_data_folder or --input-root to a root containing REGION/VARIABLE.nc. "
            "Use --variables to select the files to process."
        )

    helper_path = Path(config["grid_weights_generator"]).resolve()
    if not helper_path.is_file():
        raise FileNotFoundError(
            f"Grid-weight helper not found: {helper_path}\n"
            "Set grid_weights_generator or --gridweights-generator to an existing script."
        )

    targets = set(config["want_aggregated_targets"])
    base_version = config["olrrp_version"]
    missing_routing = set()
    for row in selected.itertuples(index=False):
        folder = Path(config["gauge_routing_product_folder"][0]) / row.Region / row.Obs_NM
        names = set()
        suffixes = set()
        if "basin" in targets:
            names.update((config["outline_shp_name"], "outline.shp"))
        if "subbasin" in targets:
            names.update((f"{config['subbasin_name_prefix']}{base_version}.shp",
                          f"finalcat_info_{base_version}.shp"))
        if targets & {"basin", "subbasin"}:
            suffixes.update(("a", "b"))
        if "hru" in targets:
            suffixes.update(normalize_hru_version(version, olrrp_version=base_version)[0]
                            for version in config["hru_versions"])
        if config["want_gridded_forcing"]:
            names.update((config["outline_shp_name"],
                          f"{config['subbasin_name_prefix']}{base_version}.shp"))
            available = {suffix for suffix in ("a", "b")
                         if (folder / f"{config['hru_name_prefix']}{base_version}{suffix}.shp").is_file()}
            suffixes.update(available or {"a", "b"})
        names.update(f"{config['hru_name_prefix']}{base_version}{suffix}.shp" for suffix in suffixes)
        for name in names:
            path = folder / name
            paths = [path]
            if path.suffix.lower() == ".shp":
                paths.extend(path.with_suffix(suffix) for suffix in (".shx", ".dbf", ".prj"))
            missing_routing.update(str(item.resolve()) for item in paths if not item.is_file())
    if missing_routing:
        paths = sorted(missing_routing)
        detail = "\n  ".join(paths[:20])
        if len(paths) > 20:
            detail += f"\n  ... and {len(paths) - 20} more missing files."
        raise FileNotFoundError(
            "Required routing files not found:\n  " + detail +
            "\nSet --routing-root to matching region/gauge folders and retain shapefile companions."
        )
    return gauge_df


def _configure_worker(config: dict) -> None:
    """Apply the same validated settings in the parent and spawned gauge workers."""
    for name in _CONFIG_DEFAULTS:
        globals()[name] = copy.deepcopy(config[name])


def _require_output(path: Optional[Path], description: str) -> Path:
    """Reject a skipped or missing requested output instead of reporting success."""
    if path is None or not Path(path).is_file():
        raise RuntimeError(f"Requested output was not produced: {description}")
    return Path(path)


def _process_one_gauge(
    gid: str,
    regional_gridded_nc: Path,
    shp_dir: Path,
    out_root: Path,
    var_file_key: str,
    precomputed: Dict[str, List[Path]],
) -> str:
    """
    Process one gauge for one input file, including its selected spatial and temporal outputs.
    """
    gauge_dir = out_root / gid
    gauge_dir.mkdir(parents=True, exist_ok=True)

    # (1) Crop gridded NETCDF (optional) + weights for GRIDDed product
    cropped_nc: Optional[Path] = None
    if want_gridded_forcing:
        cropped_nc = crop_gridded_nc_for_gauge(
            copy_region_forcing_nc_path=regional_gridded_nc,
            gauge_outline_shp=shp_dir / outline_shp_name,
            gauge_output_folder=gauge_dir,
            gauge_ID=gid,
            grid_weights_generator=Path(grid_weights_generator),
            crop_buffer_cell=crop_buffer_cell,
            use_nco=False,
            compress=True,
            compression_level=aggregate_compress_level,
            sample=False,
            sample_time_start=None,
            sample_time_end=None,
            sample_suffix="sample",
            common_header=None,
        )
        cropped_nc = _require_output(cropped_nc, f"gridded forcing for {gid}, {var_file_key}")

        weights_sentinel = gauge_dir / ".gridded_weights_done"
        if not weights_sentinel.exists() and cropped_nc is not None:
            derive_gridweights_for_gridded_products(
                nc_path=cropped_nc,
                shp_dir=shp_dir,
                out_dir=gauge_dir,
                grid_weights_generator=Path(grid_weights_generator),
                olrrp_version=olrrp_version,
                hru_name_prefix=hru_name_prefix,
                subbasin_name_prefix=subbasin_name_prefix,
                outline_shp_name=outline_shp_name,
                grid_weights_file_prefix=grid_weights_file_prefix,
                common_header=None,
                label=gid,
            )
            try:
                weights_sentinel.write_text("done", encoding="utf-8")
            except Exception:
                pass

    src_for_agg = cropped_nc if cropped_nc is not None else regional_gridded_nc

    # temperature daily context
    daily_temp_context = is_daily_freq(temporal_agg_freq) and is_temperature_var(var_file_key)

    # (2) Aggregate targets
    for tgt in want_aggregated_targets:
        tgt = tgt.strip().lower()
        if tgt not in {"basin", "subbasin", "hru"}:
            continue

        ver_idx = 1 if tgt == "basin" else 2 if tgt == "subbasin" else 3

        if tgt == "hru":
            for hv in hru_versions:
                suffix_req, _ = normalize_hru_version(hv, olrrp_version=olrrp_version)
                out_nc: Optional[Path] = None
                used_precomputed = False
                suffix_used = suffix_req

                if prefer_precomputed_regional:
                    pool = (precomputed.get(f"hru_{suffix_req}") or []) + (precomputed.get("hru") or [])
                    if pool:
                        src = max(pool, key=lambda p: p.stat().st_mtime)
                        suffix_used = suffix_req
                        try:
                            ids = unit_ids_for_gauge(
                                target="hru",
                                shp_dir=shp_dir,
                                outline_shp_name=outline_shp_name,
                                subbasin_name_prefix=subbasin_name_prefix,
                                hru_name_prefix=hru_name_prefix,
                                olrrp_version=olrrp_version,
                                hru_suffix=suffix_used,
                            )
                            out_name = f"{src.stem}_{gid}.nc"
                            out_nc = gauge_dir / out_name
                            extract_gauge_from_regional_agg_nc(
                                agg_nc_path=src,
                                unit_ids_to_keep=ids,
                                out_nc_path=out_nc,
                                compress=aggregate_compress,
                                compression_level=aggregate_compress_level,
                            )
                            used_precomputed = True
                        except Exception:
                            used_precomputed = False

                if not used_precomputed:
                    out_nc = aggregate_netcdf_forcing(
                        name_key=gid,
                        output_root=gauge_dir,
                        shapefile_root=shp_dir,
                        grid_weights_generator=Path(grid_weights_generator),
                        olrrp_version=olrrp_version,
                        nc_path=src_for_agg,
                        hru_name_prefix=hru_name_prefix,
                        subbasin_name_prefix=subbasin_name_prefix,
                        outline_shp_name=outline_shp_name,
                        hru_version=hv,
                        aggregate_target="hru",
                        compress=aggregate_compress,
                        compression_level=aggregate_compress_level,
                        overwrite=aggregate_overwrite,
                        verbose=False,
                        time_block_hours=time_block_hours,
                        output_time_chunk_hours=output_time_chunk_hours,
                        write_gridweights=False,
                    )

                out_nc = _require_output(out_nc, f"HRU {hv} forcing for {gid}, {var_file_key}")

                # PRE-temporal weights
                gw = gauge_dir / f"version{ver_idx}_gridweights_agg_hru_{olrrp_version}{suffix_used}.txt"
                if not gw.exists():
                    _write_version_gridweights_for_aggregated(
                        out_txt=gw,
                        agg_nc_path=out_nc,
                        version=ver_idx,
                        name_key=gid,
                        olrrp_version=olrrp_version,
                        hru_name_prefix=hru_name_prefix,
                        suffix_letter=suffix_used,
                        shapefile_root=shp_dir,
                        common_header=None,
                    )
                _require_output(gw, f"HRU grid weights for {gid}, {hv}")

                # Temporal aggregation
                if want_temporal_aggregation:
                    out_mean = temporal_aggregate(
                        out_nc,
                        region_name=gid,
                        var_file_key=var_file_key,
                        temporal_freq=temporal_agg_freq,
                        temporal_anchor_hour=temporal_anchor_hour,
                        temporal_tz_name=temporal_tz_name,
                        period_ending2starting=period_ending2starting,
                        compress=aggregate_compress,
                        compression_level=aggregate_compress_level,
                        output_time_chunk_hours=output_time_chunk_hours,
                        overwrite=aggregate_overwrite,
                    )
                    if daily_temp_context and (not label_avg_when_enabled_only or derive_daily_temp_minmax):
                        out_mean = rename_with_label(out_mean, gid, "AVG", overwrite=aggregate_overwrite)

                    out_mean = _require_output(out_mean, f"temporal forcing for {gid}, {var_file_key}")
                    export_metadata(out_mean, gid, gauge_dir, stat_hint=("AVG" if daily_temp_context else None))

                    # MIN/MAX
                    if derive_daily_temp_minmax and daily_temp_context:
                        pmax = temporal_aggregate_specific(
                            out_nc,
                            region_name=gid,
                            temporal_freq=temporal_agg_freq,
                            temporal_anchor_hour=temporal_anchor_hour,
                            temporal_tz_name=temporal_tz_name,
                            period_ending2starting=period_ending2starting,
                            compress=aggregate_compress,
                            compression_level=aggregate_compress_level,
                            output_time_chunk_hours=output_time_chunk_hours,
                            overwrite=aggregate_overwrite,
                            reduce_method="max",
                            label="MAX",
                        )
                        pmax = _require_output(pmax, f"daily temperature MAX for {gid}, {hv}")
                        export_metadata(pmax, gid, gauge_dir, stat_hint="MAX")

                        pmin = temporal_aggregate_specific(
                            out_nc,
                            region_name=gid,
                            temporal_freq=temporal_agg_freq,
                            temporal_anchor_hour=temporal_anchor_hour,
                            temporal_tz_name=temporal_tz_name,
                            period_ending2starting=period_ending2starting,
                            compress=aggregate_compress,
                            compression_level=aggregate_compress_level,
                            output_time_chunk_hours=output_time_chunk_hours,
                            overwrite=aggregate_overwrite,
                            reduce_method="min",
                            label="MIN",
                        )
                        pmin = _require_output(pmin, f"daily temperature MIN for {gid}, {hv}")
                        export_metadata(pmin, gid, gauge_dir, stat_hint="MIN")

        else:
            # basin/subbasin
            out_nc = aggregate_netcdf_forcing(
                name_key=gid,
                output_root=gauge_dir,
                shapefile_root=shp_dir,
                grid_weights_generator=Path(grid_weights_generator),
                olrrp_version=olrrp_version,
                nc_path=src_for_agg,
                hru_name_prefix=hru_name_prefix,
                subbasin_name_prefix=subbasin_name_prefix,
                outline_shp_name=outline_shp_name,
                hru_version="semi-distributed",
                aggregate_target=tgt,
                compress=aggregate_compress,
                compression_level=aggregate_compress_level,
                overwrite=aggregate_overwrite,
                verbose=False,
                time_block_hours=time_block_hours,
                output_time_chunk_hours=output_time_chunk_hours,
                write_gridweights=False,
            )

            out_nc = _require_output(out_nc, f"{tgt} forcing for {gid}, {var_file_key}")

            # PRE-temporal weights
            for suf in ("a", "b"):
                gw = gauge_dir / f"version{ver_idx}_gridweights_agg_{tgt}_{olrrp_version}{suf}.txt"
                if not gw.exists():
                    _write_version_gridweights_for_aggregated(
                        out_txt=gw,
                        agg_nc_path=out_nc,
                        version=ver_idx,
                        name_key=gid,
                        olrrp_version=olrrp_version,
                        hru_name_prefix=hru_name_prefix,
                        suffix_letter=suf,
                        shapefile_root=shp_dir,
                        common_header=None,
                    )
                _require_output(gw, f"{tgt} grid weights for {gid}, layout {suf}")

            if want_temporal_aggregation:
                out_mean = temporal_aggregate(
                    out_nc,
                    region_name=gid,
                    var_file_key=var_file_key,
                    temporal_freq=temporal_agg_freq,
                    temporal_anchor_hour=temporal_anchor_hour,
                    temporal_tz_name=temporal_tz_name,
                    period_ending2starting=period_ending2starting,
                    compress=aggregate_compress,
                    compression_level=aggregate_compress_level,
                    output_time_chunk_hours=output_time_chunk_hours,
                    overwrite=aggregate_overwrite,
                )
                if daily_temp_context and (not label_avg_when_enabled_only or derive_daily_temp_minmax):
                    out_mean = rename_with_label(out_mean, gid, "AVG", overwrite=aggregate_overwrite)
                out_mean = _require_output(out_mean, f"temporal forcing for {gid}, {var_file_key}")
                export_metadata(out_mean, gid, gauge_dir, stat_hint=("AVG" if daily_temp_context else None))

                if derive_daily_temp_minmax and daily_temp_context:
                    pmax = temporal_aggregate_specific(
                        out_nc, region_name=gid,
                        temporal_freq=temporal_agg_freq, temporal_anchor_hour=temporal_anchor_hour,
                        temporal_tz_name=temporal_tz_name, period_ending2starting=period_ending2starting,
                        compress=aggregate_compress, compression_level=aggregate_compress_level,
                        output_time_chunk_hours=output_time_chunk_hours, overwrite=aggregate_overwrite,
                        reduce_method="max", label="MAX",
                    )
                    pmax = _require_output(pmax, f"daily temperature MAX for {gid}, {tgt}")
                    export_metadata(pmax, gid, gauge_dir, stat_hint="MAX")

                    pmin = temporal_aggregate_specific(
                        out_nc, region_name=gid,
                        temporal_freq=temporal_agg_freq, temporal_anchor_hour=temporal_anchor_hour,
                        temporal_tz_name=temporal_tz_name, period_ending2starting=period_ending2starting,
                        compress=aggregate_compress, compression_level=aggregate_compress_level,
                        output_time_chunk_hours=output_time_chunk_hours, overwrite=aggregate_overwrite,
                        reduce_method="min", label="MIN",
                    )
                    pmin = _require_output(pmin, f"daily temperature MIN for {gid}, {tgt}")
                    export_metadata(pmin, gid, gauge_dir, stat_hint="MIN")

    return gid

##############################################################################
if __name__ == "__main__":

    runtime_config = parse_runtime_config()
    try:
        gauge_df = load_gauge_table_and_check_forcing(runtime_config)
    except (OSError, ValueError) as exc:
        build_argument_parser().error(str(exc))
    _configure_worker(runtime_config)

    print("\n====================== Running with configuration ======================")
    print(f"                     Debug mode:  {debug_mode}")
    print(f"         Gauge worker processes:  {num_cores}")
    print(f"           Want gridded forcing:  {want_gridded_forcing}")
    print(f"  Want aggregated targets (set):  {want_aggregated_targets}")
    if "hru" in want_aggregated_targets:
        print(f"    HRU versions (OWL-Hydro):  {hru_versions}")
    print(f"  Prefer precomputed region-agg:  {prefer_precomputed_regional}")
    print(f"   Temporal aggregation enabled:  {want_temporal_aggregation}")
    if want_temporal_aggregation:
        print(f"        Temporal rule / tz / ps:  {temporal_agg_freq} / {temporal_tz_name or 'UTC'} / {period_ending2starting}")
        print(f"     Daily T min/max enabled?:    {derive_daily_temp_minmax}")
    print("========================================================================\n")

    for region_i in region_for_process:
        start_time_region = time.time()
        print(f">>> Processing gauge-level forcing for region {region_i}...")

        gauge_info = gauge_df.loc[gauge_df["Region"] == region_i].reset_index(drop=True)
        if extract_gauge_list:
            gauge_info = gauge_info.loc[gauge_info["Obs_NM"].isin(extract_gauge_list)].reset_index(drop=True)
        if debug_mode:
            gauge_info = gauge_info.head(2)

        gauge_ID_all = gauge_info["Obs_NM"].tolist()
        print(f"    Total gauges: {len(gauge_ID_all)}\n")
        if not gauge_ID_all:
            continue

        input_gauge_routing_product_path = Path(gauge_routing_product_folder[0]) / region_i
        output_gauge_forcing_folder = Path(output_folder[0]) / region_i

        shutil.rmtree(output_gauge_forcing_folder, ignore_errors=True)
        output_gauge_forcing_folder.mkdir(parents=True, exist_ok=True)

        # Variables on the same grid reuse each gauge's weights and grid-info files.
        # Process files in sequence so shared weights are complete before reuse.
        for netcdf_file_i in nc_files_process:
            print(f">>> Source NETCDF (regional gridded): {netcdf_file_i}")
            start_time_var = time.time()

            regional_gridded_nc = Path(nc_data_folder[0]) / region_i / netcdf_file_i
            
            print(regional_gridded_nc)

            precomputed = {}
            if prefer_precomputed_regional:
                
                precomputed = discover_precomputed_region_agg(
                    output_region_forcing_folder=Path(output_folder[0]) / region_i,
                    var_file_key=netcdf_file_i,
                    olrrp_version=olrrp_version)

            futures: List = []
            with ProcessPoolExecutor(
                max_workers=num_cores, initializer=_configure_worker, initargs=(runtime_config,),
            ) as executor:
                for gid in gauge_ID_all:
                    shp_dir = input_gauge_routing_product_path / gid
                    futures.append(
                        executor.submit(
                            _process_one_gauge,
                            gid,
                            regional_gridded_nc,
                            shp_dir,
                            output_gauge_forcing_folder,
                            netcdf_file_i,
                            precomputed,
                        )
                    )

                failures = []
                for fu in tqdm(as_completed(futures), total=len(futures), desc="\tGauges processing", unit="gauge"):
                    try:
                        fu.result()
                    except Exception as e:
                        print(f"[ERROR] gauge task failed: {e}")
                        failures.append(str(e))
            if failures:
                raise SystemExit(
                    f"Gauge processing failed for {netcdf_file_i}; outputs were left uncompressed for inspection."
                )

            print(f"    Time used for {netcdf_file_i}: {format_elapsed_time(time.time() - start_time_var)}\n")

        # Cleanup
        print("    Cleaning intermediate files (per-gauge)...")
        for gid in gauge_ID_all:
            removed = cleanup_intermediates_in_gauge_dir(output_gauge_forcing_folder / gid)
            if removed and debug_mode:
                print(f"    - {gid}: removed {removed} file(s)")

        # Zip
        if if_zip_gauge_folder:
            print("    Zipping gauge forcing folders...")
            gauge_folders = [output_gauge_forcing_folder / gid for gid in gauge_ID_all]
            with ProcessPoolExecutor(max_workers=num_cores) as executor:
                futures = [
                    executor.submit(zip_and_cleanup_gauge_folder, folder, True, if_delete_original_gauge_folder)
                    for folder in gauge_folders
                ]
                failures = []
                for future in tqdm(as_completed(futures), total=len(futures), desc="    Zipping folders", unit="folder"):
                    try:
                        future.result()
                    except Exception as e:
                        print(f"[ERROR] ZIP task failed: {e}")
                        failures.append(str(e))
            if failures:
                raise SystemExit("ZIP packaging failed; inspect the retained gauge folders and error messages.")

        print(f"\n>>> Total time (region {region_i}): {format_elapsed_time(time.time() - start_time_region)}\n")
        gc.collect()
