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
import shutil
import time
import gc
from pathlib import Path
from typing import Dict, List, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm

dir_path   = os.path.dirname(os.path.realpath(__file__))
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
    unit_ids_for_gauge,
    crop_gridded_nc_for_gauge,
    derive_gridweights_for_gridded_products,
    extract_gauge_from_regional_agg_nc,
    _write_version_gridweights_for_aggregated,
)

##############################################################################
debug_mode                          = False
num_cores                           = 20

# config
# Assume we process region by region
# Input data are also organized region by region
region_for_process                  = ["SE"]

# Use [] to process all gauges
extract_gauge_list                  = ["02KA015", "02KF012", "02KF015"]

# Inputs
nc_data_folder                      = ['../data/casr_3.1/']
nc_files_process                    = [
                                        "CaSR_v3_1_A_PR0_SFC.nc",
                                        "CaSR_v3_1_A_RAIN_SFC.nc",
                                        "CaSR_v3_1_A_SNOW_SFC.nc",
                                        "CaSR_v3_1_A_TT_1_5m.nc",
                                    ]

# Gauge info
gauge_info_path                     = ["../data/01_Ontario_gauge_info.csv"]
gauge_routing_product_folder        = ["../data/Routing_product_simplified_extraction/"]  # gauge routing product

# OLRRP
olrrp_version                       = "v2-0-1"
subbasin_name_prefix                = "finalcat_info_"
hru_name_prefix                     = "final_hru_"
outline_shp_name                    = "outline.shp"

# Outputs
output_folder                       = ["../results"]
grid_weights_file_prefix            = "gridweights_"

# Spatial aggregation controls
want_gridded_forcing                = False           # True: produce gauge-cropped gridded *.nc
want_aggregated_targets             = ["subbasin"]     # any of: "basin", "subbasin", "hru"
hru_versions                        = ["semi-distributed", "routing-only"]  # only used when target includes "hru"
prefer_precomputed_regional         = False

# Temporal aggregation controls
want_temporal_aggregation           = False    # False: No temporal aggregation
temporal_agg_freq                   = "D"      # "D","24H","3H","M",...
temporal_anchor_hour                = 0

# Daily temperature min/max
derive_daily_temp_minmax            = False    # True: derive daily min & max temperature
label_avg_when_enabled_only         = True     # should use True when derive_daily_temp_minmax = True

# Time zone & time conventions
temporal_tz_name                    = "UTC-5"
period_ending2starting              = True

# IO / performance
aggregate_compress                  = True
aggregate_compress_level            = 1
aggregate_overwrite                 = True
time_block_hours                    = 24 * 365
output_time_chunk_hours             = 24 * 365

if_zip_gauge_folder                 = True
if_delete_original_gauge_folder     = True

# should never change
grid_weights_generator              = "../src/derive_grid_weights_modified_202312.py"
crop_buffer_cell                    = 3
##############################################################################

def _process_one_gauge(
    gid: str,
    regional_gridded_nc: Path,
    shp_dir: Path,
    out_root: Path,
    var_file_key: str,
    precomputed: Dict[str, List[Path]],
) -> str:
    """
    Per-gauge job for one variable file. For parallel processing use.
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
                suffix_req = "a" if hv.lower().startswith("semi") else "b"
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

                if out_nc is None:
                    continue

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
                        if pmax:
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
                        if pmin:
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

            if out_nc is None:
                continue

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
                    if pmax:
                        export_metadata(pmax, gid, gauge_dir, stat_hint="MAX")

                    pmin = temporal_aggregate_specific(
                        out_nc, region_name=gid,
                        temporal_freq=temporal_agg_freq, temporal_anchor_hour=temporal_anchor_hour,
                        temporal_tz_name=temporal_tz_name, period_ending2starting=period_ending2starting,
                        compress=aggregate_compress, compression_level=aggregate_compress_level,
                        output_time_chunk_hours=output_time_chunk_hours, overwrite=aggregate_overwrite,
                        reduce_method="min", label="MIN",
                    )
                    if pmin:
                        export_metadata(pmin, gid, gauge_dir, stat_hint="MIN")

    return gid

##############################################################################
if __name__ == "__main__":

    print("\n====================== Running with configuration ======================")
    print(f"                     Debug mode:  {debug_mode}")
    print(f"                Number of cores:  {num_cores}")
    print(f"           Want gridded forcing:  {want_gridded_forcing}")
    print(f"  Want aggregated targets (set):  {want_aggregated_targets}")
    print(f"  Prefer precomputed region-agg:  {prefer_precomputed_regional}")
    print(f"   Temporal aggregation enabled:  {want_temporal_aggregation}")
    if want_temporal_aggregation:
        print(f"        Temporal rule / tz / ps:  {temporal_agg_freq} / {temporal_tz_name or 'UTC'} / {period_ending2starting}")
        print(f"     Daily T min/max enabled?:    {derive_daily_temp_minmax}")
    print("========================================================================\n")

    # Gauge list
    gauge_df = pd.read_csv(gauge_info_path[0])
    gauge_df["Obs_NM"] = gauge_df["Obs_NM"].apply(lambda x: x.replace(" ", "_") if isinstance(x, str) else x)

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
            with ProcessPoolExecutor(max_workers=num_cores) as executor:
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

                for fu in tqdm(as_completed(futures), total=len(futures), desc="\tGauges processing", unit="gauge"):
                    try:
                        fu.result()
                    except Exception as e:
                        print(f"[ERROR] gauge task failed: {e}")

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
                for _ in tqdm(as_completed(futures), total=len(futures), desc="    Zipping folders", unit="folder"):
                    pass

        print(f"\n>>> Total time (region {region_i}): {format_elapsed_time(time.time() - start_time_region)}\n")
        gc.collect()