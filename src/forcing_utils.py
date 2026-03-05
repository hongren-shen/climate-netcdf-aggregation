#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2019-2026 Hongren Shen <hongren.shen@uwaterloo.ca>
# SPDX-License-Identifier: LGPL-3.0-or-later

"""
Utilities used across the aggregation workflow.

This file includes:
- Small, stable utilities (timing, zipping, naming, metadata dumping)
- Workflow helpers (temporal wrappers, metadata naming, cleanup)
- A wrapper for external gridweights generator (create_cropping_info)

"""

from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import netCDF4 as nc


def format_elapsed_time(elapsed_time_seconds: float) -> str:
    hours, remainder = divmod(float(elapsed_time_seconds), 3600)
    minutes, seconds = divmod(remainder, 60)
    hours   = int(hours)
    minutes = int(minutes)
    seconds = int(seconds)
    if hours > 0:
        return f"{hours}h {minutes}m {seconds}s"
    if minutes > 0:
        return f"{minutes}m {seconds}s"
    return f"{seconds}s"


def zip_and_cleanup_gauge_folder(
    gauge_folder: str | Path,
    do_zip: bool,
    do_delete: bool,
) -> Optional[str]:
    folder = Path(gauge_folder)
    if not folder.is_dir():
        return None
    if do_zip:
        shutil.make_archive(str(folder), "zip", root_dir=str(folder))
    if do_delete:
        shutil.rmtree(folder)
    return None


def _freq_to_token(freq_str: str) -> str:
    if not freq_str:
        return "unknown"
    s = str(freq_str).upper()

    if s in ("D", "1D"):
        return "24H"

    mD = re.fullmatch(r"(\d+)D", s)
    if mD:
        return f"{int(mD.group(1)) * 24}H"

    mH = re.fullmatch(r"(\d+)?H", s)
    if mH:
        n = int(mH.group(1) or 1)
        return f"{n}H"

    if s.startswith("W"):
        mW = re.fullmatch(r"(\d+)?W.*", s)
        n = int(mW.group(1) or 1) if mW else 1
        return f"{n}W"

    for unit in ("MS", "M", "QS", "Q", "AS", "A", "YS", "Y"):
        if s.endswith(unit):
            base = {"MS": "M", "QS": "Q", "AS": "Y", "YS": "Y"}.get(unit, unit)
            mN = re.fullmatch(r"(\d+)" + re.escape(unit), s)
            n = int(mN.group(1)) if mN else 1
            return f"{n}{base}"

    mT = re.fullmatch(r"(\d+)?T", s)  # minutes
    if mT:
        n = int(mT.group(1) or 1)
        return f"{n}T"

    mS = re.fullmatch(r"(\d+)?S", s)  # seconds
    if mS:
        n = int(mS.group(1) or 1)
        return f"{n}S"

    return s

def _infer_time_range_and_freq(ds: xr.Dataset) -> tuple[str, str]:
    if "time" in ds.coords and ds["time"].size > 0:
        times = pd.to_datetime(ds["time"].values)
        start = pd.to_datetime(times[0]).strftime("%Y%m%d")
        end   = pd.to_datetime(times[-1]).strftime("%Y%m%d")
        time_range = f"{start}_{end}"

        if ds["time"].size > 1:
            inferred = pd.infer_freq(times)
            if inferred:
                freq_token = _freq_to_token(inferred)
            else:
                deltas = pd.Series(times[1:] - times[:-1], dtype="timedelta64[ns]")
                med = deltas.median()
                if pd.isna(med):
                    freq_token = "unknown"
                else:
                    hours = float(med / pd.Timedelta(hours=1))
                    if abs(hours - round(hours)) < 1e-9:
                        freq_token = f"{int(round(hours))}H"
                    else:
                        freq_token = "unknown"
        else:
            time_range = f"{start}_{start}"
            freq_token = "fx"
    else:
        time_range = "fx"
        freq_token = "fx"
    return time_range, freq_token

def _haversine_km(lat1, lon1, lat2, lon2) -> float:
    R = 6371.0
    φ1 = math.radians(float(lat1))
    φ2 = math.radians(float(lat2))
    dφ = φ2 - φ1
    dλ = math.radians(float(lon2) - float(lon1))
    a = math.sin(dφ/2)**2 + math.cos(φ1)*math.cos(φ2)*math.sin(dλ/2)**2
    return 2 * R * math.asin(math.sqrt(a))

def _snap_km_to_nice(km: float, *, tol: float = 0.15) -> Optional[int]:
    nice = [1, 2, 2.5, 3, 4, 5, 8, 9, 10, 12, 12.5, 15, 20, 25, 50, 100]
    nice_sorted = sorted(nice, key=lambda x: (abs(km - x), -1 if x in (10,) else 0))
    for n in nice_sorted:
        if n == 0:
            continue
        if abs(km - n) / n <= tol:
            return int(round(n))
    return None

def _estimate_grid_km_from_coords(ds: xr.Dataset) -> Optional[int]:
    def _pick(obj, *names):
        for n in names:
            if n in obj:
                return obj[n]
        return None

    lat = _pick(ds.coords, "lat", "latitude")
    lon = _pick(ds.coords, "lon", "longitude")
    if lat is None or lon is None:
        lat = _pick(ds.variables, "lat", "latitude")
        lon = _pick(ds.variables, "lon", "longitude")
    if lat is None or lon is None:
        return None

    latv = np.asarray(lat.values)
    lonv = np.asarray(lon.values)

    try:
        if latv.ndim == 1 and lonv.ndim == 1:
            i = len(latv) // 2
            j = len(lonv) // 2
            dy = _haversine_km(latv[max(0, i-1)], lonv[j], latv[i], lonv[j])
            dx = _haversine_km(latv[i], lonv[max(0, j-1)], latv[i], lonv[j])
        else:
            I, J = latv.shape[-2], latv.shape[-1]
            i = max(0, min(I-2, I//2))
            j = max(0, min(J-2, J//2))
            dy = _haversine_km(latv[i, j], lonv[i, j], latv[i+1, j], lonv[i+1, j])
            dx = _haversine_km(latv[i, j], lonv[i, j], latv[i, j+1], lonv[i, j+1])

        km_raw = float(np.nanmedian([dx, dy]))
        if not np.isfinite(km_raw) or km_raw <= 0:
            return None
        snapped = _snap_km_to_nice(km_raw, tol=0.15)
        return snapped if snapped is not None else int(round(km_raw))
    except Exception:
        return None


def _detect_timezone_token(ds: xr.Dataset, tz_hint: Optional[str] = None) -> str:
    if tz_hint:
        t = tz_hint.strip().upper()
        if t == "Z":
            return "UTC"
        if t.startswith("UTC"):
            return t
        m = re.fullmatch(r"([+-])(\d{1,2})(?::?(\d{2}))?$", t)
        if m:
            hh = int(m.group(2)); mm = int(m.group(3) or 0)
            return f"UTC{m.group(1)}{hh}" + (f":{mm:02d}" if mm else "")
        return t

    def _norm(val: str) -> Optional[str]:
        if not val:
            return None
        s = str(val).strip().upper()
        if s == "Z":
            return "UTC"
        if s.startswith("UTC"):
            return s
        m = re.fullmatch(r"([+-])(\d{1,2})(?::?(\d{2}))?$", s)
        if m:
            hh = int(m.group(2)); mm = int(m.group(3) or 0)
            return f"UTC{m.group(1)}{hh}" + (f":{mm:02d}" if mm else "")
        return s

    for k in ("time_zone", "timezone", "time_zone_name", "tz_name"):
        vv = _norm(ds.attrs.get(k))
        if vv:
            return vv

    if "time" in ds.coords:
        for k in ("time_zone", "timezone", "time_zone_name", "tz_name"):
            vv = _norm(ds["time"].attrs.get(k, None))
            if vv:
                return vv

    try:
        units = str(ds["time"].attrs.get("units", "")).strip()
        m = re.search(r"\s([+-])(\d{2})(?::?(\d{2}))?$", units)
        if m:
            sign, hh, mm = m.group(1), int(m.group(2)), int(m.group(3) or 0)
            return f"UTC{sign}{hh}" + (f":{mm:02d}" if mm else "")
        if "UTC" in units.upper():
            return "UTC"
    except Exception:
        pass

    return "UTC"

def auto_generate_nc_filename(
    nc_path: str | Path,
    region_name: str,
    *,
    spatial_token: Optional[str] = None,    # e.g., "agg_subbasin_v2-0-1"
    temporal_token: Optional[str] = None,   # e.g., "24H" (timezone appended) or "24H_UTC-5"
    tz_hint: Optional[str] = None,          # e.g., "UTC", "UTC-5"
) -> str:
    """
    Build a filename from NetCDF metadata:

        <dataset_id>_<var_str>_<region>_<resolution>_<time_range>_<frequency>.nc

    Frequency always includes timezone suffix:
        "<FREQ>_<TZ>" e.g., "1H_UTC", "24H_UTC-5"
    """
    p = Path(nc_path)

    def _infer_station_agg_token(ds: xr.Dataset, path: Path) -> Optional[str]:
        attr = (ds.attrs.get("spatial_aggregation_suffix") or "").strip()
        m = re.search(r"(agg_(?:subbasin|hru|basin)_[A-Za-z0-9\-.]+[ab]?)", attr)
        if m:
            return m.group(1)
        m = re.search(r"(agg_(?:subbasin|hru|basin)_[A-Za-z0-9\-.]+[ab]?)", path.stem)
        if m:
            return m.group(1)
        return None

    def _infer_dataset_id_and_varstr(ds: xr.Dataset) -> tuple[str, str]:
        title = (ds.attrs.get("title") or "").strip()
        if ("CaSR" in title and ("v3.1" in title or "v3_1" in title)):
            dataset_id = "CaSR_v3_1"
        elif "CaSR" in title:
            dataset_id = "CaSR"
        else:
            dataset_id = ds.attrs.get("dataset_id", "unknownDataset")

        candidates = [v for v in ds.data_vars if "time" in ds[v].dims and not v.endswith("_bnds")]
        if not candidates:
            candidates = [v for v in ds.data_vars if not v.endswith("_bnds")]

        if len(candidates) > 5:
            var_str = "allvars"
        else:
            cleaned: List[str] = []
            for v in candidates:
                vv = v
                for pid in (dataset_id, dataset_id.replace("_", ".")):
                    pref = pid + "_"
                    if vv.startswith(pref):
                        vv = vv[len(pref):]
                        break
                cleaned.append(vv)
            var_str = "_".join(cleaned) if cleaned else "unknownVar"

        return dataset_id, var_str

    with xr.open_dataset(p) as ds:
        dataset_id, var_str = _infer_dataset_id_and_varstr(ds)

        time_range, base_freq = _infer_time_range_and_freq(ds)
        tz_token = _detect_timezone_token(ds, tz_hint=tz_hint)

        if temporal_token:
            tmp = str(temporal_token).strip().upper()
            if "_UTC" in tmp:
                frequency = tmp
            else:
                frequency = f"{tmp}_{tz_token}"
        else:
            frequency = f"{base_freq}_{tz_token}"

        if spatial_token:
            resolution = spatial_token
        else:
            agg_token = _infer_station_agg_token(ds, p)
            if agg_token:
                resolution = agg_token
            else:
                resolution = None
                for key in ("grid_resolution_km", "resolution_km", "nominal_resolution_km"):
                    if key in ds.attrs:
                        try:
                            km_attr = float(ds.attrs[key])
                            snapped = _snap_km_to_nice(km_attr, tol=0.20)
                            resolution = f"{snapped if snapped else int(round(km_attr))}km"
                            break
                        except Exception:
                            pass
                if not resolution:
                    km_est = _estimate_grid_km_from_coords(ds)
                    resolution = f"{km_est}km" if km_est is not None else "unknownRes"

        region_name = (region_name or "").strip().replace(" ", "_")
        return f"{dataset_id}_{var_str}_{region_name}_{resolution}_{time_range}_{frequency}.nc"

def get_variable_str(nc_path: str | Path) -> str:
    with xr.open_dataset(nc_path) as ds:
        time_vars = [v for v in ds.data_vars if "time" in ds[v].dims and not v.endswith("_bnds")]
        vars_for_name = time_vars if time_vars else [v for v in ds.data_vars if not v.endswith("_bnds")]
        return "allvars" if len(vars_for_name) > 5 else "_".join(vars_for_name)


def export_nc_metadata_to_ascii(
    nc_path: str | Path,
    output_txt_path: str | Path,
    include_encoding: bool = False,
    *,
    decode_cf: bool = False,
    show_dtypes: bool = True,
) -> None:
    with xr.open_dataset(nc_path, decode_cf=bool(decode_cf)) as ds, open(output_txt_path, "w", encoding="utf-8") as f:
        f.write("=== Global Attributes ===\n")
        for key, val in ds.attrs.items():
            f.write(f"{key}: {val}\n")
        f.write("\n")

        f.write("=== Dimensions ===\n")
        for dim, size in ds.sizes.items():
            f.write(f"{dim}: {size}\n")
        f.write("\n")

        f.write("=== Coordinates ===\n")
        for coord in ds.coords:
            da = ds[coord]
            f.write(f"{coord} {da.dims}\n")
            if show_dtypes:
                f.write(f"    dtype: {str(da.dtype)}\n")
            for attr, val in da.attrs.items():
                f.write(f"    {attr}: {val}\n")
            if include_encoding:
                enc = getattr(da, "encoding", {}) or {}
                enc_clean = {k: v for k, v in enc.items() if v is not None}
                if enc_clean:
                    f.write(f"    encoding: {enc_clean}\n")
        f.write("\n")

        f.write("=== Data Variables ===\n")
        for var in ds.data_vars:
            da = ds[var]
            f.write(f"{var} {da.dims}\n")
            if show_dtypes:
                f.write(f"    dtype: {str(da.dtype)}\n")
            for attr, val in da.attrs.items():
                f.write(f"    {attr}: {val}\n")
            if include_encoding:
                enc = getattr(da, "encoding", {}) or {}
                enc_clean = {k: v for k, v in enc.items() if v is not None}
                if enc_clean:
                    f.write(f"    encoding: {enc_clean}\n")
            f.write("\n")


def _infer_stat_label_from_attrs(ds: nc.Dataset) -> Optional[str]:
    def _map(v: str) -> Optional[str]:
        v = (v or "").strip().lower()
        if v == "mean": return "AVG"
        if v == "min":  return "MIN"
        if v == "max":  return "MAX"
        return None

    lbl = _map(getattr(ds, "temporal_aggregation_method", ""))
    if lbl:
        return lbl

    for vn in ds.variables:
        var = ds.variables[vn]
        lbl = _map(getattr(var, "temporal_aggregation_method", ""))
        if lbl:
            return lbl
        cm = (getattr(var, "cell_methods", "") or "").lower()
        if "time:" in cm:
            if " mean" in cm: return "AVG"
            if " min"  in cm: return "MIN"
            if " max"  in cm: return "MAX"
    return None


def infer_stat_label_from_name(p: Path) -> Optional[str]:
    m = re.search(r"_(AVG|MIN|MAX)_", p.name, flags=re.IGNORECASE)
    return (m.group(1).upper() if m else None)


def export_metadata(
    nc_path: Path,
    name_key: str,
    out_dir: Path,
    *,
    stat_hint: Optional[str] = None,
) -> None:
    """
    Write metadata as:
      netcdf_metadata_<name_key>_<var>[_agg_...][_STAT].txt
    """
    if not nc_path or not nc_path.exists():
        return

    try:
        var_token = get_variable_str(nc_path)
        if not var_token:
            raise ValueError("empty var_token")
    except Exception:
        stem = nc_path.stem
        m_before_gid = re.split(rf"_{re.escape(name_key)}_", stem, maxsplit=1)
        m_before_agg = re.split(r"_agg_", stem, maxsplit=1)
        if len(m_before_gid) > 1:
            var_token = m_before_gid[0]
        elif len(m_before_agg) > 1:
            var_token = m_before_agg[0]
        else:
            var_token = stem

    agg_tail = ""
    m = re.search(r"(_agg_.*)$", nc_path.stem)
    if m:
        agg_tail = m.group(1)

    stat_label: Optional[str] = stat_hint or infer_stat_label_from_name(nc_path)
    if stat_label is None:
        try:
            with nc.Dataset(str(nc_path), "r") as ds:
                stat_label = _infer_stat_label_from_attrs(ds)
        except Exception:
            stat_label = None

    meta_base = f"netcdf_metadata_{name_key}_{var_token}{agg_tail}"
    if stat_label:
        meta_base += f"_{stat_label}"
    meta_path = out_dir / (meta_base + ".txt")
    export_nc_metadata_to_ascii(nc_path, meta_path)


def create_cropping_info(
    grid_weights_generator_path: str | Path,
    shapefile_path: str | Path,
    nc_path: str | Path,
    grid_weights_path: str | Path,
    grid_info_path: str | Path | None = None,
    *,
    key_col: str = "Id",
    verbose: bool = True,
    generate_json: bool = True,
) -> None:
    """
    Calls your external gridweights generator script to produce a Raven gridweights file.
    """
    shapefile_path = str(shapefile_path)
    nc_path = str(nc_path)
    grid_weights_path = str(grid_weights_path)

    if not os.path.exists(shapefile_path):
        return

    temp_grid_info = False
    if grid_info_path is None:
        temp_grid_info = True
        grid_info_path = os.path.join(os.path.dirname(grid_weights_path), "_temp_grid_info.txt")
    else:
        grid_info_path = str(grid_info_path)

    dimname = "rlon,rlat"
    varname = "lon,lat"

    cmd = [
        "python", str(grid_weights_generator_path),
        "-i", nc_path,
        "-d", dimname,
        "-v", varname,
        "-r", shapefile_path,
        "-a", "--doall",
        "-c", key_col,
        "-o", grid_weights_path,
        "-o2", grid_info_path,
    ]
    if generate_json:
        cmd.extend(["-j", "--dojson"])

    try:
        if verbose:
            subprocess.run(cmd, check=True)
        else:
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    finally:
        if temp_grid_info and os.path.exists(grid_info_path):
            os.remove(grid_info_path)


def replicate_subbasin_coverage_to_hrus(
    subbasin_gridweights_file: str | Path,
    shapefile_path: str | Path,
    output_file: str | Path,
    subbasin_col: str = "SubId",
    hru_col: str = "HRU_ID"
) -> None:
    """
    Reads a gridweights file keyed by subbasin_col.

    - If there's exactly one unique subbasin id in it, replicate its (Cell,Weight)
      to every HRU in shapefile_path (basin-lumped coverage).
    - Otherwise, replicate each subbasin's coverage to the HRUs within that subbasin.

    The rows are sorted by HRU, then by Cell number.
    """
    shapefile_path = Path(shapefile_path)
    if not shapefile_path.exists():
        return

    coverage_data = []
    num_cells = 0
    in_block = False
    with open(subbasin_gridweights_file, encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s.startswith(":NumberGridCells"):
                try:
                    num_cells = int(s.split()[-1])
                except Exception:
                    num_cells = 0
            elif s.startswith(":GridWeights"):
                in_block = True
            elif s.startswith(":EndGridWeights"):
                in_block = False
                break
            elif in_block and s and not s.startswith("#"):
                parts = re.split(r"\s+", s)
                if len(parts) >= 3:
                    sid = int(float(parts[0]))
                    cell = int(float(parts[1]))
                    w = float(parts[2])
                    coverage_data.append((sid, cell, w))

    if not coverage_data:
        return

    df = pd.DataFrame(coverage_data, columns=[subbasin_col, "Cell", "Weight"])

    gdf = gpd.read_file(shapefile_path)
    hrus = sorted(gdf[hru_col].unique())

    new_rows = []
    unique_subs = df[subbasin_col].unique()
    if len(unique_subs) == 1:
        for _, cell, w in df.itertuples(index=False):
            for h in hrus:
                new_rows.append((h, cell, w))
    else:
        sb_to_hrus = gdf.groupby(subbasin_col)[hru_col].unique().to_dict()
        df2 = df.drop_duplicates([subbasin_col, "Cell"])
        for sid, cell, w in df2.itertuples(index=False):
            for h in sb_to_hrus.get(sid, []):
                new_rows.append((h, cell, w))

    new_rows.sort(key=lambda x: (x[0], x[1]))

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w", encoding="utf-8") as out:
        out.write(":GridWeights\n")
        out.write("   #\n")
        out.write("   # [# HRUs]\n")
        out.write(f"     :NumberHRUs       {len(hrus)}\n")
        if num_cells > 0:
            out.write(f"     :NumberGridCells  {num_cells}\n")
        out.write("   #\n")
        out.write("   # [HRU ID] [Cell #] [w_kl]\n")
        for hru, cell, w in new_rows:
            out.write(f"     {int(float(hru))}   {int(cell)}   {w}\n")
        out.write(":EndGridWeights\n")

# Workflow helpers
def is_daily_freq(freq: str) -> bool:
    s = (freq or "").upper()
    return s in ("D", "1D", "24H")

def is_temperature_var(var_file_key: str | None) -> bool:
    s = (var_file_key or "").upper()
    return any(tok in s for tok in ("TT", "T2M", "TEMP"))

def infer_hru_suffix(hru_version: str) -> str:
    return "a" if str(hru_version).strip().lower().startswith("semi") else "b"

def infer_hru_suffix_from_agg_file(nc_path: Path) -> str:
    """
    Try to infer 'a'/'b' suffix from file or global attrs.
    Returns 'a' if not detected.
    """
    try:
        with nc.Dataset(str(nc_path), "r") as ds:
            suf = (ds.__dict__.get("spatial_aggregation_suffix") or "") + " " + nc_path.name
    except Exception:
        suf = nc_path.name
    m = re.search(r"agg_hru_[^ab]*([ab])", suf)
    return (m.group(1) if m else "a")

def discover_precomputed_region_agg(
    output_region_forcing_folder: Path,
    var_file_key: str,
    *,
    olrrp_version: str,
) -> Dict[str, List[Path]]:
    """
    Search the region archive for precomputed aggregated NetCDFs for this variable.

    Returns dict with keys:
      basin, subbasin, hru, hru_a, hru_b
    """
    results: Dict[str, List[Path]] = {k: [] for k in ["basin", "subbasin", "hru", "hru_a", "hru_b"]}
    stem_hint = Path(var_file_key).stem

    pats = [
        f"*{stem_hint}*agg_basin_{olrrp_version}*.nc",
        f"*{stem_hint}*agg_subbasin_{olrrp_version}*.nc",
        f"*{stem_hint}*agg_hru_{olrrp_version}a*.nc",
        f"*{stem_hint}*agg_hru_{olrrp_version}b*.nc",
        f"*{stem_hint}*agg_hru_{olrrp_version}*.nc",
    ]

    if not output_region_forcing_folder.exists():
        return results

    for p in pats:
        for f in output_region_forcing_folder.glob(p):
            name = f.name
            if "_agg_basin_" in name:
                results["basin"].append(f)
            elif "_agg_subbasin_" in name:
                results["subbasin"].append(f)
            elif "_agg_hru_" in name:
                if name.endswith("a.nc") or f"{olrrp_version}a" in name:
                    results["hru_a"].append(f)
                elif name.endswith("b.nc") or f"{olrrp_version}b" in name:
                    results["hru_b"].append(f)
                else:
                    results["hru"].append(f)

    return results


def _infer_temporal_reduce(var_file_key: str) -> str:
    """
    Decide temporal reducer from filename/token:
      - precip-like -> sum
      - temperature-like -> mean
      - default -> mean
    """
    s = (var_file_key or "").upper()
    precip_tokens = ("PR", "PR0", "RAIN", "SNOW", "PRECIP", "PRATE")
    if any(tok in s for tok in precip_tokens):
        return "sum"
    temp_tokens = ("TT", "T2M", "TAVG", "TMEAN", "TEMP")
    if any(tok in s for tok in temp_tokens):
        return "mean"
    return "mean"


def _insert_label_after_var_token(filename: str, name_key: str, label: str) -> str:
    if f"_{label}_" in filename.upper():
        return filename
    tag = f"_{name_key}_"
    pos = filename.find(tag)
    if pos > 0:
        return filename[:pos] + f"_{label}" + filename[pos:]
    pos2 = filename.find("_agg_")
    if pos2 > 0:
        return filename[:pos2] + f"_{label}" + filename[pos2:]
    stem, ext = os.path.splitext(filename)
    return f"{stem}_{label}{ext}"


def temporal_aggregate(
    pre_nc: Path,
    *,
    region_name: str,
    var_file_key: str,
    temporal_freq: str,
    temporal_anchor_hour: int,
    temporal_tz_name: str,
    period_ending2starting: bool,
    compress: bool,
    compression_level: int,
    output_time_chunk_hours: int,
    overwrite: bool,
) -> Path:
    """
    Default temporal aggregation (reduce inferred from var_file_key).
    """
    # lazy import
    from temporal_aggregation import aggregate_netcdf_temporal

    # tz: allow IANA or fixed UTC±H
    tz_name_arg: str | None = None
    tz_offset_arg: float | None = None
    tz_cfg = temporal_tz_name or "UTC"
    if "/" in tz_cfg:
        tz_name_arg = tz_cfg
    else:
        s = tz_cfg.upper()
        if s.startswith("UTC"):
            try:
                tz_offset_arg = float(s.replace("UTC", ""))
            except Exception:
                tz_offset_arg = None

    reduce_method = _infer_temporal_reduce(var_file_key)

    freq_raw = (temporal_freq or "").upper()
    temporal_token = "24H" if freq_raw in ("D", "1D") else freq_raw

    out_name = auto_generate_nc_filename(
        pre_nc,
        region_name,
        temporal_token=temporal_token,
        tz_hint=temporal_tz_name or "UTC",
    )
    out_nc = pre_nc.with_name(out_name)

    if out_nc.exists() and not overwrite:
        return out_nc

    result_path = aggregate_netcdf_temporal(
        src_nc_path=pre_nc,
        out_nc_path=out_nc,
        temporal_reduce=reduce_method,
        temporal_freq=temporal_freq,
        temporal_anchor_hour=int(temporal_anchor_hour or 0),
        tz_name=tz_name_arg,
        tz_offset_hours=tz_offset_arg,
        period_is_ending=period_ending2starting,
        period_step_hours=None,
        time_name="time",
        compress=compress,
        compression_level=int(compression_level),
        output_time_chunk_hours=int(output_time_chunk_hours),
        overwrite=True,
        verbose=False,
    )
    return Path(result_path)


def temporal_aggregate_specific(
    pre_nc: Path,
    *,
    region_name: str,
    temporal_freq: str,
    temporal_anchor_hour: int,
    temporal_tz_name: str,
    period_ending2starting: bool,
    compress: bool,
    compression_level: int,
    output_time_chunk_hours: int,
    overwrite: bool,
    reduce_method: str,
    label: str,
) -> Optional[Path]:
    """
    Temporal aggregation with a specific reducer and an explicit label token.
    """
    if not pre_nc.exists():
        return None

    # lazy import
    from temporal_aggregation import aggregate_netcdf_temporal

    freq_raw = (temporal_freq or "").upper()
    temporal_token = "24H" if freq_raw in ("D", "1D") else freq_raw

    out_name_base = auto_generate_nc_filename(
        pre_nc,
        region_name,
        temporal_token=temporal_token,
        tz_hint=temporal_tz_name or "UTC",
    )
    out_nc_base = pre_nc.with_name(out_name_base)
    out_nc_labeled = out_nc_base.with_name(_insert_label_after_var_token(out_nc_base.name, region_name, label))

    if out_nc_labeled.exists() and not overwrite:
        return out_nc_labeled

    tz_name_arg: str | None = None
    tz_offset_arg: float | None = None
    tz_cfg = temporal_tz_name or "UTC"
    if "/" in tz_cfg:
        tz_name_arg = tz_cfg
    else:
        s = tz_cfg.upper()
        if s.startswith("UTC"):
            try:
                tz_offset_arg = float(s.replace("UTC", ""))
            except Exception:
                tz_offset_arg = None

    tmp = out_nc_labeled.with_name(out_nc_labeled.stem + "__tmp.nc")
    result_path = aggregate_netcdf_temporal(
        src_nc_path=pre_nc,
        out_nc_path=tmp,
        temporal_reduce=reduce_method,
        temporal_freq=temporal_freq,
        temporal_anchor_hour=int(temporal_anchor_hour or 0),
        tz_name=tz_name_arg,
        tz_offset_hours=tz_offset_arg,
        period_is_ending=period_ending2starting,
        period_step_hours=None,
        time_name="time",
        compress=compress,
        compression_level=int(compression_level),
        output_time_chunk_hours=int(output_time_chunk_hours),
        overwrite=True,
        verbose=False,
    )
    result_path = Path(result_path)

    if out_nc_labeled.exists() and overwrite:
        out_nc_labeled.unlink(missing_ok=True)
    result_path.replace(out_nc_labeled)

    try:
        tmp.unlink(missing_ok=True)
    except Exception:
        pass

    return out_nc_labeled


def rename_with_label(out_nc: Path, name_key: str, label: str, *, overwrite: bool) -> Path:
    if not out_nc.exists():
        return out_nc
    labeled_name = _insert_label_after_var_token(out_nc.name, name_key, label)
    if labeled_name == out_nc.name:
        return out_nc
    labeled_path = out_nc.with_name(labeled_name)
    if labeled_path.exists() and overwrite:
        labeled_path.unlink(missing_ok=True)
    if labeled_path.exists() and not overwrite:
        return labeled_path
    out_nc.replace(labeled_path)
    return labeled_path


def cleanup_intermediates_in_gauge_dir(gauge_dir: Path) -> int:
    """
    Remove per-gauge intermediate files, keeping final outputs.
    Returns count removed.
    """
    if not gauge_dir.exists():
        return 0

    patterns = [
        "_gridinfo_*.txt",
        "_gridweights_*.txt",
        "gridweights_for_crop_cached.txt",
        ".gridded_weights_done",
    ]
    removed = 0
    for patt in patterns:
        for p in gauge_dir.glob(patt):
            try:
                p.unlink(missing_ok=True)
                removed += 1
            except Exception:
                pass

    # shared crop cache
    try:
        cache_dir = gauge_dir.parent / ".crop_cache"
        cache_file = cache_dir / f"{gauge_dir.name}.txt"
        if cache_file.exists():
            cache_file.unlink()
        if cache_dir.exists():
            try:
                next(cache_dir.iterdir())
            except StopIteration:
                cache_dir.rmdir()
    except Exception:
        pass

    return removed