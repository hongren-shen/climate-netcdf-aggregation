#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2019-2026 Hongren Shen <hongren.shen@uwaterloo.ca>
# SPDX-License-Identifier: LGPL-3.0-or-later

"""
Aggregate cropped gauge NetCDFs from grid cells to a chosen target extent:
  - "basin"
  - "subbasin"
  - "hru"
using HRU layout 'a' (semi-distributed) or 'b' (routing-only).

Output filename pattern: <auto_core>_<suffix_str>.nc
  - basin    - *_agg_basin_<olrrp_version>.nc
  - subbasin - *_agg_subbasin_<olrrp_version>.nc
  - hru      - *_agg_hru_<olrrp_version><a|b>.nc

Dimensions (raven required format)
----------
- time
- stations
- char_leng_id = 64
- string30     = 30

Variables (raven required format)
-----------------
- station_id    : ('stations','char_leng_id')  fixed-width |S1
- station_names : ('stations','string30')      fixed-width |S1
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Iterable, NamedTuple, List, Optional, Tuple, Dict
import math

import geopandas as gpd
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr

from forcing_utils import (
    auto_generate_nc_filename,
    create_cropping_info,
    export_nc_metadata_to_ascii,
    replicate_subbasin_coverage_to_hrus,
    normalize_hru_version,
)

TIME_NAME = "time"

DIM_CHAR_ID   = "char_leng_id"
DIM_STRING30  = "string30"
LEN_CHAR_ID   = 64
LEN_STRING30  = 30
DEFAULT_NUMERIC_FILL = -9999.0
_CITATION_PATH = Path(__file__).resolve().parents[1] / "CITATION.cff"


def _software_creator() -> str:
    """Read the release version at write time; copied scripts can omit the CFF."""
    name = "climate-netcdf-aggregation"
    try:
        citation = _CITATION_PATH.read_text(encoding="utf-8-sig")
    except (OSError, UnicodeError):
        return name

    # Read only the top-level scalar, without requiring a YAML dependency.
    lines = [line for line in citation.splitlines() if line.startswith("version:")]
    if len(lines) != 1:
        return name
    match = re.fullmatch(
        r'''version:[ \t]*(?P<quote>["']?)'''
        r'''(?P<version>[vV]?[0-9][0-9A-Za-z._+-]*)(?P=quote)'''
        r'''(?:[ \t]+#.*)?[ \t]*''',
        lines[0],
    )
    if match is None:
        return name
    version = match.group("version")
    label = version if version.startswith(("v", "V")) else f"v{version}"
    return f"{name} {label}"


def _sanitize_fill_value(
    val,
    *,
    coerce_to_default: bool = False,
    default: float = DEFAULT_NUMERIC_FILL,
):
    """
    Return a finite float fill value or None.
    - If val is NaN/inf or non-numeric: return (default if coerce_to_default else None).
    - If val is a finite number: return float(val).
    """
    if val is None:
        return (default if coerce_to_default else None)
    try:
        f = float(val)
    except Exception:
        return (default if coerce_to_default else None)
    if math.isnan(f) or not math.isfinite(f):
        return (default if coerce_to_default else None)
    return f
# ------------------------------------------------------

class Version(NamedTuple):
    idx: int            # version number (1,2,3)
    target: str         # "basin" | "subbasin" | "hru"
    desc: str           # description
    key_col: str        # shapefile key column


VERSIONS: Tuple[Version, ...] = (
    Version(1, "basin",    "All HRUs share the same grid weights derived from the lumped basin", "ID"),
    Version(2, "subbasin", "HRUs within the same subbasin share identical grid weights",         "SubId"),
    Version(3, "hru",      "Each HRU has its own unique grid weights",                           "HRU_ID"),
)

VERSIONS_BY_IDX     = {v.idx: v for v in VERSIONS}
VERSIONS_BY_TARGET  = {v.target: v for v in VERSIONS}

__all__ = [
    "aggregate_netcdf_forcing",
    "normalize_hru_version",
    "Version",
    "VERSIONS",
    "VERSIONS_BY_IDX",
    "VERSIONS_BY_TARGET",
    "unit_ids_for_gauge",
    "derive_gridweights_for_gridded_products",
    "extract_gauge_from_regional_agg_nc",
    "crop_netcdf_by_cell_info",
    "crop_netcdf_time_range",
    "crop_gridded_nc_for_gauge",
]


def _lazy_key_cols() -> Dict[str, str]:
    try:
        return {
            "basin":    VERSIONS_BY_TARGET["basin"].key_col,
            "subbasin": VERSIONS_BY_TARGET["subbasin"].key_col,
            "hru":      VERSIONS_BY_TARGET["hru"].key_col,
        }
    except Exception:
        return {"basin": "ID", "subbasin": "SubId", "hru": "HRU_ID"}


def _hdr_kv(key: str, val, width: int = 17) -> str:
    return f"# {key:>{width}}: {'' if val is None else str(val)}"


def _gridweights_header(nc_path: Path, created_on: str) -> List[str]:
    """Describe the referenced NetCDF using its own dataset metadata."""
    dataset_label = ""
    with nc.Dataset(nc_path, "r") as dataset:
        for attribute in ("title", "source"):
            value = getattr(dataset, attribute, None)
            if isinstance(value, str) and value.strip():
                dataset_label = " ".join(value.split())
                break

    lines = [
        "# ---------------------------------------------------------------------------------------",
        "#         File type: Grid Weights for Raven NETCDF v4.n",
    ]
    if dataset_label:
        lines.append(_hdr_kv("Forcing dataset", dataset_label))
    lines.extend([
        _hdr_kv("Created by", _software_creator()),
        "#        Repository: https://github.com/hongren-shen/climate-netcdf-aggregation",
        _hdr_kv("Creation date", created_on),
        "# ---------------------------------------------------------------------------------------",
    ])
    return lines


def unit_ids_for_gauge(
    target: str,
    shp_dir: Path,
    *,
    outline_shp_name: str = "outline.shp",
    subbasin_name_prefix: str = "finalcat_info_",
    hru_name_prefix: str = "final_hru_",
    olrrp_version: str = "v2-0-1",
    key_basin: Optional[str] = None,
    key_subbasin: Optional[str] = None,
    key_hru: Optional[str] = None,
    hru_suffix: Optional[str] = None,
) -> List[int]:
    keys = {"basin": key_basin, "subbasin": key_subbasin, "hru": key_hru}
    if any(v is None for v in keys.values()):
        defaults = _lazy_key_cols()
        keys = {
            "basin": key_basin or defaults["basin"],
            "subbasin": key_subbasin or defaults["subbasin"],
            "hru": key_hru or defaults["hru"],
        }

    target = target.strip().lower()
    if target == "basin":
        shp = shp_dir / outline_shp_name
        key = keys["basin"]
    elif target == "subbasin":
        shp = shp_dir / f"{subbasin_name_prefix}{olrrp_version}.shp"
        key = keys["subbasin"]
    elif target == "hru":
        suffix = (hru_suffix or "a").strip().lower()
        if suffix not in ("a", "b"):
            suffix = "a"
        shp = shp_dir / f"{hru_name_prefix}{olrrp_version}{suffix}.shp"
        key = keys["hru"]
    else:
        raise ValueError(f"Unknown target: {target}")

    gdf = gpd.read_file(shp)
    ser = pd.to_numeric(gdf[key], errors="coerce").dropna()
    return ser.astype(int).tolist()


def _to_char_matrix(strings: Iterable[str], width: int) -> np.ndarray:
    strings = ["" if s is None else str(s) for s in strings]
    out = np.full((len(strings), width), b" ", dtype="S1")
    for i, s in enumerate(strings):
        bs = s.encode("ascii", "ignore")[:width]
        if bs:
            out[i, : len(bs)] = np.frombuffer(bs, dtype="S1")
    return out


def _select_spatial_dims(da: xr.DataArray) -> Tuple[str, str]:
    for yx in (("rlat", "rlon"), ("y", "x"), ("lat", "lon")):
        if all(d in da.dims for d in yx):
            return yx
    raise ValueError(f"No known 2D spatial dims in {da.name}: {da.dims}")


def _load_gridinfo_with_weights(grid_info_path: Path):
    df = pd.read_csv(grid_info_path, sep=r"\s+", engine="python")
    cmap = {c.lower().replace("_", ""): c for c in df.columns}

    unit_col = None
    for cand in ("hruid", "subid", "id", "unitid", "unit", "obsnm"):
        if cand in cmap:
            unit_col = cmap[cand]
            break
    if unit_col is None:
        raise ValueError(f"No unit id column found in grid-info: {grid_info_path}")

    row_col = cmap.get("row"); col_col = cmap.get("col")
    if row_col is None or col_col is None:
        raise ValueError(f"Grid-info must have Row and Col columns: {grid_info_path}")

    w_col = None
    for cand in ("wkl", "weight", "w"):
        if cand in cmap:
            w_col = cmap[cand]
            break
    if w_col is None:
        raise ValueError(f"No weight column (w_kl/Weight/w) in {grid_info_path}")

    unit_ids = pd.to_numeric(df[unit_col], errors="coerce")
    rows     = pd.to_numeric(df[row_col], errors="coerce")
    cols     = pd.to_numeric(df[col_col], errors="coerce")
    weights  = pd.to_numeric(df[w_col],   errors="coerce")

    m = unit_ids.notna() & rows.notna() & cols.notna() & weights.notna()
    unit_ids = unit_ids[m].astype(int).to_numpy()
    rows     = rows[m].astype(int).to_numpy()
    cols     = cols[m].astype(int).to_numpy()
    weights  = weights[m].astype(float).to_numpy()
    return unit_ids, rows, cols, weights


def _get_dataset_id_and_varstr(nc_path: Path) -> tuple[str, str]:
    with xr.open_dataset(nc_path) as ds:
        title = ds.attrs.get("title", "") or ""
        dataset_id = ds.attrs.get("dataset_id", "").strip()
        candidates = [v for v in ds.data_vars if TIME_NAME in ds[v].dims and not v.endswith("_bnds")]
        if not candidates:
            candidates = [v for v in ds.data_vars if not v.endswith("_bnds")]

        if not dataset_id:
            version_pattern = r"^(CaSR_v\d+(?:_\d+)*)_"
            variable_ids = {
                match.group(1)
                for name in candidates
                if (match := re.match(version_pattern, name))
            }
            filename_match = re.match(version_pattern, Path(nc_path).stem)
            if len(variable_ids) == 1:
                dataset_id = next(iter(variable_ids))
            elif not variable_ids and filename_match:
                dataset_id = filename_match.group(1)
            elif variable_ids or "casr" in title.lower():
                dataset_id = "CaSR"
            else:
                first_var = next(iter(ds.data_vars), "")
                m = re.match(r"([A-Za-z0-9_.]+)_", first_var)
                dataset_id = (m.group(1) if m else "unknownDataset")

        if len(candidates) > 5:
            var_str = "allvars"
        else:
            cleaned = []
            dids = {dataset_id, dataset_id.replace("_", ".")}
            for v in candidates:
                vv = v
                for pid in dids:
                    pref = pid + "_"
                    if vv.startswith(pref):
                        vv = vv[len(pref):]
                        break
                cleaned.append(vv)
            var_str = "_".join(cleaned)

    return dataset_id, var_str


def _read_time_from_source(nc_path: Path):
    """Return (values, units, calendar, kind) for 'time' from the source NetCDF."""
    with nc.Dataset(nc_path, "r") as src:
        if TIME_NAME not in src.variables:
            return None, None, None, "missing"
        tvar = src.variables[TIME_NAME]
        units = getattr(tvar, "units", None)
        calendar = getattr(tvar, "calendar", "proleptic_gregorian")
        if tvar.dtype.kind in ("i", "u", "f"):
            vals = np.array(tvar[:], dtype=np.float64)
            return vals, units, calendar, "numeric"
        return None, units, calendar, "datetime64"


def _append_history_line(ds, line: str) -> None:
    """
    Append a single line to the global 'history' for either xarray.Dataset or netCDF4.Dataset.
    """
    if not line:
        return
    # read current
    if hasattr(ds, "attrs"):                     # xarray.Dataset
        prev = ds.attrs.get("history", "")
        ds.attrs["history"] = prev + ("\n" if prev else "") + line
    elif hasattr(ds, "getncattr"):               # netCDF4.Dataset
        try:
            prev = ds.getncattr("history")
        except Exception:
            prev = ""
        ds.setncattr("history", prev + ("\n" if prev else "") + line)
    else:                                        # last resort
        try:
            prev = getattr(ds, "history", "")
        except Exception:
            prev = ""
        try:
            setattr(ds, "history", prev + ("\n" if prev else "") + line)
        except Exception:
            pass

try:
    from numba import njit, prange
    _HAVE_NUMBA = True
except Exception:
    _HAVE_NUMBA = False

if _HAVE_NUMBA:
    @njit(parallel=True, fastmath=True)
    def _reduce_dense_numba(V: np.ndarray, unit_idx: np.ndarray, weights: np.ndarray, U: int) -> np.ndarray:
        T, K = V.shape
        out = np.zeros((T, U), dtype=V.dtype)
        for t in prange(T):
            for k in range(K):
                out[t, unit_idx[k]] += V[t, k] * weights[k]
        return out

def _weighted_reduce_dense(V: np.ndarray, unit_idx: np.ndarray, weights: np.ndarray, U: int) -> np.ndarray:
    V = V.astype(np.float64, copy=False)
    unit_idx = unit_idx.astype(np.int64, copy=False)
    weights = weights.astype(np.float64, copy=False)
    if _HAVE_NUMBA:
        return _reduce_dense_numba(V, unit_idx, weights, U)
    T, _K = V.shape
    res = np.zeros((T, U), dtype=V.dtype)
    VW = V * weights[None, :]
    chunk = 24 * 7 * 8
    for t0 in range(0, T, chunk):
        t1 = min(T, t0 + chunk)
        for t in range(t0, t1):
            res[t, :] = np.bincount(unit_idx, weights=VW[t, :], minlength=U)
    return res


def _read_unit_ids_from_nc_netcdf4(nc_path: Path) -> np.ndarray:
    with nc.Dataset(nc_path, "r") as ds:
        if "station_id" not in ds.variables:
            raise KeyError("station_id variable missing in aggregated file.")
        ch = ds.variables["station_id"]  # stations x char_leng_id
        stations = ch.shape[0]
        ids: list[int] = []
        for i in range(stations):
            row = ch[i, :].tobytes().decode("ascii", "ignore").strip()
            if row:
                try:
                    ids.append(int(row))
                except Exception:
                    raise ValueError(f"Non-numeric station_id at index {i}: {row!r}")
        if not ids:
            raise ValueError("No station IDs decoded from station_id.")
        return np.array(ids, dtype=int)


def _write_version_gridweights_for_aggregated(
    out_txt: Path,
    agg_nc_path: Path,
    version: int,
    name_key: str | None,
    olrrp_version: str,
    hru_name_prefix: str,
    suffix_letter: str,
    shapefile_root: Path,
    *,
    common_header: str | None = None,
) -> None:
    ver = VERSIONS_BY_IDX.get(version)
    if ver is None:
        raise ValueError(f"Unsupported version: {version}")

    unit_ids = _read_unit_ids_from_nc_netcdf4(agg_nc_path)
    index_by_unitid = {int(unit_ids[i]): i for i in range(len(unit_ids))}

    outline_shp  = shapefile_root / "outline.shp"
    subbasin_shp = shapefile_root / f"finalcat_info_{olrrp_version}.shp"
    hru_shp      = shapefile_root / f"{hru_name_prefix}{olrrp_version}{suffix_letter}.shp"

    shp_by_target = {"basin": outline_shp, "subbasin": subbasin_shp, "hru": hru_shp}
    shp = shp_by_target[ver.target]
    key_col = ver.key_col
    if not shp.exists():
        raise FileNotFoundError(f"Missing shapefile: {shp}")

    gdf = gpd.read_file(hru_shp)   # still need HRU shapefile for building pairs
    cmap = {c.lower().replace("_", ""): c for c in gdf.columns}
    hru_col = cmap.get("hruid") or cmap.get("hru id") or "HRU_ID"
    if hru_col not in gdf.columns:
        raise KeyError(f"HRU_ID column not found in {hru_shp}")

    pairs: list[tuple[int, int, float]] = []
    if ver.idx == 1:  # basin
        hrus = pd.to_numeric(gdf[hru_col], errors="coerce").dropna().astype(int).unique().tolist()
        hrus.sort()
        for hru in hrus:
            pairs.append((hru, 0, 1.0))
    elif ver.idx == 2:  # subbasin
        sub_col = cmap.get("subid") or cmap.get("sub_id") or "SubId"
        if sub_col not in gdf.columns:
            raise KeyError(f"SubId column not found in {hru_shp}")
        df_sub = gdf[[hru_col, sub_col]].dropna()
        for _, row in df_sub.iterrows():
            hru = int(row[hru_col]); sub = int(row[sub_col])
            pos = index_by_unitid.get(sub)
            if pos is not None:
                pairs.append((hru, pos, 1.0))
        pairs.sort(key=lambda t: (t[0], t[1]))
    else:  # hru
        hrus = pd.to_numeric(gdf[hru_col], errors="coerce").dropna().astype(int).unique().tolist()
        hrus.sort()
        for hru in hrus:
            pos = index_by_unitid.get(hru)
            if pos is not None:
                pairs.append((hru, pos, 1.0))

    created_on = datetime.now().strftime("%Y-%m-%d")

    if common_header:
        header_lines = common_header.strip().splitlines()
    else:
        header_lines = _gridweights_header(agg_nc_path, created_on)

    header_lines += [_hdr_kv(f"Version {ver.idx}", ver.desc)]
    if name_key and str(name_key).strip():
        header_lines.append(_hdr_kv("Name key", str(name_key).strip()))

    header_lines += [
        _hdr_kv("Weights Shapefile", shp.name),
        _hdr_kv("key_col", key_col),
        _hdr_kv("NETCDF", agg_nc_path.name),
        "#              Note: This grid weights file uses pre-aggregated forcing for each HRU",
        "# ---------------------------------------------------------------------------------------",
        "",
        ":GridWeights",
        "   #",
        "   # [# HRUs]",
        f"     :NumberHRUs       {len(pairs)}",
        f"     :NumberStations   {len(unit_ids)}",
        "   #",
        "   # [HRU ID] [station index (0-based)] [weight]",
    ]

    out_txt.parent.mkdir(parents=True, exist_ok=True)
    with open(out_txt, "w", encoding="utf-8") as fp:
        fp.write("\n".join(header_lines) + "\n")
        for hru, cell_idx, w in pairs:
            fp.write(f"     {hru}   {cell_idx}   {w}\n")
        fp.write(":EndGridWeights\n")


def derive_gridweights_for_gridded_products(
    *,
    nc_path: Path,
    shp_dir: Path,
    out_dir: Path,
    grid_weights_generator: Path,
    olrrp_version: str,
    hru_name_prefix: str,
    subbasin_name_prefix: str,
    outline_shp_name: str,
    grid_weights_file_prefix: str = "gridweights_",
    common_header: Optional[str] = None,
    label: Optional[str] = None,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    outline_shp   = shp_dir / outline_shp_name
    subbasin_shp  = shp_dir / f"{subbasin_name_prefix}{olrrp_version}.shp"
    hru_shp_a     = shp_dir / f"{hru_name_prefix}{olrrp_version}a.shp"
    hru_shp_b     = shp_dir / f"{hru_name_prefix}{olrrp_version}b.shp"

    key_basin    = VERSIONS_BY_TARGET["basin"].key_col
    key_subbasin = VERSIONS_BY_TARGET["subbasin"].key_col
    key_hru      = VERSIONS_BY_TARGET["hru"].key_col

    created_on = datetime.now().strftime("%Y-%m-%d")
    default_header = [] if common_header else _gridweights_header(nc_path, created_on)

    def _prepend_header(path: Path, *, version_idx: int, version_desc: str, weights_shp_name: str, key_col: str) -> None:
        lines: List[str] = []
        if common_header:
            lines.extend(common_header.strip().splitlines())
        else:
            lines.extend(default_header)
        lines.append(_hdr_kv(f"Version {version_idx}", version_desc))
        if label and str(label).strip():
            lines.append(_hdr_kv("Name key", str(label).strip()))
        lines += [
            _hdr_kv("Weights Shapefile", weights_shp_name),
            _hdr_kv("key_col",           key_col),
            _hdr_kv("NETCDF",            nc_path.name),
            "#              Note: This grid weights file uses gridded forcing for each HRU",
            "# ---------------------------------------------------------------------------------------",
            "",
        ]
        with open(path, "r+", encoding="utf-8") as fp:
            content = fp.read()
            fp.seek(0)
            fp.write("\n".join(lines) + "\n" + (content or ""))

    for suffix, hru_shp in (("a", hru_shp_a), ("b", hru_shp_b)):
        if not hru_shp.exists():
            continue

        scratch_basin = out_dir / f"basin_{grid_weights_file_prefix}{olrrp_version}.txt"
        if not scratch_basin.exists():
            create_cropping_info(
                grid_weights_generator,
                outline_shp,
                nc_path,
                scratch_basin,
                key_col=key_basin,
                verbose=False,
                generate_json=False,
            )

        out_v1 = out_dir / f"version1_{grid_weights_file_prefix}{olrrp_version}{suffix}.txt"
        replicate_subbasin_coverage_to_hrus(
            scratch_basin, hru_shp, out_v1, subbasin_col=None, hru_col=key_hru,
        )
        _prepend_header(out_v1, version_idx=1, version_desc=VERSIONS_BY_IDX[1].desc,
                        weights_shp_name=outline_shp.name, key_col=key_basin)

        scratch_sub = out_dir / f"subbasin_{grid_weights_file_prefix}{olrrp_version}{suffix}.txt"
        if not scratch_sub.exists():
            create_cropping_info(
                grid_weights_generator,
                subbasin_shp,
                nc_path,
                scratch_sub,
                key_col=key_subbasin,
                verbose=False,
                generate_json=False,
            )

        out_v2 = out_dir / f"version2_{grid_weights_file_prefix}{olrrp_version}{suffix}.txt"
        replicate_subbasin_coverage_to_hrus(
            scratch_sub, hru_shp, out_v2, subbasin_col=key_subbasin, hru_col=key_hru,
        )
        _prepend_header(out_v2, version_idx=2, version_desc=VERSIONS_BY_IDX[2].desc,
                        weights_shp_name=subbasin_shp.name, key_col=key_subbasin)

        out_v3 = out_dir / f"version3_{grid_weights_file_prefix}{olrrp_version}{suffix}.txt"
        create_cropping_info(
            grid_weights_generator,
            hru_shp,
            nc_path,
            out_v3,
            key_col=key_hru,
            verbose=False,
            generate_json=False,
        )
        _prepend_header(out_v3, version_idx=3, version_desc=VERSIONS_BY_IDX[3].desc,
                        weights_shp_name=hru_shp.name, key_col=key_hru)

        try: scratch_basin.unlink(missing_ok=True)
        except Exception: pass
        try: scratch_sub.unlink(missing_ok=True)
        except Exception: pass

def aggregate_netcdf_forcing(
    name_key: str,
    output_root: str | Path,
    shapefile_root: str | Path,
    grid_weights_generator: str | Path,
    olrrp_version: str,
    *,
    nc_path: Path | None = None,        # explicit source netcdf (recommended)
    hru_name_prefix: str = "final_hru_",
    subbasin_name_prefix: str = "finalcat_info_",
    outline_shp_name: str = "outline.shp",
    hru_version: str = "semi-distributed",
    aggregate_target: str = "hru",
    compress: bool = True,
    compression_level: int = 1,
    overwrite: bool = True,
    verbose: bool = True,
    time_block_hours: int = 24 * 7 * 8,
    output_time_chunk_hours: int = 24 * 90,
    write_gridweights: bool = True,
    # NEW controls for CF-friendly fills:
    force_numeric_fill: bool = False,
    default_numeric_fill: float = DEFAULT_NUMERIC_FILL,
) -> Path | None:
    output_root = Path(output_root)
    shapefile_root = Path(shapefile_root)

    if nc_path is not None:
        nc_path = Path(nc_path)
        if not nc_path.exists():
            if verbose:
                print(f"[agg] skip: explicit nc_path not found: {nc_path}")
            return None
    else:
        if not output_root.exists():
            if verbose:
                print(f"[agg] skip: missing output_root {output_root}")
            return None
        nc_candidates = [p for p in output_root.glob("*.nc") if not p.name.endswith("_sample.nc")]
        if not nc_candidates:
            if verbose:
                print(f"[agg] skip: no *.nc in {output_root}")
            return None
        nc_path = max(nc_candidates, key=lambda p: p.stat().st_mtime)

    suffix_letter, _ = normalize_hru_version(hru_version, olrrp_version=olrrp_version)
    tgt = (aggregate_target or "").strip().lower()
    if tgt not in {"basin", "subbasin", "hru"}:
        raise ValueError("aggregate_target must be one of 'basin','subbasin','hru'")

    if tgt == "subbasin":
        suffix_str = f"agg_subbasin_{olrrp_version}"
        shapefile_label = f"OLRRP Subbasin {olrrp_version}"
        version_num = 2
    elif tgt == "hru":
        suffix_str = f"agg_hru_{olrrp_version}{suffix_letter}"
        shapefile_label = f"OLRRP HRU {olrrp_version}{suffix_letter}"
        version_num = 3
    else:
        suffix_str = f"agg_basin_{olrrp_version}"
        shapefile_label = f"OLRRP Basin {olrrp_version}"
        version_num = 1

    subbasin_shp = shapefile_root / f"{subbasin_name_prefix}{olrrp_version}.shp"
    hru_shp      = shapefile_root / f"{hru_name_prefix}{olrrp_version}{suffix_letter}.shp"
    outline_shp  = shapefile_root / outline_shp_name
    shp_by_target = {"basin": outline_shp, "subbasin": subbasin_shp, "hru": hru_shp}
    shp = shp_by_target[tgt]
    if not shp.exists():
        if verbose:
            print(f"[agg] skip: missing shapefile {shp}")
        return None

    try:
        key_col = VERSIONS_BY_TARGET[tgt].key_col
    except Exception:
        key_col = {"basin": "ID", "subbasin": "SubId", "hru": "HRU_ID"}[tgt]

    # ------------------------------------------------------------------ #
    # Persistent intermediates
    # ------------------------------------------------------------------ #
    gw_tmp  = output_root / f"_gridweights_{suffix_str}.txt"   # keep on disk
    gi_path = output_root / f"_gridinfo_{suffix_str}.txt"      # keep on disk

    # Build only if missing; otherwise reuse
    if gi_path.exists() and gw_tmp.exists():
        if verbose:
            print(f"[agg] using cached intermediates: {gi_path.name}, {gw_tmp.name}")
    else:
        if verbose:
            print(f"[agg] building intermediates for {tgt} ({suffix_str}) ...")
        create_cropping_info(
            grid_weights_generator,
            shp,
            nc_path,
            gw_tmp,
            gi_path,
            key_col=key_col,
            verbose=False,
            generate_json=False,
        )

    unit_ids_per_entry, rows, cols, weights = _load_gridinfo_with_weights(gi_path)
    if unit_ids_per_entry.size == 0:
        if verbose:
            print("[agg] skip: empty gridinfo")
        return None

    unit_ids = np.sort(np.unique(unit_ids_per_entry))
    unit_idx_entries = np.searchsorted(unit_ids, unit_ids_per_entry)
    weights_f64 = np.asarray(weights, dtype=np.float64)

    if verbose:
        print(f"[agg] source={nc_path.name}, entries={len(rows)}, units={len(unit_ids)}, target={tgt}, suffix={suffix_str}")

    src_time_vals_num, time_units, time_calendar, time_kind = _read_time_from_source(nc_path)
    time_vals_dt = None

    use_sparse, S = False, None
    try:
        from scipy.sparse import csr_matrix  # type: ignore
        K = rows.shape[0]; U = unit_ids.shape[0]
        S = csr_matrix(
            (weights_f64, (np.arange(K, dtype=np.int64), unit_idx_entries.astype(np.int64))),
            shape=(K, U),
            dtype=np.float64,
        )
        use_sparse = True
    except Exception:
        pass

    out_vars: dict[str, dict] = {}

    with xr.open_dataset(nc_path) as ds:
        if time_kind != "numeric" and TIME_NAME in ds:
            time_vals_dt = pd.to_datetime(ds[TIME_NAME].values).to_pydatetime().tolist()
            if not time_units:
                time_units = ds[TIME_NAME].attrs.get("units") or ds[TIME_NAME].encoding.get("units")
            if not time_calendar:
                time_calendar = (
                    ds[TIME_NAME].attrs.get("calendar")
                    or ds[TIME_NAME].encoding.get("calendar")
                    or "proleptic_gregorian"
                )

        for v in ds.data_vars:
            da = ds[v]
            if TIME_NAME not in da.dims:
                continue
            try:
                yy, xx = _select_spatial_dims(da)
            except ValueError:
                continue

            da2 = da.transpose(TIME_NAME, yy, xx, missing_dims="ignore")

            Ny = da2.sizes[yy]; Nx = da2.sizes[xx]
            if (rows.max() >= Ny) or (cols.max() >= Nx):
                raise IndexError(
                    f"Grid index out of bounds for variable '{v}': "
                    f"(rows max={rows.max()} vs Ny={Ny}, cols max={cols.max()} vs Nx={Nx})."
                )

            Tn = int(da2.sizes[TIME_NAME]); Un = int(unit_ids.size)
            res = np.empty((Tn, Un), dtype=np.float64)

            block = max(1, int(time_block_hours))
            for t0 in range(0, Tn, block):
                t1 = min(Tn, t0 + block)
                da_blk = da2.isel({TIME_NAME: slice(t0, t1)})
                Vb = np.ascontiguousarray(da_blk.values[:, rows, cols]).astype(np.float64, copy=False)

                if use_sparse and S is not None:
                    Rb = Vb @ S
                    Rb = Rb.toarray() if hasattr(Rb, "toarray") else np.asarray(Rb)
                    Rb = np.asarray(Rb, dtype=np.float64, order="C")
                else:
                    Rb = _weighted_reduce_dense(Vb, unit_idx_entries, weights_f64, Un)

                res[t0:t1, :] = Rb

            var_attrs = {
                **da.attrs,
                "spatial_aggregation_from": "grid-cells",
                "spatial_aggregation_level": tgt,
                "spatial_aggregation_shapefile": shapefile_label,
                "spatial_aggregation_note": "Area Weighted aggregation",
            }
            src_enc = getattr(da, "encoding", {}) or {}
            fill_raw = src_enc.get("_FillValue", da.attrs.get("missing_value", None))
            # If upstream has NaN/invalid fill, we default to None (omit fill attrs), unless forcing numeric.
            fill_val = _sanitize_fill_value(fill_raw, coerce_to_default=False)
            out_vars[v] = {"data": res, "attrs": var_attrs, "fill": fill_val}

        stamp = datetime.utcnow().strftime("%Y-%m-%dT%H:%M UTC")
        note  = f"{stamp}  spatial aggregation to {tgt} using {shapefile_label} (key={suffix_str})"
        _append_history_line(ds, note)
        global_attrs = dict(ds.attrs)

    key_clean = (name_key or "").strip()
    out_name = auto_generate_nc_filename(
        nc_path,
        key_clean or None,
        spatial_token=suffix_str,
    )
    out_path = output_root / out_name

    if out_path.exists() and not overwrite:
        if verbose:
            print(f"[agg] exists: {out_path}")
        return out_path

    Tn = next(iter(out_vars.values()))["data"].shape[0] if out_vars else 0
    Un = int(unit_ids.size)
    with nc.Dataset(out_path, "w", format="NETCDF4") as ds_out:
        ds_out.createDimension("time", Tn)
        ds_out.createDimension("stations", Un)
        ds_out.createDimension(DIM_CHAR_ID, LEN_CHAR_ID)
        ds_out.createDimension(DIM_STRING30, LEN_STRING30)

        for k, v in global_attrs.items():
            try: setattr(ds_out, k, v)
            except Exception: setattr(ds_out, k, str(v))

        time_var = ds_out.createVariable("time", "f8", ("time",))
        if not time_units:
            raise ValueError("Source time units missing; ensure CF-compliant time.")
        time_var.units = time_units
        time_var.calendar = time_calendar or "proleptic_gregorian"
        if src_time_vals_num is not None:
            time_var[:] = src_time_vals_num
        else:
            if not time_vals_dt:
                raise ValueError("Could not obtain time values from source dataset.")
            time_var[:] = nc.date2num(time_vals_dt, units=time_var.units, calendar=time_var.calendar)

        chunk_t = min(output_time_chunk_hours, max(Tn, 1))
        chunk_s = min(512, max(Un, 1))
        for name, info in out_vars.items():
            data = info["data"]
            attrs = info["attrs"]
            src_fill_val = info["fill"]

            if force_numeric_fill:
                fill_to_use = float(default_numeric_fill)
                data = np.where(np.isnan(data), fill_to_use, data)
            else:
                fill_to_use = src_fill_val

            kwargs = {"zlib": bool(compress), "complevel": int(compression_level), "chunksizes": (chunk_t, chunk_s)}
            if fill_to_use is not None:
                var = ds_out.createVariable(name, "f8", ("time", "stations"), fill_value=fill_to_use, **kwargs)
                var.setncattr("missing_value", fill_to_use)
            else:
                var = ds_out.createVariable(name, "f8", ("time", "stations"), **kwargs)

            for k, v in attrs.items():
                try: var.setncattr(k, v)
                except Exception: var.setncattr(k, str(v))
            var[:, :] = data

        id_strings      = [str(int(u)) for u in unit_ids]
        name_strings    = [f"RAVEN-{s}" for s in id_strings]
        sid             = ds_out.createVariable("station_id", "S1", ("stations", DIM_CHAR_ID))
        sid.long_name   = "station identifier string (unit id as text)"
        sid.cf_role     = "timeseries_id"
        sid[:, :]       = _to_char_matrix(id_strings, LEN_CHAR_ID)
        
        sname               = ds_out.createVariable("station_names", "S1", ("stations", DIM_STRING30))
        sname.long_name     = "station short name"
        sname.description   = "RAVEN-<unit_id>"
        sname[:, :]         = _to_char_matrix(name_strings, LEN_STRING30)

    dataset_id, var_str = _get_dataset_id_and_varstr(out_path)
    parts = ["netcdf_metadata"]
    if key_clean:
        parts.append(key_clean)
    parts.extend([dataset_id, var_str, suffix_str])
    meta_path = out_path.parent / ("_".join(parts) + ".txt")
    export_nc_metadata_to_ascii(out_path, meta_path, decode_cf=False, show_dtypes=True)
    if verbose:
        print(f"[agg] wrote {out_path.name}, metadata {meta_path.name}")

    # Aggregated gridweights
    vN_gw_path = out_path.parent / f"version{version_num}_gridweights_{suffix_str}.txt"
    if write_gridweights:
        if vN_gw_path.exists():
            if verbose:
                print(f"[agg] gridweights present, skip writing: {vN_gw_path.name}")
        else:
            _write_version_gridweights_for_aggregated(
                out_txt=vN_gw_path,
                agg_nc_path=out_path,
                version=version_num,
                name_key=key_clean,
                olrrp_version=olrrp_version,
                hru_name_prefix=hru_name_prefix,
                suffix_letter=suffix_letter,
                shapefile_root=shapefile_root,
            )
            if verbose:
                print(f"[agg] wrote {vN_gw_path.name}")
    else:
        if verbose:
            print(f"[agg] write_gridweights=False (skipping weights file for {suffix_str})")

    return out_path


def extract_gauge_from_regional_agg_nc(
    agg_nc_path: Path,
    unit_ids_to_keep: List[int],
    out_nc_path: Path,
    *,
    station_id_var: str = "station_id",
    station_names_var: str = "station_names",
    time_name: str = "time",
    stations_dim: str = "stations",
    char_id_dim: str = "char_leng_id",
    string30_dim: str = "string30",
    time_chunk_size: int = 24*365,
    stations_chunk_cap: int = 512,
    compress: bool = True,
    compression_level: int = 1,
    force_numeric_fill: bool = False,
    default_numeric_fill: float = DEFAULT_NUMERIC_FILL,
) -> Path:
    out_nc_path.parent.mkdir(parents=True, exist_ok=True)

    src = nc.Dataset(str(agg_nc_path), "r")
    try:
        if station_id_var not in src.variables:
            raise KeyError(f"'{station_id_var}' not found in {agg_nc_path.name}")

        sid_var = src.variables[station_id_var]
        S = sid_var.shape[0]

        src_ids = []
        for i in range(S):
            raw = sid_var[i, :].tobytes().decode("ascii", "ignore").strip()
            digits = "".join([c if c.isdigit() else " " for c in raw]).split()
            src_ids.append(int(digits[-1]) if digits else None)

        wanted = set(int(x) for x in unit_ids_to_keep)
        keep_idx = [i for i, s in enumerate(src_ids) if (s is not None) and (s in wanted)]
        if not keep_idx:
            raise RuntimeError(
                f"No matching units found in {agg_nc_path.name} for requested {len(unit_ids_to_keep)} IDs."
            )

        dst = nc.Dataset(str(out_nc_path), "w", format="NETCDF4")
        try:
            for k, v in src.__dict__.items():
                try: setattr(dst, k, v)
                except Exception: setattr(dst, k, str(v))

            for dim_name, dim in src.dimensions.items():
                if dim_name in (time_name, stations_dim, char_id_dim, string30_dim):
                    continue
                dst.createDimension(dim_name, (len(dim) if not dim.isunlimited() else None))

            T = len(src.dimensions[time_name]) if time_name in src.dimensions else 0
            dst.createDimension(time_name, T)
            dst.createDimension(stations_dim, len(keep_idx))

            if char_id_dim in src.dimensions:
                dst.createDimension(char_id_dim, len(src.dimensions[char_id_dim]))
            if string30_dim in src.dimensions:
                dst.createDimension(string30_dim, len(src.dimensions[string30_dim]))

            if time_name in src.variables:
                t_src = src.variables[time_name]
                t_dst = dst.createVariable(time_name, t_src.dtype, (time_name,))
                for k, v in t_src.__dict__.items():
                    try: setattr(t_dst, k, v)
                    except Exception: setattr(t_dst, k, str(v))
                t_dst[:] = t_src[:]

            sid_dst = dst.createVariable(station_id_var, "S1", (stations_dim, char_id_dim))
            for k, v in sid_var.__dict__.items():
                if k == "_FillValue":
                    continue
                try: setattr(sid_dst, k, v)
                except Exception: setattr(sid_dst, k, str(v))
            for j, i in enumerate(keep_idx):
                sid_dst[j, :] = sid_var[i, :]

            if station_names_var in src.variables and string30_dim in src.dimensions:
                sname_src = src.variables[station_names_var]
                sname_dst = dst.createVariable(station_names_var, "S1", (stations_dim, string30_dim))
                for k, v in sname_src.__dict__.items():
                    if k == "_FillValue":
                        continue
                    try: setattr(sname_dst, k, v)
                    except Exception: setattr(sname_dst, k, str(v))
                for j, i in enumerate(keep_idx):
                    sname_dst[j, :] = sname_src[i, :]

            def is_ts_var(v):
                dims = v.dimensions
                return (time_name in dims) and (stations_dim in dims) and (len(dims) == 2)

            st_chunk = min(stations_chunk_cap, len(keep_idx))
            for name, v in src.variables.items():
                if name in (time_name, station_id_var, station_names_var):
                    continue
                if is_ts_var(v):
                    kwargs = {}
                    if compress:
                        kwargs.update({"zlib": True, "complevel": int(compression_level)})
                    try:
                        kwargs.update({"chunksizes": (min(time_chunk_size, max(T, 1)), st_chunk)})
                    except Exception:
                        pass

                    fill_raw = getattr(v, "_FillValue", None)
                    if force_numeric_fill:
                        fill_to_use = float(default_numeric_fill)
                    else:
                        fill_to_use = _sanitize_fill_value(fill_raw, coerce_to_default=False)

                    if fill_to_use is not None:
                        var_dst = dst.createVariable(name, v.dtype, (time_name, stations_dim), fill_value=fill_to_use, **kwargs)
                        var_dst.setncattr("missing_value", fill_to_use)
                    else:
                        var_dst = dst.createVariable(name, v.dtype, (time_name, stations_dim), **kwargs)

                    for k, val in v.__dict__.items():
                        if k == "_FillValue":
                            continue
                        try: setattr(var_dst, k, val)
                        except Exception: setattr(var_dst, k, str(val))

            for t0 in range(0, T, time_chunk_size):
                t1 = min(T, t0 + time_chunk_size)
                for name, v in src.variables.items():
                    if name in (time_name, station_id_var, station_names_var):
                        continue
                    if is_ts_var(v):
                        block = v[t0:t1, keep_idx]
                        if force_numeric_fill:
                            # replace NaNs when forcing numeric
                            block = np.where(np.isnan(block), float(default_numeric_fill), block)
                        dst.variables[name][t0:t1, :] = block

            stamp = datetime.utcnow().strftime("%Y-%m-%dT%H:%M UTC")
            hist = getattr(dst, "history", "")
            note = f"{stamp}  subset {len(keep_idx)} stations from regional aggregated file"
            setattr(dst, "history", (hist + ("\n" if hist else "") + note))
        finally:
            dst.close()
    finally:
        src.close()

    return out_nc_path


def crop_netcdf_by_cell_info(
    nc_file_path: str | Path,
    grid_info_path: str | Path,
    output_nc_path: str | Path,
    crop_buffer_cell: int = 3,
    use_nco: bool = True,
    compress: bool = True,
    compression_level: int = 4,
    *,
    time_name: str = "time",
    rlat_name: str = "rlat",
    rlon_name: str = "rlon",
    copy_chunk_hours: int = 24*365,
) -> Path:
    nc_file   = Path(nc_file_path)
    grid_info = Path(grid_info_path)
    output_nc = Path(output_nc_path)

    if not grid_info.is_file():
        return output_nc

    df = pd.read_csv(grid_info, sep=r"\s+", engine="python", usecols=["Row", "Col"])
    r0 = int(df["Row"].min() - crop_buffer_cell)
    r1 = int(df["Row"].max() + crop_buffer_cell)
    c0 = int(df["Col"].min() - crop_buffer_cell)
    c1 = int(df["Col"].max() + crop_buffer_cell)

    ts = datetime.utcnow().strftime("%Y-%m-%dT%H:%M UTC")

    if use_nco:
        # Clamp indices to valid ranges for ncks
        with nc.Dataset(str(nc_file), "r") as _src:
            Ny = len(_src.dimensions.get(rlat_name, []))
            Nx = len(_src.dimensions.get(rlon_name, []))
        r0c = max(0, min(Ny - 1, r0)); r1c = max(0, min(Ny - 1, r1))
        if r1c < r0c: r0c, r1c = r1c, r0c
        c0c = max(0, min(Nx - 1, c0)); c1c = max(0, min(Nx - 1, c1))
        if c1c < c0c: c0c, c1c = c1c, c0c

        cmd = [
            "ncks", "-h",
            f"-L{compression_level if compress else 0}",
            "-d", f"{rlat_name},{r0c},{r1c}",
            "-d", f"{rlon_name},{c0c},{c1c}",
            str(nc_file),
            str(output_nc),
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        with nc.Dataset(str(output_nc), "a") as ds:
            old = ds.getncattr("history") if "history" in ds.ncattrs() else ""
            new_hist = f"{ts}  cropped extent {rlat_name}[{r0c}:{r1c}], {rlon_name}[{c0c}:{c1c}] (buffer={crop_buffer_cell})"
            ds.history = (old + "\n" + new_hist) if old else new_hist
        return output_nc

    with xr.open_dataset(nc_file, engine="netcdf4") as ds:
        cand_dims = [
            (rlat_name, rlon_name),
            ("rlat", "rlon"),
            ("y", "x"),
            ("lat", "lon"),
        ]
        y_dim = x_dim = None
        for yy, xx in cand_dims:
            if yy in ds.dims and xx in ds.dims:
                y_dim, x_dim = yy, xx
                break
        if y_dim is None or x_dim is None:
            raise ValueError(f"Could not find 2D spatial dims in {list(ds.dims)}; tried {cand_dims}")

        Ny = ds.sizes[y_dim]; Nx = ds.sizes[x_dim]
        r0c = max(0, min(Ny-1, r0)); r1c = max(0, min(Ny-1, r1))
        if r1c < r0c: r0c, r1c = r1c, r0c
        c0c = max(0, min(Nx-1, c0)); c1c = max(0, min(Nx-1, c1))
        if c1c < c0c: c0c, c1c = c1c, c0c

        sel = {y_dim: slice(r0c, r1c+1), x_dim: slice(c0c, c1c+1)}
        ds_c = ds.isel(**sel)

        for coord in (time_name, "time_bnds", "lon", "lat", y_dim, x_dim):
            if coord in ds_c.coords:
                ds_c.coords[coord].encoding.pop("_FillValue", None)
                ds_c.coords[coord].attrs.pop("_FillValue", None)

        old_hist = ds.attrs.get("history", "")
        new_hist = f"{ts}  cropped extent {y_dim}[{r0c}:{r1c}], {x_dim}[{c0c}:{c1c}] (buffer={crop_buffer_cell})"
        ds_c.attrs["history"] = (old_hist + ("\n" if old_hist else "") + new_hist)

        encoding = {}
        T = int(ds.sizes.get(time_name, 0)) if time_name in ds.sizes else 0
        t_chunk = min(copy_chunk_hours, max(T, 1)) if T else 1
        for v in ds_c.data_vars:
            da = ds_c[v]
            enc = {}
            if compress:
                enc.update(dict(zlib=True, complevel=int(compression_level)))
            chunks = []
            for d in da.dims:
                if d == time_name and T:
                    chunks.append(t_chunk)
                elif d == y_dim:
                    chunks.append(int(r1c - r0c + 1))
                elif d == x_dim:
                    chunks.append(int(c1c - c0c + 1))
                else:
                    chunks.append(int(da.sizes[d]))
            if chunks:
                enc["chunksizes"] = tuple(chunks)
            encoding[v] = enc

        ds_c.to_netcdf(output_nc, engine="netcdf4", encoding=encoding if compress else None)

    return output_nc


def crop_netcdf_time_range(
    nc_in: str | Path,
    nc_out: str | Path,
    t_start,
    t_end,
    compress: bool = True,
    compression_level: int = 4,
    *,
    time_name: str = "time",
    copy_chunk_hours: int = 24*365,
) -> Path:
    """Copy an inclusive CF time range without changing stored data or packing."""
    if Path(nc_in).resolve() == Path(nc_out).resolve():
        raise ValueError("Time-crop input and output must be different files")
    if copy_chunk_hours < 1:
        raise ValueError("copy_chunk_hours must be positive")

    def _to_ts(value):
        timestamp = pd.Timestamp(*value) if isinstance(value, (tuple, list)) else pd.Timestamp(value)
        if pd.isna(timestamp):
            raise ValueError("Time-crop bounds must be valid dates")
        if timestamp.tzinfo is not None:
            timestamp = timestamp.tz_convert("UTC").tz_localize(None)
        return timestamp

    t0, t1 = _to_ts(t_start), _to_ts(t_end)
    if t0 > t1:
        raise ValueError("Time-crop start must not be after its end")

    with nc.Dataset(str(nc_in), "r") as src:
        if time_name not in src.variables:
            raise KeyError(f"'{time_name}' coordinate missing in {nc_in}")
        tvar = src.variables[time_name]
        if tvar.dimensions != (time_name,) or not getattr(tvar, "units", None):
            raise ValueError(f"'{time_name}' must be a one-dimensional CF time coordinate with units")
        times = tvar[:]
        if np.ma.getmaskarray(times).any() or not np.isfinite(times).all():
            raise ValueError("Time coordinate contains missing or non-finite values")
        times = np.asarray(times)
        if np.any(np.diff(times) <= 0):
            raise ValueError("Time coordinate must be strictly increasing")
        bounds = nc.date2num(
            [t0.to_pydatetime(), t1.to_pydatetime()],
            units=tvar.units,
            calendar=getattr(tvar, "calendar", "standard"),
        )
        idx = np.flatnonzero((times >= bounds[0]) & (times <= bounds[1]))
        if not len(idx):
            raise ValueError(f"No times in range {t0}..{t1}")
        i0, i1 = int(idx[0]), int(idx[-1]) + 1

        # Copy raw values so scale_factor, add_offset, and fill codes remain exact.
        src.set_auto_maskandscale(False)
        with nc.Dataset(str(nc_out), "w", format="NETCDF4") as dst:
            dst.setncatts(src.__dict__)
            for dname, dim in src.dimensions.items():
                length = i1 - i0 if dname == time_name else (None if dim.isunlimited() else len(dim))
                dst.createDimension(dname, length)

            for vname, variable in src.variables.items():
                kwargs = {}
                if "_FillValue" in variable.ncattrs():
                    kwargs["fill_value"] = variable.getncattr("_FillValue")
                if compress and variable.ndim and getattr(variable.dtype, "kind", None) in ("f", "i", "u"):
                    kwargs.update(zlib=True, complevel=int(compression_level))
                if time_name in variable.dimensions:
                    kwargs["chunksizes"] = tuple(
                        min(copy_chunk_hours, i1 - i0) if dim == time_name else min(max(size, 1), 512)
                        for dim, size in zip(variable.dimensions, variable.shape)
                    )
                output = dst.createVariable(vname, variable.dtype, variable.dimensions, **kwargs)
                output.setncatts({key: value for key, value in variable.__dict__.items() if key != "_FillValue"})
                output.set_auto_maskandscale(False)
                if time_name in variable.dimensions:
                    tpos = variable.dimensions.index(time_name)
                    for start in range(i0, i1, copy_chunk_hours):
                        end = min(i1, start + copy_chunk_hours)
                        source_slice = [slice(None)] * variable.ndim
                        source_slice[tpos] = slice(start, end)
                        output_slice = list(source_slice)
                        output_slice[tpos] = slice(start - i0, end - i0)
                        output[tuple(output_slice)] = variable[tuple(source_slice)]
                else:
                    output[...] = variable[...]

            stamp = datetime.utcnow().strftime("%Y-%m-%dT%H:%M UTC")
            history = getattr(dst, "history", "")
            note = f"{stamp}  time-cropped to {t0} - {t1} (inclusive)"
            dst.history = (history + "\n" + note) if history else note

    return Path(nc_out)


def crop_gridded_nc_for_gauge(
    copy_region_forcing_nc_path: Path,
    gauge_outline_shp: Path,
    gauge_output_folder: Path,
    gauge_ID: str,
    *,
    grid_weights_generator: Path,
    crop_buffer_cell: int = 3,
    use_nco: bool = False,
    compress: bool = False,
    compression_level: int = 1,
    sample: bool = False,
    sample_time_start: Optional[Tuple[int, int, int, int, int, int]] = None,
    sample_time_end: Optional[Tuple[int, int, int, int, int, int]] = None,
    sample_suffix: str = "_sample",
    common_header: Optional[str] = None,  # not used here; kept for API symmetry
) -> Path:
    """
    Crop regional gridded NetCDF to the gauge outline extent and export metadata.
    Reuses a tiny shared crop grid-info cache (per gauge & grid domain) across variables/runs.
    Optionally writes a smaller time-sampled file.
    """
    gauge_output_folder.mkdir(parents=True, exist_ok=True)

    # Per-run cache (legacy location used by cropper)
    gridinfo = gauge_output_folder / "gridweights_for_crop_cached.txt"

    # Shared cache across variables/runs for this gauge & domain
    shared_cache_dir = gauge_output_folder.parent / ".crop_cache"
    shared_cache_dir.mkdir(parents=True, exist_ok=True)
    shared_cache = shared_cache_dir / f"{gauge_ID}.txt"

    # Try reusing shared cache first
    if not gridinfo.exists() and shared_cache.exists():
        try:
            try:
                if gridinfo.exists():
                    gridinfo.unlink()
                os.link(shared_cache, gridinfo)
            except Exception:
                try:
                    os.symlink(shared_cache, gridinfo)
                except Exception:
                    shutil.copy2(shared_cache, gridinfo)
        except Exception:
            pass

    if not gridinfo.exists():
        temp_w = gauge_output_folder / "gridweights_temp.txt"
        create_cropping_info(
            grid_weights_generator,
            gauge_outline_shp,
            copy_region_forcing_nc_path,
            temp_w,
            gridinfo,
            key_col=_lazy_key_cols()["basin"],
            verbose=False,
            generate_json=False,
        )
        temp_w.unlink(missing_ok=True)

        try:
            if shared_cache.exists():
                shared_cache.unlink()
            try:
                os.link(gridinfo, shared_cache)
            except Exception:
                try:
                    os.symlink(gridinfo, shared_cache)
                except Exception:
                    shutil.copy2(gridinfo, shared_cache)
        except Exception:
            pass

    gauge_forcing_nc_name = auto_generate_nc_filename(copy_region_forcing_nc_path, gauge_ID)
    gauge_forcing_nc_path = gauge_output_folder / gauge_forcing_nc_name
    crop_netcdf_by_cell_info(
        nc_file_path=copy_region_forcing_nc_path,
        grid_info_path=gridinfo,
        output_nc_path=gauge_forcing_nc_path,
        crop_buffer_cell=crop_buffer_cell,
        use_nco=use_nco,
        compress=compress,
        compression_level=compression_level,
    )

    if sample and sample_time_start and sample_time_end:
        tmp_sample = gauge_output_folder / "_tmp_sample.nc"
        crop_netcdf_time_range(
            nc_in=gauge_forcing_nc_path,
            nc_out=tmp_sample,
            t_start=sample_time_start,
            t_end=sample_time_end,
            compress=compress,
            compression_level=compression_level,
        )
        core_name = auto_generate_nc_filename(tmp_sample, gauge_ID)
        sample_name = core_name.replace(".nc", f"{sample_suffix}.nc")
        (gauge_output_folder / sample_name).write_bytes(tmp_sample.read_bytes())
        tmp_sample.unlink(missing_ok=True)

    export_nc_metadata_to_ascii(
        gauge_forcing_nc_path,
        gauge_output_folder / f"netcdf_metadata_{gauge_ID}_{__import__('forcing_utils').get_variable_str(gauge_forcing_nc_path)}.txt",
    )
    return gauge_forcing_nc_path
