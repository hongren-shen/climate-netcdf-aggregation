#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2019-2026 Hongren Shen <hongren.shen@uwaterloo.ca>
# SPDX-License-Identifier: LGPL-3.0-or-later

"""
Temporal aggregation utilities for gridded or station-style NetCDF forcings.

Key features
------------
- Handles gridded (time, y, x) and station (time, stations) layouts.
- Optional timezone shift BEFORE aggregation (compatiable with "America/Toronto" or fixed "UTC-5").
- Optional conversion from period-ending to period-starting timestamps (shift back by one native step).
- Aggregation methods: sum, mean, min, max, first, last.
- Robust chunked processing to keep memory bounded.
- Preserves station_id/station_names when present.
- Adds clear metadata about temporal processing.
- Daily-like grouping ("D","1D","24H") uses the **first anchor present in data** and **labels on the left edge**,
  allowing the first aggregated day to be partial (e.g., 1979-12-31 07:00-23:00 local), labeled at the first observed
  midnight (e.g., 1979-12-31 00:00).
"""

from __future__ import annotations

import re
from pathlib import Path
from datetime import datetime, timedelta
from typing import List, Tuple, Optional, Dict, Any

import numpy as np
import pandas as pd
import netCDF4 as nc

try:
    from dateutil import tz as _dtz
    _HAVE_DATEUTIL_TZ = True
except Exception:
    _HAVE_DATEUTIL_TZ = False


def _normalize_freq_for_display(freq: str) -> str:
    F = (freq or "").upper()
    return "24H" if F in ("D", "1D") else F

def _strip_source_temporal_hint(stem: str) -> str:
    s = re.sub(r'(?i)(?:^|_)(?:\d{1,2}hourly|hourly|daily|monthly)(?=_|$)', '', stem)
    s = re.sub(r'__+', '_', s).strip('_')
    return s

def _name_token(freq_display: str, tz_label: str, anchor_hour: int) -> str:
    tz = (tz_label or "UTC").replace("/", "-")
    tok = f"{freq_display}_local-{tz}"
    if int(anchor_hour) != 0:
        tok += f"_H{int(anchor_hour):02d}"
    return tok

def _decode_cf_time(tvar) -> Tuple[np.ndarray, str, str]:
    units = getattr(tvar, "units", None)
    cal   = getattr(tvar, "calendar", "proleptic_gregorian")
    if units is None:
        raise ValueError("Missing CF-compliant time 'units' attribute.")
    vals = np.asarray(tvar[:], dtype=np.float64)
    return vals, units, cal

def _cfnum_to_datetimes(vals: np.ndarray, units: str, cal: str) -> List[datetime]:
    try:
        dts = nc.num2date(vals, units=units, calendar=cal, only_use_cftime_datetimes=False)
    except TypeError:
        dts = nc.num2date(vals, units=units, calendar=cal)
    out: List[datetime] = []
    for dt in dts:
        if isinstance(dt, datetime):
            out.append(dt)
        else:
            sec_float = float(getattr(dt, "second", 0.0))
            sec_int = int(sec_float)
            usec = int(round((sec_float - sec_int) * 1e6))
            out.append(datetime(
                int(dt.year), int(dt.month), int(dt.day),
                int(dt.hour), int(dt.minute), sec_int, usec
            ))
    return out

def _apply_tz_shift(dts: List[datetime],
                    tz_name: Optional[str],
                    tz_offset_hours: Optional[float]) -> Tuple[List[datetime], str]:
    if not dts:
        return dts, "UTC"
    if tz_name and _HAVE_DATEUTIL_TZ:
        utc = _dtz.UTC
        loc = _dtz.gettz(tz_name)
        if loc is not None:
            conv = [dt.replace(tzinfo=utc).astimezone(loc).replace(tzinfo=None) for dt in dts]
            return conv, tz_name
    if tz_offset_hours is not None:
        total_minutes = int(round(float(tz_offset_hours) * 60.0))
        return [dt + timedelta(minutes=total_minutes) for dt in dts], (f"UTC{tz_offset_hours:+g}" if tz_offset_hours else "UTC")
    return dts, "UTC"

def _infer_native_step_hours(dts: List[datetime]) -> int:
    if len(dts) < 2:
        return 1
    diffs = np.diff(np.array(dts, dtype="datetime64[ns]")).astype("timedelta64[h]").astype(int)
    diffs = diffs[diffs > 0]
    if diffs.size == 0:
        return 1
    med = int(np.median(diffs))
    return max(1, med)

def _group_index(dts_local: list, *, freq: str, anchor_hour: int,
                 prefer_first_anchor_from_data: bool = False, label: str = "left"):
    if not dts_local:
        return np.zeros((0,), dtype=np.int64), []

    s = pd.DatetimeIndex(dts_local).sort_values()
    f = (freq or "").upper()
    day_like = f in ("D", "1D", "24H")

    if day_like and prefer_first_anchor_from_data:
        anchors = s[(s.hour == int(anchor_hour)) & (s.minute == 0) & (s.second == 0)]
        if len(anchors) > 0:
            t0 = anchors[0]
            delta_days = (s - t0) / np.timedelta64(1, "D")
            if (label or "").lower() == "right":
                k = np.floor(delta_days).astype(int) + 1
            else:
                k = np.floor(delta_days).astype(int)
            labels = t0 + pd.to_timedelta(k, unit="D")
            uniq, inv = np.unique(labels.values, return_inverse=True)
            group_idx = inv.astype(np.int64, copy=False)
            group_labels = pd.to_datetime(uniq).to_pydatetime().tolist()
            return group_idx, group_labels

    if day_like:
        shift = pd.to_timedelta(anchor_hour, unit="h")
        floors = (s - shift).floor("D") + shift
        if (label or "").lower() == "right":
            floors = floors + pd.to_timedelta(1, "D")
    else:
        try:
            floors = s.floor(f)
        except Exception:
            shift = pd.to_timedelta(anchor_hour, unit="h")
            floors = (s - shift).floor("D") + shift

    uniq, inv = np.unique(floors.values, return_inverse=True)
    return inv.astype(np.int64, copy=False), pd.to_datetime(uniq).to_pydatetime().tolist()

def _time_axis_pos(v, time_name: str) -> int:
    try:
        return v.dimensions.index(time_name)
    except ValueError:
        return -1

def _alloc_reducer(shape: Tuple[int, ...], method: str) -> Dict[str, np.ndarray]:
    m = method.lower()
    if m in ("sum", "mean"):
        return {"sum": np.zeros(shape, dtype=np.float64),
                "cnt": np.zeros(shape, dtype=np.int64)}
    if m == "max":
        return {"val": np.full(shape, -np.inf, dtype=np.float64)}
    if m == "min":
        return {"val": np.full(shape, +np.inf, dtype=np.float64)}
    if m == "first":
        return {"val": np.full(shape, np.nan, dtype=np.float64),
                "done": np.zeros(shape, dtype=bool)}
    if m == "last":
        return {"val": np.full(shape, np.nan, dtype=np.float64)}
    raise ValueError(f"Unknown temporal_reduce: {method}")

def _reduce_apply(R: Dict[str, np.ndarray], chunk: np.ndarray, method: str):
    m = method.lower()
    if m in ("sum", "mean"):
        R["sum"] += np.sum(chunk, axis=0, dtype=np.float64)
        R["cnt"] += chunk.shape[0]
    elif m == "max":
        R["val"] = np.maximum(R["val"], np.max(chunk, axis=0))
    elif m == "min":
        R["val"] = np.minimum(R["val"], np.min(chunk, axis=0))
    elif m == "first":
        sel = ~R["done"]
        if sel.any():
            R["val"][sel] = chunk[0].reshape(R["val"].shape)[sel]
            R["done"][sel] = True
    elif m == "last":
        R["val"] = chunk[-1].reshape(R["val"].shape)
    else:
        raise ValueError(m)

def _reduce_finalize(R: Dict[str, np.ndarray], method: str) -> np.ndarray:
    m = method.lower()
    if m == "sum":
        return R["sum"]
    if m == "mean":
        return R["sum"] / np.maximum(R["cnt"], 1)
    return R["val"]

# -------- CF cell_methods rewrite (preserve non-time clauses, replace time: …) -------- #
def _rewrite_cell_methods(prev: str, *, time_method: str, freq_display: str,
                          anchor_hour: int, tz_label: str, label_side: str,
                          period_start: bool | None) -> str:
    """
    Keep any non-time cell_methods (e.g., 'area: mean'), replace only the 'time: ...' clause.
    Append a detailed comment.
    """
    # Strip existing time: ... clause(s)
    def _strip_time(s: str) -> str:
        # Remove 'time: ...' segments up to ';' or end
        return re.sub(r'(?is)(?:^|\s)time\s*:[^;]*;?', ' ', s).strip()

    kept = _strip_time(prev or "")
    kept = kept.strip()
    if kept and not kept.endswith(";"):
        kept += ";"

    period_txt = "start" if (period_start is True) else ("end" if (period_start is False) else "as_is")
    cm_time = (f"time: {time_method.lower()} "
               f"(comment: resampled freq={freq_display}; anchor={anchor_hour}; tz={tz_label}; "
               f"label={label_side}; period={period_txt})")
    return (kept + (" " if kept else "") + cm_time).strip()

def aggregate_netcdf_temporal(
    *,
    src_nc_path: Path,
    out_nc_path: Path | None = None,
    temporal_reduce: str,                 # "sum","mean","max","min","first","last"
    temporal_freq: str,                   # "D","24H","6H","W-MON","M",...
    temporal_anchor_hour: int = 0,        # for daily-like windows
    tz_name: str | None = None,           # timezone name (preferred)
    tz_offset_hours: float | None = None, # or fixed offset
    period_is_ending: bool = False,       # if True - shift timestamps backward by source step
    period_step_hours: int | None = None, # if None, auto-infer from time deltas
    # common names/layout
    time_name: str = "time",
    stations_dim: str = "stations",
    char_id_dim: str = "char_leng_id",
    string30_dim: str = "string30",
    # IO & perf
    compress: bool = True,
    compression_level: int = 1,
    output_time_chunk_hours: int = 24 * 365,
    overwrite: bool = True,
    verbose: bool = True,
) -> Path:
    """
    Temporally aggregate a forcing NetCDF (gridded or spatially aggregated) to a coarser
    resolution (e.g., hourly to daily, hourly to monthly), with optional timezone shifting
    and handling of period-ending timestamps.

    Daily-like rules ("D","1D","24H") use:
      - first anchor present in data (e.g., first observed 00:00 local)
      - LEFT-edge labels (period_start semantics after period-ending correction)
    """
    src_nc_path = Path(src_nc_path)
    if not src_nc_path.exists():
        raise FileNotFoundError(src_nc_path)

    src = nc.Dataset(str(src_nc_path), "r")
    try:
        if time_name not in src.variables:
            raise KeyError(f"'{time_name}' variable missing in {src_nc_path.name}")
        tvar = src.variables[time_name]
        t_vals_num, t_units, t_cal = _decode_cf_time(tvar)

        dts_utc = _cfnum_to_datetimes(t_vals_num, t_units, t_cal)
        dts_local, tz_label = _apply_tz_shift(dts_utc, tz_name, tz_offset_hours)

        inferred_step_h = int(period_step_hours) if period_step_hours is not None else _infer_native_step_hours(dts_local)

        # optional period-ending -> starting shift
        pe2ps_applied = False
        if period_is_ending:
            dts_local = [dt - timedelta(hours=inferred_step_h) for dt in dts_local]
            pe2ps_applied = True

        day_like = (temporal_freq or "").upper() in ("D", "1D", "24H")
        label_side = "left"  # left-edge labeling for output coords

        group_idx, group_labels = _group_index(
            dts_local,
            freq=temporal_freq,
            anchor_hour=temporal_anchor_hour,
            prefer_first_anchor_from_data=day_like,
            label=label_side,
        )
        G = len(group_labels)
        if G == 0:
            raise ValueError("Temporal grouping produced zero groups.")

        if out_nc_path is None:
            freq_display = _normalize_freq_for_display(temporal_freq)
            stem_clean   = _strip_source_temporal_hint(src_nc_path.stem)
            token        = _name_token(freq_display, tz_label, temporal_anchor_hour)
            out_nc_path  = src_nc_path.with_name(f"{stem_clean}_{token}.nc")

        out_nc_path = Path(out_nc_path)
        if out_nc_path.exists() and not overwrite:
            if verbose:
                print(f"[t-agg] exists: {out_nc_path.name}")
            return out_nc_path

        dst = nc.Dataset(str(out_nc_path), "w", format="NETCDF4")
        try:

            for k, v in src.__dict__.items():
                try: setattr(dst, k, v)
                except Exception: setattr(dst, k, str(v))

            for dname, dim in src.dimensions.items():
                if dname == time_name:
                    continue
                dst.createDimension(dname, (len(dim) if not dim.isunlimited() else None))
            dst.createDimension(time_name, G)

            tv_out = dst.createVariable(time_name, "f8", (time_name,))
            tv_out.calendar = t_cal

            day_like = (temporal_freq or "").upper() in ("D", "1D", "24H")
            origin_dt = group_labels[0]

            if day_like:
                new_units = f"days since {origin_dt.strftime('%Y-%m-%dT%H:%M:%S')}"
                tv_out.units = new_units
                tv_out[:] = np.arange(len(group_labels), dtype=np.float64)
            else:
                new_units = f"hours since {origin_dt.strftime('%Y-%m-%dT%H:%M:%S')}"
                tv_out.units = new_units
                tv_out[:] = nc.date2num(group_labels, units=new_units, calendar=t_cal)

            def _ensure_var(_dst, name, dtype, dims):
                if name in _dst.variables:
                    return _dst.variables[name]
                return _dst.createVariable(name, dtype, dims)

            if stations_dim in src.dimensions and "station_id" in src.variables and char_id_dim in src.dimensions:
                sid_src = src.variables["station_id"]
                sid_dst = _ensure_var(dst, "station_id", "S1", (stations_dim, char_id_dim))
                sid_dst[:, :] = sid_src[:, :]
                for k, v in sid_src.__dict__.items():
                    if k == "_FillValue": continue
                    try: setattr(sid_dst, k, v)
                    except Exception: setattr(sid_dst, k, str(v))

            if stations_dim in src.dimensions and "station_names" in src.variables and string30_dim in src.dimensions:
                sname_src = src.variables["station_names"]
                sname_dst = _ensure_var(dst, "station_names", "S1", (stations_dim, string30_dim))
                sname_dst[:, :] = sname_src[:, :]
                for k, v in sname_src.__dict__.items():
                    if k == "_FillValue": continue
                    try: setattr(sname_dst, k, v)
                    except Exception: setattr(sname_dst, k, str(v))

            reducers: Dict[str, Dict] = {}

            def _dst_chunks(vsrc) -> Tuple | None:
                dims = vsrc.dimensions
                if time_name not in dims:
                    return None
                tchunk = min(max(1, output_time_chunk_hours), G)
                chunks = []
                for d in dims:
                    if d == time_name:
                        chunks.append(tchunk)
                    else:
                        length = vsrc.shape[dims.index(d)]
                        cap = 1024 if d not in (stations_dim,) else 512
                        chunks.append(min(length, cap))
                return tuple(chunks)

            freq_display = _normalize_freq_for_display(temporal_freq)

            for vname, vsrc in src.variables.items():
                if vname in dst.variables or vname == time_name:
                    continue

                if time_name not in vsrc.dimensions:
                    vdst = dst.createVariable(vname, vsrc.dtype, vsrc.dimensions)
                    vdst[:] = vsrc[:]
                    for k, v in vsrc.__dict__.items():
                        if k == "_FillValue": continue
                        try: setattr(vdst, k, v)
                        except Exception: setattr(vdst, k, str(v))
                    continue

                kwargs: Dict[str, Any] = {}
                if compress and vsrc.dtype.kind in ("f", "i", "u"):
                    ch = _dst_chunks(vsrc)
                    if ch is not None:
                        kwargs.update(dict(zlib=True, complevel=int(compression_level), chunksizes=ch))
                    else:
                        kwargs.update(dict(zlib=True, complevel=int(compression_level)))

                vdst = dst.createVariable(vname, vsrc.dtype, vsrc.dimensions, **kwargs)
                for k, v in vsrc.__dict__.items():
                    if k == "_FillValue": continue
                    try: setattr(vdst, k, v)
                    except Exception: setattr(vdst, k, str(v))

                try: vdst.temporal_source_step_hours = int(inferred_step_h)
                except Exception: pass

                prev_cm = getattr(vsrc, "cell_methods", "") or ""
                vdst.cell_methods = _rewrite_cell_methods(
                    prev_cm,
                    time_method=temporal_reduce,
                    freq_display=freq_display,
                    anchor_hour=int(temporal_anchor_hour),
                    tz_label=tz_label,
                    label_side=label_side,
                    period_start=(pe2ps_applied or day_like)
                )

                try: vdst.temporal_aggregation_method = temporal_reduce.lower()
                except Exception: pass
                try: vdst.temporal_aggregation_freq = freq_display
                except Exception: pass
                try: vdst.temporal_aggregation_anchor_hour = int(temporal_anchor_hour)
                except Exception: pass
                try: vdst.time_period_convention = ("period_start" if (pe2ps_applied or day_like) else "as_is")
                except Exception: pass
                try: vdst.time_timezone = tz_label
                except Exception: pass

                tpos = _time_axis_pos(vsrc, time_name)
                out_shape = tuple(vsrc.shape[i] for i in range(len(vsrc.shape)) if i != tpos)
                reducers[vname] = {
                    "dst": vdst,
                    "tpos": tpos,
                    "groups": [_alloc_reducer(out_shape, temporal_reduce) for _ in range(G)],
                }

            T = len(src.dimensions[time_name])
            step = max(1, output_time_chunk_hours)
            for t0 in range(0, T, step):
                t1 = min(T, t0 + step)
                gids = group_idx[t0:t1]
                if gids.size == 0:
                    continue

                for vname, slot in reducers.items():
                    vsrc = src.variables[vname]
                    tpos = slot["tpos"]
                    base_slices = [slice(None)] * vsrc.ndim
                    base_slices[tpos] = slice(t0, t1)
                    slab = vsrc[tuple(base_slices)]
                    if tpos != 0:
                        slab = np.moveaxis(slab, tpos, 0)

                    for g in np.unique(gids):
                        mask = (gids == g)
                        if not np.any(mask):
                            continue
                        chunk = slab[mask, ...]
                        _reduce_apply(slot["groups"][int(g)], chunk, temporal_reduce)

            for vname, slot in reducers.items():
                vdst = slot["dst"]
                tpos = slot["tpos"]
                stacked = [_reduce_finalize(slot["groups"][g], temporal_reduce) for g in range(G)]
                data = np.stack(stacked, axis=0)
                if tpos != 0:
                    data = np.moveaxis(data, 0, tpos)
                vdst[:] = data

            # append a history note and set time-related globals for convenience
            stamp = datetime.utcnow().strftime("%Y-%m-%dT%H:%M UTC")
            note = (f"{stamp}  temporal aggregation: method={temporal_reduce.lower()}, "
                    f"freq={freq_display}, anchor={int(temporal_anchor_hour)}, tz={tz_label}, "
                    f"label={label_side}, period={'start' if (pe2ps_applied or day_like) else 'as_is'}")
            prev_hist = getattr(dst, "history", "")
            dst.history = (prev_hist + ("\n" if prev_hist else "") + note) if note else prev_hist

            try: dst.temporal_aggregation_method = temporal_reduce.lower()
            except Exception: pass
            try: dst.temporal_aggregation_freq = freq_display
            except Exception: pass
            try: dst.temporal_aggregation_anchor_hour = int(temporal_anchor_hour)
            except Exception: pass
            try: dst.time_timezone = tz_label
            except Exception: pass
            try: dst.time_period_convention = ("period_start" if (pe2ps_applied or day_like) else "as_is")
            except Exception: pass

        finally:
            dst.close()

    finally:
        src.close()

    if verbose:
        print(f"[t-agg] wrote {out_nc_path.name}")
    return out_nc_path
