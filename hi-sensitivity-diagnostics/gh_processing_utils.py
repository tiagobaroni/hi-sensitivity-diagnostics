# gh_processing_utils.py
# Changelog:
# - Reduced to statistical utilities required by hi-sensitivity-diagnostics.
# - Standardized header under GPLv3.
# - Clarified role as supporting statistical diagnostics layer.

# =========================================================================
# PROJECT: hi-sensitivity-diagnostics
#
# DESCRIPTION:
# Statistical and numerical utility module supporting spectral
# sensitivity diagnostics for 21cm HI observations.
#
# This file provides robust statistical estimators and time-series
# diagnostics used for:
#   • Spectral line-free RMS estimation (MAD-based)
#   • Robust location and dispersion metrics
#   • Drift detection support
#   • Allan deviation computation (stability analysis)
#   • Instrumental floor modeling support
#
# The utilities here are intentionally general-purpose and contain
# no observational decision logic. All physical interpretation is
# handled by gh_sensitivity_planner.py.
#
# AUTHOR: Tiago Henrique França Baroni
# EMAIL: t.baroni (at) gmail (dot) com
# LINK: https://github.com/tiagobaroni/hi-sensitivity-diagnostics
# LICENSE: GNU General Public License v3.0 (GPL-3.0)
# =========================================================================
#
# Copyright (C) 2026 Tiago Henrique França Baroni
#
# This file is part of hi-sensitivity-diagnostics.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# =========================================================================
#
# FUNCTIONAL ROLE:
#
# Provides robust statistical primitives required by the sensitivity
# diagnostics workflow, including:
#
#   1. Robust sigma estimation via MAD scaling
#   2. Robust location-scale estimation
#   3. Drift and variability support metrics
#   4. Time-series stability analysis (Allan deviation)
#
# These utilities are:
#   • Spectrally agnostic
#   • Instrument-agnostic
#   • Free of integration-time decision logic
#
# Their purpose is to supply statistically reliable inputs to the
# physical modeling layer implemented in gh_sensitivity_planner.py.
# =========================================================================

from __future__ import annotations

from typing import Tuple

import numpy as np


def compute_mad_sigma(data: np.ndarray) -> float:
    """
    Compute a robust sigma estimate using MAD scaled to Gaussian sigma.

    Args:
        data: 1D array of samples (NaNs allowed).

    Returns:
        Robust sigma (float). Returns NaN if no finite samples exist.
    """
    x = np.asarray(data, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    return float(1.4826 * mad)


def compute_robust_location_scale(data: np.ndarray) -> Tuple[float, float]:
    """
    Compute robust location (median) and scale (MAD sigma).

    Args:
        data: 1D array of samples (NaNs allowed).

    Returns:
        (location, scale) as floats. Returns (NaN, NaN) if no finite samples exist.
    """
    x = np.asarray(data, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan"), float("nan")
    loc = float(np.median(x))
    mad = float(np.median(np.abs(x - loc)))
    scale = float(1.4826 * mad)
    return loc, scale


def compute_allan_deviation(
    series: np.ndarray,
    tau0: float,
    max_m: int | None = None,
    n_points: int = 24,
) -> dict:
    """
    Compute Allan deviation for a scalar time series sampled at tau0.

    This implementation is intended for *diagnostic* stability analysis
    (e.g., total power drift/stability), not for integration-time decisions
    by itself.

    The Allan variance for averaging factor m is:
        avar(m) = 0.5 * mean( (ȳ_{k+1}(m) - ȳ_k(m))^2 )

    where ȳ_k(m) is the average of m consecutive samples.

    Args:
        series: 1D array of samples (NaNs allowed).
        tau0: Base sampling interval in seconds (per sample).
        max_m: Maximum averaging factor. If None, chosen from data length.
        n_points: Number of averaging factors to evaluate (log-spaced).

    Returns:
        dict with:
          - taus: ndarray of averaging times (seconds)
          - adev: ndarray of Allan deviation (same units as series)
          - tau_opt_s: float, tau at minimum adev (NaN if undefined)
    """
    y = np.asarray(series, dtype=float)
    y = y[np.isfinite(y)]
    if y.size < 4 or tau0 <= 0:
        return {"taus": np.array([]), "adev": np.array([]), "tau_opt_s": float("nan")}

    n = y.size
    if max_m is None:
        # Conservative: need at least 2 cluster differences => 2*(m) <= n - m
        max_m = max(1, n // 4)
    max_m = int(max(1, min(max_m, n // 2)))

    # Log-spaced averaging factors, unique ints in [1, max_m]
    m_vals = np.unique(np.clip(np.round(np.logspace(0, np.log10(max_m), n_points)), 1, max_m).astype(int))
    taus = []
    adev = []

    for m in m_vals:
        k_max = n // m
        if k_max < 3:
            continue

        # Compute cluster averages
        y_trim = y[: k_max * m]
        y_bar = y_trim.reshape(k_max, m).mean(axis=1)

        # Two-sample Allan variance
        dy = np.diff(y_bar)
        avar = 0.5 * np.mean(dy * dy)
        if np.isfinite(avar) and avar >= 0:
            taus.append(m * tau0)
            adev.append(np.sqrt(avar))

    taus = np.asarray(taus, dtype=float)
    adev = np.asarray(adev, dtype=float)

    if taus.size == 0:
        tau_opt = float("nan")
    else:
        tau_opt = float(taus[int(np.nanargmin(adev))])

    return {"taus": taus, "adev": adev, "tau_opt_s": tau_opt}
