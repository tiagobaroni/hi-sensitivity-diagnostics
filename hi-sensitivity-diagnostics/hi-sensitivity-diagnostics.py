# hi-sensitivity-diagnostics.py
# Changelog:
# - Fixed import typos and parsing artifacts.
# - Added CLI flags: --allan and --cold_sky (implemented).
# - Added Allan deviation diagnostics (total-power time series).
# - Added cold-sky reference comparison (diagnostic only).
# - Preserved core planning logic based on spectral line-free RMS.

# =========================================================================
# PROJECT: hi-sensitivity-diagnostics
#
# DESCRIPTION:
# Diagnostic and integration-time planning module for 21cm HI observations.
# Implements a physically consistent framework for separating thermal noise
# from systematic instrumental effects using spectral RMS analysis.
#
# Designed to operate on reduced FITS spectra (e.g., SalsaSpectrum output),
# this module evaluates observational regimes, estimates instrumental floor,
# detects drift, and provides objective integration-time recommendations.
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
# SCIENTIFIC PHILOSOPHY:
#
# This module formalizes the distinction between:
#
#   1) spectral_linefree_rms  → Primary decision metric
#   2) total_power_metric     → Diagnostic only (drift/stability)
#
# WARNING:
# total_power is diagnostic only.
# Integration-time decisions must rely on spectral_linefree_rms.
#
# =========================================================================

from __future__ import annotations

import os
import argparse
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple

import numpy as np
import tkinter as tk
from tkinter import filedialog

import matplotlib.pyplot as plt
from fpdf import FPDF
from fpdf.enums import XPos, YPos

def _configure_pdf_fonts(pdf: FPDF) -> str:
    """
    Configure Unicode-capable fonts for FPDF2.

    Returns the font family name to be used (preferred: 'DejaVu').
    If DejaVu fonts are not found, falls back to core fonts and the
    caller should avoid Unicode (Greek symbols etc.) in PDF text.
    """
    # Look for fonts in common repo locations
    candidate_dirs = [
        Path(__file__).resolve().parent,
        Path(__file__).resolve().parent / "fonts",
        Path.cwd(),
        Path.cwd() / "fonts",
    ]

    def _find_font(filename: str) -> Optional[Path]:
        for d in candidate_dirs:
            p = d / filename
            if p.exists() and p.is_file():
                return p
        return None

    regular = _find_font("DejaVuSans.ttf")
    bold = _find_font("DejaVuSans-Bold.ttf")
    italic = _find_font("DejaVuSans-Oblique.ttf")
    bold_italic = _find_font("DejaVuSans-BoldOblique.ttf")

    if regular and bold:
        # 'uni=True' enables Unicode text
        pdf.add_font("DejaVu", "", str(regular))
        pdf.add_font("DejaVu", "B", str(bold))
        if italic:
            pdf.add_font("DejaVu", "I", str(italic))
        if bold_italic:
            pdf.add_font("DejaVu", "BI", str(bold_italic))
        return "DejaVu"

    # Fallback: core fonts (NOT Unicode-safe)
    return "Helvetica"


def _pdf_safe(text: str, unicode_ok: bool) -> str:
    """
    When Unicode fonts are not available, replace characters that will
    crash PDF generation (e.g., τ, σ, em-dash).
    """
    if unicode_ok:
        return text

    # Minimal transliteration for common symbols used in this project
    replacements = {
        "τ": "tau",
        "σ": "sigma",
        "Δ": "Delta",
        "≤": "<=",
        "≥": ">=",
        "≈": "~=",
        "—": "-",
        "–": "-",
        "µ": "u",
    }
    out = text
    for k, v in replacements.items():
        out = out.replace(k, v)

    # Strip any remaining non-latin1 characters defensively
    out = out.encode("latin-1", errors="ignore").decode("latin-1")
    return out



# Internal dependencies
from salsaspectrum.spectrum import SalsaSpectrum
from options import DEFAULT_PROCESS_OPTIONS, ProcessOptions, with_overrides
import gh_processing_utils as utils


def _collect_fits_files(folder_path: str) -> List[Path]:
    path = Path(folder_path)
    all_files = list(path.iterdir())
    files = sorted({f for f in all_files if f.is_file() and f.suffix.lower() == ".fits"})
    return files


def analyze_prospecting_sample(
    folder_path: str,
    tau_s_fixed: float,
    options: ProcessOptions,
    n_steps: int = 10,
) -> Dict[str, Any]:
    """
    Analyze a prospecting (pilot) sample to establish noise baseline and antenna floor.

    Outputs are designed for:
      - RMS stacking efficiency audit (RMS vs N)
      - Quadratic decomposition (sigma_obs^2 = sigma_th^2 + sigma_sys^2)
      - Optional Allan deviation (total power stability)
    """
    files = _collect_fits_files(folder_path)
    if not files:
        raise FileNotFoundError(f"No FITS files found in {folder_path}")

    individual_results: List[Dict[str, Any]] = []
    total_power_series: List[float] = []

    for f in files:
        try:
            spec = SalsaSpectrum(str(f))
            n_chan = len(spec.data)
            if n_chan < 10:
                continue

            # Identify line-free regions using 15% edge heuristic
            edge_size = max(1, int(n_chan * 0.15))
            line_free_windows = [0, edge_size, n_chan - edge_size, n_chan - 1]

            # Total power metric (diagnostic only) BEFORE baseline subtraction
            # Note: we keep this as a simple scalar per subintegration/file.
            total_power_i = float(np.nanmean(spec.data))
            if np.isfinite(total_power_i):
                total_power_series.append(total_power_i)

            # Physics metadata extraction
            t_sys = spec.get_keyword("TSYS", options.t_sys_est)
            t_int = spec.get_keyword("EXPTIME", spec.get_keyword("ONTIME", tau_s_fixed))

            spec.calc_theoretical_noise(
                t_sys=t_sys,
                integ_time=t_int,
                n_pol=options.n_pol,
                enbw_factor=options.enbw_factor,
            )

            # Baseline fitting to isolate thermal noise
            spec.fit_baseline(order=1, windows=line_free_windows, coord="pix")
            spec.subtract_baseline()

            # Define mask for robust noise measurement (line-free edges)
            is_base_mask = np.zeros(n_chan, dtype=bool)
            is_base_mask[0:edge_size] = True
            is_base_mask[(n_chan - edge_size) :] = True

            sigma_obs_i = utils.compute_mad_sigma(spec.data[is_base_mask])

            if np.isfinite(sigma_obs_i) and np.isfinite(spec.sigma_theo):
                individual_results.append(
                    {
                        "sigma_obs": float(sigma_obs_i),
                        "sigma_th": float(spec.sigma_theo),
                        "data": np.asarray(spec.data, dtype=float),
                        "mask": is_base_mask,
                    }
                )

        except Exception:
            continue

    if not individual_results:
        raise ValueError("FATAL ERROR: Informação insuficiente para garantir código funcional.")

    # Stacking audit for marginal efficiency (RMS vs N)
    n_files = len(individual_results)
    n_steps = max(2, int(n_steps))
    efficiency_report: List[Tuple[int, float]] = []
    unique_counts = sorted(list(set([max(1, int(p * n_files)) for p in np.linspace(0.1, 1.0, n_steps)])))

    for count in unique_counts:
        subset = individual_results[:count]
        stacked_data = np.mean([item["data"] for item in subset], axis=0)
        stack_rms = utils.compute_mad_sigma(stacked_data[individual_results[0]["mask"]])
        efficiency_report.append((count, float(stack_rms)))

    return {
        "n_prospecting": int(n_files),
        "sigma_prospecting_stack": float(efficiency_report[-1][1]),
        "sigma_th_med_file": float(np.nanmedian([item["sigma_th"] for item in individual_results])),
        "individual_sigmas": [float(res["sigma_obs"]) for res in individual_results],
        "efficiency": efficiency_report,
        "folder_name": Path(folder_path).name,
        "total_power_series": total_power_series,
    }


def analyze_cold_sky_reference(
    folder_path: str,
    tau_s_fixed: float,
    options: ProcessOptions,
    n_steps: int = 10,
) -> Dict[str, Any]:
    """
    Analyze a cold-sky reference dataset with the exact same pipeline as the prospecting sample.

    IMPORTANT:
    - This is diagnostic only and must not be used as the sole metric for integration-time decisions.
    - It helps detect unexpected instrumental/state changes across sessions.
    """
    results = analyze_prospecting_sample(folder_path, tau_s_fixed, options, n_steps=n_steps)
    return results


def _plot_efficiency_curve(folder_name: str, efficiency: List[Tuple[int, float]]) -> str:
    n_vals = [i[0] for i in efficiency]
    rms_vals = [i[1] for i in efficiency]
    ideal_rms = [rms_vals[0] * np.sqrt(n_vals[0] / n) for n in n_vals]

    plt.figure(figsize=(8, 5))
    plt.plot(n_vals, rms_vals, "bo-", label="Prospecting RMS (Observed)")
    plt.plot(n_vals, ideal_rms, "r--", label="Ideal 1/sqrt(N)")
    plt.title(f"Antenna Efficiency: {folder_name}")
    plt.xlabel("Sample Count (N)")
    plt.ylabel("RMS [K]")
    plt.grid(True, alpha=0.3)
    plt.legend()

    plot_path = f"plot_efficiency_{folder_name}.png"
    plt.savefig(plot_path, dpi=120, bbox_inches="tight")
    plt.close()
    return plot_path


def _plot_allan(folder_name: str, taus: np.ndarray, adev: np.ndarray) -> str:
    plt.figure(figsize=(8, 5))
    plt.loglog(taus, adev, "o-")
    plt.title(f"Allan Deviation (Total Power): {folder_name}")
    plt.xlabel("Averaging time τ [s]")
    plt.ylabel("Allan deviation [arb.]")
    plt.grid(True, which="both", alpha=0.3)

    plot_path = f"plot_allan_{folder_name}.png"
    plt.savefig(plot_path, dpi=120, bbox_inches="tight")
    plt.close()
    return plot_path


def generate_report(
    results: Dict[str, Any],
    sigma_target: float,
    n_planned: int,
    n_max: int,
    tau_s: float,
    allan_result: Optional[Dict[str, Any]] = None,
    cold_sky_result: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Generate a campaign planning report evaluating if goals are achievable.
    Predict required N or methodological changes based on the estimated noise floor.

    Decisions are based on spectral_linefree_rms only. Any total-power and cold-sky
    analysis is diagnostic and used only as cautionary information.
    """
    n_prospecting = int(results["n_prospecting"])
    sigma_obs_prospecting = float(results["sigma_prospecting_stack"])
    sigma_th_file = float(results["sigma_th_med_file"])
    r_max_criterion = 1.3

    # Fundamental law: sigma_file = sigma_prospecting * sqrt(n_prospecting)
    sigma_per_file = sigma_obs_prospecting * np.sqrt(n_prospecting)

    # Expected campaign performance with n_planned
    sigma_at_n_planned = sigma_per_file / np.sqrt(max(1, int(n_planned)))

    # Floor estimation via quadratic decomposition
    sigma_th_prospecting = sigma_th_file / np.sqrt(n_prospecting)
    sigma_sys_hat = float(np.sqrt(max(0.0, sigma_obs_prospecting**2 - sigma_th_prospecting**2)))

    # Required N to reach goal (thermal scaling assumption)
    n_req = int(np.ceil((sigma_per_file / sigma_target) ** 2))

    # Suggested outlier threshold for downstream pipeline
    med_i, mad_sigma_i = utils.compute_robust_location_scale(np.asarray(results["individual_sigmas"], dtype=float))
    suggested_max_rms = float(med_i + 3.0 * mad_sigma_i)

    plot_eff_path = _plot_efficiency_curve(results["folder_name"], results["efficiency"])

    # PDF generation (fpdf2 modern syntax)
    pdf = FPDF()
    pdf.add_page()
    font_family = _configure_pdf_fonts(pdf)
    unicode_ok = (font_family != "Helvetica")
    pdf.set_font(font_family, "B", 16)
    pdf.cell(0, 10, _pdf_safe("HI Sensitivity Diagnostics - Planning Report", unicode_ok), new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")

    pdf.set_font(font_family, "", 11)
    pdf.ln(5)
    pdf.cell(0, 8, _pdf_safe(f"Prospecting files: {n_prospecting}", unicode_ok), new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.cell(0, 8, _pdf_safe(f"Initially planned N: {n_planned}", unicode_ok), new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.cell(0, 8, _pdf_safe(f"Target RMS: {sigma_target:.4f} K", unicode_ok), new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.cell(0, 8, _pdf_safe(f"Observed RMS (stack): {sigma_obs_prospecting:.4f} K", unicode_ok), new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.cell(0, 8, _pdf_safe(f"Theoretical RMS (median file): {sigma_th_file:.4f} K", unicode_ok), new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.cell(0, 8, _pdf_safe(f"Estimated systematic floor: {sigma_sys_hat:.4f} K", unicode_ok), new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    pdf.set_font(font_family, "", 9)
    pdf.multi_cell(
        0,
        5,
        "Note: sigma_sys is estimated via quadratic decomposition. "
        f"Thermal regime criterion adopted: Rmax = {r_max_criterion}.",
    )
    pdf.ln(3)

    # Cold sky diagnostics (optional)
    if cold_sky_result is not None:
        sigma_cold = float(cold_sky_result["sigma_prospecting_stack"])
        if np.isfinite(sigma_cold) and sigma_obs_prospecting > 0:
            dev_pct = 100.0 * (sigma_cold - sigma_obs_prospecting) / sigma_obs_prospecting
        else:
            dev_pct = float("nan")

        pdf.set_font(font_family, "B", 11)
        pdf.cell(0, 8, _pdf_safe("Cold-sky reference (diagnostic):", unicode_ok), new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        pdf.set_font(font_family, "", 11)
        pdf.cell(0, 8, _pdf_safe(f"Cold-sky RMS (stack): {sigma_cold:.4f} K", unicode_ok), new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        if np.isfinite(dev_pct):
            pdf.cell(0, 8, _pdf_safe(f"Deviation vs prospecting: {dev_pct:+.1f} %", unicode_ok), new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            if abs(dev_pct) > 20.0:
                pdf.set_font(font_family, "", 10)
                pdf.multi_cell(
                    0,
                    5,
                    "WARNING: deviation > 20% suggests session-to-session changes "
                    "(instrument state, RFI, calibration or pointing).",
                )
                pdf.set_font(font_family, "", 11)
        pdf.ln(2)

    pdf.set_fill_color(230, 230, 230)
    pdf.set_font(font_family, "B", 12)
    pdf.cell(0, 10, _pdf_safe("CAMPAIGN STRATEGY ADVICE:", unicode_ok), new_x=XPos.LMARGIN, new_y=YPos.NEXT, fill=True)
    pdf.set_font(font_family, "", 11)

    # Decision matrix (spectral RMS only)
    if sigma_target <= sigma_sys_hat and np.isfinite(sigma_sys_hat):
        msg = (
            f"SYSTEMATIC LIMIT: Target ({sigma_target:.4f} K) is below the floor ({sigma_sys_hat:.4f} K).\n"
            "ACTION: Thermal scaling (increasing N) is physically ineffective in this regime.\n"
            "ADVICE: 1) Relax scientific target; 2) Reduce sigma_sys (RFI/baseline); or 3) Change strategy."
        )
    elif sigma_at_n_planned <= sigma_target:
        msg = (
            f"SUCCESS: Your planned N={n_planned} is sufficient to reach the target.\n"
            "ADVICE: Proceed with the campaign as originally scheduled."
        )
    elif n_req <= n_max:
        msg = (
            f"ADJUSTMENT REQUIRED: Your planned N={n_planned} is insufficient.\n"
            f"ADVICE: Increase sample count per point to N={n_req} for the final campaign."
        )
    else:
        tau_req_total = (n_prospecting * tau_s) * (sigma_obs_prospecting / sigma_target) ** 2
        tau_s_hybrid = tau_req_total / max(1, int(n_max))
        msg = (
            f"LOGISTICAL LIMIT: Required N ({n_req}) exceeds your limit N_max ({n_max}).\n"
            f"ADVICE: Use N={n_max} and increase integration time per sample to ~{tau_s_hybrid:.1f} s."
        )

    pdf.multi_cell(0, 8, msg)
    pdf.ln(3)

    pdf.set_font(font_family, "B", 11)
    pdf.cell(
        0,
        8,
        f"Pipeline hint: set 'max_rms' to {suggested_max_rms:.4f} K",
        new_x=XPos.LMARGIN,
        new_y=YPos.NEXT,
    )

    # Efficiency plot
    y0 = pdf.get_y() + 5
    pdf.image(plot_eff_path, x=15, y=y0, w=180)

    # Allan diagnostics (optional)
    allan_plot_path: Optional[str] = None
    if allan_result is not None and allan_result.get("taus") is not None and allan_result.get("adev") is not None:
        taus = np.asarray(allan_result["taus"], dtype=float)
        adev = np.asarray(allan_result["adev"], dtype=float)
        if taus.size > 0 and adev.size > 0:
            allan_plot_path = _plot_allan(results["folder_name"], taus, adev)

            pdf.add_page()
            pdf.set_font(font_family, "B", 14)
            pdf.cell(0, 10, _pdf_safe("Allan Deviation (Diagnostic)", unicode_ok), new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.set_font(font_family, "", 11)

            tau_opt = allan_result.get("tau_opt_s")
            if tau_opt is not None and np.isfinite(tau_opt):
                pdf.multi_cell(
                    0,
                    6,
                    "This plot is diagnostic only. It helps evaluate stability and whether increasing τ_s "
                    "is beneficial compared to increasing N.\n"
                    f"Estimated stability optimum (minimum Allan deviation): τ_s ≈ {float(tau_opt):.1f} s",
                )
            else:
                pdf.multi_cell(
                    0,
                    6,
                    "This plot is diagnostic only. It helps evaluate stability and whether increasing τ_s "
                    "is beneficial compared to increasing N.",
                )

            pdf.ln(2)
            pdf.image(allan_plot_path, x=15, y=pdf.get_y() + 5, w=180)

    output_filename = f"Sensitivity_Plan_{results['folder_name']}.pdf"
    pdf.output(output_filename)

    # Cleanup
    for p in [plot_eff_path, allan_plot_path]:
        if p and os.path.exists(p):
            try:
                os.remove(p)
            except Exception:
                pass

    print(f"\nReport generated: {output_filename}")
    print("-" * 60)
    print(msg.replace("\n", " "))
    print("-" * 60)


def main() -> None:
    parser = argparse.ArgumentParser(description="HI Prospecting & Campaign Planning Tool")
    parser.add_argument("tau_s", type=float, help="Integration time per file [s]")
    parser.add_argument("sigma_target", type=float, help="Target RMS [K]")
    parser.add_argument("n_planned", type=int, help="Initial planned sample count")
    parser.add_argument("--n_max", type=int, default=100, help="Maximum allowed sample count")
    parser.add_argument("--steps", type=int, default=10, help="Efficiency audit resolution")
    parser.add_argument("--allan", action="store_true", help="Compute Allan deviation (diagnostic) using total-power time series.",)
    parser.add_argument("--cold_sky", type=str, default=None, help="Path to cold-sky reference FITS directory (diagnostic only).",)

    args = parser.parse_args()

    # GUI: prospecting folder
    root = tk.Tk()
    root.withdraw()
    folder = filedialog.askdirectory(title="Select folder with Prospecting FITS files")
    root.destroy()

    if not folder:
        return

    planner_opts = with_overrides(DEFAULT_PROCESS_OPTIONS, integ_time_est=args.tau_s)

    try:
        results = analyze_prospecting_sample(folder, args.tau_s, planner_opts, n_steps=args.steps)

        allan_result: Optional[Dict[str, Any]] = None
        if args.allan:
            power_series = np.asarray(results.get("total_power_series", []), dtype=float)
            allan_result = utils.compute_allan_deviation(power_series, tau0=args.tau_s)

        cold_sky_result: Optional[Dict[str, Any]] = None
        if args.cold_sky is not None:
            cold_sky_result = analyze_cold_sky_reference(args.cold_sky, args.tau_s, planner_opts, n_steps=args.steps)

        generate_report(
            results=results,
            sigma_target=args.sigma_target,
            n_planned=args.n_planned,
            n_max=args.n_max,
            tau_s=args.tau_s,
            allan_result=allan_result,
            cold_sky_result=cold_sky_result,
        )

    except Exception as e:
        print(f"FATAL ERROR: {e}")


if __name__ == "__main__":
    main()