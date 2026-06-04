#!/usr/bin/env python3
"""Bagpipes backend for saguiSED.

This script is called by saguiSED::fit_region_seds().  It intentionally uses a
fast optimisation fit rather than a full posterior sampler; it is meant as a
first-look regional SED fitting backend and as a stable bridge between R and
Python.  A full posterior mode can be added later without changing the R-facing
object structure.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, least_squares

import bagpipes


FILTER_COLORS = [
    "#3F4EA1",
    "#315FB9",
    "#2577D9",
    "#1FA2D6",
    "#28B5B5",
    "#39BF75",
    "#90D94A",
    "#E4D642",
    "#F0C145",
    "#EE9A3A",
    "#E97832",
    "#D85A2A",
    "#C9432B",
    "#A9362D",
]

SAGUI_BLUE = "#213E60"
SAGUI_ORANGE = "#E68C3A"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--photometry-csv", required=True)
    parser.add_argument("--throughput-csv", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--filters", required=True)
    parser.add_argument("--unit", default="Jy", choices=["Jy", "uJy", "nJy", "10nJy"])
    parser.add_argument("--redshift", required=True, type=float)
    parser.add_argument("--systematic-frac", default=0.06, type=float)
    parser.add_argument("--sfh", default="exponential", choices=["exponential"])
    parser.add_argument("--dust", default="Calzetti", choices=["Calzetti"])
    parser.add_argument("--metallicity", default="free", choices=["free", "fixed"])
    parser.add_argument("--fixed-metallicity", default=1.0, type=float)
    parser.add_argument("--regions", default="")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def unit_to_ujy(values: np.ndarray, unit: str) -> np.ndarray:
    if unit == "Jy":
        return values * 1.0e6
    if unit == "uJy":
        return values
    if unit == "nJy":
        return values * 1.0e-3
    if unit == "10nJy":
        return values * 0.01
    raise ValueError(f"Unsupported unit: {unit}")


def infer_wavelengths(filters: list[str], throughput: pd.DataFrame) -> np.ndarray:
    waves = []
    wave_col = "wavelength_um" if "wavelength_um" in throughput.columns else "wavelength"
    for filt in filters:
        curve = throughput.loc[throughput["filter"] == filt]
        if curve.empty:
            raise ValueError(f"Cannot infer central wavelength for filter {filt}")
        if "lambda_eff_um" in curve.columns and np.isfinite(curve["lambda_eff_um"]).any():
            waves.append(float(np.nanmedian(curve["lambda_eff_um"])))
            continue
        weights = curve["throughput"].to_numpy(dtype=float)
        wave = curve[wave_col].to_numpy(dtype=float)
        if np.nanmax(weights) <= 0:
            waves.append(float(np.nanmedian(wave)))
        else:
            waves.append(float(np.average(wave, weights=np.maximum(weights, 0))))
    return np.asarray(waves, dtype=float)


def write_filter_files(filters: list[str], throughput: pd.DataFrame, out_dir: Path) -> list[str]:
    filter_dir = out_dir / "bagpipes_filters"
    filter_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    wave_col = "wavelength_um" if "wavelength_um" in throughput.columns else "wavelength"
    for filt in filters:
        curve = throughput.loc[throughput["filter"] == filt]
        if curve.empty:
            raise ValueError(f"Filter {filt} not found in throughput table.")
        path = filter_dir / f"{filt}.dat"
        np.savetxt(
            path,
            np.column_stack(
                [
                    curve[wave_col].to_numpy(dtype=float) * 1.0e4,
                    curve["throughput"].to_numpy(dtype=float),
                ]
            ),
        )
        paths.append(str(path))
    return paths


def make_components(params: np.ndarray, args: argparse.Namespace) -> dict:
    log_mass, age_gyr, log_tau_gyr, log_metallicity, av = params
    metallicity = args.fixed_metallicity if args.metallicity == "fixed" else 10.0**log_metallicity

    return {
        "redshift": args.redshift,
        "exponential": {
            "age": age_gyr,
            "tau": 10.0**log_tau_gyr,
            "massformed": log_mass,
            "metallicity": metallicity,
        },
        "dust": {"type": args.dust, "Av": av},
    }


def fit_one_region(
    row: pd.Series,
    filters: list[str],
    wave_um: np.ndarray,
    filt_paths: list[str],
    args: argparse.Namespace,
) -> tuple[dict, pd.DataFrame, np.ndarray]:
    obs = unit_to_ujy(np.asarray([row[filt] for filt in filters], dtype=float), args.unit)
    stat_err = unit_to_ujy(np.asarray([row[f"{filt}_err"] for filt in filters], dtype=float), args.unit)
    err = np.sqrt(stat_err**2 + (args.systematic_frac * np.maximum(obs, 1.0e-12)) ** 2)
    region = int(row["region"])

    spec_wavs = np.linspace(max(1000.0, 0.82 * wave_um.min() * 1.0e4), 1.12 * wave_um.max() * 1.0e4, 1400)
    start = np.array([10.0, 3.0, 0.0, np.log10(0.5), 0.4])

    model = bagpipes.model_galaxy(
        make_components(start, args),
        filt_list=filt_paths,
        spec_wavs=spec_wavs,
        spec_units="mujy",
        phot_units="mujy",
    )

    def model_photometry(params: np.ndarray) -> np.ndarray:
        try:
            model.update(make_components(params, args))
            phot = np.asarray(model.photometry, dtype=float)
            if not np.all(np.isfinite(phot)):
                return np.full_like(obs, 1.0e30)
            return phot
        except Exception:
            return np.full_like(obs, 1.0e30)

    def residuals(params: np.ndarray) -> np.ndarray:
        return (model_photometry(params) - obs) / err

    metallicity_bounds = (np.log10(0.005), np.log10(2.5))
    if args.metallicity == "fixed":
        metallicity_bounds = (np.log10(args.fixed_metallicity), np.log10(args.fixed_metallicity))

    lower = np.array([7.0, 0.08, -1.3, metallicity_bounds[0], 0.0])
    upper = np.array([12.8, 7.2, 1.0, metallicity_bounds[1], 3.5])
    if args.metallicity == "fixed":
        lower[3] -= 1.0e-8
        upper[3] += 1.0e-8

    starts = [
        start,
        [10.0, 1.0, -0.3, np.log10(max(args.fixed_metallicity, 0.005)), 0.2],
        [10.5, 5.0, 0.3, np.log10(0.2), 0.8],
        [9.8, 2.0, -0.6, np.log10(0.4), 0.0],
        [10.8, 6.5, 0.6, np.log10(1.0), 1.2],
    ]

    starts = [np.minimum(np.maximum(np.asarray(x, dtype=float), lower + 1.0e-7), upper - 1.0e-7) for x in starts]
    fits = [
        least_squares(
            residuals,
            trial,
            bounds=(lower, upper),
            max_nfev=450,
            xtol=1.0e-5,
            ftol=1.0e-5,
            gtol=1.0e-5,
        )
        for trial in starts
    ]
    best = min(fits, key=lambda fit: np.sum(fit.fun**2))

    if args.metallicity == "free":
        global_fit = differential_evolution(
            lambda x: np.sum(residuals(x) ** 2),
            list(zip(lower, upper)),
            maxiter=20,
            popsize=8,
            polish=False,
            seed=region + 2026,
            workers=1,
            updating="immediate",
        )
        polished = least_squares(
            residuals,
            global_fit.x,
            bounds=(lower, upper),
            max_nfev=550,
            xtol=1.0e-5,
            ftol=1.0e-5,
            gtol=1.0e-5,
        )
        if np.sum(polished.fun**2) < np.sum(best.fun**2):
            best = polished

    model.update(make_components(best.x, args))
    model_ujy = np.asarray(model.photometry, dtype=float)
    pulls = (obs - model_ujy) / err
    rms = float(np.sqrt(np.mean(pulls**2)))
    log_mass, age_gyr, log_tau_gyr, log_metallicity, av = best.x
    metallicity = args.fixed_metallicity if args.metallicity == "fixed" else 10.0**log_metallicity

    summary = {
        "region": region,
        "n_pix": int(row["n_pix"]) if "n_pix" in row.index and np.isfinite(row["n_pix"]) else np.nan,
        "redshift": args.redshift,
        "sfh": args.sfh,
        "dust": args.dust,
        "metallicity_mode": args.metallicity,
        "logMformed": log_mass,
        "age_gyr": age_gyr,
        "tau_gyr": 10.0**log_tau_gyr,
        "metallicity": metallicity,
        "logZ_Zsun": np.log10(max(metallicity, 1.0e-12)),
        "Av": av,
        "fit_chi2": float(np.sum(((model_ujy - obs) / err) ** 2)),
        "fit_rms_sigma": rms,
        "systematic_frac": args.systematic_frac,
    }

    photometry = pd.DataFrame(
        {
            "region": region,
            "filter": filters,
            "wave_um": wave_um,
            "obs_ujy": obs,
            "err_ujy": err,
            "model_ujy": model_ujy,
            "pull": pulls,
        }
    )
    return summary, photometry, np.asarray(model.spectrum, dtype=float)


def plot_region(
    region: int,
    throughput: pd.DataFrame,
    photometry: pd.DataFrame,
    spectrum: np.ndarray,
    summary: dict,
    out_dir: Path,
) -> None:
    plot_dir = out_dir / "sed_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    y_values = np.r_[photometry["obs_ujy"].to_numpy(), photometry["model_ujy"].to_numpy()]
    y_min = float(np.nanmin(y_values) * 0.55)
    y_max = float(np.nanmax(y_values) * 1.28)
    filter_base = y_min * 1.08
    filter_height = (y_max - y_min) * 0.20

    fig, ax = plt.subplots(figsize=(8.2, 4.7), dpi=220)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    filters = list(photometry["filter"])
    for i, filt in enumerate(filters):
        curve = throughput.loc[throughput["filter"] == filt]
        wave_col = "wavelength_um" if "wavelength_um" in curve.columns else "wavelength"
        x = curve[wave_col].to_numpy(dtype=float)
        y = curve["throughput"].to_numpy(dtype=float)
        if np.nanmax(y) > 0:
            y = y / np.nanmax(y)
        ax.fill_between(
            x,
            filter_base,
            filter_base + filter_height * y,
            color=FILTER_COLORS[i % len(FILTER_COLORS)],
            alpha=0.28,
            linewidth=0,
            zorder=0,
        )

    visible = np.isfinite(spectrum[:, 1])
    ax.plot(
        spectrum[visible, 0] / 1.0e4,
        spectrum[visible, 1],
        color=SAGUI_BLUE,
        linewidth=2.35,
        label="Bagpipes model",
        zorder=2,
    )
    ax.plot(
        photometry["wave_um"],
        photometry["model_ujy"],
        color=SAGUI_BLUE,
        linewidth=1.1,
        alpha=0.25,
        zorder=3,
    )
    ax.scatter(photometry["wave_um"], photometry["model_ujy"], s=22, color=SAGUI_BLUE, zorder=4)
    ax.errorbar(
        photometry["wave_um"],
        photometry["obs_ujy"],
        yerr=photometry["err_ujy"],
        fmt="o",
        ms=7.4,
        mfc="white",
        mec=SAGUI_ORANGE,
        mew=1.9,
        ecolor="#424242",
        elinewidth=1.25,
        capsize=2.8,
        label="Observed photometry",
        zorder=5,
    )

    ax.text(
        0.03,
        0.94,
        f"region {region}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="bold",
        color=SAGUI_BLUE,
    )
    ax.text(
        0.03,
        0.875,
        f"RMS = {summary['fit_rms_sigma']:.2f}\u03c3",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        color="#555555",
    )

    ax.set_xlim(0.92 * float(photometry["wave_um"].min()), 1.07 * float(photometry["wave_um"].max()))
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("Observed wavelength [\u00b5m]", fontsize=12)
    ax.set_ylabel("Flux density [\u00b5Jy]", fontsize=12)
    ax.legend(loc="upper right", frameon=False, fontsize=10)
    ax.tick_params(colors="#1b1b1b", labelsize=10)
    for spine in ax.spines.values():
        spine.set_color("#222222")
        spine.set_linewidth(0.9)
    fig.tight_layout(pad=0.55)
    fig.savefig(plot_dir / f"region_{region:03d}_bagpipes_fit.png", facecolor="white")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    filters = [x.strip() for x in args.filters.split(",") if x.strip()]
    throughput = pd.read_csv(args.throughput_csv)
    photometry = pd.read_csv(args.photometry_csv)
    wave_um = infer_wavelengths(filters, throughput)
    filt_paths = write_filter_files(filters, throughput, out_dir)

    if args.regions:
        region_ids = [int(x) for x in args.regions.split(",") if x.strip()]
    else:
        region_ids = [int(x) for x in photometry["region"].dropna().unique()]

    summaries = []
    photometry_rows = []
    spectrum_rows = []
    for region in region_ids:
        region_rows = photometry.loc[photometry["region"] == region]
        if region_rows.empty:
            raise ValueError(f"Region {region} not found in photometry table.")
        summary, model_photometry, spectrum = fit_one_region(
            region_rows.iloc[0],
            filters,
            wave_um,
            filt_paths,
            args,
        )
        summaries.append(summary)
        photometry_rows.append(model_photometry)
        spectrum_rows.append(
            pd.DataFrame(
                {
                    "region": region,
                    "wave_um": spectrum[:, 0] / 1.0e4,
                    "model_ujy": spectrum[:, 1],
                }
            )
        )
        plot_region(region, throughput, model_photometry, spectrum, summary, out_dir)

    summary_df = pd.DataFrame(summaries)
    model_photometry_df = pd.concat(photometry_rows, ignore_index=True)
    model_spectrum_df = pd.concat(spectrum_rows, ignore_index=True)
    summary_df.to_csv(out_dir / "region_fit_summary.csv", index=False)
    model_photometry_df.to_csv(out_dir / "model_photometry.csv", index=False)
    model_spectrum_df.to_csv(out_dir / "model_spectrum.csv", index=False)

    print(f"Wrote {out_dir / 'region_fit_summary.csv'}")
    print(f"Wrote {out_dir / 'model_photometry.csv'}")
    print(f"Wrote {out_dir / 'model_spectrum.csv'}")
    print(f"Wrote SED plots under {out_dir / 'sed_plots'}")


if __name__ == "__main__":
    main()
