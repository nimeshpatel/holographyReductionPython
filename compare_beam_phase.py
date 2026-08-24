#!/usr/bin/env python3
"""
compare_beam_phase.py - Visualize the boresight drift in the BEAM (pre-FFT) phase map.

Plots the far-field/beam phase maps with and without boresight correction, plus
their difference. The difference map isolates the applied boresight correction
(the per-row drift) because all the real beam structure cancels.

This is the pre-holis_aber2 (pre-FFT) phase, where the drift appears in its rawest
form: a per-row additive phase. Rows correspond to the slow (elevation) scan axis,
i.e. to TIME through the observation, so the abrupt end-of-track drift shows up in
the last rows.

Accepts either:
  - 1-column phaseout.dat (N*N values, row-major)            -> reshaped to (N,N)
  - 4-column rgin.dat (j k amp phase)                         -> phase scattered by (j,k)

Usage:
    python compare_beam_phase.py phaseout_nocal.dat phaseout_cal.dat -o beam_phase_compare.pdf
    python compare_beam_phase.py rgin_nocal.dat rgin_cal.dat --units deg -o cmp.pdf
"""

import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def load_beam_phase(filename):
    """Load a beam phase map from a 1-col (N*N) or 4-col (j k amp phase) file."""
    data = np.loadtxt(filename)
    if data.ndim == 1:
        n = int(round(np.sqrt(data.size)))
        if n * n != data.size:
            raise ValueError(f"{filename}: {data.size} values is not a perfect square")
        return data.reshape((n, n))
    # 2-D: assume columns (j, k, amp, phase) or at least (..., phase) in last col
    if data.shape[1] >= 4:
        j = data[:, 0].astype(int)
        k = data[:, 1].astype(int)
        ph = data[:, 3]
        n = int(max(j.max(), k.max())) + 1
        grid = np.full((n, n), np.nan)
        grid[j, k] = ph
        return grid
    # fallback: last column, square reshape
    ph = data[:, -1]
    n = int(round(np.sqrt(ph.size)))
    return ph.reshape((n, n))


def wrap180(x):
    """Wrap an angle (deg) into (-180, 180]."""
    return (x + 180.0) % 360.0 - 180.0


def main():
    ap = argparse.ArgumentParser(description="Compare beam (pre-FFT) phase maps: nocal vs cal vs difference.")
    ap.add_argument("nocal", help="Beam phase map WITHOUT boresight correction (mode 2)")
    ap.add_argument("cal", help="Beam phase map WITH boresight correction (mode 1)")
    ap.add_argument("-o", "--out", default=None, help="Output figure (PDF/PNG). If omitted, shows interactively.")
    ap.add_argument("--units", choices=["deg", "rad"], default="deg",
                    help="Units of the phase files (default: deg). Difference is wrapped accordingly.")
    ap.add_argument("--diff-clip", type=float, default=None,
                    help="Symmetric color limit for the difference panel (in file units). Default: auto (robust).")
    ap.add_argument("--mask-file", default=None,
                    help="Optional N*N mask file (1=valid,0=masked) to blank masked pixels.")
    args = ap.parse_args()

    nocal = load_beam_phase(args.nocal)
    cal = load_beam_phase(args.cal)
    if nocal.shape != cal.shape:
        raise ValueError(f"Shape mismatch: {args.nocal}{nocal.shape} vs {args.cal}{cal.shape}")
    n = nocal.shape[0]

    # Difference = applied correction (cal - nocal), wrapped to remove 2pi ambiguity
    diff = cal - nocal
    if args.units == "deg":
        diff = wrap180(diff)
        unit_label = "deg"
        to_um = (299792458.0 / 94.5e9 * 1e6) / (4 * np.pi) * (np.pi / 180.0)  # um per deg @94.5GHz
    else:
        diff = (diff + np.pi) % (2 * np.pi) - np.pi
        unit_label = "rad"
        to_um = (299792458.0 / 94.5e9 * 1e6) / (4 * np.pi)  # um per rad

    # Optional mask
    if args.mask_file:
        m = np.loadtxt(args.mask_file).reshape((n, n))
        blank = (m == 0)
        for arr in (nocal, cal, diff):
            arr[blank] = np.nan

    # Robust symmetric color limit for difference
    if args.diff_clip is not None:
        dlim = args.diff_clip
    else:
        finite = diff[np.isfinite(diff)]
        dlim = np.nanpercentile(np.abs(finite), 99) if finite.size else 1.0
        dlim = max(dlim, 1e-6)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))

    # Shared scale for the two phase panels (robust)
    both = np.concatenate([nocal[np.isfinite(nocal)].ravel(),
                           cal[np.isfinite(cal)].ravel()])
    vlo, vhi = np.nanpercentile(both, [2, 98])

    for ax, arr, title in (
        (axes[0], nocal, f"WITHOUT boresight cal (mode 2)"),
        (axes[1], cal,   f"WITH boresight cal (mode 1)"),
    ):
        im = ax.imshow(arr, origin="lower", cmap="twilight", vmin=vlo, vmax=vhi, aspect="equal")
        ax.set_title(title)
        ax.set_xlabel("azimuth pixel (fast scan)")
        ax.set_ylabel("row  ~ elevation ~ TIME")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=f"beam phase [{unit_label}]")

    rms = np.sqrt(np.nanmean(diff[np.isfinite(diff)] ** 2))
    im = axes[2].imshow(diff, origin="lower", cmap="RdBu_r", vmin=-dlim, vmax=dlim, aspect="equal")
    axes[2].set_title(f"DIFFERENCE (cal - nocal) = applied correction\n"
                      f"RMS = {rms:.2f} {unit_label}  (~{rms*to_um:.1f} um)")
    axes[2].set_xlabel("azimuth pixel (fast scan)")
    axes[2].set_ylabel("row  ~ elevation ~ TIME")
    fig.colorbar(im, ax=axes[2], fraction=0.046, pad=0.04, label=f"correction [{unit_label}]")

    fig.suptitle("Beam (pre-FFT) phase: boresight drift appears as per-row banding", y=1.02)
    fig.tight_layout()

    if args.out:
        fig.savefig(args.out, dpi=150, bbox_inches="tight")
        print(f"Saved {args.out}")
    else:
        plt.show()

    # Also print the per-row mean of the difference: this IS the drift vs row/time.
    row_mean = np.nanmean(diff, axis=1)
    print("\nPer-row mean of applied correction (row -> elevation -> time):")
    print(f"  range: {np.nanmin(row_mean):.2f} to {np.nanmax(row_mean):.2f} {unit_label}")
    print(f"  this is the boresight drift as it maps onto the map rows.")


if __name__ == "__main__":
    main()
