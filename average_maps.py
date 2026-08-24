#!/usr/bin/env python3
"""
average_maps.py

Average multiple holography surface error maps (Epr.dat files).

Each input file is a flat single-column ASCII file containing N*N phase
values in radians (or surface displacement in µm), written row-major.
The grid size N is either inferred automatically from the file length
(must be a perfect square) or specified explicitly with --ndim.

Masked/invalid pixels are marked as -9999 in Holis convention.
These are excluded from the average at each pixel; the output pixel
is set to -9999 only if ALL input maps have that pixel masked.

Usage:
    python average_maps.py Epr1.dat Epr2.dat Epr3.dat Epr4.dat

    python average_maps.py results/*_Epr.dat -o mean_Epr.dat

    python average_maps.py Epr1.dat Epr2.dat --ndim 128 --output mean.dat \\
        --verbose --mask mask128.dat

Options:
    -o / --output   Output filename (default: mean_Epr.dat)
    --ndim N        Grid size override (auto-detected if omitted)
    --mask FILE     Optional mask file; masked pixels are excluded AND the
                    output is re-masked to -9999 at those locations
    --sentinel VAL  Masked-pixel sentinel value (default: -9999.0)
    --verbose       Print per-file statistics and the averaged-map RMS
"""

import argparse
import os
import sys
import numpy as np


SENTINEL_DEFAULT = -9999.0


def infer_ndim(filename):
    """Infer N from the number of non-blank lines in a flat file."""
    with open(filename) as f:
        n_lines = sum(1 for line in f if line.strip() and not line.startswith('#'))
    n = int(round(np.sqrt(n_lines)))
    if n * n != n_lines:
        raise ValueError(
            f"{filename}: {n_lines} values is not a perfect square — "
            f"cannot infer grid size. Use --ndim to specify explicitly."
        )
    return n


def load_map(filename, ndim, sentinel):
    """
    Load a flat N*N map file.
    Returns a (N, N) float64 array with masked pixels set to NaN.
    """
    data = np.loadtxt(filename)
    if data.size != ndim * ndim:
        raise ValueError(
            f"{filename}: expected {ndim*ndim} values for {ndim}×{ndim} grid, "
            f"got {data.size}"
        )
    arr = data.reshape((ndim, ndim)).astype(np.float64)
    arr[np.abs(arr - sentinel) < 1.0] = np.nan   # sentinel → NaN
    return arr


def load_mask(mask_file, ndim):
    """Load a mask file (1=valid, 0=masked). Returns boolean array."""
    m = np.loadtxt(mask_file)
    if m.size != ndim * ndim:
        raise ValueError(
            f"Mask file {mask_file}: size {m.size} does not match {ndim}×{ndim} grid"
        )
    return m.reshape((ndim, ndim)).astype(bool)


def rms_valid(arr):
    """RMS about the mean over finite (non-NaN) pixels."""
    v = arr[np.isfinite(arr)]
    if v.size == 0:
        return np.nan
    return float(np.sqrt(np.mean((v - np.mean(v)) ** 2)))


def main():
    parser = argparse.ArgumentParser(
        description='Average multiple holography Epr.dat map files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('files', nargs='+',
                        help='Input map files (flat N*N single-column ASCII)')
    parser.add_argument('-o', '--output', default='mean_Epr.dat',
                        help='Output filename (default: mean_Epr.dat)')
    parser.add_argument('--ndim', type=int, default=None,
                        help='Grid size N (auto-detected from first file if omitted)')
    parser.add_argument('--mask', default=None,
                        help='Optional mask file (1=valid, 0=masked)')
    parser.add_argument('--sentinel', type=float, default=SENTINEL_DEFAULT,
                        help=f'Masked-pixel sentinel value (default: {SENTINEL_DEFAULT})')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print per-file and summary statistics')
    args = parser.parse_args()

    # Check input files exist
    missing = [f for f in args.files if not os.path.exists(f)]
    if missing:
        for f in missing:
            print(f"Error: file not found: {f}")
        sys.exit(1)

    n_files = len(args.files)
    if n_files < 2:
        print("Warning: only one input file — output will be identical to input.")

    # Determine grid size
    ndim = args.ndim
    if ndim is None:
        ndim = infer_ndim(args.files[0])
        if args.verbose:
            print(f"Auto-detected grid size: {ndim}×{ndim}")
    else:
        if args.verbose:
            print(f"Using specified grid size: {ndim}×{ndim}")

    # Load optional mask
    ext_mask = None
    if args.mask:
        ext_mask = load_mask(args.mask, ndim)
        if args.verbose:
            n_valid = ext_mask.sum()
            print(f"Mask: {n_valid} valid pixels out of {ndim*ndim} "
                  f"({100*n_valid/ndim/ndim:.1f}%)")

    # Load all maps into a 3-D stack (n_files, ndim, ndim)
    print(f"Loading {n_files} map(s)...")
    stack = np.full((n_files, ndim, ndim), np.nan)

    for i, fname in enumerate(args.files):
        arr = load_map(fname, ndim, args.sentinel)
        if ext_mask is not None:
            arr[~ext_mask] = np.nan
        stack[i] = arr
        if args.verbose:
            r = rms_valid(arr)
            n_valid = np.sum(np.isfinite(arr))
            print(f"  [{i+1}/{n_files}] {os.path.basename(fname):45s}  "
                  f"RMS = {r:.3f}   valid pixels = {n_valid}")

    # Compute NaN-aware mean
    mean_map = np.nanmean(stack, axis=0)

    # Pixels that were masked in ALL input maps remain NaN → restore sentinel
    all_masked = np.all(~np.isfinite(stack), axis=0)
    mean_map[all_masked] = args.sentinel

    # Also apply external mask to output
    if ext_mask is not None:
        mean_map[~ext_mask] = args.sentinel

    # Statistics
    finite_mean = mean_map[mean_map != args.sentinel]
    rms_mean = float(np.sqrt(np.mean((finite_mean - finite_mean.mean()) ** 2)))
    n_contrib = np.sum(~np.isnan(stack), axis=0)   # how many maps contributed to each pixel

    print(f"\nAveraged map:")
    print(f"  Grid size:          {ndim}×{ndim}")
    print(f"  Input maps:         {n_files}")
    print(f"  Valid output pixels:{np.sum(mean_map != args.sentinel)}")
    print(f"  Mean pixel contrib: {n_contrib[n_contrib > 0].mean():.2f} maps/pixel")
    print(f"  Output RMS:         {rms_mean:.4f}")

    # Write output
    np.savetxt(args.output, mean_map.flatten(), fmt='%15.8e')
    print(f"  Saved: {args.output}")


if __name__ == '__main__':
    main()
