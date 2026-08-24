#!/usr/bin/env python3
"""
plot_beam_maps.py

Plot the regridded beam amplitude and phase BEFORE the Fourier transform,
in the same style as unwrap2d.py's comparison figure (no panel overlays).

The input files are the flat N×N single-column files written by preprocess.py:
  ampout.dat  (or the filename listed as 'inamp'  in withphase_aber.prm)
  phaseout.dat (or the filename listed as 'inphase' in withphase_aber.prm)

Alternatively, pass --rgin to read rgin.dat (4-col: j k amp phase) directly
from the regridder output, before preprocess corrections are applied.

Usage:
    python plot_beam_maps.py                         # reads prm, auto-detects N
    python plot_beam_maps.py --prm withphase_aber.prm --output beam_maps.pdf
    python plot_beam_maps.py --amp ampout.dat --phase phaseout.dat -N 128
    python plot_beam_maps.py --rgin rgin.dat -N 128 --output beam_maps_raw.pdf
"""

import argparse
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


# ------------------------------------------------------------------ helpers

def _prm_value(line):
    """Extract value token from a Holis .prm line (col-49 first, then last token)."""
    val = line[49:].strip() if len(line) >= 50 else ""
    if not val:
        toks = line.split()
        val = toks[-1] if toks else ""
    return val


def read_prm(prm_file):
    """Read the minimal fields we need from withphase_aber.prm."""
    params = {}
    if not os.path.exists(prm_file):
        return params
    with open(prm_file) as f:
        for raw in f:
            line = raw.rstrip('\n')
            low = line.lower().strip()
            if low.startswith('!') or not low:
                continue
            try:
                if 'far field amplitude file name' in low:
                    params['inamp'] = _prm_value(line)
                elif 'far field phase file name' in low:
                    params['inphase'] = _prm_value(line)
                elif 'size n of the n by n' in low:
                    params['ndim'] = int(_prm_value(line))
                elif 'sampling interval' in low and 'far field' in low:
                    params['samp_itvl'] = float(_prm_value(line))
                elif 'observing frequency' in low:
                    params['freq'] = float(_prm_value(line))
                elif 'diameter of the primary' in low and 'secondary' not in low:
                    params['dprim'] = float(_prm_value(line))
            except ValueError:
                continue
    return params


def load_flat(filename, ndim):
    """Load a flat N*N single-column file and reshape to (N, N)."""
    data = np.loadtxt(filename)
    if data.size != ndim * ndim:
        raise ValueError(
            f"{filename}: expected {ndim*ndim} values for {ndim}×{ndim} grid, "
            f"got {data.size}"
        )
    return data.reshape((ndim, ndim))


def load_rgin(filename, ndim):
    """
    Load rgin.dat (4-col: j k amp phase) and scatter into N×N grids.
    Missing cells are NaN.
    """
    data = np.loadtxt(filename)
    amp_grid   = np.full((ndim, ndim), np.nan)
    phase_grid = np.full((ndim, ndim), np.nan)
    j = data[:, 0].astype(int)
    k = data[:, 1].astype(int)
    amp_grid[j, k]   = data[:, 2]
    phase_grid[j, k] = data[:, 3]
    return amp_grid, phase_grid


def infer_ndim(filename):
    """Guess N from the number of lines in a flat file."""
    with open(filename) as f:
        n_lines = sum(1 for line in f if line.strip() and not line.startswith('#'))
    n = int(round(np.sqrt(n_lines)))
    if n * n != n_lines:
        raise ValueError(
            f"Cannot infer grid size from {filename}: {n_lines} lines is not a perfect square."
        )
    return n


# ------------------------------------------------------------------ main plot

def make_figure(amp, phase, ndim, samp_itvl=None, prefix='', output=None,
                phase_unit='rad', amp_label='Amplitude', cmap_phase='RdBu_r',
                vmin_phase=None, vmax_phase=None, no_interactive=False):
    """
    Create side-by-side amplitude | phase figure matching unwrap2d style.

    Parameters
    ----------
    amp         : (N, N) array, amplitude (any units)
    phase       : (N, N) array, phase (rad or deg, NaN = masked)
    ndim        : grid size
    samp_itvl   : sampling interval in arcsec (for axis labelling)
    prefix      : dataset name for the title
    output      : output filename (.pdf or .png); None = display interactively
    phase_unit  : 'rad' or 'deg'
    """

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))

    # ---- amplitude ----
    ax = axes[0]
    amp_plot = np.where(amp == 0, np.nan, amp)       # zeros -> masked
    vmax_amp = np.nanpercentile(amp_plot, 99.5)
    im_amp = ax.imshow(amp_plot, origin='lower', cmap='viridis',
                       vmin=0, vmax=vmax_amp, aspect='equal')
    cb = fig.colorbar(im_amp, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(amp_label, fontsize=11)
    ax.set_title('Beam Amplitude', fontsize=13)

    # ---- phase ----
    ax = axes[1]
    ph_plot = phase.copy().astype(float)
    ph_plot[ph_plot == -9999] = np.nan               # Holis mask sentinel

    # Symmetric color scale about zero
    if vmin_phase is None or vmax_phase is None:
        pmax = np.nanpercentile(np.abs(ph_plot[np.isfinite(ph_plot)]), 99)
        vmin_phase = -pmax
        vmax_phase =  pmax

    im_ph = ax.imshow(ph_plot, origin='lower', cmap=cmap_phase,
                      vmin=vmin_phase, vmax=vmax_phase, aspect='equal')
    cb = fig.colorbar(im_ph, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(f'Phase ({phase_unit})', fontsize=11)
    ax.set_title('Beam Phase', fontsize=13)

    # Common axis labels
    if samp_itvl is not None:
        half = ndim / 2 * samp_itvl
        tick_pos  = np.linspace(0, ndim - 1, 5)
        tick_vals = np.linspace(-half, half, 5)
        for ax_ in axes:
            ax_.set_xticks(tick_pos)
            ax_.set_xticklabels([f'{v:.1f}' for v in tick_vals], fontsize=9)
            ax_.set_yticks(tick_pos)
            ax_.set_yticklabels([f'{v:.1f}' for v in tick_vals], fontsize=9)
            ax_.set_xlabel('Azimuth offset (arcsec)', fontsize=10)
            ax_.set_ylabel('Elevation offset (arcsec)', fontsize=10)
    else:
        for ax_ in axes:
            ax_.set_xlabel('Pixel (azimuth)', fontsize=10)
            ax_.set_ylabel('Pixel (elevation)', fontsize=10)

    # Overall title
    title = f'Regridded beam maps — {ndim}×{ndim} grid'
    if prefix:
        title = f'{prefix}\n{title}'
    fig.suptitle(title, fontsize=12, y=1.01)

    plt.tight_layout()

    if output:
        fmt = 'png' if output.lower().endswith('.png') else 'pdf'
        plt.savefig(output, format=fmt, dpi=150, bbox_inches='tight')
        print(f"Saved: {output}")

    if not no_interactive:
        plt.show()

    plt.close()


# ------------------------------------------------------------------ CLI

def main():
    parser = argparse.ArgumentParser(
        description='Plot pre-FFT beam amplitude and phase maps.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--prm', default='withphase_aber.prm',
                        help='Parameter file (default: withphase_aber.prm)')
    parser.add_argument('--amp', default=None,
                        help='Amplitude file (overrides prm inamp)')
    parser.add_argument('--phase', default=None,
                        help='Phase file (overrides prm inphase)')
    parser.add_argument('--rgin', default=None,
                        help='Load from rgin.dat (4-col) instead of flat files')
    parser.add_argument('-N', '--ndim', type=int, default=None,
                        help='Grid size N (auto-detected from file if omitted)')
    parser.add_argument('--prefix', default='',
                        help='Dataset name for plot title')
    parser.add_argument('--output', '-o', default=None,
                        help='Output file (.pdf or .png). Default: interactive display.')
    parser.add_argument('--phase-deg', action='store_true',
                        help='Phase file is in degrees (default: radians)')
    parser.add_argument('--vmin-phase', type=float, default=None,
                        help='Color scale minimum for phase panel')
    parser.add_argument('--vmax-phase', type=float, default=None,
                        help='Color scale maximum for phase panel')
    parser.add_argument('--no-interactive', action='store_true',
                        help='Suppress interactive display (use with --output)')
    args = parser.parse_args()

    # Read .prm for defaults
    prm = read_prm(args.prm)

    amp_file   = args.amp   or prm.get('inamp',   'ampout.dat')
    phase_file = args.phase or prm.get('inphase', 'phaseout.dat')
    ndim       = args.ndim  or prm.get('ndim')
    samp_itvl  = prm.get('samp_itvl')
    prefix     = args.prefix or prm.get('prefix', '')

    if args.rgin:
        # --- load from rgin.dat ---
        if ndim is None:
            # rgin.dat: deduce ndim from max j or k index
            data = np.loadtxt(args.rgin)
            ndim = int(max(data[:, 0].max(), data[:, 1].max())) + 1
            print(f"  Inferred grid size from rgin.dat: {ndim}×{ndim}")
        print(f"  Loading rgin.dat: {args.rgin}")
        amp, phase = load_rgin(args.rgin, ndim)
        phase_unit = 'deg' if args.phase_deg else 'rad'
        amp_label  = 'Amplitude (arb.)'
    else:
        # --- load from flat preprocess output files ---
        if not os.path.exists(amp_file):
            print(f"Error: amplitude file not found: {amp_file}")
            sys.exit(1)
        if not os.path.exists(phase_file):
            print(f"Error: phase file not found: {phase_file}")
            sys.exit(1)

        if ndim is None:
            ndim = infer_ndim(amp_file)
            print(f"  Inferred grid size from {amp_file}: {ndim}×{ndim}")

        print(f"  Loading amplitude: {amp_file}")
        print(f"  Loading phase:     {phase_file}")
        amp   = load_flat(amp_file,   ndim)
        phase = load_flat(phase_file, ndim)
        phase_unit = 'deg' if args.phase_deg else 'rad'
        amp_label  = 'Amplitude (arb.)'

    print(f"  Grid: {ndim}×{ndim},  phase unit: {phase_unit}")
    if samp_itvl:
        print(f"  Sampling interval: {samp_itvl} arcsec → FOV ≈ {ndim*samp_itvl:.1f} arcsec")

    make_figure(
        amp, phase, ndim,
        samp_itvl=samp_itvl,
        prefix=prefix,
        output=args.output,
        phase_unit=phase_unit,
        vmin_phase=args.vmin_phase,
        vmax_phase=args.vmax_phase,
        no_interactive=args.no_interactive,
    )


if __name__ == '__main__':
    main()
