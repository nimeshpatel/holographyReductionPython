#!/usr/bin/env python3
"""
panel_fit.py
============

Nimesh Patel
20 May 2026
This is a re-write of TK's original fortran program panelfit.f
which has been used for SMA and APEX antennas holography.
The panelfit.prm file is no longer needed in this python version, 
but all the relevant numbers are included below from that file, for GLT 12m
antenna.

Convert a holography Epr.dat phase map into a photogrammetry-style
residuals table that can be fed directly into plotscrews.py.

Inputs
------
* Epr.dat -- a 128x128 array (ASCII, one value per line) of residual
  phase, in radians.  Cells with no data are -9999.

The conversion:
  1. Read the 128x128 phase array.
  2. Convert phase to surface error in microns via
         err_um = phase * lambda/(4*pi) * 1e6
                  * sqrt(1 + (x^2 + y^2) / (4 * fprim^2))
     (same expression as the GLT panelfit Fortran code errorcal.f).
  3. For each panel of the dish, gather every data pixel within the
     panel's (r, theta) wedge and fit a tilted plane:
         err(x, y) ~ a + b*x + c*y
  4. Evaluate the fitted plane at five screw positions per panel:
     four corners (each inset 'inset_mm' from the panel edges) and the
     panel centre.
  5. Write out the result as P<m>_<n>  Dx  Dy  Dz  Dtot lines, matching
     the format produced by combine_photogrammetry_files.py.

Coordinate convention (matches gltDishMapSmoothed.py and the dishMap):
  * pixel spacing dx = dprim / rate / ndim = 12.0 m / 0.75 / 128 = 125 mm.
  * row 0 displays at the TOP of the image (+Y); row 127 at the bottom (-Y).
    Equivalently, y_phys = (ncent - row) * dx.
  * col 0 is at -X (left), col 127 at +X (right).

Photogrammetry n-numbering convention (matches plotscrews.py):
  * n = 1 sits at the top of the dish (+Y axis).
  * Numbering increases CLOCKWISE around the ring.
  * For centre rings (m = 2, 5, 8, ..., 23), there are Np entries and
    n = 1 is panel #1 (the panel whose centre sits closest to +Y).
  * For edge rings (m = 1, 3, 4, 6, ...), there are 2*Np entries.
    n = 2j-1 is panel j's CCW corner (high-theta side) and n = 2j is
    panel j's CW corner.

Output Dz sign:
  * Dz is given in mm, with the photogrammetry sign convention:
        Dz > 0  =>  panel sits ABOVE the ideal paraboloid.
  * plotscrews.py applies dz_displayed_um = -1000 * Dz_mm before
    plotting, so red labels mean "push panel toward subreflector"
    (positive correction).

Usage
-----
    python3 panel_fit.py panelfit/Epr.dat residuals_from_holo.txt
    python3 plotscrews.py residuals_from_holo.txt --output map.pdf
"""

import argparse
import math
import os
import sys
import numpy as np

# ----- telescope and panel geometry (panelplt.prm / panelfit.prm) -----
DPRIM_M  = 12.0                 # primary mirror diameter
FPRIM_M  =  4.8                 # primary mirror focal length
RATE     =  0.75                # Nyquist sampling rate used in errorcal
NDIM     = 128                  # phase data array size

# Panel band radii (mm). Matches plotscrews.py BASE_RADII so screws will
# fall in the same locations as the photogrammetry-mode plot expects.
BASE_RADII      = [600, 1265, 1820, 2605, 3220, 4040, 4780, 5435, 6000]
PANELS_PER_BAND = [12, 12, 24, 24, 48, 48, 48, 48]

# m-index convention used by plotscrews.py photogrammetry mode.
CENTER_MS  = [2, 5, 8, 11, 14, 17, 20, 23]
EDGE_PAIRS = [(1, 3), (4, 6), (7, 9), (10, 12),
              (13, 15), (16, 18), (19, 21), (22, 24)]

C_MPS = 2.99792458e8


# =====================================================================
#  IO
# =====================================================================
def load_epr(path):
    """Load Epr.dat as a 128x128 array of phase values (radians)."""
    vals = np.loadtxt(path)
    if vals.size != NDIM * NDIM:
        raise SystemExit(
            f"Expected {NDIM*NDIM} values in {path}, got {vals.size}.")
    return vals.reshape(NDIM, NDIM)


def phase_to_surface_um(phase_rad, x_mm, y_mm, freq_ghz, fprim_m):
    """Convert phase (radians) to surface error in microns, with the
    1/cos(grazing angle) geometric correction for the parabolic surface,
    exactly as in panelfit/errorcal.f line 28-29."""
    lam_m   = C_MPS / (freq_ghz * 1e9)
    fprim_mm = fprim_m * 1000.0
    geom = np.sqrt(1.0 + (x_mm**2 + y_mm**2) / (4.0 * fprim_mm**2))
    return phase_rad * (lam_m / (4.0 * math.pi)) * 1e6 * geom


# =====================================================================
#  Panel-by-panel plane fit
# =====================================================================
def fit_plane(x, y, z):
    """Least-squares fit of z = a + b*x + c*y.  Returns (a, b, c)."""
    A = np.column_stack([np.ones_like(x), x, y])
    coef, *_ = np.linalg.lstsq(A, z, rcond=None)
    return float(coef[0]), float(coef[1]), float(coef[2])


def panel_pixel_mask(r_grid, th_grid, rin, rout, th_min, th_max,
                     r_margin, ang_margin_deg):
    """Boolean mask of pixels inside the (r, theta) wedge of a panel,
    pulled in by `r_margin` mm radially and `ang_margin_deg` degrees on
    each tangential edge to avoid grabbing pixels from neighbours."""
    in_r = (r_grid >= rin + r_margin) & (r_grid <= rout - r_margin)
    # Handle wrap-around in theta:
    a = (th_min + ang_margin_deg) % 360.0
    b = (th_max - ang_margin_deg) % 360.0
    if a < b:
        in_th = (th_grid >= a) & (th_grid <= b)
    else:
        in_th = (th_grid >= a) | (th_grid <= b)
    return in_r & in_th


def screw_positions(rin, rout, th_c_deg, panel_width_deg,
                    radial_inset_mm, tangential_inset_mm):
    """Return the five screw (x, y) positions for one panel, plus the
    role of each (iCCW, iCW, oCCW, oCW, center)."""
    r_inner = rin  + radial_inset_mm
    r_outer = rout - radial_inset_mm
    r_centre = 0.5 * (rin + rout)

    # Tangential inset in degrees at the relevant radius
    ang_in  = math.degrees(tangential_inset_mm / r_inner)
    ang_out = math.degrees(tangential_inset_mm / r_outer)

    half = panel_width_deg / 2.0
    th_iCCW = th_c_deg + (half - ang_in)
    th_iCW  = th_c_deg - (half - ang_in)
    th_oCCW = th_c_deg + (half - ang_out)
    th_oCW  = th_c_deg - (half - ang_out)
    th_C    = th_c_deg

    def xy(r, th_deg):
        th = math.radians(th_deg)
        return r * math.cos(th), r * math.sin(th)

    return [
        ('iCCW', *xy(r_inner,  th_iCCW)),
        ('iCW',  *xy(r_inner,  th_iCW)),
        ('oCCW', *xy(r_outer,  th_oCCW)),
        ('oCW',  *xy(r_outer,  th_oCW)),
        ('C',    *xy(r_centre, th_C)),
    ]


# =====================================================================
#  Main
# =====================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Fit each dish panel to a tilted plane on the "
                    "Epr.dat phase map and emit a photogrammetry-style "
                    "residuals table consumable by plotscrews.py.")
    ap.add_argument("epr_file",
                    help="Holography phase file (128x128 ASCII, radians)")
    ap.add_argument("output_file",
                    help="Path to write the residuals table")
    ap.add_argument("--freq-ghz", type=float, default=94.5,
                    help="Holography frequency in GHz (default 94.5)")
    ap.add_argument("--radial-inset-mm", type=float, default=50.0,
                    help="Screw inset from inner/outer panel edges (mm). "
                         "Default 50 (matches panelplt.prm).")
    ap.add_argument("--tangential-inset-mm", type=float, default=50.0,
                    help="Screw inset from CCW/CW panel edges, measured "
                         "as arc length (mm). Default 50.")
    ap.add_argument("--fit-r-margin-mm", type=float, default=125.0,
                    help="Discard pixels within this many mm of the "
                         "inner/outer panel edges when fitting. "
                         "Default 125 (one pixel).")
    ap.add_argument("--fit-ang-margin-deg", type=float, default=0.25,
                    help="Discard pixels within this many degrees of "
                         "the CCW/CW panel edges when fitting (default "
                         "0.25 degrees).")
    ap.add_argument("--diagnostics", action="store_true",
                    help="Print per-panel fit statistics")
    args = ap.parse_args()

    # ---------- load & convert ----------
    phase = load_epr(args.epr_file)
    valid_mask = phase > -9000.0
    dx_mm = DPRIM_M * 1000.0 / (RATE * NDIM)
    ncent = (NDIM - 1) / 2.0   # array centre index (63.5 for ndim=128)

    rows, cols = np.indices(phase.shape)
    x_grid = (cols  - ncent) * dx_mm
    # Row 0 -> +Y (top), matching origin='upper' display convention.
    y_grid = (ncent - rows) * dx_mm
    r_grid  = np.sqrt(x_grid**2 + y_grid**2)
    th_grid = np.degrees(np.arctan2(y_grid, x_grid)) % 360.0

    lam_m = C_MPS / (args.freq_ghz * 1e9)
    err_um = phase_to_surface_um(phase, x_grid, y_grid,
                                 args.freq_ghz, FPRIM_M)
    # Cells with no data must not contaminate the fit:
    err_um = np.where(valid_mask, err_um, np.nan)

    print(f"Loaded {valid_mask.sum()}/{phase.size} valid pixels from "
          f"{args.epr_file}")
    print(f"Pixel size: dx = dprim/rate/ndim = {dx_mm:.3f} mm")
    print(f"Frequency:  {args.freq_ghz} GHz  ->  lambda = {lam_m*1e3:.4f} mm")
    print(f"Surface-error range over the dish:  "
          f"{np.nanmin(err_um):+.1f} to {np.nanmax(err_um):+.1f} microns")
    print()

    # ---------- fit each panel and emit one screw record each ----------
    lines = []
    diagnostics_rows = []

    for k in range(8):
        rin   = BASE_RADII[k]
        rout  = BASE_RADII[k + 1]
        Np    = PANELS_PER_BAND[k]
        pw    = 360.0 / Np
        m_in, m_out = EDGE_PAIRS[k]
        m_ctr       = CENTER_MS[k]

        for j in range(Np):
            # Panel index j is 0-based; the 1-based "panel number" is j+1.
            # CW numbering from +Y top: panel j has centre at
            #   th_c = 90 - (j + 0.5) * pw
            th_c = (90.0 - (j + 0.5) * pw) % 360.0
            th_min = (th_c - pw / 2.0) % 360.0
            th_max = (th_c + pw / 2.0) % 360.0

            mask = panel_pixel_mask(r_grid, th_grid, rin, rout,
                                    th_min, th_max,
                                    args.fit_r_margin_mm,
                                    args.fit_ang_margin_deg)
            mask &= valid_mask
            if mask.sum() < 4:
                if args.diagnostics:
                    print(f"  band {k+1} panel {j+1:>2}: too few pixels "
                          f"({mask.sum()}); skipping.")
                continue

            xp = x_grid[mask]
            yp = y_grid[mask]
            zp = err_um[mask]

            a, b, c = fit_plane(xp, yp, zp)
            if args.diagnostics:
                resid = zp - (a + b*xp + c*yp)
                diagnostics_rows.append(
                    (k+1, j+1, mask.sum(), zp.mean(), zp.std(), resid.std()))

            screws = screw_positions(
                rin, rout, th_c, pw,
                args.radial_inset_mm, args.tangential_inset_mm)

            # Map each screw to (m, n) and write the line.
            for role, xs, ys in screws:
                # Evaluate panel fit at this screw position
                z_um = a + b * xs + c * ys
                # Photogrammetry sign convention: Dz > 0  =>  panel above
                # ideal.  Our fit IS the (signed) surface error in microns,
                # so dz_mm = z_um / 1000.
                dz_mm = z_um / 1000.0
                dtot  = abs(dz_mm)

                # m, n indices
                if role == 'C':
                    m_idx = m_ctr
                    n_idx = j + 1
                elif role == 'iCCW':
                    m_idx = m_in
                    n_idx = 2 * (j + 1) - 1
                elif role == 'iCW':
                    m_idx = m_in
                    n_idx = 2 * (j + 1)
                elif role == 'oCCW':
                    m_idx = m_out
                    n_idx = 2 * (j + 1) - 1
                elif role == 'oCW':
                    m_idx = m_out
                    n_idx = 2 * (j + 1)
                else:
                    raise RuntimeError(f"Unknown screw role {role!r}")

                # The output file also carries x, y, z, but plotscrews.py
                # ignores them; we still emit physical x and y for debug
                # purposes. The third number per line MUST be Dz (mm).
                lines.append((m_idx, n_idx, xs, ys, dz_mm, dtot))

    # Sort lines by (m, n) for readability
    lines.sort(key=lambda t: (t[0], t[1]))

    with open(args.output_file, 'w') as f:
        for m_idx, n_idx, xs, ys, dz_mm, dtot in lines:
            # Format matches Paraboloid_Residuals*.txt: label + Dx Dy Dz Dtot
            # We use (Dx, Dy) = (0, 0); only Dz is consumed by plotscrews.py.
            label = f"P{m_idx}_{n_idx}"
            f.write(f"{label:<10}{0.0:12.3f}{0.0:12.3f}"
                    f"{dz_mm:12.6f}{dtot:12.6f}\n")
    print(f"Wrote {len(lines)} screw records to {args.output_file}")

    if args.diagnostics:
        print()
        print(f"{'band':>4} {'panel':>5} {'n_pix':>6} {'mean_um':>9} "
              f"{'std_um':>9} {'resid_std':>9}")
        for row in diagnostics_rows:
            print(f"{row[0]:>4} {row[1]:>5} {row[2]:>6} "
                  f"{row[3]:+9.2f} {row[4]:>9.2f} {row[5]:>9.2f}")


if __name__ == "__main__":
    main()
