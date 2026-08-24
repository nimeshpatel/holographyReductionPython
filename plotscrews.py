#!/usr/bin/env python3
"""
plotscrews.py
=============

Nimesh Patel

Originally written in September 2025 for photogrammetry maps. Revised on 
20 May 2026 to allow 128x128 holography maps. The Epr.dat files from holography
should first be processed by panel_fit.py to produce the adjustments table in
the same format as photogrammetry.

Plot dish-panel adjuster-screw corrections (Dz, microns) on a polar map.

Supports two input formats, selected with --mode:

  photogrammetry   (DEFAULT)
      ASCII residuals file with lines of the form
          P<m>_<n>   Dx   Dy   Dz   Dtot
      where m is the screw-ring index (1..24) and n is the entry around
      that ring. Dz is given in METRES and the screw correction shown on
      the plot is  dz_um = -Dz * 1000  (sign-flipped, per existing
      convention so that red = move panel toward subreflector).

  holography
      "correction.tabl"-style file with header lines and rows of the form
          <letter>_<n>   v1   v2   v3   v4   v5
      where <letter> is a..h (one letter per panel BAND) and n is the
      panel number within that band. The 5 numeric columns are the
      adjustments for the four CORNER screws and the CENTER screw of the
      panel, already in MICRONS, with the sign convention "positive =
      panel motion toward the subreflector" -- which matches the
      red = positive convention used by this plot, so NO sign flip is
      applied for holography by default (override with --flip-sign).

Both modes produce identical-looking output: a full-dish PDF plus four
per-quadrant zoom PDFs, with each screw labelled in microns.

Holography column / corner conventions
--------------------------------------
The angles stored in panelfit/screws.table are NOT in the same frame as
the dish-surface map. The GLT panelfit Fortran code transforms them via

    phi_physical = pi - theta_screws_table - offset

before placing each screw on the dish surface (see panelfit.f line 207
and search.f line 31). With the default offset=0 this is a reflection
about the +Y axis: x_displayed = -x_screws_table, y unchanged.

Verified empirically by sampling panelfit/Epr.dat at each screw position
in the four possible orientations and computing the correlation between
correction.tabl values and Epr.dat at that screw:
    screws-table direct                  r = -0.53
    rows inverted (image-origin upper)    r = -0.27
    x flipped (Fortran phi=pi-theta)      r = -0.80   <-- correct
    180-deg rotation                      r = -0.27

So in the dish-surface display frame, where the user reads the map:

   col 1 = inner corner on the CW  side of the panel  (iCW in display)
   col 2 = inner corner on the CCW side of the panel  (iCCW in display)
   col 3 = outer corner on the CW  side of the panel  (oCW in display)
   col 4 = outer corner on the CCW side of the panel  (oCCW in display)
   col 5 = center screw

The 5 letter rings map to the 8 panel bands in order (a innermost,
h outermost), with panel counts 12, 12, 24, 24, 48, 48, 48, 48.

In the display frame: panel n=1 has its leading edge on the -X axis
(theta=180 degrees) and numbering proceeds CLOCKWISE around the dish.

If a future file uses different conventions, override any of these:
   --holo-corner-order  comma-separated permutation of
                        {iCW, iCCW, oCW, oCCW} for the first four columns
                        (default: "iCW,iCCW,oCW,oCCW")
   --holo-col5          what column 5 represents:  "center" or "mean"
                        (default: center; "mean" hides it from the plot)
   --holo-start-deg     angle (deg, measured from +X CCW) of panel-1's
                        LEADING edge (default 180 = -X axis, matching the
                        Fortran panelfit transformation)
   --holo-ccw           numbering goes counter-clockwise instead of CW
                        (default CW in the display frame)
   --flip-sign          negate every value before plotting
"""

import argparse
import math
import os
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

# -------------------- Fixed panel-boundary radii (mm) --------------------
BASE_RADII = [600, 1265, 1820, 2605, 3220, 4040, 4780, 5435, 6000]

# Photogrammetry m-index conventions
CENTER_MS  = [2, 5, 8, 11, 14, 17, 20, 23]
EDGE_PAIRS = [(1, 3), (4, 6), (7, 9), (10, 12),
              (13, 15), (16, 18), (19, 21), (22, 24)]

# Holography ring letters -> band index (0..7)
HOLO_RING_ORDER = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

# Defaults
DEFAULT_DELTA_R_MM      = 75.0
DEFAULT_DELTA_DEG       = 1.5     # legacy fixed-angle inset (degrees)
DEFAULT_TANG_INSET_MM   = 50.0    # default tangential inset (mm of arc length)
DEFAULT_FONTSIZE        = 4
DEFAULT_FAINT_THRESH_UM = 21.0
DEFAULT_FAINT_ALPHA     = 0.15
DEFAULT_NORMAL_ALPHA    = 1.0

PMN_RE   = re.compile(r'^\s*P\s*(\d+)\s*[_]\s*(\d+)\b', re.IGNORECASE)
HOLO_RE  = re.compile(r'^\s*([a-hA-H])\s*_\s*(\d+)\b')


# ============================================================
#  Photogrammetry reader  (unchanged interface)
# ============================================================
def read_photogrammetry_table(path):
    """Read 'P<m>_<n>  Dx Dy Dz ...' lines. Returns (per_ring, nmax)."""
    per_ring = defaultdict(list)
    nmax = defaultdict(int)
    with open(path, 'r') as f:
        for line in f:
            m = PMN_RE.match(line)
            if not m:
                continue
            m_ring = int(m.group(1))
            n_idx  = int(m.group(2))
            rest   = line[m.end():].strip()
            if not rest:
                continue
            nums = []
            for tok in rest.split():
                try:
                    nums.append(float(tok))
                except ValueError:
                    break
            if len(nums) < 3:  # need at least Dx Dy Dz
                continue
            dz_mm = float(nums[2])
            per_ring[m_ring].append((n_idx, dz_mm))
            if n_idx > nmax[m_ring]:
                nmax[m_ring] = n_idx
    for m_ring in per_ring:
        per_ring[m_ring].sort(key=lambda t: t[0])
    return per_ring, nmax


def band_bounds():
    return BASE_RADII[:-1], BASE_RADII[1:]


def infer_panels_per_band_photogrammetry(per_ring, nmax):
    npanels = []
    for m_center in CENTER_MS:
        N = nmax.get(m_center, 0)
        if N <= 0:
            raise SystemExit(f"Missing/zero entries on center ring m={m_center}.")
        if len(per_ring.get(m_center, [])) != N:
            raise SystemExit(f"Center ring m={m_center} count mismatch (len!=N).")
        npanels.append(N)
    return npanels


# ============================================================
#  Holography reader
# ============================================================
def read_holography_table(path):
    """
    Read a 'correction.tabl'-style file.

    Returns
    -------
    holo : dict
        { letter (lower-case) : list of (panel_n, [v1, v2, v3, v4, v5]) sorted by n }
    nmax : dict
        { letter : max panel_n encountered }
    """
    holo = defaultdict(list)
    nmax = defaultdict(int)
    with open(path, 'r') as f:
        for line in f:
            m = HOLO_RE.match(line)
            if not m:
                continue
            letter = m.group(1).lower()
            n_idx  = int(m.group(2))
            rest   = line[m.end():].strip()
            if not rest:
                continue
            # Fortran sometimes writes "1.23E-002" -- normal float() handles it.
            nums = []
            for tok in rest.split():
                try:
                    nums.append(float(tok))
                except ValueError:
                    break
            if len(nums) < 5:
                continue  # need 5 columns
            holo[letter].append((n_idx, nums[:5]))
            if n_idx > nmax[letter]:
                nmax[letter] = n_idx
    for letter in holo:
        holo[letter].sort(key=lambda t: t[0])
    return holo, nmax


def infer_panels_per_band_holography(holo, nmax):
    npanels = []
    for letter in HOLO_RING_ORDER:
        N = nmax.get(letter, 0)
        if N <= 0:
            raise SystemExit(f"Missing/zero entries on holography ring '{letter}'.")
        if len(holo.get(letter, [])) != N:
            raise SystemExit(
                f"Holography ring '{letter}' count mismatch "
                f"(have {len(holo.get(letter, []))} rows, max n={N}).")
        npanels.append(N)
    return npanels


# ============================================================
#  Common annotation / drawing helpers
# ============================================================
def draw_panel_boundaries(ax, rmin, rmax, npanels, clip_on=False):
    for k in range(len(npanels)):
        rin, rout, n = rmin[k], rmax[k], npanels[k]
        for j in range(n):
            th1 = 2 * math.pi * j / n
            th2 = 2 * math.pi * (j + 1) / n
            poly = np.array([
                [rin  * math.cos(th1), rin  * math.sin(th1)],
                [rout * math.cos(th1), rout * math.sin(th1)],
                [rout * math.cos(th2), rout * math.sin(th2)],
                [rin  * math.cos(th2), rin  * math.sin(th2)],
            ])
            patch = Polygon(poly, closed=True, fill=False,
                            edgecolor="grey", linewidth=0.4)
            patch.set_clip_on(clip_on)
            ax.add_patch(patch)


def first_panel_clockwise_from_top(n_panels):
    """Return j (0-based) of the panel whose mid-angle is closest to +Y top,
    using the photogrammetry/dish convention of clockwise-from-top numbering."""
    width = 2.0 * math.pi / n_panels
    th_mid = [(j + 0.5) * width for j in range(n_panels)]
    top = math.pi / 2.0
    deltas = [(top - th) % (2.0 * math.pi) for th in th_mid]
    return int(np.argmin(deltas))


def theta_from_n(n_idx, n_on_ring):
    # Photogrammetry edge-ring convention:
    # Top (+Y) is pi/2; increasing n moves CLOCKWISE (decreases theta).
    return math.pi / 2.0 - 2.0 * math.pi * ((n_idx - 1) / float(n_on_ring))


def faint_alpha_for(dz_um, faint_thresh_um, faint_alpha, normal_alpha):
    return faint_alpha if abs(dz_um) <= faint_thresh_um else normal_alpha


def _label(ax, r, th, dz_um, fontsize, faint_thresh_um,
           faint_alpha, normal_alpha, clip_on):
    alpha = faint_alpha_for(dz_um, faint_thresh_um, faint_alpha, normal_alpha)
    ax.text(r * math.cos(th), r * math.sin(th),
            f"{dz_um:+.0f}",
            color=("red" if dz_um > 0 else "blue"),
            fontsize=fontsize, ha="center", va="center", alpha=alpha,
            clip_on=clip_on)


# ----- Photogrammetry annotation (panel-aware placement) --------------
# n-index convention assumed (this matches the photogrammetry data layout
# AND the output of panel_fit.py):
#   * Centre ring (m = 2, 5, 8, ...): Np entries, n = 1..Np. n = 1 is the
#     panel whose centre sits closest to the +Y axis; n increases CW.
#   * Edge rings  (m = 1, 3, 4, 6, ...): 2*Np entries. n = 2j-1 is panel
#     j's CCW corner (high-theta side) and n = 2j is panel j's CW corner.
#
# Label placement: corner labels sit at the actual screw positions in
# each band, computed from delta_r_mm (radial inset) and either
# tang_inset_mm (tangential inset; converts to an angular offset that
# adapts to the band radius) or delta_th (fixed angular inset, if
# tang_inset_mm is None).
def annotate_photogrammetry(ax, per_ring, nmax, npanels,
                            delta_r_mm, delta_deg, fontsize,
                            faint_thresh_um, faint_alpha, normal_alpha,
                            clip_on=False, tang_inset_mm=None):
    rmin, rmax = band_bounds()
    delta_th_fixed = math.radians(delta_deg)

    for k in range(8):
        rin, rout = rmin[k], rmax[k]
        Np        = npanels[k]
        panel_width = 2.0 * math.pi / Np
        m_in, m_out = EDGE_PAIRS[k]
        m_ctr       = CENTER_MS[k]

        # Radial positions for the three label rings
        r_inner = rin + delta_r_mm
        r_outer = rout - delta_r_mm
        r_centre = 0.5 * (rin + rout)

        # Per-band angular insets, in radians:
        #   - tangential mode: arc length / radius
        #   - fixed-angle  mode: delta_deg
        if tang_inset_mm is not None:
            delta_th_inner = tang_inset_mm / r_inner
            delta_th_outer = tang_inset_mm / r_outer
        else:
            delta_th_inner = delta_th_fixed
            delta_th_outer = delta_th_fixed

        def panel_center_th(j_1based):
            # Panel j (1-based, CW from +Y top): centre at
            #   pi/2 - (j - 0.5) * panel_width
            return math.pi / 2.0 - (j_1based - 0.5) * panel_width

        def corner_th(j_1based, is_ccw, delta_th):
            th_c = panel_center_th(j_1based)
            inset = panel_width / 2.0 - delta_th
            return th_c + inset if is_ccw else th_c - inset

        # CENTER ring (1 entry per panel)
        if m_ctr in per_ring and nmax.get(m_ctr, 0) == Np:
            for (n, dz_mm) in per_ring[m_ctr]:
                th_mid = panel_center_th(n)
                dz_um = -dz_mm * 1000.0
                _label(ax, r_centre, th_mid, dz_um, fontsize,
                       faint_thresh_um, faint_alpha, normal_alpha, clip_on)

        # INNER edge ring (2 entries per panel)
        if m_in in per_ring:
            for (n, dz_mm) in per_ring[m_in]:
                j = (n + 1) // 2          # 1-based panel
                th = corner_th(j, is_ccw=(n % 2 == 1),
                               delta_th=delta_th_inner)
                dz_um = -dz_mm * 1000.0
                _label(ax, r_inner, th, dz_um, fontsize,
                       faint_thresh_um, faint_alpha, normal_alpha, clip_on)

        # OUTER edge ring (2 entries per panel)
        if m_out in per_ring:
            for (n, dz_mm) in per_ring[m_out]:
                j = (n + 1) // 2
                th = corner_th(j, is_ccw=(n % 2 == 1),
                               delta_th=delta_th_outer)
                dz_um = -dz_mm * 1000.0
                _label(ax, r_outer, th, dz_um, fontsize,
                       faint_thresh_um, faint_alpha, normal_alpha, clip_on)


# ----- Holography annotation ------------------------------------------
def _holo_panel_center_angle(n_idx, Np, start_deg, ccw):
    """
    Centre angle (radians) of panel n_idx in a band of Np panels.

    start_deg : angle in degrees (from +X axis, CCW positive) of the
                LEADING edge of panel n=1.  The centre of panel n=1 sits
                half a panel width into the band, in the direction of
                numbering.
    ccw       : if True, increasing n goes CCW; if False, increasing n
                goes CW.
    """
    start_rad = math.radians(start_deg)
    width     = 2.0 * math.pi / Np
    # Panel 1 centre sits half a panel-width past the leading edge,
    # in the direction of numbering.
    if ccw:
        return start_rad + (n_idx - 1 + 0.5) * width
    else:
        return start_rad - (n_idx - 1 + 0.5) * width


def annotate_holography(ax, holo, npanels,
                        delta_r_mm, delta_deg, fontsize,
                        faint_thresh_um, faint_alpha, normal_alpha,
                        corner_order, col5_role, start_deg, ccw,
                        flip_sign, clip_on=False, tang_inset_mm=None):
    """
    corner_order : 4-tuple (or list) of strings from {'iCW','iCCW','oCW','oCCW'}
                   giving the role of columns 1..4 respectively.
    col5_role    : 'center' or 'mean'
    tang_inset_mm: if given, treat the tangential inset of corner labels as a
                   FIXED ARC LENGTH (mm) so the labels track the actual screw
                   positions in every band. If None, the fixed angular inset
                   delta_deg is used at every radius (matches old behaviour).
    """
    rmin, rmax = band_bounds()
    delta_th_fixed = math.radians(delta_deg)
    sign           = -1.0 if flip_sign else 1.0

    for k, letter in enumerate(HOLO_RING_ORDER):
        if letter not in holo:
            continue
        rin, rout = rmin[k], rmax[k]
        Np = npanels[k]
        width = 2.0 * math.pi / Np

        r_in_lbl  = rin  + delta_r_mm
        r_out_lbl = rout - delta_r_mm

        # Per-band angular insets, in radians, evaluated at the actual
        # radius where the labels are drawn.
        if tang_inset_mm is not None:
            delta_th_inner = tang_inset_mm / r_in_lbl
            delta_th_outer = tang_inset_mm / r_out_lbl
        else:
            delta_th_inner = delta_th_fixed
            delta_th_outer = delta_th_fixed

        for (n_idx, vals) in holo[letter]:
            th_mid = _holo_panel_center_angle(n_idx, Np, start_deg, ccw)

            # Tangential offsets for "CW" and "CCW" corners relative to mid.
            # In screen terms with +Y up:
            #   CW corner of a panel = smaller theta side (because clockwise
            #     numbering reduces theta in our default convention)
            #   CCW corner = larger theta side
            # But to keep CW/CCW meaningful regardless of numbering direction,
            # define:
            #   CW  side: th_mid - (width/2 - delta_th)
            #   CCW side: th_mid + (width/2 - delta_th)
            # This way the labels always sit physically on the same side of
            # the panel.
            half_in  = width / 2.0 - delta_th_inner
            half_out = width / 2.0 - delta_th_outer
            th_cw_in   = th_mid - half_in
            th_ccw_in  = th_mid + half_in
            th_cw_out  = th_mid - half_out
            th_ccw_out = th_mid + half_out

            pos_map = {
                'iCW':  (r_in_lbl,  th_cw_in),
                'iCCW': (r_in_lbl,  th_ccw_in),
                'oCW':  (r_out_lbl, th_cw_out),
                'oCCW': (r_out_lbl, th_ccw_out),
            }

            for col_i, role in enumerate(corner_order):
                if role not in pos_map:
                    continue
                r_lbl, th_lbl = pos_map[role]
                dz_um = sign * float(vals[col_i])
                _label(ax, r_lbl, th_lbl, dz_um, fontsize,
                       faint_thresh_um, faint_alpha, normal_alpha, clip_on)

            # Column 5
            if col5_role == 'center':
                r_lbl  = 0.5 * (rin + rout)
                dz_um  = sign * float(vals[4])
                _label(ax, r_lbl, th_mid, dz_um, fontsize,
                       faint_thresh_um, faint_alpha, normal_alpha, clip_on)
            # 'mean' -> not drawn


# ============================================================
#  Plotting (quadrants + full)
# ============================================================
def _do_annotate(ax, ctx, fontsize, clip_on):
    """Dispatch to the right annotator based on the loaded mode context."""
    if ctx['mode'] == 'photogrammetry':
        annotate_photogrammetry(
            ax, ctx['per_ring'], ctx['nmax'], ctx['npanels'],
            ctx['delta_r_mm'], ctx['delta_deg'], fontsize,
            ctx['faint_threshold_um'], ctx['faint_alpha'], ctx['normal_alpha'],
            clip_on=clip_on,
            tang_inset_mm=ctx.get('tang_inset_mm'))
    else:  # holography
        annotate_holography(
            ax, ctx['holo'], ctx['npanels'],
            ctx['delta_r_mm'], ctx['delta_deg'], fontsize,
            ctx['faint_threshold_um'], ctx['faint_alpha'], ctx['normal_alpha'],
            ctx['holo_corner_order'], ctx['holo_col5'],
            ctx['holo_start_deg'], ctx['holo_ccw'], ctx['flip_sign'],
            clip_on=clip_on,
            tang_inset_mm=ctx.get('tang_inset_mm'))


def create_quadrant_plot(ctx, quadrant_num, input_filename):
    rmin, rmax = band_bounds()
    dish_radius = 6100

    fig, ax = plt.subplots(figsize=(8, 8))
    zoom_fontsize = ctx['fontsize'] * 2

    ax.set_aspect("equal", adjustable="box")

    if quadrant_num == 1:
        ax.set_xlim(0, dish_radius);            ax.set_ylim(0, dish_radius)
        quad_label = "Quadrant 1 (Top-Right)"
    elif quadrant_num == 2:
        ax.set_xlim(-dish_radius, 0);           ax.set_ylim(0, dish_radius)
        quad_label = "Quadrant 2 (Top-Left)"
    elif quadrant_num == 3:
        ax.set_xlim(-dish_radius, 0);           ax.set_ylim(-dish_radius, 0)
        quad_label = "Quadrant 3 (Bottom-Left)"
    elif quadrant_num == 4:
        ax.set_xlim(0, dish_radius);            ax.set_ylim(-dish_radius, 0)
        quad_label = "Quadrant 4 (Bottom-Right)"

    draw_panel_boundaries(ax, rmin, rmax, ctx['npanels'], clip_on=True)
    _do_annotate(ax, ctx, zoom_fontsize, clip_on=True)

    ax.set_xlabel("X [mm]", fontsize=12)
    ax.set_ylabel("Y [mm]", fontsize=12)
    ax.set_title(f"Adjustments from {input_filename}\nDz (microns) - {quad_label}",
                 fontsize=14)
    ax.grid(False)
    ax.tick_params(axis='both', which='major', labelsize=10)
    return fig


# ============================================================
#  CLI / main
# ============================================================
def _parse_corner_order(s):
    parts = [p.strip() for p in s.split(',')]
    if len(parts) != 4:
        raise argparse.ArgumentTypeError(
            "--holo-corner-order must list exactly 4 entries")
    allowed = {'iCW', 'iCCW', 'oCW', 'oCCW'}
    for p in parts:
        if p not in allowed:
            raise argparse.ArgumentTypeError(
                f"Bad corner code '{p}'. Use any of {sorted(allowed)}.")
    if len(set(parts)) != 4:
        raise argparse.ArgumentTypeError("Corner codes must be unique.")
    return parts


def main():
    ap = argparse.ArgumentParser(
        description=("Plot Dz (microns) screw adjustments from either "
                     "photogrammetry residuals or holography correction tables.")
    )
    ap.add_argument("input_file")
    ap.add_argument("--mode", choices=("photogrammetry", "holography"),
                    default="photogrammetry",
                    help="Which file format / interpretation to use. "
                         "Default: photogrammetry.")
    ap.add_argument("--output", help="Save full-dish plot to this PDF; "
                                     "four _Q1.._Q4 quadrant PDFs are also produced.")
    ap.add_argument("--fontsize", type=int, default=DEFAULT_FONTSIZE)
    ap.add_argument("--delta-r-mm", type=float, default=DEFAULT_DELTA_R_MM)
    ap.add_argument("--delta-deg", type=float, default=DEFAULT_DELTA_DEG,
                    help="Fixed tangential inset of corner labels, in DEGREES. "
                         "Used when --tangential-inset-mm is not given. "
                         "Becomes less accurate in inner bands because a fixed "
                         "angle = different arc length at each radius.")
    ap.add_argument("--tangential-inset-mm", type=float,
                    default=DEFAULT_TANG_INSET_MM,
                    help="Tangential inset of corner labels, in MILLIMETRES "
                         "of arc length, evaluated at each band's label "
                         "radius. Keeps labels stuck to the screws in EVERY "
                         "band (the screws themselves sit ~50 mm in from the "
                         "panel edge). Overrides --delta-deg. Pass 0 (or any "
                         "negative number) to fall back to --delta-deg.")
    ap.add_argument("--faint-threshold-um", type=float,
                    default=DEFAULT_FAINT_THRESH_UM,
                    help="Values with |Dz| <= this appear faint (default 21 um)")
    ap.add_argument("--faint-alpha", type=float, default=DEFAULT_FAINT_ALPHA,
                    help="Opacity for faint values (default 0.15)")
    ap.add_argument("--normal-alpha", type=float, default=DEFAULT_NORMAL_ALPHA,
                    help="Opacity for normal values (default 1.0)")
    ap.add_argument("--title", type=str, default="Adjuster Screw Dz (microns)")

    # Holography-specific overrides
    ap.add_argument("--holo-corner-order", type=_parse_corner_order,
                    default=['iCW', 'iCCW', 'oCW', 'oCCW'],
                    help="Comma-separated mapping of columns 1..4 to corner "
                         "positions in the DISPLAY frame. Use iCW/iCCW/oCW/oCCW "
                         "(inner/outer, clockwise/counter-clockwise). "
                         "Default: iCW,iCCW,oCW,oCCW (after the panelfit "
                         "phi=pi-theta mirror that swaps screws.table CCW <-> "
                         "display CW).")
    ap.add_argument("--holo-col5", choices=("center", "mean"), default="center",
                    help="Meaning of column 5. 'center' plots it at the panel "
                         "centre; 'mean' hides it. Default: center.")
    ap.add_argument("--holo-start-deg", type=float, default=180.0,
                    help="Angle (deg, CCW from +X) of the LEADING edge of "
                         "panel n=1 in the DISPLAY frame. Default 180 (-X axis), "
                         "matching the Fortran panelfit phi=pi-theta transform.")
    ap.add_argument("--holo-ccw", action="store_true",
                    help="Holography panel numbering goes CCW (default is CW "
                         "in the display frame).")
    ap.add_argument("--flip-sign", action="store_true",
                    help="Negate every value before plotting. Use if your "
                         "file's sign convention is opposite to "
                         "'positive = move toward subreflector'.")
    args = ap.parse_args()

    # -------- Sniff format & complain helpfully if --mode is wrong --------
    n_pmn  = 0
    n_holo = 0
    with open(args.input_file, 'r') as f:
        for line in f:
            if PMN_RE.match(line):
                n_pmn += 1
            elif HOLO_RE.match(line):
                n_holo += 1
    if args.mode == "photogrammetry" and n_pmn == 0 and n_holo > 0:
        raise SystemExit(
            f"{args.input_file} looks like a holography correction.tabl "
            f"(found {n_holo} '<letter>_<n>' lines, no 'P<m>_<n>' lines). "
            f"Re-run with --mode holography.")
    if args.mode == "holography" and n_holo == 0 and n_pmn > 0:
        raise SystemExit(
            f"{args.input_file} looks like a photogrammetry residuals file "
            f"(found {n_pmn} 'P<m>_<n>' lines, no '<letter>_<n>' lines). "
            f"Re-run without --mode holography (the default already handles "
            f"this format -- it's what panel_fit.py emits).")

    # -------- Load the right table --------
    if args.mode == "photogrammetry":
        per_ring, nmax = read_photogrammetry_table(args.input_file)
        npanels = infer_panels_per_band_photogrammetry(per_ring, nmax)
        ctx = dict(mode="photogrammetry",
                   per_ring=per_ring, nmax=nmax, npanels=npanels)
    else:
        holo, nmax = read_holography_table(args.input_file)
        npanels = infer_panels_per_band_holography(holo, nmax)
        ctx = dict(mode="holography",
                   holo=holo, nmax=nmax, npanels=npanels,
                   holo_corner_order=args.holo_corner_order,
                   holo_col5=args.holo_col5,
                   holo_start_deg=args.holo_start_deg,
                   holo_ccw=args.holo_ccw,
                   flip_sign=args.flip_sign)

    # tang_inset_mm of None  -> use --delta-deg
    # tang_inset_mm > 0      -> use that arc length
    tang_inset_mm = args.tangential_inset_mm
    if tang_inset_mm is not None and tang_inset_mm <= 0:
        tang_inset_mm = None

    # Common parameters carried in ctx
    ctx.update(
        delta_r_mm=args.delta_r_mm,
        delta_deg=args.delta_deg,
        tang_inset_mm=tang_inset_mm,
        fontsize=args.fontsize,
        faint_threshold_um=args.faint_threshold_um,
        faint_alpha=args.faint_alpha,
        normal_alpha=args.normal_alpha,
    )
    if 'flip_sign' not in ctx:
        ctx['flip_sign'] = args.flip_sign

    # -------- Full-dish plot --------
    rmin, rmax = band_bounds()
    fig, ax = plt.subplots(figsize=(8, 8))
    draw_panel_boundaries(ax, rmin, rmax, ctx['npanels'])
    _do_annotate(ax, ctx, args.fontsize, clip_on=False)

    dish_radius = 6100
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-dish_radius, dish_radius)
    ax.set_ylim(-dish_radius, dish_radius)
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")

    input_filename = os.path.basename(args.input_file)
    ax.set_title(f"Adjustments from {input_filename}, Dz (microns)")
    ax.grid(False)

    if args.output:
        plt.savefig(args.output, bbox_inches="tight")
        print(f"Saved {args.output}")

        base_name = os.path.splitext(args.output)[0]
        for q in range(1, 5):
            quad_fig = create_quadrant_plot(ctx, q, input_filename)
            quad_output = f"{base_name}_Q{q}.pdf"
            quad_fig.savefig(quad_output, bbox_inches="tight")
            print(f"Saved {quad_output}")
            plt.close(quad_fig)
    else:
        plt.show()


if __name__ == "__main__":
    main()
