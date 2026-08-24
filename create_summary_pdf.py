#!/usr/bin/env python3
"""
create_summary_pdf.py

Nimesh Patel

Create a summary PDF for holography data reduction results.
Combines text information with embedded PDF figures.

Usage:
    python create_summary_pdf.py <data_filename> <output_prefix> [--comment "comment string"]

Arguments:
    data_filename: Original input data filename
    output_prefix: Prefix used for output files (e.g., results/holoADC-AzEl-20250928_170604)
    --comment: Optional comment string to include in summary

Output:
    Creates results/<output_prefix>_summary.pdf
"""

import argparse
import os
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from PIL import Image
import numpy as np


def extract_rms_values(output_prefix):
    """
    Extract RMS surface error values for the summary page.

    Priority order:
      1. JSON sidecar written by glt_dish_map.py  (<prefix>_rms.json)
      2. Regex parse of the holis log file         (<prefix>.log)

    Returns a dict with keys 'standard', 'area_weighted', 'illum_weighted',
    'baars_illum_weighted', 'quantity', and the fallback key 'log_rms'
    (the raw value from the log, in microns).  Missing values are None.
    """
    import json

    result = {
        'standard': None,
        'area_weighted': None,
        'illum_weighted': None,
        'baars_illum_weighted': None,
        'quantity': 'surface error',
        'log_rms': None,
    }

    # --- 1. JSON sidecar (most reliable, written by glt_dish_map.py) ---
    sidecar = output_prefix + "_rms.json"
    if os.path.exists(sidecar):
        try:
            with open(sidecar) as f:
                data = json.load(f)
            result.update(data)
            print(f"  RMS read from {sidecar}")
            return result
        except Exception as e:
            print(f"  Warning: could not read {sidecar}: {e}")

    # --- 2. Fallback: parse holis log ---
    log_file = output_prefix + ".log"
    if not os.path.exists(log_file):
        print(f"  Warning: neither {sidecar} nor {log_file} found")
        return result

    rms_microns = None
    try:
        with open(log_file) as f:
            for line in f:
                m = re.search(r'RMS after phasefit.*\(([0-9.e+\-]+)\s+microns\)', line)
                if m:
                    rms_microns = float(m.group(1))   # keep overwriting -> last value wins
        if rms_microns is not None:
            result['log_rms'] = rms_microns
            result['standard'] = rms_microns
            print(f"  RMS read from {log_file}: {rms_microns:.2f} µm")
        else:
            print(f"  Warning: no RMS line found in {log_file}")
    except Exception as e:
        print(f"  Warning: could not parse {log_file}: {e}")

    return result


def create_summary_pdf(data_filename, output_prefix, comment=None):
    """
    Create a composite PDF summary page with text and embedded figures.

    Layout (portrait-ish 16×11):
      Row 0 (top ~22%): two side-by-side inset boxes
            LEFT  — filename, RMS, comment  (small fonts)
            RIGHT — aberration fit table
      Row 1 (bottom ~78%): illumination map | surface error map
    """

    # --- file paths ---
    log_file          = output_prefix + ".log"
    illumination_pdf  = output_prefix + "_illumination.pdf"
    surface_error_pdf = output_prefix + ".pdf"

    # Auto-version summary filename
    base_summary = output_prefix + "_summary"
    summary_pdf  = base_summary + ".pdf"
    summary_png  = base_summary + ".png"
    if os.path.exists(summary_pdf):
        counter = 1
        while os.path.exists(f"{base_summary}_{counter}.pdf"):
            counter += 1
        summary_pdf = f"{base_summary}_{counter}.pdf"
        summary_png = f"{base_summary}_{counter}.png"
        print(f"  Previous summary found; saving as {os.path.basename(summary_pdf)}")

    # --- extract fit / RMS data ---
    rms  = extract_rms_values(output_prefix)
    fit  = rms.get('fit') or {}
    coef = fit.get('coef', {})
    dofit = fit.get('dofit', [])

    # RMS string: standard (pixel-based) only, with quantity label
    if rms['standard'] is not None:
        rms_str = f"Surface RMS = {rms['standard']:.2f} µm"
    else:
        rms_str = "Surface RMS = (not available)"

    # Verify input PDFs exist
    if not os.path.exists(illumination_pdf):
        print(f"Error: Illumination PDF not found: {illumination_pdf}")
        return False
    if not os.path.exists(surface_error_pdf):
        print(f"Error: Surface error PDF not found: {surface_error_pdf}")
        return False

    # ------------------------------------------------------------------ #
    #  Figure layout                                                       #
    # ------------------------------------------------------------------ #
    fig = plt.figure(figsize=(16, 11))

    # Two rows: narrow header (22%) | wide figures (78%)
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[0.22, 0.78],
        hspace=0.08, wspace=0.06,
        left=0.03, right=0.97, top=0.97, bottom=0.03,
    )

    # ---- header left: text info ----------------------------------------
    ax_info = fig.add_subplot(gs[0, 0])
    ax_info.set_facecolor('#f7f9fc')
    for spine in ax_info.spines.values():
        spine.set_edgecolor('#c0c8d8')
    ax_info.tick_params(left=False, bottom=False,
                        labelleft=False, labelbottom=False)

    # filename (strip leading path for brevity)
    fname_short = os.path.basename(data_filename)
    info_lines = [
        (fname_short,         11, 'bold',   0.88),
        (rms_str,             10, 'normal', 0.68),
    ]
    if comment:
        import textwrap
        wrapped = textwrap.fill(f"Comment: {comment}", width=72)
        info_lines.append((wrapped, 9, 'normal', 0.42))

    for text, fs, weight, y in info_lines:
        ax_info.text(0.03, y, text,
                     fontsize=fs, weight=weight,
                     va='top', ha='left',
                     transform=ax_info.transAxes,
                     wrap=True,
                     clip_on=False)

    # ---- header right: aberration fit table ----------------------------
    ax_tbl = fig.add_subplot(gs[0, 1])
    ax_tbl.set_facecolor('#f7f9fc')
    for spine in ax_tbl.spines.values():
        spine.set_edgecolor('#c0c8d8')
    ax_tbl.tick_params(left=False, bottom=False,
                       labelleft=False, labelbottom=False)

    # Pre-compute conversion factors for astig/coma -> µm pk-pk at dish edge
    # peak wavefront error = |coef| × R0^n × λ/(4π) × 1e6  (µm)
    freq_ghz  = fit.get('freq_ghz', 94.5)
    lam_mm    = fit.get('wavelength_mm', 299.792458 / freq_ghz)
    lam_m     = lam_mm * 1e-3
    um_per_rad = lam_m / (4 * np.pi) * 1e6      # µm per radian of wavefront
    R0 = 6.0   # dish radius in metres (dprim/2 = 12/2)

    def to_um_edge(raw_coef, power):
        """Convert raw rad/m^n coefficient to µm peak at dish edge."""
        return raw_coef * (R0 ** power) * um_per_rad

    term_defs = [
        # (label, json_key, display_unit, edge_power)
        # edge_power=0 -> display raw value in given unit, no conversion
        # edge_power>0 -> convert to µm pk at edge
        ('DC offset',  'DC_offset_rad',      'rad',           0),
        ('Tilt X',     'tilt_x_rad_pix',     'rad/px',        0),
        ('Tilt Y',     'tilt_y_rad_pix',     'rad/px',        0),
        ('Defocus',    'defocus_mm',          'mm',            0),
        ('Astig 45°',  'astigmatism_45_deg',  'µm pk (edge)',  2),
        ('Astig 0°',   'astigmatism_0_deg',   'µm pk (edge)',  2),
        ('Coma X',     'coma_x',              'µm pk (edge)',  3),
        ('Coma Y',     'coma_y',              'µm pk (edge)',  3),
    ]

    if coef:
        rows = []
        for idx, (label, key, unit, edge_power) in enumerate(term_defs):
            val    = coef.get(key)
            fitted = bool(dofit[idx]) if idx < len(dofit) else True
            if val is not None:
                if edge_power > 0:
                    display_val = to_um_edge(val, edge_power)
                    val_str = f"{display_val:.2f}"
                else:
                    val_str = f"{val:.4e}"
                rows.append([label, val_str, unit, '✓' if fitted else '—'])

        if rows:
            rms_bef   = fit.get('rms_before_rad')
            rms_aft_u = fit.get('rms_after_um')

            title_parts = [f"Aberration fit  ({freq_ghz:.1f} GHz, λ={lam_mm:.3f} mm)"]
            if rms_bef and rms_aft_u:
                title_parts.append(
                    f"RMS before fit: {rms_bef * um_per_rad:.1f} µm  "
                    f"→  after fit: {rms_aft_u:.1f} µm  (half-path-length, pre-projection)"
                )
            ax_tbl.set_title('\n'.join(title_parts), fontsize=8, pad=3, loc='left')

            tbl = ax_tbl.table(
                cellText=rows,
                colLabels=['Term', 'Value', 'Unit', 'Fit'],
                loc='center',
                cellLoc='center',
            )
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(8)
            tbl.scale(1, 1.15)

            # Style header row
            for j in range(4):
                tbl[0, j].set_facecolor('#d0d8e8')
                tbl[0, j].set_text_props(weight='bold')
            # Grey out non-fitted rows
            for i, row in enumerate(rows):
                if row[3] == '—':
                    for j in range(4):
                        tbl[i+1, j].set_facecolor('#efefef')
                        tbl[i+1, j].set_text_props(color='#999999')
    else:
        ax_tbl.text(0.5, 0.5, "Aberration fit data\nnot available",
                    ha='center', va='center', fontsize=9,
                    transform=ax_tbl.transAxes, color='#888888')

    # ---- figures -------------------------------------------------------
    ax_left  = fig.add_subplot(gs[1, 0])
    ax_right = fig.add_subplot(gs[1, 1])
    ax_left.axis('off')
    ax_right.axis('off')

    try:
        from pdf2image import convert_from_path
        illum_imgs  = convert_from_path(illumination_pdf,  dpi=150)
        surface_imgs = convert_from_path(surface_error_pdf, dpi=150)

        if illum_imgs:
            ax_left.imshow(illum_imgs[0])
        if surface_imgs:
            ax_right.imshow(surface_imgs[0])

    except ImportError:
        print("Warning: pdf2image not available.")
        ax_left.text(0.5, 0.5,
                     f"Illumination Map\n{os.path.basename(illumination_pdf)}",
                     ha='center', va='center', fontsize=12,
                     transform=ax_left.transAxes)
        ax_right.text(0.5, 0.5,
                      f"Surface Error Map\n{os.path.basename(surface_error_pdf)}",
                      ha='center', va='center', fontsize=12,
                      transform=ax_right.transAxes)

    # ---- save ----------------------------------------------------------
    plt.savefig(summary_pdf, format='pdf', bbox_inches='tight', dpi=150)
    plt.savefig(summary_png, format='png', bbox_inches='tight', dpi=150)
    plt.close()

    print(f"Summary PDF created: {summary_pdf}")
    print(f"Summary PNG created: {summary_png}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Create a summary PDF for holography data reduction results.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('data_filename',
                       help='Original input data filename')
    parser.add_argument('output_prefix',
                       help='Output prefix for result files (e.g., results/holoADC-AzEl-20250928_170604)')
    parser.add_argument('--comment', '-c',
                       default=None,
                       help='Optional comment string to include in summary')

    args = parser.parse_args()

    success = create_summary_pdf(args.data_filename, args.output_prefix, args.comment)

    if not success:
        exit(1)


if __name__ == "__main__":
    main()
