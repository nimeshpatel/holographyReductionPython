#!/usr/bin/env python3
"""
create_summary_pdf.py

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


def extract_rms_from_log(log_file):
    """
    Extract the final RMS surface error value in microns from the holis log file.

    The log file contains multiple RMS values. We want the final converged value
    which appears near the end with the format:
    "RMS after phasefit = X.XXXXXX radian (YY.YYYY microns)"
    """
    if not os.path.exists(log_file):
        print(f"Warning: Log file {log_file} not found")
        return None

    rms_microns = None
    with open(log_file, 'r') as f:
        for line in f:
            # Look for lines with RMS in microns
            match = re.search(r'RMS after phasefit.*\(([0-9.e+\-]+)\s+microns\)', line)
            if match:
                rms_microns = float(match.group(1))

    return rms_microns


def create_summary_pdf(data_filename, output_prefix, comment=None):
    """
    Create a composite PDF summary page with text and embedded figures.

    Parameters:
    -----------
    data_filename : str
        Original input data filename
    output_prefix : str
        Prefix for output files (e.g., results/holoADC-AzEl-20250928_170604)
    comment : str, optional
        User comment to include in the summary
    """

    # Construct file paths
    log_file = output_prefix + ".log"
    illumination_pdf = output_prefix + "_illumination.pdf"
    surface_error_pdf = output_prefix + ".pdf"
    summary_pdf = output_prefix + "_summary.pdf"

    # Extract RMS value
    rms_microns = extract_rms_from_log(log_file)

    # Verify input PDFs exist
    if not os.path.exists(illumination_pdf):
        print(f"Error: Illumination PDF not found: {illumination_pdf}")
        return False
    if not os.path.exists(surface_error_pdf):
        print(f"Error: Surface error PDF not found: {surface_error_pdf}")
        return False

    # Create figure with custom layout
    fig = plt.figure(figsize=(16, 11))

    # Create grid for layout: text at top, two images below
    # Use gridspec for better control
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 3], hspace=0.15, wspace=0.1,
                         left=0.05, right=0.95, top=0.95, bottom=0.05)

    # Top section for text (spans both columns)
    ax_text = fig.add_subplot(gs[0, :])
    ax_text.axis('off')

    # Build text content
    text_lines = []
    text_lines.append(f"Holography map data: {data_filename}")

    if rms_microns is not None:
        text_lines.append(f"Surface rms error = {rms_microns:.2f} microns")
    else:
        text_lines.append("Surface rms error = (not available)")

    if comment:
        text_lines.append(f"Comment: {comment}")

    # Position text in the text area
    y_start = 0.85
    y_step = 0.35

    for i, line in enumerate(text_lines):
        y_pos = y_start - i * y_step
        # Determine font size based on line content
        if i == 0:  # filename
            fontsize = 24
            weight = 'normal'
        elif i == 1:  # rms error
            fontsize = 18
            weight = 'normal'
        else:  # comment
            fontsize = 18
            weight = 'normal'

        ax_text.text(0.5, y_pos, line,
                    fontsize=fontsize,
                    weight=weight,
                    ha='center', va='center',
                    transform=ax_text.transAxes)

    # Bottom left: Illumination figure
    ax_left = fig.add_subplot(gs[1, 0])
    ax_left.axis('off')

    # Bottom right: Surface error figure
    ax_right = fig.add_subplot(gs[1, 1])
    ax_right.axis('off')

    # Load and display PDF images
    # We'll use PIL to convert PDF pages to images
    try:
        from pdf2image import convert_from_path

        # Convert PDFs to images
        illumination_images = convert_from_path(illumination_pdf, dpi=150)
        surface_error_images = convert_from_path(surface_error_pdf, dpi=150)

        if illumination_images:
            ax_left.imshow(illumination_images[0])
            ax_left.set_title("Illumination Map", fontsize=16, weight='normal', pad=10)

        if surface_error_images:
            ax_right.imshow(surface_error_images[0])
            ax_right.set_title("Surface Error Map", fontsize=16, weight='normal', pad=10)

    except ImportError:
        # If pdf2image is not available, add text placeholders
        print("Warning: pdf2image not available. Using matplotlib PDF embedding.")
        print("Install pdf2image with: pip install pdf2image")
        print("You may also need poppler-utils: brew install poppler (macOS) or apt-get install poppler-utils (Linux)")

        # Fallback: show file names
        ax_left.text(0.5, 0.5, f"Illumination Map\n\n{os.path.basename(illumination_pdf)}",
                    ha='center', va='center', fontsize=16, transform=ax_left.transAxes)
        ax_right.text(0.5, 0.5, f"Surface Error Map\n\n{os.path.basename(surface_error_pdf)}",
                     ha='center', va='center', fontsize=16, transform=ax_right.transAxes)

    # Save to PDF
    plt.savefig(summary_pdf, format='pdf', bbox_inches='tight', dpi=150)
    plt.close()

    print(f"Summary PDF created: {summary_pdf}")
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
