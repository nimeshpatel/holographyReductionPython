#!/usr/bin/env python3
"""
Nimesh Patel
October 2025
Python version of previous fortran/C code from Holis package.

detect_raster_start_interactive_v3.py — Interactive raster start detection

Features:
  ✓ Auto-detects raster scan start using sustained positive az slope + flat elevation
  ✓ Shows dual plots: full trajectory (left) + zoomed view (right)
  ✓ Interactive cursor selection on zoomed plot for manual adjustment
  ✓ Simplified workflow:
    - Accept detected start (y): write trimmed file immediately
    - Reject (n): enter line number or use cursor selection, replot zoomed view only
  ✓ Output file: numeric-only, columns = az, el, amp, phase
  ✓ Inspection plots (3-panel) commented out - moved to regridding program
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Tuple, List


def read_ascii(path: str, has_header: bool = False, comment: str = "#") -> Tuple[pd.DataFrame, List[str]]:
    header_lines, rows = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip().startswith(comment) or line.strip() == "":
                header_lines.append(line.rstrip("\n"))
            else:
                rows.append(line)
    if not rows:
        raise ValueError("No data rows found.")
    tmp = "".join(rows)
    if has_header:
        df = pd.read_csv(pd.io.common.StringIO(tmp), delim_whitespace=True, comment=comment)
    else:
        df = pd.read_csv(pd.io.common.StringIO(tmp), delim_whitespace=True, header=None,
                         names=["line_no", "az", "el", "amp", "phase"])
    for c in ["line_no", "az", "el", "amp", "phase"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna().reset_index(drop=True), header_lines


def rolling_median(a: pd.Series, window: int) -> pd.Series:
    return a.rolling(window=window, center=True, min_periods=max(1, window//2)).median() if window > 1 else a.copy()


def _auto_slope_threshold(series_sm: pd.Series, slope_window: int) -> float:
    span = min(len(series_sm)//2, max(20, 5*slope_window))
    diffs = series_sm.diff().iloc[:span].dropna()
    mad = np.median(np.abs(diffs - np.median(diffs))) if len(diffs) else 0.0
    return max(0.25 * mad, 1e-6)


def find_first_sustained_positive_slope_with_flat_el(az, el, smooth_window=7, slope_window=7,
                                                     min_consecutive=3, slope_threshold=None,
                                                     el_slope_threshold=None, el_min_flat_consecutive=None) -> int:
    if el_min_flat_consecutive is None:
        el_min_flat_consecutive = min_consecutive
    az_sm, el_sm = rolling_median(az, smooth_window), rolling_median(el, smooth_window)
    slope_az = (az_sm - az_sm.shift(slope_window)) / float(slope_window)
    slope_el = (el_sm - el_sm.shift(slope_window)) / float(slope_window)
    thr_az = float(slope_threshold) if slope_threshold else _auto_slope_threshold(az_sm, slope_window)
    thr_el = float(el_slope_threshold) if el_slope_threshold else _auto_slope_threshold(el_sm, slope_window)
    pos_az, flat_el = slope_az > thr_az, slope_el.abs() <= thr_el
    both = pos_az & flat_el
    count = 0
    for i in range(len(both)):
        count = count + 1 if both.iloc[i] else 0
        if count >= min(min_consecutive, el_min_flat_consecutive):
            start_idx = i - count + 1
            pre_lo, pre_hi = max(0, start_idx - 2*slope_window), max(0, start_idx - 1)
            if (not pos_az.iloc[pre_lo:pre_hi].any()) and slope_az.iloc[pre_lo:pre_hi].mean() <= thr_az / 4.0:
                return int(start_idx)
    raise RuntimeError("No sustained az-positive + el-flat segment found.")


def quick_inspection_plots(df: pd.DataFrame):
    fig, axs = plt.subplots(3, 1, figsize=(8, 10))
    axs[0].scatter(df["az"], df["el"], s=0.2)
    axs[0].set_xlabel("Azimuth")
    axs[0].set_ylabel("Elevation")
    axs[0].set_title("Azimuth vs Elevation (raster grid)")
    axs[0].set_aspect("equal", adjustable="box")
    axs[1].scatter(df["az"], df["amp"], s=0.2)
    axs[1].set_xlabel("Azimuth")
    axs[1].set_ylabel("Amplitude")
    axs[1].set_title("Amplitude vs Azimuth")
    axs[2].scatter(df["el"], df["amp"], s=0.2)
    axs[2].set_xlabel("Elevation")
    axs[2].set_ylabel("Amplitude")
    axs[2].set_title("Amplitude vs Elevation")
    plt.tight_layout()
    plt.show()


def plot_trajectory_dual(df: pd.DataFrame, start_idx: int, color="red", label_prefix="Detected start", zoom_range=500):
    """Plot full trajectory and zoomed view side by side (non-blocking)."""
    fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(16, 5))

    # Full trajectory plot (left)
    ax1.plot(df.index, df["az"], label="Azimuth", color="tab:blue", linewidth=0.8)
    ax1.set_xlabel("Index number")
    ax1.set_ylabel("Azimuth", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    ax1.set_title("Full Trajectory")

    ax2 = ax1.twinx()
    ax2.plot(df.index, df["el"], label="Elevation", color="tab:orange", alpha=0.7, linewidth=0.8)
    ax2.set_ylabel("Elevation", color="tab:orange")
    ax2.tick_params(axis="y", labelcolor="tab:orange")

    ax1.axvline(x=start_idx, color=color, linestyle="--", linewidth=2,
                label=f"{label_prefix} (line {int(start_idx)})")

    # Zoomed trajectory plot (right)
    zoom_start = max(0, start_idx - zoom_range//2)
    zoom_end = min(len(df), start_idx + zoom_range//2)
    zoom_df = df.iloc[zoom_start:zoom_end]

    ax3.plot(zoom_df.index, zoom_df["az"], label="Azimuth", color="tab:blue", linewidth=1.2)
    ax3.set_xlabel("Index number")
    ax3.set_ylabel("Azimuth", color="tab:blue")
    ax3.tick_params(axis="y", labelcolor="tab:blue")
    ax3.set_title(f"Zoomed View (±{zoom_range//2} points) - Click to select")

    ax4 = ax3.twinx()
    ax4.plot(zoom_df.index, zoom_df["el"], label="Elevation", color="tab:orange", alpha=0.7, linewidth=1.2)
    ax4.set_ylabel("Elevation", color="tab:orange")
    ax4.tick_params(axis="y", labelcolor="tab:orange")

    ax3.axvline(x=start_idx, color=color, linestyle="--", linewidth=2,
                label=f"{label_prefix} (line {int(start_idx)})")

    fig.suptitle("Raster Trajectory with Start Point", fontsize=14, fontweight='bold')
    fig.tight_layout()
    ax1.legend(loc="upper left")
    ax3.legend(loc="upper left")

    # Enable cursor interaction
    cursor_coords = {'x': start_idx}

    def onclick(event):
        if event.inaxes in [ax3]:
            cursor_coords['x'] = int(round(event.xdata)) if event.xdata is not None else start_idx
            print(f"  --> Cursor clicked at line number: {cursor_coords['x']}")

    fig.canvas.mpl_connect('button_press_event', onclick)

    # Non-blocking show - plot stays open
    plt.show(block=False)
    plt.pause(0.1)

    return fig, cursor_coords


def plot_trajectory_zoomed_only(df: pd.DataFrame, start_idx: int, color="green", label_prefix="Selected start", zoom_range=500):
    """Plot only zoomed view for manual adjustment (non-blocking)."""
    fig, ax1 = plt.subplots(figsize=(12, 5))

    zoom_start = max(0, start_idx - zoom_range//2)
    zoom_end = min(len(df), start_idx + zoom_range//2)
    zoom_df = df.iloc[zoom_start:zoom_end]

    ax1.plot(zoom_df.index, zoom_df["az"], label="Azimuth", color="tab:blue", linewidth=1.2)
    ax1.set_xlabel("Index number")
    ax1.set_ylabel("Azimuth", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    ax1.set_title(f"Zoomed View (±{zoom_range//2} points) - Click to adjust if needed")

    ax2 = ax1.twinx()
    ax2.plot(zoom_df.index, zoom_df["el"], label="Elevation", color="tab:orange", alpha=0.7, linewidth=1.2)
    ax2.set_ylabel("Elevation", color="tab:orange")
    ax2.tick_params(axis="y", labelcolor="tab:orange")

    ax1.axvline(x=start_idx, color=color, linestyle="--", linewidth=2,
                label=f"{label_prefix} (line {int(start_idx)})")

    fig.tight_layout()
    ax1.legend(loc="upper left")

    cursor_coords = {'x': start_idx}

    def onclick(event):
        if event.inaxes == ax1:
            cursor_coords['x'] = int(round(event.xdata)) if event.xdata is not None else start_idx
            print(f"  --> Cursor clicked at line number: {cursor_coords['x']}")

    fig.canvas.mpl_connect('button_press_event', onclick)

    # Non-blocking show - plot stays open
    plt.show(block=False)
    plt.pause(0.1)

    return fig, cursor_coords


def main():
    ap = argparse.ArgumentParser(description="Raster start detection and trimmed-data visualization.")
    ap.add_argument("input")
    ap.add_argument("-o", "--output", default=None)
    ap.add_argument("--smooth-window", type=int, default=7)
    ap.add_argument("--slope-window", type=int, default=7)
    ap.add_argument("--min-consecutive", type=int, default=3)
    ap.add_argument("--slope-threshold", type=float, default=None)
    ap.add_argument("--el-slope-threshold", type=float, default=None)
    ap.add_argument("--el-min-flat-consecutive", type=int, default=None)
    ap.add_argument("--assume-header", action="store_true")
    args = ap.parse_args()

    df, _ = read_ascii(args.input, has_header=args.assume_header)

    start_idx = find_first_sustained_positive_slope_with_flat_el(
        df["az"], df["el"], args.smooth_window, args.slope_window,
        args.min_consecutive, args.slope_threshold,
        args.el_slope_threshold, args.el_min_flat_consecutive)

    print(f"\nDetected start at line number = {start_idx}")
    print("Showing full trajectory (left) and zoomed view (right)...")
    print("You can click on the zoomed plot (right) to select a different start point.\n")

    # Show dual plot (full + zoomed) - non-blocking
    fig1, cursor_coords = plot_trajectory_dual(df, start_idx, "red", "Detected start")

    # Prompt immediately while plot is visible
    choice = input("Accept [y/n]: ").strip().lower()

    if choice == "y":
        # User accepted - close plot and proceed
        plt.close(fig1)
        print(f"Using detected start line number = {start_idx}")
    else:
        # User rejected - prompt for manual input immediately
        print("\nYou can:")
        print("  1. Enter a line number directly")
        print("  2. Press Enter to use the cursor-clicked value (if you clicked on the zoomed plot)")
        manual = input("\nEnter line number (or press Enter for cursor value): ").strip()

        # Close the first plot
        plt.close(fig1)

        if manual != "":
            # User typed a line number
            try:
                start_idx = int(manual)
                print(f"\nUsing manually entered line number = {start_idx}")
            except ValueError:
                print("\nInvalid input. Using cursor-selected value.")
                start_idx = cursor_coords['x']
        else:
            # Use cursor-clicked value
            start_idx = cursor_coords['x']
            print(f"\nUsing cursor-selected line number = {start_idx}")

        # Show zoomed plot only for verification
        print("Showing zoomed view for verification. Click to adjust if needed.\n")
        fig2, final_coords = plot_trajectory_zoomed_only(df, start_idx, "green", "Manual start")

        # Ask for final confirmation
        final_choice = input("Accept this start point [y/n]: ").strip().lower()
        plt.close(fig2)

        if final_choice == "y":
            start_idx = final_coords['x']
            print(f"Final start line number = {start_idx}")
        else:
            print("Aborted. No file written.")
            sys.exit(0)

    trimmed = df.iloc[start_idx:].reset_index(drop=True)

    # Commented out: inspection plots moved to regridding program
    # print("Displaying inspection plots (az/el, amp vs az, amp vs el) from trimmed data...")
    # quick_inspection_plots(trimmed)

    out_path = args.output or (args.input[:-4]+"_trimmed.txt"
                               if args.input.lower().endswith(".txt")
                               else args.input+"_trimmed.txt")

    # Fast file writing using numpy (much faster than pandas to_csv)
    print(f"\nWriting trimmed file with {len(trimmed)} data points...")
    trimmed_array = trimmed[["az", "el", "amp", "phase"]].values
    np.savetxt(out_path, trimmed_array, fmt="%.9g", delimiter=" ")
    print(f"Wrote trimmed numeric-only file: {out_path}")
    print(f"Trimmed data starts at line {start_idx}, contains {len(trimmed)} points.\n")

    # Close any remaining plots
    plt.close('all')


if __name__ == "__main__":
    main()
