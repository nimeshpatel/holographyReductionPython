#!/usr/bin/env python3
"""
detect_raster_start_interactive_v3.py — final streamlined version

Changes in v3:
  ✓ Reports start using index number, not line_no.
  ✓ Prompts: "Detected start at line number = NNNN" / "Accept [y/n]:"
  ✓ Output file contains only columns: el, az, amp, phase (az/el swapped order).
  ✓ Numeric-only file — no headers or comments.
  ✓ Uses scatter-only plots with small dots and square az/el aspect ratio.
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


def plot_trajectory(df: pd.DataFrame, start_idx: int, color="red", label_prefix="Detected start"):
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.plot(df.index, df["az"], label="Azimuth", color="tab:blue")
    ax1.set_xlabel("Index number")
    ax1.set_ylabel("Azimuth", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")

    ax2 = ax1.twinx()
    ax2.plot(df.index, df["el"], label="Elevation", color="tab:orange", alpha=0.7)
    ax2.set_ylabel("Elevation", color="tab:orange")
    ax2.tick_params(axis="y", labelcolor="tab:orange")

    ax1.axvline(x=start_idx, color=color, linestyle="--",
                label=f"{label_prefix} (line {int(start_idx)})")

    fig.suptitle("Raster Trajectory with Start Point")
    fig.tight_layout()
    fig.legend(loc="upper right")
    plt.show()


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

    print(f"Detected start at line number = {start_idx}")
    plot_trajectory(df, start_idx, "red", "Detected start")

    choice = input("Accept [y/n]: ").strip().lower()
    if choice != "y":
        print("Replotting for manual inspection...")
        plot_trajectory(df, start_idx, "red", "Detected start")
        manual = input("Enter line number to trim from (or press Enter to cancel): ").strip()
        if manual == "":
            print("Aborted. No file written.")
            sys.exit(0)
        try:
            manual_val = int(manual)
        except ValueError:
            print("Invalid input.")
            sys.exit(1)
        start_idx = manual_val
        print(f"Using manually selected line number = {start_idx}")
        plot_trajectory(df, start_idx, "green", "Manual start (confirmed)")

    trimmed = df.iloc[start_idx:].reset_index(drop=True)

    print("Displaying inspection plots (az/el, amp vs az, amp vs el) from trimmed data...")
    quick_inspection_plots(trimmed)

    out_path = args.output or (args.input[:-4]+"_trimmed.txt"
                               if args.input.lower().endswith(".txt")
                               else args.input+"_trimmed.txt")

    trimmed_out = trimmed[["az", "el", "amp", "phase"]]
    trimmed_out.to_csv(out_path, sep=" ", index=False, header=False, float_format="%.9g")
    print(f"Wrote trimmed numeric-only file: {out_path}")


if __name__ == "__main__":
    main()
