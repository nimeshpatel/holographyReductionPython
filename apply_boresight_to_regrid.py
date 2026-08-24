#!/usr/bin/env python3
"""
apply_boresight_to_regrid.py - Apply boresight phase correction to regridded data

This script applies phase drift correction AFTER regridding, using the mean
timestamp of each grid cell. This may be more effective than pre-regrid
correction because:
1. Each cell gets a single, clean correction value
2. The regridded phase is a cleaner average of many samples
3. No within-cell correction fighting

Usage:
    python apply_boresight_to_regrid.py rgin_with_time.dat boresight_data.txt -o rgin.dat

Input files:
    rgin_with_time.dat: Regridded data with timestamps (j k amp phase mean_time)
    boresight_data.txt: Boresight measurements from boresight_cal.py --save-boresight-data

Output:
    rgin.dat: Corrected regridded data (j k amp corrected_phase) - no timestamps
"""

import argparse
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import median_filter
import sys


def robust_unwrap_phase(phases):
    """
    Robustly unwrap a sequence of phase values.
    """
    if len(phases) == 0:
        return phases.copy()
    
    unwrapped = np.zeros_like(phases)
    unwrapped[0] = phases[0]
    
    for i in range(1, len(phases)):
        diff = phases[i] - phases[i-1]
        # Unwrap: if jump > 180, subtract 360; if < -180, add 360
        while diff > 180:
            diff -= 360
        while diff < -180:
            diff += 360
        unwrapped[i] = unwrapped[i-1] + diff
    
    return unwrapped


def load_boresight_data(filepath, verbose=False):
    """
    Load boresight measurements from file saved by boresight_cal.py.
    
    Expected format: time amplitude phase (3 columns)
    """
    if verbose:
        print(f"Loading boresight data from {filepath}...")
    
    data = np.loadtxt(filepath)
    if len(data.shape) == 1:
        data = data.reshape(1, -1)
    
    times = data[:, 0]
    amps = data[:, 1]
    phases = data[:, 2]
    
    if verbose:
        print(f"  Loaded {len(times)} boresight measurements")
        print(f"  Time range: {times.min():.1f} to {times.max():.1f}")
        print(f"  Phase range: {phases.min():.1f} to {phases.max():.1f} deg")
    
    return times, amps, phases


def fit_phase_drift(bore_times, bore_phases, smoothing=None, median_window=5, verbose=False):
    """
    Fit phase drift from boresight measurements.
    
    Returns a function that gives phase drift at any time.
    """
    if verbose:
        print("Fitting phase drift...")
    
    # Unwrap phases
    unwrapped = robust_unwrap_phase(bore_phases)
    
    if verbose:
        print(f"  Unwrapped phase range: {unwrapped.min():.1f} to {unwrapped.max():.1f} deg")
    
    # Apply median filter to reject outliers
    if median_window is not None and median_window > 1:
        filtered = median_filter(unwrapped, size=median_window, mode='nearest')
        if verbose:
            diff = np.abs(unwrapped - filtered)
            n_changed = np.sum(diff > 1.0)
            print(f"  Median filter (window={median_window}): {n_changed} points adjusted")
        unwrapped = filtered
    
    # Reference to first measurement
    ref_phase = unwrapped[0]
    drift_at_bore = unwrapped - ref_phase
    
    if verbose:
        print(f"  Drift range: {drift_at_bore.min():.2f} to {drift_at_bore.max():.2f} deg")
    
    # Fit spline
    if smoothing is None:
        smoothing = len(bore_times) * 4.0
    
    if len(bore_times) < 3:
        # Linear interpolation
        def drift_func(t):
            return np.interp(t, bore_times, drift_at_bore)
    else:
        try:
            spline = UnivariateSpline(bore_times, drift_at_bore, s=smoothing, k=3)
            
            if verbose:
                fit_at_bore = spline(bore_times)
                residuals = drift_at_bore - fit_at_bore
                print(f"  Spline fit RMS residual: {np.std(residuals):.2f} deg")
            
            drift_func = spline
        except Exception as e:
            print(f"  Warning: Spline fit failed ({str(e)}), using linear interpolation")
            def drift_func(t):
                return np.interp(t, bore_times, drift_at_bore)
    
    return drift_func, ref_phase


def apply_correction(regrid_file, boresight_file, output_file, smoothing=None, 
                     median_window=5, verbose=False):
    """
    Apply boresight correction to regridded data.
    """
    # Load boresight data
    bore_times, bore_amps, bore_phases = load_boresight_data(boresight_file, verbose)
    
    # Fit phase drift
    drift_func, ref_phase = fit_phase_drift(bore_times, bore_phases, 
                                            smoothing=smoothing,
                                            median_window=median_window,
                                            verbose=verbose)
    
    # Load regridded data
    if verbose:
        print(f"\nLoading regridded data from {regrid_file}...")
    
    data = np.loadtxt(regrid_file)
    n_cells = len(data)
    
    if data.shape[1] < 5:
        print(f"Error: Expected 5 columns (j k amp phase time), got {data.shape[1]}")
        print("  Make sure to use regrid_holo_with_time.py to generate the input file")
        sys.exit(1)
    
    j_idx = data[:, 0].astype(int)
    k_idx = data[:, 1].astype(int)
    amp = data[:, 2]
    phase = data[:, 3]  # In radians
    mean_time = data[:, 4]
    
    if verbose:
        valid = amp > 0.01 * amp.max()
        print(f"  Loaded {n_cells} grid cells")
        print(f"  High-amplitude cells: {np.sum(valid)}")
        print(f"  Time range of cells: {mean_time[valid].min():.1f} to {mean_time[valid].max():.1f}")
    
    # Calculate correction for each cell
    corrections = drift_func(mean_time)  # In degrees
    
    if verbose:
        print(f"\nApplying corrections...")
        print(f"  Correction range: {corrections.min():.2f} to {corrections.max():.2f} deg")
    
    # Apply correction (phase is in radians, correction is in degrees)
    phase_corrected = phase - np.radians(corrections)
    
    # Wrap to [-pi, pi]
    phase_corrected = np.mod(phase_corrected + np.pi, 2*np.pi) - np.pi
    
    # Write output (4 columns: j k amp phase, no timestamps)
    if verbose:
        print(f"\nWriting corrected data to {output_file}...")
    
    with open(output_file, 'w') as f:
        for i in range(n_cells):
            f.write(f"{j_idx[i]} {k_idx[i]} {amp[i]:.16e} {phase_corrected[i]:.16e}\n")
    
    if verbose:
        print(f"  Wrote {n_cells} cells")
        
        # Show correction statistics for high-amplitude cells
        valid = amp > 0.01 * amp.max()
        corr_valid = corrections[valid]
        print(f"\nCorrection statistics (high-amplitude cells):")
        print(f"  Mean: {np.mean(corr_valid):.2f} deg")
        print(f"  Std:  {np.std(corr_valid):.2f} deg")
        print(f"  Range: {corr_valid.min():.2f} to {corr_valid.max():.2f} deg")
    
    return n_cells


def main():
    parser = argparse.ArgumentParser(
        description="Apply boresight phase correction to regridded holography data"
    )
    parser.add_argument("regrid_file", type=str,
                        help="Input regridded data with timestamps (j k amp phase time)")
    parser.add_argument("boresight_file", type=str,
                        help="Boresight data file (time amp phase)")
    parser.add_argument("-o", "--output", type=str, default="rgin.dat",
                        help="Output file (default: rgin.dat)")
    parser.add_argument("--smoothing", type=float, default=None,
                        help="Spline smoothing factor (default: auto)")
    parser.add_argument("--median-window", type=int, default=5,
                        help="Median filter window size (default: 5)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")
    
    args = parser.parse_args()
    
    n_cells = apply_correction(
        args.regrid_file,
        args.boresight_file,
        args.output,
        smoothing=args.smoothing,
        median_window=args.median_window,
        verbose=args.verbose
    )
    
    print(f"\nPost-regrid boresight correction complete.")
    print(f"  Input:  {args.regrid_file}")
    print(f"  Output: {args.output}")


if __name__ == "__main__":
    main()
