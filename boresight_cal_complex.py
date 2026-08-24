#!/usr/bin/env python3
"""
boresight_cal_complex.py - Boresight drift calibration for GLT holography data

COMPLEX CORRECTION VERSION - applies both amplitude and phase correction.

The correction is:
    S_corrected = S_measured × (A_ref / A_bore(t)) × exp(-i × phase_drift(t))

This properly handles the case where both amplitude and phase drift over time.

Usage:
    python boresight_cal_complex.py input_raw.txt -o output_cal.txt [options]
    
Modes:
    Default:              Apply complex (amp + phase) calibration
    --no-correction:      Remove boresight points WITHOUT applying correction
"""

import argparse
import numpy as np
from scipy.interpolate import UnivariateSpline
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calibrate holography data using boresight drift measurements"
    )
    parser.add_argument("input_file", type=str,
                        help="Input raw data file (with timestamps)")
    parser.add_argument("-o", "--output", type=str, default="calibrated.txt",
                        help="Output calibrated data file (default: calibrated.txt)")
    parser.add_argument("--center-az", type=float, default=230.179,
                        help="Beacon center azimuth in degrees (default: 230.179)")
    parser.add_argument("--center-el", type=float, default=2.8724,
                        help="Beacon center elevation in degrees (default: 2.8724)")
    parser.add_argument("--tolerance", type=float, default=0.01,
                        help="Position tolerance for boresight detection in degrees (default: 0.01)")
    parser.add_argument("--min-samples", type=int, default=100,
                        help="Minimum samples to identify a segment (default: 100)")
    parser.add_argument("--min-bore-samples", type=int, default=5000,
                        help="Minimum samples for a REAL boresight visit (default: 5000)")
    parser.add_argument("--min-amp", type=float, default=0.5,
                        help="Minimum amplitude for valid boresight segment (default: 0.5)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Print detailed progress information")
    parser.add_argument("--remove-boresight", action="store_true",
                        help="Remove boresight data points from output")
    parser.add_argument("--slew-margin", type=float, default=0.5,
                        help="Time margin (seconds) around boresight segments to remove (default: 0.5)")
    # Note: --correct-amplitude is removed; this version always does complex correction
    parser.add_argument("--smoothing", type=float, default=None,
                        help="Spline smoothing factor (default: auto)")
    parser.add_argument("--linear", action="store_true",
                        help="Use piecewise linear interpolation instead of spline")
    parser.add_argument("--median-window", type=int, default=None,
                        help="Median filter window size to reject outlier jumps (e.g., 5)")
    parser.add_argument("--median-threshold", type=float, default=None,
                        help="Only filter points deviating by more than this many degrees")
    parser.add_argument("--plot", action="store_true",
                        help="Show diagnostic plots")
    parser.add_argument("--no-correction", action="store_true",
                        help="Remove boresight points but do NOT apply any correction")
    parser.add_argument("--save-prefix", type=str, default=None,
                        help="Prefix for saving diagnostic files (default: auto from input filename)")
    parser.add_argument("--save-boresight-data", action="store_true",
                        help="Save boresight measurements to a text file")
    return parser.parse_args()


def parse_timestamp(ts_str):
    """Parse timestamp string to seconds."""
    ts_str = str(ts_str).strip()
    if ":" in ts_str:
        parts = ts_str.split(":")
        hours = float(parts[0])
        minutes = float(parts[1])
        seconds = float(parts[2]) if len(parts) > 2 else 0.0
        return hours * 3600 + minutes * 60 + seconds
    else:
        return float(ts_str)


def read_raw_data(filepath, verbose=False):
    """Read raw holography data file."""
    if verbose:
        print("Reading " + str(filepath) + "...")
    
    header = None
    data_lines = []
    
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                header = line
                continue
            data_lines.append(line)
    
    if verbose:
        print("  Found " + str(len(data_lines)) + " data lines")
    
    n = len(data_lines)
    timestamps = np.zeros(n)
    az = np.zeros(n)
    el = np.zeros(n)
    amplitude = np.zeros(n)
    phase = np.zeros(n)
    
    for i, line in enumerate(data_lines):
        parts = line.split()
        if len(parts) >= 5:
            timestamps[i] = parse_timestamp(parts[0])
            az[i] = float(parts[1])
            el[i] = float(parts[2])
            amplitude[i] = float(parts[3])
            phase[i] = float(parts[4])
    
    return timestamps, az, el, amplitude, phase, header


def identify_boresight_segments(timestamps, az, el, center_az, center_el, 
                                 tolerance, min_samples, verbose=False):
    """Identify boresight measurement segments."""
    if verbose:
        print("Identifying boresight segments...")
    
    dist = np.sqrt((az - center_az)**2 + (el - center_el)**2)
    is_boresight = dist < tolerance
    
    segments = []
    in_segment = False
    start_idx = 0
    
    for i in range(len(is_boresight)):
        if is_boresight[i] and not in_segment:
            in_segment = True
            start_idx = i
        elif not is_boresight[i] and in_segment:
            in_segment = False
            end_idx = i
            if end_idx - start_idx >= min_samples:
                segments.append((start_idx, end_idx))
    
    if in_segment and len(timestamps) - start_idx >= min_samples:
        segments.append((start_idx, len(timestamps)))
    
    if verbose:
        print("  Found " + str(len(segments)) + " segments")
    
    return segments


def compute_boresight_values(timestamps, amplitude, phase, segments, min_bore_samples, min_amp=0.5, verbose=False):
    """Compute mean amplitude and phase for valid boresight segments."""
    if verbose:
        print("Analyzing boresight segments...")
    
    bore_times = []
    bore_amps = []
    bore_phases = []
    valid_segments = []
    
    for i, (start, end) in enumerate(segments):
        n_samples = end - start
        mean_time = np.mean(timestamps[start:end])
        mean_amp = np.mean(amplitude[start:end])
        
        # Robust phase average using complex mean
        phase_rad = np.deg2rad(phase[start:end])
        complex_avg = np.mean(np.exp(1j * phase_rad))
        mean_phase = np.rad2deg(np.angle(complex_avg))
        
        is_valid = n_samples >= min_bore_samples and mean_amp > min_amp
        
        if verbose:
            status = "VALID" if is_valid else "SKIP"
            print("  Seg {:2d}: t={:.1f}, amp={:.3f}, phase={:7.2f} deg, n={:5d} [{}]".format(
                i+1, mean_time, mean_amp, mean_phase, n_samples, status))
        
        if is_valid:
            bore_times.append(mean_time)
            bore_amps.append(mean_amp)
            bore_phases.append(mean_phase)
            valid_segments.append((start, end))
    
    if verbose:
        print("  Valid segments: " + str(len(valid_segments)))
    
    return np.array(bore_times), np.array(bore_amps), np.array(bore_phases), valid_segments


def robust_unwrap_phase(phases):
    """
    Robustly unwrap a sequence of phase values.
    
    Handles cases where phases span the ±180° boundary.
    Uses sequential unwrapping with 180° threshold.
    """
    if len(phases) == 0:
        return phases.copy()
    
    unwrapped = np.zeros_like(phases)
    unwrapped[0] = phases[0]
    
    for i in range(1, len(phases)):
        diff = phases[i] - phases[i-1]
        
        # Wrap difference to [-180, 180]
        while diff > 180:
            diff -= 360
        while diff < -180:
            diff += 360
        
        unwrapped[i] = unwrapped[i-1] + diff
    
    return unwrapped


def median_filter_phases(phases, window_size=5, threshold=None, verbose=False):
    """
    Apply median filtering to reject outlier phase jumps.
    
    Args:
        phases: Array of phase measurements
        window_size: Size of median filter window (must be odd)
        threshold: If provided, only replace points that deviate by more than
                   this many degrees from the local median. If None, always use median.
        verbose: Print progress
    
    Returns:
        filtered_phases: Median-filtered phases
        n_replaced: Number of points that were modified
    """
    from scipy.ndimage import median_filter
    
    if window_size % 2 == 0:
        window_size += 1  # Ensure odd window size
    
    # Apply median filter
    median_phases = median_filter(phases, size=window_size, mode='nearest')
    
    if threshold is not None:
        # Only replace points that deviate significantly from local median
        deviation = np.abs(phases - median_phases)
        outliers = deviation > threshold
        filtered_phases = phases.copy()
        filtered_phases[outliers] = median_phases[outliers]
        n_replaced = np.sum(outliers)
    else:
        # Replace all with median
        filtered_phases = median_phases
        n_replaced = len(phases)
    
    if verbose:
        print("  Median filter: window={}, threshold={}°".format(
            window_size, threshold if threshold else "none"))
        print("  Points filtered: {} of {}".format(n_replaced, len(phases)))
    
    return filtered_phases, n_replaced


def fit_phase_drift(bore_times, bore_phases, timestamps, smoothing=None, use_linear=False, 
                    median_window=None, median_threshold=None, verbose=False):
    """
    Fit phase drift using either piecewise linear interpolation or spline.
    
    Args:
        bore_times: Array of boresight measurement times
        bore_phases: Array of boresight phase measurements (degrees)
        timestamps: Array of all data timestamps to interpolate to
        smoothing: Spline smoothing factor (ignored if use_linear=True)
        use_linear: If True, use piecewise linear interpolation
        median_window: If provided, apply median filter with this window size
        median_threshold: Only filter points deviating by more than this (degrees)
        verbose: Print progress
    
    Returns:
        phase_drift: Interpolated drift at each timestamp
        ref_phase: Reference phase (first measurement)
        unwrapped_phases: Unwrapped boresight phases (after filtering if applied)
    """
    method = "piecewise linear" if use_linear else "spline"
    if verbose:
        print("Fitting phase drift ({})...".format(method))
        print("  Raw phase range: {:.2f} to {:.2f} deg".format(
            bore_phases.min(), bore_phases.max()))
    
    # Unwrap phase to handle wraparound
    unwrapped_phases = robust_unwrap_phase(bore_phases)
    
    if verbose:
        print("  Unwrapped phase range: {:.2f} to {:.2f} deg".format(
            unwrapped_phases.min(), unwrapped_phases.max()))
    
    # Apply median filtering if requested
    if median_window is not None:
        unwrapped_phases, n_filtered = median_filter_phases(
            unwrapped_phases, window_size=median_window, 
            threshold=median_threshold, verbose=verbose
        )
        if verbose:
            print("  Filtered phase range: {:.2f} to {:.2f} deg".format(
                unwrapped_phases.min(), unwrapped_phases.max()))
    
    # Reference to first measurement
    ref_phase = unwrapped_phases[0]
    drift_at_bore = unwrapped_phases - ref_phase
    
    if verbose:
        print("  Drift range: {:.2f} to {:.2f} deg".format(
            drift_at_bore.min(), drift_at_bore.max()))
    
    if len(bore_times) < 2:
        print("  Warning: Need at least 2 points for interpolation")
        print("  Using constant (zero) drift")
        phase_drift = np.zeros_like(timestamps)
        return phase_drift, ref_phase, unwrapped_phases
    
    if use_linear:
        # Piecewise linear interpolation
        phase_drift = np.interp(timestamps, bore_times, drift_at_bore)
        
        if verbose:
            print("  Using piecewise linear interpolation between {} boresight points".format(
                len(bore_times)))
    else:
        # Spline fit
        if len(bore_times) < 3:
            print("  Warning: Need at least 3 points for spline fit")
            print("  Falling back to linear interpolation")
            phase_drift = np.interp(timestamps, bore_times, drift_at_bore)
        else:
            if smoothing is None:
                smoothing = len(bore_times) * 4.0
            
            try:
                spline = UnivariateSpline(bore_times, drift_at_bore, s=smoothing, k=3)
                phase_drift = spline(timestamps)
                
                if verbose:
                    fit_at_bore = spline(bore_times)
                    residuals = drift_at_bore - fit_at_bore
                    print("  Spline fit RMS residual: {:.2f} deg".format(np.std(residuals)))
            except Exception as e:
                print("  Warning: Spline fit failed ({}), using linear interpolation".format(str(e)))
                phase_drift = np.interp(timestamps, bore_times, drift_at_bore)
    
    if verbose:
        print("  Phase drift correction range: {:.2f} to {:.2f} deg".format(
            phase_drift.min(), phase_drift.max()))
        # Show correction at different time percentiles
        print("  Correction at time percentiles:")
        for pct in [0, 25, 50, 75, 100]:
            idx = int(len(phase_drift) * pct / 100)
            if idx >= len(phase_drift):
                idx = len(phase_drift) - 1
            t = timestamps[idx] if len(timestamps) > idx else 0
            print(f"    {pct:3d}%: t={t:.1f}, correction={phase_drift[idx]:+.2f} deg")
    
    return phase_drift, ref_phase, unwrapped_phases


# Keep old name as alias for compatibility
def fit_phase_drift_spline(bore_times, bore_phases, timestamps, smoothing=None, verbose=False):
    """Alias for fit_phase_drift with spline method (for compatibility)."""
    return fit_phase_drift(bore_times, bore_phases, timestamps, smoothing=smoothing, 
                          use_linear=False, verbose=verbose)


def apply_phase_correction(phase, phase_drift, verbose=False):
    """Apply phase drift correction."""
    if verbose:
        print("Applying phase correction...")
    
    cal_phase = phase - phase_drift
    
    # Wrap back to [-180, 180]
    cal_phase = np.mod(cal_phase + 180, 360) - 180
    
    return cal_phase


def apply_complex_correction(amplitude, phase, amp_drift, phase_drift, verbose=False):
    """
    Apply complex (amplitude + phase) drift correction.
    
    The correction is:
        S_corrected = S_measured × (A_ref / A_bore(t)) × exp(-i × phase_drift(t))
    
    Which gives:
        A_corrected = A_measured × (A_ref / A_bore(t))
        φ_corrected = φ_measured - phase_drift(t)
    
    Where amp_drift = A_bore(t) / A_ref - 1, so A_bore(t) / A_ref = 1 + amp_drift
    
    Args:
        amplitude: Array of measured amplitudes
        phase: Array of measured phases (degrees)
        amp_drift: Array of amplitude drift values (as fraction of reference)
        phase_drift: Array of phase drift values (degrees)
        verbose: Print progress
        
    Returns:
        cal_amplitude: Corrected amplitudes
        cal_phase: Corrected phases (degrees)
    """
    if verbose:
        print("Applying COMPLEX correction (amplitude + phase)...")
        print(f"  Amplitude drift range: {amp_drift.min():.4f} to {amp_drift.max():.4f}")
        print(f"  Phase drift range: {phase_drift.min():.2f} to {phase_drift.max():.2f} deg")
    
    # Amplitude correction: divide by the relative amplitude drift
    # amp_drift = (A_bore / A_ref) - 1, so A_bore / A_ref = 1 + amp_drift
    # A_corrected = A_measured / (A_bore / A_ref) = A_measured / (1 + amp_drift)
    cal_amplitude = amplitude / (1.0 + amp_drift)
    
    # Phase correction (same as before)
    cal_phase = phase - phase_drift
    
    # Wrap phase back to [-180, 180]
    cal_phase = np.mod(cal_phase + 180, 360) - 180
    
    if verbose:
        amp_correction_pct = (amplitude / cal_amplitude - 1) * 100
        print(f"  Amplitude correction range: {amp_correction_pct.min():.2f}% to {amp_correction_pct.max():.2f}%")
    
    return cal_amplitude, cal_phase


def fit_amplitude_drift(bore_times, bore_amps, timestamps, smoothing=None, verbose=False):
    """
    Fit amplitude drift using spline interpolation.
    
    Returns the relative amplitude drift: (A_bore(t) / A_ref) - 1
    So that A_bore(t) / A_ref = 1 + amp_drift
    """
    if verbose:
        print("Fitting amplitude drift...")
        print(f"  Raw amplitude range: {bore_amps.min():.6f} to {bore_amps.max():.6f}")
    
    ref_amp = bore_amps[0]
    amp_ratio = bore_amps / ref_amp  # A_bore / A_ref
    
    if verbose:
        print(f"  Reference amplitude: {ref_amp:.6f}")
        print(f"  Amplitude ratio range: {amp_ratio.min():.4f} to {amp_ratio.max():.4f}")
    
    if len(bore_times) < 2:
        print("  Warning: Need at least 2 points for interpolation")
        amp_drift = np.zeros_like(timestamps)
        return amp_drift, ref_amp
    
    if len(bore_times) < 3:
        # Linear interpolation
        amp_drift = np.interp(timestamps, bore_times, amp_ratio - 1)
    else:
        # Spline fit
        if smoothing is None:
            smoothing = len(bore_times) * 4.0
        
        try:
            spline = UnivariateSpline(bore_times, amp_ratio - 1, s=smoothing, k=3)
            amp_drift = spline(timestamps)
            
            if verbose:
                fit_at_bore = spline(bore_times)
                residuals = (amp_ratio - 1) - fit_at_bore
                print(f"  Spline fit RMS residual: {np.std(residuals):.6f}")
        except Exception as e:
            print(f"  Warning: Spline fit failed ({str(e)}), using linear interpolation")
            amp_drift = np.interp(timestamps, bore_times, amp_ratio - 1)
    
    if verbose:
        print(f"  Amplitude drift range: {amp_drift.min():.4f} to {amp_drift.max():.4f}")
    
    return amp_drift, ref_amp


def write_output_data(filepath, timestamps, az, el, amplitude, phase, 
                          header=None, verbose=False, mask=None, calibrated=True):
    """Write data to file."""
    if verbose:
        if calibrated:
            print("Writing calibrated data to " + str(filepath) + "...")
        else:
            print("Writing data (uncalibrated) to " + str(filepath) + "...")
    
    n_written = 0
    with open(filepath, "w") as f:
        if header:
            f.write(header + "\n")
        
        for i in range(len(timestamps)):
            if mask is not None and not mask[i]:
                continue
            f.write("{:15.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(
                timestamps[i], az[i], el[i], amplitude[i], phase[i]))
            n_written += 1
    
    if verbose:
        print("  Wrote " + str(n_written) + " data points")
    
    return n_written


def main():
    args = parse_args()
    
    # Determine save prefix
    if args.save_prefix:
        prefix = args.save_prefix
    else:
        # Extract prefix from input filename
        import os
        basename = os.path.basename(args.input_file)
        # Remove common suffixes
        for suffix in ['_trimmed', '_raw', '.txt', '_with_time']:
            if basename.endswith(suffix):
                basename = basename[:-len(suffix)]
            elif suffix in basename:
                basename = basename.replace(suffix, '')
        prefix = basename.replace('.txt', '')
    
    timestamps, az, el, amplitude, phase, header = read_raw_data(
        args.input_file, verbose=args.verbose
    )
    
    if len(timestamps) == 0:
        print("Error: No data found in input file")
        sys.exit(1)
    
    # Identify boresight segments
    all_segments = identify_boresight_segments(
        timestamps, az, el, 
        args.center_az, args.center_el,
        args.tolerance, args.min_samples,
        verbose=args.verbose
    )
    
    # Get valid boresight measurements
    bore_times, bore_amps, bore_phases, valid_segments = compute_boresight_values(
        timestamps, amplitude, phase, all_segments, args.min_bore_samples,
        min_amp=args.min_amp, verbose=args.verbose
    )
    
    # Save boresight data if requested
    if args.save_boresight_data and len(bore_times) > 0:
        bore_data_file = f"{prefix}_boresight_data.txt"
        with open(bore_data_file, 'w') as f:
            f.write("# Boresight measurements\n")
            f.write("# time(s)  amplitude  phase(deg)  phase_unwrapped(deg)\n")
            unwrapped = robust_unwrap_phase(bore_phases) if len(bore_phases) > 0 else bore_phases
            for i in range(len(bore_times)):
                f.write(f"{bore_times[i]:.3f}  {bore_amps[i]:.6f}  {bore_phases[i]:.3f}  {unwrapped[i]:.3f}\n")
        print(f"  Saved boresight data: {bore_data_file}")
    
    # Create mask to remove boresight segments (used in both modes)
    if args.remove_boresight or args.no_correction:
        if args.verbose:
            print("Creating mask to exclude boresight visits...")
        
        keep_mask = np.ones(len(timestamps), dtype=bool)
        
        for start_idx, end_idx in valid_segments:
            seg_start_time = timestamps[start_idx] - args.slew_margin
            seg_end_time = timestamps[end_idx - 1] + args.slew_margin
            
            time_mask = (timestamps >= seg_start_time) & (timestamps <= seg_end_time)
            keep_mask[time_mask] = False
        
        if args.verbose:
            print("  Removed " + str(np.sum(~keep_mask)) + " boresight points")
    else:
        keep_mask = None
    
    # Mode: no-correction (just remove boresight, no calibration)
    if args.no_correction:
        if args.verbose:
            print("\n*** NO-CORRECTION MODE: Skipping phase calibration ***")
        
        n_written = write_output_data(
            args.output, timestamps, az, el, amplitude, phase,
            header, verbose=args.verbose, mask=keep_mask, calibrated=False
        )
        
        print("")
        print("Output (no calibration applied):")
        print("  Input:  " + args.input_file)
        print("  Output: " + args.output)
        print("  Boresight points removed: " + str(np.sum(~keep_mask)))
        print("  Data points written: " + str(n_written))
        sys.exit(0)
    
    # Normal calibration mode
    if len(bore_times) < 2:
        print("Warning: Insufficient valid boresight segments for calibration")
        print("         Copying input to output without calibration")
        write_output_data(args.output, timestamps, az, el, amplitude, phase,
                              header, verbose=args.verbose, mask=keep_mask, calibrated=False)
        sys.exit(0)
    
    # Fit phase drift using selected method
    phase_drift, ref_phase, unwrapped_bore_phases = fit_phase_drift(
        bore_times, bore_phases, timestamps, 
        smoothing=args.smoothing, use_linear=args.linear,
        median_window=args.median_window, median_threshold=args.median_threshold,
        verbose=args.verbose
    )
    
    # Fit amplitude drift
    amp_drift, ref_amp = fit_amplitude_drift(
        bore_times, bore_amps, timestamps,
        smoothing=args.smoothing, verbose=args.verbose
    )
    
    # Apply COMPLEX correction (amplitude + phase)
    cal_amplitude, cal_phase = apply_complex_correction(
        amplitude, phase, amp_drift, phase_drift, verbose=args.verbose
    )
    
    # Write output
    n_written = write_output_data(
        args.output, timestamps, az, el, cal_amplitude, cal_phase,
        header, verbose=args.verbose, mask=keep_mask, calibrated=True
    )
    
    # Summary
    method_name = "piecewise linear" if args.linear else "spline"
    print("")
    print("COMPLEX Calibration complete:")
    print("  Input:  " + args.input_file)
    print("  Output: " + args.output)
    print("  Valid boresight segments: " + str(len(valid_segments)))
    print("  Method: " + method_name)
    print("  Phase drift corrected: {:.2f} deg".format(
        phase_drift.max() - phase_drift.min()))
    print("  Amplitude drift corrected: {:.4f} (relative)".format(
        amp_drift.max() - amp_drift.min()))
    if args.remove_boresight:
        print("  Boresight points removed: " + str(np.sum(~keep_mask)))
    print("  Data points written: " + str(n_written))
    
    # Optional diagnostic plot
    if args.plot:
        try:
            import matplotlib.pyplot as plt
            
            fig, axes = plt.subplots(3, 1, figsize=(12, 10))
            
            # Use unwrapped phases for plotting
            drift_bore = unwrapped_bore_phases - unwrapped_bore_phases[0]
            
            # Panel 1: Phase drift with fit
            ax = axes[0]
            
            # If median filtering was applied, show raw phases too
            if args.median_window is not None:
                raw_unwrapped = robust_unwrap_phase(bore_phases)
                ax.scatter(bore_times, raw_unwrapped, c='lightgray', s=30, zorder=3, 
                           label='Raw measurements', alpha=0.7)
                ax.scatter(bore_times, unwrapped_bore_phases, c='red', s=50, zorder=5, 
                           label='After median filter')
            else:
                ax.scatter(bore_times, unwrapped_bore_phases, c='red', s=50, zorder=5, 
                           label='Boresight measurements (unwrapped)')
            
            t_plot = np.linspace(bore_times.min(), bore_times.max(), 500)
            if args.linear:
                # Show piecewise linear interpolation
                fit_plot = np.interp(t_plot, bore_times, drift_bore) + unwrapped_bore_phases[0]
                ax.plot(t_plot, fit_plot, 'b-', lw=2, label='Piecewise linear fit')
            elif len(bore_times) >= 3:
                # Show spline fit
                s = args.smoothing if args.smoothing else len(bore_times)*4
                spline = UnivariateSpline(bore_times, drift_bore, s=s, k=3)
                ax.plot(t_plot, spline(t_plot) + unwrapped_bore_phases[0], 'b-', lw=2, label='Spline fit')
            
            ax.set_xlabel('Time (seconds)')
            ax.set_ylabel('Phase (degrees)')
            title = 'Boresight Phase vs Time ({})'.format(method_name)
            if args.median_window:
                title += ' [median filtered, window={}]'.format(args.median_window)
            ax.set_title(title)
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            # Panel 2: Amplitude drift with fit
            ax = axes[1]
            amp_ratio_bore = bore_amps / bore_amps[0]
            ax.scatter(bore_times, amp_ratio_bore, c='green', s=50, zorder=5, 
                       label='Boresight amplitude (relative)')
            
            # Show spline fit for amplitude
            if len(bore_times) >= 3:
                s = args.smoothing if args.smoothing else len(bore_times)*4
                try:
                    amp_spline = UnivariateSpline(bore_times, amp_ratio_bore, s=s, k=3)
                    ax.plot(t_plot, amp_spline(t_plot), 'g-', lw=2, label='Spline fit')
                except:
                    pass
            
            ax.set_xlabel('Time (seconds)')
            ax.set_ylabel('Relative amplitude')
            ax.set_title('Boresight Amplitude vs Time (normalized to first)')
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.axhline(y=1.0, color='k', linestyle='--', alpha=0.3)
            
            # Panel 3: Applied corrections
            ax = axes[2]
            step = max(1, len(timestamps) // 2000)
            ax.plot(timestamps[::step], phase_drift[::step], 'b-', lw=1, alpha=0.7, label='Phase correction (deg)')
            ax.set_xlabel('Time (seconds)')
            ax.set_ylabel('Phase correction (degrees)', color='blue')
            ax.tick_params(axis='y', labelcolor='blue')
            ax.grid(True, alpha=0.3)
            
            # Secondary y-axis for amplitude correction
            ax2 = ax.twinx()
            ax2.plot(timestamps[::step], (1 + amp_drift[::step]) * 100, 'g-', lw=1, alpha=0.7, label='Amplitude factor (%)')
            ax2.set_ylabel('Amplitude factor (%)', color='green')
            ax2.tick_params(axis='y', labelcolor='green')
            ax2.axhline(y=100, color='green', linestyle='--', alpha=0.3)
            
            ax.set_title('Applied COMPLEX Correction ({})'.format(method_name))
            
            plt.tight_layout()
            plot_filename = f"{prefix}_boresight_cal_complex_diagnostic.png"
            plt.savefig(plot_filename, dpi=150)
            print(f"  Saved diagnostic plot: {plot_filename}")
            plt.show()
        except ImportError:
            print("  (matplotlib not available for plotting)")


if __name__ == "__main__":
    main()
