#!/usr/bin/env python3
"""
boresight_diagnostics.py - Diagnose boresight calibration data quality

This script analyzes boresight calibration data to:
1. Plot boresight amplitude and phase vs time
2. Show statistics and identify outliers
3. Compare calibrated vs uncalibrated data distributions
4. Check for data corruption

Usage:
    python boresight_diagnostics.py trimmed_with_time.txt [calibrated.txt]
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys


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


def read_data(filepath):
    """Read holography data file."""
    print("Reading " + filepath + "...")
    
    data_lines = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            data_lines.append(line)
    
    n = len(data_lines)
    print("  Found " + str(n) + " data lines")
    
    # Detect format: 4 columns (no timestamp) or 5 columns (with timestamp)
    first_parts = data_lines[0].split()
    has_timestamp = len(first_parts) >= 5
    
    timestamps = np.zeros(n)
    az = np.zeros(n)
    el = np.zeros(n)
    amplitude = np.zeros(n)
    phase = np.zeros(n)
    
    for i, line in enumerate(data_lines):
        parts = line.split()
        if has_timestamp:
            timestamps[i] = parse_timestamp(parts[0])
            az[i] = float(parts[1])
            el[i] = float(parts[2])
            amplitude[i] = float(parts[3])
            phase[i] = float(parts[4])
        else:
            timestamps[i] = i  # Use index as pseudo-timestamp
            az[i] = float(parts[0])
            el[i] = float(parts[1])
            amplitude[i] = float(parts[2])
            phase[i] = float(parts[3])
    
    return timestamps, az, el, amplitude, phase, has_timestamp


def identify_boresight_segments(timestamps, az, el, center_az=230.179, center_el=2.8724, 
                                 tolerance=0.01, min_samples=100):
    """Identify boresight measurement segments."""
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
    
    return segments, is_boresight


def compute_segment_stats(timestamps, amplitude, phase, segments, min_samples=5000):
    """Compute statistics for each boresight segment."""
    stats = []
    
    for i, (start, end) in enumerate(segments):
        n = end - start
        t_mean = np.mean(timestamps[start:end])
        amp_mean = np.mean(amplitude[start:end])
        amp_std = np.std(amplitude[start:end])
        
        # Complex average for phase
        phase_rad = np.deg2rad(phase[start:end])
        complex_avg = np.mean(np.exp(1j * phase_rad))
        phase_mean = np.rad2deg(np.angle(complex_avg))
        phase_std = np.std(phase[start:end])
        
        is_valid = n >= min_samples and amp_mean > 0.6
        
        stats.append({
            'segment': i + 1,
            'start': start,
            'end': end,
            'n': n,
            't_mean': t_mean,
            'amp_mean': amp_mean,
            'amp_std': amp_std,
            'phase_mean': phase_mean,
            'phase_std': phase_std,
            'is_valid': is_valid
        })
    
    return stats


def plot_boresight_timeseries(timestamps, amplitude, phase, segments, stats, output_prefix="bore_diag"):
    """Plot boresight amplitude and phase vs time."""
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    
    # Extract boresight data points
    bore_times = []
    bore_amps = []
    bore_phases = []
    
    for start, end in segments:
        bore_times.extend(timestamps[start:end])
        bore_amps.extend(amplitude[start:end])
        bore_phases.extend(phase[start:end])
    
    bore_times = np.array(bore_times)
    bore_amps = np.array(bore_amps)
    bore_phases = np.array(bore_phases)
    
    # Plot amplitude
    ax = axes[0]
    ax.scatter(bore_times, bore_amps, s=0.5, alpha=0.3, c='blue')
    
    # Overlay segment means
    for s in stats:
        color = 'green' if s['is_valid'] else 'red'
        ax.axhline(y=s['amp_mean'], xmin=0, xmax=1, alpha=0.3)
        ax.scatter([s['t_mean']], [s['amp_mean']], s=100, c=color, marker='o', 
                   edgecolors='black', zorder=5)
    
    ax.set_ylabel('Amplitude')
    ax.set_title('Boresight Amplitude vs Time')
    ax.grid(True, alpha=0.3)
    
    # Plot phase
    ax = axes[1]
    ax.scatter(bore_times, bore_phases, s=0.5, alpha=0.3, c='blue')
    
    # Overlay segment means
    for s in stats:
        color = 'green' if s['is_valid'] else 'red'
        ax.scatter([s['t_mean']], [s['phase_mean']], s=100, c=color, marker='o',
                   edgecolors='black', zorder=5)
    
    ax.set_xlabel('Time (seconds)')
    ax.set_ylabel('Phase (degrees)')
    ax.set_title('Boresight Phase vs Time')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_prefix + '_timeseries.png', dpi=150)
    print("Saved: " + output_prefix + "_timeseries.png")
    plt.show()


def plot_data_histograms(amp_orig, phase_orig, amp_cal, phase_cal, output_prefix="bore_diag"):
    """Plot histograms comparing original vs calibrated data."""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Amplitude histograms
    ax = axes[0, 0]
    ax.hist(amp_orig, bins=100, alpha=0.7, label='Original', density=True)
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Density')
    ax.set_title('Original Amplitude Distribution')
    ax.legend()
    
    ax = axes[0, 1]
    ax.hist(amp_cal, bins=100, alpha=0.7, label='Calibrated', color='orange', density=True)
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Density')
    ax.set_title('Calibrated Amplitude Distribution')
    ax.legend()
    
    # Phase histograms
    ax = axes[1, 0]
    ax.hist(phase_orig, bins=100, alpha=0.7, label='Original', density=True)
    ax.set_xlabel('Phase (degrees)')
    ax.set_ylabel('Density')
    ax.set_title('Original Phase Distribution')
    ax.legend()
    
    ax = axes[1, 1]
    ax.hist(phase_cal, bins=100, alpha=0.7, label='Calibrated', color='orange', density=True)
    ax.set_xlabel('Phase (degrees)')
    ax.set_ylabel('Density')
    ax.set_title('Calibrated Phase Distribution')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_prefix + '_histograms.png', dpi=150)
    print("Saved: " + output_prefix + "_histograms.png")
    plt.show()


def detect_outliers(amplitude, phase, amp_threshold=3.0, phase_threshold=3.0):
    """Detect outliers using sigma clipping."""
    
    amp_mean = np.mean(amplitude)
    amp_std = np.std(amplitude)
    phase_mean = np.mean(phase)
    phase_std = np.std(phase)
    
    amp_outliers = np.abs(amplitude - amp_mean) > amp_threshold * amp_std
    phase_outliers = np.abs(phase - phase_mean) > phase_threshold * phase_std
    
    return amp_outliers, phase_outliers, amp_mean, amp_std, phase_mean, phase_std


def plot_outlier_map(az, el, amplitude, phase, amp_outliers, phase_outliers, output_prefix="bore_diag"):
    """Plot spatial distribution of outliers."""
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Amplitude outliers
    ax = axes[0]
    ax.scatter(az[~amp_outliers], el[~amp_outliers], s=0.1, alpha=0.1, c='blue', label='Normal')
    ax.scatter(az[amp_outliers], el[amp_outliers], s=5, c='red', label='Outliers')
    ax.set_xlabel('Azimuth (deg)')
    ax.set_ylabel('Elevation (deg)')
    ax.set_title('Amplitude Outliers (n=' + str(np.sum(amp_outliers)) + ')')
    ax.legend()
    
    # Phase outliers
    ax = axes[1]
    ax.scatter(az[~phase_outliers], el[~phase_outliers], s=0.1, alpha=0.1, c='blue', label='Normal')
    ax.scatter(az[phase_outliers], el[phase_outliers], s=5, c='red', label='Outliers')
    ax.set_xlabel('Azimuth (deg)')
    ax.set_ylabel('Elevation (deg)')
    ax.set_title('Phase Outliers (n=' + str(np.sum(phase_outliers)) + ')')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_prefix + '_outliers.png', dpi=150)
    print("Saved: " + output_prefix + "_outliers.png")
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Diagnose boresight calibration data")
    parser.add_argument("input_file", help="Input file with timestamps (trimmed_with_time.txt)")
    parser.add_argument("calibrated_file", nargs='?', default=None,
                        help="Calibrated file (calibrated.txt) for comparison")
    parser.add_argument("--center-az", type=float, default=230.179)
    parser.add_argument("--center-el", type=float, default=2.8724)
    parser.add_argument("--tolerance", type=float, default=0.01)
    parser.add_argument("--output-prefix", "-o", default="bore_diag")
    args = parser.parse_args()
    
    # Read original data
    timestamps, az, el, amplitude, phase, has_ts = read_data(args.input_file)
    
    # Identify boresight segments
    print("\nIdentifying boresight segments...")
    segments, is_bore = identify_boresight_segments(
        timestamps, az, el, args.center_az, args.center_el, args.tolerance
    )
    print("  Found " + str(len(segments)) + " segments")
    
    # Compute segment statistics
    print("\nBoresight segment statistics:")
    stats = compute_segment_stats(timestamps, amplitude, phase, segments)
    
    print("\n{:>4} {:>8} {:>10} {:>10} {:>10} {:>10} {:>8}".format(
        "Seg", "Samples", "Amp Mean", "Amp Std", "Ph Mean", "Ph Std", "Valid"))
    print("-" * 70)
    
    for s in stats:
        valid_str = "YES" if s['is_valid'] else "NO"
        print("{:>4} {:>8} {:>10.4f} {:>10.4f} {:>10.2f} {:>10.2f} {:>8}".format(
            s['segment'], s['n'], s['amp_mean'], s['amp_std'], 
            s['phase_mean'], s['phase_std'], valid_str))
    
    # Plot boresight timeseries
    print("\nPlotting boresight timeseries...")
    plot_boresight_timeseries(timestamps, amplitude, phase, segments, stats, args.output_prefix)
    
    # Detect outliers in original data
    print("\nDetecting outliers in original data...")
    amp_out, phase_out, amp_mean, amp_std, phase_mean, phase_std = detect_outliers(amplitude, phase)
    print("  Amplitude: mean={:.4f}, std={:.4f}, outliers={}".format(amp_mean, amp_std, np.sum(amp_out)))
    print("  Phase: mean={:.2f}, std={:.2f}, outliers={}".format(phase_mean, phase_std, np.sum(phase_out)))
    
    # Plot outlier map
    plot_outlier_map(az, el, amplitude, phase, amp_out, phase_out, args.output_prefix)
    
    # If calibrated file provided, compare
    if args.calibrated_file:
        print("\nReading calibrated data...")
        ts_cal, az_cal, el_cal, amp_cal, phase_cal, _ = read_data(args.calibrated_file)
        
        print("\nDetecting outliers in calibrated data...")
        amp_out_cal, phase_out_cal, amp_mean_cal, amp_std_cal, phase_mean_cal, phase_std_cal = detect_outliers(amp_cal, phase_cal)
        print("  Amplitude: mean={:.4f}, std={:.4f}, outliers={}".format(amp_mean_cal, amp_std_cal, np.sum(amp_out_cal)))
        print("  Phase: mean={:.2f}, std={:.2f}, outliers={}".format(phase_mean_cal, phase_std_cal, np.sum(phase_out_cal)))
        
        # Compare phase distributions
        print("\nPhase comparison:")
        print("  Original range: {:.2f} to {:.2f}".format(phase.min(), phase.max()))
        print("  Calibrated range: {:.2f} to {:.2f}".format(phase_cal.min(), phase_cal.max()))
        
        # Check for phase wrapping issues
        phase_diff = phase_cal.max() - phase_cal.min()
        if phase_diff > 360:
            print("  WARNING: Calibrated phase range > 360 degrees - possible wrapping issue!")
        
        # Plot histograms
        print("\nPlotting histograms...")
        plot_data_histograms(amplitude, phase, amp_cal, phase_cal, args.output_prefix)
        
        # Plot calibrated outliers
        plot_outlier_map(az_cal, el_cal, amp_cal, phase_cal, amp_out_cal, phase_out_cal, 
                        args.output_prefix + "_cal")
    
    print("\nDiagnostics complete.")


if __name__ == "__main__":
    main()
