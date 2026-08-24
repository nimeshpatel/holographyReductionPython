#!/usr/bin/env python3
"""
plot_boresight_phase.py - Quick visualization of boresight phase drift

This script extracts and plots boresight phase measurements from raw holography
data without running the full reduction pipeline. Useful for quickly assessing
whether boresight calibration is likely to help.

Usage:
    python plot_boresight_phase.py <raw_data_file> [options]

Examples:
    python plot_boresight_phase.py ../data/holoADC-AzEl-20260522_123456.txt
    python plot_boresight_phase.py ../data/holoADC-AzEl-20260522_123456.txt --center-az 230.18 --center-el 2.87
    python plot_boresight_phase.py ../data/holoADC-AzEl-20260522_123456.txt --output boresight_phase.png
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys


def parse_timestamp(ts_str):
    """Parse timestamp string to float seconds."""
    ts_str = str(ts_str).strip()
    if ':' in ts_str:
        parts = ts_str.split(':')
        hours = float(parts[0])
        minutes = float(parts[1])
        seconds = float(parts[2]) if len(parts) > 2 else 0.0
        return hours * 3600 + minutes * 60 + seconds
    else:
        return float(ts_str)


def read_raw_data(filename, verbose=False):
    """Read raw holography data file."""
    if verbose:
        print(f"Reading {filename}...")
    
    timestamps = []
    az = []
    el = []
    amplitude = []
    phase = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('!'):
                continue
            
            parts = line.split()
            if len(parts) >= 5:
                # Format: timestamp az el amp phase
                timestamps.append(parse_timestamp(parts[0]))
                az.append(float(parts[1]))
                el.append(float(parts[2]))
                amplitude.append(float(parts[3]))
                phase.append(float(parts[4]))
            elif len(parts) == 4:
                # Format: az el amp phase (no timestamp)
                az.append(float(parts[0]))
                el.append(float(parts[1]))
                amplitude.append(float(parts[2]))
                phase.append(float(parts[3]))
    
    if verbose:
        print(f"  Read {len(az)} data points")
    
    has_timestamps = len(timestamps) == len(az)
    if not has_timestamps:
        # Create synthetic timestamps based on sample index
        timestamps = list(range(len(az)))
    
    return (np.array(timestamps), np.array(az), np.array(el), 
            np.array(amplitude), np.array(phase), has_timestamps)


def identify_boresight_segments(timestamps, az, el, center_az, center_el, 
                                 tolerance=0.01, min_samples=100, verbose=False):
    """Identify contiguous segments where telescope is at boresight position."""
    if verbose:
        print(f"Identifying boresight segments...")
        print(f"  Center: az={center_az:.4f}, el={center_el:.4f}")
        print(f"  Tolerance: {tolerance} deg")
    
    at_boresight = (np.abs(az - center_az) < tolerance) & (np.abs(el - center_el) < tolerance)
    
    segments = []
    in_segment = False
    start_idx = 0
    
    for i in range(len(at_boresight)):
        if at_boresight[i] and not in_segment:
            in_segment = True
            start_idx = i
        elif not at_boresight[i] and in_segment:
            in_segment = False
            if i - start_idx >= min_samples:
                segments.append((start_idx, i))
    
    if in_segment and len(timestamps) - start_idx >= min_samples:
        segments.append((start_idx, len(timestamps)))
    
    if verbose:
        print(f"  Found {len(segments)} segments with >= {min_samples} samples")
    
    return segments


def compute_boresight_values(timestamps, amplitude, phase, segments, 
                              min_bore_samples=5000, min_amp=0.3, verbose=False):
    """Compute mean amplitude and phase for valid boresight segments."""
    if verbose:
        print("Analyzing boresight segments...")
    
    bore_times = []
    bore_amps = []
    bore_phases = []
    
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
            print(f"  Seg {i+1:3d}: t={mean_time:.1f}, amp={mean_amp:.3f}, "
                  f"phase={mean_phase:7.2f} deg, n={n_samples:5d} [{status}]")
        
        if is_valid:
            bore_times.append(mean_time)
            bore_amps.append(mean_amp)
            bore_phases.append(mean_phase)
    
    if verbose:
        print(f"  Valid segments: {len(bore_times)}")
    
    return np.array(bore_times), np.array(bore_amps), np.array(bore_phases)


def robust_unwrap_phase(phases):
    """Robustly unwrap a sequence of phase values."""
    if len(phases) == 0:
        return phases.copy()
    
    unwrapped = np.zeros_like(phases)
    unwrapped[0] = phases[0]
    
    for i in range(1, len(phases)):
        diff = phases[i] - phases[i-1]
        while diff > 180:
            diff -= 360
        while diff < -180:
            diff += 360
        unwrapped[i] = unwrapped[i-1] + diff
    
    return unwrapped


def plot_boresight_phase(bore_times, bore_phases, bore_amps, output=None, 
                         title=None, show_amp=False):
    """Create diagnostic plot of boresight phase vs time."""
    
    if len(bore_times) == 0:
        print("No valid boresight segments found!")
        return
    
    # Unwrap phases
    unwrapped = robust_unwrap_phase(bore_phases)
    
    # Calculate statistics
    phase_range = unwrapped.max() - unwrapped.min()
    phase_std = np.std(unwrapped)
    duration = bore_times.max() - bore_times.min()
    
    # Estimate surface error impact
    # wavelength at 94.5 GHz = 3.17 mm
    # surface error = wavelength / (4*pi) * phase_in_radians
    error_scaling = 299792458.0 / (94.5e9) * 1e6 / (4.0 * np.pi)  # microns per radian
    phase_range_rad = np.deg2rad(phase_range)
    surface_error_range = phase_range_rad * error_scaling
    
    # Create plot
    if show_amp:
        fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    else:
        fig, axes = plt.subplots(2, 1, figsize=(12, 6))
        axes = [axes[0], axes[1]] if hasattr(axes, '__len__') else [axes]
    
    # Panel 1: Phase vs time
    ax = axes[0]
    ax.scatter(bore_times, bore_phases, c='gray', s=30, alpha=0.5, label='Raw phase')
    ax.scatter(bore_times, unwrapped, c='red', s=50, label='Unwrapped phase')
    ax.plot(bore_times, unwrapped, 'b-', alpha=0.5, lw=1)
    ax.set_xlabel('Time (seconds)')
    ax.set_ylabel('Phase (degrees)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    if title:
        ax.set_title(title)
    else:
        ax.set_title('Boresight Phase vs Time')
    
    # Add statistics text box
    stats_text = (f"Duration: {duration:.0f} s ({duration/60:.1f} min)\n"
                  f"Phase range: {phase_range:.1f}°\n"
                  f"Phase std: {phase_std:.1f}°\n"
                  f"Surface error range: {surface_error_range:.1f} µm\n"
                  f"Valid segments: {len(bore_times)}")
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Panel 2: Phase derivative (rate of change)
    ax2 = axes[1]
    if len(bore_times) > 1:
        dt = np.diff(bore_times)
        dphi = np.diff(unwrapped)
        rate = dphi / dt * 60  # degrees per minute
        mid_times = (bore_times[:-1] + bore_times[1:]) / 2
        
        ax2.plot(mid_times, rate, 'g.-', lw=1, markersize=4)
        ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
        ax2.set_xlabel('Time (seconds)')
        ax2.set_ylabel('Phase rate (deg/min)')
        ax2.set_title('Phase Drift Rate')
        ax2.grid(True, alpha=0.3)
        
        # Add rate statistics
        rate_text = f"Mean rate: {np.mean(rate):.2f} deg/min\nMax rate: {np.max(np.abs(rate)):.2f} deg/min"
        ax2.text(0.02, 0.98, rate_text, transform=ax2.transAxes, fontsize=9,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    plt.tight_layout()
    
    if output:
        plt.savefig(output, dpi=150)
        print(f"Saved plot to {output}")
    
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Quick visualization of boresight phase drift from raw holography data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s ../data/holoADC-AzEl-20260522_123456.txt
  %(prog)s data.txt --center-az 230.18 --center-el 2.87
  %(prog)s data.txt --output phase_drift.png --verbose
        """
    )
    parser.add_argument("input_file", help="Raw holography data file")
    parser.add_argument("--center-az", type=float, default=230.179,
                        help="Beacon center azimuth in degrees (default: 230.179)")
    parser.add_argument("--center-el", type=float, default=2.8724,
                        help="Beacon center elevation in degrees (default: 2.8724)")
    parser.add_argument("--tolerance", type=float, default=0.01,
                        help="Position tolerance for boresight detection (default: 0.01)")
    parser.add_argument("--min-samples", type=int, default=100,
                        help="Minimum samples to identify a segment (default: 100)")
    parser.add_argument("--min-bore-samples", type=int, default=5000,
                        help="Minimum samples for valid boresight visit (default: 5000)")
    parser.add_argument("--min-amp", type=float, default=0.3,
                        help="Minimum amplitude for valid segment (default: 0.3)")
    parser.add_argument("--output", "-o", help="Save plot to file (PNG, PDF)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Print detailed progress")
    parser.add_argument("--title", help="Custom plot title")
    
    args = parser.parse_args()
    
    # Read data
    timestamps, az, el, amplitude, phase, has_timestamps = read_raw_data(
        args.input_file, verbose=args.verbose
    )
    
    if not has_timestamps:
        print("Warning: No timestamps found in data file")
        print("         Using sample index as time axis")
    
    # Find boresight segments
    segments = identify_boresight_segments(
        timestamps, az, el,
        args.center_az, args.center_el,
        args.tolerance, args.min_samples,
        verbose=args.verbose
    )
    
    if len(segments) == 0:
        print("\nNo boresight segments found!")
        print("Check --center-az and --center-el values, or --tolerance")
        
        # Show data range to help debug
        print(f"\nData ranges:")
        print(f"  AZ: {az.min():.4f} to {az.max():.4f}")
        print(f"  EL: {el.min():.4f} to {el.max():.4f}")
        sys.exit(1)
    
    # Compute boresight values
    bore_times, bore_amps, bore_phases = compute_boresight_values(
        timestamps, amplitude, phase, segments,
        min_bore_samples=args.min_bore_samples,
        min_amp=args.min_amp,
        verbose=args.verbose
    )
    
    if len(bore_times) == 0:
        print("\nNo valid boresight segments found!")
        print("Try lowering --min-bore-samples or --min-amp")
        sys.exit(1)
    
    # Generate title from filename if not provided
    title = args.title
    if not title:
        import os
        basename = os.path.basename(args.input_file)
        title = f"Boresight Phase: {basename}"
    
    # Plot
    plot_boresight_phase(bore_times, bore_phases, bore_amps, 
                         output=args.output, title=title)
    
    # Print summary
    unwrapped = robust_unwrap_phase(bore_phases)
    phase_range = unwrapped.max() - unwrapped.min()
    error_scaling = 299792458.0 / (94.5e9) * 1e6 / (4.0 * np.pi)
    surface_error_range = np.deg2rad(phase_range) * error_scaling
    
    print(f"\n{'='*50}")
    print("SUMMARY")
    print(f"{'='*50}")
    print(f"Valid boresight segments: {len(bore_times)}")
    print(f"Phase range: {phase_range:.1f} degrees")
    print(f"Equivalent surface error: {surface_error_range:.1f} µm")
    print()
    
    if phase_range < 20:
        print("RECOMMENDATION: Phase drift is small (<20°)")
        print("                Boresight calibration unlikely to help significantly")
    elif phase_range < 50:
        print("RECOMMENDATION: Moderate phase drift (20-50°)")
        print("                Boresight calibration may provide modest improvement")
    else:
        print("RECOMMENDATION: Large phase drift (>50°)")
        print("                Boresight calibration recommended")


if __name__ == "__main__":
    main()
