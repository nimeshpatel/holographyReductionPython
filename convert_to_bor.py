#!/usr/bin/env python3
"""
convert_to_bor.py - Convert boresight_cal.py output to original Holis bor format

The original Holis preprocess.py expects a boresight file with N+1 lines
(one per row boundary), each containing: amplitude phase

This script converts the boresight_data.txt from boresight_cal.py to that format.

Usage:
    python convert_to_bor.py holoADC-AzEl-20260509_175932_boresight_data.txt 32 -o bor32

Input format (from boresight_cal.py --save-boresight-data):
    # time(s)  amplitude  phase(deg)  phase_unwrapped(deg)
    64796.500  0.785189  -144.200  -144.200
    ...
    
Output format (for preprocess.py with ido_bore=1):
    amplitude phase
"""

import argparse
import numpy as np
import sys


def convert_boresight(input_file, grid_size, output_file, use_unwrapped=True, 
                      negate_phase=False, verbose=False):
    """
    Convert boresight data to original Holis format.
    
    Args:
        input_file: Path to boresight_data.txt (time, amp, phase, phase_unwrapped)
        grid_size: Grid size (32, 64, or 128)
        output_file: Output bor file path
        use_unwrapped: Use unwrapped phase (column 4) instead of wrapped (column 3)
        negate_phase: Negate phases (needed for preprocess.py scale_ph compatibility)
        verbose: Print details
    """
    # Expected number of boresight measurements
    n_expected = grid_size + 1
    
    if verbose:
        print(f"Converting boresight data for {grid_size}x{grid_size} map")
        print(f"Expected measurements: {n_expected}")
    
    # Read input data (skip comment lines)
    data = []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 3:
                    data.append([float(x) for x in parts])
    
    data = np.array(data)
    n_measurements = len(data)
    
    if verbose:
        print(f"Found {n_measurements} measurements")
    
    if data.shape[1] >= 4:
        times = data[:, 0]
        amps = data[:, 1]
        phases_wrapped = data[:, 2]
        phases_unwrapped = data[:, 3]
        phases = phases_unwrapped if use_unwrapped else phases_wrapped
        if verbose:
            print(f"  Using {'unwrapped' if use_unwrapped else 'wrapped'} phases")
    else:
        times = data[:, 0]
        amps = data[:, 1]
        phases = data[:, 2]
    
    if verbose:
        print(f"  Time range: {times.min():.1f} to {times.max():.1f} s")
        print(f"  Phase range: {phases.min():.1f} to {phases.max():.1f} deg")
        print(f"  Amplitude range: {amps.min():.4f} to {amps.max():.4f}")
    
    # Check if we have enough measurements
    if n_measurements < n_expected:
        print(f"Warning: Only {n_measurements} measurements, expected {n_expected}")
        print(f"  Will pad with last measurement")
        n_to_use = n_measurements
    elif n_measurements > n_expected:
        print(f"Note: {n_measurements} measurements, expected {n_expected}")
        print(f"  Will use first {n_expected} measurements")
        n_to_use = n_expected
    else:
        n_to_use = n_expected
    
    # Select or pad to get exactly N+1 measurements
    selected_amps = np.zeros(n_expected)
    selected_phases = np.zeros(n_expected)
    
    # Use first n_to_use measurements
    selected_amps[:n_to_use] = amps[:n_to_use]
    selected_phases[:n_to_use] = phases[:n_to_use]
    
    # Pad remaining with last value if needed
    if n_to_use < n_expected:
        selected_amps[n_to_use:] = amps[-1]
        selected_phases[n_to_use:] = phases[-1]
    
    # Negate phases if requested (needed for preprocess.py compatibility)
    if negate_phase:
        if verbose:
            print(f"  Negating phases (for preprocess.py scale_ph compatibility)")
        selected_phases = -selected_phases
    
    # Write output in original format
    # Format: amplitude phase (space separated)
    with open(output_file, 'w') as f:
        for i in range(n_expected):
            f.write(f"  {selected_amps[i]:.6f}  {selected_phases[i]:.6f}\n")
    
    if verbose:
        print(f"\nWrote {n_expected} lines to {output_file}")
        print(f"  First: amp={selected_amps[0]:.4f}, phase={selected_phases[0]:.1f} deg")
        print(f"  Last:  amp={selected_amps[-1]:.4f}, phase={selected_phases[-1]:.1f} deg")
        
        # Show drift
        phase_drift = selected_phases[-1] - selected_phases[0]
        print(f"  Total phase drift: {phase_drift:.1f} deg")
    
    return n_expected


def main():
    parser = argparse.ArgumentParser(
        description="Convert boresight_cal.py output to original Holis bor format"
    )
    parser.add_argument("input_file", type=str,
                        help="Input boresight_data.txt file")
    parser.add_argument("grid_size", type=int, choices=[32, 64, 128],
                        help="Grid size (32, 64, or 128)")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output file (default: bor<size>)")
    parser.add_argument("--wrapped", action="store_true",
                        help="Use wrapped phases instead of unwrapped")
    parser.add_argument("--negate", action="store_true",
                        help="Negate phases (needed for preprocess.py scale_ph compatibility)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")
    
    args = parser.parse_args()
    
    if args.output is None:
        args.output = f"bor{args.grid_size}"
    
    convert_boresight(args.input_file, args.grid_size, args.output, 
                      use_unwrapped=not args.wrapped,
                      negate_phase=args.negate,
                      verbose=args.verbose)
    
    print(f"\nCreated {args.output}")
    print(f"\nTo use with original approach:")
    print(f"  1. In preprocess.prm, set:")
    print(f"       Do bore-sight drift correction (1/0)....         1")
    print(f"       Bore-sight data file name...............         {args.output}")
    print(f"  2. Run: python preprocess.py")


if __name__ == "__main__":
    main()
