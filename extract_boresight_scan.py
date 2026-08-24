#!/usr/bin/env python3
"""
Extract the middle scan at boresight elevation from holography map data.

Input: holoADC-AzEl-YYYMMDD_HHMMSS.txt
Output: holoADC-AzEl-YYYMMDD_HHMMSS_boresight.txt

The program finds the scan row that passes through the boresight elevation
and writes it to a separate file.
"""

import sys
import numpy as np
from pathlib import Path


def find_map_start(data, boresight_el, map_range=0.729):
    """
    Find where the actual map data starts (after initial positioning).

    The map should start from the bottom elevation (boresight - map_range)
    and scan left to right in azimuth.
    """
    min_el = boresight_el - map_range

    # Find first index where elevation is near the minimum
    for i, row in enumerate(data):
        if abs(row[2] - min_el) < 0.05:  # Within 0.05 deg tolerance
            return i

    # If not found by elevation, look for where elevation starts changing
    # significantly from the initial positioning values
    elevations = data[:, 2]
    for i in range(100, len(elevations) - 100):
        # Check if we're at a stable elevation (scanning in azimuth)
        el_window = elevations[i:i+50]
        if np.std(el_window) < 0.01 and abs(np.mean(el_window) - boresight_el) < map_range:
            return i

    # Default: skip first 1000 lines as mentioned
    return 1000


def extract_boresight_scan(filename, boresight_el=2.8724, map_range=0.729):
    """
    Extract the scan at boresight elevation from holography map data.

    Parameters:
    -----------
    filename : str
        Input filename (holoADC-AzEl-YYYMMDD_HHMMSS.txt)
    boresight_el : float
        Boresight elevation in degrees
    map_range : float
        Map extent in degrees (±range around boresight)
    """

    # Read the data file
    print(f"Reading {filename}...")
    data = []
    header_lines = []

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                # Parse data line: time, az, el, amp, phase
                values = line.split()
                if len(values) == 5:
                    data.append([float(v) for v in values])

    data = np.array(data)
    print(f"Read {len(data)} data points")

    # Find where the map actually starts
    map_start = find_map_start(data, boresight_el, map_range)
    print(f"Map data starts at line {map_start}")

    # Get elevations from the map data
    map_data = data[map_start:]
    elevations = map_data[:, 2]

    # Find unique elevation values (rows)
    # Round elevations to 0.001 deg precision to group them efficiently
    rounded_elevations = np.round(elevations / 0.001) * 0.001
    unique_elevations = sorted(np.unique(rounded_elevations))
    print(f"Found {len(unique_elevations)} elevation rows")
    print(f"Elevation range: {min(unique_elevations):.4f} to {max(unique_elevations):.4f} deg")

    # Find the elevation row closest to boresight
    closest_el = min(unique_elevations, key=lambda x: abs(x - boresight_el))
    print(f"Boresight elevation: {boresight_el:.4f} deg")
    print(f"Closest row elevation: {closest_el:.4f} deg")

    # Extract only ONE continuous scan at this elevation
    # Find the first occurrence of boresight elevation
    tolerance_el = 0.003  # 3 millidegrees tolerance for elevation
    boresight_data = []
    in_scan = False
    scan_start_idx = None

    for i, row in enumerate(map_data):
        el = row[2]
        at_boresight = abs(el - closest_el) < tolerance_el

        if at_boresight and not in_scan:
            # Start of the boresight scan
            in_scan = True
            scan_start_idx = i
            boresight_data.append(row)
        elif at_boresight and in_scan:
            # Continue the scan
            boresight_data.append(row)
        elif not at_boresight and in_scan:
            # We've moved away from boresight elevation, stop extracting
            break

    boresight_data = np.array(boresight_data)
    print(f"Extracted {len(boresight_data)} points from one continuous scan at boresight elevation")
    if scan_start_idx is not None:
        print(f"Scan starts at index {scan_start_idx} in map data")

    # Create output filename
    input_path = Path(filename)
    output_filename = input_path.stem + "_boresight" + input_path.suffix
    output_path = input_path.parent / output_filename

    # Write output file
    print(f"Writing to {output_path}...")
    with open(output_path, 'w') as f:
        # Write header
        for header_line in header_lines:
            f.write(header_line)

        # Write boresight scan data
        for row in boresight_data:
            f.write(f"{row[0]:15.6f} {row[1]:15.6f} {row[2]:15.6f} {row[3]:15.6f} {row[4]:15.6f}\n")

    print(f"Done! Boresight scan saved to {output_filename}")

    return output_path


def main():
    if len(sys.argv) < 2:
        print("Usage: python extract_boresight_scan.py <input_file> [boresight_elevation]")
        print("Example: python extract_boresight_scan.py holoADC-AzEl-20231215_143022.txt 2.8724")
        sys.exit(1)

    input_file = sys.argv[1]
    boresight_el = float(sys.argv[2]) if len(sys.argv) > 2 else 2.8724

    extract_boresight_scan(input_file, boresight_el)


if __name__ == "__main__":
    main()
