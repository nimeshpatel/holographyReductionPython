#!/usr/bin/env python3
"""
fix_missing_cell.py - Fix missing grid cells in rgin.dat

Adds missing cells by interpolating from neighboring cells.
"""

import numpy as np
import sys

def fix_missing_cells(input_file='rgin.dat', output_file=None, grid_size=32):
    """Fix missing cells in regridded holography data."""
    
    if output_file is None:
        output_file = input_file
    
    # Read existing data
    data = np.loadtxt(input_file)
    rows = data[:,0].astype(int)
    cols = data[:,1].astype(int)
    
    print(f"Read {len(data)} lines from {input_file}")
    
    # Find all missing cells
    expected = set((r, c) for r in range(grid_size) for c in range(grid_size))
    actual = set(zip(rows, cols))
    missing = expected - actual
    
    if not missing:
        print("No missing cells found.")
        return
    
    print(f"Found {len(missing)} missing cell(s): {sorted(missing)}")
    
    # Add missing cells with interpolated values
    new_rows = []
    for (mr, mc) in missing:
        # Get neighboring cells for interpolation
        neighbors = []
        for dr, dc in [(-1,0), (1,0), (0,-1), (0,1), (-1,-1), (-1,1), (1,-1), (1,1)]:
            nr, nc = mr + dr, mc + dc
            mask = (rows == nr) & (cols == nc)
            if any(mask):
                neighbors.append(data[mask][0])
        
        if neighbors:
            neighbors = np.array(neighbors)
            avg_amp = np.mean(neighbors[:, 2])
            avg_phase = np.mean(neighbors[:, 3])
            print(f"  Cell ({mr}, {mc}): interpolated from {len(neighbors)} neighbors")
            print(f"    amp = {avg_amp:.6e}, phase = {avg_phase:.6f}")
        else:
            avg_amp, avg_phase = 0.0, 0.0
            print(f"  Cell ({mr}, {mc}): no neighbors, using zeros")
        
        new_rows.append([mr, mc, avg_amp, avg_phase])
    
    # Append new rows
    new_rows = np.array(new_rows)
    data = np.vstack([data, new_rows])
    
    # Sort by row, then column
    idx = np.lexsort((data[:,1], data[:,0]))
    data = data[idx]
    
    # Write back
    with open(output_file, 'w') as f:
        for row in data:
            r = int(row[0])
            c = int(row[1])
            amp = row[2]
            phase = row[3]
            f.write(f"{r} {c} {amp:.16e} {phase:.16e}\n")
    
    print(f"Wrote {len(data)} lines to {output_file}")


if __name__ == '__main__':
    if len(sys.argv) > 1:
        grid_size = int(sys.argv[1])
    else:
        grid_size = 32
    
    fix_missing_cells(grid_size=grid_size)
