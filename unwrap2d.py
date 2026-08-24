#!/usr/bin/env python3
"""
unwrap2d.py - 2D Phase unwrapping for holography aperture field

Handles disconnected regions (quadrants separated by quadrupod legs) by:
1. Finding all connected regions
2. Unwrapping each region independently using flood-fill
3. Stitching regions together using diagonal neighbors across the gaps

Usage:
    python unwrap2d.py -i Ep.dat -o tk.dat -d 64
    python unwrap2d.py -i Ep.dat -o tk.dat -d 64 -m quality --plot
"""

import argparse
import numpy as np
from pathlib import Path
from collections import deque
import heapq


def read_ndim_from_prm(prm_file='withphase_aber.prm'):
    """Read the grid dimension from parameter file."""
    prm_path = Path(prm_file)
    if not prm_path.exists():
        return 128
    
    try:
        with open(prm_path, 'r') as f:
            for line in f:
                if 'Size N of the N by N data file' in line:
                    value_str = line[49:].strip()
                    return int(value_str)
    except Exception:
        pass
    
    return 128


def find_connected_regions(phase_2d):
    """
    Find all connected regions of valid pixels.
    Uses 4-connectivity (up, down, left, right).
    
    Returns:
        List of regions, where each region is a list of (i, j) coordinates
    """
    n = phase_2d.shape[0]
    visited = np.zeros((n, n), dtype=bool)
    regions = []
    
    for start_i in range(n):
        for start_j in range(n):
            if phase_2d[start_i, start_j] != -9999.0 and not visited[start_i, start_j]:
                # BFS to find all connected pixels
                region = []
                queue = deque([(start_i, start_j)])
                visited[start_i, start_j] = True
                
                while queue:
                    i, j = queue.popleft()
                    region.append((i, j))
                    
                    for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                        ni, nj = i + di, j + dj
                        if 0 <= ni < n and 0 <= nj < n:
                            if phase_2d[ni, nj] != -9999.0 and not visited[ni, nj]:
                                visited[ni, nj] = True
                                queue.append((ni, nj))
                
                regions.append(region)
    
    return regions


def unwrap_region(phase_2d, region, verbose=False):
    """
    Unwrap a single connected region using flood-fill.
    
    Args:
        phase_2d: Full 2D phase array
        region: List of (i, j) coordinates in this region
    
    Returns:
        Dictionary mapping (i, j) -> unwrapped phase value
    """
    if not region:
        return {}
    
    n = phase_2d.shape[0]
    result = {}
    region_set = set(region)
    
    # Find seed: pixel closest to center of region
    center_i = np.mean([p[0] for p in region])
    center_j = np.mean([p[1] for p in region])
    
    seed = min(region, key=lambda p: (p[0] - center_i)**2 + (p[1] - center_j)**2)
    
    # Initialize seed
    result[seed] = phase_2d[seed]
    
    # BFS from seed
    queue = deque()
    processed = {seed}
    
    i0, j0 = seed
    for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        ni, nj = i0 + di, j0 + dj
        if (ni, nj) in region_set and (ni, nj) not in processed:
            queue.append((ni, nj))
            processed.add((ni, nj))
    
    while queue:
        i, j = queue.popleft()
        
        # Find processed neighbors
        neighbors = []
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            ni, nj = i + di, j + dj
            if (ni, nj) in result:
                neighbors.append(result[(ni, nj)])
        
        if not neighbors:
            result[(i, j)] = phase_2d[i, j]
        else:
            neighbor_mean = np.mean(neighbors)
            raw_phase = phase_2d[i, j]
            diff = raw_phase - neighbor_mean
            n_wraps = np.round(diff / (2 * np.pi))
            result[(i, j)] = raw_phase - n_wraps * 2 * np.pi
        
        # Add unprocessed neighbors in this region
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            ni, nj = i + di, j + dj
            if (ni, nj) in region_set and (ni, nj) not in processed:
                queue.append((ni, nj))
                processed.add((ni, nj))
    
    return result


def find_diagonal_neighbors(phase_2d, region1, region2):
    """
    Find pairs of pixels from two regions that are diagonal neighbors.
    This allows stitching across quadrupod gaps.
    
    Returns:
        List of ((i1, j1), (i2, j2)) pairs where pixels are diagonally adjacent
    """
    set1 = set(region1)
    set2 = set(region2)
    pairs = []
    
    # Check all 8 directions including diagonals
    for i1, j1 in region1:
        for di, dj in [(-1, -1), (-1, 1), (1, -1), (1, 1),  # diagonals
                       (-2, 0), (2, 0), (0, -2), (0, 2),     # 2-step cardinal
                       (-1, -2), (-1, 2), (1, -2), (1, 2),   # knight moves
                       (-2, -1), (-2, 1), (2, -1), (2, 1)]:
            i2, j2 = i1 + di, j1 + dj
            if (i2, j2) in set2:
                pairs.append(((i1, j1), (i2, j2)))
    
    return pairs


def stitch_regions(phase_2d, regions, unwrapped_regions, verbose=False):
    """
    Stitch multiple unwrapped regions together.
    
    Uses diagonal neighbors to find the best 2*pi offset between regions.
    """
    if len(regions) <= 1:
        return unwrapped_regions[0] if unwrapped_regions else {}
    
    if verbose:
        print(f"Stitching {len(regions)} regions together...")
    
    # Start with the largest region as reference
    region_sizes = [len(r) for r in regions]
    ref_idx = np.argmax(region_sizes)
    
    final_result = dict(unwrapped_regions[ref_idx])
    stitched = {ref_idx}
    
    # Iteratively stitch remaining regions
    while len(stitched) < len(regions):
        best_offset = None
        best_region_idx = None
        best_n_pairs = 0
        
        for idx in range(len(regions)):
            if idx in stitched:
                continue
            
            # Find diagonal neighbors between this region and any stitched region
            all_pairs = []
            for stitched_idx in stitched:
                pairs = find_diagonal_neighbors(phase_2d, regions[stitched_idx], regions[idx])
                all_pairs.extend([(p, stitched_idx) for p in pairs])
            
            if not all_pairs:
                continue
            
            # Calculate offset using all pairs
            offsets = []
            for (p1, p2), stitched_idx in all_pairs:
                ref_phase = final_result[p1]
                new_phase = unwrapped_regions[idx][p2]
                diff = new_phase - ref_phase
                n_wraps = np.round(diff / (2 * np.pi))
                offsets.append(-n_wraps * 2 * np.pi)
            
            # Use median offset (robust to outliers)
            offset = np.median(offsets)
            
            if len(all_pairs) > best_n_pairs:
                best_n_pairs = len(all_pairs)
                best_offset = offset
                best_region_idx = idx
        
        if best_region_idx is None:
            # No more regions can be stitched via neighbors
            # Just add remaining regions with zero offset
            if verbose:
                print("  Warning: Some regions have no diagonal neighbors, using zero offset")
            for idx in range(len(regions)):
                if idx not in stitched:
                    for pos, phase in unwrapped_regions[idx].items():
                        final_result[pos] = phase
                    stitched.add(idx)
        else:
            # Apply offset and add region
            if verbose:
                print(f"  Region {best_region_idx}: {len(regions[best_region_idx])} pixels, "
                      f"offset = {best_offset:.2f} rad ({best_offset/(2*np.pi):.1f} wraps), "
                      f"using {best_n_pairs} neighbor pairs")
            
            for pos, phase in unwrapped_regions[best_region_idx].items():
                final_result[pos] = phase + best_offset
            
            stitched.add(best_region_idx)
    
    return final_result


def unwrap_2d_multiregion(phase_2d, verbose=False):
    """
    Perform 2D phase unwrapping handling multiple disconnected regions.
    
    1. Find all connected regions
    2. Unwrap each region independently
    3. Stitch regions together using diagonal neighbors
    """
    n = phase_2d.shape[0]
    
    # Find connected regions
    regions = find_connected_regions(phase_2d)
    
    if verbose:
        print(f"Found {len(regions)} connected regions:")
        for i, r in enumerate(regions):
            print(f"  Region {i}: {len(r)} pixels")
    
    if len(regions) == 0:
        return phase_2d
    
    # Unwrap each region independently
    unwrapped_regions = []
    for i, region in enumerate(regions):
        if verbose:
            print(f"Unwrapping region {i}...")
        unwrapped = unwrap_region(phase_2d, region, verbose=verbose)
        unwrapped_regions.append(unwrapped)
    
    # Stitch regions together
    final_dict = stitch_regions(phase_2d, regions, unwrapped_regions, verbose=verbose)
    
    # Convert back to 2D array
    result = np.full_like(phase_2d, -9999.0)
    for (i, j), phase in final_dict.items():
        result[i, j] = phase
    
    if verbose:
        valid_count = np.sum(phase_2d != -9999.0)
        unwrapped_count = np.sum(result != -9999.0)
        print(f"Unwrapped {unwrapped_count} of {valid_count} valid pixels")
    
    return result


def unwrap_2d_floodfill(phase_2d, verbose=False):
    """
    Perform 2D phase unwrapping using flood-fill (BFS) algorithm.
    Now uses multi-region handling for disconnected areas.
    """
    return unwrap_2d_multiregion(phase_2d, verbose=verbose)


def unwrap_2d_quality_guided(phase_2d, verbose=False):
    """
    Quality-guided 2D phase unwrapping.
    Now uses multi-region handling for disconnected areas.
    """
    # For now, use the same multi-region approach
    # Quality guidance could be added within each region if needed
    return unwrap_2d_multiregion(phase_2d, verbose=verbose)


def main():
    parser = argparse.ArgumentParser(
        description='2D Phase unwrapping for holography data'
    )
    parser.add_argument('-i', '--input', default='Ep.dat',
                        help='Input file (default: Ep.dat)')
    parser.add_argument('-o', '--output', default='tk.dat',
                        help='Output file (default: tk.dat)')
    parser.add_argument('-d', '--dim', type=int, default=None,
                        help='Grid dimension (default: read from withphase_aber.prm)')
    parser.add_argument('-m', '--method', choices=['flood', 'quality'], default='quality',
                        help='Unwrapping method: flood (BFS) or quality (quality-guided)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print detailed progress')
    parser.add_argument('--plot', action='store_true',
                        help='Show before/after plots')
    
    args = parser.parse_args()
    
    if args.dim is None:
        args.dim = read_ndim_from_prm()
    
    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    print(f"Dimension: {args.dim}")
    print(f"Method: {args.method}")
    
    data = np.loadtxt(args.input)
    n = args.dim
    
    if len(data) != n * n:
        print(f"Warning: Expected {n*n} values, got {len(data)}")
        n = int(np.sqrt(len(data)))
        print(f"Using detected dimension: {n}")
    
    phase_2d = data.reshape(n, n)
    
    valid_before = phase_2d[phase_2d != -9999]
    range_before = valid_before.max() - valid_before.min()
    
    if args.method == 'flood':
        unwrapped = unwrap_2d_floodfill(phase_2d, verbose=args.verbose)
    else:
        unwrapped = unwrap_2d_quality_guided(phase_2d, verbose=args.verbose)
    
    valid_after = unwrapped[unwrapped != -9999]
    range_after = valid_after.max() - valid_after.min()
    
    print(f"\nPhase range before: {range_before:.2f} rad ({range_before/(2*np.pi):.1f} wraps)")
    print(f"Phase range after:  {range_after:.2f} rad ({range_after/(2*np.pi):.1f} wraps)")
    
    with open(args.output, 'w') as f:
        for j in range(n):
            for i in range(n):
                if unwrapped[j, i] == -9999.0:
                    f.write("-9999.0000\n")
                else:
                    f.write(f"{unwrapped[j, i]:.6f}\n")
    
    print(f"\nWrote unwrapped phase to {args.output}")
    
    if args.plot:
        try:
            import matplotlib.pyplot as plt
            
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            
            ax = axes[0]
            phase_plot = np.where(phase_2d == -9999, np.nan, phase_2d)
            im = ax.imshow(phase_plot, cmap='RdBu_r')
            ax.set_title('Before Unwrapping')
            plt.colorbar(im, ax=ax, label='Phase (rad)')
            
            ax = axes[1]
            unwrap_plot = np.where(unwrapped == -9999, np.nan, unwrapped)
            im = ax.imshow(unwrap_plot, cmap='RdBu_r')
            ax.set_title(f'After 2D Unwrapping ({args.method})')
            plt.colorbar(im, ax=ax, label='Phase (rad)')
            
            plt.tight_layout()
            plt.savefig('unwrap_comparison.png', dpi=150)
            print("Saved: unwrap_comparison.png")
            plt.show()
        except ImportError:
            print("(matplotlib not available for plotting)")


if __name__ == "__main__":
    main()
