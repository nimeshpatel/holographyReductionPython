#!/usr/bin/env python3
"""
generate_mask.py - Generate aperture mask for GLT holography

Creates a mask file defining valid pixels on the dish aperture.
Masks out: outside primary edge, inside secondary obstruction, quadrupod legs.

Usage:
    python generate_mask.py 64 -o mask64.dat
    python generate_mask.py 128 -o mask128.dat
"""

import argparse
import numpy as np


def generate_mask(grid_size, outer_dia=11.9, inner_dia=0.75, quad_half_width=0.0375,
                  dish_dia=12.0, verbose=False):
    """
    Generate aperture mask for holography.
    
    Parameters:
        grid_size: Size of the NxN grid
        outer_dia: Outer diameter for masking (primary edge) in meters
        inner_dia: Inner diameter for masking (subreflector) in meters
        quad_half_width: Half-width of quadrupod legs in meters
        dish_dia: Full dish diameter in meters (for scaling)
    
    Returns:
        mask: NxN array of 0s and 1s (1 = valid pixel)
    """
    
    # Pixel spacing in meters
    # The aperture field spans the dish diameter
    pixel_size = dish_dia / grid_size
    
    if verbose:
        print("Generating {}x{} mask".format(grid_size, grid_size))
        print("  Pixel size: {:.4f} m".format(pixel_size))
        print("  Outer diameter: {} m".format(outer_dia))
        print("  Inner diameter: {} m".format(inner_dia))
        print("  Quadrupod half-width: {} m".format(quad_half_width))
    
    # Create coordinate grid centered on dish
    # Pixel centers go from -dish_dia/2 to +dish_dia/2
    half_size = dish_dia / 2.0
    x = np.linspace(-half_size + pixel_size/2, half_size - pixel_size/2, grid_size)
    y = np.linspace(-half_size + pixel_size/2, half_size - pixel_size/2, grid_size)
    X, Y = np.meshgrid(x, y)
    
    # Radial distance from center
    R = np.sqrt(X**2 + Y**2)
    
    # Start with all pixels valid
    mask = np.ones((grid_size, grid_size), dtype=int)
    
    # Mask outside primary edge
    outer_radius = outer_dia / 2.0
    mask[R > outer_radius] = 0
    
    # Mask inside secondary obstruction
    inner_radius = inner_dia / 2.0
    mask[R < inner_radius] = 0
    
    # Mask quadrupod legs (4 legs at 45 degree angles)
    # Legs run diagonally from center to edge
    # A point is on a leg if it's close to one of the diagonal lines
    
    # For 45-degree legs: |x - y| < width or |x + y| < width
    # But we need to account for the leg width perpendicular to the leg direction
    # For a 45-degree line, perpendicular distance = |x ± y| / sqrt(2)
    
    leg_width = quad_half_width * 2  # Full width
    perp_dist_threshold = leg_width / np.sqrt(2) * 2  # Perpendicular distance threshold
    
    # Leg 1: y = x (upper right to lower left)
    leg1 = np.abs(X - Y) < perp_dist_threshold
    # Leg 2: y = -x (upper left to lower right)  
    leg2 = np.abs(X + Y) < perp_dist_threshold
    
    # Only mask the legs outside the secondary
    legs = (leg1 | leg2) & (R > inner_radius)
    mask[legs] = 0
    
    n_valid = np.sum(mask)
    n_total = grid_size * grid_size
    
    if verbose:
        print("  Valid pixels: {} of {} ({:.1f}%)".format(
            n_valid, n_total, 100.0 * n_valid / n_total))
    
    return mask


def write_mask(mask, filepath, verbose=False):
    """Write mask to file (one value per line)."""
    grid_size = mask.shape[0]
    
    with open(filepath, 'w') as f:
        for i in range(grid_size):
            for j in range(grid_size):
                f.write("{}\n".format(mask[i, j]))
    
    if verbose:
        print("Wrote mask to {}".format(filepath))


def main():
    parser = argparse.ArgumentParser(description="Generate aperture mask for GLT holography")
    parser.add_argument("grid_size", type=int, help="Grid size (e.g., 32, 64, 128)")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output filename (default: mask<N>.dat)")
    parser.add_argument("--outer-dia", type=float, default=11.9,
                        help="Outer diameter for masking (default: 11.9 m)")
    parser.add_argument("--inner-dia", type=float, default=0.75,
                        help="Inner diameter for masking (default: 0.75 m)")
    parser.add_argument("--quad-width", type=float, default=0.0375,
                        help="Quadrupod half-width (default: 0.0375 m)")
    parser.add_argument("--dish-dia", type=float, default=12.0,
                        help="Full dish diameter (default: 12.0 m)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print detailed information")
    args = parser.parse_args()
    
    # Generate mask
    mask = generate_mask(
        args.grid_size,
        outer_dia=args.outer_dia,
        inner_dia=args.inner_dia,
        quad_half_width=args.quad_width,
        dish_dia=args.dish_dia,
        verbose=True
    )
    
    # Output filename
    if args.output is None:
        output_path = "mask{}.dat".format(args.grid_size)
    else:
        output_path = args.output
    
    # Write mask
    write_mask(mask, output_path, verbose=True)
    
    print("\nTo use this mask, update your withphase_aber_{}.prm file:".format(args.grid_size))
    print("  name of the output mask file                     {}".format(output_path))


if __name__ == "__main__":
    main()
