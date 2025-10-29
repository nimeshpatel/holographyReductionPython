"""
HOLIS_MASK.PY - Masking operations for holography code

This module contains functions for applying masks to aperture data,
including dish/subreflector masking and quadrupod leg masking.

Translated from FORTRAN by Python conversion 2025
"""

import numpy as np
import logging
from typing import Tuple

logger = logging.getLogger(__name__)

# Masking flag for invalid pixels
MASK_VALUE = complex(-9999.0, -9999.0)
MASK_VALUE_REAL = -9999.0


def apply_mask2(c: np.ndarray, ndim: int, ncent: int, ndish: float, nsub: float) -> np.ndarray:
    """
    Apply simple circular mask to aperture data.
    Points outside main reflector and inside subreflector are set to zero.

    Args:
        c: Complex aperture array
        ndim: Array dimension
        ncent: Center pixel
        ndish: Dish radius in pixels
        nsub: Subreflector radius in pixels

    Returns:
        Masked complex array
    """
    # Create coordinate arrays
    i_coords, j_coords = np.meshgrid(np.arange(1, ndim + 1), np.arange(1, ndim + 1), indexing='ij')

    # Calculate radius from center
    radius = np.sqrt((i_coords - ncent)**2 + (j_coords - ncent)**2)

    # Apply mask: set to zero if outside dish or inside subreflector
    mask = (radius <= ndish) & (radius >= nsub)
    c_masked = c.copy()
    c_masked[~mask] = 0.0 + 0.0j

    return c_masked


def apply_mask4(c: np.ndarray, ndim: int, ncent: int, pixbylen: float,
                out_mask: float, ain_mask: float, quad_hw: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply detailed mask including dish, subreflector, and quadrupod legs.
    Masked pixels are set to -9999-9999j.

    Args:
        c: Complex aperture array
        ndim: Array dimension
        ncent: Center pixel
        pixbylen: Pixels per unit length
        out_mask: Outer mask diameter (in physical units)
        ain_mask: Inner mask diameter (in physical units)
        quad_hw: Half-width of quadrupod legs (in physical units)

    Returns:
        Tuple of (masked complex array, mask pattern array)
    """
    # Initialize mask pattern (1=valid, 0=masked)
    mask = np.ones((ndim, ndim), dtype=np.int32)

    # Calculate mask parameters in pixels
    pi = np.pi
    theta = 90.0 * pi / 180.0  # Quadrupod orientation
    slope1 = np.tan(theta)
    slope2 = -1.0 / slope1 if slope1 != 0 else 0.0

    n_out = int(np.round(out_mask * pixbylen / 2.0))
    n_in = int(np.round(ain_mask * pixbylen / 2.0))
    n_quad_hw = int(np.round(quad_hw * pixbylen))

    # Create coordinate arrays
    i_coords, j_coords = np.meshgrid(np.arange(1, ndim + 1), np.arange(1, ndim + 1), indexing='ij')

    # Calculate x, y from center
    x = j_coords - ncent
    y = i_coords - ncent

    # Calculate radius
    radius = np.sqrt((i_coords - ncent)**2 + (j_coords - ncent)**2)

    # Apply dish and subreflector mask
    mask[(radius >= n_out) | (radius <= n_in)] = 0

    # Apply quadrupod legs mask
    # This creates a cross pattern (vertical and horizontal)
    quad_mask = (np.abs(x) <= n_quad_hw) | (np.abs(y) <= n_quad_hw)
    mask[quad_mask] = 0

    # Apply mask to complex array
    c_masked = c.copy()
    c_masked[mask == 0] = MASK_VALUE

    return c_masked, mask


def apply_mask5(ph: np.ndarray, maskfile: str) -> np.ndarray:
    """
    Apply mask pattern from file to phase array.
    Masked pixels are set to -9999.

    Args:
        ph: Phase array
        maskfile: Mask file name

    Returns:
        Masked phase array
    """
    # Read mask from file
    mask = np.loadtxt(maskfile, dtype=np.int32)

    # Reshape mask to match phase array dimensions if needed
    if mask.ndim == 1:
        mask = mask.reshape(ph.shape)

    # Apply mask to phase
    ph_masked = ph.copy()
    ph_masked[mask == 0] = MASK_VALUE_REAL

    return ph_masked


def onedim(ph: np.ndarray, mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Convert 2D phase array to 1D array including only unmasked pixels.
    Also generate index arrays for mapping back to 2D.

    Args:
        ph: Phase array (2D)
        mask: Mask array (1=valid, 0=masked)

    Returns:
        Tuple of (1D phase array, i indices, j indices, number of valid points)
    """
    ndim = ph.shape[0]

    # Also mark pixels with -9999 as masked
    mask_combined = mask.copy()
    mask_combined[ph == MASK_VALUE_REAL] = 0

    # Find valid pixels
    valid = mask_combined == 1
    nmasked = np.sum(valid)

    # Create 1D arrays
    ph1 = ph[valid]

    # Get indices of valid pixels (convert to 1-based for consistency with Fortran)
    ii_coords, ij_coords = np.where(valid)
    ii = ii_coords + 1  # Convert to 1-based indexing
    ij = ij_coords + 1

    return ph1, ii, ij, nmasked


def rms1d(ph1: np.ndarray) -> float:
    """
    Calculate RMS of 1D phase array.

    Args:
        ph1: 1D phase array

    Returns:
        RMS value in radians
    """
    if len(ph1) == 0:
        return 0.0

    # Calculate mean and mean square
    mean_ph = np.mean(ph1)
    mean_sq = np.mean(ph1**2)

    # Calculate RMS
    rms = np.sqrt(mean_sq - mean_ph**2)

    return rms
