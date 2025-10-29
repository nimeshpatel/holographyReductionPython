"""
HOLIS_MISELL.PY - Misell algorithm operations for holography code

This module contains functions specific to the Misell phase retrieval algorithm,
including initial guess generation and convergence checking.

Translated from FORTRAN by Python conversion 2025
"""

import numpy as np
import logging
from typing import Tuple

logger = logging.getLogger(__name__)


def create_guess(ndim: int, ncent: int, ndish: float, nsub: float,
                 taper: float, radphase: float, seed: int = None) -> np.ndarray:
    """
    Create initial guess of complex aperture pattern for Misell algorithm.
    Uses specified edge taper and random phase fluctuation.

    Args:
        ndim: Array dimension
        ncent: Center pixel
        ndish: Dish radius in pixels
        nsub: Subreflector radius in pixels
        taper: Edge taper in dB (0 = uniform)
        radphase: Random phase fluctuation amplitude in radians
        seed: Random seed (optional)

    Returns:
        Initial guess complex aperture array
    """
    if seed is not None:
        np.random.seed(seed)

    c = np.zeros((ndim, ndim), dtype=np.complex128)

    aloge = 0.4342944  # log10(e)
    ndish2 = ndish * ndish
    nsub2 = nsub * nsub

    # Calculate taper sigma if taper is specified
    if taper != 0.0:
        sig2 = 10.0 * ndish2 * aloge / taper

    # Create coordinate arrays
    i_coords, j_coords = np.meshgrid(np.arange(1, ndim + 1), np.arange(1, ndim + 1), indexing='ij')

    nx = j_coords - ncent
    ny = i_coords - ncent
    n2 = nx**2 + ny**2

    # Generate random phases
    ph = radphase * np.random.standard_normal((ndim, ndim))

    # Create complex field with random phase
    c_phase = np.exp(1j * ph)

    # Apply amplitude taper if specified
    if taper != 0.0:
        amplitude = np.exp(-n2 / (2.0 * sig2))
    else:
        amplitude = np.ones((ndim, ndim), dtype=np.float64)

    # Apply dish and subreflector mask
    mask = (n2 < ndish2) & (n2 > nsub2)
    c[mask] = amplitude[mask] * c_phase[mask]

    return c


def check_converge(c: np.ndarray, cold: np.ndarray, ndim: int, ncent: int,
                   ndish: float, nsub: float, maxiter: int,
                   niter: int) -> Tuple[float, bool]:
    """
    Check convergence of Misell algorithm.
    Currently uses maximum iteration number check and phase difference RMS.

    Args:
        c: Current complex aperture array
        cold: Previous complex aperture array
        ndim: Array dimension
        ncent: Center pixel
        ndish: Dish radius in pixels
        nsub: Subreflector radius in pixels
        maxiter: Maximum number of iterations
        niter: Current iteration number

    Returns:
        Tuple of (RMS error in radians, convergence flag)
    """
    nsub2 = nsub * nsub
    ndish2 = ndish * ndish

    err = 0.0
    npoint = 0

    # Create coordinate arrays
    i_coords, j_coords = np.meshgrid(np.arange(1, ndim + 1), np.arange(1, ndim + 1), indexing='ij')

    n2 = (i_coords - ncent)**2 + (j_coords - ncent)**2

    # Calculate amplitudes
    am = np.abs(c)
    am1 = np.abs(cold)

    # Calculate phases
    phase = np.angle(c)
    phase1 = np.angle(cold)

    # Select points within dish and outside subreflector with non-zero amplitude
    mask = (n2 >= nsub2) & (n2 <= ndish2) & (am != 0) & (am1 != 0)

    if np.sum(mask) > 0:
        # Calculate phase difference
        phase_diff = phase[mask] - phase1[mask]
        err = np.sqrt(np.mean(phase_diff**2))
        npoint = np.sum(mask)

    # Check if maximum iterations reached
    converged = (niter >= maxiter)

    return err, converged
