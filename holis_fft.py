"""
Nimesh Patel
October 2025
Python version of previous fortran/C code from Holis package.

HOLIS_FFT.PY - FFT and field manipulation functions for holography code

This module contains FFT operations, field transformations, and
phase corrections for the holography aperture field calculation program.

Translated from FORTRAN by Python conversion 2025
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftshift
import logging
from typing import Tuple

logger = logging.getLogger(__name__)


def do_fft(c: np.ndarray, inverse: bool = False) -> np.ndarray:
    """
    Perform 2D FFT on complex array.

    Sign convention (opposite to Numerical Recipes, agrees with MATLAB):
    - From aperture to far field: inverse=True (inverse FFT)
    - From far field to aperture: inverse=False (forward FFT)

    Args:
        c: Complex 2D array
        inverse: If True, do inverse FFT; if False, do forward FFT

    Returns:
        Transformed complex array
    """
    if inverse:
        # Inverse FFT: aperture to far field
        result = ifft2(c) * c.shape[0] * c.shape[1]
    else:
        # Forward FFT: far field to aperture
        result = fft2(c) / (c.shape[0] * c.shape[1])

    return result


def do_fftshift(c: np.ndarray) -> np.ndarray:
    """
    Shift array values after FFT to conform to DFT convention.
    This is a wrapper around numpy's fftshift for consistency.

    Args:
        c: Complex 2D array

    Returns:
        Shifted array
    """
    return fftshift(c)


def defocus1(radius1: float, alambda: float, fprim: float, fm: float,
             df: float, pi: float) -> float:
    """
    Calculate defocus phase using Ruze formula for Cassegrain telescope.

    Args:
        radius1: Radial distance from center (meters)
        alambda: Wavelength (meters)
        fprim: Primary focal length (meters)
        fm: Magnified focal length (meters)
        df: Defocus amount (meters)
        pi: Pi constant

    Returns:
        Phase in radians
    """
    # Ruze formula for defocus phase
    term1 = (radius1 / (2 * fprim))**2 / (1 + (radius1 / (2 * fprim))**2)
    term2 = (radius1 / (2 * fm))**2 / (1 + (radius1 / (2 * fm))**2)

    phi = (4 * pi / alambda) * df * (term1 + term2)

    return phi


def defocus(c: np.ndarray, ndim: int, ncent: int, fprim: float, fm: float,
            pi: float, alambda: float, df: float, dprim: float, rate: float) -> np.ndarray:
    """
    Create defocused complex aperture pattern.

    Args:
        c: Complex aperture array
        ndim: Array dimension
        ncent: Center pixel
        fprim: Primary focal length
        fm: Magnified focal length
        pi: Pi constant
        alambda: Wavelength
        df: Defocus amount
        dprim: Primary diameter
        rate: Sampling rate

    Returns:
        Defocused complex array
    """
    dx = dprim / rate / ndim

    # Create coordinate arrays
    i_coords, j_coords = np.meshgrid(np.arange(ndim), np.arange(ndim), indexing='ij')

    # Calculate radius at each point
    radius1 = np.sqrt((i_coords - ncent + 1)**2 + (j_coords - ncent + 1)**2) * dx

    # Calculate phase at each point (vectorized)
    term1 = (radius1 / (2 * fprim))**2 / (1 + (radius1 / (2 * fprim))**2)
    term2 = (radius1 / (2 * fm))**2 / (1 + (radius1 / (2 * fm))**2)
    phi = (4 * pi / alambda) * df * (term1 + term2)

    # Apply phase correction
    c_defocused = c * np.exp(1j * phi)

    return c_defocused


def undefocus(c: np.ndarray, ndim: int, ncent: int, fprim: float, fm: float,
              pi: float, alambda: float, df: float, dprim: float, rate: float) -> np.ndarray:
    """
    Undo the defocus on aperture pattern (inverse of defocus operation).

    Args:
        c: Complex aperture array
        ndim: Array dimension
        ncent: Center pixel
        fprim: Primary focal length
        fm: Magnified focal length
        pi: Pi constant
        alambda: Wavelength
        df: Defocus amount
        dprim: Primary diameter
        rate: Sampling rate

    Returns:
        Undefocused complex array
    """
    dx = dprim / rate / ndim

    # Create coordinate arrays
    i_coords, j_coords = np.meshgrid(np.arange(ndim), np.arange(ndim), indexing='ij')

    # Calculate radius at each point
    radius1 = np.sqrt((i_coords - ncent + 1)**2 + (j_coords - ncent + 1)**2) * dx

    # Calculate negative phase (opposite of defocus)
    term1 = (radius1 / (2 * fprim))**2 / (1 + (radius1 / (2 * fprim))**2)
    term2 = (radius1 / (2 * fm))**2 / (1 + (radius1 / (2 * fm))**2)
    phi = -(4 * pi / alambda) * df * (term1 + term2)

    # Apply phase correction
    c_undefocused = c * np.exp(1j * phi)

    return c_undefocused


def nearfield(c: np.ndarray, ndim: int, ncent: int, ndish: float, nsub: float,
              dprim: float, distan: float, rate: float, freq: float,
              pi: float, clight: float) -> np.ndarray:
    """
    Perform second-order correction for near field effect.

    Args:
        c: Complex aperture array
        ndim: Array dimension
        ncent: Center pixel
        ndish: Dish radius in pixels
        nsub: Subreflector radius in pixels
        dprim: Primary diameter
        distan: Distance to transmitter (can be negative for correction direction)
        rate: Sampling rate
        freq: Frequency in GHz
        pi: Pi constant
        clight: Speed of light

    Returns:
        Corrected complex array
    """
    freq1 = freq * 1e9  # Convert GHz to Hz
    alambda = clight / freq1
    resol = dprim / rate / ndim
    wk = 2 * pi / alambda

    # Create coordinate arrays
    n_coords, m_coords = np.meshgrid(np.arange(ndim), np.arange(ndim), indexing='ij')

    # Calculate x, y positions
    x = (n_coords - ncent + 1) * resol
    y = (m_coords - ncent + 1) * resol
    radius = np.sqrt(x**2 + y**2)

    # Calculate phase correction
    phi = wk * radius**2 / (2 * distan)

    # Apply correction
    c_corrected = c * np.exp(1j * phi)

    return c_corrected


def separate_ap(c: np.ndarray, mask: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Separate amplitude and phase from complex aperture array.

    Args:
        c: Complex aperture array
        mask: Optional mask array (1=valid, 0=invalid)

    Returns:
        Tuple of (amplitude array, phase array)
    """
    # Calculate amplitude
    am = np.abs(c)

    # Calculate phase
    ph = np.angle(c)

    # Handle special values (masked pixels marked as -9999-9999j)
    if mask is not None:
        masked = (c == complex(-9999.0, -9999.0))
        ph[masked] = -9999.0

    # Handle zero amplitude
    ph[am == 0] = 0.0

    return am, ph


def change_amplitude(c: np.ndarray, am_new: np.ndarray) -> np.ndarray:
    """
    Change the amplitude of complex array while preserving phase.

    Args:
        c: Complex array
        am_new: New amplitude array

    Returns:
        Complex array with new amplitude
    """
    # Get current phase
    phase = np.angle(c)

    # Create new complex array with new amplitude and old phase
    c_new = am_new * np.exp(1j * phase)

    return c_new


def c_pass(c_src: np.ndarray) -> np.ndarray:
    """
    Copy complex array (pass values from source to destination).

    Args:
        c_src: Source complex array

    Returns:
        Copy of source array
    """
    return c_src.copy()


def initialize(ndim: int) -> np.ndarray:
    """
    Initialize cumulative array to zeros.

    Args:
        ndim: Array dimension

    Returns:
        Zero complex array
    """
    return np.zeros((ndim, ndim), dtype=np.complex128)


def cumulate(c: np.ndarray, ccum: np.ndarray) -> np.ndarray:
    """
    Cumulate (add) current array to cumulative array.

    Args:
        c: Current complex array
        ccum: Cumulative complex array

    Returns:
        Updated cumulative array
    """
    return ccum + c


def divide(ccum: np.ndarray, nmap: int) -> np.ndarray:
    """
    Divide cumulative array by number of maps to get average.

    Args:
        ccum: Cumulative complex array
        nmap: Number of maps

    Returns:
        Averaged complex array
    """
    return ccum / float(nmap)
