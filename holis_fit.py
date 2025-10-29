"""
HOLIS_FIT.PY - Phase fitting operations for holography code

This module contains functions for least-squares fitting of large-scale
phase errors including constant offset, tilt, defocus, astigmatism, and coma.

Translated from FORTRAN by Python conversion 2025
"""

import numpy as np
import logging
from typing import Tuple, Callable

logger = logging.getLogger(__name__)


def func2_aber2(i: int, ii: np.ndarray, ij: np.ndarray, ndim: int, ncent: int,
                fprim: float, fmag: float, diffp: np.ndarray, dprim: float,
                rate: float, freq: float, pi: float, clight: float,
                dofit: np.ndarray) -> np.ndarray:
    """
    Calculate basis functions for phase fitting (aberration version 2).

    Fitting functions:
    1. Constant offset
    2. Tilt in x direction
    3. Tilt in y direction
    4. Defocus (using Ruze formula)
    5. Astigmatism at 45 degrees
    6. Astigmatism along x
    7. Coma in x
    8. Coma in y

    Args:
        i: Index into 1D arrays (0-based in Python)
        ii: Array of i indices (1-based from Fortran)
        ij: Array of j indices (1-based from Fortran)
        ndim: Array dimension
        ncent: Center pixel
        fprim: Primary focal length
        fmag: Focal magnification
        diffp: Diffraction pattern (not used in aber2 version)
        dprim: Primary diameter
        rate: Sampling rate
        freq: Frequency in GHz
        pi: Pi constant
        clight: Speed of light
        dofit: Array indicating which functions to fit (1=on, 0=off)

    Returns:
        Array of basis function values
    """
    freq1 = freq * 1e9  # Convert GHz to Hz
    alambda = clight / freq1  # Wavelength in m
    fm = fprim * fmag

    # Get 2D indices (convert from 1-based to 0-based)
    # Note: In FORTRAN, ii contains rows, ij contains columns
    # The FORTRAN code confusingly names them n=ii(i) and m=ij(i)
    row_idx = ii[i] - 1
    col_idx = ij[i] - 1

    # Calculate coordinates from center
    # x corresponds to column (horizontal), y to row (vertical)
    x = float(col_idx - (ncent - 1))
    y = float(row_idx - (ncent - 1))
    radius = np.sqrt(x**2 + y**2)
    radius1 = radius * dprim / rate / ndim

    # Initialize basis functions
    afunc = np.zeros(8, dtype=np.float64)

    # Function 1: Constant offset
    afunc[0] = dofit[0] * 1.0

    # Function 2: Tilt in x direction
    afunc[1] = dofit[1] * x

    # Function 3: Tilt in y direction
    afunc[2] = dofit[2] * y

    # Function 4: Defocus (Ruze formula)
    term1 = ((radius1 / 2 / fprim)**2) / (1 + (radius1 / 2 / fprim)**2)
    term2 = ((radius1 / 2 / fm)**2) / (1 + (radius1 / 2 / fm)**2)
    afunc[3] = dofit[3] * (4 * pi / alambda * (term1 + term2))

    # Function 5: Astigmatism at 45 degrees (replaces diffraction)
    if radius > 0:
        afunc[4] = dofit[4] * (radius1**2) * (1.414 * (x + y) / radius)**2
    else:
        afunc[4] = 0.0

    # Function 6: Astigmatism along x
    if radius > 0:
        afunc[5] = dofit[5] * (radius1**2) * (x / radius)**2
    else:
        afunc[5] = 0.0

    # Function 7: Coma in x
    if radius > 0:
        afunc[6] = dofit[6] * (radius1**3) * (x / radius)
    else:
        afunc[6] = 0.0

    # Function 8: Coma in y
    if radius > 0:
        afunc[7] = dofit[7] * (radius1**3) * (y / radius)
    else:
        afunc[7] = 0.0

    return afunc


def svdfit1_aber2(ph1: np.ndarray, sig: np.ndarray, ii: np.ndarray,
                  ij: np.ndarray, ndim: int, ncent: int, fprim: float,
                  fmag: float, diffp: np.ndarray, dprim: float, rate: float,
                  freq: float, pi: float, clight: float,
                  dofit: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Perform SVD-based least squares fit for phase data.

    Args:
        ph1: 1D phase data array
        sig: 1D error array (typically all 1.0)
        ii: i indices of valid points
        ij: j indices of valid points
        ndim: Array dimension
        ncent: Center pixel
        fprim: Primary focal length
        fmag: Focal magnification
        diffp: Diffraction pattern
        dprim: Primary diameter
        rate: Sampling rate
        freq: Frequency in GHz
        pi: Pi constant
        clight: Speed of light
        dofit: Array indicating which functions to fit

    Returns:
        Tuple of (coefficient array, chi-squared value)
    """
    ndata = len(ph1)
    ma = len(dofit)

    # Build design matrix
    u = np.zeros((ndata, ma), dtype=np.float64)
    b = np.zeros(ndata, dtype=np.float64)

    for i in range(ndata):
        afunc = func2_aber2(i, ii, ij, ndim, ncent, fprim, fmag, diffp,
                            dprim, rate, freq, pi, clight, dofit)
        tmp = 1.0 / sig[i]
        u[i, :] = afunc * tmp
        b[i] = ph1[i] * tmp

    # Perform SVD
    U, w, Vt = np.linalg.svd(u, full_matrices=False)

    # Threshold small singular values
    TOL = 1.e-10
    wmax = np.max(w)
    thresh = TOL * wmax
    w[w < thresh] = 0.0

    # Back-substitution
    # coef = V @ (w^-1 @ (U^T @ b))
    w_inv = np.zeros_like(w)
    w_inv[w != 0] = 1.0 / w[w != 0]

    coef = Vt.T @ (w_inv * (U.T @ b))

    # Log coefficients
    logger.info(f"coef(1) = {coef[0]:15.8e}   ! constant offset")
    logger.info(f"coef(2) = {coef[1]:15.8e}   ! tilt in x direction")
    logger.info(f"coef(3) = {coef[2]:15.8e}   ! tilt in y direction")
    logger.info(f"coef(4) = {coef[3]*1000:15.8e}   ! defocus in mm")
    logger.info(f"coef(5) = {coef[4]:15.8e}   ! astigmatism-45 term")
    logger.info(f"coef(6) = {coef[5]:15.8e}   ! astigmatism term")
    logger.info(f"coef(7) = {coef[6]:15.8e}   ! coma-x term")
    logger.info(f"coef(8) = {coef[7]:15.8e}   ! coma-y term")

    # Calculate chi-squared
    chisq = 0.0
    for i in range(ndata):
        afunc = func2_aber2(i, ii, ij, ndim, ncent, fprim, fmag, diffp,
                            dprim, rate, freq, pi, clight, dofit)
        fitted = np.dot(coef, afunc)
        chisq += ((ph1[i] - fitted) / sig[i])**2

    return coef, chisq


def phasefit_aber2(ph: np.ndarray, ph_um: np.ndarray, ph1: np.ndarray,
                   ii: np.ndarray, ij: np.ndarray, ndim: int, ncent: int,
                   fprim: float, fmag: float, diffp: np.ndarray, dprim: float,
                   rate: float, freq: float, pi: float, clight: float,
                   dofit: np.ndarray, mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Least squares fitting on aperture phase large scale error terms.

    Args:
        ph: 2D masked phase array
        ph_um: 2D unmasked phase array
        ph1: 1D phase array of valid points
        ii: i indices of valid points
        ij: j indices of valid points
        ndim: Array dimension
        ncent: Center pixel
        fprim: Primary focal length
        fmag: Focal magnification
        diffp: Diffraction pattern
        dprim: Primary diameter
        rate: Sampling rate
        freq: Frequency in GHz
        pi: Pi constant
        clight: Speed of light
        dofit: Array indicating which functions to fit
        mask: Mask array

    Returns:
        Tuple of (corrected masked phase, corrected unmasked phase, chi-squared)
    """
    nmasked = len(ph1)

    # Create uniform error array
    sig = np.ones(nmasked, dtype=np.float64)

    # Perform SVD fit
    coef, chisq = svdfit1_aber2(ph1, sig, ii, ij, ndim, ncent, fprim, fmag,
                                diffp, dprim, rate, freq, pi, clight, dofit)

    logger.info(f"chisq = {chisq:15.8e}")

    # Calculate correction and apply to phase arrays
    freq1 = freq * 1e9
    alambda = clight / freq1
    fm = fprim * fmag

    # Open file for phase correction values
    with open('phasefit.Ep', 'w') as f:
        for m in range(ndim):
            for n in range(ndim):
                x = float(n - ncent + 1)
                y = float(m - ncent + 1)
                radius = np.sqrt(x**2 + y**2)
                radius1 = radius * dprim / rate / ndim

                # Calculate correction
                correction = (
                    -coef[0]
                    - coef[1] * x
                    - coef[2] * y
                    - coef[4] * (radius1**2) * (1.414 * (x + y) / radius)**2 if radius > 0 else 0.0
                )

                # Defocus term
                term1 = ((radius1 / 2 / fprim)**2) / (1 + (radius1 / 2 / fprim)**2)
                term2 = ((radius1 / 2 / fm)**2) / (1 + (radius1 / 2 / fm)**2)
                correction -= coef[3] * 4 * pi / alambda * (term1 + term2)

                # Astigmatism and coma terms
                if radius > 0:
                    correction -= coef[5] * (radius1**2) * (x / radius)**2
                    correction -= coef[6] * (radius1**3) * (x / radius)
                    correction -= coef[7] * (radius1**3) * (y / radius)

                # Apply correction to unmasked phase
                ph_um[m, n] += correction

                # Apply correction to masked phase
                if mask[m, n] == 1 or ph[m, n] != -9999.0:
                    ph[m, n] += correction
                else:
                    ph[m, n] = -9999.0

                # Write correction to file
                f.write(f"{correction:15.8e}\n")

    return ph, ph_um, chisq
