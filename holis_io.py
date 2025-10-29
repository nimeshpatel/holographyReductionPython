"""
HOLIS_IO.PY - Input/Output functions for holography code

This module contains all file reading and writing functions for the
holography aperture field calculation program.

Translated from FORTRAN by Python conversion 2025
"""

import numpy as np
import logging
from typing import Tuple, List

logger = logging.getLogger(__name__)


def read_prm1(filename: str = 'withphase_aber.prm') -> dict:
    """
    Read parameters for phase coherent holography (Algorithm 1).

    Args:
        filename: Parameter file name

    Returns:
        Dictionary containing all parameters
    """
    logger.info("")
    logger.info("-- Reading parameters from '%s' --", filename)
    logger.info("")

    params = {}
    dofit = []

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Skip comment lines (starting with '!')
    def get_next_line(lines_iter):
        for line in lines_iter:
            if not line.strip().startswith('!') and line.strip():
                return line.strip()
        return None

    lines_iter = iter(lines)

    # Read file names and parameters in order
    params['inamp'] = get_next_line(lines_iter).split()[-1]
    logger.info(f"Far field amplitude file name = {params['inamp']}")

    params['inphase'] = get_next_line(lines_iter).split()[-1]
    logger.info(f"Far field phase file name = {params['inphase']}")

    params['outamp'] = get_next_line(lines_iter).split()[-1]
    logger.info(f"Aperture field amplitude file name = {params['outamp']}")

    params['outphase'] = get_next_line(lines_iter).split()[-1]
    logger.info(f"Aperture field phase file name = {params['outphase']}")

    params['resiphase'] = get_next_line(lines_iter).split()[-1]
    logger.info(f"Residual aperture phase file name = {params['resiphase']}")

    params['resiph_um'] = get_next_line(lines_iter).split()[-1]
    logger.info(f"Unmasked residual aperture phase file name = {params['resiph_um']}")

    params['fam_um'] = get_next_line(lines_iter).split()[-1]
    logger.info(f"Unmasked aperture amplitude file name = {params['fam_um']}")

    params['fph_um'] = get_next_line(lines_iter).split()[-1]
    logger.info(f"Unmasked aperture phase file name = {params['fph_um']}")

    # Numeric parameters
    params['ndim'] = int(get_next_line(lines_iter).split()[-1])
    logger.info(f"Size N of the N by N array = {params['ndim']}")

    params['rate'] = float(get_next_line(lines_iter).split()[-1])
    logger.info(f"Nyquist sampling rate = {params['rate']}")

    params['distan'] = float(get_next_line(lines_iter).split()[-1])
    logger.info(f"Distance of the transmitter (meters) = {params['distan']}")

    params['freq'] = float(get_next_line(lines_iter).split()[-1])
    logger.info(f"Observing Frequency (GHz) = {params['freq']}")

    params['samp_itvl'] = float(get_next_line(lines_iter).split()[-1])
    logger.info(f"Sampling interval in the far field (arcsec) = {params['samp_itvl']}")

    params['dprim'] = float(get_next_line(lines_iter).split()[-1])
    logger.info(f"Diameter of the primary mirror (m) = {params['dprim']}")

    params['dsec'] = float(get_next_line(lines_iter).split()[-1])
    logger.info(f"Diameter of the secondary mirror (m) = {params['dsec']}")

    params['fprim'] = float(get_next_line(lines_iter).split()[-1])
    logger.info(f"Focal length of the primary (m) = {params['fprim']}")

    params['fmag'] = float(get_next_line(lines_iter).split()[-1])
    logger.info(f"Magnification of the Cassegrain system = {params['fmag']}")

    params['difffile'] = get_next_line(lines_iter).split()[-1]
    logger.info(f"Name of phase diffraction file = {params['difffile']}")

    params['nfit'] = int(get_next_line(lines_iter).split()[-1])
    logger.info(f"Number of large scale fitting parameters = {params['nfit']}")

    params['near_f'] = int(get_next_line(lines_iter).split()[-1])
    logger.info(f"Near field correction = {params['near_f']}")

    params['defoc'] = float(get_next_line(lines_iter).split()[-1])
    logger.info(f"Defocus correction = {params['defoc']}")

    # Derived parameters
    params['ncent'] = params['ndim'] // 2 + 1  # Center of dish in mesh
    params['ndish'] = params['ndim'] // 2 * params['rate']  # Radius of dish in pixels
    ratio = params['dsec'] / params['dprim']  # Ratio of secondary to primary dia
    params['nsub'] = params['ndish'] * ratio  # Radius of subreflector in pixels
    params['ndimsq'] = params['ndim'] * params['ndim']  # Dimension of 1D array

    params['out_mask'] = float(get_next_line(lines_iter).split()[-1])
    logger.info(f"Outer dia for masking = {params['out_mask']}")

    params['ain_mask'] = float(get_next_line(lines_iter).split()[-1])
    logger.info(f"Inner dia for masking = {params['ain_mask']}")

    params['quad_hw'] = float(get_next_line(lines_iter).split()[-1])
    logger.info(f"Half-width for quadrupod masking = {params['quad_hw']}")

    params['maskfile'] = get_next_line(lines_iter).split()[-1]
    logger.info(f"Maskfile name = {params['maskfile']}")

    # Read fitting function flags
    for i in range(params['nfit']):
        dofit_val = int(get_next_line(lines_iter).split()[-1])
        dofit.append(dofit_val)
        logger.info(f"Fitting function {i+1} ON/OFF = {dofit_val}")

    params['dofit'] = np.array(dofit, dtype=np.int32)

    return params


def read_amph(amp_file: str, phase_file: str, ndim: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read amplitude and phase files and create complex array.

    Args:
        amp_file: Amplitude file name
        phase_file: Phase file name
        ndim: Array dimension

    Returns:
        Tuple of (amplitude array, phase array, complex array)
    """
    # Read amplitude
    am = np.loadtxt(amp_file, dtype=np.float64).reshape((ndim, ndim))

    # Read phase
    ph = np.loadtxt(phase_file, dtype=np.float64).reshape((ndim, ndim))

    # Create complex array
    c = am * np.exp(1j * ph)

    return am, ph, c


def read_amm(files: List[str], ndim: int, namp: int) -> List[np.ndarray]:
    """
    Read multiple amplitude maps for Misell algorithm.

    Args:
        files: List of amplitude file names
        ndim: Array dimension
        namp: Number of amplitude maps (2, 3, or 4)

    Returns:
        List of amplitude arrays
    """
    amps = []
    for i in range(namp):
        am = np.loadtxt(files[i], dtype=np.float64).reshape((ndim, ndim))
        amps.append(am)

    return amps


def read_diff(filename: str, ndim: int) -> np.ndarray:
    """
    Read unit diffraction phase map.

    Args:
        filename: Diffraction file name
        ndim: Array dimension

    Returns:
        Diffraction phase array
    """
    diffp = np.loadtxt(filename, dtype=np.float64).reshape((ndim, ndim))
    return diffp


def read_tk(filename: str, ndim: int) -> np.ndarray:
    """
    Read unwrapped phase file (tk.dat).

    Args:
        filename: Unwrapped phase file name
        ndim: Array dimension

    Returns:
        Phase array
    """
    ph = np.loadtxt(filename, dtype=np.float64).reshape((ndim, ndim))
    return ph


def write_amph(am: np.ndarray, ph: np.ndarray, amp_file: str, phase_file: str):
    """
    Write amplitude and phase arrays to files.
    Flattens the 2D arrays to write one value per line (Fortran-style).

    Args:
        am: Amplitude array
        ph: Phase array
        amp_file: Amplitude output file name
        phase_file: Phase output file name
    """
    np.savetxt(amp_file, am.flatten(), fmt='%15.8e')
    np.savetxt(phase_file, ph.flatten(), fmt='%15.8e')


def write_resiph(ph: np.ndarray, filename: str):
    """
    Write residual phase to file.
    Flattens the 2D array to write one value per line (Fortran-style).

    Args:
        ph: Phase array
        filename: Output file name
    """
    np.savetxt(filename, ph.flatten(), fmt='%15.8e')


def write_mask(mask: np.ndarray, filename: str):
    """
    Write mask pattern to file.
    Flattens the 2D array to write one value per line (Fortran-style).

    Args:
        mask: Mask array (integer array with 0s and 1s)
        filename: Output file name
    """
    np.savetxt(filename, mask.flatten(), fmt='%d')
