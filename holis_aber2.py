#!/usr/bin/env python3
"""
HOLIS_ABER2.PY - Holography Aperture Field Calculation

Calculates the aperture field of an antenna from the measured radiation patterns,
using different types of "holography" techniques.

Three types of algorithms are implemented:
1) Standard phase-coherent holography, which requires 1 far-field amplitude map
   and 1 far-field phase map.
2) Misell's phase-retrieval algorithm, which requires 2 far-field amplitude maps
   obtained at different focus settings.
3) Global chi-squared fitting approach (not yet implemented in this version).

Radiation patterns obtained at the near field can be accommodated.
Aperture masking and diffraction patterns are incorporated into the initial guess
for the phase-retrieval and fitting approaches.

Original FORTRAN version: Xiaolei Zhang 6/25/93
Last FORTRAN revision: 11/20/96 X.Z.
Modified for nearfield and defocus: 01/Mar/97 TK
Python translation: 2025

SIGN CONVENTION:
- Beam data is read as increasing azimuth values at a series of (increasing) elevations
- Data files are of dimension (NDIM, NDIM), where NDIM is a power of 2
- Center of map is always at (NDIM/2+1, NDIM/2+1)
- Positive defocus denotes movement of secondary away from main reflector
- From aperture to far field is inverse FFT
- From far field to aperture is forward FFT
"""

import numpy as np
import sys
import logging
import argparse
from typing import Tuple, Optional

# Import helper modules
import holis_io
import holis_fft
import holis_mask
import holis_fit
import holis_misell

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s',
    handlers=[
        logging.FileHandler('holis.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Constants from holis.inc
NMAX = 256
NMAXSQ = NMAX * NMAX
NFITMAX = 8  # Updated to 8 to match aber2 version (includes astigmatism and coma terms)

# Physical constants
PI = np.pi
CLIGHT = 2.99792e8  # Speed of light in m/s


class HolisData:
    """Container for holography data arrays and parameters."""

    def __init__(self, ndim: int):
        """
        Initialize data arrays.

        Args:
            ndim: Dimension of data arrays (must be power of 2)
        """
        self.ndim = ndim

        # Complex arrays
        self.c = np.zeros((ndim, ndim), dtype=np.complex128)  # FFT data array
        self.ccum = np.zeros((ndim, ndim), dtype=np.complex128)  # Cumulative array
        self.cmask = np.zeros((ndim, ndim), dtype=np.complex128)  # Masking array
        self.cold = np.zeros((ndim, ndim), dtype=np.complex128)  # Old c data array

        # Real arrays
        self.am = np.zeros((ndim, ndim), dtype=np.float64)  # Amplitude
        self.ph = np.zeros((ndim, ndim), dtype=np.float64)  # Phase
        self.am_um = np.zeros((ndim, ndim), dtype=np.float64)  # Unmasked amplitude
        self.ph_um = np.zeros((ndim, ndim), dtype=np.float64)  # Unmasked phase

        # Misell amplitude arrays
        self.am1 = np.zeros((ndim, ndim), dtype=np.float64)
        self.am2 = np.zeros((ndim, ndim), dtype=np.float64)
        self.am3 = np.zeros((ndim, ndim), dtype=np.float64)
        self.am4 = np.zeros((ndim, ndim), dtype=np.float64)

        # Other arrays
        self.diffp = np.zeros((ndim, ndim), dtype=np.float64)  # Unit diffraction phase map
        self.mask = np.zeros((ndim, ndim), dtype=np.int32)  # Mask pattern

        # 1D arrays for fitting (size determined later)
        self.ph1 = None
        self.sig = None
        self.ii = None  # i indices for valid masked points
        self.ij = None  # j indices for valid masked points

        # SVD arrays (size determined later)
        self.u = None
        self.v = None
        self.w = None
        self.coef = None
        self.dofit = None


class HolisParameters:
    """Container for telescope and algorithm parameters."""

    def __init__(self):
        # File names
        self.inamp = ""
        self.inphase = ""
        self.outamp = ""
        self.outphase = ""
        self.resiphase = ""
        self.resiph_um = ""
        self.fam_um = ""
        self.fph_um = ""
        self.difffile = ""
        self.maskfile = ""

        # Misell file names
        self.inamp1 = ""
        self.inamp2 = ""
        self.inamp3 = ""
        self.inamp4 = ""

        # Telescope parameters
        self.ndim = 0  # Array dimension
        self.rate = 0.0  # Sampling rate
        self.distan = 0.0  # Distance for near field
        self.freq = 0.0  # Frequency in GHz
        self.samp_itvl = 0.0  # Sampling interval
        self.dprim = 0.0  # Primary dish diameter
        self.dsec = 0.0  # Secondary diameter
        self.fprim = 0.0  # Primary focal length
        self.fmag = 0.0  # Focal magnification

        # Algorithm parameters
        self.nfit = 0  # Number of fitting parameters
        self.ncent = 0  # Center pixel
        self.ndish = 0  # Dish radius in pixels
        self.nsub = 0  # Subreflector radius in pixels
        self.ndimsq = 0  # ndim squared
        self.near_f = 0  # Near field flag
        self.defoc = 0.0  # Defocus value

        # Masking parameters
        self.out_mask = 0.0  # Outer mask radius
        self.ain_mask = 0.0  # Inner mask radius
        self.quad_hw = 0.0  # Quadrupod half-width

        # Misell parameters
        self.namp = 0  # Number of amplitude maps
        self.defo1 = 0.0  # Defocus values
        self.defo2 = 0.0
        self.defo3 = 0.0
        self.defo4 = 0.0
        self.maxiter = 0  # Maximum iterations
        self.taper = 0.0  # Edge taper
        self.radphase = 0.0  # Random phase
        self.nmap = 0  # Number of maps with different seeds
        self.idum = 0  # Random seed

        # Derived parameters
        self.freq1 = 0.0  # Frequency in Hz
        self.alambda = 0.0  # Wavelength
        self.pixbylen = 0.0  # Pixels per length


def parse_arguments():
    """
    Parse command line arguments.

    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description='Holography Aperture Field Calculation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Algorithms:
  1 - Phase coherent holography (default)
  2 - Misell phase recovery
  3 - Global fitting (not implemented)

Examples:
  %(prog)s                     # Run with defaults (algorithm 1, no unwrap)
  %(prog)s -a 2                # Use Misell algorithm
  %(prog)s --unwrap            # Use unwrapped phase from tk.dat
  %(prog)s -a 1 -u             # Explicit algorithm 1 with unwrap
        """
    )

    parser.add_argument('-a', '--algorithm', type=int, default=1,
                        choices=[1, 2, 3],
                        help='Holography algorithm (default: 1)')
    parser.add_argument('-u', '--unwrap', action='store_true',
                        help='Read unwrapped phase from tk.dat')

    return parser.parse_args()


def post_process(c: np.ndarray, params: dict, ndim: int, ncent: int,
                 ndish: float, nsub: float, pixbylen: float, use_unwrap: bool = False):
    """
    Common post-processing for all algorithms.
    Includes masking, phase fitting, and output writing.

    Args:
        c: Complex aperture array
        params: Parameter dictionary
        ndim: Array dimension
        ncent: Center pixel
        ndish: Dish radius in pixels
        nsub: Subreflector radius in pixels
        pixbylen: Pixels per unit length
        use_unwrap: If True, read unwrapped phase from tk.dat
    """
    # Separate amplitude and phase
    am, ph = holis_fft.separate_ap(c)

    # Save unmasked amplitude and phase
    am_um = am.copy()
    ph_um = ph.copy()

    # Write unmasked files
    logger.info(f"Writing unmasked amplitude to {params['fam_um']}...")
    holis_io.write_resiph(am_um, params['fam_um'])

    logger.info(f"Writing unmasked phase to {params['fph_um']}...")
    holis_io.write_resiph(ph_um, params['fph_um'])

    # Apply mask
    logger.info("Applying mask to aperture field...")
    c, mask = holis_mask.apply_mask4(c, ndim, ncent, pixbylen,
                                     params['out_mask'], params['ain_mask'],
                                     params['quad_hw'])

    # Write mask
    logger.info(f"Writing mask to {params['maskfile']}...")
    holis_io.write_mask(mask, params['maskfile'])

    # Separate amplitude and phase after masking
    am, ph = holis_fft.separate_ap(c, mask)

    # Write aperture amplitude and phase
    logger.info(f"Writing aperture amplitude to {params['outamp']}...")
    logger.info(f"Writing aperture phase to {params['outphase']}...")
    holis_io.write_amph(am, ph, params['outamp'], params['outphase'])

    # Read unwrapped phase file if requested
    if use_unwrap:
        logger.info("Reading unwrapped phase from tk.dat...")
        ph = holis_io.read_tk('tk.dat', ndim)

    # Convert to 1D and calculate RMS before fit
    ph1, ii, ij, nmasked = holis_mask.onedim(ph, mask)
    rms1 = holis_mask.rms1d(ph1)
    logger.info("")
    logger.info(f"RMS before phasefit = {rms1:.8e} radian")

    # Read diffraction file
    logger.info(f"Reading diffraction pattern from {params['difffile']}...")
    diffp = holis_io.read_diff(params['difffile'], ndim)

    # Perform phase fit
    logger.info("Performing phase fit...")
    logger.info("phasefit_aber2 call")
    ph, ph_um, chisq = holis_fit.phasefit_aber2(
        ph, ph_um, ph1, ii, ij, ndim, ncent,
        params['fprim'], params['fmag'], diffp, params['dprim'],
        params['rate'], params['freq'], PI, CLIGHT,
        params['dofit'], mask
    )

    # Calculate RMS after fit
    ph1, ii, ij, nmasked = holis_mask.onedim(ph, mask)
    rms1 = holis_mask.rms1d(ph1)

    # Convert RMS to microns: surface_error = phase * lambda / (4*pi)
    freq1 = params['freq'] * 1e9  # Convert GHz to Hz
    alambda = CLIGHT / freq1  # Wavelength in meters
    rms1_microns = rms1 * alambda / (4 * PI) * 1e6  # Convert to microns

    logger.info(f"RMS after phasefit = {rms1:.8e} radian ({rms1_microns:.8e} microns)")

    # Apply mask to residual phase
    ph = holis_mask.apply_mask5(ph, params['maskfile'])

    # Write residual phase
    logger.info(f"Writing residual phase to {params['resiphase']}...")
    holis_io.write_resiph(ph, params['resiphase'])

    logger.info(f"Writing unmasked residual phase to {params['resiph_um']}...")
    holis_io.write_resiph(ph_um, params['resiph_um'])

    logger.info("")
    logger.info("Processing complete!")



def main():
    """Main program entry point."""

    # Parse command line arguments
    args = parse_arguments()

    logger.info("="*70)
    logger.info("HOLIS_ABER2 - Holography Aperture Field Calculation")
    logger.info("="*70)
    logger.info("")

    if args.algorithm == 1:
        logger.info("Use with-phase holography algorithm")
        logger.info("(You must first fill out 'withphase_aber.prm')")
        logger.info("")
        run_phase_coherent(args.unwrap)

    elif args.algorithm == 2:
        logger.info("Use Misell phase retrieval algorithm")
        logger.info("(You must first fill out 'misell.prm')")
        logger.info("")
        run_misell(args.unwrap)

    elif args.algorithm == 3:
        logger.info("Use global fitting holography algorithm")
        logger.info("(You must first fill out 'fitting.prm')")
        logger.info("")
        logger.error("This algorithm not implemented yet")
        sys.exit(1)


def run_phase_coherent(use_unwrap: bool = False):
    """
    Implementation of phase coherent holography (Algorithm 1).

    This section performs the inverse Fast Fourier Transform (FFT) to obtain
    the complex aperture field distribution from the input complex far-field
    distribution.

    Args:
        use_unwrap: If True, read unwrapped phase from tk.dat
    """
    logger.info("Starting phase coherent holography...")

    # Read parameters
    params = holis_io.read_prm1('withphase_aber.prm')

    # Extract key parameters
    ndim = params['ndim']
    ncent = params['ncent']
    ndish = params['ndish']
    nsub = params['nsub']
    rate = params['rate']
    distan = params['distan']
    freq = params['freq']
    dprim = params['dprim']
    dsec = params['dsec']
    fprim = params['fprim']
    fmag = params['fmag']
    near_f = params['near_f']
    defoc = params['defoc']
    pixbylen = 2 * ndish / dprim

    # Calculate derived constants
    freq1 = freq * 1e9  # Convert GHz to Hz
    alambda = CLIGHT / freq1

    # Read far-field amplitude and phase
    logger.info(f"Reading far-field data from {params['inamp']} and {params['inphase']}...")
    am, ph, c = holis_io.read_amph(params['inamp'], params['inphase'], ndim)

    # Do origin shift before FFT
    c = holis_fft.do_fftshift(c)

    # Forward Fourier transform from far field to aperture
    logger.info("Performing FFT from far field to aperture...")
    c = holis_fft.do_fft(c, inverse=False)

    # Do origin shift after FFT
    c = holis_fft.do_fftshift(c)

    # Perform second-order correction for near field effect
    if near_f == 1:
        logger.info("Applying near-field correction...")
        c = holis_fft.nearfield(c, ndim, ncent, ndish, nsub, dprim, distan,
                                rate, freq, PI, CLIGHT)

    # Take out a rough defocus if needed
    if defoc != 0.0:
        logger.info(f"Applying defocus correction: {defoc} mm...")
        fm = fprim * fmag
        c = holis_fft.undefocus(c, ndim, ncent, fprim, fm, PI, alambda,
                                defoc, dprim, rate)

    # Common post-processing
    post_process(c, params, ndim, ncent, ndish, nsub, pixbylen, use_unwrap)


def run_misell(use_unwrap: bool = False):
    """
    Implementation of Misell's phase recovery holography (Algorithm 2).

    This algorithm uses 2-4 far-field amplitude maps obtained at different
    focus settings to recover the aperture phase through iterative refinement.

    Args:
        use_unwrap: If True, read unwrapped phase from tk.dat
    """
    logger.info("Starting Misell phase retrieval...")
    logger.info("NOTE: Misell algorithm translation requires parameter file reader.")
    logger.info("Please create a read_prm2 function in holis_io.py based on read_prm1")
    logger.info("as a template, then uncomment the code below.")

    # NOTE: This implementation requires read_prm2 to be added to holis_io.py
    # Uncomment below after implementing read_prm2

    raise NotImplementedError(
        "Misell algorithm requires read_prm2() function in holis_io.py.\n"
        "Please implement it based on read_prm2.f, following the pattern of read_prm1()."
    )

    # # Read parameters
    # params = holis_io.read_prm2('misell.prm')
    #
    # # Extract parameters
    # ndim = params['ndim']
    # ncent = params['ncent']
    # ...
    #
    # # Read amplitude maps
    # logger.info("Reading amplitude maps...")
    # amp_files = [params['inamp1'], params['inamp2']]
    # if params['namp'] >= 3:
    #     amp_files.append(params['inamp3'])
    # if params['namp'] == 4:
    #     amp_files.append(params['inamp4'])
    #
    # amps = holis_io.read_amm(amp_files, ndim, params['namp'])
    #
    # # Initialize for multiple maps with different seeds
    # ccum = holis_fft.initialize(ndim)
    # niter = 0
    # imap = 0
    #
    # # Iteration on maps with different initial seeds
    # while imap < params['nmap']:
    #     imap += 1
    #     logger.info(f"\nGenerating map {imap} of {params['nmap']}...")
    #
    #     # Create initial guess
    #     c = holis_misell.create_guess(ndim, ncent, ndish, nsub,
    #                                   params['taper'], params['radphase'],
    #                                   params['idum'])
    #     cold = c.copy()
    #
    #     # Misell iteration loop
    #     niter = 0
    #     converged = False
    #
    #     while not converged:
    #         # Iterate through all defocus/amplitude pairs
    #         for i, (defo, am) in enumerate(zip(defocus_values, amps)):
    #             # Defocus aperture
    #             c = holis_fft.defocus(c, ndim, ncent, fprim, fm, PI,
    #                                  alambda, defo, dprim, rate)
    #
    #             # Near field correction (going to far field)
    #             c = holis_fft.nearfield(c, ndim, ncent, ndish, nsub,
    #                                    dprim, -distan, rate, freq, PI, CLIGHT)
    #
    #             # Propagate to far field
    #             c = holis_fft.do_fftshift(c)
    #             c = holis_fft.do_fft(c, inverse=True)
    #             c = holis_fft.do_fftshift(c)
    #
    #             # Replace amplitude with measured
    #             c = holis_fft.change_amplitude(c, am)
    #
    #             # Propagate back to aperture
    #             c = holis_fft.do_fftshift(c)
    #             c = holis_fft.do_fft(c, inverse=False)
    #             c = holis_fft.do_fftshift(c)
    #
    #             # Undefocus
    #             c = holis_fft.undefocus(c, ndim, ncent, fprim, fm, PI,
    #                                    alambda, defo, dprim, rate)
    #
    #             # Near field correction (back to aperture)
    #             c = holis_fft.nearfield(c, ndim, ncent, ndish, nsub,
    #                                    dprim, distan, rate, freq, PI, CLIGHT)
    #
    #         # Apply mask
    #         c = holis_mask.apply_mask2(c, ndim, ncent, ndish, nsub)
    #
    #         # Check convergence
    #         niter += 1
    #         err, converged = holis_misell.check_converge(
    #             c, cold, ndim, ncent, ndish, nsub, params['maxiter'], niter
    #         )
    #
    #         logger.info(f"Iteration {niter}: RMS convergence error = {err:.8e} radian")
    #
    #         # Update cold
    #         cold = c.copy()
    #
    #     logger.info(f"Map {imap} converged after {niter} iterations")
    #
    #     # Accumulate result
    #     ccum = holis_fft.cumulate(c, ccum)
    #
    #     # Update seed for next map
    #     params['idum'] += 95
    #     niter = 0
    #
    # # Average over all maps
    # c = holis_fft.divide(ccum, params['nmap'])
    #
    # # Common post-processing
    # post_process(c, params, ndim, ncent, ndish, nsub, pixbylen, use_unwrap)


if __name__ == "__main__":
    main()
