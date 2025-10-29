# HOLIS_ABER2 Python Translation

## Overview

This directory contains the Python translation of the FORTRAN holography aperture field calculation program `holis_aber2.f`. The translation preserves the original algorithms while leveraging modern Python libraries for improved readability and maintainability.

## Translation Summary

### Original FORTRAN Files (37 needed files identified)
The translation is based on `holis_aber2.f` and its dependencies. We identified and removed 26 obsolete FORTRAN files, keeping only the 37 files actually used by the main program.

### Translated Python Modules

1. **holis_aber2.py** - Main program
   - User interface and algorithm selection
   - Phase coherent holography (Algorithm 1) - **COMPLETE**
   - Misell phase recovery (Algorithm 2) - **FRAMEWORK PROVIDED**
   - Common post-processing routines

2. **holis_io.py** - Input/Output operations
   - `read_prm1()` - Read parameters for phase coherent algorithm
   - `read_amph()` - Read amplitude and phase files
   - `read_amm()` - Read multiple amplitude maps
   - `read_diff()` - Read diffraction pattern
   - `read_tk()` - Read unwrapped phase
   - `write_amph()` - Write amplitude and phase
   - `write_resiph()` - Write residual phase
   - `write_mask()` - Write mask pattern
   - **TODO**: `read_prm2()` - Needs implementation for Misell algorithm

3. **holis_fft.py** - FFT and field operations
   - `do_fft()` - 2D FFT (uses numpy.fft)
   - `do_fftshift()` - FFT shift (uses numpy.fft.fftshift)
   - `defocus()` - Apply defocus using Ruze formula
   - `undefocus()` - Remove defocus
   - `nearfield()` - Near-field corrections
   - `separate_ap()` - Separate amplitude and phase
   - `change_amplitude()` - Replace amplitude preserving phase
   - `initialize()`, `cumulate()`, `divide()` - Array operations

4. **holis_mask.py** - Masking operations
   - `apply_mask2()` - Simple circular mask
   - `apply_mask4()` - Detailed mask with quadrupod
   - `apply_mask5()` - Apply mask from file
   - `onedim()` - Convert 2D to 1D with indexing
   - `rms1d()` - Calculate RMS of 1D array

5. **holis_fit.py** - Phase fitting operations
   - `func2_aber2()` - Basis functions for fitting (8 terms):
     1. Constant offset
     2. Tilt in x
     3. Tilt in y
     4. Defocus (Ruze formula)
     5. Astigmatism at 45°
     6. Astigmatism along x
     7. Coma in x
     8. Coma in y
   - `svdfit1_aber2()` - SVD least squares fit (uses numpy.linalg.svd)
   - `phasefit_aber2()` - Main phase fitting routine

6. **holis_misell.py** - Misell-specific operations
   - `create_guess()` - Generate initial guess with random phase
   - `check_converge()` - Check convergence criteria

### Key Translation Decisions

1. **FFT Library**: Used numpy.fft instead of translating custom FFT
   - Maintains sign convention: aperture→far field = inverse FFT

2. **SVD**: Used numpy.linalg.svd instead of translating svdcmp1/svbksb1

3. **Random Numbers**: Use numpy.random instead of ran1/gasdev

4. **Array Indexing**:
   - FORTRAN uses 1-based indexing
   - Python uses 0-based indexing
   - Careful conversion in coordinate calculations

5. **File I/O**: Simplified using numpy.loadtxt/savetxt

6. **Logging**: Added logging module for better output management

## Usage

### Phase Coherent Holography (Algorithm 1)

1. Prepare parameter file `withphase_aber.prm`
2. Run the program:
   ```bash
   python holis_aber2.py
   ```
3. Select option 1 when prompted
4. The program will:
   - Read far-field amplitude and phase
   - Perform FFT to get aperture field
   - Apply near-field and defocus corrections (if specified)
   - Apply masking
   - Perform phase fitting
   - Write output files

### Misell Phase Recovery (Algorithm 2)

**STATUS**: Framework provided, requires completing `read_prm2()` function

1. Implement `read_prm2()` in `holis_io.py` following the pattern of `read_prm1()`
2. Uncomment the Misell implementation in `run_misell()` function
3. Prepare parameter file `misell.prm`
4. Run and select option 2

## Files and Dependencies

### Required Python Packages
```
numpy
scipy (optional, currently using numpy for all operations)
```

### Input Files Required (Phase Coherent)
- Parameter file: `withphase_aber.prm`
- Far-field amplitude: specified in .prm file
- Far-field phase: specified in .prm file
- Diffraction pattern: specified in .prm file

### Output Files
- Aperture amplitude and phase (masked and unmasked)
- Residual phase after fitting
- Mask pattern
- Fitting coefficients (logged)
- `phasefit.Ep` - Phase corrections applied
- `holis.log` - Full execution log

## Testing and Validation

### Recommended Testing Approach

1. **Compare with FORTRAN**:
   - Run identical input through both versions
   - Compare output files numerically
   - Expected differences:
     - Minor numerical differences due to different FFT implementations
     - Possible differences in random number generation (Misell)

2. **Test Cases**:
   - Start with Phase Coherent algorithm (fully implemented)
   - Use small test datasets (e.g., 64×64 or 128×128)
   - Verify RMS values before/after fitting
   - Check mask pattern correctness
   - Validate fitting coefficients

3. **Known Considerations**:
   - Coordinate system conversions (1-based vs 0-based indexing)
   - FFT normalization conventions
   - Phase wrapping/unwrapping

## Remaining Work

### High Priority
1. **Implement read_prm2()** in holis_io.py for Misell algorithm
2. **Test Phase Coherent algorithm** with real data
3. **Validate numerical accuracy** against FORTRAN version

### Medium Priority
1. Complete and test Misell algorithm after read_prm2() is implemented
2. Add error handling for file I/O
3. Add input validation for parameters

### Low Priority
1. Implement Algorithm 3 (Global Fitting) if needed
2. Optimize performance (vectorization, etc.)
3. Add unit tests for individual functions

## Code Structure

```
holisCodes/
├── holis_aber2.py      # Main program
├── holis_io.py         # I/O operations
├── holis_fft.py        # FFT and field ops
├── holis_mask.py       # Masking operations
├── holis_fit.py        # Phase fitting
├── holis_misell.py     # Misell operations
├── *.f                 # Original FORTRAN files (37 files kept)
└── README_TRANSLATION.md  # This file
```

## Algorithm Details

### Phase Coherent Holography (Algorithm 1)

```
1. Read far-field amplitude and phase
2. Combine into complex array
3. FFT shift
4. Forward FFT (far field → aperture)
5. FFT shift
6. Near-field correction (if enabled)
7. Defocus correction (if needed)
8. Apply mask (dish, subreflector, quadrupod)
9. Fit and remove large-scale phase errors:
   - Constant offset
   - Tilts (x, y)
   - Defocus
   - Astigmatism (2 terms)
   - Coma (x, y)
10. Write outputs
```

### Misell Phase Recovery (Algorithm 2 - Framework)

```
For each random seed (nmap times):
    1. Create initial guess (random phase, tapered amplitude)
    2. Iterate until convergence:
        For each amplitude map:
            a. Defocus to measurement plane
            b. Near-field correction
            c. FFT to far field
            d. Replace amplitude with measured
            e. FFT back to aperture
            f. Undefocus
            g. Near-field correction back
        h. Apply mask
        i. Check convergence
    3. Accumulate result
4. Average over all maps
5. Post-process (same as Algorithm 1)
```

## Contact and History

**Original FORTRAN Author**: Xiaolei Zhang (1993-1996)
**FORTRAN Modifications**: TK (1997)
**Python Translation**: 2025

For questions about the Python translation or to report issues, please document them in your project repository.

## License

[Specify license if applicable]
