#!/usr/bin/env python3
"""
Nimesh Patel
27 October 2025

Regrid raw on-the-fly holography scan data onto a regular grid.

This code bins irregularly sampled (el, az, amplitude, phase) data 
into a regular 2D grid, averaging samples that fall within each grid bin.
Phase averaging is done using vector averaging to handle phase wrapping correctly.

CORRECTED VERSION: This version fixes a bug in the original C code where amplitude
was incorrectly computed from the last sample instead of the vector average.

Usage: python regrid.py rawfile regrid.prm
"""

from pathlib import Path
import math
import sys


def read_regrid_params(param_file: Path):
    """
    Read regridding parameters from parameter file.
    
    Expected format (single line):
    naz nel azstrt azstep0 elstrt elstep az_thresh el_thresh
    
    Returns:
        dict with parameters
    """
    with open(param_file, 'r') as f:
        line = f.readline().strip()
        parts = line.split()
        
    if len(parts) < 8:
        raise ValueError(f"Parameter file must have 8 values, got {len(parts)}")
    
    params = {
        'naz': int(parts[0]),           # Number of azimuth bins
        'nel': int(parts[1]),           # Number of elevation bins
        'azstrt': float(parts[2]),      # Starting azimuth
        'azstep0': float(parts[3]),     # Azimuth step at reference elevation
        'elstrt': float(parts[4]),      # Starting elevation
        'elstep': float(parts[5]),      # Elevation step
        'az_thresh': float(parts[6]),   # Azimuth threshold for bin assignment
        'el_thresh': float(parts[7])    # Elevation threshold for new row
    }
    
    print(f"Regridding parameters:")
    print(f"  Grid size: {params['nel']} × {params['naz']}")
    print(f"  Azimuth: start={params['azstrt']}, step={params['azstep0']}")
    print(f"  Elevation: start={params['elstrt']}, step={params['elstep']}")
    print(f"  Thresholds: az={params['az_thresh']}, el={params['el_thresh']}")
    
    return params


def regrid_data(rawfile: Path, param_file: Path, outfile: Path):
    """
    Regrid on-the-fly scan data onto a regular grid.
    
    The algorithm:
    1. For each elevation row, calculate azimuth step (corrected for elevation)
    2. For each azimuth bin, collect all data points within threshold
    3. Average amplitude and phase using complex vector averaging
    4. Handle missing bins by interpolation or setting to zero
    
    Complex vector averaging:
    - Convert each (amp, phase) to complex: z = amp * exp(i*phase)
    - Average the complex values: z_avg = mean(z)
    - Extract: amp_avg = |z_avg|, phase_avg = arg(z_avg)
    
    This properly handles phase wrapping and gives the correct amplitude
    of the averaged complex signal.
    
    Args:
        rawfile: Input file with columns: el az amplitude phase
        param_file: Parameter file with grid specifications
        outfile: Output file for regridded data (j k amplitude phase)
    """
    params = read_regrid_params(param_file)
    
    naz = params['naz']
    nel = params['nel']
    azstrt = params['azstrt']
    azstep0 = params['azstep0']
    elstrt = params['elstrt']
    elstep = params['elstep']
    az_thresh = params['az_thresh']
    el_thresh = params['el_thresh']
    
    pi = math.pi
    
    # Initialize output grids
    amp_grd = [[0.0 for _ in range(naz)] for _ in range(nel)]
    ph_grd = [[0.0 for _ in range(naz)] for _ in range(nel)]
    
    # Open files
    print(f"Reading raw data from {rawfile}...")
    fin = open(rawfile, 'r')
    fout = open(outfile, 'w')
    
    # State variables
    el_prev = elstrt
    el = elstrt
    elgrd = elstrt
    new_el = 0
    wait_for_row = 0
    end_flag = 0
    line_num = 0
    
    # Loop over elevation rows
    j = 0
    while j <= nel - 1:
        azgrd = azstrt
        azstep = azstep0 / math.cos(el * pi / 180.0)
        
        i = 0  # Number of points in current az bin
        start_flag = 0
        if new_el == 1:
            end_flag = 0
        flag = 0  # Flag for missing data in previous bin
        k = 0  # Azimuth grid point counter
        
        # Initialize first bin
        amp_grd[j][k] = 0.0
        ph_grd[j][k] = 0.0
        sin_ph_grd = 0.0
        cos_ph_grd = 0.0
        xave = 0.0  # Running average of real part (amp*cos(phase))
        yave = 0.0  # Running average of imaginary part (amp*sin(phase))
        
        # Process data points for this elevation row
        while True:
            line = fin.readline()
            if not line:  # EOF
                break
            
            parts = line.strip().split()
            if len(parts) < 4:
                continue
                
            az = float(parts[0])
            el = float(parts[1])
            amp = float(parts[2])
            ph = float(parts[3])
            line_num += 1
            
            # Check if we've moved to a new elevation row
            if new_el != 1:
                if (el - elgrd) > el_thresh:
                    if el >= elstrt:
                        el_prev = el
                    new_el = 1
                    wait_for_row = 1
                    elgrd = elgrd + elstep
                    if j == 0 and k == 0 and end_flag != 1:
                        j -= 1
                    end_flag = 0
                    break
            
            # Check if point is within current azimuth bin
            if abs(az - azgrd) < az_thresh:
                # Inside current bin - accumulate
                new_el = 0
                wait_for_row = 0
                end_of_bin = 0
                
                # Convert phase to radians (input is in degrees)
                ph = ph * (pi / 180.0)
                
                start_flag = 1
                
                # Accumulate for simple averaging (used for diagnostic)
                amp_grd[j][k] = amp_grd[j][k] + amp
                
                # Accumulate phase using vector averaging
                sin_ph = math.sin(ph)
                cos_ph = math.cos(ph)
                sin_ph_grd = sin_ph_grd + sin_ph
                cos_ph_grd = cos_ph_grd + cos_ph
                
                # Complex vector average - properly average the complex quantity
                x = amp * cos_ph
                y = amp * sin_ph
                # Running average: new_avg = (old_avg * n + new_value) / (n+1)
                xave = (xave * i + x) / (i + 1)
                yave = (yave * i + y) / (i + 1)
                
                i += 1
            
            # End of current bin
            elif end_flag == 1:
                pass  # Already at end of row
            
            elif (azgrd - az) < 0:
                # Data point is past current bin - finalize this bin
                if wait_for_row != 1:
                    # Process the completed bin
                    if i == 0:
                        # No data in this bin
                        if k == 0:
                            amp_grd[j][k] = 0.0
                            ph_grd[j][k] = 0.0
                            ph = 0.0
                            print(f" :00 {line_num} {i} {j} {k} {el:.6f} {azgrd:.6f} {amp_grd[j][k]:.6f} {ph:.6f}")
                        else:
                            if flag == 1:
                                # Two successive bins without data
                                print(f"two successive points without data, using prev. values")
                                ph = 0.0
                                print(f" :prev {line_num} {i} {j} {k-1} {el:.6f} {azgrd:.6f} {amp_grd[j][k-2]:.6f} {ph:.6f}")
                                ph = 0.0
                                print(f" :prev {line_num} {i} {j} {k} {el:.6f} {azgrd:.6f} {amp_grd[j][k-2]:.6f} {ph:.6f}")
                                flag = 0
                            else:
                                flag = 1
                    else:
                        # Average the accumulated data using complex vector average
                        # This is the CORRECTED version - uses xave, yave (the averages)
                        # not x, y (the last sample) as in the buggy C code
                        amp_grd[j][k] = math.sqrt(xave**2 + yave**2)
                        ph_grd[j][k] = math.atan2(yave, xave)
                        
                        # Handle previous missing bin by interpolation
                        if flag == 1:
                            ph_grd[j][k-1] = (ph_grd[j][k-2] + ph_grd[j][k]) / 2.0
                            ph = ph_grd[j][k-1]
                            amp_grd[j][k-1] = (amp_grd[j][k-2] + amp_grd[j][k]) / 2.0
                            print(f" :intr {line_num} {i} {j} {k-1} {el:.6f} {azgrd:.6f} {amp_grd[j][k-1]:.6f} {ph:.6f}")
                            flag = 0
                        
                        ph = ph_grd[j][k]
                        print(f" : {line_num} {i} {j} {k} {el:.6f} {azgrd:.6f} {amp_grd[j][k]:.6f} {ph:.6f}")
                        fout.write(f"{j} {k} {amp_grd[j][k]:.16e} {ph:.16e}\n")
                    
                    # Move to next azimuth bin
                    azgrd = azgrd - azstep
                    i = 0
                    start_flag = 0
                    
                    k += 1
                    if k > naz - 1:
                        # Completed this elevation row
                        j -= 1
                        end_flag = 1
                        break
                    else:
                        # Initialize next bin
                        amp_grd[j][k] = 0.0
                        ph_grd[j][k] = 0.0
                        sin_ph_grd = 0.0
                        cos_ph_grd = 0.0
                        xave = 0.0
                        yave = 0.0
            
            if el >= elstrt:
                el_prev = el
        
        # Move to next elevation row
        j += 1
    
    fin.close()
    fout.close()
    
    print(f"\nRegridding complete. Output written to {outfile}")
    print(f"Grid dimensions: {nel} × {naz}")
    print(f"Processed {line_num} data points")
    print(f"\nNote: This version uses corrected complex vector averaging.")
    print(f"Amplitude values will differ from the original C code which had a bug.")


def main():
    if len(sys.argv) < 3:
        print("Usage: python regrid.py rawfile regrid.prm")
        print("  rawfile: Input file with columns (el az amplitude phase)")
        print("  regrid.prm: Parameter file with grid specifications")
        sys.exit(1)
    
    rawfile = Path(sys.argv[1])
    param_file = Path(sys.argv[2])
    outfile = Path("rgin.dat")
    
    print(f"Regridding: {sys.argv[0]} {rawfile} {param_file}")
    print(f"Using CORRECTED complex vector averaging (fixes C code bug)\n")
    
    if not rawfile.exists():
        print(f"Error: Raw data file '{rawfile}' not found")
        sys.exit(1)
    
    if not param_file.exists():
        print(f"Error: Parameter file '{param_file}' not found")
        sys.exit(1)
    
    regrid_data(rawfile, param_file, outfile)


if __name__ == '__main__':
    main()
