#!/usr/bin/env python3
"""
Nimesh Patel
October 2025
Python version of previous fortran/C code from Holis package.

Phase unwrapping code for aperture electric field.

Logic:
For every line (n points) of data,
- if phase = -9999.0 (missing data), leave it.
- if phase jumps by > limit radians, set phase = phase + 2*pi

Original C code: 04/dec/96, tks
Python translation with improvements
"""

import sys
import math
import argparse


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Unwrap phase data from input file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                              # Run with all defaults
  %(prog)s -i mydata.dat -o output.dat  # Custom input/output files
  %(prog)s -d 256 -l 4.0                # Custom dimension and limit
  %(prog)s --validate                   # Interactive validation mode
        """
    )
    parser.add_argument('-i', '--input', dest='inputfile', default='Ep.dat',
                        help='Input file path (default: Ep.dat)')
    parser.add_argument('-o', '--output', dest='outputfile', default='tk.dat',
                        help='Output file path (default: tk.dat)')
    parser.add_argument('-d', '--dim', dest='ndim', type=int, default=128,
                        help='Dimension of the data grid (n x n) (default: 128)')
    parser.add_argument('-l', '--limit', type=float, default=3.5,
                        help='Limit for phase jump detection in radians (default: 3.5)')
    parser.add_argument('-v', '--validate', action='store_true',
                        help='Enable interactive validation mode (default: False)')

    args = parser.parse_args()

    # Use math.pi for precise value instead of hardcoding
    pi = math.pi

    print(f"Input file: {args.inputfile}")
    print(f"Output file: {args.outputfile}")
    print(f"Dimension: {args.ndim}")
    print(f"Limit: {args.limit}")
    print(f"Validation mode: {args.validate}")

    limit = args.limit
    ans1 = 1 if args.validate else 0

    n = args.ndim
    phase_begin = -9999.0
    flag1 = 0

    # Open files
    with open(args.inputfile, 'r') as fpi, open(args.outputfile, 'w') as fpo:
        # Process n x n grid
        for j in range(n):
            step = 0.0
            ans = 0
            phase_old = phase_begin
            flag = 1

            for i in range(n):
                # Read phase value
                line = fpi.readline().strip()
                if not line:
                    break

                phase = float(line)

                # Handle missing data
                if phase == -9999.0:
                    fpo.write("-9999.0000\n")
                else:
                    if flag1 == 0:
                        flag1 = 1
                    else:
                        change = phase - phase_old
                        abschange = abs(change)

                        # Check if phase jump exceeds limit
                        if abschange > limit:
                            ans = 1
                            if ans1 == 1:
                                # Interactive mode: ask user
                                user_input = input(f"\n {i} {j} Want to step (1/0)? ")
                                ans = int(user_input)

                            if ans == 1:
                                # Add or subtract 2*pi based on direction of jump
                                step = step - (abschange / change) * 2 * pi
                                ans = 0

                    phase_old = phase
                    fpo.write(f"{phase + step:.6f}\n")

                    if flag == 1:
                        phase_begin = phase + step
                        flag = 0


if __name__ == "__main__":
    main()
