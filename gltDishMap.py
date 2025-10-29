# gltDishMap.py
#
# Nimesh Patel, original code from 2023-Mar-09
#
# This is a modified version of Dick Plambeck's plot_surface.py code
# for making a surface error map based on photogrammetry measurements.
# This version is for the GLT antenna. Panel geometry and dimensions are
# from Vertex drawings provided by Philippe Raffin. 
# The sample input file Parajob10rotated.txt is provided by Teddy Huang, from
# a previous photogrammetry measurement at Thule. (A rotation correction
# of 1.46 deg is applied to the x,y positions of targets).
#
# Revised on 2025-Jun-07, updated for GLT holography maps.
# This same code can be used for plotting both photogrammetry and
# holography surface error maps (difference is recognized by shape of
# input file. For holography maps, the data are assumed to be from a
# 128x128 map. 
# usage: gltDishMap.py [-h] [--output OUTPUT] [--vmin VMIN] [--vmax VMAX]
#                     [--cmap CMAP] [--label LABEL]
#                     filename

#Plot dish surface error map with optional color scale and output file.
#
#positional arguments:
#  filename         Input data file (1-column or 4-column)
#
#optional arguments:
#  -h, --help       show this help message and exit
#  --output OUTPUT  Save figure to this PDF file
#  --vmin VMIN      Minimum value for color scale (in microns)
#  --vmax VMAX      Maximum value for color scale (in microns)
#  --cmap CMAP      Color table (use seismic, bwr, or similar with opposite  colors around 0).
#  --label LABEL    Text for color wedge
#  --holo-smooth    nearest , bilinear (default), bicubic 
#  --holo-filterrad 0.8 (default 1)
#  --holo-filternorm True (default False); Omit to keep interpolation gentler


import numpy
import matplotlib.pyplot as pyplot
from scipy.interpolate import CloughTocher2DInterpolator
from matplotlib.patches import Polygon
import argparse
import os

def read_holography_params(prm_file='withphase_aber.prm'):
    """Read parameters from holography parameter file."""
    params = {
        'rate': 0.7521,
        'dprim': 12.0,
        'freq': 94.5,
        'grid_size': 128
    }

    if not os.path.exists(prm_file):
        print(f"  Warning: {prm_file} not found, using default parameters")
        return params

    try:
        with open(prm_file, 'r') as f:
            for line in f:
                if 'Nyquist sampling rate' in line:
                    params['rate'] = float(line[49:].strip())
                elif 'Diameter of primary' in line and 'secondary' not in line:
                    params['dprim'] = float(line[49:].strip())
                elif 'Observing frequency' in line:
                    params['freq'] = float(line[49:].strip())
                elif 'Size N of the N by N' in line:
                    params['grid_size'] = int(line[49:].strip())
    except Exception as e:
        print(f"  Warning: Error reading {prm_file}: {e}")
        print(f"  Using default parameters")

    return params

def load_surface_file(filename, args):
    data = numpy.loadtxt(filename)
    if data.ndim == 1:
        # Read holography parameters from withphase_aber.prm
        prm_file = args.prm_file if hasattr(args, 'prm_file') and args.prm_file else 'withphase_aber.prm'
        params = read_holography_params(prm_file)

        grid_size = params['grid_size']
        rate = params['rate']
        dprim = params['dprim']
        freq_ghz = params['freq']

        # Calculate correct pixel spacing from holography parameters
        # Field of view = dprim / rate, then divide by grid_size
        fov_meters = dprim / rate
        pixel_spacing_mm = (fov_meters * 1000.0) / grid_size  # Convert to mm

        # Calculate error scaling: wavelength/2 in microns
        # wavelength = c/f = 299792458 / (94.5e9) = 3.171 mm = 3171 microns
        # Surface error = wavelength / (4*pi) * phase_error_radians
        wavelength_microns = 299792458.0 / (freq_ghz * 1e9) * 1e6
        error_scaling = wavelength_microns / (4.0 * numpy.pi)

        print(f"Holography parameters:")
        print(f"  Grid size: {grid_size} x {grid_size}")
        print(f"  Pixel spacing: {pixel_spacing_mm:.3f} mm")
        print(f"  Field of view: {fov_meters:.3f} m")
        print(f"  Error scaling: {error_scaling:.3f} microns/radian")

        error_grid = data.reshape((grid_size, grid_size))
        half_size = (grid_size * pixel_spacing_mm) / 2
        x = numpy.linspace(-half_size, half_size, grid_size) + args.x_shift
        y = numpy.linspace(-half_size, half_size, grid_size) + args.y_shift

        # Try to read the mask file from holography processing
        mask_file = args.mask_file if hasattr(args, 'mask_file') and args.mask_file else "mask32.dat"
        try:
            mask_data = numpy.loadtxt(mask_file, dtype=int)
            mask_grid = mask_data.reshape((grid_size, grid_size))
            # In holography mask: 1=valid, 0=masked
            # For matplotlib masked array: True=masked, False=valid
            mask = (mask_grid == 0)
            print(f"  Using mask from: {mask_file}")
            print(f"  Masked pixels: {numpy.sum(mask)} of {grid_size*grid_size}")
        except (FileNotFoundError, ValueError) as e:
            print(f"  Warning: Could not read {mask_file}, using simple circular mask")
            # Fallback to simple circular mask
            xtarg, ytarg = numpy.meshgrid(x, y)
            r = numpy.sqrt(xtarg**2 + ytarg**2)
            mask = (r > 6000) | (r < 375)  # Mask outside 6m and inside subreflector

        error = error_grid * error_scaling
        error_masked = numpy.ma.array(error, mask=mask)

        ztarg = numpy.zeros((grid_size, grid_size))
        return None, None, ztarg, None, True, error_masked, x, y
    else:
        xtarg, ytarg, ztarg, error = numpy.loadtxt(filename, unpack=True)
        # --- recentre PG points so the dish is centered on (0,0) like the panel overlay ---
        if args.pg_recenter == "mean":
            cx, cy = float(numpy.mean(xtarg)), float(numpy.mean(ytarg))
        elif args.pg_recenter == "median":
            cx, cy = float(numpy.median(xtarg)), float(numpy.median(ytarg))
        else:  # "none"
            cx, cy = 0.0, 0.0

        xtarg = xtarg - cx + args.x_shift
        ytarg = ytarg - cy + args.y_shift

        error *= 1000.  # mm to microns (only for photogrammetry)
        return xtarg, ytarg, ztarg, error, False, None, None, None

def parse_panel_geometry(filename="panelplt.prm"):
    rmin = []
    rmax = []
    npanels = []
    radii = []
    panel_counts = []
    with open(filename, "r") as f:
        lines = f.readlines()
    for line in lines:
        if "Radius of the" in line:
            try:
                value = float(line[49:].strip()) * 1000
                radii.append(value)
            except ValueError:
                continue
        elif "Number of panels in the" in line:
            try:
                value = int(line[49:].strip())
                panel_counts.append(value)
            except ValueError:
                continue
    if len(radii) >= 2:
        rmin = radii[:-1]
        rmax = radii[1:]
    else:
        raise ValueError("Not enough radii extracted to form rmin/rmax.")
    npanels = panel_counts
    return rmin, rmax, npanels

def draw_panel_boundaries(ax, rmin, rmax, npanels):
    for i in range(len(npanels)):
        inner = rmin[i]
        outer = rmax[i]
        n = npanels[i]
        for j in range(n):
            theta1 = 2 * numpy.pi * j / n
            theta2 = 2 * numpy.pi * (j + 1) / n
            x0 = inner * numpy.cos(theta1)
            y0 = inner * numpy.sin(theta1)
            x1 = outer * numpy.cos(theta1)
            y1 = outer * numpy.sin(theta1)
            x2 = outer * numpy.cos(theta2)
            y2 = outer * numpy.sin(theta2)
            x3 = inner * numpy.cos(theta2)
            y3 = inner * numpy.sin(theta2)
            panel = Polygon([[x0, y0], [x1, y1], [x2, y2], [x3, y3]], closed=True, fill=False, edgecolor='black', linewidth=0.3)
            ax.add_patch(panel)

def main():
    parser = argparse.ArgumentParser(description="Plot dish surface error map with optional color scale and output file.")

    parser.add_argument('--holo-smooth', nargs='?', const='bilinear', choices=['nearest','bilinear','bicubic'],
                        help='Holography display interpolation. If provided without value, uses bilinear. If nearest: no interpolation; bicubic: too smooth')
    parser.add_argument('--holo-filterrad', type=float, default=1.0,
                        help='Interpolation filter radius (imshow filterrad); lower is less smoothing, try 0.8 (default 1.0).')
    parser.add_argument('--holo-filternorm', action='store_true',
                        help='Use filternorm=True for imshow. Omit to keep it gentler (False by default).')
    parser.add_argument('--mask-file', default='mask32.dat',
                        help='Mask file from holography processing (default: mask32.dat)')
    parser.add_argument('--prm-file', default='withphase_aber.prm',
                        help='Parameter file from holography processing (default: withphase_aber.prm)')
    parser.add_argument("filename", help="Input data file (1-column or 4-column)")
    parser.add_argument("--output", help="Save figure to this PDF file")
    parser.add_argument("--vmin", type=float, help="Minimum value for color scale (in microns)")
    parser.add_argument("--vmax", type=float, help="Maximum value for color scale (in microns)")
    parser.add_argument("--cmap",  help="Color table ")
    parser.add_argument("--label",  help="Text for color wedge ")
    parser.add_argument("--pg-recenter", choices=["none", "mean", "median"], default="median",
                help="How to recenter photogrammetry x,y before interpolation (default: median).")
    parser.add_argument("--x-shift", type=float, default=0.0,
                help="Extra x shift in mm (positive shifts the map right).")
    parser.add_argument("--y-shift", type=float, default=0.0,
                help="Extra y shift in mm (positive shifts the map up).")

    args = parser.parse_args()

    xtarg, ytarg, ztarg, error, is_grid_data, error_grid, xgrid_raw, ygrid_raw = load_surface_file(args.filename,args)

    fig, ax = pyplot.subplots(figsize=(8, 7))
    cmap = pyplot.get_cmap(args.cmap)

    if is_grid_data:
        extent = (xgrid_raw.min(), xgrid_raw.max(), ygrid_raw.min(), ygrid_raw.max())
        im = ax.imshow(error_grid, extent=extent, origin='upper', cmap=cmap, aspect='equal', vmin=args.vmin, vmax=args.vmax, interpolation=(args.holo_smooth or 'nearest'), filternorm=args.holo_filternorm, filterrad=args.holo_filterrad)
    else:
        xp = xtarg
        yp = ytarg
        zp = error
        interp = CloughTocher2DInterpolator(list(zip(xp, yp)), zp)
        grid_x = numpy.linspace(min(xp), max(xp), 300)
        grid_y = numpy.linspace(min(yp), max(yp), 300)
        xgrid, ygrid = numpy.meshgrid(grid_x, grid_y)
        zgrid = interp(xgrid, ygrid)
        im = ax.imshow(zgrid, extent=(grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()),
                       origin='lower', cmap=cmap, aspect='equal', vmin=args.vmin, vmax=args.vmax)

    pyplot.colorbar(im, label=args.label)
    rmin, rmax, npanels = parse_panel_geometry("panelplt.prm")
    draw_panel_boundaries(ax, rmin, rmax, npanels)

    ax.set_title(args.filename)
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_xlim(-6100, 6100)
    ax.set_ylim(-6100, 6100)
    fig.tight_layout()

    if args.output:
        pyplot.savefig(args.output, format='pdf')
        print(f"Saved figure to {args.output}")
    else:
        pyplot.show()

if __name__ == "__main__":
    main()
