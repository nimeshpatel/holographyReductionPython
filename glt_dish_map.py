# gltDishMap.py
#
# Nimesh Patel
# October 2025
# Python version of previous fortran/C code from Holis package.
# Original code from 2023-Mar-09
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
#  --cmap CMAP      Color table (default: coolwarm; can use seismic, bwr, or similar).
#  --label LABEL    Text for color wedge
#  --holo-smooth    nearest , bilinear (default), bicubic 
#  --holo-filterrad 0.8 (default 1)
#  --holo-filternorm True (default False); Omit to keep interpolation gentler
#
# 2026-May: Added interactive cursor tracking to display surface error values


import numpy
import matplotlib.pyplot as pyplot
from scipy.interpolate import CloughTocher2DInterpolator
from matplotlib.patches import Polygon
import argparse
import os


def cos_alpha_factor(r_m, fprim):
    """
    Baars Eq. 57 projection factor 1/cos(alpha) that converts the half-path-length
    error e (= phase * lambda/4pi) into the normal surface displacement d:

        d = e / cos(alpha),   cos(alpha) = 1 / sqrt(1 + r^2 / (4 f^2))

    so  1/cos(alpha) = sqrt(1 + r^2 / (4 f^2)).

    Args:
        r_m: radius from dish axis in METERS (array or scalar)
        fprim: primary focal length in METERS

    Returns:
        Projection factor 1/cos(alpha) (>= 1), same shape as r_m.
    """
    return numpy.sqrt(1.0 + (r_m ** 2) / (4.0 * fprim ** 2))


def baars_illum_weight(r_m, dprim, taper_db):
    """
    Parabolic illumination amplitude weight, generalized from Baars Eq. 56.

    Baars Eq. 56 (ALMA, 12 dB power taper, 6 m radius):
        w(r) = 1 - (1 - 10^(-0.6)) (r/6)^2      -> w(edge)=0.251

    Generalized for an arbitrary power taper T (dB) and dish radius R0 = dprim/2:
        w(r) = 1 - (1 - 10^(-T/20)) (r/R0)^2

    (T is a power taper in dB; amplitude edge value is 10^(-T/20).)

    Args:
        r_m: radius in METERS
        dprim: primary diameter in METERS (R0 = dprim/2)
        taper_db: power taper at the dish edge, in dB (e.g. 12.0)

    Returns:
        Amplitude weight array, clipped at 0.
    """
    R0 = dprim / 2.0
    edge_amp = 10.0 ** (-taper_db / 20.0)
    w = 1.0 - (1.0 - edge_amp) * (r_m / R0) ** 2
    return numpy.clip(w, 0.0, None)


def _weighted_rms(values, weights):
    """Weighted RMS about the weighted mean."""
    wsum = numpy.sum(weights)
    if wsum <= 0:
        return float(numpy.std(values))
    mean = numpy.average(values, weights=weights)
    return float(numpy.sqrt(numpy.average((values - mean) ** 2, weights=weights)))


def calculate_holography_rms(error_grid, xgrid, ygrid, illum_grid=None,
                             illum_power=2, dprim=12.0,
                             baars_taper_db=None):
    """
    Calculate multiple RMS values for holography surface error data.

    NOTE on the quantity being summed: if the caller has already applied the
    cos-alpha projection (see load_surface_file / cos_alpha_factor), then
    error_grid holds the normal surface displacement d_i and all RMS values
    below are displacement RMS. If not, they are half-path-length (e_i) RMS.

    Args:
        error_grid: 2D array of surface errors in microns (may be masked array)
        xgrid: 2D array of x coordinates (mm), centered at dish center
        ygrid: 2D array of y coordinates (mm), centered at dish center
        illum_grid: Optional 2D array of illumination amplitude (e.g. Ea_um.dat)
        illum_power: Exponent on the illumination amplitude for the file-based
                     illum weighting. 2 = power weighting (original behavior),
                     1 = amplitude weighting (Baars Eq. 54 convention).
        dprim: primary diameter in meters (for the analytic Baars weight)
        baars_taper_db: if not None, also compute the Baars analytic
                        illumination-weighted RMS at this power taper (dB)

    Returns:
        Dictionary with RMS values in microns (missing ones are None)
    """
    # Handle both masked arrays and regular arrays with NaN
    if hasattr(error_grid, 'mask'):
        data = numpy.ma.filled(error_grid, numpy.nan)
        valid = ~numpy.isnan(data) & ~error_grid.mask
    else:
        data = error_grid
        valid = ~numpy.isnan(data)

    out = {'standard': 0, 'area_weighted': 0,
           'illum_weighted': None, 'baars_illum_weighted': None}

    if numpy.sum(valid) == 0:
        return out

    error_valid = data[valid]
    x_valid = xgrid[valid]
    y_valid = ygrid[valid]
    r_mm = numpy.sqrt(x_valid ** 2 + y_valid ** 2)
    r_m = r_mm / 1000.0

    # 1. Standard RMS (per-pixel equal weight)
    out['standard'] = float(numpy.std(error_valid))

    # 2. Area-weighted RMS (weight proportional to r^2)
    r2_weights = r_mm ** 2
    out['area_weighted'] = (_weighted_rms(error_valid, r2_weights)
                            if numpy.sum(r2_weights) > 0 else out['standard'])

    # 3. Illumination-weighted RMS from a supplied amplitude map (Ea_um.dat)
    if illum_grid is not None:
        illum_valid = illum_grid[valid]
        illum_weights = numpy.abs(illum_valid) ** illum_power
        if numpy.sum(illum_weights) > 0:
            out['illum_weighted'] = _weighted_rms(error_valid, illum_weights)
        else:
            out['illum_weighted'] = out['standard']

    # 4. Baars analytic illumination-weighted RMS (Eq. 54 + Eq. 56),
    #    amplitude weighting per Baars.
    if baars_taper_db is not None:
        w = baars_illum_weight(r_m, dprim, baars_taper_db)
        out['baars_illum_weighted'] = _weighted_rms(error_valid, w)

    return out

def _prm_value(line):
    """
    Extract the numeric/string value from a holography .prm line.

    The Holis convention puts a dotted label in the first ~49 columns and the
    value afterwards, but spacing varies between files. Try the fixed column
    first, then fall back to the last whitespace-separated token.
    """
    val = line[49:].strip() if len(line) >= 50 else ""
    if not val:
        toks = line.split()
        val = toks[-1] if toks else ""
    return val


def read_holography_params(prm_file='withphase_aber.prm'):
    """Read parameters from holography parameter file."""
    params = {
        'rate': 0.7521,
        'dprim': 12.0,
        'freq': 94.5,
        'grid_size': 128,
        'fprim': None,   # Primary focal length (m) - needed for cos-alpha projection
        'fmag': None,    # Cassegrain magnification (informational)
    }

    if not os.path.exists(prm_file):
        print(f"  Warning: {prm_file} not found, using default parameters")
        return params

    try:
        with open(prm_file, 'r') as f:
            for raw in f:
                line = raw.rstrip('\n')
                low = line.lower()
                if low.strip().startswith('!'):
                    continue
                try:
                    if 'nyquist sampling rate' in low:
                        params['rate'] = float(_prm_value(line))
                    elif 'diameter of the primary' in low or ('diameter of primary' in low):
                        params['dprim'] = float(_prm_value(line))
                    elif 'observing frequency' in low:
                        params['freq'] = float(_prm_value(line))
                    elif 'size n of the n by n' in low:
                        params['grid_size'] = int(_prm_value(line))
                    elif 'focal length' in low and 'primary' in low:
                        params['fprim'] = float(_prm_value(line))
                    elif 'magnification' in low:
                        params['fmag'] = float(_prm_value(line))
                except ValueError:
                    continue
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

        # Auto-detect grid size from data if possible, otherwise use prm file
        data_size = len(data)
        detected_grid_size = int(numpy.sqrt(data_size))
        if detected_grid_size * detected_grid_size == data_size:
            grid_size = detected_grid_size
            if grid_size != params['grid_size']:
                print(f"  Auto-detected grid size from data: {grid_size} x {grid_size}")
        else:
            grid_size = params['grid_size']
        rate = params['rate']
        dprim = params['dprim']
        freq_ghz = params['freq']
        # Make params available to the caller (for dprim / fprim / taper in RMS)
        try:
            args._holo_params = params
        except Exception:
            pass

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
        
        # Calculate shifts - auto-shift scales with pixel size
        # Reference: -65mm x-shift and +65mm y-shift work for 128x128 (124.65 mm pixels)
        # This corresponds to ~0.52 pixel offset
        if hasattr(args, 'auto_shift') and args.auto_shift:
            pixel_offset = 0.52  # pixels (calibrated for 128x128 at -65mm)
            x_shift = -pixel_offset * pixel_spacing_mm
            y_shift = pixel_offset * pixel_spacing_mm
            print(f"  Auto-shift enabled: x={x_shift:.1f} mm, y={y_shift:.1f} mm")
        else:
            x_shift = args.x_shift
            y_shift = args.y_shift
            if x_shift != 0 or y_shift != 0:
                print(f"  Manual shift: x={x_shift:.1f} mm, y={y_shift:.1f} mm")
        
        x = numpy.linspace(-half_size, half_size, grid_size) + x_shift
        y = numpy.linspace(-half_size, half_size, grid_size) + y_shift

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

        # Also mask any -9999 flag values (invalid data marker)
        flag_mask = (error_grid < -9000)
        if numpy.sum(flag_mask) > 0:
            print(f"  Flagged pixels (-9999): {numpy.sum(flag_mask)}")
            mask = mask | flag_mask

        error = error_grid * error_scaling

        # --- cos-alpha projection (Baars Eq. 55/57): path-length error e -> normal displacement d ---
        # d = e / cos(alpha),  1/cos(alpha) = sqrt(1 + r^2/(4 f^2))
        apply_cos = getattr(args, 'cos_alpha', True)
        fprim = getattr(args, 'fprim', None) or params.get('fprim', None)
        if apply_cos:
            if fprim and fprim > 0:
                # Unshifted, dish-centered radius in meters (display shift must NOT enter here)
                xg0 = numpy.linspace(-half_size, half_size, grid_size)  # mm
                yg0 = numpy.linspace(-half_size, half_size, grid_size)  # mm
                Xg0, Yg0 = numpy.meshgrid(xg0, yg0)
                r_m = numpy.sqrt(Xg0 ** 2 + Yg0 ** 2) / 1000.0  # mm -> m
                proj = cos_alpha_factor(r_m, fprim)             # = 1/cos(alpha)
                error = error * proj
                edge_factor = float(cos_alpha_factor(dprim / 2.0, fprim))
                print(f"  cos-alpha projection: ON  (fprim = {fprim:.3f} m, "
                      f"dish-edge 1/cos(alpha) = {edge_factor:.4f})")
            else:
                print("  cos-alpha projection: SKIPPED (no focal length; "
                      "add 'Focal length of the primary' to the .prm or pass --fprim)")
        else:
            print("  cos-alpha projection: OFF (reporting half-path-length error e)")

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


class InteractiveCursor:
    """
    Interactive cursor that displays surface error value at mouse position.
    """
    def __init__(self, ax, im, error_data, x_coords, y_coords, is_grid_data, interp_func=None):
        self.ax = ax
        self.im = im
        self.error_data = error_data
        self.x_coords = x_coords
        self.y_coords = y_coords
        self.is_grid_data = is_grid_data
        self.interp_func = interp_func
        
        # Create text annotation for displaying values
        self.text = ax.text(0.02, 0.98, '', transform=ax.transAxes,
                           fontsize=10, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        # Connect mouse motion event
        self.cid = ax.figure.canvas.mpl_connect('motion_notify_event', self.on_move)
    
    def on_move(self, event):
        if event.inaxes != self.ax:
            self.text.set_text('')
            self.ax.figure.canvas.draw_idle()
            return
        
        x, y = event.xdata, event.ydata
        
        # Calculate radius from center
        r = numpy.sqrt(x**2 + y**2)
        
        # Get the surface error value at this position
        if self.is_grid_data:
            # For grid data, find nearest pixel
            if self.x_coords is not None and self.y_coords is not None:
                # Find nearest indices
                ix = numpy.argmin(numpy.abs(self.x_coords - x))
                
                # For origin='upper', y-axis is inverted:
                # - y_coords goes from -half_size to +half_size (bottom to top in world coords)
                # - But imshow origin='upper' puts array row 0 at top
                # - So we need to invert: row index = (len-1) - iy
                iy_world = numpy.argmin(numpy.abs(self.y_coords - y))
                iy_array = len(self.y_coords) - 1 - iy_world
                
                # Check bounds
                if 0 <= ix < len(self.x_coords) and 0 <= iy_array < len(self.y_coords):
                    value = self.error_data[iy_array, ix]
                    
                    if numpy.ma.is_masked(value):
                        self.text.set_text(f'X: {x:.0f} mm\nY: {y:.0f} mm\nR: {r:.0f} mm\n(masked)')
                    else:
                        self.text.set_text(f'X: {x:.0f} mm\nY: {y:.0f} mm\nR: {r:.0f} mm\nError: {value:.1f} μm')
                else:
                    self.text.set_text(f'X: {x:.0f} mm\nY: {y:.0f} mm\nR: {r:.0f} mm')
            else:
                self.text.set_text(f'X: {x:.0f} mm\nY: {y:.0f} mm\nR: {r:.0f} mm')
        else:
            # For photogrammetry data, use interpolator
            if self.interp_func is not None:
                value = self.interp_func(x, y)
                if numpy.isnan(value):
                    self.text.set_text(f'X: {x:.0f} mm\nY: {y:.0f} mm\nR: {r:.0f} mm\n(outside data)')
                else:
                    self.text.set_text(f'X: {x:.0f} mm\nY: {y:.0f} mm\nR: {r:.0f} mm\nError: {value:.1f} μm')
            else:
                self.text.set_text(f'X: {x:.0f} mm\nY: {y:.0f} mm\nR: {r:.0f} mm')
        
        self.ax.figure.canvas.draw_idle()


def main():
    parser = argparse.ArgumentParser(description="Plot dish surface error map with optional color scale and output file.")

    parser.add_argument('--no-panels', action='store_true',
                        help='Suppress panel boundary overlays on the dish map.')
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
    parser.add_argument("--cmap", default="coolwarm", help="Color table (default: coolwarm)")
    parser.add_argument("--label",  help="Text for color wedge ")
    parser.add_argument("--pg-recenter", choices=["none", "mean", "median"], default="median",
                help="How to recenter photogrammetry x,y before interpolation (default: median).")
    parser.add_argument("--x-shift", type=float, default=0.0,
                help="Extra x shift in mm (positive shifts the map right).")
    parser.add_argument("--y-shift", type=float, default=0.0,
                help="Extra y shift in mm (positive shifts the map up).")
    parser.add_argument("--auto-shift", action="store_true",
                help="Automatically calculate shift based on grid size (overrides --x-shift and --y-shift).")
    parser.add_argument("--no-interactive", action="store_true",
                help="Disable interactive cursor tracking (useful for batch processing).")
    parser.add_argument("--illum-file", 
                help="Illumination amplitude file (Ea_um.dat) for illumination-weighted RMS.")
    parser.add_argument("--no-cos-alpha", dest="cos_alpha", action="store_false",
                help="Disable the cos-alpha projection (Baars Eq.55/57). "
                     "Default: projection ON, reporting normal surface displacement.")
    parser.set_defaults(cos_alpha=True)
    parser.add_argument("--fprim", type=float, default=None,
                help="Primary focal length (m) override for cos-alpha projection "
                     "(otherwise read from the .prm file).")
    parser.add_argument("--illum-power", type=int, choices=[1, 2], default=2,
                help="Exponent on illumination amplitude for --illum-file weighting: "
                     "2 = power (default, original), 1 = amplitude (Baars Eq.54).")
    parser.add_argument("--baars-illum", action="store_true",
                help="Also report the Baars analytic illumination-weighted RMS "
                     "(Eq.54 + Eq.56), using the parabolic taper model.")
    parser.add_argument("--taper-db", type=float, default=12.0,
                help="Edge power taper (dB) for the Baars analytic illumination weight "
                     "(default 12.0, matching ALMA).")

    args = parser.parse_args()

    xtarg, ytarg, ztarg, error, is_grid_data, error_grid, xgrid_raw, ygrid_raw = load_surface_file(args.filename,args)

    fig, ax = pyplot.subplots(figsize=(8, 7))
    cmap = pyplot.get_cmap(args.cmap)

    interp_func = None  # For photogrammetry data interpolation

    if is_grid_data:
        extent = (xgrid_raw.min(), xgrid_raw.max(), ygrid_raw.min(), ygrid_raw.max())
        im = ax.imshow(error_grid, extent=extent, origin='upper', cmap=cmap, aspect='equal', vmin=args.vmin, vmax=args.vmax, interpolation=(args.holo_smooth or 'nearest'), filternorm=args.holo_filternorm, filterrad=args.holo_filterrad)
    else:
        xp = xtarg
        yp = ytarg
        zp = error
        interp_func = CloughTocher2DInterpolator(list(zip(xp, yp)), zp)
        grid_x = numpy.linspace(min(xp), max(xp), 300)
        grid_y = numpy.linspace(min(yp), max(yp), 300)
        xgrid, ygrid = numpy.meshgrid(grid_x, grid_y)
        zgrid = interp_func(xgrid, ygrid)
        im = ax.imshow(zgrid, extent=(grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()),
                       origin='lower', cmap=cmap, aspect='equal', vmin=args.vmin, vmax=args.vmax)

    pyplot.colorbar(im, label=args.label)
    rmin, rmax, npanels = parse_panel_geometry("panelplt.prm")
    if not args.no_panels:
        draw_panel_boundaries(ax, rmin, rmax, npanels)

    # Calculate and display RMS
    title_text = args.filename
    
    if is_grid_data:
        # Holography data - calculate multiple RMS values
        # Load illumination data if provided
        illum_grid = None
        if args.illum_file and os.path.exists(args.illum_file):
            try:
                illum_data = numpy.loadtxt(args.illum_file)
                grid_size = int(numpy.sqrt(len(illum_data)))
                illum_grid = illum_data.reshape((grid_size, grid_size))
                # Apply same masking as error_grid
                illum_grid = numpy.where(numpy.isnan(error_grid), numpy.nan, illum_grid)
            except Exception as e:
                print(f"Warning: Could not load illumination file: {e}")
        
        # Create coordinate grids for RMS calculation
        # IMPORTANT: Use coordinates centered at (0,0), NOT shifted coordinates
        # The x_shift/y_shift are for display alignment, not for calculating r²
        grid_size = error_grid.shape[0]
        half_extent = (xgrid_raw.max() - xgrid_raw.min()) / 2
        x_centered = numpy.linspace(-half_extent, half_extent, grid_size)
        y_centered = numpy.linspace(half_extent, -half_extent, grid_size)  # Flipped for image coords
        xgrid_2d, ygrid_2d = numpy.meshgrid(x_centered, y_centered)
        
        rms_values = calculate_holography_rms(
            error_grid, xgrid_2d, ygrid_2d, illum_grid,
            illum_power=args.illum_power,
            dprim=getattr(args, '_holo_params', {}).get('dprim', 12.0) if hasattr(args, '_holo_params') else 12.0,
            baars_taper_db=(args.taper_db if args.baars_illum else None),
        )

        # Label reflects whether values are normal displacement (projection on) or path-length
        qty = "normal displacement" if getattr(args, 'cos_alpha', True) else "half-path-length"

        # Print RMS to console (standard/pixel-based only)
        print(f"\nSurface RMS:  {rms_values['standard']:.2f} µm")

        # Write JSON sidecar so create_summary_pdf.py can read the correct RMS
        if args.output:
            import json
            sidecar = os.path.splitext(args.output)[0] + "_rms.json"
            rms_out = {'standard': round(rms_values['standard'], 4)}
            # Merge fit results if available (written by holis_aber2.py)
            fit_json = os.path.splitext(args.output)[0] + "_fit.json"
            if not os.path.exists(fit_json):
                base = os.path.splitext(args.output)[0]
                fit_json = base.replace('_Epr', '') + "_fit.json" if '_Epr' in base else fit_json
            if os.path.exists(fit_json):
                try:
                    with open(fit_json) as _f:
                        fit_data = json.load(_f)
                    rms_out['fit'] = fit_data
                except Exception:
                    pass
            with open(sidecar, 'w') as _f:
                json.dump(rms_out, _f, indent=2)
            print(f"  RMS sidecar saved to {sidecar}")

        # Title on the plot: standard RMS only
        title_text += f"\nRMS: {rms_values['standard']:.1f} µm"
    else:
        # Photogrammetry data - standard RMS only
        rms_error = numpy.std(error)
        title_text += f'\nRMS = {rms_error:.1f} μm'
    
    ax.set_title(title_text)
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_xlim(-6100, 6100)
    ax.set_ylim(-6100, 6100)
    fig.tight_layout()

    # Set up interactive cursor (unless disabled or saving to file)
    if not args.no_interactive and not args.output:
        cursor = InteractiveCursor(
            ax, im, 
            error_grid if is_grid_data else None,
            xgrid_raw if is_grid_data else None,
            ygrid_raw if is_grid_data else None,
            is_grid_data,
            interp_func
        )
        print("\nInteractive mode: Move cursor over the plot to see surface error values.")
        print("Close the window to exit.\n")

    if args.output:
        # Auto-detect format from file extension
        if args.output.lower().endswith('.png'):
            pyplot.savefig(args.output, format='png', dpi=150)
        elif args.output.lower().endswith('.pdf'):
            pyplot.savefig(args.output, format='pdf')
        else:
            # Default to the extension or PDF if unknown
            pyplot.savefig(args.output)
        print(f"Saved figure to {args.output}")
    else:
        pyplot.show()

if __name__ == "__main__":
    main()
