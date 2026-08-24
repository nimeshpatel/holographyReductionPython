#!/bin/bash
#
# holoMap.sh - GLT Holography Data Reduction Pipeline
#
# Usage: ./holoMap.sh <input_file> <size> [boresight_mode] [options]
#
#   input_file:     Raw holography data file
#   size:           Map size (32, 64, or 128)
#   boresight_mode: 0 = no boresight cal data in file (default)
#                   1 = apply boresight calibration (pre-regrid, boresight_cal.py)
#                   2 = boresight data present, but skip correction (just remove boresight points)
#                   3 = apply boresight calibration (post-regrid, original ido_bore=1 method)
#   options:        --no-plot     Skip interactive plots (faster over slow network)
#                   --compare     Run both with and without calibration, show comparison
#
# Examples:
#   ./holoMap.sh ../data/holoADC-AzEl-20260522_123456.txt 64
#   ./holoMap.sh ../data/holoADC-AzEl-20260522_123456.txt 128 1
#   ./holoMap.sh ../data/holoADC-AzEl-20260522_123456.txt 128 2      # Skip correction
#   ./holoMap.sh ../data/holoADC-AzEl-20260522_123456.txt 128 3      # Original ido_bore method
#   ./holoMap.sh ../data/holoADC-AzEl-20260522_123456.txt 128 1 --no-plot
#   ./holoMap.sh ../data/holoADC-AzEl-20260522_123456.txt 128 1 --compare
#   ./holoMap.sh ../data/holoADC-AzEl-20260522_123456.txt 128 1 --no-plot --compare
#

set -e  # Exit on error

# Activate the nimesh_holo environment
source $HOME/.mamba_rc
mamba activate nimesh_holo

# Parse arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_file> <size> [boresight_mode] [--no-plot] [--compare]"
    echo ""
    echo "  boresight_mode: 0 = no boresight data (default)"
    echo "                  1 = apply boresight calibration (pre-regrid)"
    echo "                  2 = boresight data present, skip correction"
    echo "                  3 = apply boresight calibration (post-regrid, ido_bore)"
    echo ""
    echo "  --no-plot    Skip interactive plots (faster over slow network)"
    echo "  --compare    Compare with and without calibration side-by-side"
    exit 1
fi

INPUTFILE="$1"
SIZE="$2"
DO_BORESIGHT_CAL="${3:-0}"

# Parse optional flags
NO_PLOT=0
DO_COMPARE=0
shift 2  # Remove first two positional args
if [ $# -gt 0 ] && [[ "$1" =~ ^[0-9]+$ ]]; then
    shift  # Remove boresight_mode if it was provided
fi
while [ $# -gt 0 ]; do
    case "$1" in
        --no-plot)
            NO_PLOT=1
            shift
            ;;
        --compare)
            DO_COMPARE=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Extract prefix from input filename
BASENAME=$(basename "$INPUTFILE" .txt)
PREFIX="${BASENAME}"

echo "=============================================="
echo "GLT Holography Reduction Pipeline"
echo "=============================================="
echo "Input file: $INPUTFILE"
echo "Map size: ${SIZE}x${SIZE}"
case "$DO_BORESIGHT_CAL" in
    0) echo "Boresight calibration: OFF (no boresight data)" ;;
    1) echo "Boresight calibration: pre-regrid (boresight_cal.py)" ;;
    2) echo "Boresight calibration: skip correction (remove points only)" ;;
    3) echo "Boresight calibration: post-regrid (original ido_bore method)" ;;
    *) echo "Boresight calibration mode: $DO_BORESIGHT_CAL" ;;
esac
echo "Skip plots: $NO_PLOT"
echo "Compare mode: $DO_COMPARE"
echo "Output prefix: $PREFIX"
echo ""

# Set up parameter file symlinks
ln -sf "regrid_${SIZE}x${SIZE}.prm" regrid.prm
ln -sf "preprocess_${SIZE}.prm" preprocess.prm
ln -sf "withphase_aber_${SIZE}.prm" withphase_aber.prm

# Build detect_raster_start options
DETECT_OPTS=""
if [ "$NO_PLOT" == "1" ]; then
    DETECT_OPTS="--no-plot"
fi

#----------------------------------------------
# Step 1: Detect raster start and trim data
#----------------------------------------------
echo "Step 1: Detecting raster start..."

if [ "$DO_BORESIGHT_CAL" == "1" ] || [ "$DO_BORESIGHT_CAL" == "2" ] || [ "$DO_BORESIGHT_CAL" == "3" ] || [ "$DO_COMPARE" == "1" ]; then
    python detect_raster_start.py "$INPUTFILE" -o trimmed_with_time.txt --keep-timestamps $DETECT_OPTS
else
    python detect_raster_start.py "$INPUTFILE" -o trimmed.txt $DETECT_OPTS
fi

#----------------------------------------------
# Step 2: Boresight calibration (if applicable)
#----------------------------------------------
# Build boresight_cal.py options
BORE_PLOT_OPT=""
if [ "$NO_PLOT" == "0" ]; then
    BORE_PLOT_OPT="--plot"
fi

if [ "$DO_BORESIGHT_CAL" == "1" ]; then
    echo ""
    echo "Step 2: Applying boresight calibration (pre-regrid method)..."
    python boresight_cal.py trimmed_with_time.txt -o calibrated.txt --verbose --remove-boresight \
        --min-amp 0.35 --slew-margin 2.0 --median-window 5 --smoothing 128 \
        --save-boresight-data --save-prefix "$PREFIX" $BORE_PLOT_OPT
    
    awk '!/^#/ {print $2, $3, $4, $5}' calibrated.txt > trimmed.txt
    
elif [ "$DO_BORESIGHT_CAL" == "2" ]; then
    echo ""
    echo "Step 2: Removing boresight points (no phase correction)..."
    python boresight_cal.py trimmed_with_time.txt -o calibrated.txt --verbose --remove-boresight \
        --min-amp 0.35 --slew-margin 2.0 --no-correction \
        --save-boresight-data --save-prefix "$PREFIX" $BORE_PLOT_OPT
    
    awk '!/^#/ {print $2, $3, $4, $5}' calibrated.txt > trimmed.txt

elif [ "$DO_BORESIGHT_CAL" == "3" ]; then
    echo ""
    echo "Step 2: Preparing for post-regrid boresight correction (ido_bore method)..."
    
    # Extract boresight data without applying correction
    python boresight_cal.py trimmed_with_time.txt -o calibrated.txt --verbose --remove-boresight \
        --min-amp 0.35 --slew-margin 2.0 --no-correction \
        --save-boresight-data --save-prefix "$PREFIX" $BORE_PLOT_OPT
    
    # Convert boresight data to original bor format
    # Use --negate because preprocess.py applies scale_ph to main data but not to boresight
    BORFILE="bor${SIZE}"
    echo "Converting boresight data to ${BORFILE}..."
    python convert_to_bor.py "${PREFIX}_boresight_data.txt" "$SIZE" -o "$BORFILE" --negate -v
    
    # Extract trimmed.txt without timestamps
    awk '!/^#/ {print $2, $3, $4, $5}' calibrated.txt > trimmed.txt
    
    # Update preprocess.prm to use ido_bore=1
    echo "Updating preprocess.prm for ido_bore=1..."
    sed -i "s/^Bore-sight data file name.*/Bore-sight data file name...............         ${BORFILE}/" preprocess.prm
    sed -i "s/^Do bore-sight drift correction.*/Do bore-sight drift correction (1\/0)....         1/" preprocess.prm

elif [ "$DO_COMPARE" == "1" ]; then
    echo ""
    echo "Step 2: Comparison mode - will run both calibrated and uncalibrated"
    # We'll handle this after the main pipeline
    # For now, run WITH calibration as the primary path
    python boresight_cal.py trimmed_with_time.txt -o calibrated.txt --verbose --remove-boresight \
        --min-amp 0.35 --slew-margin 2.0 --median-window 5 --smoothing 128 \
        --save-boresight-data --save-prefix "$PREFIX" $BORE_PLOT_OPT
    
    awk '!/^#/ {print $2, $3, $4, $5}' calibrated.txt > trimmed.txt
fi

#----------------------------------------------
# Step 3: Regridding
#----------------------------------------------
echo ""
echo "Step 3: Regridding data to a ${SIZE}x${SIZE} map..."
if [ "$NO_PLOT" == "1" ]; then
    python regrid_holo.py trimmed.txt regrid.prm --no-plot
else
    python regrid_holo.py trimmed.txt regrid.prm
fi

#----------------------------------------------
# Step 4: Fix missing grid cells (if boresight cal)
#----------------------------------------------
if [ "$DO_BORESIGHT_CAL" == "1" ] || [ "$DO_BORESIGHT_CAL" == "2" ] || [ "$DO_BORESIGHT_CAL" == "3" ] || [ "$DO_COMPARE" == "1" ]; then
    echo ""
    echo "Step 4: Fixing any missing grid cells..."
    python fix_missing_cell.py "$SIZE" --force
fi

#----------------------------------------------
# Step 5: Preprocessing
#----------------------------------------------
echo ""
echo "Step 5: Preprocessing..."
python preprocess.py

#----------------------------------------------
# Step 5b: Plot pre-FFT beam maps (amp + phase)
#----------------------------------------------
echo ""
echo "Step 5b: Plotting pre-FFT beam maps..."
if [ "$NO_PLOT" == "0" ]; then
    python plot_beam_maps.py \
        --prefix "$PREFIX" \
        --output "results/${PREFIX}_beam_maps.pdf"
else
    python plot_beam_maps.py \
        --prefix "$PREFIX" \
        --output "results/${PREFIX}_beam_maps.pdf" \
        --no-interactive
fi

#----------------------------------------------
# Step 6: FFT and aberration fitting
#----------------------------------------------
echo ""
echo "Step 6: FFT and aberration fitting..."
python holis_aber2.py --prefix "$PREFIX"

#----------------------------------------------
# Step 7: Phase unwrapping
#----------------------------------------------
echo ""
echo "Step 7: Phase unwrapping..."
if [ "$NO_PLOT" == "0" ]; then
    python unwrap2d.py -d "$SIZE" --plot
else
    python unwrap2d.py -d "$SIZE"
fi

#----------------------------------------------
# Step 8: Final FFT with unwrapped phase
#----------------------------------------------
echo ""
echo "Step 8: Final FFT with unwrapped phase..."
python holis_aber2.py --unwrap --prefix "$PREFIX"

#----------------------------------------------
# Step 9: Copy results with prefix
#----------------------------------------------
echo ""
echo "Step 9: Saving results..."

# Create results directory if it doesn't exist
mkdir -p results

# Copy output files to results directory
cp ampout.dat results/$PREFIX.ampout
cp phaseout.dat results/$PREFIX.phaseout
cp rgin.dat results/$PREFIX.rgrd
cp Epr.dat results/${PREFIX}_Epr.dat
cp holis.log results/$PREFIX.log
cp holis_fit.json results/${PREFIX}_fit.json 2>/dev/null || true
cp Ep_um.dat results/$PREFIX.Ep_um.dat 2>/dev/null || true
cp Ea_um.dat results/$PREFIX.Ea_um.dat 2>/dev/null || true

# If in compare mode, also copy to _cal suffix
if [ "$DO_COMPARE" == "1" ]; then
    cp Epr.dat results/${PREFIX}_cal_Epr.dat
fi

echo ""
echo "WITH calibration results saved to results/ directory"

#----------------------------------------------
# Step 10: Comparison mode - run without calibration
#----------------------------------------------
if [ "$DO_COMPARE" == "1" ]; then
    echo ""
    echo "=============================================="
    echo "Comparison: Running WITHOUT boresight calibration"
    echo "=============================================="
    
    # Save calibrated results
    mkdir -p compare_backup
    cp rgin.dat compare_backup/rgin_cal.dat
    cp Epr.dat compare_backup/Epr_cal.dat
    
    # Run boresight_cal with --no-correction
    echo ""
    echo "Removing boresight points (no phase correction)..."
    python boresight_cal.py trimmed_with_time.txt -o calibrated_nocal.txt --verbose --no-correction \
        --min-amp 0.35 --slew-margin 2.0 --save-boresight-data --save-prefix "${PREFIX}_nocal"
    
    awk '!/^#/ {print $2, $3, $4, $5}' calibrated_nocal.txt > trimmed.txt
    
    echo "Regridding..."
    if [ "$NO_PLOT" == "1" ]; then
        python regrid_holo.py trimmed.txt regrid.prm --no-plot
    else
        python regrid_holo.py trimmed.txt regrid.prm
    fi
    
    echo "Fixing missing cells..."
    python fix_missing_cell.py "$SIZE" --force
    
    echo "Preprocessing..."
    python preprocess.py
    
    echo "FFT and aberration fitting..."
    python holis_aber2.py --prefix "${PREFIX}_nocal"
    
    echo "Phase unwrapping..."
    python unwrap2d.py -d "$SIZE"
    
    echo "Final FFT with unwrapped phase..."
    python holis_aber2.py --unwrap --prefix "${PREFIX}_nocal"
    
    # Save uncalibrated results
    cp Epr.dat results/${PREFIX}_nocal_Epr.dat
    cp holis.log results/${PREFIX}_nocal.log
    cp holis_fit.json results/${PREFIX}_nocal_fit.json 2>/dev/null || true
    cp Ep_um.dat results/${PREFIX}_nocal.Ep_um.dat 2>/dev/null || true
    cp Ea_um.dat results/${PREFIX}_nocal.Ea_um.dat 2>/dev/null || true
    
    echo ""
    echo "WITHOUT calibration results saved to results/ directory"
    
    #----------------------------------------------
    # Generate comparison plot
    #----------------------------------------------
    echo ""
    echo "Generating comparison plot..."
    
    python << EOF
import numpy as np
import matplotlib.pyplot as plt
import os

prefix = "${PREFIX}"
cal_file = f"results/{prefix}_cal_Epr.dat"
nocal_file = f"results/{prefix}_nocal_Epr.dat"

if not os.path.exists(cal_file) or not os.path.exists(nocal_file):
    print(f"Error: Could not find result files")
    print(f"  Looking for: {cal_file}")
    print(f"  Looking for: {nocal_file}")
    exit(1)

cal_data = np.loadtxt(cal_file)
nocal_data = np.loadtxt(nocal_file)

# Determine grid size
grid_size = int(np.sqrt(len(cal_data)))

# Reshape to grids
cal_grid = cal_data.reshape((grid_size, grid_size))
nocal_grid = nocal_data.reshape((grid_size, grid_size))

# Mask invalid values
cal_grid = np.ma.masked_where(cal_grid < -9000, cal_grid)
nocal_grid = np.ma.masked_where(nocal_grid < -9000, nocal_grid)

# Convert to microns (data is in radians)
error_scaling = 299792458.0 / (94.5e9) * 1e6 / (4.0 * np.pi)
cal_grid = cal_grid * error_scaling
nocal_grid = nocal_grid * error_scaling

# Calculate RMS for valid pixels
cal_rms = np.std(cal_grid.compressed())
nocal_rms = np.std(nocal_grid.compressed())

# Calculate difference
diff_grid = cal_grid - nocal_grid

# Find common color scale
vmax = max(np.abs(cal_grid).max(), np.abs(nocal_grid).max())
vmin = -vmax

# Create figure with 3 subplots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# WITH calibration
im0 = axes[0].imshow(cal_grid, cmap='coolwarm', vmin=vmin, vmax=vmax, origin='upper')
axes[0].set_title(f'WITH boresight cal\nRMS: {cal_rms:.1f} µm')
axes[0].set_xlabel('X pixel')
axes[0].set_ylabel('Y pixel')
plt.colorbar(im0, ax=axes[0], label='Surface error (µm)')

# WITHOUT calibration
im1 = axes[1].imshow(nocal_grid, cmap='coolwarm', vmin=vmin, vmax=vmax, origin='upper')
axes[1].set_title(f'WITHOUT boresight cal\nRMS: {nocal_rms:.1f} µm')
axes[1].set_xlabel('X pixel')
axes[1].set_ylabel('Y pixel')
plt.colorbar(im1, ax=axes[1], label='Surface error (µm)')

# Difference
diff_max = np.abs(diff_grid).max()
im2 = axes[2].imshow(diff_grid, cmap='coolwarm', vmin=-diff_max, vmax=diff_max, origin='upper')
axes[2].set_title(f'Difference (CAL - NOCAL)\nMax diff: {diff_max:.1f} µm')
axes[2].set_xlabel('X pixel')
axes[2].set_ylabel('Y pixel')
plt.colorbar(im2, ax=axes[2], label='Difference (µm)')

plt.tight_layout()
plt.savefig(f"results/{prefix}_boresight_comparison.png", dpi=150)
print(f"Saved comparison plot: results/{prefix}_boresight_comparison.png")
plt.show()
EOF

fi

#----------------------------------------------
# Step 10: Generate surface map visualizations
#----------------------------------------------
echo ""
echo "Step 10: Generating surface maps..."

# Illumination map
if [ -f "results/$PREFIX.Ea_um.dat" ]; then
    echo "Generating illumination map..."
    if [ "$NO_PLOT" == "0" ]; then
        # Show interactively first
        python glt_dish_map.py results/$PREFIX.Ea_um.dat --holo-smooth bilinear --holo-filterrad 0.8 --auto-shift
    fi
    # Then save to PDF
    python glt_dish_map.py results/$PREFIX.Ea_um.dat \
           --holo-smooth bilinear --holo-filterrad 0.8\
            --auto-shift --output results/${PREFIX}_illumination.pdf
fi

# Phase residual map
echo "Generating phase residual map..."
if [ "$NO_PLOT" == "0" ]; then
    # Show interactively first
    python glt_dish_map.py results/${PREFIX}_Epr.dat --vmin -180 --vmax 180 --auto-shift \
        --holo-smooth bilinear --holo-filterrad 0.8\
        --mask-file "mask${SIZE}.dat" --prm-file withphase_aber.prm
fi
# Then save to PDF
python glt_dish_map.py results/${PREFIX}_Epr.dat --vmin -180 --vmax 180 --auto-shift \
     --holo-smooth bilinear --holo-filterrad 0.8\
    --mask-file "mask${SIZE}.dat" --prm-file withphase_aber.prm \
    --output results/${PREFIX}.pdf

#----------------------------------------------
# Step 11: Generate summary PDF
#----------------------------------------------
echo ""
if [ "$NO_PLOT" == "0" ]; then
    read -p "Any comments to add on the summary page? (press Enter to skip): " user_comment
    
    echo ""
    echo "Generating summary PDF..."
    if [ -z "$user_comment" ]; then
        python create_summary_pdf.py "$INPUTFILE" "results/${PREFIX}"
    else
        python create_summary_pdf.py "$INPUTFILE" "results/${PREFIX}" --comment "$user_comment"
    fi
    echo ""
    echo "Summary PDF saved to results/${PREFIX}_summary.pdf"
else
    echo "Generating summary PDF (no comments in --no-plot mode)..."
    python create_summary_pdf.py "$INPUTFILE" "results/${PREFIX}"
    echo "Summary PDF saved to results/${PREFIX}_summary.pdf"
fi

echo ""
echo "=============================================="
echo "REDUCTION COMPLETE"
echo "=============================================="
echo ""
echo "Results saved to results/ directory:"
echo "  Phase residual map: results/${PREFIX}_Epr.dat"
echo "  Surface map PDF: results/${PREFIX}.pdf"
echo "  Illumination map: results/${PREFIX}_illumination.pdf"
if [ "$DO_COMPARE" == "1" ]; then
    echo ""
    echo "Comparison results:"
    echo "  WITH cal: results/${PREFIX}_cal_Epr.dat"
    echo "  WITHOUT cal: results/${PREFIX}_nocal_Epr.dat"
    echo "  Comparison plot: results/${PREFIX}_boresight_comparison.png"
fi

# Reset ido_bore=0 after mode 3 to avoid affecting subsequent runs
if [ "$DO_BORESIGHT_CAL" == "3" ]; then
    echo ""
    echo "Resetting preprocess.prm: ido_bore=0"
    sed -i "s/^Do bore-sight drift correction.*/Do bore-sight drift correction (1\/0)....         0/" preprocess.prm
fi

echo ""

