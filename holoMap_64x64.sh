#!/bin/bash
# Pipeline script for holography maps reduction (64x64 version)
# Nimesh Patel
# 16 December 2025
#
# Usage: ./holoMap_64x64.sh <input_file.txt>
# The output file prefix is automatically derived from the input filename.
#
# This version produces a 64x64 map for better comparison with photogrammetry.
# To be run on gltobscon, under holo/holoReducePy area.
# If you want to run this on another machine, copy over
# all the .py codes, .prm and diff512.dat files (and this script).
#

set -e  # Exit immediately if any command fails
set -x  # Print commands as they execute

# Activate the nimesh_holo environment (has all the needed packages)
#source ~/locutus/mambaforge/bin/activate_mamba nimesh_holo
source $HOME/.mamba_rc
mamba activate nimesh_holo



# Extract prefix from input filename (remove .txt extension)
INPUTFILE="$1"
PREFIX=$(basename "$INPUTFILE" .txt)

# Create results directory if it doesn't exist
mkdir -p results
echo "Results will be saved to the 'results' subdirectory."

# Swap parameter files for 64x64 processing
echo "Setting up 64x64 parameter files..."
cp preprocess.prm preprocess.prm.bak
cp withphase_aber.prm withphase_aber.prm.bak
cp preprocess_64x64.prm preprocess.prm
cp withphase_aber_64x64.prm withphase_aber.prm

# Function to restore original parameter files on exit
cleanup() {
    echo "Restoring original parameter files..."
    mv preprocess.prm.bak preprocess.prm
    mv withphase_aber.prm.bak withphase_aber.prm
}
trap cleanup EXIT

echo "Starting with the first line detection of start of the map,"
echo "and writing out trimmed.txt file, deleting all the earlier lines."

python detect_raster_start.py $1 -o trimmed.txt

echo "Regridding data to a 64 x 64 map..."
python regrid_holo.py trimmed.txt regrid_64x64.prm

echo "Preprocessing..."
python preprocess.py

echo "Fourier transform...running holis..."
python holis_aber2.py

echo "Unwrapping phase..."
python unwrap.py

echo "Fourier transform...repeating holis..."
python holis_aber2.py --unwrap

echo "Copying relevant output files to results directory..."
cp ampout.dat results/$PREFIX.ampout
cp phaseout.dat results/$PREFIX.phaseout
cp rgin.dat results/$PREFIX.rgrd
cp Epr.dat results/${PREFIX}_Epr.dat
cp holis.log results/$PREFIX.log
cp Ep_um.dat results/$PREFIX.Ep_um.dat
cp Ea_um.dat results/$PREFIX.Ea_um.dat

echo "Generating illumination map..."
python glt_dish_map.py results/$PREFIX.Ea_um.dat --x-shift -65 --y-shift 65
python glt_dish_map.py results/$PREFIX.Ea_um.dat --x-shift -65 --y-shift 65 --output results/${PREFIX}_illumination.pdf
python glt_dish_map.py results/$PREFIX.Ea_um.dat --x-shift -65 --y-shift 65 --output results/${PREFIX}_illumination.png

echo "Generating phase residual map..."
python glt_dish_map.py results/${PREFIX}_Epr.dat --vmin -180 --vmax 180 --x-shift -65 --y-shift 65
python glt_dish_map.py results/${PREFIX}_Epr.dat --vmin -180 --vmax 180 --x-shift -65 --y-shift 65 --output results/$PREFIX.pdf
python glt_dish_map.py results/${PREFIX}_Epr.dat --vmin -180 --vmax 180 --x-shift -65 --y-shift 65 --output results/${PREFIX}.png

echo "Holography data reduction completed."
echo "Phase residual map saved to results/$PREFIX.pdf"
echo "Illumination map saved to results/${PREFIX}_illumination.pdf"

# Prompt for user comments
echo ""
read -p "Any comments to add on the summary page? (press Enter to skip): " user_comment

# Generate summary PDF
echo ""
echo "Generating summary PDF..."
if [ -z "$user_comment" ]; then
    python create_summary_pdf.py "$INPUTFILE" "results/$PREFIX"
else
    python create_summary_pdf.py "$INPUTFILE" "results/$PREFIX" --comment "$user_comment"
fi

echo ""
echo "Summary PDF saved to results/${PREFIX}_summary.pdf"
