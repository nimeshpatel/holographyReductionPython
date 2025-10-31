#!/bin/bash
# Pipeline script for holography maps reduction
# Nimesh Patel
# 29 October 2025
#
# To be run on gltobscon, under holo/holoReducePy area.
# If you want to run this on another machine, copy over 
# all the .py codes, .prm and diff512.dat files (and this script).
# 

set -e  # Exit immediately if any command fails
set -x  # Print commands as they execute

# Create results directory if it doesn't exist
mkdir -p results
echo "Results will be saved to the 'results' subdirectory."

echo "Starting with the first line detection of start of the map,"
echo "and writing out trimmed.txt file, deleting all the earlier lines."

python detect_raster_start.py $1 -o trimmed.txt

echo "Regridding data to a 128 x 128 map..."
python regrid_holo.py trimmed.txt regrid_128x128_29oct25.prm

echo "Preprocessing..."
python preprocess.py

echo "Fourier transform...running holis..."
python holis_aber2.py

echo "Unwrapping phase..."
python unwrap.py

echo "Fourier transform...repeating holis..."
python holis_aber2.py --unwrap

echo "Copying relevant output files to results directory..."
cp ampout.dat results/$2.ampout
cp phaseout.dat results/$2.phaseout
cp rgin.dat results/$2.rgrd
cp Epr.dat results/$2_Epr.dat
cp holis.log results/$2.log
cp Ep_um.dat results/$2.Ep_um.dat
cp Ea_um.dat results/$2.Ea_um.dat

echo "Generating illumination map..."
python glt_dish_map.py results/$2.Ea_um.dat --x-shift -65 --y-shift 65 
python glt_dish_map.py results/$2.Ea_um.dat --x-shift -65 --y-shift 65 --output results/$2_illumination.pdf

echo "Generating phase residual map..."
python glt_dish_map.py results/$2_Epr.dat --vmin -180 --vmax 180 --x-shift -65 --y-shift 65
python glt_dish_map.py results/$2_Epr.dat --vmin -180 --vmax 180 --x-shift -65 --y-shift 65 --output results/$2.pdf

echo "Holography data reduction completed."
echo "Phase residual map saved to results/$2.pdf"
echo "Illumination map saved to results/$2_illumination.pdf"

# Prompt for user comments
echo ""
read -p "Any comments to add on the summary page? (press Enter to skip): " user_comment

# Generate summary PDF
echo ""
echo "Generating summary PDF..."
if [ -z "$user_comment" ]; then
    python create_summary_pdf.py "$1" "results/$2"
else
    python create_summary_pdf.py "$1" "results/$2" --comment "$user_comment"
fi

echo ""
echo "Summary PDF saved to results/$2_summary.pdf"
