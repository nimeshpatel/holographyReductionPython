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

echo "Starting with the first line detection of start of the map, and writing
out trimmed.txt file, deleting all the earlier lines."

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


python glt_dish_map.py $2_Epr.dat --cmap coolwarm --vmin -180 --vmax 180 --x-shift -65 --y-shift 65
python glt_dish_map.py $2_Epr.dat --cmap coolwarm --vmin -180 --vmax 180 --x-shift -65 --y-shift 65 --output $2.pdf

echo "Copying relevant output files to save them for this run.
cp ampout.dat $2.ampout
cp phaseout.dat $2.phaseout
cp rgin.dat $2.rgrd
cp Epr.dat $2_Epr.dat
cp holis.log $2.log
cp Ep_um.dat $2.Ep_um.dat
cp Ea_um.dat $2.Ea_um.dat

echo "Holography data reduction completed. Output saved to $2.pdf."
