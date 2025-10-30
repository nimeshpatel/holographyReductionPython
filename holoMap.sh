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

python detect_raster_start_interactive_v3.py $1 -o trimmed.txt
python regrid_holo.py trimmed.txt regrid_128x128_29oct25.prm
python preprocess.py
python holis_aber2.py
python unwrap.py
python holis_aber2.py --unwrap
cp Epr.dat $2_Epr.dat
python glt_dish_map.py $2_Epr.dat --cmap coolwarm --vmin -180 --vmax 180 --x-shift -65 --y-shift 65
python glt_dish_map.py $2_Epr.dat --cmap coolwarm --vmin -180 --vmax 180 --x-shift -65 --y-shift 65 --output $2.pdf

echo "Holography data reduction complete! Output saved to $2.pdf"
