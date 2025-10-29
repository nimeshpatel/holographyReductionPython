set -x
python detect_raster_start_interactive_v3.py $1 -o trimmed.txt
python regrid_holo.py trimmed.txt regrid_128x128_15sep25.prm
python preprocess.py
python holis_aber2.py 
python unwrap.py 
python holis_aber2.py --unwrap
python gltDishMap.py Epr.dat --cmap coolwarm --vmin -180 --vmax 180
