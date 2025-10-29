set -x
python detect_raster_start_interactive_v3.py $1 -o trimmed.txt
python regrid_holo.py trimmed.txt regrid_128x128_15sep25.prm
python preprocess.py
python holis_aber2.py 
python unwrap.py 
python holis_aber2.py --unwrap
cp Epr.dat $2_Epr.dat
python gltDishMap.py $2_Epr.dat --cmap coolwarm --vmin -180 --vmax 180 
python gltDishMap.py $2_Epr.dat --cmap coolwarm --vmin -180 --vmax 180 --output $2.pdf
