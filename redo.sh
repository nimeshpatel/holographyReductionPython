cpy3 preprocess.py
cpy3 holis_aber2.py
cpy3 unwrap2d.py -d 128 --plot
cpy3 holis_aber2.py --unwrap
cpy3 glt_dish_map.py Epr.dat --vmin -180 --vmax 180 --x-shift -65 --y-shift 65 --mask-file mask128.dat --prm-file withphase_aber.prm
