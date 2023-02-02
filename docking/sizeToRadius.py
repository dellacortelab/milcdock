import sys
import numpy as np

if len(sys.argv) == 1:
    print("Script to convert size of box coordinates to a radius.")
    print("Args: xsize, ysize, zsize are first 3 arguments")
    print("Use --plants flag after size values to run Plants-specific radius conversion.")
    print("Output: radius")
    sys.exit()

do_plants = False
size = np.array(sys.argv[1:4],dtype=float)
for i,arg in enumerate(sys.argv):
    if arg == '--plants':
        do_plants = True

if do_plants:
    radius = np.linalg.norm(size/2)/2
else:
    radius = np.linalg.norm(size/2)/3

print(f'{radius:3.3f}')

