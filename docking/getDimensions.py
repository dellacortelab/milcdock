#script to get dimensions of a protein for both full blind docking and iterative blind docking. 

import MDAnalysis as mda
import numpy as np
import sys
import warnings

ledock = False
rdock = False
plants = False
out_path_given = False

if len(sys.argv) > 1:
    for i,arg in enumerate(sys.argv):
        if arg == "-r":
            rec_path = sys.argv[i+1]
        elif arg == "--ledock":
            ledock = True
        elif arg == "--rdock":
            rdock = True
        elif arg == "--plants":
            plants = True
        elif arg == "-o":
            out_path = sys.argv[i+1]
            out_path_given = True
else:
    print("\nScript to get the dimensions of a protein for automated blind docking.\n")
    print("[-r PATH]: Path to protein or crystal ligand file.")
    print("[-o PATH]: Path to file to write coordinates to. If not given, will print to standard out.")
    print("[--ledock]: Output for Ledock simulation (in min/max form, rather than center/size")
    print("[--rdock]: Output for rDock simulation (in center, radius form)")
    print("[--plants]: Output for PLANTS simulation (in center, radius form)")
    sys.exit()

#create MDA universe object

with warnings.catch_warnings():
    warnings.simplefilter('ignore',category=UserWarning)
    uni = mda.Universe(rec_path)

all_positions = uni.atoms.positions

#Defines box by widest dimensions of molecule
center_coords = []
xyz_size = []
for i in range(3): 
    dim_max = max(all_positions[:,i])
    dim_min = min(all_positions[:,i])
    center_coords.append((dim_max + dim_min) / 2.)
    xyz_size.append(abs(dim_max - dim_min)+15) #add 15 Angst. to create wider binding region

xyz_size = np.array(xyz_size)
center_coords = np.array(center_coords)
#add length to each dimension if less than 22.5x22.5x22.5 angst.^3
if np.average(xyz_size) < 25: #22.5 Angs.^3 is the Autodock standard I remember.
    size_increase = 25 - np.average(xyz_size)
    xyz_size = xyz_size + size_increase

if ledock:
    str1 = f"Mins: {center_coords[0]-xyz_size[0]/2:3.3f} {center_coords[1]-xyz_size[1]/2:3.3f} {center_coords[2]-xyz_size[2]/2:3.3f}"
    str2 = f"Maxs: {center_coords[0]+xyz_size[0]/2:3.3f} {center_coords[1]+xyz_size[1]/2:3.3f} {center_coords[2]+xyz_size[2]/2:3.3f}"
elif rdock:
    str1 = f"Center: {center_coords[0]:3.3f} {center_coords[1]:3.3f} {center_coords[2]:3.3f}"
    str2 = f"Radius: {np.linalg.norm(xyz_size/2)/3}"
elif plants:
    str1 = f"Center: {center_coords[0]:3.3f} {center_coords[1]:3.3f} {center_coords[2]:3.3f}"
    str2 = f"Radius: {np.linalg.norm(xyz_size/2)/2}"
else:
    str1 = f"Center: {center_coords[0]:3.3f} {center_coords[1]:3.3f} {center_coords[2]:3.3f}"
    str2 = f"Size: {xyz_size[0]:3.3f} {xyz_size[1]:3.3f} {xyz_size[2]:3.3f}"
if not out_path_given:
    print(str1)
    print(str2)
else:
    file_lines = [str1 + "\n", str2 + "\n"]
    with open(out_path, "w") as f:
        print(f"Writing coordinates to {out_path}")
        f.writelines(file_lines)
