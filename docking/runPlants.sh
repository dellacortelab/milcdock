#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Script to run PLANTS"
    echo "[-l PATH]: Path to all ligands in mol2 format."
    echo "[-r PATH]: Path to all receptors (pdb or mol2)"
    echo "[-o PATH]: Path to out directory"
    echo "[--center X Y Z]: Coordinates of the center of the binding site"
    echo "[--size R]: Radius of the binding site" 
    exit
fi

center_given=false
size_given=false

IFS=' ' read -r -a all_args <<< "$@"

for ((idx=0; idx<$#; idx++)); do
    if [ ${all_args[${idx}]} == "-l" ]; then
        lig_path=${all_args[$(($idx+1))]} 
    elif [ ${all_args[${idx}]} == "-r" ]; then
        rec_path=${all_args[$(($idx+1))]}
    elif [ ${all_args[${idx}]} == "-o" ]; then
        out_path=${all_args[$(($idx+1))]}
        if ! [ -d "${out_path}" ]; then
            mkdir ${out_path}
        fi
    elif [ ${all_args[${idx}]} == "--center" ]; then
        x_center=${all_args[$(($idx+1))]}
        y_center=${all_args[$(($idx+2))]}
        z_center=${all_args[$(($idx+3))]}
        center_given=true
    elif [ ${all_args[${idx}]} == "--size" ]; then
        radius=${all_args[$(($idx+1))]}
        size_given=true
    fi
done 

script_path=$(dirname "$0")
script_path=$(realpath ${script_path})
plants_config="${script_path}/../data/plantsconfig"

lig_path=$(realpath $lig_path)
rec_path=$(realpath $rec_path)

#convert receptor to mol2 if not already
if [ -d "$rec_path" ]; then
    for rec in ${rec_path}/*pdb; do
        base_name=${rec%.*}
        mol2_name=${base_name}.mol2
        if ! [ -f $mol2_name ]; then
            echo "Converting ${rec} to mol2"
            obabel -ipdb $rec -omol2 -O$mol2_name -h
        fi
    done
    all_recs=(${rec_path}/*mol2)
else
    extension="${rec_path##*.}"
    base_name=${rec_path%.*}
    mol2_name=${base_name}.mol2
    if [ $extension != "mol2" ] && ! [ -f $mol2_name ]; then
        echo "Converting ${rec_path} to mol2"
        obabel -ipdb $rec_path -omol2 -O$mol2_name -h
    fi  
    rec_path=${rec_path%.*}.mol2
    all_recs=(${rec_path})
fi

#convert all lig files to mol2 using SPORES if they are not already in mol2 format, and create ligand list in a .txt file
if [ -d "$lig_path" ]; then
    all_ligs=(${lig_path}/*mol2)
else
    all_ligs=(${lig_path})
fi

for rec in ${all_recs[@]}; do
    rec_basename=$(basename ${rec} .mol2)
    rec_name=$(basename ${rec})
    rec_config="${rec_basename}_plantsconfig"
    rec_out_path="${out_path}/${rec_basename}"
    if ! [ -d ${rec_out_path} ]; then
        mkdir ${rec_out_path}
    fi
    full_out_path="${rec_out_path}/plants"
    if ! [ -d ${full_out_path} ]; then
        mkdir ${full_out_path}
    fi
    cd $full_out_path
    lig_list="./ligands.txt"
    if ! [ -f $lig_list ]; then
        for lig in ${all_ligs[@]}; do
            echo $lig >> $lig_list
        done
    fi
    cp $plants_config ./$rec_config
    cp $rec .

    #get binding site coordinates for blind docking if none are supplied (Blind docking not yet implemented, but included here to avoid errors)
    if [ "$center_given" = false ] || [ "$size_given" = false ]; then
        echo "Getting coordinates for blind docking on ${rec_name}"

        raw_coords=$(python ${script_path}/getDimensions.py -r ${rec_name} --plants)
        raw_coords=($raw_coords)
        x_center=${raw_coords[1]}
        y_center=${raw_coords[2]}
        z_center=${raw_coords[3]}
        radius=${raw_coords[5]}
    fi 

    #Modify config file
    sed -i "s/RECEPTOR/${rec_name}/" $rec_config
    sed -i "s/LIGANDS/ligands.txt/" $rec_config
    sed -i "s/OUTPUT/./" $rec_config
    sed -i "s/XCOORD/${x_center}/" $rec_config
    sed -i "s/YCOORD/${y_center}/" $rec_config
    sed -i "s/ZCOORD/${z_center}/" $rec_config
    sed -i "s/RADIUS/${radius}/" $rec_config

    #Perform docking with PLANTS
    echo "Starting docking calculations on ${rec_basename}"
    echo $rec_config
    exit
    plants --mode screen $rec_config
    cd -
done
echo "All docking runs complete"


