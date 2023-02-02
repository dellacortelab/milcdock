if [ $# -eq 0 ]; then
    echo -e "\nScript to run Ledock"
    echo "[-l PATH]: Path to all ligands (mol2 format)"
    echo "[-r PATH]: Path to all receptors"
    echo "[-o PATH]: Path to out directory"
    echo "[--blind]: Perform blind docking (automatically generate box around entire receptor, default if no coordinates given)"
    echo "[--center X Y Z]: Coordinates of the center of the binding site (only needed if not performing blind docking. Assumes same coordinates for all receptors)"
    echo -e "[--size X Y Z]: Size of the binding site (only needed if not performing blind docking)\n" 
    exit
fi


center_given=false
size_given=false
blind=false

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
        x_size=${all_args[$(($idx+1))]}
        y_size=${all_args[$(($idx+2))]}
        z_size=${all_args[$(($idx+3))]}
        size_given=true
    elif [ ${all_args[${idx}]} == "--blind" ]; then
        blind=true
    fi
done

script_path=$(dirname $0)
script_path=$(realpath ${script_path})

if [ "${center_given}" = true ] && [ "${size_given}" = true ]; then
    xvals=($(python ${script_path}/centerSizeToMinMax.py ${x_center} ${x_size}))
    yvals=($(python ${script_path}/centerSizeToMinMax.py ${y_center} ${y_size}))
    zvals=($(python ${script_path}/centerSizeToMinMax.py ${z_center} ${z_size}))
    xmin=${xvals[0]}
    xmax=${xvals[1]}
    ymin=${yvals[0]}
    ymax=${yvals[1]}
    zmin=${zvals[0]}
    zmax=${zvals[1]}
fi

if [ -d "$lig_path" ]; then
    all_lig=(${lig_path}/*mol2)
else
    all_lig=(${lig_path})
fi

#create array of receptor file paths
if [ -d ${rec_path} ]; then
    all_rec=(${rec_path}/*pdb)
else
    all_rec=(${rec_path})
    if [ $(basename ${rec_path} .pdb) == $(basename $rec_path) ]; then
        echo "For Ledock, receptor must be in PDB format. Cancelling calculation."
        exit
    fi
fi

current_dir=$(pwd)


for rec in "${all_rec[@]}"; do
    rec_name=$(basename $rec .pdb)
    if ! [ -d "${out_path}/${rec_name}" ]; then
        mkdir ${out_path}/${rec_name} 
    fi
    rec_out="${out_path}/${rec_name}/ledock"
    if ! [ -d "${rec_out}" ]; then
        mkdir ${rec_out}
    fi
    cp ${rec} ${rec_out}
    for lig in ${all_lig[@]}; do
        cp ${lig} $rec_out
    done 
    cd ${rec_out}
    ls *mol2 > ligands

    # remove H from receptor, then change HIS residue names to match CHARMM convention
    awk '$3 !~ /^H/ { print $0 }' ${rec_name}.pdb > receptor_noh.pdb
    sed -i 's/HID/HSD/g' receptor_noh.pdb
    sed -i 's/HIE/HSE/g' receptor_noh.pdb
    sed -i 's/HIP/HSP/g' receptor_noh.pdb
    sed -i 's/ZN5/ ZN/g' receptor_noh.pdb
    mv receptor_noh.pdb ${rec_name}.pdb

    echo "Removed H from $rec_name and changed HIS residue names to match CHARMM format."

    lepro ${rec_name}.pdb

    if [ "${center_given}" = false ] || [ "${size_given}" = false ]; then
        #get coordinates from python script
        echo "Getting coordinates for: $rec_name"
        raw_coords=$(python ${script_path}/getDimensions.py -r ${rec_name}.pdb --ledock)
        raw_coords=($raw_coords)
        xmin=${raw_coords[1]}
        ymin=${raw_coords[2]}
        zmin=${raw_coords[3]}
        xmax=${raw_coords[5]}
        ymax=${raw_coords[6]}
        zmax=${raw_coords[7]}
    fi
   
    #convert center/size coords to min/max 
    sed -i "s/xmin/${xmin}/" dock.in
    sed -i "s/xmax/${xmax}/" dock.in
    sed -i "s/ymin/${ymin}/" dock.in
    sed -i "s/ymax/${ymax}/" dock.in
    sed -i "s/zmin/${zmin}/" dock.in
    sed -i "s/zmax/${zmax}/" dock.in

    echo "Running Ledock on ${rec_name}" 
    ledock dock.in
    
    cd $current_dir

done
