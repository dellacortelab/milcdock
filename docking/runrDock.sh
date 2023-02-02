if [ $# -eq 0 ]; then
    echo -e "\nScript to run rDock"
    echo "[-l PATH]: Path to all ligands (pdb, mol2, or sd format. If not sd, requires a working install of Open Babel to convert to mol2).)"
    echo "[-r PATH]: Path to all receptors (pdb or mol2 format. If pdb, requires Open Babel to convert to mol2. If mol2, requires hydrogens and charges)"
    echo "[-o PATH]: Path to out directory"
    echo "[--blind]: Perform blind docking (automatically generate box around entire receptor, default if no coordinates given)"
    echo "[--center X Y Z]: Coordinates of the center of the binding site (only needed if not performing blind docking. Assumes same coordinates for all receptors)"
    echo -e "[--size R]: Radius of the binding site (only needed if not performing blind docking)\n" 
    exit
fi

center_given=false
size_given=false
blind=false

IFS=' ' read -r -a allArgs <<< "$@"

for ((idx=0; idx<$#; idx++)); do
    if [ ${allArgs[${idx}]} == "-l" ]; then
        lig_path=${allArgs[$(($idx+1))]} 
    elif [ ${allArgs[${idx}]} == "-r" ]; then
        rec_path=${allArgs[$(($idx+1))]}
    elif [ ${allArgs[${idx}]} == "-o" ]; then
        out_path=${allArgs[$(($idx+1))]}
        if ! [ -d "${out_path}" ]; then
            mkdir ${out_path}
        fi
    elif [ ${allArgs[${idx}]} == "--center" ]; then
        x_center=${allArgs[$(($idx+1))]}
        y_center=${allArgs[$(($idx+2))]}
        z_center=${allArgs[$(($idx+3))]}
        center_given=true
    elif [ ${allArgs[${idx}]} == "--size" ]; then
        radius=${allArgs[$(($idx+1))]}
        size_given=true
    elif [ ${allArgs[${idx}]} == "--blind" ]; then
        blind=true
    fi
done 

script_dir=$(dirname $0)
script_dir=$(realpath ${script_dir})

#convert receptor to mol2 if not already
if [ -d "$rec_path" ]; then
    for rec in ${rec_path}/*pdb; do
        base_name=${rec%.*}
        mol2_name=${base_name}.mol2
        if ! [ -f $mol2_name ]; then
            echo "Converting ${rec} to mol2"
            obabel -ipdb $rec -omol2 -O$mol2_name -h --partialcharge gasteiger #documentation recommends adding charges, but didn't have them added for training.
        fi
    done
    all_recs=(${rec_path}/*mol2)
else
    extension="${rec_path##*.}"
    base_name=${rec_path%.*}
    mol2_name=${base_name}.mol2
    if [ $extension != "mol2" ] && ! [ -f $mol2_name ]; then
        echo "Converting ${rec_path} to mol2"
        obabel -ipdb $rec_path -omol2 -O$mol2_name -h --partialcharge gasteiger #documentation recommends adding charges, but didn't have them added for training.
    fi  
    rec_path=${rec_path%.*}.mol2
    all_recs=(${rec_path})
fi

lig_path=$(realpath ${lig_path})
if [ -d "$lig_path" ]; then
    for lig in ${lig_path}/*mol2; do
        base_name=${lig%.*}
        sd_name=${base_name}.sd
        if ! [ -f $sd_name ]; then
            echo "Converting ${lig} to sd"
            obabel -imol2 $lig -osd -O $sd_name -h
        fi
    done
    all_ligs=(${lig_path}/*sd)
else
    extension="${lig_path##*.}"
    base_name=${lig_path%.*}
    sd_name=${base_name}.sd
    if [ $extension != "sd" ] && ! [ -f $sd_name ]; then
        echo "Converting ${lig_path} to sd"
        if [ "${extension}" == "mol2" ]; then
            obabel -imol2 $lig_path -osd -O $sd_name -h
        fi
    fi
    lig_path=${lig_path%.*}.sd
    all_ligs=(${lig_path})
fi

for rec in ${all_recs[@]}; do
    rec_basename=$(basename ${rec} .mol2)
    rec_name=$(basename ${rec})
    new_dir=$out_path/${rec_basename}
    if ! [ -d $new_dir ]; then
        mkdir ${new_dir}
    fi
    full_dir="${new_dir}/rdock"
    if ! [ -d $full_dir ]; then
        mkdir ${full_dir}
    fi
    cp ${rec} ${full_dir}/

    cd $full_dir

    #get binding site coordinates for blind docking if none are supplied
    if [ "$center_given" = false ] || [ "$size_given" = false ]; then
        echo "Getting coordinates for blind docking on ${rec_name}"

        raw_coords=$(python ${script_dir}/getDimensions.py -r ${rec_name} --rdock)
        raw_coords=($raw_coords)
        x_center=${raw_coords[1]}
        y_center=${raw_coords[2]}
        z_center=${raw_coords[3]}
        radius=${raw_coords[5]}
    fi 
    
    #prep rbcavity input
    cp ${script_dir}/../data/rdock_prep.prm .
    sed -i "s/xcoord/${x_center}/" rdock_prep.prm
    sed -i "s/ycoord/${y_center}/" rdock_prep.prm
    sed -i "s/zcoord/${z_center}/" rdock_prep.prm
    sed -i "s/binding_radius/${radius}/" rdock_prep.prm
    sed -i "s/receptor_path/${rec_name}/" rdock_prep.prm
    if [ "$blind" = true ]; then
       sed -i "s/cavities/99/" rdock_prep.prm
    else
       sed -i "s/cavities/1/" rdock_prep.prm
    fi 

    #calculate docking cavities
    rbcavity -W -d -r rdock_prep.prm

    #loop through all ligands, run docking calculations 
    cp ${script_dir}/../data/rdock_dock.prm .
    for lig in ${all_ligs[@]}; do
        cp $lig .
        cp ${lig%.*}.mol2 .
        lig_name=$(basename $lig)
        lig_base=$(basename $lig .sd)
        echo "Starting docking calculations of $lig_base on $rec_basename"
        rbdock -i ${lig_name} -o ${lig_base}_out -r rdock_prep.prm -p dock.prm -n 20 -t 8 -C
    done
    cd -
done
