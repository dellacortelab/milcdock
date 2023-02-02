if [ $# -eq 0 ]; then
    echo "Script to run Vina"
    echo "[-l PATH]: Path to all ligands (mol2 or pdbqt format)"
    echo "[-r PATH]: Path to all receptors"
    echo "[-o PATH]: Path to out directory"
    echo "[--blind]: Perform blind docking (automatically generate box around entire receptor)"
    echo "If the structure files are not in PDBQT format, the path to AutodockTools Utilities24 must be provided to convert them."
    echo "[--center X Y Z]: Coordinates of the center of the binding site (only needed if not performing blind docking. Assumes same coordinates for all receptors)"
    echo "[--size X Y Z]: Size of the binding site (only needed if not performing blind docking)" 
    echo "[-e EXHAUSTIVENESS]: set exhaustiveness parameter (default=8)"
    exit
fi

#read in command line arguments
exhaustiveness=8
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
        x_coord=${all_args[$(($idx+1))]}
        y_coord=${all_args[$(($idx+2))]}
        z_coord=${all_args[$(($idx+3))]}
        center_given=true
    elif [ ${all_args[${idx}]} == "--size" ]; then
        x_size=${all_args[$(($idx+1))]}
        y_size=${all_args[$(($idx+2))]}
        z_size=${all_args[$(($idx+3))]}
        size_given=true
    elif [ ${all_args[${idx}]} == "-e" ]; then
        exhaustiveness=${all_args[$(($idx+1))]}
    elif [ ${all_args[${idx}]} == "--blind" ]; then
        blind=true
    fi
done

script_path=$(dirname "$0")
script_path=$(realpath ${script_path})
autodock_tools_path="/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24"

#convert all lig files to pdbqt if they are not already
if [ -d "$lig_path" ]; then
    cd $lig_path # must be in lig directory for prepare_ligand4.py to work
    for lig in *mol2; do
        base_name=$(basename ${lig} .mol2)
        pdbqt_name=${base_name}.pdbqt
        if ! [ -f $pdbqt_name ]; then
            echo "Converting ${lig} to pdbqt"
            python ${script_path}/fix_atom_names.py ${lig}
            updated_lig_name="${base_name}_renamed.mol2"
            pythonsh ${autodock_tools_path}/prepare_ligand4.py -l ${updated_lig_name} -o ${pdbqt_name}
            rm ${updated_lig_name}
        fi
    done
    cd -
    all_ligs=(${lig_path}/*pdbqt)
else
    extension="${lig_path##*.}"
    dir=$(dirname ${lig_path})
    base_name=$(basename ${lig_path} .mol2)
    pdbqt_name=${base_name}.pdbqt
    if [ $extension != "pdbqt" ] && ! [ -f ${pdbqt_name} ]; then
        echo "Converting ${lig_path} to pdbqt"
        python ${script_path}/fix_atom_names.py ${lig_path}
        updated_lig_name="${base_name}_renamed.mol2"
        cd $dir
        pythonsh ${autodock_tools_path}/prepare_ligand4.py -l ${updated_lig_name} -o ${pdbqt_name}
        rm ${updated_lig_name}
        cd -
    fi
    lig_path=${lig_path%.*}.pdbqt
    all_ligs=(${lig_path})
fi

# Use obabel to add atom type for HN hydrogens because prepare_receptor4.py cannot identify atom type of HN hydrogens
obabel ${rec_path} -O ${rec_path} ---errorlevel 0

#convert receptors to pdbqt if they are not already
if [ -d "$rec_path" ]; then
    for rec in ${rec_path}/*pdb; do
        base_name=${rec%.*}
        pdbqt_name=${base_name}.pdbqt
        if ! [ -f $pdbqt_name ]; then
            echo "Converting ${rec} to pdbqt"
            fixed_rec=${base_name}_amber.pdb
            #convert HIS residues to AMBER convention
            cp $rec $fixed_rec
            sed -i 's/HSD/HID/g' $fixed_rec
            sed -i 's/HSE/HIE/g' $fixed_rec
            sed -i 's/HSP/HIP/g' $fixed_rec
            pythonsh ${autodock_tools_path}/prepare_receptor4.py -r ${fixed_rec} -o ${pdbqt_name} -A hydrogens
            rm $fixed_rec
        fi
    done
    all_recs=(${rec_path}/*pdbqt)
else
    extension="${rec_path##*.}"
    base_name=${rec_path%.*}
    pdbqt_name=${base_name}.pdbqt
    if [ $extension != "pdbqt" ] && ! [ -f ${pdbqt_name} ]; then
        echo "Converting ${rec_path} to pdbqt"
        fixed_rec=${base_name}_amber.pdb
        cp $rec_path $fixed_rec
        #convert HIS residues to AMBER convention
        sed -i 's/HSD/HID/g' $fixed_rec
        sed -i 's/HSE/HIE/g' $fixed_rec
        sed -i 's/HSP/HIP/g' $fixed_rec
        pythonsh ${autodock_tools_path}/prepare_receptor4.py -r ${fixed_rec} -o ${pdbqt_name} -A hydrogens
        rm $fixed_rec
    fi
    rec_path=${rec_path%.*}.pdbqt
    all_recs=(${rec_path})
fi

for rec in "${all_recs[@]}"; do
    rec_name=$(basename $rec .pdbqt)
    config_name="${rec_name}_config"
    rec_out_path="${out_path}/${rec_name}"
    if ! [ -d ${rec_out_path} ]; then
        mkdir $rec_out_path
    fi
    full_out_path="${rec_out_path}/vina"
    if ! [ -d ${full_out_path} ]; then
        mkdir $full_out_path
    fi
    cp ${script_path}/../data/vinaConfig ${full_out_path}/${config_name}

    if [ "${center_given}" = false ] || [ "${size_given}" = false ]; then
        #get coordinates from python script
        echo "Getting coordinates for: $rec"
        raw_coords=$(python ${script_path}/getDimensions.py -r ${rec})
        raw_coords=($raw_coords)
        x_coord=${raw_coords[1]}
        y_coord=${raw_coords[2]}
        z_coord=${raw_coords[3]}
        x_size=${raw_coords[5]}
        y_size=${raw_coords[6]}
        z_size=${raw_coords[7]}
    fi

    sed -i "s/XCOORD/${x_coord}/" ${full_out_path}/${config_name}
    sed -i "s/YCOORD/${y_coord}/" ${full_out_path}/${config_name}
    sed -i "s/ZCOORD/${z_coord}/" ${full_out_path}/${config_name}
    sed -i "s/XSIZE/${x_size}/" ${full_out_path}/${config_name}
    sed -i "s/YSIZE/${y_size}/" ${full_out_path}/${config_name}
    sed -i "s/ZSIZE/${z_size}/" ${full_out_path}/${config_name}
    sed -i "s/EXHAUST/${exhaustiveness}/" ${full_out_path}/${config_name}

    for lig in ${all_ligs[@]}; do
        lig_name=$(basename $lig .pdbqt)
        cmd="vina --config ${full_out_path}/${config_name} --receptor ${rec} --ligand ${lig} --out ${full_out_path}/${lig_name}.pdbqt"
        echo -e "Running command:\n${cmd}"
        ${cmd}
        echo "Done"
    done
done 
echo "All Vina runs complete."
