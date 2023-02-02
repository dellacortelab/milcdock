if [ $# -eq 0 ]; then
    echo "Script to run AutoDock"
    echo "[-l PATH]: Path to all ligands (mol2 or pdbqt format)"
    echo "[-r PATH]: Path to all receptors"
    echo "[-o PATH]: Path to out directory"
    echo "[--center X Y Z]: Coordinates of the center of the binding site (only needed if not performing blind docking. Assumes same coordinates for all receptors)"
    exit
fi

IFS=' ' read -r -a all_args <<< "$@"

for ((idx=0; idx<$#; idx++)); do
    if [ ${all_args[${idx}]} == "-l" ]; then
        lig_path=${all_args[$(($idx+1))]}
    elif [ ${all_args[${idx}]} == "-c" ]; then
        crystal_path=${all_args[$(($idx+1))]}
        crystal_base_name=${crystal_path%.*}
        crystal_pdbqt_name=${crystal_base_name}.pdbqt
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
    fi
done

lig_path=$(realpath ${lig_path})
rec_path=$(realpath ${rec_path})
script_path=$(dirname "$0")
script_path=$(realpath ${script_path})
autodock_tools_path="/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24"

#convert all lig files to pdbqt if they are not already
if [ -d "$lig_path" ]; then
    cd $lig_path
    for lig in *mol2; do
        base_name=${lig%.*}
        pdbqt_name=${base_name}.pdbqt
        if ! [ -f $pdbqt_name ]; then
            echo "Converting ${lig} to pdbqt"
            python ${script_path}/fix_atom_names.py ${lig}
            updated_lig_name="${base_name}_renamed.mol2"
            pythonsh ${autodock_tools_path}/prepare_ligand4.py -l ${updated_lig_name} -o ${pdbqt_name}
            rm $updated_lig_name
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

# Convert receptors to pdbqt if they are not already
if [ -d "$rec_path" ]; then
    for rec in ${rec_path}/*pdb; do
        base_name=${rec%.*}
        pdbqt_name=${base_name}.pdbqt
        if ! [ -f $pdbqt_name ]; then
            echo "Converting ${rec} to pdbqt"
	        fixedRec=${base_name}_amber.pdb
            cp $rec $fixedRec
            #convert HIS residues to AMBER convention
            sed -i 's/HSD/HID/g' $fixedRec
            sed -i 's/HSE/HIE/g' $fixedRec
            sed -i 's/HSP/HIP/g' $fixedRec
	        pythonsh ${autodock_tools_path}/prepare_receptor4.py -r ${fixedRec} -o ${pdbqt_name} -A hydrogens
            rm $fixedRec
        fi
    done
    allRecs=(${rec_path}/*pdbqt)
else
    extension="${rec_path##*.}"
    base_name=${rec_path%.*}
    pdbqt_name=${base_name}.pdbqt
    if [ $extension != "pdbqt" ] && ! [ -f ${pdbqt_name} ]; then
        echo "Converting ${rec_path} to pdbqt"
        fixedRec=${base_name}_amber.pdb
        cp $rec_path $fixedRec
        #convert HIS residues to AMBER convention
        sed -i 's/HSD/HID/g' $fixedRec
        sed -i 's/HSE/HIE/g' $fixedRec
        sed -i 's/HSP/HIP/g' $fixedRec
        pythonsh ${autodock_tools_path}/prepare_receptor4.py -r ${fixedRec} -o ${pdbqt_name} -A hydrogens
        rm $fixedRec
    fi
    rec_path=${rec_path%.*}.pdbqt
    allRecs=(${rec_path})
fi

for rec in ${allRecs[@]}; do
    recName=$(basename $rec .pdbqt)
    recOutPath="${out_path}/${recName}"
    if ! [ -d $recOutPath ]; then
        mkdir $recOutPath
    fi
    recOutPath="${recOutPath}/autodock"
    if ! [ -d $recOutPath ]; then
        mkdir $recOutPath
    fi
    cp $rec $recOutPath
    cd $recOutPath

    if [ "${center_given}" = false ]; then
        #get coordinates from python script
        echo "Getting coordinates for: $rec"
        rawCoords=$(python getDimensions.py -r ${rec})
        rawCoords=($rawCoords)
        x_coord=${rawCoords[1]}
        y_coord=${rawCoords[2]}
        z_coord=${rawCoords[3]}
    fi
    
    # Create gpf file and run autogrid
    pythonsh ${autodock_tools_path}/prepare_gpf4.py -r ${rec} -d ${all_ligs%/*} -p gridcenter=CENTER -o ${recName}.gpf
    
    #To implement blind docking, modify this sed command to change CENTER to auto
    sed -i "s/CENTER/${x_coord} ${y_coord} ${z_coord}/" ${recName}.gpf

    autogrid4 -p ${recName}.gpf -l ${recName}_grid.glg
    
    echo "Performing docking on ${recName}"

    # Run Docking for each ligand
    for lig in ${all_ligs[@]}; do
        ligName=$(basename $lig .pdbqt)
        cp $lig .
        echo "Docking ${ligName}"

        pythonsh ${autodock_tools_path}/prepare_dpf42.py -r ${rec} -l ${lig} -p ga_run=20 -p tran0=CENTER -o ${ligName}_${recName}.dpf

        #To implement blind docking, modify this sed command to change CENTER to random
        sed -i "s/CENTER/${x_coord} ${y_coord} ${z_coord}/" ${ligName}_${recName}.dpf

        autodock4 -p ${ligName}_${recName}.dpf -l ${ligName}_${recName}.dlg
        rm ${ligName}_${recName}.dpf
        rm ./${ligName}.pdbqt
    done
    echo "${recName} docking finished"
    cd -
done
echo "All runs complete"


