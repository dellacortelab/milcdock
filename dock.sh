#Script to run Autodock Vina, Ledock, Dock6, and rDock for consensus and ensemble docking.

function print_help {
    echo -e "\nScript to run Autodock Vina, Ledock, rDock, Autodock4, and Plants for consensus and ensemble docking.\n"
    echo "[-r PATH]: Path to hydrogenated receptor files in PDB format. Can be a path to a single file or a path to a directory containing multiple structures."
    echo "[-l PATH]: Path to hydrogenated ligand files in mol2 format. Can be a path to a single file or a path to a directory containing multiple files."
    echo "[-o PATH]: Path to out directory."
    echo "[-c PATH]: Path to crystal ligand structure to define binding site (negates blind docking)."
    echo "--vina run docking with Vina."
    echo "Vina parameters:"
    echo "      [-e EXHAUSTIVENESS]: set exhaustiveness parameter (default=8)."
    echo "--ledock run docking with Ledock."
    echo "--rdock: run docking with rDock."
    echo "--autodock: run docking with autodock."
    echo "--plants: run docking with PLANTS"
    echo "--blind: Run blind docking (automatic box generation)."
    echo "To give custom binding site coordinates (negates blind docking and coordinate generation from ligand structure):"
    echo "  [--center X Y Z]: set center of binding site"
    echo "  [--size X Y Z]: set size of binding site."
    echo -e "-h: print this help text.\n"
    echo "Each docking program must be installed and path to executable directory added to docking/config.sh."
    echo -e "See dependencies.txt for details on how to install docking programs and version information.\n" 
    exit
}

if [ $# -eq 0 ]; then
    print_help
    exit
fi

vina=false
ledock=false
rdock=false
autodock=false
plants=false
exhaustiveness=8
exhaust_given=false
blind=false
center_given=false
size_given=false
crystal_given=false



#read in command line arguments
IFS=' ' read -r -a allArgs <<< "$@"
for ((idx=0; idx<$#; idx++)); do
    if [ ${allArgs[${idx}]} == "-h" ]; then
        print_help
    elif [ ${allArgs[${idx}]} == "-r" ]; then
        rec_path=${allArgs[$(($idx+1))]}
        rec_path=$(realpath $rec_path)
    elif [ ${allArgs[${idx}]} == "-l" ]; then
        lig_path=${allArgs[$(($idx+1))]}
        lig_path=$(realpath $lig_path)
    elif [ ${allArgs[${idx}]} == "-c" ]; then
        crystal_path=${allArgs[$(($idx+1))]}
        crystal_path=$(realpath $crystal_path)
        crystal_given=true
    elif [ ${allArgs[${idx}]} == "-o" ]; then
        out_path=${allArgs[$(($idx+1))]}
        out_path=$(realpath $out_path)
        if ! [ -d "${out_path}" ]; then
            mkdir $out_path
            echo $out_path
        fi
    elif [ ${allArgs[${idx}]} == "--vina" ]; then
        vina=true
        echo "Will run Autodock Vina docking."
    elif [ ${allArgs[${idx}]} == "-e" ]; then
        exhaustiveness=${allArgs[$(($idx+1))]}
        exhaust_given=true
    elif [ ${allArgs[${idx}]} == "--ledock" ]; then
        ledock=true
        echo "Will run Ledock docking."
    elif [ ${allArgs[${idx}]} == "--rdock" ]; then
        rdock=true
        echo "Will run rDock docking."
    elif [ ${allArgs[${idx}]} == "--autodock" ]; then
        autodock=true
        echo "Will run Autodock docking."
    elif [ ${allArgs[${idx}]} == "--plants" ]; then
        plants=true
        echo "Will run PLANTS docking."
    elif [ ${allArgs[${idx}]} == "--blind" ]; then
        blind=true
    elif [ ${allArgs[${idx}]} == "--center" ]; then
        x_coord=${allArgs[$(($idx+1))]}
        y_coord=${allArgs[$(($idx+2))]}
        z_coord=${allArgs[$(($idx+3))]}
        center_given=true
    elif [ ${allArgs[${idx}]} == "--size" ]; then
        x_size=${allArgs[$(($idx+1))]}
        y_size=${allArgs[$(($idx+2))]}
        z_size=${allArgs[$(($idx+3))]}
        size_given=true
    fi
done

script_path=$(dirname "$0")

if [ "$crystal_given" = true ] && ! ([ "${center_given}" = true ] && [ "${size_given}" = true ]); then
    echo "Generating box from $crystal_path"
    #write something smarter for generating the temp file name. Maybe add an idx that iterates up if the file of that idx exists.
    coords_file="coords_${RANDOM}.txt"
    python ./docking/getDimensions.py -r $crystal_path -o $coords_file
    raw_coords=$(cat $coords_file)
    rm $coords_file
    raw_coords=($raw_coords)
    x_coord=${raw_coords[1]}
    y_coord=${raw_coords[2]}
    z_coord=${raw_coords[3]}
    x_size=${raw_coords[5]}
    y_size=${raw_coords[6]}
    z_size=${raw_coords[7]}
    center_given=true
    size_given=true
fi

#start vina runs

if [ "$vina" = true ]; then
    echo "Running Vina..."
    vinacmd="./docking/runVina.sh -l ${lig_path} -r ${rec_path} -o ${out_path}"
    
    if [ "${exhaust_given}" = true ]; then
        vinacmd="${vinacmd} -e ${exhaustiveness}"
    fi
    
    #add coordinates to command
    if [ "$center_given" = true ] && [ "$size_given" = true ]; then
        vinacmd="${vinacmd} --center $x_coord $y_coord $z_coord --size $x_size $y_size $z_size" 
    elif [ "$blind" = true ]; then
        vinacmd="${vinacmd} --blind"
    else
        echo "Search space coordinates or --blind specification required for Vina to run."
        vinacmd="" 
    fi
    
    $vinacmd
fi

if [ "$ledock" = true ]; then
    echo "Running Ledock..."
    ledockcmd="./docking/runLedock.sh -l ${lig_path} -r ${rec_path} -o ${out_path}"
    
    #add coordinates to command
    if [ "$center_given" = true ] && [ "$size_given" = true ]; then
        ledockcmd="${ledockcmd} --center $x_coord $y_coord $z_coord --size $x_size $y_size $z_size" 
    elif [ "$blind" = true ]; then
        ledockcmd="${ledockcmd} --blind"
    else
        echo "Search space coordinates or --blind specification required for Ledock to run."
        ledockcmd="" 
    fi
    
    $ledockcmd
fi

if [ "$rdock" = true ]; then
    echo "Running rDock..."
    rdockcmd="./docking/runrDock.sh -l ${lig_path} -r ${rec_path} -o ${out_path}"
    #add coordinates to command
    if [ "$center_given" = true ] && [ "$size_given" = true ]; then
        radius=$(python ./docking/sizeToRadius.py $x_size $y_size $z_size)
        rdockcmd="${rdockcmd} --center $x_coord $y_coord $z_coord --size $radius"
    elif [ "$blind" = true ]; then
        rdockcmd="${rdockcmd} --blind"
    else
        echo "Search space coordinates or --blind specification required for rDock to run."
        rdockcmd=""
    fi
    $rdockcmd
fi

if [ "$autodock" = true ]; then
    echo "Running Autodock..."
    autodockcmd="./docking/runAutodock.sh -l ${lig_path} -r ${rec_path} -o ${out_path}"
    
    #add coordinates to command
    if [ "$center_given" = true ]; then
        autodockcmd="${autodockcmd} --center $x_coord $y_coord $z_coord" 
    else
        echo "Search space coordinates required for Autodock to run."
        autodockcmd="" 
    fi
    
    $autodockcmd
fi

if [ "$plants" = true ]; then
    echo "Running PLANTS..."
    plantscmd="./docking/runPlants.sh -l ${lig_path} -r ${rec_path} -o ${out_path}"
    
    #add coordinates to command
    if [ "$center_given" = true ] && [ "$size_given" = true ]; then
        radius=$(python ./docking/sizeToRadius.py $x_size $y_size $z_size --plants)
        plantscmd="${plantscmd} --center $x_coord $y_coord $z_coord --size $radius" 
    else
        echo "Search space coordinates required for rDock to run."
        plantscmd="" 
    fi
    
    $plantscmd
fi
