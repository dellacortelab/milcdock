import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import os
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import sys

def convertUniToDict(uni,newNames = [],add_idx = True):
    allCoords = dict()
    atom_types = dict()
    for i,atom in enumerate(uni.atoms):
        if not newNames:
            name = atom.name
        elif newNames[i].startswith(atom.name[0]):
            name = newNames[i]
        else:
            print('Atom names improperly mapped.')
            sys.exit(1)
        if add_idx:
            if name not in atom_types:
                atom_types[name] = 1
            name_type = name
            if not name.endswith('*') and not name == 'OXT' and not name[-1].isdigit():
                name = name + str(atom_types[name])
            atom_types[name_type] += 1
        coord = atom.position
        allCoords[name] = coord
    return allCoords

def modifyPDBFiles(filePath):
    #Use openbabel to convert autodock .out files from pdbqt to pdb format before using this function.
    #Keep the .out and the .pdb files in the same directory under the same name except for the file extension.
    #This function is also used to modify dock6 files that have been converted from mol2 to pdb format for
    #input into MDAnalysis.Universe
    with open(filePath, "r") as f:
        fileLines = f.readlines()
    fileToWrite = open(filePath, "w")
    for line in fileLines:
        if not line.startswith("MODEL") and not line.startswith("CONECT") and not line.startswith("END"):
            fileToWrite.write(line)
    fileToWrite.close()
    return

def make_names_unique(names):
    name_types = dict()
    unique_names = []
    for name in names:
        if name not in name_types:
            name_types[name] = 1
        name_type = name
        if not name.endswith('*') and not name == 'OXT' and not name[-1].isdigit():
            name = name + str(name_types[name])
        unique_names.append(name)
        name_types[name_type] += 1
    return unique_names

#helper functions to build dataframe
def appendLedockDataToDataframe(ledock_path, dataFrame, proteinModel):
    dockingTool = "ledock"
    receptorName = proteinModel
    df_list = []
    failed_ligs = 0
    for file in os.listdir(ledock_path):
        if not file.endswith("dok"):
            continue
        #get length of base ligand to extract single ligands from mda universe
        with open(os.path.join(ledock_path, file)) as f:
            f.readline()
            file_lines = f.readlines()
        ligandName = file[:file.rfind(".")]
        atom_names = []
        for line in file_lines:
            if line.startswith('ATOM'):
                atom_names.append(line.split()[2])
            elif line.startswith('END'):
                break
        atom_names = make_names_unique(atom_names)
        scoresList = []
        poseAmount = []
        lig_idx = 0
        poses = []
        for line in file_lines:
            if line.startswith("REMARK"):
                scoresList.append(float(line.split()[7]))
                poseAmount.append(int(line.split()[5]))
            elif line.startswith('ATOM'):
                line_split = line.split()
                coords = np.array(line_split[5:8],dtype=float)
                if len(poses) == lig_idx:
                    poses.append(dict())
                    atom_idx = 0
                atom_name = atom_names[atom_idx]
                if not atom_name.startswith('H') and not atom_name.startswith('1H') and not atom_name.startswith('2H'):
                    pose_dict = poses[lig_idx]
                    pose_dict[atom_name] = coords
                atom_idx += 1
            elif line.startswith('END'):
                lig_idx += 1
        if len(poses) == 0:
            failed_ligs += 1
            continue
        for i, score in enumerate(scoresList):
            df_list.append([score,ligandName,receptorName,poses[i],dockingTool,poseAmount[i]])
    ledock_frame = pd.DataFrame(df_list,columns=[i for i in dataFrame.columns])
    dataFrame = pd.concat([dataFrame,ledock_frame],ignore_index=True)
    if failed_ligs != 0:
        print(f'{failed_ligs} ligands failed.')
    return dataFrame

def extractAutodockScores(path_to_out_file):
    #This function gets scores from autodock output files.
    scores_list = []
    with open(path_to_out_file) as f:
        file_lines = f.readlines()
    for line in file_lines:
        if line.split()[0] == "REMARK" and line.split()[1] == "VINA":
            scores_list.append(line.split()[3])
    return scores_list

def appendVinaDataToDataframe(vina_path, dataFrame, proteinModel):
    dockingTool = "vina"
    receptorName = proteinModel
    df_list = []
    for file in sorted(os.listdir(vina_path)):
        if file.endswith(".pdbqt"):
            scores = extractAutodockScores(os.path.join(vina_path, file))
            modifyPDBFiles(os.path.join(vina_path, file[:-6] + ".pdb"))
            try:
                ligandUniverse = mda.Universe(os.path.join(vina_path, file[:-6] + ".pdb"))
            except:
                print(f"Unable to create MDA universe from: {file}")
            ligandUniverse = ligandUniverse.select_atoms('not type H')
            atomNames = []
            hasLen = False
            for i,atom in enumerate(ligandUniverse.atoms):
                name = atom.name
                if name in atomNames:
                    ligandLength = i
                    hasLen = True
                    break
                atomNames.append(name)
            if not hasLen:
                #when this error pops up it is due to having only 1 pose in the file.
                ligandLength = len(ligandUniverse.atoms)
            ligandName = file[:file.rfind(".")]
            for i in range(len(scores)):
                #convert to dict to allow pickling
                uniCoords = convertUniToDict(ligandUniverse.atoms[ligandLength*i:ligandLength*(i+1)],add_idx=False)
                df_list.append([float(scores[i]),ligandName,receptorName,uniCoords,dockingTool])
    vina_frame = pd.DataFrame(df_list,columns=[i for i in dataFrame.columns])
    dataFrame = pd.concat([dataFrame,vina_frame],ignore_index=True)
    return dataFrame

def appendrDockDataToDataframe(rdock_path, dataFrame, proteinModel):
    #depends on original mol2 files being copied over into the same directory.
    rdockDataFrame = pd.DataFrame(columns=[i for i in dataFrame.columns])
    dockingTool = "rdock"
    receptorName = proteinModel
    df_list = []
    for file in sorted(os.listdir(rdock_path)):
        if file.endswith("_out.sd"):
            modifyPDBFiles(os.path.join(rdock_path, file[:file.rfind(".")] + ".pdb"))
            universe = mda.Universe(os.path.join(rdock_path, file[:file.rfind(".")] + ".pdb"),format="PDB")
            universe = universe.select_atoms('not type H')
            origLigFile = file[:file.find('_out.sd')] + '.mol2'
            origAtoms = mda.Universe(os.path.join(rdock_path,origLigFile))
            origAtoms = origAtoms.select_atoms('not type H')
            with open(os.path.join(rdock_path, file), "r") as f:
                fileLines = f.readlines()

            for index,line in enumerate(fileLines):
                if line.find("<Name>") != -1:
                    realLigandName = fileLines[index + 1][:-1] #real name of specific ligand from inside file
                    break

            #get length of ligand, create universes with each ligand type
            ligandLength = len(origAtoms)
            atomNames = [atom.name for atom in origAtoms.atoms]
            ligandName = file[:file.rfind('_')]
            ligandScores = []
            score_inter = []
            score_inter_vdw = []
            score_inter_norm = []
            score_intra = []
            score_norm = []
            for index, line in enumerate(fileLines):
                if line.find("<SCORE>") != -1:
                    ligandScore = float(fileLines[index + 1])
                    ligandScores.append(ligandScore)
                elif line.find("<SCORE.INTER>") != -1:
                    score_inter.append(float(fileLines[index+1]))
                elif line.find("<SCORE.INTER.VDW>") != -1:
                    score_inter_vdw.append(float(fileLines[index+1]))
                elif line.find("<SCORE.INTER.norm>") != -1:
                    score_inter_norm.append(float(fileLines[index+1]))
                elif line.find("<SCORE.INTRA>") != -1:
                    score_intra.append(float(fileLines[index+1]))
                elif line.find("<SCORE.norm>") != -1:
                    score_norm.append(float(fileLines[index+1]))
            for i in range(len(ligandScores)):
                uniCoords = convertUniToDict(universe.atoms[ligandLength*i:ligandLength*(i+1)],atomNames)
                df_list.append([ligandScores[i],score_inter[i],score_inter_vdw[i],score_inter_norm[i],score_intra[i],score_norm[i],ligandName,receptorName,uniCoords,dockingTool,realLigandName])
    rdock_frame = pd.DataFrame(df_list,columns=[i for i in dataFrame.columns])
    dataFrame = pd.concat([dataFrame,rdock_frame],ignore_index=True)
    return dataFrame

def get_lig_length(uni):
    atomNames = []
    for i,atom in enumerate(uni.atoms):
        name = atom.name
        if name in atomNames:
            ligandLength = i
            return ligandLength
        atomNames.append(name)
    print('Unable to find length of ligand.')
    return None

def appendPLANTSDataToDataframe(plants_path, dataFrame, proteinModel):
    dockingTool = "plants"
    df_list = []
    receptorName = proteinModel
    with open(os.path.join(plants_path,'ligands.txt'),'r') as f:
        lig_file_paths = f.readlines()
    lig_names = []
    for path in lig_file_paths:
        lig_file = os.path.basename(path)
        lig_names.append(lig_file[:lig_file.rfind('.mol2')])
    features = pd.read_csv(os.path.join(plants_path,'features.csv'))
    #use openbabel to convert docked_ligands.mol2 to pdb before running.
    modifyPDBFiles(os.path.join(plants_path,'docked_ligands.pdb'))
    uni = mda.Universe(os.path.join(plants_path,'docked_ligands.pdb'))
    uni = uni.atoms
    #assume order in features.csv matches order in docking_ligands.mol2
    with open(os.path.join(plants_path,'docked_ligands.mol2'),'r') as f:
        mol2_lines = f.readlines()
    lig_lengths = []
    for i,line in enumerate(mol2_lines):
        if line.startswith('@<TRIPOS>ATOM'):
            atoms_idx = i
        elif line.startswith('@<TRIPOS>BOND'):
            bond_idx = i
            lig_length = bond_idx - (atoms_idx+1)
            lig_lengths.append(lig_length)
    lig_idx = 0
    lig_conf = -1
    i = 0
    for idx,row in features.iterrows():
        lig_length = lig_lengths[i]
        i += 1
        new_lig_conf = int(idx[idx.rfind('_')+1:])
        if new_lig_conf < lig_conf:
            lig_idx += 1
        coords = convertUniToDict(uni[:lig_length].select_atoms('not type H'))
        uni = uni[lig_length:]
        real_name = idx[:idx.find('_')]
        name = lig_names[lig_idx]
        new_row = [name,real_name,proteinModel,dockingTool,coords]
        new_row.extend(row.iloc[:5])
        new_row.append(row['SCORE_NORM_CONTACT'])
        new_row.append(row['PLPtotal'])
        new_row.append(row['CHEMparthbond'])
        df_list.append(new_row)
        lig_conf = new_lig_conf
    plants_frame = pd.DataFrame(df_list,columns=[i for i in dataFrame.columns])
    dataFrame = pd.concat([dataFrame,plants_frame],ignore_index=True)
    return dataFrame

def appendAutodockDataToDataframe(dir_path,df,receptor_name):
    #consider reading in charge to normalize scores by charge
    auto_df = pd.DataFrame(columns = [i for i in df.columns])
    docking_tool = 'autodock4'
    for file_name in sorted(os.listdir(dir_path)):
        if file_name.endswith(".dlg"):
            #read in scores and atom coordinates
            #I don't think we'll be able to use a universe to read in coords with this file type.
            ligand_name = file_name[:file_name.rfind('_')]
            with open(os.path.join(dir_path,file_name)) as f:
                f_lines = f.readlines()
            # atom_map_file = file_name[:file_name.rfind('.')] + '.map'
            # with open(os.path.join(dir_path,atom_map_file),'r') as f:
            #     atom_map = f.readlines()
            # for i,atom in enumerate(atom_map):
            #     atom_map[i] = atom.strip()
            #initialize lists to save poses and score values
            lig_idx = 0
            poses = []
            scores = []
            #inhibition_scores = [] # don't use inhibition constant for now
            final_intermol_energy = []
            vdw_hbond_desolv_energy = []
            electro_energy = []
            final_internal_energy = []
            tors_energy = []
            for line in f_lines:
                #don't use any clustered data for now. Just stick with standard setup. Consider using some later.
                if line.startswith('DOCKED: USER'):
                    line_split = line.split()
                    if len(line_split) > 5:
                        if line_split[2] == 'Estimated' and line_split[3] == 'Free':
                            scores.append(float(line_split[8]))
                        elif line_split[2] == '(1)' and line_split[4] == 'Intermolecular':
                            final_intermol_energy.append(float(line_split[7]))
                        elif line_split[2] == 'vdW' and line_split[3] == '+':
                            vdw_hbond_desolv_energy.append(float(line_split[9]))
                        elif line_split[2] == 'Electrostatic' and line_split[3] == 'Energy':
                            electro_energy.append(float(line_split[5]))
                        elif line_split[2] == '(2)' and line_split[4] == 'Total':
                            final_internal_energy.append(float(line_split[8]))
                        elif line_split[2] == '(3)' and line_split[3] == 'Torsional':
                            tors_energy.append(float(line_split[7]))
                elif line.startswith('DOCKED: ATOM'):
                    line_split = line.split()
                    atom_name = line_split[3]
                    firstCoordIdx = 6
                    lastCoordIdx = 9
                    if atom_name.startswith('H'):
                        continue
                    elif atom_name.startswith("1H") or atom_name.startswith("2H"):
                        continue
                    try:
                        if float(line_split[firstCoordIdx]) == 1:
                            firstCoordIdx += 1
                            lastCoordIdx += 1
                        atom_pos = np.array(line_split[firstCoordIdx:lastCoordIdx],dtype=float)
                    except ValueError:
                        #fix error where z coordinate combines with the next column's terms.
                        atom_pos = line_split[firstCoordIdx:lastCoordIdx]
                        i = 2 #idx for z coordinate
                        val = atom_pos[i]
                        for j,char in enumerate(val):
                            if j > 0 and not char.isnumeric() and char != '.':
                                val = val[:j]
                                atom_pos[i] = val
                                break
                        atom_pos = np.array(atom_pos,dtype=float)
                    if len(poses) == lig_idx:
                        poses.append(dict())
                        atom_idx = 0
                    # try:
                    #     atom_name = atom_map[atom_idx]
                    # except IndexError:
                    #     print("INDEX ERROR! Ligand: ", ligand_name)
                    #     quit()
                    pose_dict = poses[lig_idx]
                    pose_dict[atom_name] = atom_pos #may need to reassign dict to list. Not sure if pass by reference or value.
                    atom_idx += 1
                elif line.startswith('DOCKED: ENDMDL'):
                    lig_idx += 1
            for i,score in enumerate(scores):
                new_row = [ligand_name,receptor_name,docking_tool,poses[i],score,final_intermol_energy[i],vdw_hbond_desolv_energy[i],electro_energy[i],final_internal_energy[i],tors_energy[i]]
                auto_df.loc[len(auto_df)] = new_row
    df = df.append(auto_df,ignore_index=True)
    return df

def getRMSD(uni1,uni2):
    """Calculate RMSD between two ligands based on atom names. Ligand length and atom names must match. DOESN'T WORK WITH CURRENT SETUP SWITCH TO DICT COMPATABILITY"""
    uni1Names = uni1.atoms.names
    uni2Names = uni2.atoms.names
    uni1Positions = uni1.atoms.positions
    uni2Positions = uni2.atoms.positions
    squareSum = 0
    for uni1Idx,name in enumerate(uni1Names):
        uni2Idx = np.where(uni2Names == name)
        squareSum += (np.linalg.norm(uni1Positions[uni1Idx] - uni2Positions[uni2Idx]))**2
    RMSD = np.sqrt(1/len(uni1Names)*squareSum)
    return RMSD


def convert_docking_output_files_to_pdb(path, vina, rdock, plants):
    if vina:
        os.system("cd {}/vina && obabel -ipdbqt *.pdbqt -opdb -m".format(path))
    if rdock:
        os.system("cd {}/rdock && obabel -isd *_out.sd -opdb -m".format(path))
    if plants:
        file = "{}/plants/docked_ligands.mol2".format(path)
        basename = os.path.basename(file).replace('.mol2', '')
        out_path = "{}/plants/{}.pdb".format(path, basename)
        os.system("obabel -imol2 {} -opdb -O {}".format(file, out_path))


def convert_docking_to_pkl(docking_output_dir, docking_pkl_dir, ledock=False, rdock=False, vina=False, plants=False, autodock=False):
    """
    Convert docking data into a dataframe then save it as a pkl.

    Parameters:
    docking_output_dir (str): path to docking data
    docking_pkl_dir (str): path to save pkl
    ledock (bool): generate pickle from Ledock data
    dock6 (bool): generate pickle from Dock6 data
    rdock (bool): generate pickle from rDock data
    vina (bool): generate pickle from Vina data
    plants (bool): generate pickle from Plants data
    autodock (bool): generate pickle from Autodock4 data
    all (bool): generate pickles from all available docking data

    Returns:
    None
    """
    targetName = os.path.basename(os.path.normpath(docking_output_dir))

    print(f"Creating dataframes for {targetName}")
    vinaFrame = pd.DataFrame(columns=['score','ligand name','protein model','universe','docking tool'])
    ledockFrame = pd.DataFrame(columns=['score','ligand name','protein model','universe','docking tool','num_poses'])
    rdockFrame = pd.DataFrame(columns=['score','score_inter','score_inter_vdw','score_inter_norm','score_intra','score_norm','ligand name','protein model','universe','docking tool','real_lig_name'])
    plantsFrame = pd.DataFrame(columns=['ligand name','real_lig_name','protein model','docking tool','universe','TOTAL_SCORE','SCORE_RB_PEN','SCORE_NORM_HEVATOMS','SCORE_NORM_CRT_HEVATOMS','SCORE_NORM_WEIGHT','SCORE_NORM_CONTACT','PLPtotal','CHEMparthbond'])
    autoFrame = pd.DataFrame(columns=['ligand name','protein model','docking tool','universe','score','final_intermol_energy','vdw_hbond_desolv_energy','electro_energy','final_internal_energy','tors_energy'])

    print(f'Reading in data from {docking_output_dir}')
    if vina:
        vinaFrame = appendVinaDataToDataframe(os.path.join(docking_output_dir, 'vina/'),vinaFrame,targetName)
    if ledock:
        ledockFrame = appendLedockDataToDataframe(os.path.join(docking_output_dir, 'ledock/'),ledockFrame,targetName)
    if rdock:
        rdockFrame = appendrDockDataToDataframe(os.path.join(docking_output_dir, 'rdock/'),rdockFrame,targetName)
    if plants:
        plantsFrame = appendPLANTSDataToDataframe(os.path.join(docking_output_dir, 'plants/'),plantsFrame,targetName)
    if autodock:
        autoFrame = appendAutodockDataToDataframe(os.path.join(docking_output_dir, 'autodock/'),autoFrame,targetName)

    if vina:
        fileName = f'{targetName}_vina.pkl'
        vinaFrame.to_pickle(os.path.join(docking_pkl_dir,fileName))
    if ledock:
        fileName = f'{targetName}_ledock.pkl'
        ledockFrame.to_pickle(os.path.join(docking_pkl_dir,fileName))
    if rdock:
        fileName = f'{targetName}_rdock.pkl'
        rdockFrame.to_pickle(os.path.join(docking_pkl_dir,fileName))
    if plants:
        fileName = f'{targetName}_plants.pkl'
        plantsFrame.to_pickle(os.path.join(docking_pkl_dir,fileName))
    if autodock:
        fileName = f'{targetName}_autodock.pkl'
        autoFrame.to_pickle(os.path.join(docking_pkl_dir,fileName))
    

def remove_extra_files(docking_output_dir) -> None:
    """
    Remove extra files from target directories after pickling to prep for tarring.

    Parameters:
    dir (str): path to the docking output target directory

    Returns:
    None
    """
    if not os.path.isdir(docking_output_dir):
        raise ValueError(f"{docking_output_dir} is not a valid directory.")

    if os.path.isdir(f"{docking_output_dir}/rdock"):
        os.system(f"rm {docking_output_dir}/rdock/*out.pdb")
        os.system(f"mkdir {docking_output_dir}/rdock/temp_save")
        os.system(f"mv {docking_output_dir}/rdock/*out.sd {docking_output_dir}/rdock/temp_save")
        os.system(f"rm {docking_output_dir}/rdock/*sd")
        os.system(f"mv {docking_output_dir}/rdock/temp_save/* {docking_output_dir}/rdock/")
        os.system(f"rm -r {docking_output_dir}/rdock/temp_save/")

    if os.path.isdir(f"{docking_output_dir}/vina"):
        os.system(f"rm {docking_output_dir}/vina/*pdb")

    if os.path.isdir(f"{docking_output_dir}/plants"):
        os.system(f"rm {docking_output_dir}/plants/*dat")
        os.system(f"rm {docking_output_dir}/plants/docked_ligands.pdb")
        os.system(f"rm {docking_output_dir}/plants/docked_proteins.mol2")
        os.system(f"rm {docking_output_dir}/plants/skippedligands.csv")
        os.system(f"rm {docking_output_dir}/plants/ranking.csv")
        os.system(f"rm {docking_output_dir}/plants/bestranking.csv")
        os.system(f"rm {docking_output_dir}/plants/correspondingNames.csv")
        os.system(f"rm {docking_output_dir}/plants/constraints.csv")