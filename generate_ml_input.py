import os
import pandas as pd
import numpy as np
import sys
from argparse import ArgumentParser

from docking.docking_to_pkl import convert_docking_output_files_to_pdb, convert_docking_to_pkl, remove_extra_files


parser = ArgumentParser()
parser.add_argument('--docking_output_dir', help='Path to a directory containing docking output files.', default='./data/1_docked_output/')
parser.add_argument('--docking_pkl_dir', help='Path to a directory containing only the pickles containing docking data. Each pickle file in the specified directory should be named using the format: TARGETNAME_DOCKINGTOOL.pkl', default='./data/2_docking_pkls/')
parser.add_argument('--ml_pkl_path', help='file to write machine learning inputs to (should have .pkl extension)', default='./data/3_ml_inputs/ml_inputs.pkl')
parser.add_argument('--dock_tools', help='Space-separated list of docking tools to use. Options are: vina, rdock, ledock, plants, autodock', nargs='+', default=['vina', 'rdock', 'ledock', 'plants', 'autodock'], choices=['vina', 'rdock', 'ledock', 'plants', 'autodock'])

def removePoses(df,maxLigs=3,col='score'):
    '''Remove extra poses from a df containing poses from only one ligand type'''
    if len(df) > maxLigs:
        scores = np.array(df[col],dtype=float)
        allidxs = np.argsort(scores)
        badidxs = allidxs[maxLigs:]
        removeIdxs = df.iloc[badidxs,:].index
        return df.drop(index = removeIdxs)
    else:
        return df

def getRMSD(coords_dict_1,coords_dict_2):
    atomNames = coords_dict_1.keys()
    coords1 = np.zeros((len(atomNames),3),dtype=float)
    coords2 = np.zeros((len(atomNames),3),dtype=float)
    for i,atom in enumerate(atomNames):
        coords1[i,:] = coords_dict_1[atom]
        coords2[i,:] = coords_dict_2[atom]
    rmsd = np.sqrt(np.mean(np.linalg.norm(coords1 - coords2,axis=1)**2))
    return rmsd
    
def avgRMSDFromTopScore(df,col='score'):
    '''get the avg RMSD of all poses from the top pose.'''
    topScore = min(df[col])
    topScoreDF = df[df[col] == topScore]
    topIdx = topScoreDF.index[0]
    topScoreCoords = list(topScoreDF['universe'])[0]
    rmsd = 0
    for idx in df.index:
        if idx != topIdx:
            rmsd += getRMSD(topScoreCoords,df.loc[idx,'universe'])
    if rmsd != 0:
        rmsd = rmsd / (len(df)-1) #-1 because # of RMSD is len(df)-1
    return rmsd

def getPairwiseRMSD(df1,df2,df1_col='score',df2_col='score'):
    '''returns the RMSD between the top scoring (most negative) pose from df1 and df2'''
    df1_topscore = min(df1[df1_col])
    topdf1 = df1[df1[df1_col] == df1_topscore]
    df1_coords = list(topdf1['universe'])[0]
    df2_topscore = min(df2[df2_col])
    topdf2 = df2[df2[df2_col] == df2_topscore]
    df2_coords = list(topdf2['universe'])[0]
    rmsd = getRMSD(df1_coords,df2_coords)
    return rmsd

def determineVal(row):
    if row['name'].startswith('active'):
        return 1
    else:
        return 0

def addActiveCol(df):
    actives = df.apply(lambda row: determineVal(row), axis = 1)
    # Only create this column if there are actives
    if actives.sum() > 0:
        df['active'] = actives

def generate_ml_input(docking_pkl_dir, ml_pkl_path, dockingTools):

    allTargets = []
    allFiles = []

    for file in os.listdir(docking_pkl_dir):
        if not file.endswith('.pkl'):
            continue
        allFiles.append(file)
        target = file[:file.rfind('_')]
        if target not in allTargets:
            allTargets.append(target)

    #generate columns for the final df
    df_columns = ['name','realname','receptor']
    for tool in dockingTools:
        df_columns.extend([tool + '_top_score'])
        df_columns.extend([tool + '_avg_3_score'])
        if tool == 'rdock':
            df_columns.extend([tool + '_top_score_inter'])
            df_columns.extend([tool + '_avg_3_score_inter'])
            df_columns.extend([tool + '_top_score_inter_vdw'])
            df_columns.extend([tool + '_avg_3_score_inter_vdw'])
            df_columns.extend([tool + '_top_score_inter_norm'])
            df_columns.extend([tool + '_avg_3_score_inter_norm'])
            df_columns.extend([tool + '_top_score_intra'])
            df_columns.extend([tool + '_avg_3_score_intra'])
            df_columns.extend([tool + '_top_score_norm'])
            df_columns.extend([tool + '_avg_3_score_norm'])
        elif tool == 'plants':
            df_columns.extend([tool + '_SCORE_RB_PEN'])
            df_columns.extend([tool + '_SCORE_NORM_HEVATOMS'])
            df_columns.extend([tool + '_SCORE_NORM_CRT_HEVATOMS'])
            df_columns.extend([tool + '_SCORE_NORM_WEIGHT'])
            df_columns.extend([tool + '_SCORE_NORM_CONTACT'])
            df_columns.extend([tool + '_PLPtotal'])
            df_columns.extend([tool + '_CHEMparthbond'])
        elif tool == 'autodock':
            df_columns.extend([tool + '_final_intermol_energy'])
            df_columns.extend([tool + '_vdw_hbond_desolv_desolv_energy'])
            df_columns.extend([tool + '_electro_energy'])
            df_columns.extend([tool + '_final_internal_energy'])
            df_columns.extend([tool + '_tors_energy'])
        df_columns.extend([tool + '_rng_3_score'])
        df_columns.extend([tool + '_avg_3_rmsd'])



    for i in range(len(dockingTools)):
        for j in range(len(dockingTools)):
            if i < j:
                name = f'RMSD_{dockingTools[i]}-{dockingTools[j]}'
                df_columns.extend([name])

    all_target_dfs = pd.DataFrame(columns=df_columns)

    #write df specific to each target, concatenate later.
    for i,target in enumerate(allTargets):
        targetDF = pd.DataFrame(columns=df_columns)
        targetFiles = [file for file in allFiles if file.startswith(target)]
        print(f"Processing {target}")
        hasVina = False
        hasrDock = False
        hasLedock = False
        hasPlants= False
        hasAutodock = False
        for file in targetFiles:
            if file.endswith('_vina.pkl'):
                dockingTool = 'vina'
            elif file.endswith('_ledock.pkl'):
                dockingTool = 'ledock'
            elif file.endswith('_rdock.pkl'):
                dockingTool = 'rdock'
            elif file.endswith('_plants.pkl'):
                dockingTool = 'plants'
            elif file.endswith('_autodock.pkl'):
                dockingTool = 'autodock'
            else:
                dockingTool = 'none'
                continue
            df = pd.read_pickle(os.path.join(docking_pkl_dir,file))
            if dockingTool == 'vina':
                vinaDF = df
                hasVina = True
            elif dockingTool == 'ledock':
                ledockDF = df
                hasLedock = True
            elif dockingTool == 'rdock':
                rdockDF = df
                hasrDock = True
            elif dockingTool == 'plants':
                plantsDF = df
                hasPlants= True
            elif dockingTool == 'autodock':
                autodockDF = df
                hasAutodock = True
        if not (hasVina and hasrDock and hasLedock and hasPlants and hasAutodock):
            print(f"Skipping {target} due to missing pkl.")
            continue
        vinaLigs = vinaDF['ligand name'].unique()
        ledockLigs = ledockDF['ligand name'].unique()
        rdockLigs = rdockDF['ligand name'].unique()
        plantsLigs = plantsDF['ligand name'].unique()
        autodockLigs = autodockDF['ligand name'].unique()
        #grab only the intersection of ligands that succeeded for all programs
        allLigs = [lig for lig in vinaLigs if lig in ledockLigs and lig in rdockLigs and lig in plantsLigs and lig in autodockLigs]
        print(f'allLigs has {len(allLigs)} ligands.')
        num_errors = 0
        vina_errors = 0
        rdock_errors = 0
        pairwise_rmsd_errors = 0
        va_errors = 0
        vl_errors = 0
        vp_errors = 0
        for j,lig in enumerate(allLigs):
            vdf = vinaDF[vinaDF['ligand name'] == lig]
            ldf = ledockDF[ledockDF['ligand name'] == lig]
            rdf = rdockDF[rdockDF['ligand name'] == lig]
            pdf = plantsDF[plantsDF['ligand name'] == lig]
            adf = autodockDF[autodockDF['ligand name'] == lig]
            vdf = removePoses(vdf)
            ldf = removePoses(ldf)
            rdf = removePoses(rdf)
            pdf = removePoses(pdf,col='TOTAL_SCORE')
            adf = removePoses(adf)
            vinaTopScore = min(vdf['score'])
            vinaAvgScore = np.mean(vdf['score'])
            if abs(vinaAvgScore) <= 1e-5:
                vinaRngScore = abs(max(vdf['score'])-min(vdf['score'])) #can't divide by 0
            else:
                vinaRngScore = abs((max(vdf['score'])-min(vdf['score']))/vinaAvgScore) #normalized as % of avg score
            try:
                vinaAvgRMSD = avgRMSDFromTopScore(vdf)
            except KeyError:
                num_errors += 1
                vina_errors += 1
                continue
            #read in scores for each docking tool
            rdockTopScore = min(rdf['score'])
            rdockAvgScore = np.mean(rdf['score'])
            rdockRngScore = abs((max(rdf['score'])-min(rdf['score']))/rdockAvgScore) #normalized as % of avg score
            #consider getting the scores here that correspond with the top scoring pose in general, rather than just the lowest of the three
            rdockTopScoreInter = min(rdf['score_inter'])
            rdockAvgScoreInter = np.mean(rdf['score_inter'])
            rdockTopScoreInterVDW = min(rdf['score_inter_vdw'])
            rdockAvgScoreInterVDW = np.mean(rdf['score_inter_vdw'])
            rdockTopScoreInterNorm = min(rdf['score_inter_norm'])
            rdockAvgScoreInterNorm = np.mean(rdf['score_inter_norm'])
            rdockTopScoreIntra = min(rdf['score_intra'])
            rdockAvgScoreIntra = np.mean(rdf['score_intra'])
            rdockTopScoreNorm = min(rdf['score_norm'])
            rdockAvgScoreNorm = np.mean(rdf['score_norm'])
            try:
                rdockAvgRMSD = avgRMSDFromTopScore(rdf)
            except IndexError:
                num_errors += 1
                rdock_errors += 1
                continue
            realLigName = list(rdf['real_lig_name'])[0]
            ledockTopScore = min(ldf['score'])
            ledockAvgScore = np.mean(ldf['score'])
            ledockRngScore = abs((max(ldf['score'])-min(ldf['score']))/ledockAvgScore) #normalized as % of avg score
            ledockAvgRMSD = avgRMSDFromTopScore(ldf)

            plantsTopScore = min(pdf['TOTAL_SCORE'])
            plantsAvgScore = np.mean(pdf['TOTAL_SCORE'])
            plantsRngScore = abs((max(pdf['TOTAL_SCORE'])-min(pdf['TOTAL_SCORE']))/plantsAvgScore) #normalized as % of avg score
            plants_score_rb_pen = min(pdf['SCORE_RB_PEN'])
            plants_score_norm_hevatoms= min(pdf['SCORE_NORM_HEVATOMS'])
            plants_score_norm_crt_hevatoms= min(pdf['SCORE_NORM_CRT_HEVATOMS'])
            plants_score_norm_weight = min(pdf['SCORE_NORM_WEIGHT'])
            plants_score_norm_contact = min(pdf['SCORE_NORM_CONTACT'])
            plants_plptotal = min(pdf['PLPtotal'])
            plants_chemparthbond = min(pdf['CHEMparthbond'])
            plantsAvgRMSD = avgRMSDFromTopScore(pdf,col='TOTAL_SCORE')

            autodockTopScore = min(adf['score'])
            autodockAvgScore = np.mean(adf['score'])
            autodockRngScore = abs((max(adf['score'])-min(adf['score']))/autodockAvgScore) #normalized as % of avg score
            autodockAvgRMSD = avgRMSDFromTopScore(adf)
            autodock_final_intermol_energy = min(adf['final_intermol_energy'])
            autodock_vdw_hbond_desolv = min(adf['vdw_hbond_desolv_energy'])
            autodock_electro_energy = min(adf['electro_energy'])
            autodock_final_internal_energy = min(adf['final_internal_energy'])
            autodock_tors_energy = min(adf['tors_energy'])

            try:
                vinaledock = getPairwiseRMSD(vdf,ldf)
            except KeyError:
                num_errors += 1
                vl_errors += 1
                continue
            vinardock = getPairwiseRMSD(vdf,rdf)
            try:
                vinaplants = getPairwiseRMSD(vdf,pdf,df2_col='TOTAL_SCORE')
            except KeyError:
                num_errors += 1
                vp_errors += 1
                continue
            try:
                vinaautodock = getPairwiseRMSD(vdf,adf)
            except KeyError:
                num_errors += 1
                va_errors += 1
                continue
            ledockrdock = getPairwiseRMSD(ldf,rdf)
            ledockplants = getPairwiseRMSD(ldf,pdf,df2_col='TOTAL_SCORE')
            ledockautodock = getPairwiseRMSD(ldf,adf)
            rdockplants = getPairwiseRMSD(rdf,pdf,df2_col='TOTAL_SCORE')
            rdockautodock = getPairwiseRMSD(rdf,adf)
            autodockplants = getPairwiseRMSD(adf,pdf,df2_col='TOTAL_SCORE')
            newRow = [lig,
                    realLigName,
                    target,
                    vinaTopScore,
                    vinaAvgScore,
                    vinaRngScore,
                    vinaAvgRMSD,
                    rdockTopScore,
                    rdockAvgScore,
                    rdockTopScoreInter,
                    rdockAvgScoreInter,
                    rdockTopScoreInterVDW,
                    rdockAvgScoreInterVDW,
                    rdockTopScoreInterNorm,
                    rdockAvgScoreInterNorm,
                    rdockTopScoreIntra,
                    rdockAvgScoreIntra,
                    rdockTopScoreNorm,
                    rdockAvgScoreNorm,
                    rdockRngScore,
                    rdockAvgRMSD,
                    ledockTopScore,
                    ledockAvgScore,
                    ledockRngScore,
                    ledockAvgRMSD,
                    plantsTopScore,
                    plantsAvgScore,
                    plants_score_rb_pen,
                    plants_score_norm_hevatoms,
                    plants_score_norm_crt_hevatoms,
                    plants_score_norm_weight,
                    plants_score_norm_contact,
                    plants_plptotal,
                    plants_chemparthbond,
                    plantsRngScore,
                    plantsAvgRMSD,
                    autodockTopScore,
                    autodockAvgScore,
                    autodock_final_intermol_energy,
                    autodock_vdw_hbond_desolv,
                    autodock_electro_energy,
                    autodock_final_internal_energy,
                    autodock_tors_energy,
                    autodockRngScore,
                    autodockAvgRMSD,
                    vinardock,
                    vinaledock,
                    vinaplants,
                    vinaautodock,
                    ledockrdock,
                    rdockplants,
                    rdockautodock,
                    ledockplants,
                    ledockautodock,
                    autodockplants]
            targetDF.loc[len(targetDF)] = newRow
            if j % 1000 == 0:
                print(j)

        print(f"Finished with {vina_errors} vina errors")
        print(f"Finished with {rdock_errors} rdock errors")
        print(f"Finished with {vl_errors} vina-ledock errors")
        print(f"Finished with {vp_errors} vina-plants errors")
        print(f"Finished with {va_errors} vina-autodock errors")
        print(f"Finished with {num_errors} ligand errors")

        concat_frames = [all_target_dfs,targetDF]
        all_target_dfs = pd.concat(concat_frames)

    addActiveCol(all_target_dfs)
    print(all_target_dfs.head())
    print(f"Writing pkl to {ml_pkl_path}")
    all_target_dfs.to_pickle(ml_pkl_path)


if __name__ == "__main__":
    args = parser.parse_args()

    # Create boolean flags for each docking tool
    vina = False
    rdock = False
    ledock = False
    plants = False
    autodock = False
    if 'vina' in args.dock_tools:
        vina = True
    if 'rdock' in args.dock_tools:
        rdock = True
    if 'ledock' in args.dock_tools:
        ledock = True
    if 'plants' in args.dock_tools:
        plants = True
    if 'autodock' in args.dock_tools:
        autodock = True
    
    # Convert all docking file types to PDB
    # convert_docking_output_files_to_pdb(args.docking_output_dir, vina=args.vina, rdock=args.rdock, plants=args.plants)
    convert_docking_output_files_to_pdb(args.docking_output_dir, vina=vina, rdock=rdock, plants=plants)

    # Convert docking files to pkl format
    convert_docking_to_pkl(args.docking_output_dir, args.docking_pkl_dir, vina=vina, rdock=rdock, ledock=ledock, plants=plants, autodock=autodock)

    # Remove extra files
    remove_extra_files(args.docking_output_dir)

    # Convert pkl docking files to a pkl of compiled inputs to an ML program
    generate_ml_input(docking_pkl_dir=args.docking_pkl_dir, ml_pkl_path=args.ml_pkl_path, dockingTools=args.dock_tools)
