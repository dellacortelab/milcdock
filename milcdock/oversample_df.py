#script that takes a df as input and saves a version with all targets and ligands oversampled. 
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--input',help='Path to input pkled DF',required=True)
parser.add_argument('-o','--out',help='Path to write oversampled df to (should end in .pkl)',required=True)
args = parser.parse_args()
path = args.input
out_path = args.out

def oversample_df(df):
    '''Oversamples active/decoy ratio and total ligands per target ratio
        Returns oversampled dataframe'''
    new_df = pd.DataFrame(columns=df.columns)
    #find # of times to multiply actives
    #use df.append(active_df) that many times to expand tdf
    print('Starting ligand oversampling')
    for target in df['receptor'].unique():
        print(target)
        tdf = df[df['receptor'] == target]
        tdf = tdf.reset_index(drop=True)
        actives = tdf[tdf['active'] == 1]
        num_active = len(actives)
        decoys = tdf[tdf['active'] == 0]
        num_decoy = len(decoys)
        active_oversample_rate = num_decoy // num_active
        if active_oversample_rate <= 1:
            print('Active oversample rate only 1, skipping.')
            continue
        for i in range(active_oversample_rate-1):
            tdf = pd.concat([tdf,actives])
        new_df = pd.concat([new_df,tdf])
    max_ligands = 0
    for target in new_df['receptor'].unique():
        tdf = new_df[new_df['receptor'] == target]
        t_ligands = len(tdf)
        if t_ligands > max_ligands:
            max_ligands = t_ligands
    
    print('Starting target oversampling')
    for target in new_df['receptor'].unique():
        print(target)
        tdf = new_df[new_df['receptor'] == target]
        num_ligs = len(tdf)
        oversample_rate = max_ligands // num_ligs
        if oversample_rate > 1:
            for i in range(oversample_rate-1):
                new_df = pd.concat([new_df,tdf])
            
    new_df = new_df.reset_index(drop=True)
    return new_df

def reduce_decoys(df,max_decoys = 80000):
    for rec in df['receptor'].unique():
        rdf = df[df['receptor'] == rec]
        num_decoy = len(rdf[rdf['active'] == 0])
        if num_decoy > max_decoys:
            drop_num = num_decoy - max_decoys
            decoy_idxs = rdf[rdf['active'] == 0].index
            drop_idxs = np.random.choice(decoy_idxs,drop_num,replace=False)
            df = df.drop(drop_idxs)
            print(f'Dropped {drop_num} decoys on {rec}')
            rdf = df[df['receptor'] == rec]
            reduced_num_decoy = len(rdf[rdf['active'] == 0])
            print(f'Went from {num_decoy} to {reduced_num_decoy} decoys')
    return df
            

df = pd.read_pickle(path)

#some failed jobs put inf values in plants_SCORE_NORM_CONTACT column. Remove these.
bad_idxs = []
df.reset_index(inplace=True,drop=True)
print('Removing bad ligands...')
for i,val in enumerate(df['plants_SCORE_NORM_CONTACT']):
    if np.isposinf(val):
        bad_idxs.append(i)
df = df.drop(bad_idxs)
df.reset_index(inplace=True,drop=True)

#limit each target to max 80k ligands
print('Reducing dataframe...')
df = reduce_decoys(df)

print('Oversampling dataframe...')
df = oversample_df(df)

df.to_pickle(out_path)