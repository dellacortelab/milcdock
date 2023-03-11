import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader, Sampler

dataset_dict = {'dude':'./data/3_ml_inputs/dude.pkl', 'lit-pcba':'./data/3_ml_inputs/lit-pcba.pkl', 'dude-lit-pcba':'./data/3_ml_inputs/dude-lit-pcba.pkl', 'dude-mini':'./data/3_ml_inputs/dude-mini.pkl', 'cache_1':'./data/3_ml_inputs/cache_combined.pkl', 'mcule':'./data/3_ml_inputs/mcule_all_ml_inputs.pkl', 'dude-2':'./data/3_ml_inputs/dude-2.pkl', 'dude-lit-pcba-xgboost-oversampled':'./data/3_ml_inputs/dude-lit-pcba_xgboost_80k.pkl', 'cache_1_part_2':'./data/3_ml_inputs/cache_merged_selection.pkl', 'cache_2':'./data/3_ml_inputs/cache_2_merged.pkl'}

dude_train_receptors = [['hdac2', 'cp3a4', 'ada17', 'pygm', 'mmp13', 'lkha4', 'bace1', 'dpp4', 'fnta', 'pde5a', 'hdac8'], ['hivrt', 'pnph', 'aces', 'gria2', 'hxk4', 'fpps', 'grik1'], ['pparg', 'fa10', 'try1', 'urok'], ['ptn1', 'pa2ga', 'fa7', 'tryb1', 'thrb', 'fabp4', 'ital', 'def', 'kif11'], ['hivint', 'hivpr', 'kith', 'parp1', 'nram', 'xiap', 'casp3', 'ampc', 'hs90a', 'cah2']]
lit_pcba_train_receptors = [['FEN1'], ['PKM2'], ['PPARG'], ['KAT2A', 'TP53'], []]
dude_val_receptors = [['ppard', 'ppara', 'andr', 'prgr', 'gcr', 'fkb1a'], ['esr1', 'thb', 'esr2', 'rxra', 'mcr', 'adrb1', 'cxcr4', 'drd3'], ['cdk2', 'mk10', 'rock1', 'akt2', 'akt1', 'abl1', 'mk01', 'braf']]
lit_pcba_val_receptors = [['MTORC1'], ['ESR1_ago', 'ESR1_ant', 'VDR', 'OPRK1'], ['MAPK1']]
dude_test_receptors = ['mk14', 'lck', 'igf1r', 'wee1', 'pur2', 'vgfr2', 'kpcb', 'fgfr1', 'hmdh', 'reni', 'glcm', 'adrb2', 'egfr', 'jak2', 'tgfr1', 'met', 'mapk2', 'src', 'fak1', 'plk1', 'csf1r']
lit_pcba_test_receptors = ['GBA', 'ALDH1', 'IDH1', 'ADRB2']

dude_excluded_test_receptors = ['aldr', 'aofb', 'comt', 'cp2c9', 'dhi1', 'dyr', 'inha', 'mp2k1', 'nos1', 'pgh1', 'pgh2', 'pyrd', 'sahh', 'tysy']

unbiased_features = ['vina_avg_3_rmsd', 'rdock_top_score_inter_norm', 'rdock_avg_3_score_inter_norm', 'rdock_top_score_norm', 'rdock_avg_3_score_norm', 'rdock_avg_3_rmsd', 'ledock_avg_3_rmsd', 'plants_SCORE_NORM_HEVATOMS', 'plants_SCORE_NORM_CRT_HEVATOMS', 'plants_SCORE_NORM_WEIGHT', 'plants_SCORE_NORM_CONTACT', 'plants_avg_3_rmsd', 'autodock_avg_3_rmsd', 'RMSD_vina-rdock', 'RMSD_vina-ledock', 'RMSD_vina-plants', 'RMSD_vina-autodock', 'RMSD_rdock-ledock', 'RMSD_rdock-plants', 'RMSD_rdock-autodock', 'RMSD_ledock-plants', 'RMSD_ledock-autodock', 'RMSD_plants-autodock']

def flatten(l):
    return [item for sublist in l for item in sublist]


class DockingDataset(Dataset):
    def __init__(self, dataset=None, dataset_path=None, mode='train', fold=None, target=None, dock_tools=['vina', 'autodock', 'plants', 'ledock', 'rdock']):
        """
        Args:
            dataset (string): Name of the dataset to use. Must provide either this argument or the dataset_path argument.
            dataset_path (string): Path to a dataset. Alternate way to specify dataset.
            mode (string): One of 'train', 'val', 'test'. Only necessary if specifying known dataset with dataset argument.
        """
        if dataset_path is not None:
            input_path = dataset_path
        else:
            input_path = dataset_dict[dataset]
        
        print(f'Reading data from {input_path}')
        df = pd.read_pickle(input_path)

        if dataset_path is not None:
            self.train_receptors = df['receptor'].unique().tolist()
            self.val_receptors = df['receptor'].unique().tolist()
            self.test_receptors = df['receptor'].unique().tolist()
        elif dataset == 'dude':
            if fold is not None:
                self.train_receptors = flatten([cluster for i, cluster in enumerate(dude_train_receptors + dude_val_receptors) if i != fold])
                self.val_receptors = (dude_train_receptors + dude_val_receptors)[fold]
                self.test_receptors = dude_test_receptors
            else:
                self.train_receptors = flatten(dude_train_receptors)
                self.val_receptors = flatten(dude_val_receptors)
                self.test_receptors = dude_test_receptors
        elif dataset == 'lit-pcba':
            if fold is not None:
                self.train_receptors = flatten([cluster for i, cluster in enumerate(lit_pcba_train_receptors + lit_pcba_val_receptors) if i != fold])
                self.val_receptors = (lit_pcba_train_receptors + lit_pcba_val_receptors)[fold]
                self.test_receptors = lit_pcba_test_receptors
            else:
                self.train_receptors = flatten(lit_pcba_train_receptors)
                self.val_receptors = flatten(lit_pcba_val_receptors)
                self.test_receptors = lit_pcba_test_receptors
        elif dataset == 'dude-lit-pcba':
            if fold is not None:
                dude_train_receptors_fold = flatten([cluster for i, cluster in enumerate(dude_train_receptors + dude_val_receptors) if i != fold])
                lit_pcba_train_receptors_fold = flatten([cluster for i, cluster in enumerate(lit_pcba_train_receptors + lit_pcba_val_receptors) if i != fold])
                self.train_receptors = dude_train_receptors_fold + lit_pcba_train_receptors_fold
                dude_val_receptors_fold = (dude_train_receptors + dude_val_receptors)[fold]
                lit_pcba_val_receptors_fold = (lit_pcba_train_receptors + lit_pcba_val_receptors)[fold]
                self.val_receptors = dude_val_receptors_fold + lit_pcba_val_receptors_fold
                self.test_receptors = dude_test_receptors + lit_pcba_test_receptors
            else:
                self.train_receptors = flatten(dude_train_receptors + lit_pcba_train_receptors)
                self.val_receptors = flatten(dude_val_receptors + lit_pcba_val_receptors)
                self.test_receptors = dude_test_receptors + lit_pcba_test_receptors
        elif dataset == 'dude-lit-pcba-xgboost-oversampled':
            self.train_receptors = flatten(dude_train_receptors + lit_pcba_train_receptors)
            self.val_receptors = flatten(dude_val_receptors + lit_pcba_val_receptors)
            self.test_receptors = dude_test_receptors + lit_pcba_test_receptors
        elif dataset == 'cache_1' or dataset == 'cache_1_part_2' or dataset == 'mcule':
            self.train_receptors = ['6DLO_A']
            self.val_receptors = ['6DLO_A']
            self.test_receptors = ['6DLO_A']
        elif dataset == 'dude-2':
            self.train_receptors = dude_excluded_test_receptors
            self.val_receptors = dude_excluded_test_receptors
            self.test_receptors = dude_excluded_test_receptors
        elif dataset == 'dude-mini':
            self.train_receptors = ['hivrt', 'ptn1']
            self.val_receptors = ['hivrt', 'ptn1']
            self.test_receptors = ['hivrt', 'ptn1']
        elif dataset == 'cache_2':
            self.train_receptors = ['5RLZ']
            self.val_receptors = ['5RLZ']
            self.test_receptors = ['5RLZ']            

        self.train_receptors = sorted(self.train_receptors)
        self.val_receptors = sorted(self.val_receptors)
        self.test_receptors = sorted(self.test_receptors)
        
        self.target = target
        
        # some failed jobs put inf values in plants_SCORE_NORM_CONTACT column. Remove these. 
        bad_idxs = []
        df.reset_index(inplace=True,drop=True)
        print('Removing bad ligands...')
        for i,val in enumerate(df['plants_SCORE_NORM_CONTACT']):
            if np.isposinf(val):
                bad_idxs.append(i)
        df = df.drop(bad_idxs)
        df.reset_index(inplace=True, drop=True)

        # Drop tools not specified in dock_tools argument
        all_tools = ['vina','ledock','rdock','plants','autodock']
        drop_tools_list = []
        for tool in all_tools:
            if tool not in dock_tools:
                drop_tools_list.append(tool)
        del_cols = []
        orig_len = len(df.columns)
        for i,col in enumerate(df.columns):
            for tool in drop_tools_list:
                if tool in col:
                    del_cols.append(col)
        del_cols += ['RMSD_avg', 'avg_RMSD_<_2']
        orig_len = len(df.columns)
        self.df = df.drop(del_cols, axis=1).sort_values('receptor')
        new_len = len(self.df.columns)
        print(f'{orig_len - new_len} columns dropped')


        if 'active' not in self.df.columns:
            self.label_available = False
        else:
            self.label_available = True

        print('Prepping ML inputs')
        self.train_df, self.val_df, self.test_df = self.gen_train_val_test_df()
        print('Generated train/val/test split')

        self.train_df = self.train_df.reset_index(drop=True)
        
        self.set_mode(mode=mode, target=target)

    def gen_train_val_test_df(self):
        '''
        Function that converts a df into a training df and a validation df based on receptor names.
        Inputs:
        df, training set targets, and list of validation set targets
        Returns:
        training df, validation df
        '''
        train_df = self.df[self.df['receptor'].isin(self.train_receptors)]
        val_df = self.df[self.df['receptor'].isin(self.val_receptors)]
        test_df = self.df[self.df['receptor'].isin(self.test_receptors)]
        return train_df, val_df, test_df

    def set_mode(self, mode=None, target=None):
        """Set the dataset to train/validation/test or a specific target"""
        self.mode = mode
        if mode == 'train':
            self.current_df = self.train_df
            self.receptors = self.train_receptors
        elif mode == 'val':
            self.current_df = self.val_df
            self.receptors = self.val_receptors
        elif mode == 'test':
            self.current_df = self.test_df
            self.receptors = self.test_receptors
        elif mode == 'single_target':
            self.target = target
            self.current_df = self.df.loc[self.df['receptor'] == target]
            self.receptors = [target]
        
        self.current_df = self.current_df.reset_index(drop=True)

    def __getitem__(self, i):
        if self.label_available:
            item = self.current_df.iloc[i, 3:-1]
            item = np.array(item).astype(np.float32)
            label = self.current_df.iloc[i, -1]
            label = np.array(label).astype(np.float32)
            return item, label

        else:
            item = self.current_df.iloc[i, 3:]
            item = np.array(item).astype(np.float32)
            
            return item

    def __len__(self):
        return len(self.current_df.index)


class OverSampler(Sampler):
    """OverSampler that balances data by oversampling active ligands and oversampling receptors with fewer ligands."""
    def __init__(self, dataset, shuffle=True):
        self.dataset = dataset
        self.shuffle = shuffle

    def __iter__(self):
        # Find max_n_items
        max_n_items = 0
        for receptor in self.dataset.current_df['receptor'].unique():
            n_items = sum((self.dataset.current_df['receptor'] == receptor) & (self.dataset.current_df['active'] == False))
            if n_items >= max_n_items:
                max_n_items = n_items

        # Oversample active ligands to make them equal with the inactive ligands
        # Oversample both actives and inactives to be equal with the receptor containing the most inactives
        self.all_indices = []
        for receptor in self.dataset.current_df['receptor'].unique():
            inactive_indices = self.dataset.current_df.index[(self.dataset.current_df['receptor'] == receptor) & (self.dataset.current_df['active'] == False)].tolist()
            active_indices = self.dataset.current_df.index[(self.dataset.current_df['receptor'] == receptor) & (self.dataset.current_df['active'] == True)].tolist()
            active_oversample_rate = len(inactive_indices) // len(active_indices)
            target_oversample_rate = max_n_items // len(inactive_indices)
            self.all_indices.extend(inactive_indices*target_oversample_rate)
            self.all_indices.extend(active_indices*active_oversample_rate*target_oversample_rate)
            
        if self.shuffle:
            self.all_indices = [self.all_indices[i] for i in torch.randperm(len(self.all_indices))]

        return iter(self.all_indices)

    def __len__(self):
        return len(self.all_indices)
        
        
def get_dataloader(batch_size=32, balance_data=True, shuffle=True, drop_last=False, **kwargs):

    dataset = DockingDataset(**kwargs)
    if balance_data:
        loader = DataLoader(
            dataset,
            batch_size=batch_size,
            sampler=OverSampler(dataset, shuffle=shuffle),
            num_workers=16,
            pin_memory=True,
            drop_last=drop_last
        )
    else:
        loader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=shuffle,
            # num_workers=16,
            pin_memory=True,
            drop_last=drop_last
        )
    return loader