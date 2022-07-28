import pandas as pd
import numpy as np
from torch.utils.data import Dataset, DataLoader
import os

class DockingDataset(Dataset):
    def __init__(self, data_path):
        if not os.path.exists(data_path):
            os.makedirs(os.path.dirname(data_path), exist_ok=True)
            if os.path.basename(data_path) == 'test_data.pkl':
                uri = https://byu.box.com/shared/static/90iyu2y8416pwyxqcec5il63w90pz45x
            elif os.path.basename(data_path) == 'test_data_2.pkl':
                uri = https://byu.box.com/shared/static/5u4vaqtregam8sfidjij4h7tm7x35j62
            os.system(f'curl -L {uri} --output {data_path}')
        self.df = pd.read_pickle(data_path)
        self.mode = 'full_data'
        self.target = None

    def set_mode(self, mode=None, target=None):
        """Set the dataset to train/validation/test or a specific target"""
        self.mode = mode
        if mode == 'full_data':
            self.current_df = self.df
        elif mode == 'single_target':
            self.target = target
            self.current_df = self.df.loc[self.df['receptor'] == target]
        
        self.current_df = self.current_df.reset_index(drop=True)

    def __getitem__(self, i):
        item = self.current_df.iloc[i, 3:-1]
        label = self.current_df.iloc[i, -1]
        item, label = np.array(item).astype(np.float32), np.array(label).astype(np.float32)
        return item, label

    def __len__(self):
        return len(self.current_df.index)

        
def get_dataloader(batch_size=32, **kwargs):

    dataset = DockingDataset(**kwargs)
    loader = DataLoader(
        dataset,
        batch_size=batch_size,
        num_workers=2,
        pin_memory=True
    )
    return loader
