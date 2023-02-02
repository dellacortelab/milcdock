import os
import pandas as pd
from sklearn.cluster import SpectralClustering
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_pickle('./alignment_matrix.pkl')

# df = df.drop(columns=['mk14', 'hmdh', 'vgfr2', 'pur2', 'met', 'kpcb', 'igf1r', 'egfr', 'mapk2', 'plk1', 'csf1r', 'reni', 'src', 'fak1', 'lck', 'glcm', 'jak2', 'adrb2', 'fgfr1', 'tgfr1', 'wee1', 'GBA', 'ALDH1', 'IDH1', 'ADRB2'])
# df = df.drop(index=['mk14', 'hmdh', 'vgfr2', 'pur2', 'met', 'kpcb', 'igf1r', 'egfr', 'mapk2', 'plk1', 'csf1r', 'reni', 'src', 'fak1', 'lck', 'glcm', 'jak2', 'adrb2', 'fgfr1', 'tgfr1', 'wee1', 'GBA', 'ALDH1', 'IDH1', 'ADRB2'])
for i in range(len(df)):
    for j in range(len(df)):
        if j < i:
            df.iloc[j,i] = df.iloc[i,j]

for i in range(len(df)):
    df.iloc[i,i] = 1000

score_mat = df.to_numpy(dtype=float)
delta = 12 #this has to be appropriately weighted
test = np.exp(- score_mat ** 2 / (2. * delta ** 2))
plt.imshow(test)
plt.show()
clustering = SpectralClustering(n_clusters=25, assign_labels='discretize', random_state=0).fit(test)
plt.hist(clustering.labels_)
# clustering.labels_
# train_idxs = np.where((clustering.labels_==0) | (clustering.labels_==5) )
# train_targets = np.array(df.columns)[train_idxs]
# valid_idxs = np.where((clustering.labels_==2) | (clustering.labels_==3) | (clustering.labels_==7) | (clustering.labels_==8) | (clustering.labels_==9))
# valid_targets = np.array(df.columns)[valid_idxs]
# test_idxs = np.where((clustering.labels_==1) | (clustering.labels_==4) | (clustering.labels_==6))
# test_targets = np.array(df.columns)[test_idxs]


for i in range(30):
    print(i)
    i_idxs = np.where(clustering.labels_==i)
    i_targets = np.array(df.columns)[i_idxs]
    print(i_targets)
# print(train_targets)
# print(valid_targets)
# print(test_targets)

# 0                                                                                                                                                                                                
# ['hdac2' 'cp3a4' 'ada17' 'pygm' 'mmp13' 'lkha4' 'bace1' 'dpp4'                                                                                                                             
#  'fnta' 'pde5a' 'hdac8' 'FEN1']                                                                                                                                                                  

# 2                                                                                                                                                                                                
# ['hivrt' 'pur2' 'pnph' 'gria2' 'hxk4' 'fpps' 'grik1' 'PKM2']                                                                                                                                     

# 3                                                                                                                                                                                                
# ['pparg' 'PPARG']                                                                                                                                                                            
# 11
# ['fa10' 'try1' 'urok']

# 12
# ['ptn1' 'pa2ga']
# 15
# ['fa7' 'tryb1' 'thrb']
# 17
# ['fabp4' 'ital' 'def' 'kif11' 'KAT2A' 'TP53']

# 24
# ['hivint' 'hivpr' 'kith']
# 27
# ['parp1' 'xiap']
# 29
# ['casp3' 'ampc' 'hs90a' 'cah2']



                                                                                                                                                                                                                                                                                                                                      
# 8                                                                                                                                                                                                
# ['ppard' 'ppara']      
# 14
# ['andr' 'prgr' 'gcr']       
# 19
# ['fkb1a' 'MTORC1']


# 26
# ['esr1' 'ESR1_ago' 'ESR1_ant']
# 18
# ['thb' 'esr2' 'rxra' 'mcr' 'VDR']
# 5                                                                                                                                                                                                
# ['OPRK1']     
# 20
# ['adrb1' 'cxcr4' 'drd3']

# 7                                                                                                                                                                                                
# ['mk14' 'cdk2' 'mk10' 'lck']
# 16
# ['rock1' 'akt2' 'akt1'] 
# 21
# ['abl1' 'mk01' 'MAPK1']





                                                                                                                                                                
# 1                                                                                                                                                                                                
# ['igf1r' 'wee1']                                                                                                                                                                                
# 4                                                                                                                                                                                                
# ['vgfr2' 'kpcb' 'fgfr1']          
# 6                                                                                                                                                                                                
# ['hmdh' 'aces' 'nram' 'reni' 'glcm' 'GBA' 'ALDH1' 'IDH1']                                                                                                                                                                                                                                                                                                                  
# 10                                                                                                                                                                                               
# ['adrb2' 'ADRB2']      
# 13
# ['egfr' 'jak2' 'tgfr1']
# 23
# ['met' 'kit' 'mapk2' 'src' 'fak1']                                                                                                                                                                          
# 9                                                                                                                                                                                                
# ['plk1' 'csf1r' 'braf']