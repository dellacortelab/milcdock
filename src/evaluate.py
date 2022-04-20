import torch
from torch import nn
from torch.utils.data import DataLoader
import os
import numpy as np
import warnings
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

def evaluate(model, dataloader, device, batch_size=32):

    prev_mode = dataloader.dataset.mode
    prev_target = dataloader.dataset.target

    all_preds = []
    all_labels = []
    model.eval() #use to make predictions
    with torch.no_grad():
        all_t_auc = dict()
        all_t_bedroc = dict()
        all_t_ef = dict()

        for target in dataloader.dataset.df['receptor'].unique():
            dataloader.dataset.set_mode(mode='single_target', target=target)
            # Reset dataloader iterator with new sampler
            dataloader = DataLoader(
                dataloader.dataset,
                batch_size=batch_size,
                pin_memory=True,
                num_workers=2
            )
            pred_vals = []
            label_vals = []
            for X_batch, y_batch in dataloader:
                
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                y_pred = model(X_batch)
                y_pred = y_pred.cpu().detach().numpy()[:,0]
                y_batch = y_batch.cpu().numpy()
                pred_vals.extend(y_pred)
                label_vals.extend(y_batch)


            pred_vals = np.array(pred_vals)
            label_vals = np.array(label_vals)
            t_auc = calc_auc(pred_vals, label_vals)
            all_t_auc.update({target:t_auc})
            t_bedroc = calc_bedroc(pred_vals, label_vals)
            all_t_bedroc.update({target:t_bedroc})
            t_ef = calc_ef(pred_vals, label_vals)
            all_t_ef.update({target:t_ef})
            print("Finished evaluating: ", target)

        all_preds.extend(pred_vals)
        all_labels.extend(label_vals)
        
    dataloader.dataset.set_mode(mode=prev_mode, target=prev_target)

    return {'all_preds':all_preds, 'all_labels':all_labels, 'all_t_auc':all_t_auc, 'all_t_bedroc':all_t_bedroc, 'all_t_ef':all_t_ef}


def calc_ef(y_pred, y_label, percent=.01):
    '''Returns enrichment factor based on predicted probability of being active'''
    num_active = len(np.where(y_label == 1)[0])
    tot_molecules = len(y_label)
    sort_idxs = np.argsort(y_pred)[::-1] #in descending order
    sorted_label = y_label[sort_idxs]
    len_top_percent = int(percent*tot_molecules)
    num_top_actives = len(np.where(sorted_label[:len_top_percent]==1)[0])
    ef = (num_top_actives/len_top_percent)/(num_active/tot_molecules)
    return ef
    
def calc_auc(y_pred, y_label, plot_roc=False):
    fpr, tpr, _ = roc_curve(y_label,y_pred)
    roc_auc = auc(fpr,tpr)
    if plot_roc:
        plt.plot(fpr,tpr,label=f"ROC curve, area = {roc_auc:.2f}")
        plt.legend()
        plt.plot([0,1],[0,1],color='r',linestyle='--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Curve')
        plt.show()
    return roc_auc

def calc_bedroc(y_pred, y_label, alpha=160.9):
    '''
    Calculates and returns BEDROC metric from labels and predicted probabilities.
    Formula comes from: https://pubs.acs.org/doi/full/10.1021/ci600426e
    '''
    assert isinstance(y_pred, np.ndarray)
    assert isinstance(y_label, np.ndarray)
    n = len(np.where(y_label == 1)[0]) # # of actives
    N = len(y_label) # total # of molecules
    sort_idxs = np.argsort(y_pred)[::-1] #in descending order
    sorted_label = y_label[sort_idxs]
    rie_num_sum = []
    for i in range(N):
        if sorted_label[i] == 1:
            rie_num_sum.append(np.exp(-alpha*(i+1)/N))
    rie_num = np.sum(rie_num_sum)
    rie_denom = n/N*((1-np.exp(-alpha))/(np.exp(alpha/N)-1))
    rie = rie_num / rie_denom
    rie_min = (1-np.exp(alpha*n/N))/(n/N*(1-np.exp(alpha)))
    rie_max = (1-np.exp(-alpha*n/N))/(n/N*(1-np.exp(-alpha)))
    bedroc = (rie - rie_min) / (rie_max - rie_min)
    return bedroc