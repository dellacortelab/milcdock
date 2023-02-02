import torch
from torch import nn
import numpy as np
import warnings
from sklearn.metrics import roc_curve, auc
from torch.utils.data import DataLoader


from .model import get_saved_model, DockingClassifier

warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

def evaluate(model_save_path, dataloader, batch_size, device):
    model = get_saved_model(model_path=model_save_path, device=device)
    
    criterion = nn.BCEWithLogitsLoss()
    model.eval() #use to make predictions
    with torch.no_grad():
        all_t_auc = dict()
        all_t_bedroc = dict()
        all_t_ef = dict()

        prev_mode = dataloader.dataset.mode
        prev_target = dataloader.dataset.target
        loss = 0
        acc = 0
        num_examples = 0
        for target in dataloader.dataset.receptors:
            if target not in dataloader.dataset.df['receptor'].unique():
                raise ValueError(f"Receptor {target} not in dataset")
            dataloader.dataset.set_mode(mode='single_target', target=target)
            # Reset dataloader iterator with new sampler
            dataloader = DataLoader(
                dataloader.dataset,
                batch_size=batch_size,
                pin_memory=True,
                num_workers=16
            )
            pred_vals = []
            label_vals = []
            for X_batch, y_batch in dataloader:
                
                num_examples += len(X_batch)
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                y_pred = model(X_batch)

                if type(model) == (DockingClassifier):
                    y_pred = torch.sigmoid(y_pred)
                
                batch_loss = criterion(y_pred.squeeze(-1), y_batch)
                batch_acc = binary_acc(y_pred.squeeze(-1), y_batch)
                loss += batch_loss.item() * len(X_batch)
                acc += batch_acc.item() * len(X_batch)
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
            _,t_ef = calc_ef(pred_vals, label_vals, return_fraction=True)
            all_t_ef.update({target:t_ef})

        if num_examples == 0:
            return {'loss':0., 'acc':0., 'all_t_auc':{'tmp':0.}, 'all_t_bedroc':{'tmp':0.}, 'all_t_ef':{'tmp':0.}}

        loss = loss / num_examples
        acc = acc / num_examples
        
        dataloader.dataset.set_mode(mode=prev_mode, target=prev_target)
    
    model.train()


    return {'loss':loss, 'acc':acc, 'all_t_auc':all_t_auc, 'all_t_bedroc':all_t_bedroc, 'all_t_ef':all_t_ef}

# Metrics
def binary_acc(y_pred, y_test):
    y_pred_tag = torch.round(y_pred)

    correct_results_sum = (y_pred_tag == y_test).sum().float()
    acc = correct_results_sum/y_test.shape[0]
    acc = torch.round(acc * 100)
    
    return acc

def calc_ef(y_pred, y_label, percent=.01,return_fraction=False):
    '''Returns enrichment factor based on predicted probability of being active'''
    num_active = len(np.where(y_label == 1)[0])
    tot_molecules = len(y_label)
    sort_idxs = np.argsort(y_pred)[::-1] #in descending order
    sorted_label = y_label[sort_idxs]
    len_top_percent = int(percent*tot_molecules)
    if num_active == 0 or len_top_percent == 0:
        if return_fraction:
            return np.inf, np.inf
        else:
            return np.inf
    num_top_actives = len(np.where(sorted_label[:len_top_percent]==1)[0])
    ef = (num_top_actives/len_top_percent)/(num_active/tot_molecules)

    if return_fraction:
        #calc theoretical max
        if num_active >= len_top_percent:
            theoretical_max = (1)/(num_active/tot_molecules)
        else:
            theoretical_max = (num_active/len_top_percent)/(num_active/tot_molecules)
        fractional_ef = ef / theoretical_max
        return ef, fractional_ef
    else:
        return ef
    
def calc_auc(y_pred, y_label):
    fpr, tpr, _ = roc_curve(y_label,y_pred)
    roc_auc = auc(fpr,tpr)
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
    if len(rie_num_sum) != n:
        print('Error. Wrong sum length.')
        return
    rie_num = np.sum(rie_num_sum)
    rie_denom = n/N*((1-np.exp(-alpha))/(np.exp(alpha/N)-1))
    rie = rie_num / rie_denom
    rie_min = (1-np.exp(alpha*n/N))/(n/N*(1-np.exp(alpha)))
    rie_max = (1-np.exp(-alpha*n/N))/(n/N*(1-np.exp(-alpha)))
    bedroc = (rie - rie_min) / (rie_max - rie_min)
    return bedroc