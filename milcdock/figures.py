from matplotlib import pyplot as plt
import numpy as np
from sklearn.metrics import precision_score, recall_score

def get_prec_rec(scores, labels):
    precisions = []
    recalls = []
    sort_idx = np.argsort(scores)
    scores, labels = scores[sort_idx], labels[sort_idx]
    n_thresholds = 200
    n_scores = len(scores)
    interval = n_scores // (n_thresholds - 1)
    thresholds = scores[::interval]
    for threshold in thresholds:
        pred = scores >= threshold
        precisions.append(precision_score(labels, pred))
        recalls.append(recall_score(labels, pred))
    thresholds = np.append(thresholds, [1.])
    precisions.append(precision_score(labels, scores >= max(scores)))
    recalls.append(0.)
    return thresholds, precisions, recalls

def make_hist(preds, labels, receptors):
    preds, labels = np.array(preds), np.array(labels)

    model_name = 'dude-lit-pcba_ensemble'
    target_list = ['reni', 'igf1r', 'ALDH1', 'GBA']

    # results_df, plot_dict = load_files(plot_dict_path, dataset_scores_path)

    fig, axs = plt.subplots(nrows=2, ncols=len(target_list), figsize=(22.5, 7))

    colors = ['blue', 'orange']

    for i, targ in enumerate(target_list):
            
        targ_preds, targ_labels = preds[receptors == targ], labels[receptors == targ]
        
        label_1 = 'Decoys'
        label_2 = 'Actives'
        label_3 = 'Precision'
        label_4 = 'Recall'
            
        num_bin = 30
        bin_lims = np.linspace(0,1,num_bin+1)
        bin_centers = 0.5*(bin_lims[:-1]+bin_lims[1:])
        bin_widths = (bin_lims[1:]-bin_lims[:-1])/2
        true_neg_scores = targ_preds[targ_labels == False]
        true_pos_scores = targ_preds[targ_labels == True]
        neg_hist, _ = np.histogram(true_neg_scores, bins=bin_lims)
        pos_hist, _ = np.histogram(true_pos_scores, bins=bin_lims)

        axs[0, i].bar(bin_centers - bin_widths/2, neg_hist, width=bin_widths, align='center', color='C0', label=label_1)
        ax2 = axs[0, i].twinx()
        ax2.bar(bin_centers + bin_widths/2, pos_hist, width=bin_widths, align='center', color='C1', label=label_2)
        axs[0, i].set_xlabel('MLC Score')

        axs[0, i].set_title(targ)
        
        thresholds, precisions, recalls = get_prec_rec(scores=targ_preds, labels=targ_labels)
        axs[1, i].plot(thresholds, precisions, color='C2', label=label_3)
        axs[1, i].plot(thresholds, recalls, color='C3', label=label_4)
        axs[1, i].set_xlabel('Score Threshold')
        
        axs[0, i].spines['top'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        
        ax2.spines['right'].set_color('C1')
        ax2.spines['left'].set_color('C0')
        axs[0, i].tick_params(axis='y', colors='C0')
        ax2.tick_params(axis='y', colors='C1')

        axs[1, i].spines['top'].set_visible(False)
        axs[1, i].spines['right'].set_visible(False)
        
        if i == 2:
            axs[1,i].legend()

    fig.tight_layout(pad=1)
    plt.subplots_adjust(top=0.85)
    plt.show()
