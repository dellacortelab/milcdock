import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import precision_score, recall_score
import os
import torch
import pandas as pd
import wandb
import torch.nn.functional as F
import numpy as np

from .dataset import dataset_dict
from .evaluate import calc_auc, calc_bedroc, calc_ef, binary_acc


def make_fig(args):
    results_dfs = dict()
    if args.no_wandb:
        for dataset in args.dataset:
            base_path = os.path.splitext(dataset_dict[dataset])[0]
            mode = args.eval_targets
            pkl_path = base_path + f'-no_wandb-{mode}-results.pkl'
            results_df = pd.read_pickle(pkl_path)
            results_dfs.update({dataset:results_df})
    else:
        for dataset in args.dataset:
            base_path = os.path.splitext(dataset_dict[dataset])[0]
            mode = args.eval_targets
            pkl_path = base_path + f'-{mode}-results.pkl'
            os.system('git config --global --add safe.directory /code')
            os.environ['WANDB_API_KEY'] = '42d41437d35eaa5242dbd4a5b77f64ae09a00e90'
            run = wandb.init(dir='./results/wandb/', project='consensus-docking', entity='connormorris', resume=True, name='model_test_pkl')
            artifact = run.use_artifact(artifact_or_name=os.path.basename(pkl_path) + ':latest', type='test_results')
            artifact.download(os.path.dirname(pkl_path))
            results_df = pd.read_pickle(pkl_path)
            results_dfs.update({dataset:results_df})
            
    exp_dict[args.test_fig](args, results_dfs)

def make_fig_1(args, results_dfs):
    """Log scores for all metrics for each model on each dataset"""
    for model_name in args.model_name:
        eval_output = get_metrics(model_name, results_dfs[args.dataset[0]])

        if 'vina' in model_name or 'rdock' in model_name or 'ledock' in model_name or 'plants' in model_name or 'autodock' in model_name:
            print(
                "Model: " + model_name, \
                "Dataset: " + args.dataset[0], \
                "Mean Per-target AUC: " + str(np.mean(list(eval_output['all_t_auc'].values()))), \
                "Mean Per-target BEDROC: " + str(np.mean(list(eval_output['all_t_bedroc'].values()))), \
                "Mean Per-target EF: " + str(np.mean(list(eval_output['all_t_ef'].values())))
            )
        else:
            print(
                "Model: " + model_name, \
                "Dataset: " + args.dataset[0], \
                "Mean Per-target loss: " + str(np.mean(list(eval_output['all_t_loss'].values()))), \
                "Mean Per-target acc: " + str(np.mean(list(eval_output['all_t_acc'].values()))), \
                "Mean Per-target AUC: " + str(np.mean(list(eval_output['all_t_auc'].values()))), \
                "Mean Per-target BEDROC: " + str(np.mean(list(eval_output['all_t_bedroc'].values()))), \
                "Mean Per-target EF: " + str(np.mean(list(eval_output['all_t_ef'].values())))
            )

def make_fig_2():
    pass


def make_fig_3(args, results_dfs):
    """Log per-target bedroc scores for each model on each dataset"""
    results = pd.DataFrame()
    for model_name in args.model_name:
        print(model_name)
        eval_output = get_metrics(model_name, results_dfs[args.dataset[0]])
        results = pd.concat((results, pd.Series(eval_output['all_t_auc'], name=model_name)), axis=1)
    csv_path = os.path.join(args.figs_dir, args.dataset[0] + ".csv")
    results.to_csv(path_or_buf=csv_path)



def make_fig_4(args, results_dfs):
    
    # fig, axs = plt.subplots(nrows=2, ncols=len(args.dataset), figsize=(22.5, 7))
    fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(22.5, 7))

    # for i, model_name in enumerate(args.model_name):
    for j, dataset in enumerate(results_dfs.keys()):
        if j == 0:
            label_1 = 'Decoys'
            label_2 = 'Actives'
            label_3 = 'Precision'
            label_4 = 'Recall'
        else:
            label_1 = None
            label_2 = None
            label_3 = None
            label_4 = None

        results_df = results_dfs[dataset]
        target_list = ['ALDH1', 'GBA']
        for i, targ in enumerate(target_list):
            scores, labels = np.array(results_df[results_df['receptor']==targ][args.model_name[0]]), np.array(results_df[results_df['receptor']==targ]['active'])
            # scores, labels = np.array(results_df[args.model_name[0]]), np.array(results_df['active'])

            num_bin = 30
            bin_lims = np.linspace(0,1,num_bin+1)
            bin_centers = 0.5*(bin_lims[:-1]+bin_lims[1:])
            bin_widths = (bin_lims[1:]-bin_lims[:-1])/2

            true_neg_scores = scores[labels == False]
            true_pos_scores = scores[labels == True]
            neg_hist, _ = np.histogram(true_neg_scores, bins=bin_lims)
            pos_hist, _ = np.histogram(true_pos_scores, bins=bin_lims)

            axs[0, i].bar(bin_centers - bin_widths/2, neg_hist, width=bin_widths, align='center', color='C0', label=label_1)
            ax2 = axs[0, i].twinx()
            ax2.bar(bin_centers + bin_widths/2, pos_hist, width=bin_widths, align='center', color='C1', label=label_2)
            axs[0, i].set_xlabel('Model Score')
            axs[0, i].set_ylabel('Number of Negative Examples')
            ax2.set_ylabel('Number of Positive Examples')

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
            axs[1, i].plot(thresholds, precisions, color='C2', label=label_3)
            axs[1, i].plot(thresholds, recalls, color='C3', label=label_4)
            axs[1, i].set_xlabel('Score Threshold')

    for ax, col in zip(axs[0], target_list):
        ax.set_title(col, fontsize=20)

    fig.tight_layout(pad=2)
    plt.subplots_adjust(top=0.85)
    fig.suptitle('Model Score Distributions', fontsize=24)

    fig_path = os.path.join(args.figs_dir, f'{args.model_name[0]}_double_hist.png')
    plt.savefig(fig_path)
    

def make_fig_5(model_name, dataset_name, eval_output):
    pass


exp_dict = {'fig_1':make_fig_1, 'fig_2':make_fig_2, 'fig_3':make_fig_3, 'fig_4':make_fig_4, 'fig_5':make_fig_5}


def get_metrics(model_name, results_df):
    all_t_bedroc = {}
    all_t_auc = {}
    all_t_ef = {}
    all_t_acc = {}
    all_t_loss = {}
    for target in results_df['receptor'].unique():
        t_df = results_df[results_df['receptor'] == target]
        scores = np.array(t_df[model_name])
        labels = np.array(t_df['active'])
        
        if model_name == 'vina_top_score' or model_name == 'rdock_top_score' or model_name == 'ledock_top_score' or model_name == 'plants_top_score' or model_name == 'autodock_top_score':
            scores = -scores
            
        if model_name == 'vina_top_score' or model_name == 'rdock_top_score' or model_name == 'ledock_top_score' or model_name == 'plants_top_score' or model_name == 'autodock_top_score' or 'naive_consensus' in model_name:
            t_acc = 0.
            t_loss = 0.
        else:
            t_acc = binary_acc(torch.Tensor(scores), torch.Tensor(labels)).item()
            t_loss = F.binary_cross_entropy(torch.Tensor(scores), torch.Tensor(labels)).item()
        t_bedroc = calc_bedroc(scores, labels)
        t_ef, t_ef_frac = calc_ef(scores, labels, return_fraction=True)
        t_auc = calc_auc(scores, labels)
        
        all_t_bedroc[target] = t_bedroc
        all_t_auc[target] = t_auc
        all_t_ef[target] = t_ef
        all_t_acc[target] = t_acc
        all_t_loss[target] = t_loss
        
    return {'all_t_bedroc':all_t_bedroc, 'all_t_auc':all_t_auc, 'all_t_ef':all_t_ef, 'all_t_acc':all_t_acc, 'all_t_loss':all_t_loss}
