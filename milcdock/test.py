import os
import torch
from torch.nn import functional as F
import numpy as np
import pandas as pd

from .dataset import DockingDataset, get_dataloader
from .model import DockingClassifier, get_saved_model
from .evaluate import calc_auc, calc_bedroc, calc_ef, binary_acc
    


def score_new(model, dataloader, device):
    """Get model scores for ligands in a dataset without active/inactive labels
    Args:
        model (torch.nn.Module): The model to score
        dataloader (torch.utils.data.DataLoader): The dataloader for the dataset
        device (torch.device): The device to use for scoring
    Returns:
        (dict): A dictionary containing the scores and labels
    """
    model.eval()
    with torch.no_grad():
        pred = []

        for X_batch in dataloader:
            X_batch = X_batch.to(device)
            y_pred = model(X_batch)
            if type(model) == DockingClassifier:
                y_pred = torch.sigmoid(y_pred)
            y_pred = y_pred.squeeze().tolist()
            pred.extend(y_pred)

    model.train()
    return {'scores':pred}


def score(model, dataloader, device):
    """Score a model on a dataset
    Args:
        model (torch.nn.Module): The model to score
        dataloader (torch.utils.data.DataLoader): The dataloader for the dataset
        device (torch.device): The device to use for scoring
    Returns:
        (dict): A dictionary containing the scores and labels
    """
    model.eval()
    with torch.no_grad():
        pred = []
        label = []

        for X_batch, y_batch in dataloader:
            label.extend(y_batch.tolist())
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            y_pred = model(X_batch)
            if type(model) == DockingClassifier:
                y_pred = torch.sigmoid(y_pred)
            y_pred = y_pred.squeeze().tolist()
            pred.extend(y_pred)

    model.train()
    return {'scores':pred, 'labels':label}

def report_metrics(model_name, results_df):
    """Report metrics for a model
    Args:
        model_name (str): The name of the model
        results_df (pd.DataFrame): The dataframe containing the results
    """
    if 'active' not in results_df.columns:
        print("No label column found in results_df. Skipping metrics calculation.")
        return

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
        
    print(f'BEDROC: {all_t_bedroc}', f'AUC: {all_t_auc}', f'EF: {all_t_ef}', f'Accuracy: {all_t_acc}', f'Loss: {all_t_loss}', sep='\n')


def predict_new(args):
    """Make predictions on ligands in a dataset without active/inactive labels
    Args:
        args (argparse.Namespace): Arguments from command line
    """
    dataset = DockingDataset(dataset_path=args.ml_pkl_path, mode='test')
    results_df = dataset.current_df
    results_df = results_df.loc[:,results_df.columns.isin(['name', 'receptor', 'vina_top_score', 'autodock_top_score', 'plants_top_score', 'rdock_top_score', 'ledock_top_score', 'active'])]

    device = torch.device('cuda:' + args.device)
    test_loader = get_dataloader(dataset_path=args.ml_pkl_path, balance_data=False, shuffle=False, drop_last=False, batch_size=args.batch_size, mode='test', dock_tools=args.dock_tools)
    
    for model_load_path in args.model_load_path:
        model_suffix = '_' + '_'.join(args.model_name_suffix) if len(args.model_name_suffix) > 0 else ''
        base_model_name = os.path.basename(os.path.splitext(model_load_path)[0])
        model_name = base_model_name + model_suffix
        model = get_saved_model(model_path=model_load_path, device=device)
        eval_output = score_new(model=model, dataloader=test_loader, device=device)
        results_df = pd.concat((results_df.reset_index(drop=True), pd.Series(eval_output['scores'], name=model_name)), axis=1)

    results_df.to_pickle(args.score_path)

def predict_test(args):
    """Make predictions on ligands in a dataset with active/inactive labels
    Args:
        args (argparse.Namespace): Arguments from command line
    """
    dataset = DockingDataset(dataset_path=args.ml_pkl_path, dataset=args.dataset, mode='test')
    results_df = dataset.current_df
    results_df = results_df.loc[:,results_df.columns.isin(['name', 'receptor', 'vina_top_score', 'autodock_top_score', 'plants_top_score', 'rdock_top_score', 'ledock_top_score', 'active'])]

    device = torch.device('cuda:' + args.device)
    test_loader = get_dataloader(dataset_path=args.ml_pkl_path, dataset=args.dataset, balance_data=False, shuffle=False, drop_last=False, batch_size=args.batch_size, mode='test', dock_tools=args.dock_tools)
    
    for model_load_path in args.model_load_path:
        model_suffix = '_' + '_'.join(args.model_name_suffix) if len(args.model_name_suffix) > 0 else ''
        base_model_name = os.path.basename(os.path.splitext(model_load_path)[0])
        model_name = base_model_name + model_suffix
        model = get_saved_model(model_path=model_load_path, device=device)
        eval_output = score(model=model, dataloader=test_loader, device=device)
        results_df = pd.concat((results_df.reset_index(drop=True), pd.Series(eval_output['scores'], name=model_name)), axis=1)
        report_metrics(model_name=model_name, results_df=results_df)

    results_df.to_pickle(args.score_path)