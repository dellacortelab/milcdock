import numpy as np
import os
import torch
import torch.nn as nn
import torch.optim as optim
# import wandb
import copy
import collections

from .evaluate import calc_auc, calc_bedroc, calc_ef, binary_acc, evaluate
from .dataset import get_dataloader
from .model import DockingClassifier

def xavier_weight_gain_wrapper(weight_gain, bias):

    def xavier_weight_gain(m):
        if isinstance(m,nn.Linear):
            torch.nn.init.xavier_uniform_(m.weight, gain=weight_gain)
            torch.nn.init.constant_(m.bias, bias)

    return xavier_weight_gain

weights_init_dict = {'xavier_weight_gain':xavier_weight_gain_wrapper}

def save_model(model_save_path, model, optimizer, args, e, i):
    os.makedirs(os.path.dirname(model_save_path), exist_ok=True)
    save_dict = {}
    save_dict.update({'model':copy.deepcopy(model.state_dict())})
    save_dict.update({'optimizer':copy.deepcopy(optimizer.state_dict())})
    save_dict.update({'settings':args})
    save_dict.update({'best_epoch':e})
    save_dict.update({'best_iter':i})
    torch.save(save_dict, model_save_path)

def train(args):

    np.random.seed(0)
    torch.manual_seed(0)
    
    train_loader = get_dataloader(dataset=args.train_dataset, mode='train', fold=args.fold, batch_size=args.batch_size, balance_data=(args.no_balance_data == False), shuffle=True, drop_last=True, dock_tools=args.dock_tools)
    X = train_loader.dataset.train_df.iloc[:,3:-1]
    train_dataset_statistics = {'mean':torch.tensor(X.mean(axis=0).values.astype(np.float32)), 'std':torch.tensor(X.std(axis=0).values.astype(np.float32))}
    device = torch.device('cuda:' + args.device)
    model = DockingClassifier(depth=args.depth, layer_4_dim=args.layer_4_dim, activation=args.activation, dropout=args.dropout, dock_tools=args.dock_tools, train_dataset_statistics=train_dataset_statistics)
    model.to(device)
    
    if args.weights_init is not None:
        weights_init_fn = weights_init_dict[args.weights_init](weight_gain=args.weight_gain, bias=args.bias)
        model.apply(weights_init_fn)
    optimizer = optim.Adam(model.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay) #calculates gradient from loss function value

    if args.model_load_path != []:
        checkpoint = torch.load(args.model_load_path[0], map_location=device)
        model.load_state_dict(checkpoint['model'])
        optimizer.load_state_dict(checkpoint['optimizer'])
        # Don't load the train dataset statistics during training - we want the network to be able to use the
        # statistics for the new dataset if fine-tuning
        del checkpoint

    # We use the same train dataset statistics for each validation dataset, because the train dataset statistics 
    # correspond to the model, even if the validation dataset is different from the train dataset
    val_loaders = [get_dataloader(dataset=val_dataset, mode='val', fold=args.fold, batch_size=args.batch_size, balance_data=False, drop_last=False, shuffle=False, dock_tools=args.dock_tools) for val_dataset in args.val_dataset]

    # loss function - this one includes sigmoid activator on output automatically in training (not for predictions though)
    criterion = nn.BCEWithLogitsLoss()

    # os.system('git config --global --add safe.directory /code')
    # os.environ['WANDB_API_KEY'] = '<API_KEY>'
    # wandb.init(dir='./results', project='consensus-docking', entity='connormorris',name=args.run_name)
    # # Test positive and negative samples, receptor samples, make sure equal
    # # Log the arguments
    # config = wandb.config
    # arguments = [a for a in dir(args) if not a.startswith('_')]
    # for argument in arguments:
    #     config[argument] = getattr(args, argument)
        
    # wandb.watch(model)

    print('Beginning training')
    best_bedroc = -1
    best_epoch = -1
    # Only track the last 5k predictions, otherwise we get memory error
    train_pred = collections.deque(maxlen=5000)
    train_label = collections.deque(maxlen=5000)
    train_losses = collections.deque(maxlen=5000)
    train_accs = collections.deque(maxlen=5000)
    log_cnt = 0
    val_cnt = 0
    current_model_save_path = os.path.splitext(args.model_save_path)[0] + '_current.pt'

    for e in range(args.n_epochs):
        model.train()
        for i, (X_batch, y_batch) in enumerate(train_loader):
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            optimizer.zero_grad()
            y_pred = model(X_batch)
            loss = criterion(y_pred.squeeze(), y_batch)
            acc = binary_acc(y_pred.squeeze(), y_batch)

            loss.backward()
            optimizer.step()

            train_losses.append(loss.item())
            train_accs.append(acc.item())

            y_fornp = torch.sigmoid(y_pred)
            train_pred.extend(y_fornp.squeeze().tolist())
            train_label.extend(y_batch.tolist())
            if i % 10000 == 0:
                print(i)

            # Evaluate 5 times per epoch. The i != 0 is necessary, otherwise you log two batches in a row (at beginning and end)
            if i % (len(train_loader) // 5) == 0:# and i != 0:
                train_loss = np.mean(train_losses)
                train_acc = np.mean(train_accs)
                train_auc = calc_auc(np.array(train_pred), np.array(train_label))
                # print(len(train_pred), len(train_label))
                # import pdb; pdb.set_trace()
                train_bedroc = calc_bedroc(np.array(train_pred), np.array(train_label))
                _,train_ef = calc_ef(np.array(train_pred), np.array(train_label), return_fraction=True)

                # wandb.log({'train_acc':train_acc,'train_loss':train_loss, 'train_auc':train_auc,'train_bedroc':train_bedroc,'train_ef':train_ef}, step=log_cnt)
                for j, (dataset_name, val_loader) in enumerate(zip(args.val_dataset, val_loaders)):
                    save_model(current_model_save_path, model, optimizer, args, e, i)
                    eval_output = evaluate(model_save_path=current_model_save_path, dataloader=val_loader, batch_size=args.batch_size, device=device)
                    val_loss, val_acc = eval_output['loss'], eval_output['acc']
                    all_t_val_auc_mean = np.mean(list(eval_output['all_t_auc'].values()))
                    all_t_val_bedroc_mean = np.mean(list(eval_output['all_t_bedroc'].values()))
                    all_t_val_ef_mean = np.mean(list(eval_output['all_t_ef'].values()))

                    # wandb.log({dataset_name + '_val_acc':val_acc, dataset_name + '_val_loss':val_loss, \
                    #     dataset_name + '_target_val_auc':all_t_val_auc_mean, dataset_name + '_target_val_bedroc':all_t_val_bedroc_mean, dataset_name + '_target_val_ef':all_t_val_ef_mean}, step=log_cnt)

                    # val_dataset 0 is the dataset used for model performance criterion
                    if j == 0:

                        print(f'Epoch {e+0:03}: | Loss: {train_loss:.3f} | Acc: {train_acc:.3f} | Val Loss: {val_loss:.3f} | Val Acc: {val_acc:.3f} | Val Bedroc: {all_t_val_bedroc_mean:.3f}')

                        if (args.model_save_path is not None) and (all_t_val_bedroc_mean > best_bedroc):
                            best_bedroc = all_t_val_bedroc_mean
                            best_epoch = e
                            print(f'New best epoch: {best_epoch}')
                            print(f'New best bedroc: {best_bedroc}')

                            os.rename(current_model_save_path, args.model_save_path)

                            onnx_path = args.model_save_path[:args.model_save_path.rfind('.')] + '.onnx'
                            torch.onnx.export(copy.deepcopy(model), X_batch, onnx_path)
                            # wandb.save(onnx_path)
                            val_cnt = 0

                # wandb.log({}, step=log_cnt, commit=True)

                log_cnt += 1
                val_cnt += 1

                if val_cnt >= 20:
                    return