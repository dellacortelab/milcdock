import argparse
import time
import numpy as np
import os
from milcdock.train import train

parser = argparse.ArgumentParser()

# Train/Test settings
parser.add_argument('--mode', help='Whether to train or test (default: %(default)s)', default='train', choices=['train', 'test', 'train_simple'])
parser.add_argument('--model_save_path', help='The file path for the saved model', default=None, required=True)
parser.add_argument('--model_load_path', help='The file path to load one or more models (only one allowed if --mode==train)', nargs='+', default=[])
parser.add_argument('--run_name', help='The run name for WandB (default: %(default)s)', default='tmp')
parser.add_argument('--model_id', type=int, default=None)
parser.add_argument('--device', help='The index for the desired CUDA device (default: %(default)s)', default='0')
parser.add_argument('--batch_size', help='The batch size (default: %(default)d)', type=int, default=32)
parser.add_argument('--n_epochs', help='The number of epochs to train for (default: %(default)d)', default=200, type=int)
parser.add_argument('--learning_rate', help='The learning rate (default: %(default)f)', default=.0001, type=float)
parser.add_argument('--weight_decay', help='The weight regularization parameter (default: %(default)f)', default=0.0, type=float)

# Data settings
parser.add_argument('--train_dataset', help='Which dataset to train on (default: %(default)s)', choices=['dude', 'lit-pcba', 'dude-lit-pcba', 'dude-mini', 'dude-orig', 'dude-orig-oversampled', 'dude-lit-pcba-xgboost-oversampled'], default='dude-lit-pcba')
parser.add_argument('--val_dataset', help='Which datasets to validate on. test_dataset 0 is the dataset used as the model performance criterion. (default: %(default)s)', nargs='+', choices=['dude', 'lit-pcba', 'dude-lit-pcba', 'dude-mini', 'dude-orig', 'dude-orig-oversampled'], default='dude-lit-pcba')
parser.add_argument('--fold', help='The fold to train on, if doing cross-validation or bagging', type=int, default=None)
parser.add_argument('--dock_tools', help='The docking tools to include in the ensemble (default: %(default)s)', nargs='+', default=['vina', 'autodock', 'plants', 'rdock', 'ledock'])
parser.add_argument('--no_balance_data', help='Whether to NOT rebalance train data', action='store_true')

# Model settings
parser.add_argument('--activation', help='The activation function for the network (default: %(default)s)', default='relu', choices=['relu', 'tanh', 'sigmoid', 'selu'])
parser.add_argument('--depth', help='Number of layers (default: %(default)d)', default=4, type=int)
parser.add_argument('--layer_4_dim', help='Width of the 3rd hidden layer (default: %(default)d)', default=64, type=int)
parser.add_argument('--weight_gain', help='The weight gain of each weight layer (default: %(default)f)', default=np.sqrt(2), type=float)
parser.add_argument('--bias', help='The bias of each weight layer (default: %(default)f)', default=0.01, type=float)
parser.add_argument('--dropout', help='The proportion of neurons to drop out during training (default: %(default)f)', default=.4, type=float)
parser.add_argument('--residual', help='Whether to include skip connections', action='store_true')
parser.add_argument('--weights_init', help='The type of weight initialization (default: %(default)s)', default='xavier_weight_gain', choices=['xavier_weight_gain', 'std_normal'])


def update_args(args):
    """Update arguments for a test_id corresponding to an ablation test.
    Defaults: activation: relu, weight_decay: 0., depth: 4, width: 64, dropout: .4, 
    """
    cnt = 0
    model_list = []
    for dropout in [.2, .4, .6]:
        for weight_decay in [0, 1e-6, 1e-5, 1e-4]:    
            for activation in ['relu', 'selu']:
                for layer_4_dim in [32, 64, 128, 256]:
                    for weight_gain in [1., np.sqrt(2)]:
                        model_list.append({'dropout':dropout, 'weight_decay':weight_decay, 'activation':activation, 'layer_4_dim':layer_4_dim, 'weight_gain':weight_gain})

                        if cnt == 90:
                            print(model_list[-1])
                        cnt += 1
    args.n_epochs = 40
    filename, ext = os.path.splitext(args.model_save_path)
    args.model_save_path = filename + '_' + str(args.model_id) + ext
    settings = model_list[args.model_id]
    print(settings)

    for characteristic, option in settings.items():
        setattr(args, characteristic, option)
        
    setattr(args, 'run_name', '-'.join([characteristic + ':' + str(option) for characteristic, option in settings.items()]))

    return args


if __name__ == "__main__":

    args = parser.parse_args()

    if args.model_id is not None:
        args = update_args(args)
    
    if args.run_name is None:
        args.run_name = os.path.basename(args.model_save_path)
        
    t = time.time()
    train(args)
    print('total time for ' + args.run_name + ' ' + str(time.time()-t))
