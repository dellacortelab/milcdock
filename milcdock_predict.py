import argparse
from milcdock.test import predict_test, predict_new


parser = argparse.ArgumentParser()

# Test settings
parser.add_argument('--model_load_path', help='The base file path to load one or more models. Provide a single file to \
    specify a single model, provide a directory to specify an ensemble model', nargs='+', default=[])
parser.add_argument('--model_name_suffix', help='Any number of suffixes (ensemble, no_rdock, etc.) to be added to the'+
    ' model name in the results_df', nargs='+', default=[])
parser.add_argument('--report_metrics', help='Whether to report metrics for the model', action='store_true')
parser.add_argument('--device', help='The index for the desired CUDA device', default='0')
parser.add_argument('--batch_size', help='The batch size', default=128)
parser.add_argument('--ml_pkl_path', help='The path to the .pkl file containing inputs for the ML model', default=None)
parser.add_argument('--dataset', help='Alternative to ml_pkl_path. Provide a dataset name whose path is hard-coded in dataset.py', default=None)
parser.add_argument('--score_path', help='The path to the .pkl file containing inputs for the ML model', default='./data/4_milcdock_scores/scores.pkl')
parser.add_argument('--dock_tools', help='The docking tools to include in the ensemble', nargs='+', default=['vina', 'autodock', 'plants', 'rdock', 'ledock'])

if __name__ == "__main__":
    args = parser.parse_args()
    
    if args.report_metrics:
        # Evaluate model, write a pkl with the results
        predict_test(args)
    else:
        # Make predictions on unlabeled data, write a pkl with the results
        predict_new(args)
