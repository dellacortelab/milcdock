import torch

from src.model import get_model
from src.dataset import get_dataloader
from src.evaluate import evaluate
from src.figures import make_hist

#model_path = './data/models/ensemble_mlp'
model_path = './data/models/mlp.pt'
data_path = './data/test_data/test_data_2.pkl'

device = 'cuda' if torch.cuda.is_available() else 'cpu'
model = get_model(model_path=model_path, device=device)
dataloader = get_dataloader(data_path=data_path)
scores = evaluate(model, dataloader, device)
print("AUCs: ", scores['all_t_auc'])
print("BEDROCs: ", scores['all_t_bedroc'])
print("EFs: ", scores['all_t_ef'])

#make_hist(scores['all_preds'], scores['all_labels'], dataloader.dataset.df['receptor'])
