
import torch.nn as nn
import torch
import os

def get_model(model_path, device):

    if os.path.isdir(model_path):
        models = []
        for model_load_path in os.listdir(model_path):
            checkpoint = torch.load(os.path.join(model_path, model_load_path), map_location=device)
            model_settings = checkpoint['settings']
            model = DockingClassifier(depth=model_settings.depth, layer_4_dim=model_settings.layer_4_dim, activation=model_settings.activation, dropout=model_settings.dropout)
            model.load_state_dict(checkpoint['model'])
            model.to(device)

            models.append(model)
            
        model = EnsembleDockingClassifier(models)

    elif os.path.isfile(model_path):
        checkpoint = torch.load(model_path, map_location=device)
        model_settings = checkpoint['settings']
        model = DockingClassifier(depth=model_settings.depth, layer_4_dim=model_settings.layer_4_dim, activation=model_settings.activation, dropout=model_settings.dropout)
        model.load_state_dict(checkpoint['model'])
        model.to(device)
        
    else:
        model = DockingClassifier()

    return model

class EnsembleDockingClassifier(nn.Module):
    def __init__(self, models):
        super().__init__()
        self.models = nn.ModuleList(models)

    def forward(self, x):
        preds = []
        for model in self.models:
            pred = torch.sigmoid(model(x))
            preds.append(pred)
        return torch.mean(torch.stack(preds), dim=0)

class DockingClassifier(nn.Module):
    def __init__(self, output_dim=1, depth=90, layer_4_dim=64, width=64, activation='relu', dropout=.4, dock_tools=['vina', 'autodock', 'plants', 'ledock', 'rdock'], train_dataset_statistics={'mean':torch.zeros(52), 'std':torch.ones(52)}):
        super().__init__()

        vec_len_dict = {'vina':4, 'rdock':14, 'ledock':4, 'plants':11, 'autodock':9}
        n = len(dock_tools)
        num_combinations = (n*(n-1)) // 2
        input_vec_len = sum([vec_len_dict[tool] for tool in dock_tools]) + num_combinations

        self.dropout_prob = dropout
        
        if activation == 'relu':
            self.act_func = nn.ReLU()
        elif activation == 'tanh':
            self.act_func = nn.Tanh()
        elif activation == 'sigmoid':
            self.act_func = nn.Sigmoid()
        elif activation == 'selu':
            self.act_func = nn.SELU()
            
        self.layer_1 = nn.Sequential(
            nn.Linear(input_vec_len, width),
            self.act_func,
            nn.Dropout(p=self.dropout_prob)
        )
        self.hiddens = []
        for i in range(depth):
            in_dim = width
            out_dim = width
            if i == 2:
                out_dim = layer_4_dim
            if i == 3:
                in_dim = layer_4_dim
            self.hiddens.append(nn.Sequential(
                nn.Linear(in_dim, out_dim),
                self.act_func,
                nn.Dropout(p=self.dropout_prob)
            ))
            self.add_module('hidden_layer' + str(i), self.hiddens[-1])
        self.layer_out = nn.Sequential(
            nn.Linear(width, output_dim)
        )

        # train_dataset_statistics['mean'], train_dataset_statistics['std'] = torch.tensor(train_dataset_statistics['mean']), torch.tensor(train_dataset_statistics['std'])
        # self.train_dataset_statistics = train_dataset_statistics
        self.register_buffer('train_data_mean', train_dataset_statistics['mean'])
        self.register_buffer('train_data_std', train_dataset_statistics['std'])

    def forward(self, inputs):
        inputs = (inputs - self.train_data_mean) / self.train_data_std
        x = self.layer_1(inputs)
        for layer in self.hiddens:
            x = layer(x)
        x = self.layer_out(x)
        return x