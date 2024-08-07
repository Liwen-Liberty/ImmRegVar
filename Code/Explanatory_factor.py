import torch
import torch.nn as nn
import numpy as np
from sklearn.metrics import r2_score
from Transformer_model import Model
from args import Config
import matplotlib.pyplot as plt
from load_data import load_data, load_data2
from tqdm import tqdm
import os
import pandas as pd

def explanatory_factor(config):
    # Load configuration
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Initialize the model
    model = Model(config)
    model.to(config.device)
    
    # Load saved model weights
    model.load_state_dict(torch.load(script_dir + '/Results/model_state_dict.pt'))
    
    # Prepare test data
    train_load, targets = load_data2(config)
    x_train_load = torch.from_numpy(train_load[:]).to(torch.long).to(config.device)
    y_train = torch.from_numpy(targets[:]).to(torch.float).to(config.device)
    
    # Set model to evaluation mode
    model.eval()
    
    # Prediction
    with torch.no_grad():
        predictions = []
        actuals = []
        adj_train = []
        for batch in range(0, x_train_load.shape[0], config.batch_size):
            train_pre, _, aw = model(x_train_load[batch:batch + config.batch_size])
            predictions.extend(train_pre.view(-1).cpu().numpy())
            actuals.extend(y_train[batch:batch + config.batch_size].cpu().numpy())
            adj_train.extend(_.cpu().numpy().tolist())
    
    
    # Initialize matrices for storing predictions
    predictions_hat_matrix = []
    
    # Set model to evaluation mode
    model.eval()
    
    # Iterate through each feature for perturbation analysis
    for i in tqdm(range(x_train_load.size()[1])):
        ones_vector = torch.ones(x_train_load.size()[1], dtype=torch.long).to(config.device)
        ones_vector[i] = 0
        x_train_new = x_train_load * ones_vector
    
        # Prediction with perturbed input
        with torch.no_grad():
            predictions_hat = []
            for batch in range(0, x_train_new.shape[0], config.batch_size):
                train_pre_hat, _hat, aw_hat = model(x_train_new[batch:batch + config.batch_size])
                predictions_hat.extend(train_pre_hat.view(-1).cpu().numpy())
    
        predictions_hat_matrix.append(predictions_hat)
    
    # Compute explanatory factor
    difference_matrix = [[hat - act for hat, act in zip(row, actuals)] for row in predictions_hat_matrix]
    difference_vector = [pre - act for pre, act in zip(predictions, actuals)]
    ratio_matrix = [[abs(hat) / abs(pre) for hat, pre in zip(row, difference_vector)] for row in difference_matrix]
    C_Delta_xi = np.array(list(zip(*ratio_matrix)))
    C_Delta_xi_log = np.log2(C_Delta_xi)
    
    # Load attention weights
    attention_weight = np.load(script_dir + "/Results/attention_weights.npy")
    attention_weight_torch = torch.from_numpy(attention_weight)
    A_sum = torch.sum(attention_weight_torch, dim=1)
    A_sum_np = A_sum.numpy()
    
    # Weighted C_Delta_xi_log
    weighted_C_Delta_xi_log = C_Delta_xi_log * A_sum_np
    
    # Compute average of log values
    C_Delta_xi_log_avg = np.nanmean(C_Delta_xi_log, axis=0)
    sorted_indices_C_log = np.argsort(C_Delta_xi_log_avg)[::-1]
    
    # Compute average of weighted log values
    weighted_log_avg = np.nanmean(weighted_C_Delta_xi_log, axis=0)
    sorted_indices_weight_log = np.argsort(weighted_log_avg)[::-1]
    
    return sorted_indices_C_log, sorted_indices_weight_log


config = Config()
unweighted_rank, weighted_rank = explanatory_factor(config)














