import pandas as pd
import numpy as np
import torch
from tqdm import tqdm
from data_process import load_data_save, load_data_save2
import os

#-------------------------- Data load ----------------------------
def load_data2(config):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    load_data_save()
    load_data_save2()

    score_gene_num_vector = np.load(script_dir + "/data/score_mutation.npy")

    train = np.array(score_gene_num_vector[:, 1:], dtype=int)
    targets = np.array(score_gene_num_vector[:, 0], dtype=float)

    return train, targets