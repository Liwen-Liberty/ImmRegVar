import pandas as pd
import numpy as np
import os

def load_data_save():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    df = pd.read_csv(script_dir + "/data/Solid_tumor_CYT_score_df.csv")
    rf = pd.read_csv(script_dir + "/data/Solid_tumor_mutation_status_mat.csv")
    
    final=[]
    p_list = list(set(rf.columns) & set(df["PatientID"]))
    for i in range(len(df)):
        pair=[]
        if df.loc[i]["PatientID"] in p_list:
            PatientID=df.loc[i]["PatientID"]
            CYT_score=df.loc[i]["CYT_score"]
            PatientID_feature=rf.loc[:,PatientID].tolist() 
            pair.append(CYT_score)
            pair.extend(PatientID_feature)
            final.append(pair)
    
    final_array=np.array(final)
    np.save(script_dir + "/data/score_mutation.npy",final_array)

def load_data_save2():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    df = pd.read_csv(script_dir + "/data/Solid_tumor_CYT_score_df.csv")
    rf = pd.read_csv(script_dir + "/data/Solid_tumor_mutation_status_mat.csv")
    
    n_gene=691
    final=[]
    for i in range(len(df)):
        pair=[]
        CYT_score=df.loc[i]["CYT_score"]
        pair.append(CYT_score)
        for j in range(n_gene):
            pair.append(j)
        final.append(pair)
    
    
    final_array=np.array(final)
    np.save(script_dir + "/data/score_gene_num_vector.npy",final_array)

    
