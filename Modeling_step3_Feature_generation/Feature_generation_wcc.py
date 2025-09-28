# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 11:14:02 2025

@author: PC
"""

import os
import numpy as np
import pandas as pd
from numba import njit
from joblib import Parallel, delayed
from tqdm.auto import tqdm
import warnings
warnings.filterwarnings('ignore')

# 1. Fast weight correlation
@njit(cache=True)
def numba_wcc(x, y, tr_width=10):
    n = len(x)
    max_lag = min(tr_width, n)
    weights = 1 - np.arange(max_lag) / tr_width
    ac1 = 0.0
    ac2 = 0.0
    cross = 0.0
    for lag in range(max_lag):
        x_lag = x[lag:]
        y_lag = y[lag:]
        x_lead = x[:n-lag]
        y_lead = y[:n-lag]        
        ac1 += weights[lag] * np.sum(x_lead * x_lag)
        ac2 += weights[lag] * np.sum(y_lead * y_lag)
        cross += weights[lag] * np.sum(x_lead * y_lag)    
    denominator = np.sqrt(ac1 * ac2)
    return cross / denominator if denominator > 1e-10 else np.nan

# 2. De noise - default is 0.001
def preprocess_profile(profile, impute_NA=True, noise_floor=0.001):
    # clean profile - each protein
    profile = profile.astype(np.float64).copy()
    profile[profile < noise_floor] = noise_floor    
    # if there is NA - imputate 0 or col_mean (here is 0 as example)
    if impute_NA:
        profile[np.isnan(profile)] = 0.0    
    return profile

# 3. Sort input protein matrix
def load_and_preprocess(tpc_path, pairs_path):
    tpc_chunks = pd.read_csv(tpc_path, sep='\t', chunksize=10000) # size = 10,000    
    first_chunk = next(tpc_chunks)    
    # sorted the column - make sure from 1 to 60
    numeric_cols = [col for col in first_chunk.columns if col != 'ID' and str(col).isdigit()]
    numeric_cols_sorted = sorted(numeric_cols, key=lambda x: int(x))
    sorted_cols = ['ID'] + numeric_cols_sorted
    # reload
    tpc_data = pd.concat([first_chunk[sorted_cols]] + [chunk[sorted_cols] for chunk in tpc_chunks])    
    # make sure the ID is "protein ID"
    tpc_data = tpc_data.set_index('ID')    
    # clean profile - each protein
    tpc_cleaned = np.array([
        preprocess_profile(row.values) 
        for _, row in tqdm(tpc_data.iterrows(), desc="Preprocessing profiles")])    
    # build dict for ID
    id_dict = {pid: i for i, pid in enumerate(tpc_data.index)}
    # load the pair matrix
    pairs_df = pd.read_csv(pairs_path, sep='\t')    
    return tpc_cleaned, id_dict, pairs_df

# 3. Parallel
def calculate_wcc_for_pair(pair, tpc_data, id_dict, tr_width=10):
    try:
        idx1 = id_dict[pair[0]]
        idx2 = id_dict[pair[1]]
        return {
            'Protein_A': pair[0], # column1 name of pair matrix
            'Protein_B': pair[1], # column1 name of pair matrix
            'WCC_score': numba_wcc(tpc_data[idx1], tpc_data[idx2], tr_width)}
    except KeyError:
        return None

# 4. Chunk
def compute_wcc_parallel(tpc_data, id_dict, pairs_df, chunk_size=10000, n_jobs=-1):
    # pair matrix
    pairs = pairs_df[['protein_A', 'protein_B']].values    
    n_chunks = len(pairs) // chunk_size + 1  # chunk = 10,000
    chunks = np.array_split(pairs, n_chunks)    
    # parallel
    results = Parallel(n_jobs=n_jobs)(
        delayed(lambda chunk: [
            calculate_wcc_for_pair(pair, tpc_data, id_dict) 
            for pair in chunk])(chunk)
        for chunk in tqdm(chunks, desc="Processing chunks"))   
    return [item for sublist in results for item in sublist if item is not None] # comb the results and drop NA

# 5. Main function
if __name__ == "__main__":
    tpc_path = "TPC3_df.tsv" # Protein matrix - protein trace source
    pairs_path = "TPC_1v5_FP_95w_MGL_20250603.tsv" # Pair matrix - check the column names
    output_path = "TPC_95wFP_wcc_scores_optimized_1v5_no3_20250608.tsv" # Output files    
    print("Loading and preprocessing data...") # For large size data - loading
    tpc_data, id_dict, pairs_df = load_and_preprocess(tpc_path, pairs_path)    
    print("Calculating WCC scores...") # Calculate wcc - calculating
    results = compute_wcc_parallel(
        tpc_data, 
        id_dict, 
        pairs_df,
        chunk_size=10000, # chunk size
        n_jobs=min(16, os.cpu_count()-1)) # CPU core - 16 is not the maximum
    result_df = pd.DataFrame(results)
    result_df.to_csv(output_path, sep='\t', index=False)    
    # # Optional
    # valid_scores = result_df['WCC_score'].dropna()
    # print(f"\nFinish WCC score without NA: {len(valid_scores)}/{len(result_df)}") # check the results
    # print("Top 5 WCC scores:\n", valid_scores.nlargest(5))