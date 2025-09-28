# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 09:14:23 2025

@author: PC
"""

import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import euclidean
from numba import njit
from joblib import Parallel, delayed
from tqdm.auto import tqdm
import warnings

warnings.filterwarnings('ignore')

os.chdir("E:\\Feature_generation")

# 1. Fast pearson correlation
@njit(cache=True)
def numba_pearson(x, y):    
    # collect n_valid
    x_vals = []
    y_vals = []
    for i in range(len(x)):
        if not np.isnan(x[i]) and not np.isnan(y[i]):
            x_vals.append(x[i])
            y_vals.append(y[i])    
    n = len(x_vals)
    if n < 2:  # at least n_valid
        return np.nan, np.nan    
    # mean
    sum_x = sum(x_vals)
    sum_y = sum(y_vals)
    mean_x = sum_x / n
    mean_y = sum_y / n   
    # cv and std
    cov = 0.0
    std_x = 0.0
    std_y = 0.0
    for i in range(n):
        x_diff = x_vals[i] - mean_x
        y_diff = y_vals[i] - mean_y
        cov += x_diff * y_diff
        std_x += x_diff ** 2
        std_y += y_diff ** 2    
    std_x = np.sqrt(std_x / n)
    std_y = np.sqrt(std_y / n)    
    if std_x == 0 or std_y == 0:
        return np.nan, np.nan    
    r = (cov / n) / (std_x * std_y)    
    # P value of pearson corr (t distribution)
    if np.abs(r) == 1.0:
        p = 0.0
    else:
        t = r * np.sqrt((n - 2) / (1 - r**2))
        # two-way t test
        p = 2 * np.exp(-0.717 * np.abs(t) - 0.416 * t**2)
        p = min(p, 1.0)      
    return r, p

# 2. Gaussian smoothing
def gaussian_smoothing(profile, window_size=5, std=None):   
    # default std = window_size/4，consistent with R code    
    if len(profile) < window_size:
        return profile    
    if std is None:
        std = window_size / 4 # window_size setting    
    x = np.linspace(-(window_size-1)/2, (window_size-1)/2, window_size)
    window = np.exp(-x**2/(2*std**2))
    window /= window.sum()    
    smoothed = np.convolve(profile, window, mode='same') # boundary setting
    return smoothed

# 3. De noise
def clean_profile(profile, impute_NA=True, smooth=True, smooth_width=4, noise_floor=0.001):
    # each protein - clean profile
    profile = profile.astype(np.float64).copy()    
    # setting mise floor
    profile[profile < noise_floor] = noise_floor    
    # missing imputation - col_mean  
    if impute_NA:
        col_mean = np.nanmean(profile)
        profile[np.isnan(profile)] = col_mean    
    # smooth
    if smooth and smooth_width > 1:
        profile = gaussian_smoothing(profile, window_size=smooth_width)    
    return profile

def clean_profiles(profile_matrix, impute_NA=True, smooth=True, smooth_width=4, noise_floor=0.001):
    # each matrix - clean profile
    cleaned = np.apply_along_axis(
        lambda x: clean_profile(x, impute_NA, smooth, smooth_width, noise_floor),
        axis=1,
        arr=profile_matrix.values if isinstance(profile_matrix, pd.DataFrame) else profile_matrix)    
    if isinstance(profile_matrix, pd.DataFrame):
        cleaned = pd.DataFrame(cleaned, index=profile_matrix.index, columns=profile_matrix.columns)    
    return cleaned

# 4. Feature generation
def calculate_features_for_pair(pair, tpc_data, gaussians_dict=None, smooth_width=4):
    # features per pairs
    prot_a, prot_b = pair['protein_A'], pair['protein_B'] # pair matrix column name  
    try:
        # input data - here is the example tpc-1 matrix
        vec_a_raw = tpc_data.loc[prot_a].values
        vec_b_raw = tpc_data.loc[prot_b].values        
        # clean_profile - each matrix
        smoothed_matrix = clean_profiles(
            tpc_data.loc[[prot_a, prot_b]], 
            impute_NA=True, 
            smooth=True, 
            smooth_width=smooth_width)
        vec_a_smoothed = smoothed_matrix.loc[prot_a].values
        vec_b_smoothed = smoothed_matrix.loc[prot_b].values        
    except KeyError:
        return None  # if the ID is not existed in input data    
    # Raw pearson and P value
    pearson_raw, pearson_raw_p = numba_pearson(vec_a_raw, vec_b_raw)    
    # Smoothed pearson (no P value)
    pearson_smoothed, _ = numba_pearson(vec_a_smoothed, vec_b_smoothed)   
    # n_valid RAW
    valid_idx = ~(np.isnan(vec_a_raw) | np.isnan(vec_b_raw))
    n_pairs = valid_idx.sum()    
    if n_pairs < 3 or np.isnan(pearson_raw):
        return None   
    # n_valid SMOOTH
    valid_data_a = vec_a_smoothed[valid_idx]
    valid_data_b = vec_b_smoothed[valid_idx]    
    # Feature generation
    features = {
        'perturbation_A': prot_a,
        'perturbation_B': prot_b,
        'n_pairs': int(n_pairs),
        'pearson_R_raw': float(1-pearson_raw),
        # 'pearson_pvalue_raw': float(pearson_raw_p), 
        'pearson_R_smoothed': float(1-pearson_smoothed),
        'euclidean_distance': euclidean(valid_data_a, valid_data_b),
        'co_peak': abs(np.nanargmax(valid_data_a) - np.nanargmax(valid_data_b))}    
    # Co peak (if have gaussian model)
    # gaussians_dict = {
    # 'protein_A': {'means': [1.0, 2.0], 'stds': [0.1, 0.2]},
    # 'protein_B': {'means': [1.5, 2.5], 'stds': [0.15, 0.25]}}
    if gaussians_dict and prot_a in gaussians_dict and prot_b in gaussians_dict:
        gmm_a = gaussians_dict[prot_a]
        gmm_b = gaussians_dict[prot_b]
        if gmm_a and gmm_b:
            features['co_apex'] = min(
                abs(ma - mb) 
                for ma in gmm_a['means'] 
                for mb in gmm_b['means'])    
    return features

# 5. Chunk
def process_chunk(chunk, tpc_data, gaussians_dict=None, smooth_width=4):
    results = []
    for _, row in chunk.iterrows():
        result = calculate_features_for_pair(row, tpc_data, gaussians_dict, smooth_width)
        if result:
            results.append(result)
    return results
# parallel
def compute_features_parallel(pairs_df, tpc_data, gaussians_dict=None, 
                            chunk_size=10000, n_jobs=-1, smooth_width=4):
    # 10,000 per chunk
    n_chunks = len(pairs_df) // chunk_size + 1
    chunks = np.array_split(pairs_df, n_chunks)    
    # Parallel
    all_results = Parallel(n_jobs=n_jobs)(
        delayed(process_chunk)(chunk, tpc_data, gaussians_dict, smooth_width)
        for chunk in tqdm(chunks, desc="Processing chunks")
    )   
    # Combine results
    features_list = [item for sublist in all_results for item in sublist]
    return pd.DataFrame(features_list)
# load protein matrix - TPC-1 as example
def load_tpc_data(tpc_path):
    tpc = pd.read_csv(tpc_path, sep='\t')
    # re sort the colnames from 1 to 60 (should with out character)
    numeric_cols = [col for col in tpc.columns if col not in ['ID'] and str(col).isdigit()]
    numeric_cols_sorted = sorted(numeric_cols, key=lambda x: int(x))
    # ID + NO.1-60
    tpc = tpc[['ID'] + numeric_cols_sorted]
    return tpc.set_index('ID')
# load protein pair matrix - two columns (protein_A and protein_B, check the column name)
def load_pairs_data(pairs_path):
    return pd.read_csv(pairs_path, sep='\t')

# 6. Data input
if __name__ == "__main__":
    tpc_data = load_tpc_data("TPC_1v5_perturbation_MGL_20250529.tsv")
    pairs_df = load_pairs_data("TPC_1v5_perturbationtwocolumns_20250603.tsv")    
    # smooth_width- consistent with R
    smooth_width = 4    
    # feature generation
    features_df = compute_features_parallel(
        pairs_df, 
        tpc_data, 
        gaussians_dict=None,  # default is None (if there is no gaussian model), gaussians_dict refers to use model
        chunk_size=10000,
        n_jobs=18,
        smooth_width=smooth_width)    
    features_df.to_csv("Computed_features_with_smoothed.tsv", sep='\t', index=False)