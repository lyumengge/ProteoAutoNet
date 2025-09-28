# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 23:07:00 2025

@author: PC
"""

import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import euclidean
from tqdm import tqdm
import warnings
from numba import njit
from joblib import Parallel, delayed
import multiprocessing
from datetime import datetime

warnings.filterwarnings('ignore') 

os.chdir("D:/TP_perturbation/")

# 1. Config setting
class Config:
    N_JOBS = max(multiprocessing.cpu_count() - 2, 1)  # core number total-2
    BATCH_SIZE = 500  # 500 per run                               
    MIN_VALID_POINTS = 4  # min number of one trace                          
    TOP_N = 50  # choose 50 perturbations per pair                                    
    REP_RATIO = 0.8  # choose 50*0.8 = 40 pairs with high similarity per trace

# 2. Fast pairwise Pearson correlation
@njit(fastmath=True, cache=True)
def fast_pearson(x, y):
    """Fast Pearson corr"""
    n = len(x)
    sum_x = sum_y = sum_xy = sum_x2 = sum_y2 = 0.0
    valid = 0    
    for i in range(n):
        xi, yi = x[i], y[i]
        if not (np.isnan(xi) or np.isnan(yi)):
            sum_x += xi
            sum_y += yi
            sum_xy += xi * yi
            sum_x2 += xi * xi
            sum_y2 += yi * yi
            valid += 1    
    if valid < 2:
        return np.nan    
    tmp = (valid * sum_xy - sum_x * sum_y)
    denom = np.sqrt((valid * sum_x2 - sum_x**2) * (valid * sum_y2 - sum_y**2))
    return tmp / denom if denom != 0 else np.nan

# 3. Framework per batch
def process_batch(pairs_batch, perturbations):
    # 500 per batch
    batch_results = []
    protein_set = set(perturbations['protein'].unique())    
    for _, pair in pairs_batch.iterrows():
        protein_A, protein_B = pair['protein_A'], pair['protein_B']        
        # skip the protein is not existed
        if protein_A not in protein_set or protein_B not in protein_set:
            continue            
        # get perturbation
        mask_A = perturbations['protein'] == protein_A
        mask_B = perturbations['protein'] == protein_B
        A_data = perturbations.loc[mask_A].iloc[:, :-2].values.astype(np.float32)  # float 32 
        B_data = perturbations.loc[mask_B].iloc[:, :-2].values.astype(np.float32)
        A_ids = perturbations.loc[mask_A, 'perturbation_id'].values
        B_ids = perturbations.loc[mask_B, 'perturbation_id'].values        
        # generation pairs score
        scores = []
        for i in range(len(A_data)):
            dataA = A_data[i]
            for j in range(len(B_data)):
                dataB = B_data[j]                
                # get the mask of per pairs
                valid_A = ~np.isnan(dataA)
                valid_B = ~np.isnan(dataB)
                valid_mask = valid_A & valid_B
                n_valid = valid_mask.sum()                
                if n_valid < Config.MIN_VALID_POINTS:
                    continue                
                # mask_A and mask _B
                valid_A = dataA[valid_mask]
                valid_B = dataB[valid_mask]
                pearson_r = fast_pearson(valid_A, valid_B)
                # euclidean distance              
                if not np.isnan(pearson_r):
                    euc_dist = euclidean(valid_A, valid_B) / n_valid
                    # generated score
                    score = n_valid * 0.5 + abs(pearson_r) * 30 - euc_dist * 0.01
                    scores.append((
                        protein_A, protein_B, 
                        A_ids[i], B_ids[j],
                        score, pearson_r, euc_dist, n_valid))        
        if len(scores) >= 10:
            # fast align and select
            scores_arr = np.array(
                scores, 
                dtype=[('protein_A', 'U20'), ('protein_B', 'U20'),
                       ('perturbation_A', 'U50'), ('perturbation_B', 'U50'),
                       ('score', float), ('pearson_r', float),
                       ('euclidean_dist', float), ('n_valid', int)])
            scores_arr.sort(order='score')
            scores_arr = scores_arr[::-1]  #descending the score             
            # choose top 40 similarity - can change the percentage in config
            n_total = min(Config.TOP_N, len(scores_arr))
            n_rep = int(n_total * Config.REP_RATIO)            
            # store results
            rep = scores_arr[:n_rep]
            non_rep = scores_arr[n_rep:n_total]
            # make the similarity label in case need
            for item in rep:
                batch_results.append((*item, 'representative'))
            for item in non_rep:
                batch_results.append((*item, 'non_representative'))
    
    return batch_results

# 4. Running time and data loading
def preprocess_perturbations(perturbations):
    # match protein_pairs and perturbations
    perturbations['protein'] = perturbations.index.str.split('_').str[0]
    perturbations['perturbation_id'] = perturbations.index
    return perturbations

def main():
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")    
    # data load - all original positive pairs(scpre>0.5 from RF)
    protein_pairs = pd.read_csv("Nthy_1v5_pairs_20250603.tsv", sep='\t', usecols=['protein_A', 'protein_B'])
    perturbations = pd.read_csv("Nthy_1v5_perturbation_20250529.tsv", sep='\t', index_col=0,
        dtype={col: np.float32 for col in pd.read_csv("Nthy_1v5_perturbation_20250529.tsv", sep='\t', nrows=1).columns[1:]}) # float 32 input
    perturbations = preprocess_perturbations(perturbations)    
    # loading chunks
    batch_size = 1000  # load pairs from all original positive pairs
    batches = [protein_pairs.iloc[i:i + batch_size] for i in range(0, len(protein_pairs), batch_size)]    
    # parallel - use tqdm as processing
    all_results = []
    with tqdm(total=len(batches), desc="Processes") as pbar:
        for batch in batches:
            results = Parallel(n_jobs=Config.N_JOBS, batch_size=Config.BATCH_SIZE)(
                delayed(process_batch)(sub_batch, perturbations)
                for _, sub_batch in batch.groupby(np.arange(len(batch)) // (batch_size//Config.N_JOBS)))
            all_results.extend([item for sublist in results for item in sublist])
            pbar.update(1)    
    # save results
    if all_results:
        cols = ['protein_A', 'protein_B', 'perturbation_A', 'perturbation_B',
                'score', 'pearson_r', 'euclidean_dist', 'n_valid', 'type']
        final_df = pd.DataFrame(all_results, columns=cols)
        final_df.to_csv("Nthy_1v5_final_optimized_50per_perturbations_20250602.tsv", sep='\t', index=False)
        
        print(f"\nFinish time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Process number of pairs: {len(protein_pairs)}")
        print(f"Output the number of perturbations: {len(final_df)}")
        print(f"Representative ratio: {len(final_df[final_df['type']=='representative'])/len(final_df):.1%}")


if __name__ == "__main__":
    main()