# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 10:10:15 2025

@author: PC
"""

import pandas as pd
import numpy as np

# Support 3 replicates
ori1 = pd.read_csv("rep1.tsv", sep="\t")
ori2 = pd.read_csv("rep2.tsv", sep="\t")
ori3 = pd.read_csv("rep3.tsv", sep="\t")
all_data = pd.concat([ori1, ori2, ori3], keys=['ori1', 'ori2', 'ori3'])
common_keys = all_data.groupby(['perturbation_A', 'perturbation_B']).filter(lambda x: len(x) == 3) # exsited in 3 reps
two_common_keys = all_data.groupby(['perturbation_A', 'perturbation_B']).filter(lambda x: len(x) == 2) # exsited in 2 reps
unique_keys = all_data.groupby(['perturbation_A', 'perturbation_B']).filter(lambda x: len(x) == 1) # exsited in 1 reps

# process 3 reps
if not common_keys.empty:
    common_keys = common_keys.groupby(['perturbation_A', 'perturbation_B']).mean().reset_index()
    common_keys['n_pairs'] = common_keys['n_pairs'].astype(int)
    common_keys['co_peak'] = np.ceil(common_keys['co_peak']).astype(int) # transfer int type
# process 2 reps
if not two_common_keys.empty:
    two_common_keys = two_common_keys.groupby(['perturbation_A', 'perturbation_B']).mean().reset_index()
    two_common_keys['n_pairs'] = two_common_keys['n_pairs'].astype(int)
    two_common_keys['co_peak'] = np.ceil(two_common_keys['co_peak']).astype(int) # transfer int type
# process 1 reps
if not unique_keys.empty:
    unique_keys = unique_keys.reset_index(level=[0], drop=True)
    unique_keys['n_pairs'] = unique_keys['n_pairs'].astype(int)
    unique_keys['co_peak'] = np.ceil(unique_keys['co_peak']).astype(int) # transfer int type

# Comb results
result = pd.concat([common_keys, two_common_keys, unique_keys]).drop_duplicates(subset=['perturbation_A', 'perturbation_B'], keep='first')
result = result.reset_index(drop=True) # reset index
result.to_csv("Three_rep_comb.tsv", sep = "\t", index=False)

'''
If you only have 2 replicates
Please use the following codes
'''

# Support 2 reps
ori1 = pd.read_csv("rep1.tsv", sep="\t")
ori2 = pd.read_csv("rep2.tsv", sep="\t")
all_data = pd.concat([ori1, ori2], keys=['ori1', 'ori2'])
two_common_keys = all_data.groupby(['perturbation_A', 'perturbation_B']).filter(lambda x: len(x) == 2) # exsited in 2 reps
unique_keys = all_data.groupby(['perturbation_A', 'perturbation_B']).filter(lambda x: len(x) == 1) # exsited in 1 reps

# process 2 reps
if not two_common_keys.empty:
    two_common_keys = two_common_keys.groupby(['perturbation_A', 'perturbation_B']).mean().reset_index()
    two_common_keys['n_pairs'] = two_common_keys['n_pairs'].astype(int)
    two_common_keys['co_peak'] = np.ceil(two_common_keys['co_peak']).astype(int)
# process 1 reps
if not unique_keys.empty:
    unique_keys = unique_keys.reset_index(level=[0], drop=True)
    unique_keys['n_pairs'] = unique_keys['n_pairs'].astype(int)
    unique_keys['co_peak'] = np.ceil(unique_keys['co_peak']).astype(int)

# Comb results
result = pd.concat([two_common_keys, unique_keys]).drop_duplicates(subset=['perturbation_A', 'perturbation_B'], keep='first')
result = result.reset_index(drop=True) # reset index
result.to_csv("Two_rep_comb.tsv", sep = "\t", index=False)