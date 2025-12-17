#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import numpy as np
from scipy.integrate import simps
import re

os.chdir("/path/visualization")

# input files
data_files = [
    'Nthy1_gaussian_withoutlog2.tsv',
    'Nthy2_gaussian_withoutlog2.tsv',
    'Nthy3_gaussian_withoutlog2.tsv'
]

# output dir
output_dir = "normalized_matrices"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# creat storage
normalized_matrices = {}

for data_file in data_files:
    print(f"Processing: {data_file}")    
    # read input
    data = pd.read_csv(data_file, sep='\t', index_col=0)
    print(f"Original data shape: {data.shape}")
    print(f"First few columns: {data.columns[:5].tolist()}")
    # validate the normalization
    normalized_data = data.div(data.sum(axis=0), axis=1)
    col_sums = normalized_data.sum(axis=0)
    print(f"Column sums after normalization: {col_sums.head().values.round(1)}")
    print(f"All columns sum to 1: {np.allclose(col_sums, 1.0, atol=1e-10)}")    
    # save to files
    base_name = os.path.splitext(os.path.basename(data_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}_normalized.tsv")
    normalized_data.to_csv(output_file, sep='\t', float_format='%.10e')
    print(f"Saved to: {output_file}")
    print(f"File size: {os.path.getsize(output_file)} bytes\n")
    # save to storage
    normalized_matrices[base_name] = normalized_data
print("All files processed successfully!")
print("Available normalized matrices:", list(normalized_matrices.keys()))

# details of each matrix
for matrix_name, matrix_data in normalized_matrices.items():
    print(f"\n{matrix_name} matrix details:")
    print(f"Shape: {matrix_data.shape}")
    print(f"Total proteins: {len(matrix_data)}")
    print(f"Total fractions: {len(matrix_data.columns)}")
    print("First 5 proteins and first 5 fractions:")
    print(matrix_data.iloc[:5, :5])
    
    stats_file = os.path.join(output_dir, f"{matrix_name}_stats.txt")
    with open(stats_file, 'w') as f:
        f.write(f"Matrix: {matrix_name}\n")
        f.write(f"Shape: {matrix_data.shape}\n")
        f.write(f"Total proteins: {len(matrix_data)}\n")
        f.write(f"Total fractions: {len(matrix_data.columns)}\n")
        f.write("\nColumn sums (should all be 1.0):\n")
        f.write(str(matrix_data.sum(axis=0).round(6).to_dict()))
        f.write("\n\nFirst 10 proteins:\n")
        f.write(str(matrix_data.index[:10].tolist()))
        f.write("\n\nAll fractions:\n")
        f.write(str(matrix_data.columns.tolist()))

print("\nDetailed statistics files have been saved in the normalized_matrices folder.")



