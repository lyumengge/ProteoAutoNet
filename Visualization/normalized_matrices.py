#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os
import numpy as np
from scipy.integrate import simps
import re


# In[3]:


# 设置工作目录
os.chdir("G:/Tower3/NC_rebuttal_t3/New_feature_RF/ComplexXGB/DEPPS")

# 数据文件列表
data_files = [
    'G:/Tower3/Figure_NC/Nthy_202410/Nthy1_gaussian_withoutlog2.tsv',
    'G:/Tower3/Figure_NC/Nthy_202410/Nthy2_gaussian_withoutlog2.tsv',
    'G:/Tower3/Figure_NC/Nthy_202410/Nthy3_gaussian_withoutlog2.tsv'
]

# 创建输出目录
output_dir = "normalized_matrices"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 存储所有normalized矩阵的字典
normalized_matrices = {}

for data_file in data_files:
    print(f"Processing: {data_file}")
    
    # 读取数据
    data = pd.read_csv(data_file, sep='\t', index_col=0)
    
    print(f"Original data shape: {data.shape}")
    print(f"First few columns: {data.columns[:5].tolist()}")
    
    # NORMALIZATION: 每列除以该列的总和
    normalized_data = data.div(data.sum(axis=0), axis=1)
    
    # 验证normalization
    col_sums = normalized_data.sum(axis=0)
    print(f"Column sums after normalization: {col_sums.head().values.round(1)}")
    print(f"All columns sum to 1: {np.allclose(col_sums, 1.0, atol=1e-10)}")
    
    # 保存到文件
    base_name = os.path.splitext(os.path.basename(data_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}_normalized.tsv")
    
    # 写入完整的矩阵内容
    normalized_data.to_csv(output_file, sep='\t', float_format='%.10e')
    print(f"Saved to: {output_file}")
    print(f"File size: {os.path.getsize(output_file)} bytes\n")
    
    # 存储到字典
    normalized_matrices[base_name] = normalized_data

print("All files processed successfully!")
print("Available normalized matrices:", list(normalized_matrices.keys()))

# 显示每个矩阵的详细信息
for matrix_name, matrix_data in normalized_matrices.items():
    print(f"\n{matrix_name} matrix details:")
    print(f"Shape: {matrix_data.shape}")
    print(f"Total proteins: {len(matrix_data)}")
    print(f"Total fractions: {len(matrix_data.columns)}")
    print("First 5 proteins and first 5 fractions:")
    print(matrix_data.iloc[:5, :5])
    
    # 也可以选择保存更详细的统计信息
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


# In[5]:


# 设置工作目录
os.chdir("G:/Tower3/NC_rebuttal_t3/New_feature_RF/ComplexXGB/DEPPS")

# 数据文件列表
data_files = [
    'G:/Tower3/Figure_NC/FTC238_202410/FTC238NO1_gaussian_withoutlog2.tsv',
    'G:/Tower3/Figure_NC/FTC238_202410/FTC238NO2_gaussian_withoutlog2.tsv',
    'G:/Tower3/Figure_NC/FTC238_202410/FTC238NO3_gaussian_withoutlog2.tsv'
]

# 创建输出目录
output_dir = "normalized_matrices"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 存储所有normalized矩阵的字典
normalized_matrices = {}

for data_file in data_files:
    print(f"Processing: {data_file}")
    
    # 读取数据
    data = pd.read_csv(data_file, sep='\t', index_col=0)
    
    print(f"Original data shape: {data.shape}")
    print(f"First few columns: {data.columns[:5].tolist()}")
    
    # NORMALIZATION: 每列除以该列的总和
    normalized_data = data.div(data.sum(axis=0), axis=1)
    
    # 验证normalization
    col_sums = normalized_data.sum(axis=0)
    print(f"Column sums after normalization: {col_sums.head().values.round(1)}")
    print(f"All columns sum to 1: {np.allclose(col_sums, 1.0, atol=1e-10)}")
    
    # 保存到文件
    base_name = os.path.splitext(os.path.basename(data_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}_normalized.tsv")
    
    # 写入完整的矩阵内容
    normalized_data.to_csv(output_file, sep='\t', float_format='%.10e')
    print(f"Saved to: {output_file}")
    print(f"File size: {os.path.getsize(output_file)} bytes\n")
    
    # 存储到字典
    normalized_matrices[base_name] = normalized_data

print("All files processed successfully!")
print("Available normalized matrices:", list(normalized_matrices.keys()))

# 显示每个矩阵的详细信息
for matrix_name, matrix_data in normalized_matrices.items():
    print(f"\n{matrix_name} matrix details:")
    print(f"Shape: {matrix_data.shape}")
    print(f"Total proteins: {len(matrix_data)}")
    print(f"Total fractions: {len(matrix_data.columns)}")
    print("First 5 proteins and first 5 fractions:")
    print(matrix_data.iloc[:5, :5])
    
    # 也可以选择保存更详细的统计信息
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


# In[6]:


# 设置工作目录
os.chdir("G:/Tower3/NC_rebuttal_t3/New_feature_RF/ComplexXGB/DEPPS")

# 数据文件列表
data_files = [
    'G:/Tower3/Figure_NC/TPC_202410/TPC1_gaussian_withoutlog2.tsv',
    'G:/Tower3/Figure_NC/TPC_202410/TPC2_gaussian_withoutlog2.tsv',
    'G:/Tower3/Figure_NC/TPC_202410/TPC3_gaussian_withoutlog2.tsv'
]

# 创建输出目录
output_dir = "normalized_matrices"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 存储所有normalized矩阵的字典
normalized_matrices = {}

for data_file in data_files:
    print(f"Processing: {data_file}")
    
    # 读取数据
    data = pd.read_csv(data_file, sep='\t', index_col=0)
    
    print(f"Original data shape: {data.shape}")
    print(f"First few columns: {data.columns[:5].tolist()}")
    
    # NORMALIZATION: 每列除以该列的总和
    normalized_data = data.div(data.sum(axis=0), axis=1)
    
    # 验证normalization
    col_sums = normalized_data.sum(axis=0)
    print(f"Column sums after normalization: {col_sums.head().values.round(1)}")
    print(f"All columns sum to 1: {np.allclose(col_sums, 1.0, atol=1e-10)}")
    
    # 保存到文件
    base_name = os.path.splitext(os.path.basename(data_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}_normalized.tsv")
    
    # 写入完整的矩阵内容
    normalized_data.to_csv(output_file, sep='\t', float_format='%.10e')
    print(f"Saved to: {output_file}")
    print(f"File size: {os.path.getsize(output_file)} bytes\n")
    
    # 存储到字典
    normalized_matrices[base_name] = normalized_data

print("All files processed successfully!")
print("Available normalized matrices:", list(normalized_matrices.keys()))

# 显示每个矩阵的详细信息
for matrix_name, matrix_data in normalized_matrices.items():
    print(f"\n{matrix_name} matrix details:")
    print(f"Shape: {matrix_data.shape}")
    print(f"Total proteins: {len(matrix_data)}")
    print(f"Total fractions: {len(matrix_data.columns)}")
    print("First 5 proteins and first 5 fractions:")
    print(matrix_data.iloc[:5, :5])
    
    # 也可以选择保存更详细的统计信息
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


# In[ ]:




