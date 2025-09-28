# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 20:13:03 2024

@author: PC
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import re  # For extracting numbers from column names

os.chdir("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_F_0710")

# Step 1: Read the real protein list from the TSV file
protein_list_file_path = 'E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\NthyvsFTC238_20250710.tsv'  # replace with actual file path
protein_list = pd.read_csv(protein_list_file_path, sep='\t')

# Extract individual protein IDs from the "combined" column
protein_list_split = protein_list['combine'].str.split('_', expand=True)
protein_list_split.columns = ['protein1', 'protein2']

# Initialize a dataframe to store the results

result_df = pd.DataFrame(columns=['group', 'protein_pair', 'auc_protein1', 'auc_protein2', 'max_col_protein1', 'max_col_protein2'])

# Define data file groups
tpc_files = [
    'E:\\MenggeLYU\\Figure_NC\\FTC238_202410\\FTC238NO1_gaussian_withoutlog2.tsv',
    'E:\\MenggeLYU\\Figure_NC\\FTC238_202410\\FTC238NO2_gaussian_withoutlog2.tsv',
    'E:\\MenggeLYU\\Figure_NC\\FTC238_202410\\FTC238NO3_gaussian_withoutlog2.tsv'
]

nthy_files = [
    'E:\\MenggeLYU\\Figure_NC\\Nthy_202410\\Nthy1_gaussian_withoutlog2.tsv',
    'E:\\MenggeLYU\\Figure_NC\\Nthy_202410\\Nthy2_gaussian_withoutlog2.tsv',
    'E:\\MenggeLYU\\Figure_NC\\Nthy_202410\\Nthy3_gaussian_withoutlog2.tsv'
]

groups = {'FTC238': tpc_files, 'Nthy': nthy_files}

# Function to extract numeric part from column names
def extract_numeric_part(colname):
    match = re.search(r'_(\d+)', colname)
    return int(match.group(1)) if match else None

# Step 2: Process each protein pair and generate curves for both groups
for idx, row in protein_list_split.iterrows():
    protein1 = row['protein1']
    protein2 = row['protein2']

    plt.figure(figsize=(10, 6))  # Create a new figure for each protein pair
    
    for group_name, data_files in groups.items():
        # Initialize lists to store protein data across files in the current group
        protein1_data = []
        protein2_data = []
        numeric_x_axis = None
        
        for data_file_path in data_files:
            # Step 3: Read each data file
            data = pd.read_csv(data_file_path, sep='\t', index_col=0)

            # Normalize each sample by dividing by the sum of all protein intensities in that sample
            data = data.div(data.sum(axis=0), axis=1)

            # Check if the proteins exist in the data
            if protein1 in data.index and protein2 in data.index:
                # Extract the numeric part of the column names for the x-axis
                numeric_x_axis = [extract_numeric_part(col) for col in data.columns]

                # Collect the data for averaging
                protein1_data.append(data.loc[protein1].to_numpy())
                protein2_data.append(data.loc[protein2].to_numpy())
            else:
                print(f'{protein1} or {protein2} not found in {data_file_path}.')

        # Step 4: Average the data for each protein in the current group
        if protein1_data and protein2_data:
            avg_protein1 = sum(protein1_data) / len(protein1_data)
            avg_protein2 = sum(protein2_data) / len(protein2_data)

            # Plot the average curves
            color = 'blue' if group_name == 'FTC238' else 'red'
            plt.plot(numeric_x_axis, avg_protein1, label=f'{protein1} ({group_name})', color=color, linestyle='--' if group_name == 'Nthy' else '-')
            plt.plot(numeric_x_axis, avg_protein2, label=f'{protein2} ({group_name})', color=color, linestyle=':' if group_name == 'Nthy' else '-.')

            # Step 5: Calculate the area under the curve (AUC) for both proteins
            auc_protein1 = simpson(avg_protein1, dx=1)  # numerical integration
            auc_protein2 = simpson(avg_protein2, dx=1)

            # Step 6: Find the column where the maximum occurs
            max_col_protein1 = numeric_x_axis[avg_protein1.argmax()]
            max_col_protein2 = numeric_x_axis[avg_protein2.argmax()]

            # Store the results in the dataframe
            result_df = pd.concat([
                result_df,
                pd.DataFrame({
                    'group': [group_name],
                    'protein_pair': [f'{protein1}_{protein2}'],
                    'auc_protein1': [auc_protein1],
                    'auc_protein2': [auc_protein2],
                    'max_col_protein1': [max_col_protein1],
                    'max_col_protein2': [max_col_protein2]
                })
            ], ignore_index=True)

    # Plot settings for the current protein pair
    plt.xlabel('The number of fractions')
    plt.ylabel('Normalized Intensity')
    plt.title(f'Curves for {protein1} and {protein2} in FTC238 and Nthy Groups')

    plt.legend()
    plt.tight_layout()

    # Automatically save the plot as a PDF
    pdf_filename = f'{protein1}_{protein2}_group_comparison.pdf'
    plt.savefig(pdf_filename, format='pdf')
    plt.close()  # Close the plot to avoid displaying it in the console

# Step 7: Display the resulting dataframe
print(result_df)

# Optionally, save the results to a CSV file
result_df.to_csv('NthyandFTC238_maxcol_forMWfiltered_20250710.csv', index=False)
