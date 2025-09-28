import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.integrate import simps
import re  # For extracting numbers from column names


os.chdir("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_F_auc_0710\\")
# Step 1: Read the real protein list from the TSV file
protein_list_file_path = 'E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_F_pair_list_0710.tsv'  # replace with actual file path
protein_list = pd.read_csv(protein_list_file_path, sep='\t')

# Extract individual protein IDs from the "combined" column
protein_list_split = protein_list['combine'].str.split('_', expand=True)
protein_list_split.columns = ['protein1', 'protein2']


result_df = pd.DataFrame(columns=['protein_pair', 'data_file', 'protein_label', 'auc', 'max_col'])

# List of data files to pick protein pairs from
data_files = [
    'E:\\MenggeLYU\\Figure_NC\\Nthy_202410\\Nthy1_gaussian_withoutlog2.tsv',
    'E:\\MenggeLYU\\Figure_NC\\Nthy_202410\\Nthy2_gaussian_withoutlog2.tsv',
    'E:\\MenggeLYU\\Figure_NC\\Nthy_202410\\Nthy3_gaussian_withoutlog2.tsv',
    'E:\\MenggeLYU\\Figure_NC\\FTC238_202410\\FTC238NO1_gaussian_withoutlog2.tsv',
    'E:\\MenggeLYU\\Figure_NC\\FTC238_202410\\FTC238NO2_gaussian_withoutlog2.tsv',
    'E:\\MenggeLYU\\Figure_NC\\FTC238_202410\\FTC238NO3_gaussian_withoutlog2.tsv'
]

# Function to extract numeric part from column names
def extract_numeric_part(colname):
    return int(re.search(r'_(\d+)', colname).group(1))

# Step 2: Pick the rows corresponding to the proteins and plot curves from multiple files
for idx, row in protein_list_split.iterrows():
    protein1 = row['protein1']
    protein2 = row['protein2']
    protein_pair = f'{protein1}_{protein2}'

    plt.figure(figsize=(10, 6))  # Create a new figure for each protein pair
    for data_file_path in data_files:
        # Step 3: Read each data file
        data = pd.read_csv(data_file_path, sep='\t', index_col=0)
        
        # Normalize each sample column by dividing by the total intensity of the sample
        data = data.div(data.sum(axis=0), axis=1)

        # Check if the proteins exist in the data
        if protein1 in data.index and protein2 in data.index:
            # Step 4: Plot the curves for protein1 and protein2 from each file
            base_file_name = os.path.splitext(os.path.basename(data_file_path))[0].rsplit('_', 1)[0]  # Get the base name without '_gaussian_withoutlog2.tsv'
            label1 = f'{protein1} ({base_file_name})'
            label2 = f'{protein2} ({base_file_name})'

            # Extract the numeric part of the column names for the x-axis
            numeric_x_axis = [extract_numeric_part(col) for col in data.columns]

            # Convert to numpy array to avoid FutureWarning
            plt.plot(numeric_x_axis, data.loc[protein1].to_numpy(), label=label1)
            plt.plot(numeric_x_axis, data.loc[protein2].to_numpy(), label=label2)
            
            # Step 5: Calculate the area under the curve (AUC) for both proteins
            auc_protein1 = simps(data.loc[protein1], dx=1)  # numerical integration
            auc_protein2 = simps(data.loc[protein2], dx=1)

            # Step 6: Find the column where the maximum occurs
            max_col_protein1 = data.loc[protein1].idxmax()
            max_col_protein2 = data.loc[protein2].idxmax()

            # Store the results for both proteins in the dataframe
            new_row = pd.DataFrame([{
                    'protein_pair': protein_pair,
                    'data_file': base_file_name,
                    'protein_label': protein1,
                    'auc': auc_protein1,
                    'max_col': max_col_protein1}])
            result_df = pd.concat([result_df, new_row], ignore_index=True)

            
            new_row = pd.DataFrame([{
                    'protein_pair': protein_pair,
                    'data_file': base_file_name,
                    'protein_label': protein2,
                    'auc': auc_protein2,
                    'max_col': max_col_protein2}])
            result_df = pd.concat([result_df, new_row], ignore_index=True)
        else:
            print(f'{protein1} or {protein2} not found in {data_file_path}.')

    # Remove all but one X-axis label
    plt.xlabel('The number of fractions')
    plt.ylabel('Normalized Intensity')
    plt.title(f'Curves for {protein1} and {protein2} from Multiple Data Files')
    
    plt.legend()
    plt.tight_layout()

    # Automatically save the plot as a PDF
    pdf_filename = f'{protein1}_{protein2}_combined.pdf'
    plt.savefig(pdf_filename, format='pdf')
    plt.close()  # Close the plot to avoid displaying it in the console

# Step 7: Display the resulting dataframe
print(result_df)

# Optionally, save the results to a CSV file
result_df.to_csv('N_F_AUC_forDEPPS20250712.csv', index=False)
