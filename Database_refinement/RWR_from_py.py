# -*- coding: utf-8 -*-
import polars as pl
import networkx as nx
import numpy as np
import os

# Set file path
input_folder = "/storage/generateBinaryNetwork/Part2/Ten"
output_folder = "/storage/generateBinaryNetwork/Part2/Ten"

# Make sure the output path
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Random Walk with Restart 
def random_walk_with_restart(G, start_node, p, r, num_iterations):
    probabilities = {node: 0.0 for node in G.nodes()}
    probabilities[start_node] = 1.0

    for _ in range(num_iterations):
        new_probabilities = probabilities.copy()
        for node, prob in probabilities.items():
            if np.random.rand() < r or len(list(G.successors(node))) == 0:
                new_probabilities[node] += prob * r
            else:
                successors = list(G.successors(node))
                if successors:
                    prob_to_successors = prob * p / len(successors)
                    for successor in successors:
                        new_probabilities[successor] += prob_to_successors
                new_probabilities[node] *= (1 - p)

        probabilities = new_probabilities

    # Normalized probability
    total_prob = sum(probabilities.values())
    normalized_probabilities = {k: v / total_prob for k, v in probabilities.items()}
    return normalized_probabilities

# Iterate through all TSV files in the folder
for filename in os.listdir(input_folder):
    if filename.endswith(".tsv"):
        input_path = os.path.join(input_folder, filename)
        
        # Read the first row of tsv as start_node
        first_row_df = pl.read_csv(input_path, separator='\t').head(1)  # read first row
        if not first_row_df.is_empty():
            start_node = first_row_df[0, "a"]  # get the first column

            # read full tsv
            df = pl.read_csv(input_path, separator='\t')
            edges = [tuple(row[:2]) for row in df.iter_rows()]  # the eadg is starting the start of two rows
            G = nx.DiGraph()
            G.add_edges_from(edges)
            
            # set RWR 
            p = 0.5  # prob
            r = 0.5  # retart
            num_iterations = 1000  # iteration

            # RWR function
            distribution = random_walk_with_restart(G, start_node, p, r, num_iterations)
            
            # output fold
            output_path = os.path.join(output_folder, "{}_rwr.tsv".format(os.path.splitext(filename)[0]))
            
            # write to the output file
            with open(output_path, "w", encoding="utf-8") as f:
                f.write("Node\tValue\n")
                for node, value in distribution.items():
                    f.write("{}\t{:.4f}\n".format(node, value))
            
            print("Random Walk with Restart distribution has been written to '{}'".format(output_path))
        else:
            print("File '{}' is empty, skipping...".format(input_path))

