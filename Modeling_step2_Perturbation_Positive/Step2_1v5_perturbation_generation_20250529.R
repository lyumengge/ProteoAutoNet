##################################################
## Project: NC_rebuttal
## Script purpose: True sample preturbation
## Date: 2025-05-27
## Author: MenggeLYU
## Version: 1.6.1
##################################################

# path and library --------------------------------------------------------
setwd("E:\\final_ratio\\TP_perturbation")
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(purrr)
library(zoo)
library(gridExtra)
library(datawizard)
library(tidyverse)


# protein_pairs_over05 ----------------------------------------------------
tpc15 <- rio::import("TPC_1v5_score05_proteinpairs_features_MGL.tsv")
tpc_prot1 <- read.delim("TPC1_protmatrix_20240918.tsv", row.names=1, quote="")
tpc_prot2 <- read.delim("TPC2_protmatrix_20240918.tsv", row.names=1, quote="")
tpc_prot3 <- read.delim("TPC3_protmatrix_20240918.tsv", row.names=1, quote="")
colnames(tpc_prot1) <- gsub("TPC.1_NO1_30minDIA_","", colnames(tpc_prot1))
colnames(tpc_prot2) <- gsub("TPC.1_NO2_30minDIA_","", colnames(tpc_prot2))
colnames(tpc_prot3) <- gsub("TPC.1_NO3_30minDIA_","", colnames(tpc_prot3))
tpc_prot1 <- tpc_prot1[, order(as.numeric(colnames(tpc_prot1)))]
tpc_prot2 <- tpc_prot2[, order(as.numeric(colnames(tpc_prot2)))]
tpc_prot3 <- tpc_prot3[, order(as.numeric(colnames(tpc_prot3)))]
tpc_prot1 <- tpc_prot1[rowSums(is.na(tpc_prot1)) != ncol(tpc_prot1), ]
tpc_prot2 <- tpc_prot2[rowSums(is.na(tpc_prot2)) != ncol(tpc_prot2), ]
tpc_prot3 <- tpc_prot3[rowSums(is.na(tpc_prot3)) != ncol(tpc_prot3), ]


# sample_rows <- sample(nrow(tpc15), 10) # test row: first time need to use
tpc15_sample <- tpc15

prot1 <- tpc_prot1
prot2 <- tpc_prot2
prot3 <- tpc_prot3

# gaussian smoth function
gaussian_smooth <- function(x, window_size = 3) {
  if (length(x) < window_size) return(x)
  weights <- dnorm(seq(-2, 2, length.out = window_size))
  weights <- weights / sum(weights)
  rollapply(x, width = window_size, FUN = function(x) sum(x * weights), 
            fill = NA, align = "center")}

# perturbation
generate_perturbed_sample <- function(protA_orig, protB_orig, matrix_name) {
  # gaussian function
  protA_smooth <- gaussian_smooth(protA_orig)
  protB_smooth <- gaussian_smooth(protB_orig)

  # value perturbation 
  perturb_value <- function(v) {
    if (is.na(v)) return(NA)
    runif(1, 0.9 * v, 1.1 * v)}
  
  protA_perturbed <- sapply(protA_smooth, perturb_value)
  protB_perturbed <- sapply(protB_smooth, perturb_value)
  
  # NA three choice
  handle_na_perturbation <- function(orig, perturbed) {
    na_positions <- which(is.na(orig))
    for (pos in na_positions) {
      action <- sample(1:3, 1)
      if (action == 1) {
        perturbed[pos] <- NA
      } else if (action == 2 && pos > 1) {
        perturbed[pos - 1] <- NA
      } else if (action == 3 && pos < length(perturbed)) {
        perturbed[pos + 1] <- NA
      }
    }
    return(perturbed)}
  
  list(protA = handle_na_perturbation(protA_orig, protA_perturbed), protB = handle_na_perturbation(protB_orig, protB_perturbed),matrix = matrix_name)}

# each pair generate 100 perturbations from one rep
perturbed_samples <- list()

for (i in 1:nrow(tpc15_sample)) {
  protA_name <- tpc15_sample[i, "protein_A"]
  protB_name <- tpc15_sample[i, "protein_B"]
  
  # original 
  matrices <- list(prot1 = prot1, prot2 = prot2, prot3 = prot3)
  
  # storage
  matrix_data <- lapply(names(matrices), function(mat_name) {
    mat <- matrices[[mat_name]]
    protA_orig <- as.numeric(mat[protA_name, ])
    protB_orig <- as.numeric(mat[protB_name, ])
    
    # 100 perturbation generation
    samples <- lapply(1:100, function(x) {
      generate_perturbed_sample(protA_orig, protB_orig, mat_name)
    })
    
    list(
      original_A = protA_orig,
      original_B = protB_orig,
      perturbed = samples
    )
  })
  names(matrix_data) <- names(matrices)
  
  perturbed_samples[[i]] <- list(
    protein_A = protA_name,
    protein_B = protB_name,
    matrices = matrix_data
  )
  # process output
  cat(sprintf("Generated samples for pair %d: %s vs %s across all matrices\n", 
              i, protA_name, protB_name))
}

# list to data frame
perturbed_df <- map_dfr(1:length(perturbed_samples), function(pair_idx) {
  pair_data <- perturbed_samples[[pair_idx]]
  
  # three matrices
  map_dfr(names(pair_data$matrices), function(mat_name) {
    mat_data <- pair_data$matrices[[mat_name]]
    
    # pick original firstly
    orig_df <- data.frame(
      pair_id = pair_idx,
      protein_A = pair_data$protein_A,
      protein_B = pair_data$protein_B,
      matrix = mat_name,
      sample_type = "original",
      sample_id = 0,  # original == 0
      timepoint = seq_along(mat_data$original_A),
      value_A = mat_data$original_A,
      value_B = mat_data$original_B)
    
    # pick perturbation secondly
    perturbed_df <- map_dfr(1:length(mat_data$perturbed), function(sample_idx) {
      perturbed <- mat_data$perturbed[[sample_idx]]
      data.frame(
        pair_id = pair_idx,
        protein_A = pair_data$protein_A,
        protein_B = pair_data$protein_B,
        matrix = mat_name,
        sample_type = "perturbed",
        sample_id = sample_idx,  # perturbation ID
        timepoint = seq_along(perturbed$protA),
        value_A = perturbed$protA,
        value_B = perturbed$protB
      )
    })
    
    # comb original and perturbation
    bind_rows(orig_df, perturbed_df)
  })
})

# check data format
str(perturbed_df)
head(perturbed_df)

# # long format to wide format
# perturbed_wide <- perturbed_df %>%
#   pivot_longer(
#     cols = c(value_A, value_B),
#     names_to = "protein_role",
#     values_to = "intensity"
#   ) %>%
#   mutate(protein = ifelse(protein_role == "value_A", protein_A, protein_B)) %>%
#   select(-protein_role, -protein_A, -protein_B)
# 
# # pick 
perturbed_only <- perturbed_df %>%
  filter(sample_type == "perturbed")

# protein_A and protein_B fill up value_A and value_B
# protein_A == value_A
proteinA_matrix <- perturbed_only %>%
  mutate(
    row_name = paste(protein_A, sample_id, matrix, sep = "_"),
    timepoint = as.character(timepoint)
  ) %>%
  select(row_name, timepoint, intensity = value_A) %>%
  # check the timepoint (fraction 1 to 60)
  group_by(row_name, timepoint) %>%
  summarise(intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  # long to wide
  pivot_wider(
    names_from = timepoint,
    values_from = intensity
  ) %>%
  column_to_rownames("row_name")

# protein_B == value_B
proteinB_matrix <- perturbed_only %>%
  mutate(
    row_name = paste(protein_B, sample_id, matrix, sep = "_"),
    timepoint = as.character(timepoint)
  ) %>%
  select(row_name, timepoint, intensity = value_B) %>%
  group_by(row_name, timepoint) %>%
  summarise(intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = timepoint,
    values_from = intensity
  ) %>%
  column_to_rownames("row_name")

# comb protein_A and protein_B (protein ID need to grep"_")
final_matrix <- rbind(proteinA_matrix, proteinB_matrix)
# manipulate NaN to NA
final_matrix2 <- replace_nan_inf(final_matrix)
final_matrix2_with_id <- cbind(ID = rownames(final_matrix2), final_matrix2)
rio::export(perturbed_df,"TPC_1v5_allpairs_20250529.tsv")
rio::export(final_matrix2_with_id ,"TPC_1v5_perturbation_20250529.tsv")
