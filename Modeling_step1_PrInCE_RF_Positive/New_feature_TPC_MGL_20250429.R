##################################################
## Project: NC_rebuttal
## Script purpose: eatures_generation_TPC
## Date: 2025-04-22
## Author: Mengge LYU
## Version: 1.4.1
##################################################

# library load --------------------------------------------------------------------
setwd("E:\\rebuttal_t3")
library(data.table)
library(dplyr)
library(PrInCE)
library(ggplot2)
library(jsonlite)
library(naivebayes)
# data --------------------------------------------------------------------
test1 <- read.delim("TPC1_protmatrix_20240918.tsv", row.names=1, quote="")
test2 <- read.delim("TPC2_protmatrix_20240918.tsv", row.names=1, quote="")
test3 <- read.delim("TPC3_protmatrix_20240918.tsv", row.names=1, quote="")

# gold standard -----------------------------------------------------------
load("three_comb_adj_matrix_20250409.RData")
gold_standard <- adj_df
colnames(test1) <- gsub("TPC.1_NO1_30minDIA_","", colnames(test1))
colnames(test2) <- gsub("TPC.1_NO2_30minDIA_","", colnames(test2))
colnames(test3) <- gsub("TPC.1_NO3_30minDIA_","", colnames(test3))
test1 <- test1[, order(as.numeric(colnames(test1)))]
test2 <- test2[, order(as.numeric(colnames(test2)))]
test3 <- test3[, order(as.numeric(colnames(test3)))]
test1 <- test1[,-c(1:3)]
test2 <- test2[,-c(1:3)]
test3 <- test3[,-c(1:3)]
test1 <- test1[rowSums(is.na(test1)) != ncol(test1), ]
test2 <- test2[rowSums(is.na(test2)) != ncol(test2), ]
test3 <- test3[rowSums(is.na(test3)) != ncol(test3), ]
matrices <- list(test1, test2, test3)
gaussian <- lapply(matrices, build_gaussians)
feature1 <- calculate_features(matrices[[1]], gaussian[[1]])
feature2 <- calculate_features(matrices[[2]], gaussian[[2]])
feature3 <- calculate_features(matrices[[3]], gaussian[[3]])
features <- concatenate_features(list(feature1, feature2, feature3)) %>% replace_missing_data()
label <- make_labels(gold_standard,features)
save.image("New_feature_TPC_20250429.Rdata")