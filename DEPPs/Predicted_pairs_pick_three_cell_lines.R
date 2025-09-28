##################################################
## Project: NC rebuttal
## Script purpose: pick NA protein pairs in PrInCE
## Date: 2025-06-03
## Author: MenggeLYU
## Version: 1.0
##################################################


# path and library --------------------------------------------------------
library(rio)
library(data.table)
library(dplyr)
library(tidyverse)
library(tidyr)
library(PrInCE)

# label NA
## TPC_predicted pairs -----------------------------------------------------------
load("E:/MenggeLYU/NC_rebuttal_t3/New_feature_RF/final_ratio/TPC_basevalues_MGL_20250509.RData")
Pair <- TPC_features[,1:2]
label <- as.data.frame(TPC_label)
Pairwithlabel <- cbind(Pair,label)
Pairwithlabel_NA <- Pairwithlabel[-which(Pairwithlabel$TPC_label ==1),]
Pairwithlabel_NA <- Pairwithlabel_NA[-which(Pairwithlabel_NA$TPC_label ==0),]
rio::export(Pairwithlabel_NA, "E:/MenggeLYU/NC_rebuttal_t3/New_feature_RF/ComplexXGB/TPC_predicted_pairs_20250616.tsv")
rm(list=ls())
gc()

## FTC238_predicted pairs --------------------------------------------------
load("E:/MenggeLYU/NC_rebuttal_t3/New_feature_RF/final_ratio/FTC238_basevalues_MGL_20250509.RData")
Pair <- FTC238_features[,1:2]
label <- as.data.frame(FTC238_label)
Pairwithlabel <- cbind(Pair,label)
Pairwithlabel_NA <- Pairwithlabel[-which(Pairwithlabel$FTC238_label ==1),]
Pairwithlabel_NA <- Pairwithlabel_NA[-which(Pairwithlabel_NA$FTC238_label ==0),]
rio::export(Pairwithlabel_NA, "E:/MenggeLYU/NC_rebuttal_t3/New_feature_RF/ComplexXGB/FTC238_predicted_pairs_20250616.tsv")
rm(list=ls())
gc()

## Nthy_predicted pairs --------------------------------------------------
load("E:/MenggeLYU/NC_rebuttal_t3/New_feature_RF/final_ratio/Nthy_basevalues_MGL_20250509.RData")
Pair <- Nthy_features[,1:2]
label <- as.data.frame(Nthy_label)
Pairwithlabel <- cbind(Pair,label)
Pairwithlabel_NA <- Pairwithlabel[-which(Pairwithlabel$Nthy_label ==1),]
Pairwithlabel_NA <- Pairwithlabel_NA[-which(Pairwithlabel_NA$Nthy_label ==0),]
rio::export(Pairwithlabel_NA, "E:/MenggeLYU/NC_rebuttal_t3/New_feature_RF/ComplexXGB/Nthy_predicted_pairs_20250616.tsv")
gc()

# label = 1 and 0
## TPC_predicted pairs -----------------------------------------------------------
load("E:/MenggeLYU/NC_rebuttal_t3/New_feature_RF/final_ratio/TPC_basevalues_MGL_20250509.RData")
Pair <- TPC_features[,1:2]
label <- as.data.frame(TPC_label)
Pairwithlabel <- cbind(Pair,label)
Pairwithlabel_01 <- Pairwithlabel[which(Pairwithlabel$TPC_label == 1|Pairwithlabel$TPC_label == 0),]

rio::export(Pairwithlabel_01, "E:/MenggeLYU/NC_rebuttal_t3/New_feature_RF/ComplexXGB/TPC_label01_pairs_20250616.tsv")
rm(list=ls())
gc()

## FTC238_predicted pairs --------------------------------------------------
load("E:/MenggeLYU/NC_rebuttal_t3/New_feature_RF/final_ratio/FTC238_basevalues_MGL_20250509.RData")
Pair <- FTC238_features[,1:2]
label <- as.data.frame(FTC238_label)
Pairwithlabel <- cbind(Pair,label)
Pairwithlabel_01 <- Pairwithlabel[which(Pairwithlabel$FTC238_label == 1|Pairwithlabel$FTC238_label == 0),]
rio::export(Pairwithlabel_01, "E:/MenggeLYU/NC_rebuttal_t3/New_feature_RF/ComplexXGB/FTC238_label01_pairs_20250616.tsv")
rm(list=ls())
gc()

## Nthy_predicted pairs --------------------------------------------------
load("E:/MenggeLYU/NC_rebuttal_t3/New_feature_RF/final_ratio/Nthy_basevalues_MGL_20250509.RData")
Pair <- Nthy_features[,1:2]
label <- as.data.frame(Nthy_label)
Pairwithlabel <- cbind(Pair,label)
Pairwithlabel_01 <- Pairwithlabel[which(Pairwithlabel$Nthy_label == 1|Pairwithlabel$Nthy_label == 0),]
rio::export(Pairwithlabel_01, "E:/MenggeLYU/NC_rebuttal_t3/New_feature_RF/ComplexXGB/Nthy_label01_pairs_20250616.tsv")
rm(list=ls())
gc()
