##################################################
## Project: NC_rebuttal
## Script purpose: pick true samples
## Date: 2025-05-27
## Author: MenggeLYU
## Version: 1.5.0
##################################################


# path and library --------------------------------------------------------
setwd("E:\\final_ratio\\TP_distribution")
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)

# Nthy --------------------------------------------------------------------
nthy_results <- fread("E:\\Nthy\\final_predictions_Nthy.csv")
onefive <- nthy_results[which(nthy_results$ratio=="1:5"),]
oneten <- nthy_results[which(nthy_results$ratio=="1:10"),]
onefivetrue <- onefive[which(onefive$true_label=="1"),]
onetentrue <- oneten[which(oneten$true_label=="1"),]
combined_data <- rbind(
  transform(onetentrue, ratio = "1:10"),
  transform(onefivetrue, ratio = "1:15")
)


ggplot(combined_data, aes(x = predicted_score, fill = ratio)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Predicted Score",
       x = "Predicted Score",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("1:10" = "lightblue", "1:15" = "pink"),
                    labels = c("1:10" = "1:10", "1:15" = "1:15"),
                    name = "Ratio")
ggsave("Nthy_tp_distribution.pdf", width = 8, height = 6)

onefivepick <- onefivetrue[which(onefivetrue$predicted_score>0.5),]
onetenpick <- onetentrue[which(onetentrue$predicted_score>0.5),]

load("E:\\Nthy\\optimized_analysis_final_Nthy.RData")

onefiverows <- new_features[onefivepick$sample_id, ]
onetenrows <- new_features[onetenpick$sample_id, ]

rio::export(onefiverows, "Nthy_1v5_score05_proteinpairs_features_MGL.tsv")
rio::export(onetenrows, "Nthy_1v10_score05_proteinpairs_features_MGL.tsv")
rm(list=ls())
gc()

# TPC --------------------------------------------------------------------
tpc_results <- fread("E:\\TPC\\final_predictions_TPC.csv")
onefive <- tpc_results[which(tpc_results$ratio=="1:5"),]
oneten <- tpc_results[which(tpc_results$ratio=="1:10"),]
onefivetrue <- onefive[which(onefive$true_label=="1"),]
onetentrue <- oneten[which(oneten$true_label=="1"),]
combined_data <- rbind(
  transform(onetentrue, ratio = "1:10"),
  transform(onefivetrue, ratio = "1:15")
)


ggplot(combined_data, aes(x = predicted_score, fill = ratio)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Predicted Score",
       x = "Predicted Score",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("1:10" = "lightblue", "1:15" = "pink"),
                    labels = c("1:10" = "1:10", "1:15" = "1:15"),
                    name = "Ratio")
ggsave("TPC_tp_distribution.pdf", width = 8, height = 6)

onefivepick <- onefivetrue[which(onefivetrue$predicted_score>0.5),]
onetenpick <- onetentrue[which(onetentrue$predicted_score>0.5),]

load("E:\\TPC\\optimized_analysis_final_TPC.RData")

onefiverows <- new_features[onefivepick$sample_id, ]
onetenrows <- new_features[onetenpick$sample_id, ]

rio::export(onefiverows, "TPC_1v5_score05_proteinpairs_features_MGL.tsv")
rio::export(onetenrows, "TPC_1v10_score05_proteinpairs_features_MGL.tsv")
rm(list=ls())
gc()

# FTC238 --------------------------------------------------------------------
ftc238_results <- fread("E:\\FTC238\\final_predictions_FTC238.csv")
onefive <- ftc238_results[which(ftc238_results$ratio=="1:5"),]
oneten <- ftc238_results[which(ftc238_results$ratio=="1:10"),]
onefivetrue <- onefive[which(onefive$true_label=="1"),]
onetentrue <- oneten[which(oneten$true_label=="1"),]
combined_data <- rbind(
  transform(onetentrue, ratio = "1:10"),
  transform(onefivetrue, ratio = "1:5")
)


ggplot(combined_data, aes(x = predicted_score, fill = ratio)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Predicted Score",
       x = "Predicted Score",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("1:10" = "lightblue", "1:5" = "pink"),
                    labels = c("1:10" = "1:10", "1:5" = "1:5"),
                    name = "Ratio")
ggsave("FTC238_tp_distribution.pdf", width = 8, height = 6)

onefivepick <- onefivetrue[which(onefivetrue$predicted_score>0.5),]
onetenpick <- onetentrue[which(onetentrue$predicted_score>0.5),]

load("E:\\FTC238\\optimized_analysis_final_FTC238.RData")

onefiverows <- new_features[onefivepick$sample_id, ]
onetenrows <- new_features[onetenpick$sample_id, ]

rio::export(onefiverows, "FTC238_1v5_score05_proteinpairs_features_MGL.tsv")
rio::export(onetenrows, "FTC238_1v10_score05_proteinpairs_features_MGL.tsv")
rm(list=ls())
gc()

