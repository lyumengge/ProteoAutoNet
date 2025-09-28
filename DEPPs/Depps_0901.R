##################################################
## Project:NC_rebuttal
## Script purpose:TPC-1 vs Nthy-ori 3-1 results-depps
## Date: 2025-09-06
## Author: MenggeLYU
## Version: 1.1
##################################################



# path and library --------------------------------------------------------
setwd("G://Tower3//NC_rebuttal_t3//New_feature_RF//ComplexXGB")
library(data.table)
library(rio)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(pROC)

# load data ---------------------------------------------------------------

nf <- rio::import("G://Tower3//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//N_F_auc_0710//N_F_AUC_forDEPPS20250712.csv")
nt <- rio::import("G://Tower3//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//N_T_auc_0710//N_T_AUC_forDEPPS20250712.csv")

# Nthy vS TPC-1 - each replicate as a unique value
comb <- nt
comb2 <- unique(comb) # 80202

# protein auc average
filter1 <- comb2 %>%
  group_by(protein_pair, data_file, protein_label) %>%
  summarise(mean_sum = mean(auc, na.rm = TRUE), .groups = "drop") 
# protein pair auc sum
filter1 <- unique(filter1) # 80202
filter2 <- filter1 %>%
  group_by(protein_pair, data_file) %>%
  summarise(sum = sum(mean_sum, na.rm = TRUE), .groups = "drop") # 40101

filter3 <- unique(filter2)
# group
filter3$data_file <- gsub("_gaussian", "", filter3$data_file)

filter4 <- filter3 %>%
  mutate(cell = case_when(
    # data_file %in% c("TPCNO1", "TPCNO2", "TPCNO3") ~ "TPC",
    data_file %in% c("Nthy1", "Nthy2", "Nthy3") ~ "Nthy",
    # data_file %in% c("FTC133NO1", "FTC133NO2", "FTC133NO3") ~ "FTC133",
    data_file %in% c("TPC1", "TPC2", "TPC3") ~ "TPC",
    TRUE ~ data_file))

nt_results <- filter4 %>%
  group_by(protein_pair) %>%
  summarise(
    # Subset data for the two cell types
    Nthy_values = list(sum[cell == "Nthy"]),
    TPC_values = list(sum[cell == "TPC"]),
    # Calculate medians
    Nthy_median = median(unlist(Nthy_values), na.rm = TRUE),
    TPC_median = median(unlist(TPC_values), na.rm = TRUE),
    # Fold change (log2 scale)
    fold_change = log2(TPC_median / Nthy_median),
    # Likelihood ratio test
    p_value = ifelse(length(unlist(Nthy_values)) > 1 & length(unlist(TPC_values)) > 1,
                     {
                       # Get the current protein pair
                       current_pair <- unique(protein_pair)
                       # Filter data for the current protein pair
                       data_subset <- filter(filter3, protein_pair == current_pair)
                       # Fit models
                       null_model <- lm(sum ~ 1, data = data_subset)
                       alt_model <- lm(sum ~ cell, data = data_subset)
                       # Perform likelihood ratio test
                       anova(null_model, alt_model, test = "LRT")$`Pr(>Chi)`[2]
                     },
                     NA_real_)
  ) %>%
  # Adjust p-values using BH method
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  # Filter results with adjusted p-value < 0.05
  filter(!is.na(adj_p_value) & adj_p_value < 0.05) # 3211
results1 <- nt_results[which(nt_results$p_value<0.05),] # 3211
results2 <- results1[which(results1$fold_change>log2(1.5)|results1$fold_change < -log2(1.5)),] # 2198
rio::export(results2, "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\TPC_Nthy_depp_p005_log21_5_20250714.tsv")
results2pp <- as.data.frame(unique(results2$protein_pair)) # 2198
colnames(results2pp)[1] <- "combine"
rio::export(results2pp, "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\TPC_Nthy_depp_pair_list_20250714.tsv")

# ## TPC vs Nthy
# results <- nt_results[-which(nt_results$adj_p_value==0),]
# all three pairs ---------------------------------------------------------
N_Tcomb_matched <- rio::import("G://Tower3//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//N_T_pair_listwithlabel_0710.tsv")
N_Fcomb_matched <- rio::import("G://Tower3//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//N_F_pair_listwithlabel_0710.tsv")

allthree <- rbind(N_Tcomb_matched,N_Fcomb_matched)
allthree <- unique(allthree) # 10571
table(allthree$label) # 1150

allthree_withlabel <- rbind(nt,nf)
allthree_withlabel <- allthree_withlabel[,1:3]
allthree_withlabel <- unique(allthree_withlabel)
colnames(allthree_withlabel)[1] <- "combine"
allthree_withlabel2 <- left_join(allthree_withlabel,allthree)
allthree_withlabel2$cell <- gsub("[1-3]_gaussian","",allthree_withlabel2$data_file)
allthree_withlabel2$cell <- gsub("NO","",allthree_withlabel2$cell)
table(allthree_withlabel2$cell)
head(allthree_withlabel2)

comb <- rbind(nf,nt)
head(comb)
pairlength <- unique(comb$protein_pair) # 10571
head(allthree_withlabel2)
colnames(allthree_withlabel2)[1] <- "protein_pair"
allthree_withlabel3 <- left_join(allthree_withlabel2, comb, by = c("protein_pair","data_file","protein_label"))
allthree_withlabel3$col_num <- gsub('(.*)_','',allthree_withlabel3$max_col)
allthree_withlabel3 <- allthree_withlabel3 %>%
  mutate(col_num = as.numeric(col_num))  # mutate to numeric
allthree_withlabel4 <- allthree_withlabel3 %>%
  group_by(protein_pair, protein_label, cell) %>%
  filter(max(col_num) - min(col_num) <= 2) %>%
  ungroup()
pairlength <- unique(allthree_withlabel4$protein_pair) # 8984

# depps filter with 8000 pairs --------------------------------------------
depps_nt <-  nt_results[-which(nt_results$adj_p_value==0),]
depps_nt <- depps_nt[which(depps_nt$adj_p_value<0.05),]
depps_nt <- depps_nt[which(depps_nt$fold_change > log2(1.5) |depps_nt$fold_change< -log2(1.5) ),] #2183
depps_nt2 <- depps_nt[which(depps_nt$protein_pair %in% allthree_withlabel4$protein_pair==TRUE),] # 2014


# nf depps network --------------------------------------------------------
depps_nt2  <- depps_nt2  %>%
  separate(protein_pair, into = c("protein_1", "protein_2"), sep = "_")
protein_nt2 <- as.data.frame(unique(c(depps_nt2$protein_1,depps_nt2$protein_2)))
colnames(protein_nt2) <- 'prot'
protein_nt2 <- unique(protein_nt2) #1061
rio::export(protein_nt2, "G://Tower3//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//tpc_nthy_network_final_0714.tsv")
nt2_gene <- rio::import("G://Tower3//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//tpc_nthy_network_string_mapping_gene_0714.tsv")
protein_nt2$queryItem <- gsub("-(.*)",'', protein_nt2$prot)
protein_nt2 <- left_join(protein_nt2,nt2_gene, by = "queryItem")
protein_nt2 <- protein_nt2[,c(1,2,4,5)]
protein_nt2 <- unique(protein_nt2)
## string
nt2_egde <- rio::import('G://Tower3//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//tpc_nthy_network_string_mapping_pairs_0714.tsv')
nt2_edge_string <- unique(nt2_egde[,c(1:2)])
colnames(nt2_edge_string) <- c("node1","node2")
nt2_edge_string$edge <- "STRING"
nt2_edge_string$fold_change <- NA
## CF-MS
colnames(protein_nt2)[1] <- "protein_1"
depps_nt3 <- left_join(depps_nt2, protein_nt2, by = "protein_1")
colnames(depps_nt3)[12] <- "node1"
depps_nt3 <- depps_nt3[,c(1:9,12)]
colnames(protein_nt2)[1] <- "protein_2"
depps_nt4 <- left_join(depps_nt3, protein_nt2, by = "protein_2")
colnames(depps_nt4)[13] <- "node2"

depps_nt5 <- depps_nt4[,c(1:10,13)]
depps_nt5 <- depps_nt5[,c(10,11,7)]
depps_nt5$edge <- "CF-DIA-MS"

## nT CF-DIA-MS cluster
Nthy1 <- rio::import("G:\\Tower3\\Figure_NC\\Nthy_202410\\Nthy1_gaussian_withoutlog2.tsv")
Nthy2 <- rio::import("G:\\Tower3\\Figure_NC\\Nthy_202410\\Nthy2_gaussian_withoutlog2.tsv")
Nthy3 <- rio::import("G:\\Tower3\\Figure_NC\\Nthy_202410\\Nthy3_gaussian_withoutlog2.tsv")

TPC1 <- rio::import("G:\\Tower3\\Figure_NC\\TPC_202410\\TPC1_gaussian_withoutlog2.tsv")
TPC2 <- rio::import("G:\\Tower3\\Figure_NC\\TPC_202410\\TPC2_gaussian_withoutlog2.tsv")
TPC3 <- rio::import("G:\\Tower3\\Figure_NC\\TPC_202410\\TPC3_gaussian_withoutlog2.tsv")

combined <- bind_rows(Nthy1, Nthy2, Nthy3)
prot <- unique(c(depps_nt4$protein_1,depps_nt4$protein_2))

Nthy_mean <- combined %>%
  group_by(V1) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  

combined2 <- bind_rows(TPC1, TPC2, TPC3)
#prot <- unique(c(depps_nt4$protein_1,depps_nt4$protein_2))

TPC_mean <- combined2 %>%
  group_by(V1) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  

heatmap <- Nthy_mean[Nthy_mean$V1%in%prot,] # 1061
heatmap_clean <- heatmap[,2:61]
rownames(heatmap_clean) <- heatmap$V1
map1 <- pheatmap::pheatmap(heatmap_clean, scale = "row",
                           clustering_distance_rows = "euclidean", 
                           clustering_method = "ward.D2", 
                           cluster_rows = T,cluster_cols = F,
                           treeheight_row = 0,
                           show_rownames = FALSE)

heatmap2 <- TPC_mean[TPC_mean$V1%in%prot,] # 1061
heatmap_clean2 <- heatmap2[,2:61]
rownames(heatmap_clean2) <- heatmap2$V1
map2 <- pheatmap::pheatmap(heatmap_clean2, scale = "row",
                           clustering_distance_rows = "euclidean", 
                           clustering_method = "ward.D2",
                           cluster_rows = T,cluster_cols = F,
                           treeheight_row = 0,
                           show_rownames = FALSE)

## dendogram similarity
dend1 <- as.dendrogram(map1$tree_row)
dend2 <- as.dendrogram(map2$tree_row)

hclust1 <- as.hclust(map1$tree_row)
hclust2 <- as.hclust(map2$tree_row)
# plot(hclust1, cex = 0.6, labels = FALSE)  
# abline(h = 50, col = "red")  
# plot(hclust2, cex = 0.6)  
# abline(h = 50, col = "red")  

plot(hclust1, 
     labels = FALSE,  
     hang = -1,       
     cex = 0.6,
     main = "Nthy-ori 3-1 Cluster Dendrogram (Compared with TPC-1)")
abline(h = 50, col = "red")

plot(hclust2, 
     labels = FALSE,  
     hang = -1,       
     cex = 0.6,
     main = "TPC-1 Cluster Dendrogram (Compared with Nthy-ori 3-1)")
abline(h = 50, col = "red")

# fviz_nbclust(heatmap_clean, FUN = hcut, method = "wss")
# fviz_nbclust(heatmap_clean, FUN = hcut, method = "silhouette")

## cluster number
k <- 10
clusters1 <- cutree(hclust1, k = k)
clusters2 <- cutree(hclust2, k = k)
head(clusters1)
head(clusters2)

df_clusters <- data.frame(
  protein = names(clusters1),
  cluster_map1 = clusters1,
  cluster_map2 = clusters2)

overlap <- df_clusters %>%
  group_by(cluster_map1, cluster_map2) %>%
  summarise(
    n_proteins = n(),
    proteins = list(protein),
    .groups = "drop"
  ) %>%
  arrange(desc(n_proteins))

overlap_unnested <- overlap %>%
  unnest(proteins) %>%
  group_by(across(-proteins)) %>%  
  summarise(proteins = paste(proteins, collapse = ","), .groups = "drop")

colnames(overlap_unnested)[1] <- "Nthy-ori 3-1 clusters"
colnames(overlap_unnested)[2] <- "TPC-1 clusters"
colnames(overlap_unnested)[4] <- "protein list"

rio::export(overlap_unnested,"TPC-1vsNthy_overlapped_cluster_results.xlsx")
## tpc-1 cluster 6 VS nthy cluster 1 - nthy-------------------------------------------------------
current_row <- overlap_unnested[3,]
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster1vsTPC_cluster6_87prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster1vsTPC_cluster6_87prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("Nthy_cluster1vsTPC_cluster6_87prots_nthy_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot$comb <- paste(cluster_prot$preferredName, " (", cluster_prot$Prot, ")", sep = "")
colnames(cluster_prot)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot, by = "V1")

subset_data2 <- subset_data1[,2:61]
rownames(subset_data2) <- subset_data1$comb
p1 <- pheatmap(
  subset_data2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE) # 25-26 cols
row_order <- p1$tree_row$order
ordered_rows <- rownames(subset_data2)[row_order]
## A4 heatmap standard
pdf("Nthy_cluster1vsTPC_cluster6_87prots_nthy_ori.pdf", 
    width = 8.27,         
    height = 11.69)   

pheatmap(
  subset_data2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 23,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("Nthy_cluster1vsTPC_cluster6_87prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot$preferredName & node2 %in% cluster_prot$preferredName)
string_normalized <- string %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")))

cfms_marked <- cfms_interactions %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")),
    edge = ifelse(pair %in% string_normalized$pair, "Both", edge)
  ) %>%
  select(-pair)  
rio::export(cfms_marked,"Nthy_cluster1vsTPC_cluster6_87prots_nthy_final.tsv")

## tpc-1 cluster 6 VS nthy cluster 1 - tpc-------------------------------------------------------
# current_row <- overlap[3,]
# proteins_to_plot <- unlist(current_row$proteins)
# subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
# 
# current_row <- overlap_unnested[1,]
# proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster1vsTPC_cluster6_87prots_tpc.tsv")
rio::export(cluster_prot,"Nthy_cluster1vsTPC_cluster6_87prots_tpc.xlsx")



## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("Nthy_cluster1vsTPC_cluster6_87prots_nthy_gene.tsv") # tpc same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot$comb <- paste(cluster_prot$preferredName, " (", cluster_prot$Prot, ")", sep = "")
colnames(cluster_prot)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot, by = "V1")

subset_data2 <- subset_data1[,2:61]
rownames(subset_data2) <- subset_data1$comb
subset_data2_ordered <- subset_data2[ordered_rows, ]
rownames(subset_data2_ordered) <- ordered_rows
pheatmap(
  subset_data2_ordered,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE) # 10 cols

## A4 heatmap standard
pdf("Nthy_cluster1vsTPC_cluster6_87prots_tpc_ori.pdf", 
    width = 8.27,         
    height = 11.69)   

pheatmap(
  subset_data2_ordered,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 23,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - ftc238 is similar with Nthy
string <- rio::import("Nthy_cluster1vsTPC_cluster6_87prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot$preferredName & node2 %in% cluster_prot$preferredName)
string_normalized <- string %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")))

cfms_marked <- cfms_interactions %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")),
    edge = ifelse(pair %in% string_normalized$pair, "Both", edge)
  ) %>%
  select(-pair)  
rio::export(cfms_marked,"Nthy_cluster1vsTPC_cluster6_87prots_tpc_final.tsv")

### tpc vs nthy nodes mean intensity ----------------------------------------------
node1 <- cfms_marked[,1]
node2 <- cfms_marked[,2]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node <- unique(cf_node) # 32
colnames(cf_node)[1] <- "preferredName" 
cf_node <- left_join(cf_node, cluster_prot, by = "preferredName") # 45

calculate_curve_auc <- function(protein_matrix, protein_list) {
  protein_names <- protein_matrix$V1
  expr_matrix <- as.matrix(protein_matrix[, -1])
  rownames(expr_matrix) <- protein_names
  target_proteins <- expr_matrix[rownames(expr_matrix) %in% protein_list, , drop = FALSE]
  if (nrow(target_proteins) == 0) {
    return(data.frame(protein = character(), auc = numeric()))
  }
  auc_results <- apply(target_proteins, 1, function(expr_values) {
    x <- 1:length(expr_values)
    y <- as.numeric(expr_values)
    auc_value <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    return(auc_value)
  })
  return(data.frame(protein = names(auc_results), auc = auc_results))
}

protein_list <- cf_node$V1

Nthy1_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy1_gaussian_withoutlog2_normalized.tsv")
Nthy2_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy2_gaussian_withoutlog2_normalized.tsv")
Nthy3_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy3_gaussian_withoutlog2_normalized.tsv")

TPC1_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC1_gaussian_withoutlog2_normalized.tsv")
TPC2_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC2_gaussian_withoutlog2_normalized.tsv")
TPC3_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC3_gaussian_withoutlog2_normalized.tsv")

combined_nor <- bind_rows(Nthy1_nor, Nthy2_nor, Nthy3_nor)
Nthy_nor_mean <- combined_nor %>%
  group_by(V1) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  

combined2_nor <- bind_rows(TPC1_nor, TPC2_nor, TPC3_nor)
#prot <- unique(c(depps_nt4$protein_1,depps_nt4$protein_2))

TPC_nor_mean <- combined2_nor %>%
  group_by(V1) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  


auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list)
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") # 78
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster1vsTPC_cluster6_87prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(TPC_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","TPC")
cf_node_tpc <- left_join(cf_node, auc_results2, by = "V1") #45
#cf_node_tpc$tpc <- log10(cf_node_tpc$tpc)
colnames(cf_node_tpc)[1] <- "name"
rio::export(cf_node_tpc, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster1vsTPC_cluster6_87prots_tpc_intensity.tsv")

## tpc-1 cluster 2 VS nthy cluster 2 - nthy-------------------------------------------------------
current_row <- overlap_unnested[6,]
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster2vsTPC_cluster2_46prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster2vsTPC_cluster2_46prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("Nthy_cluster2vsTPC_cluster2_46prots_nthy_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot$comb <- paste(cluster_prot$preferredName, " (", cluster_prot$Prot, ")", sep = "")
colnames(cluster_prot)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot, by = "V1")

subset_data2 <- subset_data1[,2:61]
rownames(subset_data2) <- subset_data1$comb
p1 <- pheatmap(
  subset_data2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE) # 25-26 cols
row_order <- p1$tree_row$order
ordered_rows <- rownames(subset_data2)[row_order]
## A4 heatmap standard
pdf("Nthy_cluster2vsTPC_cluster2_46prots_nthy_ori.pdf", 
    width = 8.27,         
    height = 11.69)   

pheatmap(
  subset_data2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 23,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("Nthy_cluster2vsTPC_cluster2_46prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot$preferredName & node2 %in% cluster_prot$preferredName)
string_normalized <- string %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")))

cfms_marked <- cfms_interactions %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")),
    edge = ifelse(pair %in% string_normalized$pair, "Both", edge)
  ) %>%
  select(-pair)  
rio::export(cfms_marked,"Nthy_cluster2vsTPC_cluster2_46prots_nthy_final.tsv")

## tpc-1 cluster 2 VS nthy cluster 2 - tpc-------------------------------------------------------
# current_row <- overlap[3,]
# proteins_to_plot <- unlist(current_row$proteins)
# subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
# 
# current_row <- overlap_unnested[1,]
# proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster2vsTPC_cluster2_46prots_tpc.tsv")
rio::export(cluster_prot,"Nthy_cluster2vsTPC_cluster2_46prots_tpc.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("Nthy_cluster2vsTPC_cluster2_46prots_nthy_gene.tsv") # tpc same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot$comb <- paste(cluster_prot$preferredName, " (", cluster_prot$Prot, ")", sep = "")
colnames(cluster_prot)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot, by = "V1")

subset_data2 <- subset_data1[,2:61]
rownames(subset_data2) <- subset_data1$comb
subset_data2_ordered <- subset_data2[ordered_rows, ]
rownames(subset_data2_ordered) <- ordered_rows
pheatmap(
  subset_data2_ordered,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE) # 10 cols

## A4 heatmap standard
pdf("Nthy_cluster2vsTPC_cluster2_46prots_tpc_ori.pdf", 
    width = 8.27,         
    height = 11.69)   

pheatmap(
  subset_data2_ordered,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 23,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - ftc238 is similar with Nthy
string <- rio::import("Nthy_cluster2vsTPC_cluster2_46prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot$preferredName & node2 %in% cluster_prot$preferredName)
string_normalized <- string %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")))

cfms_marked <- cfms_interactions %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")),
    edge = ifelse(pair %in% string_normalized$pair, "Both", edge)
  ) %>%
  select(-pair)  
rio::export(cfms_marked,"Nthy_cluster2vsTPC_cluster2_46prots_tpc_final.tsv")

### tpc vs nthy nodes mean intensity ----------------------------------------------
node1 <- cfms_marked[,1]
node2 <- cfms_marked[,2]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node <- unique(cf_node) # 32
colnames(cf_node)[1] <- "preferredName" 
cf_node <- left_join(cf_node, cluster_prot, by = "preferredName") # 45

calculate_curve_auc <- function(protein_matrix, protein_list) {
  protein_names <- protein_matrix$V1
  expr_matrix <- as.matrix(protein_matrix[, -1])
  rownames(expr_matrix) <- protein_names
  target_proteins <- expr_matrix[rownames(expr_matrix) %in% protein_list, , drop = FALSE]
  if (nrow(target_proteins) == 0) {
    return(data.frame(protein = character(), auc = numeric()))
  }
  auc_results <- apply(target_proteins, 1, function(expr_values) {
    x <- 1:length(expr_values)
    y <- as.numeric(expr_values)
    auc_value <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    return(auc_value)
  })
  return(data.frame(protein = names(auc_results), auc = auc_results))
}

protein_list <- cf_node$V1

# Nthy1_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy1_gaussian_withoutlog2_normalized.tsv")
# Nthy2_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy2_gaussian_withoutlog2_normalized.tsv")
# Nthy3_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy3_gaussian_withoutlog2_normalized.tsv")
# 
# TPC1_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC1_gaussian_withoutlog2_normalized.tsv")
# TPC2_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC2_gaussian_withoutlog2_normalized.tsv")
# TPC3_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC3_gaussian_withoutlog2_normalized.tsv")
# 
# combined_nor <- bind_rows(Nthy1_nor, Nthy2_nor, Nthy3_nor)
# Nthy_nor_mean <- combined_nor %>%
#   group_by(V1) %>%
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  
# 
# combined2_nor <- bind_rows(TPC1_nor, TPC2_nor, TPC3_nor)
# #prot <- unique(c(depps_nt4$protein_1,depps_nt4$protein_2))
# 
# TPC_nor_mean <- combined2_nor %>%
#   group_by(V1) %>%
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  


auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list)
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") # 78
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster2vsTPC_cluster2_46prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(TPC_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","TPC")
cf_node_tpc <- left_join(cf_node, auc_results2, by = "V1") #45
#cf_node_tpc$tpc <- log10(cf_node_tpc$tpc)
colnames(cf_node_tpc)[1] <- "name"
rio::export(cf_node_tpc, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster2vsTPC_cluster2_46prots_tpc_intensity.tsv")



## tpc-1 cluster 2 VS nthy cluster 2 - nthy-------------------------------------------------------
current_row <- overlap_unnested[6,]
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster2vsTPC_cluster2_46prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster2vsTPC_cluster2_46prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("Nthy_cluster2vsTPC_cluster2_46prots_nthy_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot$comb <- paste(cluster_prot$preferredName, " (", cluster_prot$Prot, ")", sep = "")
colnames(cluster_prot)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot, by = "V1")

subset_data2 <- subset_data1[,2:61]
rownames(subset_data2) <- subset_data1$comb
p1 <- pheatmap(
  subset_data2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE) # 25-26 cols
row_order <- p1$tree_row$order
ordered_rows <- rownames(subset_data2)[row_order]
## A4 heatmap standard
pdf("Nthy_cluster2vsTPC_cluster2_46prots_nthy_ori.pdf", 
    width = 8.27,         
    height = 11.69)   

pheatmap(
  subset_data2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 23,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("Nthy_cluster2vsTPC_cluster2_46prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot$preferredName & node2 %in% cluster_prot$preferredName)
string_normalized <- string %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")))

cfms_marked <- cfms_interactions %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")),
    edge = ifelse(pair %in% string_normalized$pair, "Both", edge)
  ) %>%
  select(-pair)  
rio::export(cfms_marked,"Nthy_cluster2vsTPC_cluster2_46prots_nthy_final.tsv")

## tpc-1 cluster 2 VS nthy cluster 2 - tpc-------------------------------------------------------
# current_row <- overlap[3,]
# proteins_to_plot <- unlist(current_row$proteins)
# subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
# 
# current_row <- overlap_unnested[1,]
# proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster2vsTPC_cluster2_46prots_tpc.tsv")
rio::export(cluster_prot,"Nthy_cluster2vsTPC_cluster2_46prots_tpc.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("Nthy_cluster2vsTPC_cluster2_46prots_nthy_gene.tsv") # tpc same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot$comb <- paste(cluster_prot$preferredName, " (", cluster_prot$Prot, ")", sep = "")
colnames(cluster_prot)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot, by = "V1")

subset_data2 <- subset_data1[,2:61]
rownames(subset_data2) <- subset_data1$comb
subset_data2_ordered <- subset_data2[ordered_rows, ]
rownames(subset_data2_ordered) <- ordered_rows
pheatmap(
  subset_data2_ordered,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE) # 10 cols

## A4 heatmap standard
pdf("Nthy_cluster2vsTPC_cluster2_46prots_tpc_ori.pdf", 
    width = 8.27,         
    height = 11.69)   

pheatmap(
  subset_data2_ordered,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 23,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - ftc238 is similar with Nthy
string <- rio::import("Nthy_cluster2vsTPC_cluster2_46prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot$preferredName & node2 %in% cluster_prot$preferredName)
string_normalized <- string %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")))

cfms_marked <- cfms_interactions %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")),
    edge = ifelse(pair %in% string_normalized$pair, "Both", edge)
  ) %>%
  select(-pair)  
rio::export(cfms_marked,"Nthy_cluster2vsTPC_cluster2_46prots_tpc_final.tsv")

### tpc vs nthy nodes mean intensity ----------------------------------------------
node1 <- cfms_marked[,1]
node2 <- cfms_marked[,2]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node <- unique(cf_node) # 32
colnames(cf_node)[1] <- "preferredName" 
cf_node <- left_join(cf_node, cluster_prot, by = "preferredName") # 45

calculate_curve_auc <- function(protein_matrix, protein_list) {
  protein_names <- protein_matrix$V1
  expr_matrix <- as.matrix(protein_matrix[, -1])
  rownames(expr_matrix) <- protein_names
  target_proteins <- expr_matrix[rownames(expr_matrix) %in% protein_list, , drop = FALSE]
  if (nrow(target_proteins) == 0) {
    return(data.frame(protein = character(), auc = numeric()))
  }
  auc_results <- apply(target_proteins, 1, function(expr_values) {
    x <- 1:length(expr_values)
    y <- as.numeric(expr_values)
    auc_value <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    return(auc_value)
  })
  return(data.frame(protein = names(auc_results), auc = auc_results))
}

protein_list <- cf_node$V1

# Nthy1_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy1_gaussian_withoutlog2_normalized.tsv")
# Nthy2_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy2_gaussian_withoutlog2_normalized.tsv")
# Nthy3_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy3_gaussian_withoutlog2_normalized.tsv")
# 
# TPC1_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC1_gaussian_withoutlog2_normalized.tsv")
# TPC2_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC2_gaussian_withoutlog2_normalized.tsv")
# TPC3_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC3_gaussian_withoutlog2_normalized.tsv")
# 
# combined_nor <- bind_rows(Nthy1_nor, Nthy2_nor, Nthy3_nor)
# Nthy_nor_mean <- combined_nor %>%
#   group_by(V1) %>%
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  
# 
# combined2_nor <- bind_rows(TPC1_nor, TPC2_nor, TPC3_nor)
# #prot <- unique(c(depps_nt4$protein_1,depps_nt4$protein_2))
# 
# TPC_nor_mean <- combined2_nor %>%
#   group_by(V1) %>%
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  


auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list)
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") # 78
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster2vsTPC_cluster2_46prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(TPC_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","TPC")
cf_node_tpc <- left_join(cf_node, auc_results2, by = "V1") #45
#cf_node_tpc$tpc <- log10(cf_node_tpc$tpc)
colnames(cf_node_tpc)[1] <- "name"
rio::export(cf_node_tpc, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster2vsTPC_cluster2_46prots_tpc_intensity.tsv")

# tpc-1 cluster 1 VS nthy cluster 5 - nthy-------------------------------------------------------
current_row <- overlap_unnested[22,]
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster5vsTPC_cluster1_80prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster5vsTPC_cluster1_80prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("Nthy_cluster5vsTPC_cluster1_80prots_nthy_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot$comb <- paste(cluster_prot$preferredName, " (", cluster_prot$Prot, ")", sep = "")
colnames(cluster_prot)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot, by = "V1")

subset_data2 <- subset_data1[,2:61]
rownames(subset_data2) <- subset_data1$comb
p1 <- pheatmap(
  subset_data2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE) # 25-26 cols
row_order <- p1$tree_row$order
ordered_rows <- rownames(subset_data2)[row_order]
## A4 heatmap standard
pdf("Nthy_cluster5vsTPC_cluster1_80prots_nthy_ori.pdf", 
    width = 8.27,         
    height = 11.69)   

pheatmap(
  subset_data2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 23,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("Nthy_cluster5vsTPC_cluster1_80prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot$preferredName & node2 %in% cluster_prot$preferredName)
string_normalized <- string %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")))

cfms_marked <- cfms_interactions %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")),
    edge = ifelse(pair %in% string_normalized$pair, "Both", edge)
  ) %>%
  select(-pair)  
rio::export(cfms_marked,"Nthy_cluster5vsTPC_cluster1_80prots_nthy_final.tsv")

## tpc-1 cluster 1 VS nthy cluster 5 - tpc-------------------------------------------------------
# current_row <- overlap[3,]
# proteins_to_plot <- unlist(current_row$proteins)
# subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
# 
# current_row <- overlap_unnested[1,]
# proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster5vsTPC_cluster1_80prots_tpc.tsv")
rio::export(cluster_prot,"Nthy_cluster5vsTPC_cluster1_80prots_tpc.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("Nthy_cluster5vsTPC_cluster1_80prots_nthy_gene.tsv") # tpc same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot$comb <- paste(cluster_prot$preferredName, " (", cluster_prot$Prot, ")", sep = "")
colnames(cluster_prot)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot, by = "V1")

subset_data2 <- subset_data1[,2:61]
rownames(subset_data2) <- subset_data1$comb
subset_data2_ordered <- subset_data2[ordered_rows, ]
rownames(subset_data2_ordered) <- ordered_rows
pheatmap(
  subset_data2_ordered,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE) # 10 cols

## A4 heatmap standard
pdf("Nthy_cluster5vsTPC_cluster1_80prots_tpc_ori.pdf", 
    width = 8.27,         
    height = 11.69)   

pheatmap(
  subset_data2_ordered,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 23,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - ftc238 is similar with Nthy
string <- rio::import("Nthy_cluster5vsTPC_cluster1_80prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot$preferredName & node2 %in% cluster_prot$preferredName)
string_normalized <- string %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")))

cfms_marked <- cfms_interactions %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")),
    edge = ifelse(pair %in% string_normalized$pair, "Both", edge)
  ) %>%
  select(-pair)  
rio::export(cfms_marked,"Nthy_cluster5vsTPC_cluster1_80prots_tpc_final.tsv")

### tpc vs nthy nodes mean intensity ----------------------------------------------
node1 <- cfms_marked[,1]
node2 <- cfms_marked[,2]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node <- unique(cf_node) # 32
colnames(cf_node)[1] <- "preferredName" 
cf_node <- left_join(cf_node, cluster_prot, by = "preferredName") # 45

calculate_curve_auc <- function(protein_matrix, protein_list) {
  protein_names <- protein_matrix$V1
  expr_matrix <- as.matrix(protein_matrix[, -1])
  rownames(expr_matrix) <- protein_names
  target_proteins <- expr_matrix[rownames(expr_matrix) %in% protein_list, , drop = FALSE]
  if (nrow(target_proteins) == 0) {
    return(data.frame(protein = character(), auc = numeric()))
  }
  auc_results <- apply(target_proteins, 1, function(expr_values) {
    x <- 1:length(expr_values)
    y <- as.numeric(expr_values)
    auc_value <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    return(auc_value)
  })
  return(data.frame(protein = names(auc_results), auc = auc_results))
}

protein_list <- cf_node$V1

# Nthy1_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy1_gaussian_withoutlog2_normalized.tsv")
# Nthy2_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy2_gaussian_withoutlog2_normalized.tsv")
# Nthy3_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy3_gaussian_withoutlog2_normalized.tsv")
# 
# TPC1_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC1_gaussian_withoutlog2_normalized.tsv")
# TPC2_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC2_gaussian_withoutlog2_normalized.tsv")
# TPC3_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC3_gaussian_withoutlog2_normalized.tsv")
# 
# combined_nor <- bind_rows(Nthy1_nor, Nthy2_nor, Nthy3_nor)
# Nthy_nor_mean <- combined_nor %>%
#   group_by(V1) %>%
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  
# 
# combined2_nor <- bind_rows(TPC1_nor, TPC2_nor, TPC3_nor)
# #prot <- unique(c(depps_nt4$protein_1,depps_nt4$protein_2))
# 
# TPC_nor_mean <- combined2_nor %>%
#   group_by(V1) %>%
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  


auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list)
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") # 78
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster5vsTPC_cluster1_80prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(TPC_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","TPC")
cf_node_tpc <- left_join(cf_node, auc_results2, by = "V1") #45
#cf_node_tpc$tpc <- log10(cf_node_tpc$tpc)
colnames(cf_node_tpc)[1] <- "name"
rio::export(cf_node_tpc, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster5vsTPC_cluster1_80prots_tpc_intensity.tsv")

# tpc-1 cluster 4 VS nthy cluster 6 - nthy-------------------------------------------------------
current_row <- overlap_unnested[29,]
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster6vsTPC_cluster4_9prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster6vsTPC_cluster4_9prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("Nthy_cluster6vsTPC_cluster4_9prots_nthy_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot$comb <- paste(cluster_prot$preferredName, " (", cluster_prot$Prot, ")", sep = "")
colnames(cluster_prot)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot, by = "V1")

subset_data2 <- subset_data1[,2:61]
rownames(subset_data2) <- subset_data1$comb
p1 <- pheatmap(
  subset_data2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE) # 25-26 cols
row_order <- p1$tree_row$order
ordered_rows <- rownames(subset_data2)[row_order]
## A4 heatmap standard
pdf("Nthy_cluster6vsTPC_cluster4_9prots_nthy_ori.pdf", 
    width = 8.27,         
    height = 11.69)   

pheatmap(
  subset_data2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 5,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("Nthy_cluster6vsTPC_cluster4_9prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot$preferredName & node2 %in% cluster_prot$preferredName)
string_normalized <- string %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")))

cfms_marked <- cfms_interactions %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")),
    edge = ifelse(pair %in% string_normalized$pair, "Both", edge)
  ) %>%
  select(-pair)  
rio::export(cfms_marked,"Nthy_cluster6vsTPC_cluster4_9prots_nthy_final.tsv")

## tpc-1 cluster 4 VS nthy cluster 6 - tpc-------------------------------------------------------
# current_row <- overlap[3,]
# proteins_to_plot <- unlist(current_row$proteins)
# subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
# 
# current_row <- overlap_unnested[1,]
# proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster6vsTPC_cluster4_9prots_tpc.tsv")
rio::export(cluster_prot,"Nthy_cluster6vsTPC_cluster4_9prots_tpc.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("Nthy_cluster6vsTPC_cluster4_9prots_nthy_gene.tsv") # tpc same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot$comb <- paste(cluster_prot$preferredName, " (", cluster_prot$Prot, ")", sep = "")
colnames(cluster_prot)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot, by = "V1")

subset_data2 <- subset_data1[,2:61]
rownames(subset_data2) <- subset_data1$comb
subset_data2_ordered <- subset_data2[ordered_rows, ]
rownames(subset_data2_ordered) <- ordered_rows
pheatmap(
  subset_data2_ordered,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE) # 10 cols

## A4 heatmap standard
pdf("Nthy_cluster6vsTPC_cluster4_9prots_tpc_ori.pdf", 
    width = 8.27,         
    height = 11.69)   

pheatmap(
  subset_data2_ordered,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 5,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - ftc238 is similar with Nthy
string <- rio::import("Nthy_cluster6vsTPC_cluster4_9prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot$preferredName & node2 %in% cluster_prot$preferredName)
string_normalized <- string %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")))

cfms_marked <- cfms_interactions %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")),
    edge = ifelse(pair %in% string_normalized$pair, "Both", edge)
  ) %>%
  select(-pair)  
rio::export(cfms_marked,"Nthy_cluster6vsTPC_cluster4_9prots_tpc_final.tsv")

### tpc vs nthy nodes mean intensity ----------------------------------------------
node1 <- cfms_marked[,1]
node2 <- cfms_marked[,2]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node <- unique(cf_node) # 32
colnames(cf_node)[1] <- "preferredName" 
cf_node <- left_join(cf_node, cluster_prot, by = "preferredName") # 45

calculate_curve_auc <- function(protein_matrix, protein_list) {
  protein_names <- protein_matrix$V1
  expr_matrix <- as.matrix(protein_matrix[, -1])
  rownames(expr_matrix) <- protein_names
  target_proteins <- expr_matrix[rownames(expr_matrix) %in% protein_list, , drop = FALSE]
  if (nrow(target_proteins) == 0) {
    return(data.frame(protein = character(), auc = numeric()))
  }
  auc_results <- apply(target_proteins, 1, function(expr_values) {
    x <- 1:length(expr_values)
    y <- as.numeric(expr_values)
    auc_value <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    return(auc_value)
  })
  return(data.frame(protein = names(auc_results), auc = auc_results))
}

protein_list <- cf_node$V1

# Nthy1_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy1_gaussian_withoutlog2_normalized.tsv")
# Nthy2_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy2_gaussian_withoutlog2_normalized.tsv")
# Nthy3_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy3_gaussian_withoutlog2_normalized.tsv")
# 
# TPC1_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC1_gaussian_withoutlog2_normalized.tsv")
# TPC2_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC2_gaussian_withoutlog2_normalized.tsv")
# TPC3_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC3_gaussian_withoutlog2_normalized.tsv")
# 
# combined_nor <- bind_rows(Nthy1_nor, Nthy2_nor, Nthy3_nor)
# Nthy_nor_mean <- combined_nor %>%
#   group_by(V1) %>%
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  
# 
# combined2_nor <- bind_rows(TPC1_nor, TPC2_nor, TPC3_nor)
# #prot <- unique(c(depps_nt4$protein_1,depps_nt4$protein_2))
# 
# TPC_nor_mean <- combined2_nor %>%
#   group_by(V1) %>%
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  


auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list)
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") # 78
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster6vsTPC_cluster4_9prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(TPC_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","TPC")
cf_node_tpc <- left_join(cf_node, auc_results2, by = "V1") #45
#cf_node_tpc$tpc <- log10(cf_node_tpc$tpc)
colnames(cf_node_tpc)[1] <- "name"
rio::export(cf_node_tpc, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster6vsTPC_cluster4_9prots_tpc_intensity.tsv")


# tpc-1 cluster 2 VS nthy cluster 4 - nthy-------------------------------------------------------
current_row <- overlap_unnested[17,]
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster4vsTPC_cluster2_28prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster4vsTPC_cluster2_28prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("Nthy_cluster4vsTPC_cluster2_28prots_nthy_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot$comb <- paste(cluster_prot$preferredName, " (", cluster_prot$Prot, ")", sep = "")
colnames(cluster_prot)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot, by = "V1")

subset_data2 <- subset_data1[,2:61]
rownames(subset_data2) <- subset_data1$comb
p1 <- pheatmap(
  subset_data2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE) # 25-26 cols
row_order <- p1$tree_row$order
ordered_rows <- rownames(subset_data2)[row_order]
## A4 heatmap standard
pdf("Nthy_cluster4vsTPC_cluster2_28prots_nthy_ori.pdf", 
    width = 8.27,         
    height = 11.69)   

pheatmap(
  subset_data2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 5,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("Nthy_cluster4vsTPC_cluster2_28prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot$preferredName & node2 %in% cluster_prot$preferredName)
string_normalized <- string %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")))

cfms_marked <- cfms_interactions %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")),
    edge = ifelse(pair %in% string_normalized$pair, "Both", edge)
  ) %>%
  select(-pair)  
rio::export(cfms_marked,"Nthy_cluster4vsTPC_cluster2_28prots_nthy_final.tsv")

## tpc-1 cluster 2 VS nthy cluster 4 - tpc-------------------------------------------------------
# current_row <- overlap[3,]
# proteins_to_plot <- unlist(current_row$proteins)
# subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
# 
# current_row <- overlap_unnested[1,]
# proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster4vsTPC_cluster2_28prots_tpc.tsv")
rio::export(cluster_prot,"Nthy_cluster4vsTPC_cluster2_28prots_tpc.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("Nthy_cluster4vsTPC_cluster2_28prots_nthy_gene.tsv") # tpc same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot$comb <- paste(cluster_prot$preferredName, " (", cluster_prot$Prot, ")", sep = "")
colnames(cluster_prot)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot, by = "V1")

subset_data2 <- subset_data1[,2:61]
rownames(subset_data2) <- subset_data1$comb
subset_data2_ordered <- subset_data2[ordered_rows, ]
rownames(subset_data2_ordered) <- ordered_rows
pheatmap(
  subset_data2_ordered,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE) # 10 cols

## A4 heatmap standard
pdf("Nthy_cluster4vsTPC_cluster2_28prots_tpc_ori.pdf", 
    width = 8.27,         
    height = 11.69)   

pheatmap(
  subset_data2_ordered,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 5,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - ftc238 is similar with Nthy
string <- rio::import("Nthy_cluster4vsTPC_cluster2_28prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot$preferredName & node2 %in% cluster_prot$preferredName)
string_normalized <- string %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")))

cfms_marked <- cfms_interactions %>%
  mutate(
    pair = ifelse(node1 < node2, 
                  paste(node1, node2, sep = "|"), 
                  paste(node2, node1, sep = "|")),
    edge = ifelse(pair %in% string_normalized$pair, "Both", edge)
  ) %>%
  select(-pair)  
rio::export(cfms_marked,"Nthy_cluster4vsTPC_cluster2_28prots_tpc_final.tsv")

### tpc vs nthy nodes mean intensity ----------------------------------------------
node1 <- cfms_marked[,1]
node2 <- cfms_marked[,2]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node <- unique(cf_node) # 32
colnames(cf_node)[1] <- "preferredName" 
cf_node <- left_join(cf_node, cluster_prot, by = "preferredName") # 45

calculate_curve_auc <- function(protein_matrix, protein_list) {
  protein_names <- protein_matrix$V1
  expr_matrix <- as.matrix(protein_matrix[, -1])
  rownames(expr_matrix) <- protein_names
  target_proteins <- expr_matrix[rownames(expr_matrix) %in% protein_list, , drop = FALSE]
  if (nrow(target_proteins) == 0) {
    return(data.frame(protein = character(), auc = numeric()))
  }
  auc_results <- apply(target_proteins, 1, function(expr_values) {
    x <- 1:length(expr_values)
    y <- as.numeric(expr_values)
    auc_value <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    return(auc_value)
  })
  return(data.frame(protein = names(auc_results), auc = auc_results))
}

protein_list <- cf_node$V1

# Nthy1_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy1_gaussian_withoutlog2_normalized.tsv")
# Nthy2_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy2_gaussian_withoutlog2_normalized.tsv")
# Nthy3_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\Nthy3_gaussian_withoutlog2_normalized.tsv")
# 
# TPC1_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC1_gaussian_withoutlog2_normalized.tsv")
# TPC2_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC2_gaussian_withoutlog2_normalized.tsv")
# TPC3_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\TPC3_gaussian_withoutlog2_normalized.tsv")
# 
# combined_nor <- bind_rows(Nthy1_nor, Nthy2_nor, Nthy3_nor)
# Nthy_nor_mean <- combined_nor %>%
#   group_by(V1) %>%
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  
# 
# combined2_nor <- bind_rows(TPC1_nor, TPC2_nor, TPC3_nor)
# #prot <- unique(c(depps_nt4$protein_1,depps_nt4$protein_2))
# 
# TPC_nor_mean <- combined2_nor %>%
#   group_by(V1) %>%
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  


auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list)
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") # 78
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster4vsTPC_cluster2_28prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(TPC_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","TPC")
cf_node_tpc <- left_join(cf_node, auc_results2, by = "V1") #45
#cf_node_tpc$tpc <- log10(cf_node_tpc$tpc)
colnames(cf_node_tpc)[1] <- "name"
rio::export(cf_node_tpc, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster4vsTPC_cluster2_28prots_tpc_intensity.tsv")