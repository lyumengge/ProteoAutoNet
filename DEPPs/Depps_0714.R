##################################################
## Project:NC_rebuttal
## Script purpose:Thyroid cell line predicted results-depps
## Date: 2025-06-06
## Author: MenggeLYU
## Version: 1.0
##################################################



# path and library --------------------------------------------------------
setwd("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\Network_fig")
library(data.table)
library(rio)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggvenn)
library(ggrepel)
library(pheatmap)
library(factoextra)
library(pROC)

# load data ---------------------------------------------------------------
# high precision plot
TPC <- rio::import("TPC_prob095_final_20250628.tsv")
FTC238 <- rio::import("FTC238_prob095_final_20250628.tsv")
Nthy <- rio::import("Nthy_prob095_final_20250628.tsv")
plot(1, type = "n", 
     xlim = c(1, max(length(TPC$weighted_precision), 
                     length(Nthy$weighted_precision),
                     length(FTC238$weighted_precision))),
     ylim = c(0, 1),
     xlab = "Protein Pair Rank", 
     ylab = "Weighted Precision",
     main = "Comparison of Weighted Precision Across Cell Lines")
lines(TPC$weighted_precision, col = "blue", lwd = 2)
lines(Nthy$weighted_precision, col = "red", lwd = 2)
lines(FTC238$weighted_precision, col = "green", lwd = 2)
# legend
legend("topright", 
       legend = c("TPC", "Nthy", "FTC238"),
       col = c("blue", "red", "green"),
       lty = 1, lwd = 2,
       cex = 0.8)

FTC238_twocols <- FTC238[,1:2]
FTC238_twocols <- unique(FTC238_twocols)
FTC238_twocols$label <- "FTC238"

TPC_twocols <- TPC[,1:2]
TPC_twocols <- unique(TPC_twocols)
TPC_twocols$label <- "TPC"

Nthy_twocols <- Nthy[,1:2]
Nthy_twocols <- unique(Nthy_twocols)
Nthy_twocols$label <- "Nthy"

TPC <- TPC[TPC$weighted_precision>0.6,]
FTC238<- FTC238[FTC238$weighted_precision>0.6,]
Nthy <- Nthy[Nthy$weighted_precision>0.6,]

TPC$pair <- apply(TPC[, c("protein_A", "protein_B")], 1, function(x) paste(sort(x), collapse = "_"))
Nthy$pair <- apply(Nthy[, c("protein_A", "protein_B")], 1, function(x) paste(sort(x), collapse = "_"))
FTC238$pair <- apply(FTC238[, c("protein_A", "protein_B")], 1, function(x) paste(sort(x), collapse = "_"))

pair_list <- list(
  TPC = unique(TPC$pair),
  Nthy = unique(Nthy$pair),
  FTC238 = unique(FTC238$pair))

ggvenn(pair_list,
       fill_color = c("blue", "red", "green"),
       stroke_size = 0.5,
       set_name_size = 4)

load("E:\\MenggeLYU\\NC_rebuttal_t3\\Comb_threedb\\three_comb_adj_matrix_MGL_20250409.RData")

# test_matrix <- matrix(c(0,1,0,1,0,1,0,1,0), nrow=3, 
#                       dimnames=list(c("A","B","C"), c("A","B","C")))
# print(test_matrix)

result <- which(adj_df== 1, arr.ind = TRUE) %>%
  as.data.frame() %>%
  mutate(
    protein_A = rownames(adj_df)[row],
    protein_B = colnames(adj_df)[col]) %>%
  filter(protein_A != protein_B) %>%  # delete the inner interactions
  select(protein_A, protein_B)

# the number of pairs
overlap_counts <- list(
  "TPC_only" = setdiff(pair_list$TPC, union(pair_list$Nthy, pair_list$FTC238)) %>% length(),
  "Nthy_only" = setdiff(pair_list$Nthy, union(pair_list$TPC, pair_list$FTC238)) %>% length(),
  "FTC238_only" = setdiff(pair_list$FTC238, union(pair_list$TPC, pair_list$Nthy)) %>% length(),   # the unique pairs of each zone
  "TPC_Nthy" = length(intersect(pair_list$TPC, pair_list$Nthy)) - length(Reduce(intersect, pair_list)),
  "TPC_FTC238" = length(intersect(pair_list$TPC, pair_list$FTC238)) - length(Reduce(intersect, pair_list)),
  "Nthy_FTC238" = length(intersect(pair_list$Nthy, pair_list$FTC238)) - length(Reduce(intersect, pair_list)),   # overlapped of two zones
  "All_three" = length(Reduce(intersect, pair_list)))   # all three zones
overlap_df <- data.frame(
  Region = names(overlap_counts),
  Count = unlist(overlap_counts))
# the name of pairs
overlap_pairs <- list(
  "TPC_only" = setdiff(pair_list$TPC, union(pair_list$Nthy, pair_list$FTC238)),
  "Nthy_only" = setdiff(pair_list$Nthy, union(pair_list$TPC, pair_list$FTC238)),
  "FTC238_only" = setdiff(pair_list$FTC238, union(pair_list$TPC, pair_list$Nthy)),
  "TPC_Nthy_only" = setdiff(intersect(pair_list$TPC, pair_list$Nthy), pair_list$FTC238),
  "TPC_FTC238_only" = setdiff(intersect(pair_list$TPC, pair_list$FTC238), pair_list$Nthy),
  "Nthy_FTC238_only" = setdiff(intersect(pair_list$Nthy, pair_list$FTC238), pair_list$TPC),
  "All_three" = Reduce(intersect, pair_list))

ggvenn(pair_list,
       fill_color = c("blue", "red", "green"),
       stroke_size = 0.5,
       set_name_size = 4,
       show_percentage = TRUE,  
       text_size = 5,            
       label_sep = "\n") +       
  labs(title = "Protein Pair Overlap") +
  theme(plot.title = element_text(hjust = 0.5))

result$pair <- apply(result[, c("protein_A", "protein_B")], 1, 
                      function(x) paste(sort(x), collapse = "_"))

venn_regions <- list(
  "TPC_only" = setdiff(pair_list$TPC, union(pair_list$Nthy, pair_list$FTC238)),
  "Nthy_only" = setdiff(pair_list$Nthy, union(pair_list$TPC, pair_list$FTC238)),
  "FTC238_only" = setdiff(pair_list$FTC238, union(pair_list$TPC, pair_list$Nthy)),
  "TPC_Nthy_only" = setdiff(intersect(pair_list$TPC, pair_list$Nthy), pair_list$FTC238),
  "TPC_FTC238_only" = setdiff(intersect(pair_list$TPC, pair_list$FTC238), pair_list$Nthy),
  "Nthy_FTC238_only" = setdiff(intersect(pair_list$Nthy, pair_list$FTC238), pair_list$TPC),
  "All_three" = Reduce(intersect, pair_list))

overlap_results <- sapply(venn_regions, function(x) {
  sum(x %in% result$pair)})

overlap_df2 <- data.frame(
  Region = names(venn_regions),
  Total_in_Venn = lengths(venn_regions),
  In_results = overlap_results,
  Percentage = round(overlap_results / lengths(venn_regions) * 100, 1)) # identification vs three databases

overlap_df2 <- overlap_df2 %>% arrange(desc(In_results))

ggplot(overlap_df2, aes(x = reorder(Region, -In_results), y = In_results)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(In_results, "\n(", Percentage, "%)")), 
            vjust = -0.5, size = 3) +
  labs(x = "Venn Diagram Region", 
       y = "Number of Pairs in Results",
       title = "Overlap between Venn Regions and Results") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# string_g2p <- rio::import("E:\\MenggeLYU\\CHINDA\\STRING\\final\\STRING_protein_list_withlabel_MGL_20240611.tsv")
# protein_df2g <- left_join(protein_df,string_g2p, by ="Prot")
# match <- protein_df2g[which(is.na(protein_df2g$preferred_name)),]
# rio::export(protein_df, "string_protein.tsv")
# prot2 <- as.data.frame(match$Prot)
# rio::export(prot2, "string_protein_nolabel.tsv")
# stringmap <- rio::import("string_mapping (1).tsv")
# second <- stringmap[,c(2,4)]
# colnames(second) <- c("Prot", "Gene")
# first <- protein_df2g[-which(is.na(protein_df2g$preferred_name)),]
# first <- first[,c(1,3)]
# colnames(first) <- c("Prot", "Gene")
# g2p <- rbind(first,second)
# gene <- unique(g2p$Gene)
# rio::export(g2p,"string_gene.tsv")
# edge <- rio::import("ensp_id.xlsx")
# node <- rio::import("STRING network default node.csv")
# node <- node[,c(15,19)]
# colnames(node)[2] <- "A"
# edge1 <- left_join(edge,node, by = "A")
# colnames(edge1)[3]<-'gene_A' 
# colnames(node)[2] <- "B"
# edge2 <- left_join(edge1,node, by = "B")
# colnames(edge2)[4]<-'gene_B' 
# colnames(g2p)[2] <-"gene_A" 
# edge3 <- left_join(edge2,g2p, by = "gene_A")
# colnames(edge3)[5]<-'protein_A' 
# colnames(g2p)[2] <-"gene_B" 
# edge4 <- left_join(edge3,g2p, by = "gene_B")
# colnames(edge4)[6]<-'protein_B' 
# STRING_FINAL <- edge4[,5:6]
# STRING_FINAL <- unique(STRING_FINAL) # 44658
# STRING_FINAL$label <- 'STRING'
# 
# 
# COMB_ALL <- rbind(comb,STRING_FINAL)
# rio::export(COMB_ALL,"network_pair_all.tsv")

# MW ----------------------------------------------------------------------
comb <- rbind(FTC238,TPC,Nthy) # weighted precision > 0.6
comb_pairs <- as.data.frame(unique(comb[,1:2]))
comb_pairs$comb <- apply(comb_pairs[, c("protein_A", "protein_B")], 1, function(x) paste(sort(x), collapse = "_"))
lengthcheck <- unique(comb_pairs$comb) # 25173 pairs
protein <- unique(c(comb_pairs$protein_A,comb_pairs$protein_B))
protein_df <- as.data.frame(protein)
colnames(protein_df) <- "Prot"
protein_iso <-protein_df[grep("-", protein_df$Prot),] 
protein_iso <- as.data.frame(protein_iso)
colnames(protein_iso)[1] <- "protein"
rio::export(protein_iso, "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\MW\\Isoform_paris_all_weightedprecision06.tsv")
protein_cano <- protein_df[-grep("-", protein_df$Prot),] 
protein_cano <- as.data.frame(protein_cano)
rio::export(protein_cano, "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\MW\\Cano_paris_all_weightedprecision06.tsv")
iso_input <- rio::import("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\MW\\uniprot_molecular_weights_iso_0710.csv")
iso_input <- iso_input[,c(1,3)]
iso_input$KDa <- iso_input$molWeight/1000
iso_input$frac <- -0.0312*iso_input$KDa
iso_input$frac <- iso_input$frac+40.033
iso_input$Frac_Num <- ceiling(iso_input$frac)
cano_input <- rio::import("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\MW\\idmapping_model_organism_9606_2025_07_10.tsv")
cano_input <- cano_input[c(2,7)]
cano_input$KDa <- cano_input$Mass/1000
cano_input$frac <- -0.0312*cano_input$KDa
cano_input$frac <- cano_input$frac+40.033
cano_input$Frac_Num <- ceiling(cano_input$frac)

colnames(cano_input)[1:2] <- colnames(iso_input)[1:2] 
prot_mw <- rbind(cano_input,iso_input)

# DEPPS candidate-------------------------------------------------------------------
# comb <- rbind(FTC238,TPC,Nthy) # weighted precision > 0.6
FTC238$cell <- "FTC238"
TPC$cell <- "TPC"
Nthy$cell <- "Nthy"
three_all <- rbind(FTC238,TPC,Nthy)
colnames(three_all)[17] <- "protein_pair"
N_F <- three_all[-which(three_all$cell=="TPC"),]
# ftc238_rows <- three_all %>%
#   filter(cell == "FTC238")
# nthy_rows <- three_all %>%
#   filter(cell == "Nthy")
# common_protein_pairs <- intersect(ftc238_rows$protein_pair, nthy_rows$protein_pair)
# result <- three_all %>%
#   filter(protein_pair %in% common_protein_pairs)
N_F <- as.data.frame(unique(N_F$protein_pair))
colnames(N_F)[1] <- "combine"
rio::export(N_F, "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\NthyvsFTC238_20250710.tsv")
N_T <- three_all[-which(three_all$cell=="FTC238"),]
# TPC_rows <- three_all %>%
#   filter(cell == "TPC")
# nthy_rows <- three_all %>%
#   filter(cell == "Nthy")
# common_protein_pairs <- intersect(TPC_rows$protein_pair, nthy_rows$protein_pair)
# result <- three_all %>%
#   filter(protein_pair %in% common_protein_pairs)
N_T <- as.data.frame(unique(N_T$protein_pair))
colnames(N_T)[1] <- "combine"
rio::export(N_T, "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\NthyvsTPC_20250710.tsv")
# filtered MW -------------------------------------------------------------
# import from python results
N_T <- rio::import("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_T_0710\\NthyandTPC_maxcol_forMWfiltered_20250710.csv")
N_F <- rio::import("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_F_0710\\NthyandFTC238_maxcol_forMWfiltered_20250710.csv")
N_T  <- N_T  %>%
  separate(protein_pair, into = c("protein_1", "protein_2"), sep = "_")
N_F <- N_F  %>%
  separate(protein_pair, into = c("protein_1", "protein_2"), sep = "_")
## TPC-1
colnames(prot_mw)[1] <- "protein_1"
N_T1 <- left_join(N_T,prot_mw, by = "protein_1")
colnames(N_T1)[8] <- "Frac1"
colnames(prot_mw)[1] <- "protein_2"
N_T2 <- left_join(N_T1,prot_mw, by = "protein_2")
colnames(N_T2)[9] <- "Frac2"
N_T3 <- N_T2 %>%
  filter(max_col_protein1 < Frac1)
N_T4 <- N_T3 %>%
  filter(max_col_protein2 < Frac2)
## FTC238
colnames(prot_mw)[1] <- "protein_1"
N_F1 <- left_join(N_F,prot_mw, by = "protein_1")
colnames(N_F1)[8] <- "Frac1"
colnames(prot_mw)[1] <- "protein_2"
N_F2 <- left_join(N_F1,prot_mw, by = "protein_2")
colnames(N_F2)[9] <- "Frac2"
N_F3 <- N_F2 %>%
  filter(max_col_protein1 < Frac1)
N_F4 <- N_F3 %>%
  filter(max_col_protein2 < Frac2)

# combined MW filtered pairs ------------------------------------------------------------------
N_T4$combine <- apply(N_T4[, c("protein_1", "protein_2")], 1, 
                      function(x) paste(sort(x), collapse = "_"))
N_F4$combine <- apply(N_F4[, c("protein_1", "protein_2")], 1, 
                      function(x) paste(sort(x), collapse = "_"))
rio::export(N_T4,"E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_T_pair_all_0710.tsv")
rio::export(N_F4,"E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_F_pair_all_0710.tsv")

N_Tcomb <- as.data.frame(N_T4$combine)
colnames(N_Tcomb)[1] <- "combine"
N_Tcomb <- unique(N_Tcomb)

N_Fcomb <- as.data.frame(N_F4$combine)
colnames(N_Fcomb)[1] <- "combine"
N_Fcomb <- unique(N_Fcomb)


## filtered 1 and NA - delete 0
comb_twocol <- comb[,c(17,9)]
comb_twocol <- comb_twocol[-which(comb_twocol$label==0),]

N_Tcomb_matched <- N_Tcomb %>%
  filter(combine %in% comb_twocol$pair) %>%
  left_join(comb_twocol, by = c("combine" = "pair"))
N_Fcomb_matched <- N_Fcomb %>%
  filter(combine %in% comb_twocol$pair) %>%
  left_join(comb_twocol, by = c("combine" = "pair"))

rio::export(N_Tcomb_matched,"E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_T_pair_listwithlabel_0710.tsv")
rio::export(N_Fcomb_matched,"E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_F_pair_listwithlabel_0710.tsv")

N_Tcomb_new <- N_Tcomb[N_Tcomb$combine %in% comb_twocol$pair, ]
N_Fcomb_new <- N_Fcomb[N_Fcomb$combine %in% comb_twocol$pair, ]

N_Tcomb_new <- as.data.frame(N_Tcomb_new)
colnames(N_Tcomb_new)[1] <- "combine"

N_Fcomb_new <- as.data.frame(N_Fcomb_new)
colnames(N_Fcomb_new)[1] <- "combine"

rio::export(N_Tcomb_new,"E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_T_pair_list_0710.tsv")
rio::export(N_Fcomb_new,"E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_F_pair_list_0710.tsv")

# combined pairs ----------------------------------------------------------
nf <- rio::import("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_F_auc_0710\\N_F_AUC_forDEPPS20250712.csv")
nt <- rio::import("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\N_T_auc_0710\\N_T_AUC_forDEPPS20250712.csv")


## Nthy vs FTC238 - each replicate as a unique value
comb <- nf
comb2 <- unique(comb)

## protein auc average
filter1 <- comb2 %>%
  group_by(protein_pair, data_file, protein_label) %>%
  summarise(mean_sum = mean(auc, na.rm = TRUE), .groups = "drop") 
## protein pair auc sum
filter1 <- unique(filter1) # 65470
filter2 <- filter1 %>%
  group_by(protein_pair, data_file) %>%
  summarise(sum = sum(mean_sum, na.rm = TRUE), .groups = "drop") # 32735

filter3 <- unique(filter2)
## group
filter3$data_file <- gsub("_gaussian", "", filter3$data_file)

filter4 <- filter3 %>%
  mutate(cell = case_when(
    # data_file %in% c("FTC238NO1", "FTC238NO2", "FTC238NO3") ~ "FTC238",
    data_file %in% c("Nthy1", "Nthy2", "Nthy3") ~ "Nthy",
    # data_file %in% c("FTC133NO1", "FTC133NO2", "FTC133NO3") ~ "FTC133",
    data_file %in% c("FTC238NO1", "FTC238NO2", "FTC238NO3") ~ "FTC238",
    TRUE ~ data_file))

# # pick the table to do the likelihood and fold change
# filter4 <- filter3[,c(1,3,4)]

# auc distribution between different cell types
ggplot(filter4, aes(x = sum, fill = cell)) +
  geom_density(alpha = 0.5) + 
  labs(
    title = "Density Plot of Protein Pair Sums by Cell Type",
    x = "Sum of Protein Pair",
    y = "Density"
  ) +
  theme_minimal() + 
  scale_fill_brewer(palette = "Set2") 

# # Nthy vs FTC238
# data_filtered <- filter4 %>%
#   filter(cell %in% c("FTC238", "Nthy"))
# Group by protein_pair and calculate likelihood ratio test and fold change
nf_results <- filter4 %>%
  group_by(protein_pair) %>%
  summarise(
    # Subset data for the two cell types
    Nthy_values = list(sum[cell == "Nthy"]),
    FTC238_values = list(sum[cell == "FTC238"]),
    # Calculate medians
    Nthy_median = median(unlist(Nthy_values), na.rm = TRUE),
    FTC238_median = median(unlist(FTC238_values), na.rm = TRUE),
    # Fold change (log2 scale)
    fold_change = log2(FTC238_median / Nthy_median),
    # Likelihood ratio test
    p_value = ifelse(length(unlist(Nthy_values)) > 1 & length(unlist(FTC238_values)) > 1,
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
  filter(!is.na(adj_p_value) & adj_p_value < 0.05)
results1 <- nf_results[which(nf_results$p_value<0.05),]
results2 <- results1[which(results1$fold_change>log2(1.5)|results1$fold_change < -log2(1.5)),] ## 2628
rio::export(results2, "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\FTC238_Nthy_depp_p005_log21_5_20250714.tsv")
results2pp <- as.data.frame(unique(results2$protein_pair)) # 2628
colnames(results2pp)[1] <- "combine"
rio::export(results2pp, "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\FTC238_Nthy_depp_pair_list_20250714.tsv")


## Nthy vs TPC-1 - each replicate as a unique value
comb <- nt
comb2 <- unique(comb)

## protein auc average
filter1 <- comb2 %>%
  group_by(protein_pair, data_file, protein_label) %>%
  summarise(mean_sum = mean(auc, na.rm = TRUE), .groups = "drop") 
## protein pair auc sum
filter1 <- unique(filter1) # 80202
filter2 <- filter1 %>%
  group_by(protein_pair, data_file) %>%
  summarise(sum = sum(mean_sum, na.rm = TRUE), .groups = "drop") # 40101

filter3 <- unique(filter2)
## group
filter3$data_file <- gsub("_gaussian", "", filter3$data_file)

filter4 <- filter3 %>%
  mutate(cell = case_when(
    # data_file %in% c("FTC238NO1", "FTC238NO2", "FTC238NO3") ~ "FTC238",
    data_file %in% c("Nthy1", "Nthy2", "Nthy3") ~ "Nthy",
    # data_file %in% c("FTC133NO1", "FTC133NO2", "FTC133NO3") ~ "FTC133",
    data_file %in% c("TPC1", "TPC2", "TPC3") ~ "TPC",
    TRUE ~ data_file))

# # pick the table to do the likelihood and fold change
# filter4 <- filter3[,c(1,3,4)]

# auc distribution between different cell types
ggplot(filter4, aes(x = sum, fill = cell)) +
  geom_density(alpha = 0.5) + 
  labs(
    title = "Density Plot of Protein Pair Sums by Cell Type",
    x = "Sum of Protein Pair",
    y = "Density"
  ) +
  theme_minimal() + 
  scale_fill_brewer(palette = "Set2") 

# # Nthy vs FTC133
# data_filtered <- filter4 %>%
#   filter(cell %in% c("FTC238", "Nthy"))
# Group by protein_pair and calculate likelihood ratio test and fold change
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
  filter(!is.na(adj_p_value) & adj_p_value < 0.05)
results1 <- nt_results[which(nt_results$p_value<0.05),]
results2 <- results1[which(results1$fold_change>log2(1.5)|results1$fold_change < -log2(1.5)),] ## 2198
rio::export(results2, "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\TPC_Nthy_depp_p005_log21_5_20250714.tsv")
results2pp <- as.data.frame(unique(results2$protein_pair)) # 2198
colnames(results2pp)[1] <- "combine"
rio::export(results2pp, "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\TPC_Nthy_depp_pair_list_20250714.tsv")


# volcano plot ------------------------------------------------------------
## TPC vs Nthy
results <- nt_results[-which(nt_results$adj_p_value==0),]
results_vol <- results %>%
  mutate(
    significance = case_when(
      fold_change > log2(1.5) & adj_p_value < 0.05 ~ "Upregulated",
      fold_change < -log2(1.5) & adj_p_value < 0.05 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

top_up <- results_vol %>%
  filter(significance == "Upregulated") %>%
  arrange(adj_p_value) %>%
  slice(1:10)

top_down <- results_vol %>%
  filter(significance == "Downregulated") %>%
  arrange(adj_p_value) %>%
  slice(1:10)

top_labeled <- bind_rows(top_up, top_down)

volcano_plot <- ggplot(results_vol, aes(x = fold_change, y = -log10(adj_p_value), color = significance)) +
  geom_point(alpha = 0.7, size = 2) + 
  scale_color_manual(values = c("Upregulated" = "#ea9e53", "Downregulated" = "#5353ea", "Not Significant" = "grey")) +
  labs(
    title = "Volcano Plot of Differentially Expressed Protein Pairs",
    x = "Log2 Fold Change (FTC238 vs TPC)",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Threshold for significance
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "black") +  # Threshold for fold change
  geom_text_repel(
    data = top_labeled,
    aes(label = protein_pair),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf)

ggsave("TPC_Nthy_volcano.pdf", 
       volcano_plot, 
       width = 8,
       height = 6, 
       units = "in",  
       dpi = 300, 
       bg = "white" ) 

# ggplot(results_vol, aes(x = fold_change, y = -log10(adj_p_value), color = significance)) +
#   geom_point(alpha = 0.7, size = 2) + 
#   scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
#   labs(
#     title = "Volcano Plot of Differentially Expressed Protein Pairs",
#     x = "Log2 Fold Change (FTC133 vs Nthy)",
#     y = "-Log10 Adjusted P-Value"
#   ) +
#   theme_minimal() +
#   theme(
#     legend.title = element_blank(),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Threshold for significance
#   geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "black")  # Threshold for fold change


## FTC238 vs Nthy
results <- nf_results[-which(nf_results$adj_p_value==0),]
results_vol <- results %>%
  mutate(
    significance = case_when(
      fold_change > log2(1.5) & adj_p_value < 0.05 ~ "Upregulated",
      fold_change < -log2(1.5) & adj_p_value < 0.05 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

top_up <- results_vol %>%
  filter(significance == "Upregulated") %>%
  arrange(adj_p_value) %>%
  slice(1:10)

top_down <- results_vol %>%
  filter(significance == "Downregulated") %>%
  arrange(adj_p_value) %>%
  slice(1:10)

top_labeled <- bind_rows(top_up, top_down)

volcano_plot <- ggplot(results_vol, aes(x = fold_change, y = -log10(adj_p_value), color = significance)) +
  geom_point(alpha = 0.7, size = 2) + 
  scale_color_manual(values = c("Upregulated" = "#ea9e53", "Downregulated" = "#5353ea", "Not Significant" = "grey")) +
  labs(
    title = "Volcano Plot of Differentially Expressed Protein Pairs",
    x = "Log2 Fold Change (FTC238 vs TPC)",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Threshold for significance
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "black") +  # Threshold for fold change
  geom_text_repel(
    data = top_labeled,
    aes(label = protein_pair),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf)

ggsave("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\FTC238_Nthy_volcano.pdf", 
       volcano_plot, 
       width = 8,
       height = 6, 
       units = "in",  
       dpi = 300,
       bg = "white")

# all three pairs ---------------------------------------------------------

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

unique_pairs <- allthree_withlabel2 %>%
  distinct(combine, cell) %>% 
  split(.$cell) %>%            
  map(~ unique(.x$combine))    

map_int(unique_pairs, length)

nthy_only <- setdiff(unique_pairs$Nthy, union(unique_pairs$FTC238, unique_pairs$TPC))
ftc238_only <- setdiff(unique_pairs$FTC238, union(unique_pairs$Nthy, unique_pairs$TPC))
tpc_only <- setdiff(unique_pairs$TPC, union(unique_pairs$Nthy, unique_pairs$FTC238))

venn_data <- list(
  Nthy = unique_pairs$Nthy,
  FTC238 = unique_pairs$FTC238,
  TPC = unique_pairs$TPC)

ggvenn(
  venn_data, 
  fill_color = c("#1b9e77", "#d95f02", "#7570b3"),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 3) +
  labs(title = "Overlap of Protein Pairs Across Cell Types")

rep_counts <- allthree_withlabel2 %>%
  group_by(combine, cell) %>%
  summarise(
    n_replicates = n_distinct(data_file),  
    .groups = 'drop')

ggplot(rep_counts, aes(x = cell, fill = factor(n_replicates))) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set2", name = "Number of Replicates") +
  labs(
    x = "Cell Type",
    y = "Count of Protein Pairs",
    title = "Protein Pair Replication Across Biological Replicates") +
  theme_minimal()

consistent_pairs <- rep_counts %>%
  filter(n_replicates == 3) %>%
  distinct(combine, cell)

partial_pairs <- rep_counts %>%
  filter(n_replicates == 2) %>%
  distinct(combine, cell)

rep_counts <- allthree_withlabel2 %>%
  mutate(
    comparison = case_when(
      cell %in% c("Nthy", "FTC238") ~ "Nthy_vs_FTC238",
      cell %in% c("Nthy", "TPC") ~ "Nthy_vs_TPC",
      TRUE ~ "Other")) %>%
  filter(comparison != "Other") %>%
  group_by(combine, comparison) %>%
  summarise(
    n_replicates = n_distinct(data_file),
    .groups = 'drop')
high_rep_pairs <- rep_counts %>%
  filter(n_replicates >= 5) %>%
  arrange(desc(n_replicates))

head(high_rep_pairs)

ggplot(rep_counts, aes(x = n_replicates, fill = comparison)) +
  geom_bar(position = "dodge", alpha = 0.8) +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("#1f78b4", "#33a02c")) +
  labs(
    x = "Number of Biological Replicates",
    y = "Count of Protein Pairs",
    title = "Protein Pair Replication Frequency",
    subtitle = "Dashed line indicates 5+ replicates cutoff") +
  theme_minimal() +
  theme(legend.position = "top")

# check the fraction - max cols -------------------------------------------
## left_ join the max cols
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

unique_pairs <- allthree_withlabel4 %>%
  distinct(protein_pair, cell) %>% 
  split(.$cell) %>%            
  map(~ unique(.x$protein_pair))    

map_int(unique_pairs, length)

nthy_only <- setdiff(unique_pairs$Nthy, union(unique_pairs$FTC238, unique_pairs$TPC))
ftc238_only <- setdiff(unique_pairs$FTC238, union(unique_pairs$Nthy, unique_pairs$TPC))
tpc_only <- setdiff(unique_pairs$TPC, union(unique_pairs$Nthy, unique_pairs$FTC238))

venn_data <- list(
  Nthy = unique_pairs$Nthy,
  FTC238 = unique_pairs$FTC238,
  TPC = unique_pairs$TPC)

ggvenn(
  venn_data, 
  fill_color = c("#1b9e77", "#d95f02", "#7570b3"),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 3) +
  labs(title = "Overlap of Protein Pairs Across Cell Types")

# depps filter with 8000 pairs --------------------------------------------
depps_nt <-  nt_results[-which(nt_results$adj_p_value==0),]
depps_nt <- depps_nt[which(depps_nt$adj_p_value<0.05),]
depps_nt <- depps_nt[which(depps_nt$fold_change > log2(1.5) |depps_nt$fold_change< -log2(1.5) ),] #2183
depps_nt2 <- depps_nt[which(depps_nt$protein_pair %in% allthree_withlabel4$protein_pair==TRUE),] # 2014

depps_nf <-  nf_results[-which(nf_results$adj_p_value==0),]
depps_nf <- depps_nf[which(depps_nf$adj_p_value<0.05),]
depps_nf <- depps_nf[which(depps_nf$fold_change > log2(1.5) |depps_nf$fold_change< -log2(1.5) ),] #2607
depps_nf2 <- depps_nf[which(depps_nf$protein_pair %in% allthree_withlabel4$protein_pair==TRUE),] # 1967


# nt depps network --------------------------------------------------------
depps_nt2  <- depps_nt2  %>%
  separate(protein_pair, into = c("protein_1", "protein_2"), sep = "_")
protein_nt2 <- as.data.frame(unique(c(depps_nt2$protein_1,depps_nt2$protein_2)))
colnames(protein_nt2) <- 'prot'
rio::export(protein_nt2, "E://MenggeLYU//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//tpc_nthy_network_final_0714.tsv")
nt2_gene <- rio::import("E://MenggeLYU//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//tpc_nthy_network_string_mapping_gene_0714.tsv")
protein_nt2$queryItem <- gsub("-(.*)",'', protein_nt2$prot)
protein_nt2 <- left_join(protein_nt2,nt2_gene, by = "queryItem")
protein_nt2 <- protein_nt2[,c(1,2,4,5)]
## string
nt2_egde <- rio::import('E://MenggeLYU//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//tpc_nthy_network_string_mapping_pairs_0714.tsv')
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

## nt CF-DIA-MS cluster
Nthy1 <- rio::import("E:\\MenggeLYU\\Figure_NC\\Nthy_202410\\Nthy1_gaussian_withoutlog2.tsv")
Nthy2 <- rio::import("E:\\MenggeLYU\\Figure_NC\\Nthy_202410\\Nthy2_gaussian_withoutlog2.tsv")
Nthy3 <- rio::import("E:\\MenggeLYU\\Figure_NC\\Nthy_202410\\Nthy3_gaussian_withoutlog2.tsv")

TPC1 <- rio::import("E:\\MenggeLYU\\Figure_NC\\TPC_202410\\TPC1_gaussian_withoutlog2.tsv")
TPC2 <- rio::import("E:\\MenggeLYU\\Figure_NC\\TPC_202410\\TPC2_gaussian_withoutlog2.tsv")
TPC3 <- rio::import("E:\\MenggeLYU\\Figure_NC\\TPC_202410\\TPC3_gaussian_withoutlog2.tsv")

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

## cluster1-172 prot nthy-------------------------------------------------------
current_row <- overlap[1,]
proteins_to_plot <- unlist(current_row$proteins)
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot172 <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("cluster_172prot_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot172)[1] <- "Prot"
rio::export(cluster_prot172,"cluster_172prot.tsv")
cluster_prot172$P2G <- gsub("-(.*)","",cluster_prot172$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot172 <- left_join(cluster_prot172,gene, by = "P2G")
cluster_prot172$comb <- paste(cluster_prot172$preferredName, " (", cluster_prot172$Prot, ")", sep = "")
colnames(cluster_prot172)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot172, by = "V1")

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
    show_rownames = TRUE) # 10 cols
row_order <- p1$tree_row$order
ordered_rows <- rownames(subset_data2)[row_order]
## A4 heatmap standard
pdf("NvsT_N172prots_ori.pdf", 
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
  gaps_col = 10,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("NvsT_N172prots_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot172$preferredName & node2 %in% cluster_prot172$preferredName)
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
rio::export(cfms_marked,"cluster172prot_nthy.tsv")

## cluster1-172 prot tpc-------------------------------------------------------
current_row <- overlap[1,]
proteins_to_plot <- unlist(current_row$proteins)
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot172 <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("cluster_172prot_gene.tsv") # tpc-1 same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot172)[1] <- "Prot"
rio::export(cluster_prot172,"cluster_172prot.tsv")
cluster_prot172$P2G <- gsub("-(.*)","",cluster_prot172$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot172 <- left_join(cluster_prot172,gene, by = "P2G")
cluster_prot172$comb <- paste(cluster_prot172$preferredName, " (", cluster_prot172$Prot, ")", sep = "")
colnames(cluster_prot172)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot172, by = "V1")

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
pdf("NvsT_T172prots_ori.pdf", 
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
  gaps_col = 10,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - tpc-1 is similar with Nthy
string <- rio::import("NvsT_N172prots_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot172$preferredName & node2 %in% cluster_prot172$preferredName)
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
rio::export(cfms_marked,"cluster172prot_nthy.tsv")
node1 <- cfms_marked[,c(1,3)]
node2 <- cfms_marked[,c(2,3)]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node_filtered <- cf_node %>%
  group_by(node1) %>%                          
  arrange(desc(abs(fold_change)), .by_group = TRUE) %>%  
  slice_head(n = 1) %>%                        
  ungroup() 
colnames(cf_node_filtered)[1] <- "name"
rio::export(cf_node_filtered,"cluster172prot_foldchange.tsv")


## cluster2-87 prot nthy-------------------------------------------------------
current_row <- overlap[2,]
proteins_to_plot <- unlist(current_row$proteins)
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot87 <- as.data.frame(subset_data$V1)
rio::export(cluster_prot87,"cluster_87prot.tsv")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("cluster_87prot_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot87)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot87$P2G <- gsub("-(.*)","",cluster_prot87$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot87 <- left_join(cluster_prot87,gene, by = "P2G")
cluster_prot87$comb <- paste(cluster_prot87$preferredName, " (", cluster_prot87$Prot, ")", sep = "")
colnames(cluster_prot87)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot87, by = "V1")

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
pdf("NvsT_N87prots_ori.pdf", 
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
  gaps_col = 25,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("NvsT_N87prots_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot87$preferredName & node2 %in% cluster_prot87$preferredName)
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
rio::export(cfms_marked,"cluster87prot_nthy.tsv")

## cluster2-87 prot tpc-------------------------------------------------------
current_row <- overlap[2,]
proteins_to_plot <- unlist(current_row$proteins)
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot87 <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("cluster_87prot_gene.tsv") # tpc-1 same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot87)[1] <- "Prot"
rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot87$P2G <- gsub("-(.*)","",cluster_prot87$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot87 <- left_join(cluster_prot87,gene, by = "P2G")
cluster_prot87$comb <- paste(cluster_prot87$preferredName, " (", cluster_prot87$Prot, ")", sep = "")
colnames(cluster_prot87)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot87, by = "V1")

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
pdf("NvsT_T87prots_ori.pdf", 
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
  gaps_col = 26,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - tpc-1 is similar with Nthy
string <- rio::import("NvsT_N87prots_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot87$preferredName & node2 %in% cluster_prot87$preferredName)
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
rio::export(cfms_marked,"cluster87prot_nthy.tsv")
node1 <- cfms_marked[,c(1,3)]
node2 <- cfms_marked[,c(2,3)]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node_filtered <- cf_node %>%
  group_by(node1) %>%                          
  arrange(desc(abs(fold_change)), .by_group = TRUE) %>%  
  slice_head(n = 1) %>%                        
  ungroup() 
colnames(cf_node_filtered)[1] <- "name"
rio::export(cf_node_filtered,"cluster87prot_foldchange.tsv")

## cluster3-84 prot nthy-------------------------------------------------------
current_row <- overlap[3,]
proteins_to_plot <- unlist(current_row$proteins)
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot84 <- as.data.frame(subset_data$V1)
rio::export(cluster_prot84,"cluster_84prot.tsv")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("cluster_84prot_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot84)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot84$P2G <- gsub("-(.*)","",cluster_prot84$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot84 <- left_join(cluster_prot84,gene, by = "P2G")
cluster_prot84$comb <- paste(cluster_prot84$preferredName, " (", cluster_prot84$Prot, ")", sep = "")
colnames(cluster_prot84)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot84, by = "V1")

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
pdf("NvsT_N84prots_ori.pdf", 
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
  gaps_col = 25,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("NvsT_N84prots_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot84$preferredName & node2 %in% cluster_prot84$preferredName)
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
rio::export(cfms_marked,"cluster84prot_nthy.tsv")

## cluster3-84 prot tpc-------------------------------------------------------
current_row <- overlap[3,]
proteins_to_plot <- unlist(current_row$proteins)
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot84 <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("cluster_84prot_gene.tsv") # tpc-1 same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot84)[1] <- "Prot"
rio::export(cluster_prot84,"cluster_84prot.tsv")
cluster_prot84$P2G <- gsub("-(.*)","",cluster_prot84$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot84 <- left_join(cluster_prot84,gene, by = "P2G")
cluster_prot84$comb <- paste(cluster_prot84$preferredName, " (", cluster_prot84$Prot, ")", sep = "")
colnames(cluster_prot84)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot84, by = "V1")

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
pdf("NvsT_T84prots_ori.pdf", 
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
  gaps_col = 26,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - tpc-1 is similar with Nthy
string <- rio::import("NvsT_N84prots_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot87$preferredName & node2 %in% cluster_prot87$preferredName)
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
rio::export(cfms_marked,"cluster84prot_nthy.tsv")
node1 <- cfms_marked[,c(1,3)]
node2 <- cfms_marked[,c(2,3)]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node_filtered <- cf_node %>%
  group_by(node1) %>%                          
  arrange(desc(abs(fold_change)), .by_group = TRUE) %>%  
  slice_head(n = 1) %>%                        
  ungroup() 
colnames(cf_node_filtered)[1] <- "name"
rio::export(cf_node_filtered,"cluster84prot_foldchange.tsv")
## cluster4-80 prot nthy-------------------------------------------------------
current_row <- overlap[4,]
proteins_to_plot <- unlist(current_row$proteins)
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot80 <- as.data.frame(subset_data$V1)
rio::export(cluster_prot80,"cluster_80prot.tsv")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("cluster_80prot_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot80)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot80$P2G <- gsub("-(.*)","",cluster_prot80$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot80 <- left_join(cluster_prot80,gene, by = "P2G")
cluster_prot80$comb <- paste(cluster_prot80$preferredName, " (", cluster_prot80$Prot, ")", sep = "")
colnames(cluster_prot80)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot80, by = "V1")

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
pdf("NvsT_N80prots_ori.pdf", 
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
string <- rio::import("NvsT_N80prots_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot80$preferredName & node2 %in% cluster_prot80$preferredName)
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
rio::export(cfms_marked,"cluster80prot_nthy.tsv")

## cluster4-80 prot tpc-------------------------------------------------------
current_row <- overlap[4,]
proteins_to_plot <- unlist(current_row$proteins)
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot80 <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("cluster_80prot_gene.tsv") # tpc-1 same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot80)[1] <- "Prot"
rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot80$P2G <- gsub("-(.*)","",cluster_prot80$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot80 <- left_join(cluster_prot80,gene, by = "P2G")
cluster_prot80$comb <- paste(cluster_prot80$preferredName, " (", cluster_prot80$Prot, ")", sep = "")
colnames(cluster_prot80)[1] <- "V1"
subset_data1 <- left_join(subset_data,cluster_prot80, by = "V1")

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
pdf("NvsT_T80prots_ori.pdf", 
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

## string mapping results - tpc-1 is similar with Nthy
string <- rio::import("NvsT_N80prots_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nt5 %>%
  filter(node1 %in% cluster_prot80$preferredName & node2 %in% cluster_prot80$preferredName)
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
rio::export(cfms_marked,"cluster80prot_nthy.tsv")
node1 <- cfms_marked[,c(1,3)]
node2 <- cfms_marked[,c(2,3)]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node_filtered <- cf_node %>%
  group_by(node1) %>%                          
  arrange(desc(abs(fold_change)), .by_group = TRUE) %>%  
  slice_head(n = 1) %>%                        
  ungroup() 
colnames(cf_node_filtered)[1] <- "name"
rio::export(cf_node_filtered,"cluster80prot_foldchange.tsv")

## cluster 5 string ---------------------------------------------------------
current_row <- overlap[5:58,]
proteins_to_plot <- unlist(current_row$proteins)
cluster_prot <- as.data.frame(unique(proteins_to_plot))
colnames(cluster_prot)[1] <- "Prot"
rio::export(cluster_prot,"cluster5toall_prot.tsv")

gene <- rio::import("cluster5toall_prot_gene.tsv")
gene <- gene[,c(2,4)]
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")

string <- rio::import("cluster5toall_prot_string.tsv")
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
cfms_marked2 <- cfms_marked[,c(1,2,4)]
cfms_marked2 <- unique(cfms_marked2)
rio::export(cfms_marked2,"cluster_5toall_prot_nthy.tsv")
node1 <- cfms_marked[,c(1,3)]
node2 <- cfms_marked[,c(2,3)]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node_filtered <- cf_node %>%
  group_by(node1) %>%                          
  arrange(desc(abs(fold_change)), .by_group = TRUE) %>%  
  slice_head(n = 1) %>%                        
  ungroup() 
colnames(cf_node_filtered)[1] <- "name"
rio::export(cf_node_filtered,"cluster_5toall_foldchange.tsv")

# nf depps network --------------------------------------------------------
depps_nf2  <- depps_nf2  %>%
  separate(protein_pair, into = c("protein_1", "protein_2"), sep = "_")
protein_nf2 <- as.data.frame(unique(c(depps_nf2$protein_1,depps_nf2$protein_2)))
colnames(protein_nf2) <- 'prot'
rio::export(protein_nf2, "E://MenggeLYU//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//ftc238_nthy_network_final_0714.tsv")
nf2_gene <- rio::import("E://MenggeLYU//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//ftc238_nthy_network_string_mapping_gene_0714.tsv")
protein_nf2$queryItem <- gsub("-(.*)",'', protein_nf2$prot)
protein_nf2 <- left_join(protein_nf2,nf2_gene, by = "queryItem")
protein_nf2 <- protein_nf2[,c(1,2,4,5)]
## string
nf2_egde <- rio::import('E://MenggeLYU//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//ftc238_nthy_network_string_mapping_pairs_0714.tsv')
nf2_edge_string <- unique(nf2_egde[,c(1:2)])
colnames(nf2_edge_string) <- c("node1","node2")
nf2_edge_string$edge <- "STRING"
nf2_edge_string$fold_change <- NA
## CF-MS
colnames(protein_nf2)[1] <- "protein_1"
depps_nf3 <- left_join(depps_nf2, protein_nf2, by = "protein_1")
colnames(depps_nf3)[12] <- "node1"
depps_nf3 <- depps_nf3[,c(1:9,12)]
colnames(protein_nf2)[1] <- "protein_2"
depps_nf4 <- left_join(depps_nf3, protein_nf2, by = "protein_2")
colnames(depps_nf4)[13] <- "node2"

depps_nf5 <- depps_nf4[,c(1:10,13)]
depps_nf5 <- depps_nf5[,c(10,11,7)]
depps_nf5$edge <- "CF-DIA-MS"

## nf CF-DIA-MS cluster
Nthy1 <- rio::import("E:\\MenggeLYU\\Figure_NC\\Nthy_202410\\Nthy1_gaussian_withoutlog2.tsv")
Nthy2 <- rio::import("E:\\MenggeLYU\\Figure_NC\\Nthy_202410\\Nthy2_gaussian_withoutlog2.tsv")
Nthy3 <- rio::import("E:\\MenggeLYU\\Figure_NC\\Nthy_202410\\Nthy3_gaussian_withoutlog2.tsv")

FTC2381 <- rio::import("E:\\MenggeLYU\\Figure_NC\\FTC238_202410\\FTC238NO1_gaussian_withoutlog2.tsv")
FTC2382 <- rio::import("E:\\MenggeLYU\\Figure_NC\\FTC238_202410\\FTC238NO2_gaussian_withoutlog2.tsv")
FTC2383 <- rio::import("E:\\MenggeLYU\\Figure_NC\\FTC238_202410\\FTC238NO3_gaussian_withoutlog2.tsv")

combined <- bind_rows(Nthy1, Nthy2, Nthy3)
prot <- unique(c(depps_nf4$protein_1,depps_nf4$protein_2))

Nthy_mean <- combined %>%
  group_by(V1) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  

combined2 <- bind_rows(FTC2381, FTC2382, FTC2383)
#prot <- unique(c(depps_nt4$protein_1,depps_nt4$protein_2))

FTC238_mean <- combined2 %>%
  group_by(V1) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  

heatmap <- Nthy_mean[Nthy_mean$V1%in%prot,] # 959
heatmap_clean <- heatmap[,2:61]
rownames(heatmap_clean) <- heatmap$V1
map1 <- pheatmap::pheatmap(heatmap_clean, scale = "row",
                           clustering_distance_rows = "euclidean", 
                           clustering_method = "ward.D2", 
                           cluster_rows = T,cluster_cols = F,
                           treeheight_row = 0,
                           show_rownames = FALSE)

heatmap2 <- FTC238_mean[FTC238_mean$V1%in%prot,] # 959
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
     main = "FTC238 Cluster Dendrogram (Compared with Nthy-ori 3-1)")
abline(h = 50, col = "red")

# fviz_nbclust(heatmap_clean, FUN = hcut, method = "wss")
# fviz_nbclust(heatmap_clean, FUN = hcut, method = "silhouette")

## cluster number
k <- 9
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
colnames(overlap_unnested)[2] <- "FTC238 clusters"
colnames(overlap_unnested)[4] <- "protein list"


## ftc238 cluster1 prot nthy-------------------------------------------------------
current_row <- overlap_unnested[1,]
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1)
rio::export(cluster_prot,"Nthy_cluster1vsFTC238_cluster1_52prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster1vsFTC238_cluster1_52prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238vsnthy_cluster1_0830\\Nthy_cluster1vsFTC238_cluster1_52prots_nthy_gene.tsv")
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
pdf("Nthy_cluster1vsFTC238_cluster1_52prots_nthy_ori.pdf", 
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
string <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238vsnthy_cluster1_0830\\Nthy_cluster1vsFTC238_cluster1_52prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
rio::export(cfms_marked,"Nthy_cluster1vsFTC238_cluster1_52prots_nthy_final.tsv")

## ftc238 cluster1 prot ftc238-------------------------------------------------------
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
rio::export(cluster_prot,"Nthy_cluster1vsFTC238_cluster1_52prots_ftc238.tsv")
rio::export(cluster_prot,"Nthy_cluster1vsFTC238_cluster1_52prots_ftc238.xlsx")



## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238vsnthy_cluster1_0830\\Nthy_cluster1vsFTC238_cluster1_52prots_nthy_gene.tsv") # ftc238 same with nthy
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
pdf("Nthy_cluster1vsFTC238_cluster1_52prots_ftc238_ori.pdf", 
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
string <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238vsnthy_cluster1_0830\\Nthy_cluster1vsFTC238_cluster1_52prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
rio::export(cfms_marked,"Nthy_cluster1vsFTC238_cluster1_52prots_ftc238_final.tsv")

### cluster1 ftc238 vs nthy nodes mean intensity ----------------------------------------------
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

FTC2381_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\FTC238NO1_gaussian_withoutlog2_normalized.tsv")
FTC2382_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\FTC238NO2_gaussian_withoutlog2_normalized.tsv")
FTC2383_nor <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\normalized_matrices\\FTC238NO3_gaussian_withoutlog2_normalized.tsv")

combined_nor <- bind_rows(Nthy1_nor, Nthy2_nor, Nthy3_nor)
Nthy_nor_mean <- combined_nor %>%
  group_by(V1) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  

combined2_nor <- bind_rows(FTC2381_nor, FTC2382_nor, FTC2383_nor)
#prot <- unique(c(depps_nt4$protein_1,depps_nt4$protein_2))

FTC238_nor_mean <- combined2_nor %>%
  group_by(V1) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  


auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list)
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") #45
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238vsnthy_cluster1_0830\\Nthy_cluster1vsFTC238_cluster1_52prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(FTC238_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","FTC238")
cf_node_ftc238 <- left_join(cf_node, auc_results2, by = "V1") #45
#cf_node_ftc238$FTC238 <- log10(cf_node_ftc238$FTC238)
colnames(cf_node_ftc238)[1] <- "name"
rio::export(cf_node_ftc238, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238vsnthy_cluster1_0830\\Nthy_cluster1vsFTC238_cluster1_52prots_ftc238_intensity.tsv")

# cf_node_filtered <- cf_node %>%
#   group_by(node1) %>%                          
#   arrange(desc(abs(fold_change)), .by_group = TRUE) %>%  
#   slice_head(n = 1) %>%                        
#   ungroup() 
# colnames(cf_node_filtered)[1] <- "name"
# rio::export(cf_node_filtered,"cluster3_72protinNF_foldchange.tsv")


## ftc238 cluster3 prot nthy-------------------------------------------------------
setwd("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238vsnthy_cluster3_0830")
current_row <- overlap_unnested[15,] # line 15 - cluster3
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1) #210
rio::export(cluster_prot,"Nthy_cluster3vsFTC238_cluster3_210prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster3vsFTC238_cluster3_210prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("Nthy_cluster3vsFTC238_cluster3_210prots_nthy_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
# cluster_prot[duplicated(cluster_prot$Prot),]
cluster_prot[duplicated(cluster_prot$P2G),]
# > cluster_prot[duplicated(cluster_prot$P2G),]
# Prot    P2G preferredName
# 105 Q14980-2 Q14980         NUMA1
# 106 Q14980-3 Q14980         NUMA1
# 107 Q14980-3 Q14980         NUMA1
# 163   Q99996 Q99996         AKAP9
# 164 Q99996-3 Q99996         AKAP9
# 165 Q99996-3 Q99996         AKAP9
cluster_prot <- unique(cluster_prot) # 214 to 210
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
pdf("Nthy_cluster3vsFTC238_cluster3_210prots_nthy_ori2.pdf", 
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
  cellheight = 3.7,  
  # fontfamily = "Arial",       
  fontsize = 3.7,              
  fontsize_row = 4,           
  fontsize_col = 4)

dev.off()

## string mapping results
string <- rio::import("Nthy_cluster3vsFTC238_cluster3_210prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
cfms_marked <- unique(cfms_marked)
rio::export(cfms_marked,"Nthy_cluster3vsFTC238_cluster3_210prots_nthy_final.tsv")

## ftc238 cluster3 prot ftc238-------------------------------------------------------
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
rio::export(cluster_prot,"Nthy_cluster3vsFTC238_cluster3_210prots_ftc238.tsv")
rio::export(cluster_prot,"Nthy_cluster3vsFTC238_cluster3_210prots_ftc238.xlsx")


## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("Nthy_cluster3vsFTC238_cluster3_210prots_nthy_gene.tsv") # ftc238 same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot[duplicated(cluster_prot$P2G),]
# Prot    P2G preferredName
# 105 Q14980-2 Q14980         NUMA1
# 106 Q14980-3 Q14980         NUMA1
# 107 Q14980-3 Q14980         NUMA1
# 163   Q99996 Q99996         AKAP9
# 164 Q99996-3 Q99996         AKAP9
# 165 Q99996-3 Q99996         AKAP9
cluster_prot <- unique(cluster_prot)
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
pdf("Nthy_cluster3vsFTC238_cluster3_210prots_ftc238_ori2.pdf", 
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
  cellheight = 3.7,  
  # fontfamily = "Arial",       
  fontsize = 3.7,              
  fontsize_row = 4,           
  fontsize_col = 4)

dev.off()

## string mapping results - ftc238 is similar with Nthy
string <- rio::import("Nthy_cluster3vsFTC238_cluster3_210prots_ftc238_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
cfms_marked <- unique(cfms_marked)
rio::export(cfms_marked,"Nthy_cluster3vsFTC238_cluster3_210prots_ftc238_final.tsv")

### cluster3 ftc238 vs nthy nodes mean intensity ----------------------------------------------
node1 <- cfms_marked[,1]
node2 <- cfms_marked[,2]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node <- unique(cf_node) # 198
colnames(cf_node)[1] <- "preferredName" 
cf_node <- left_join(cf_node, cluster_prot, by = "preferredName") # 200
cf_node <- unique(cf_node)
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

protein_list <- cf_node$V1 # 200
auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list)
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") # 200
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster3vsFTC238_cluster3_210prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(FTC238_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","FTC238")
cf_node_ftc238 <- left_join(cf_node, auc_results2, by = "V1") # 200
#cf_node_ftc238$FTC238 <- log10(cf_node_ftc238$FTC238)
colnames(cf_node_ftc238)[1] <- "name"
cf_node_ftc238 <- unique(cf_node_ftc238)
rio::export(cf_node_ftc238, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster3vsFTC238_cluster3_210prots_ftc238_intensity.tsv")


## nthy cluster5 ftc238 cluster 8 prot nthy-------------------------------------------------------
current_row <- overlap_unnested[32,] # line 32 - cluster5
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1) # 26
rio::export(cluster_prot,"Nthy_cluster5vsFTC238_cluster8_26prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster5vsFTC238_cluster8_26prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238_cluster8vsnthy_cluster5_0830\\Nthy_cluster5vsFTC238_cluster8_26prots_nthy_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
# cluster_prot[duplicated(cluster_prot$Prot),]
cluster_prot[duplicated(cluster_prot$P2G),]
# > cluster_prot[duplicated(cluster_prot$P2G),]
# Prot    P2G preferredName
# 105 Q14980-2 Q14980         NUMA1
# 106 Q14980-3 Q14980         NUMA1
# 107 Q14980-3 Q14980         NUMA1
# 163   Q99996 Q99996         AKAP9
# 164 Q99996-3 Q99996         AKAP9
# 165 Q99996-3 Q99996         AKAP9
cluster_prot <- unique(cluster_prot) # 26
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
pdf("Nthy_cluster5vsFTC238_cluster8_26prots_nthy_ori.pdf", 
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
string <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238_cluster8vsnthy_cluster5_0830\\Nthy_cluster5vsFTC238_cluster8_26prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
cfms_marked <- unique(cfms_marked)
rio::export(cfms_marked,"Nthy_cluster5vsFTC238_cluster8_26prots_nthy_final.tsv")

## nthy cluster5 ftc238 cluster 8 prot ftc238-------------------------------------------------------
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
rio::export(cluster_prot,"Nthy_cluster5vsFTC238_cluster8_26prots_ftc238.tsv")
rio::export(cluster_prot,"Nthy_cluster5vsFTC238_cluster8_26prots_ftc238.xlsx")


## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238_cluster8vsnthy_cluster5_0830\\Nthy_cluster5vsFTC238_cluster8_26prots_nthy_gene.tsv") # ftc238 same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot[duplicated(cluster_prot$P2G),]
# Prot    P2G preferredName
# 105 Q14980-2 Q14980         NUMA1
# 106 Q14980-3 Q14980         NUMA1
# 107 Q14980-3 Q14980         NUMA1
# 163   Q99996 Q99996         AKAP9
# 164 Q99996-3 Q99996         AKAP9
# 165 Q99996-3 Q99996         AKAP9
cluster_prot <- unique(cluster_prot)
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
pdf("Nthy_cluster5vsFTC238_cluster8_26prots_ftc238_ori.pdf", 
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
string <- rio::import("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238_cluster8vsnthy_cluster5_0830\\Nthy_cluster5vsFTC238_cluster8_26prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
cfms_marked <- unique(cfms_marked)
rio::export(cfms_marked,"Nthy_cluster5vsFTC238_cluster8_26prots_ftc238_final.tsv")

### nthy cluster5 ftc238 cluster 8 nodes mean intensity ----------------------------------------------
node1 <- cfms_marked[,1]
node2 <- cfms_marked[,2]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node <- unique(cf_node) # 25
colnames(cf_node)[1] <- "preferredName" 
cf_node <- left_join(cf_node, cluster_prot, by = "preferredName") # 45
cf_node <- unique(cf_node)
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

protein_list <- cf_node$V1 # 25
auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list) # 25
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") # 25
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster5vsFTC238_cluster8_26prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(FTC238_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","FTC238")
cf_node_ftc238 <- left_join(cf_node, auc_results2, by = "V1") # 25
#cf_node_ftc238$FTC238 <- log10(cf_node_ftc238$FTC238)
colnames(cf_node_ftc238)[1] <- "name"
cf_node_ftc238 <- unique(cf_node_ftc238)
rio::export(cf_node_ftc238, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster5vsFTC238_cluster8_26prots_ftc238_intensity.tsv")

## nthy ftc238 cluster 8 prot nthy-------------------------------------------------------
setwd("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238vsnthy_cluster8_0830")
current_row <- overlap_unnested[54,] # line 54 - cluster8
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1) # 58
rio::export(cluster_prot,"Nthy_cluster8vsFTC238_cluster8_58prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster8vsFTC238_cluster8_58prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("Nthy_cluster8vsFTC238_cluster8_58prots_nthy_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
# cluster_prot[duplicated(cluster_prot$Prot),]
cluster_prot[duplicated(cluster_prot$P2G),]
# > cluster_prot[duplicated(cluster_prot$P2G),]
# Prot    P2G preferredName
# 35   P50990 P50990          CCT8
# 36 P50990-2 P50990          CCT8
# 37 P50990-2 P50990          CCT8
cluster_prot <- unique(cluster_prot) # 58
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
pdf("Nthy_cluster8vsFTC238_cluster8_58prots_nthy_ori.pdf", 
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
  gaps_col = 28,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("Nthy_cluster8vsFTC238_cluster8_58prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
cfms_marked <- unique(cfms_marked) #140
rio::export(cfms_marked,"Nthy_cluster8vsFTC238_cluster8_58prots_nthy_final.tsv")

## nthy cluster5 ftc238 cluster 8 prot ftc238-------------------------------------------------------
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
rio::export(cluster_prot,"Nthy_cluster8vsFTC238_cluster8_58prots_ftc238.tsv")
rio::export(cluster_prot,"Nthy_cluster8vsFTC238_cluster8_58prots_ftc238.xlsx")


## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("Nthy_cluster8vsFTC238_cluster8_58prots_nthy_gene.tsv") # ftc238 same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot[duplicated(cluster_prot$P2G),]
# Prot    P2G preferredName
# 35   P50990 P50990          CCT8
# 36 P50990-2 P50990          CCT8
# 37 P50990-2 P50990          CCT8
cluster_prot <- unique(cluster_prot) # 58
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
pdf("Nthy_cluster8vsFTC238_cluster8_58prots_ftc238_ori.pdf", 
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
  gaps_col = 28,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - ftc238 is similar with Nthy
string <- rio::import("Nthy_cluster8vsFTC238_cluster8_58prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
cfms_marked <- unique(cfms_marked) # 140
rio::export(cfms_marked,"Nthy_cluster8vsFTC238_cluster8_58prots_ftc238_final.tsv")

### nthy cluster5 ftc238 cluster 8 nodes mean intensity ----------------------------------------------
node1 <- cfms_marked[,1]
node2 <- cfms_marked[,2]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node <- unique(cf_node) # 25
colnames(cf_node)[1] <- "preferredName" 
cf_node <- left_join(cf_node, cluster_prot, by = "preferredName") # 45
cf_node <- unique(cf_node)
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

protein_list <- cf_node$V1 # 58
auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list) # 58
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") # 58
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster8vsFTC238_cluster8_58prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(FTC238_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","FTC238")
cf_node_ftc238 <- left_join(cf_node, auc_results2, by = "V1") # 58
#cf_node_ftc238$FTC238 <- log10(cf_node_ftc238$FTC238)
colnames(cf_node_ftc238)[1] <- "name"
cf_node_ftc238 <- unique(cf_node_ftc238)
rio::export(cf_node_ftc238, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster8vsFTC238_cluster8_58prots_ftc238_intensity.tsv")

## nthy cluster 8 vs ftc238 cluster 7 prot nthy-------------------------------------------------------
setwd("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238_cluster7vsnthy_cluster8_0830")
current_row <- overlap_unnested[53,] # line 53 
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ",")) # 13
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1) # 13
rio::export(cluster_prot,"Nthy_cluster8vsFTC238_cluster7_13prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster8vsFTC238_cluster7_13prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("Nthy_cluster8vsFTC238_cluster7_13prots_nthy_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
# cluster_prot[duplicated(cluster_prot$Prot),]
cluster_prot[duplicated(cluster_prot$P2G),]
# > cluster_prot[duplicated(cluster_prot$P2G),]
# Prot    P2G preferredName
# 35   P50990 P50990          CCT8
# 36 P50990-2 P50990          CCT8
# 37 P50990-2 P50990          CCT8
cluster_prot <- unique(cluster_prot) #13
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
pdf("Nthy_cluster8vsFTC238_cluster7_13prots_nthy_ori.pdf", 
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
  gaps_col = 28,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("Nthy_cluster8vsFTC238_cluster7_13prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
cfms_marked <- unique(cfms_marked) # 4
rio::export(cfms_marked,"Nthy_cluster8vsFTC238_cluster7_13prots_nthy_final.tsv")

## nthy cluster 8 ftc238 cluster 7 prot ftc238-------------------------------------------------------
# current_row <- overlap[3,]
# proteins_to_plot <- unlist(current_row$proteins)
# subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
# 
# current_row <- overlap_unnested[1,]
# proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1) #13
rio::export(cluster_prot,"Nthy_cluster8vsFTC238_cluster7_13prots_ftc238.tsv")
rio::export(cluster_prot,"Nthy_cluster8vsFTC238_cluster7_13prots_ftc238.tsv")


## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("Nthy_cluster8vsFTC238_cluster7_13prots_nthy_gene.tsv") # ftc238 same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot[duplicated(cluster_prot$P2G),]
# Prot    P2G preferredName
# 35   P50990 P50990          CCT8
# 36 P50990-2 P50990          CCT8
# 37 P50990-2 P50990          CCT8
cluster_prot <- unique(cluster_prot) # 13
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
pdf("Nthy_cluster8vsFTC238_cluster7_13prots_ftc238_ori.pdf", 
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
  gaps_col = 35,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - ftc238 is similar with Nthy
string <- rio::import("Nthy_cluster8vsFTC238_cluster7_13prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
cfms_marked <- unique(cfms_marked) # 140
rio::export(cfms_marked,"Nthy_cluster8vsFTC238_cluster7_13prots_ftc238_final.tsv")

### nthy cluster 8 ftc238 cluster 7 nodes mean intensity ----------------------------------------------
node1 <- cfms_marked[,1]
node2 <- cfms_marked[,2]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node <- unique(cf_node) 
colnames(cf_node)[1] <- "preferredName" 
cf_node <- left_join(cf_node, cluster_prot, by = "preferredName") 
cf_node <- unique(cf_node)
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

protein_list <- cf_node$V1 # 6
auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list) # 6
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") # 6
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster8vsFTC238_cluster7_13prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(FTC238_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","FTC238")
cf_node_ftc238 <- left_join(cf_node, auc_results2, by = "V1") # 6
#cf_node_ftc238$FTC238 <- log10(cf_node_ftc238$FTC238)
colnames(cf_node_ftc238)[1] <- "name"
cf_node_ftc238 <- unique(cf_node_ftc238)
rio::export(cf_node_ftc238, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster8vsFTC238_cluster7_13prots_ftc238_intensity.tsv")


## nthy cluster 2vs ftc238 cluster 2 prot nthy-------------------------------------------------------
setwd("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238vsnthy_cluster2_0830")
current_row <- overlap_unnested[8,] # line 8
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ",")) # 58
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1) # 58
rio::export(cluster_prot,"Nthy_cluster2vsFTC238_cluster2_58prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster2vsFTC238_cluster2_58prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("Nthy_cluster2vsFTC238_cluster2_58prots_nthy_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
# cluster_prot[duplicated(cluster_prot$Prot),]
cluster_prot[duplicated(cluster_prot$P2G),]
# > cluster_prot[duplicated(cluster_prot$P2G),]
# Prot    P2G preferredName
# 35   P50990 P50990          CCT8
# 36 P50990-2 P50990          CCT8
# 37 P50990-2 P50990          CCT8
cluster_prot <- unique(cluster_prot) # 58
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
pdf("Nthy_cluster2vsFTC238_cluster2_58prots_nthy_ori.pdf", 
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
  gaps_col = 20,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("Nthy_cluster2vsFTC238_cluster2_58prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
cfms_marked <- unique(cfms_marked) # 112
rio::export(cfms_marked,"Nthy_cluster2vsFTC238_cluster2_58prots_nthy_final.tsv")

## nthy cluster 2 ftc238 cluster 2 prot ftc238-------------------------------------------------------
# current_row <- overlap[3,]
# proteins_to_plot <- unlist(current_row$proteins)
# subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
# 
# current_row <- overlap_unnested[1,]
# proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1) # 58
rio::export(cluster_prot,"Nthy_cluster2vsFTC238_cluster2_58prots_ftc238.tsv")
rio::export(cluster_prot,"Nthy_cluster2vsFTC238_cluster2_58prots_ftc238.tsv")


## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("Nthy_cluster2vsFTC238_cluster2_58prots_nthy_gene.tsv") # ftc238 same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot[duplicated(cluster_prot$P2G),]

cluster_prot <- unique(cluster_prot) # 58
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
pdf("Nthy_cluster2vsFTC238_cluster2_58prots_ftc238_ori.pdf", 
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
  gaps_col = 20,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - ftc238 is similar with Nthy
string <- rio::import("Nthy_cluster2vsFTC238_cluster2_58prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
cfms_marked <- unique(cfms_marked) # 112
rio::export(cfms_marked,"Nthy_cluster2vsFTC238_cluster2_58prots_ftc238_final.tsv")

### nthy cluster 2 ftc238 cluster 2 nodes mean intensity ----------------------------------------------
node1 <- cfms_marked[,1]
node2 <- cfms_marked[,2]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node <- unique(cf_node) 
colnames(cf_node)[1] <- "preferredName" 
cf_node <- left_join(cf_node, cluster_prot, by = "preferredName") 
cf_node <- unique(cf_node) # 49
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

protein_list <- cf_node$V1 #  49
auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list) # 49
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") 
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster2vsFTC238_cluster2_58prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(FTC238_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","FTC238")
cf_node_ftc238 <- left_join(cf_node, auc_results2, by = "V1") # 49
#cf_node_ftc238$FTC238 <- log10(cf_node_ftc238$FTC238)
colnames(cf_node_ftc238)[1] <- "name"
cf_node_ftc238 <- unique(cf_node_ftc238)
rio::export(cf_node_ftc238, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster2vsFTC238_cluster2_58prots_ftc238_intensity.tsv")


## nthy cluster 5 vs ftc238 cluster 6 prot nthy-------------------------------------------------------
setwd("G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\ftc238_cluster6vsnthy_cluster5_0830")
current_row <- overlap_unnested[30,] # line 30
proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ",")) # 32
subset_data <- heatmap[heatmap$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1) # 32
rio::export(cluster_prot,"Nthy_cluster5vsFTC238_cluster6_32prots_nthy.tsv")
rio::export(cluster_prot,"Nthy_cluster5vsFTC238_cluster6_32prots_nthy.xlsx")

## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
## string p2g
gene <- rio::import("Nthy_cluster5vsFTC238_cluster6_32prots_nthy_gene.tsv")
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot87,"cluster_87prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
# cluster_prot[duplicated(cluster_prot$Prot),]
cluster_prot[duplicated(cluster_prot$P2G),]
# > cluster_prot[duplicated(cluster_prot$P2G),]
# Prot    P2G preferredName
# 35   P50990 P50990          CCT8
# 36 P50990-2 P50990          CCT8
# 37 P50990-2 P50990          CCT8
cluster_prot <- unique(cluster_prot) # 32
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
pdf("Nthy_cluster5vsFTC238_cluster6_32prots_nthy_ori.pdf", 
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
  gaps_col = 20,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results
string <- rio::import("Nthy_cluster5vsFTC238_cluster6_32prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
cfms_marked <- unique(cfms_marked) # 59
rio::export(cfms_marked,"Nthy_cluster5vsFTC238_cluster6_32prots_nthy_final.tsv")

## nthy cluster 5 ftc238 cluster 6 prot ftc238-------------------------------------------------------
# current_row <- overlap[3,]
# proteins_to_plot <- unlist(current_row$proteins)
# subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
# 
# current_row <- overlap_unnested[1,]
# proteins_to_plot  <- unlist(strsplit(current_row$`protein list`, ","))
subset_data <- heatmap2[heatmap2$V1 %in% proteins_to_plot, ]
subset_data1 <- subset_data[,2:61]
rownames(subset_data1) <- subset_data$V1
cluster_prot <- as.data.frame(subset_data$V1) # 59
rio::export(cluster_prot,"Nthy_cluster5vsFTC238_cluster6_32prots_ftc238.tsv")
rio::export(cluster_prot,"Nthy_cluster5vsFTC238_cluster6_32prots_ftc238.xlsx")


## rename col and row
minutes <- 9
seconds <- seq(0, by = 19, length.out = 60)  

total_seconds <- minutes * 60 + seconds
new_row_names <- sprintf("%d min %02ds", total_seconds %/% 60, total_seconds %% 60)

colnames(subset_data)[2:61] <- new_row_names
cluster_prot <- as.data.frame(subset_data$V1)
## string p2g
gene <- rio::import("Nthy_cluster5vsFTC238_cluster6_32prots_nthy_gene.tsv") # ftc238 same with nthy
gene <- gene[,c(2,4)]
colnames(cluster_prot)[1] <- "Prot"
#rio::export(cluster_prot80,"cluster_80prot.tsv")
cluster_prot$P2G <- gsub("-(.*)","",cluster_prot$Prot)
colnames(gene)[1] <- "P2G"
cluster_prot <- left_join(cluster_prot,gene, by = "P2G")
cluster_prot[duplicated(cluster_prot$P2G),]

cluster_prot <- unique(cluster_prot) # 32
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
pdf("Nthy_cluster5vsFTC238_cluster6_32prots_ftc238_ori.pdf", 
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
  gaps_col = 10,
  cellheight = 4.6,  
  # fontfamily = "Arial",       
  fontsize = 5.5,              
  fontsize_row = 5,           
  fontsize_col = 5.5)

dev.off()

## string mapping results - ftc238 is similar with Nthy
string <- rio::import("Nthy_cluster5vsFTC238_cluster6_32prots_nthy_string.tsv")
string <- string[,1:2]
colnames(string)[1:2] <- c("node1","node2")
cfms_interactions <- depps_nf5 %>%
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
cfms_marked <- unique(cfms_marked) # 59
rio::export(cfms_marked,"Nthy_cluster5vsFTC238_cluster6_32prots_ftc238_final.tsv")

### nthy cluster 5 ftc238 cluster 6 nodes mean intensity ----------------------------------------------
node1 <- cfms_marked[,1]
node2 <- cfms_marked[,2]
colnames(node2)[1] <- "node1"
cf_node <- unique(rbind(node1,node2))
cf_node <- unique(cf_node) 
colnames(cf_node)[1] <- "preferredName" 
cf_node <- left_join(cf_node, cluster_prot, by = "preferredName") 
cf_node <- unique(cf_node) # 49
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

protein_list <- cf_node$V1 #  49
auc_results <- calculate_curve_auc(Nthy_nor_mean, protein_list) # 49
colnames(auc_results) <- c("V1","Nthy")
cf_node_nthy <- left_join(cf_node, auc_results, by = "V1") 
colnames(cf_node_nthy)[1] <- "name"
rio::export(cf_node_nthy, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster5vsFTC238_cluster6_32prots_nthy_intensity.tsv")

auc_results2 <- calculate_curve_auc(FTC238_nor_mean, protein_list)
colnames(auc_results2) <-  c("V1","FTC238")
cf_node_ftc238 <- left_join(cf_node, auc_results2, by = "V1") # 49
#cf_node_ftc238$FTC238 <- log10(cf_node_ftc238$FTC238)
colnames(cf_node_ftc238)[1] <- "name"
cf_node_ftc238 <- unique(cf_node_ftc238)
rio::export(cf_node_ftc238, "G:\\Tower3\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\DEPPS\\Nthy_cluster5vsFTC238_cluster6_32prots_ftc238_intensity.tsv")

# for (i in 1:nrow(overlap)) {
#   current_row <- overlap[i, ]
#   proteins_to_plot <- unlist(current_row$proteins)  # protein list
#   subset_data <- heatmap_clean[rownames(heatmap_clean) %in% proteins_to_plot, ]
#   
#   p <- pheatmap(
#     subset_data,
#     scale = "row",
#     clustering_distance_rows = "euclidean",
#     clustering_method = "ward.D2",
#     cluster_rows = TRUE,
#     cluster_cols = FALSE,
#     treeheight_row = 0,
#     show_rownames = TRUE,
#     main = paste(
#       "Overlap: Map1 Cluster", current_row$cluster_map1,
#       "vs Map2 Cluster", current_row$cluster_map2,
#       "| Proteins =", current_row$n_proteins
#     )
#   )
#   
#  
#   filename <- paste0(
#     "heatmap_Map1_Cluster", current_row$cluster_map1,
#     "_Map2_Cluster", current_row$cluster_map2, ".pdf"
#   )
#   png(filename, width = 800, height = 600)
#   print(p)
#   dev.off()
# }

# ## nt--Nthy max col check
# nt_sep <- nt  %>%
#   separate(protein_pair, into = c("protein_1", "protein_2"), sep = "_")
# colnames(nt_sep)
# colnames(nt_sep)[2] <- "protein_2_nt_sep"
# 
# depps_6 <- left_join(depps_nt4,nt_sep, by = "protein_1")
# colnames(depps_6)[18] <- "protein_1_max_col"
# depps_6 <- depps_6[,c("protein_1","protein_2","fold_change","node1","node2","data_file","protein_1_max_col")]
# depps_6 <- unique(depps_6)
# 
# ## amend protein1
# depps_6$protein_1_col_num <- gsub("(.*)_","",depps_6$protein_1_max_col)
# depps_6$cell <- gsub("_(.*)","",depps_6$protein_1_max_col)
# depps_6$cell <- gsub(".1","",depps_6$cell)
# depps_6_nthy <- depps_6[which(depps_6$cell=="Nthy"),]
# 
# 
# colnames(nt_sep)[1] <- "protein_1_nt_sep"
# colnames(nt_sep)[2] <- "protein_2"
# 
# depps_6 <- left_join(depps_6,nt_sep, by = "protein_2")
# colnames(depps_6)[13] <- "protein_2_max_col"
# depps_6$protein_2_col_num <- gsub("(.*)_","",depps_6$protein_2_max_col)
# depps_6 <- depps_6[,c("protein_1","protein_2","fold_change","node1","node2","data_file.x","data_file.y","protein_label","protein_1_max_col",
#                       "protein_1_col_num","protein_2_max_col","protein_2_col_num")]
# depps_6 <- unique(depps_6)
## nt comb
# nt_comb_final <- rbind(nt2_edge_string,depps_nt5)
# nt_comb_final <- nt_comb_final[,c(1:3)]
# 
# nt_comb_final <- nt_comb_final %>%
#   mutate(
#     node_min = pmin(node1, node2),
#     node_max = pmax(node1, node2),
#     edge_id = paste(node_min, node_max, sep = "|"))
# 
# common_edges <- nt_comb_final %>%
#   group_by(edge_id) %>%
#   filter(n_distinct(edge) == 2) %>%  # BOTH
#   ungroup() %>%
#   distinct(edge_id)
# 
# both_edges <- nt_comb_final %>%
#   filter(edge_id %in% common_edges$edge_id) %>%
#   group_by(edge_id) %>%
#   summarise(
#     node1 = first(node_min),
#     node2 = first(node_max),
#     edge = "Both") %>%
#   select(-edge_id)
# 
# other_edges <- nt_comb_final %>%
#   filter(!edge_id %in% common_edges$edge_id) %>%
#   select(node1, node2, edge)
# 
# final_df <- bind_rows(both_edges, other_edges) %>%
#   arrange(node1, node2)
# 
# final_df <- final_df %>% select(node1, node2, edge)
# rio::export(final_df,"E://MenggeLYU//NC_rebuttal_t3//New_feature_RF//ComplexXGB//DEPPS//tpcvsnthy_edge_20250714-2.tsv")
# 
# final_node1 <- depps_nt5[,c(1,3)]
# colnames(final_node1)[1] <- "node"
# final_node2 <- depps_nt5[,c(2,3)]
# colnames(final_node2)[1] <- "node"
# final_node <- rbind(final_node1,final_node2)
# final_node <- unique(final_node)
# final_node3 <- final_node %>%
#   group_by(node) %>%
#   filter(!all(is.na(fold_change))) %>%  
#   mutate(abs_fc = ifelse(is.na(fold_change), NA, abs(fold_change))) %>%
#   arrange(desc(abs_fc), .by_group = TRUE) %>%
#   slice(1) %>%
#   select(-abs_fc) %>%
#   ungroup()
# rio::export(final_node3,"tpcvsnthy_node_20250714-3.tsv")
