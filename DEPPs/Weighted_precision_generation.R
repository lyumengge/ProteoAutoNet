##################################################
## Project: ProteoAutoNet
## Script purpose: Weighted_precision_NC-rebuttal
## Date: 2025-06-30
## Author: MenggeLYU
## Version: 1.0
##################################################

# path and library --------------------------------------------------------
setwd("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\Network_fig")
library(data.table)
library(rio)
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggvenn)
library(ggrepel)
library(tidyverse)
# load data ---------------------------------------------------------------
TPC <- rio::import("TPC_prob095_final_20250628.tsv")
FTC238 <- rio::import("FTC238_prob095_final_20250628.tsv")
Nthy <- rio::import("Nthy_prob095_final_20250628.tsv")


# Weighted cut off  -------------------------------------------------------
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

legend("topright", 
       legend = c("TPC", "Nthy", "FTC238"),
       col = c("blue", "red", "green"),
       lty = 1, lwd = 2,
       cex = 0.8)

#FTC238 <- FTC238[,1:2]
FTC238  <- unique(FTC238 )
FTC238 $cell <- "FTC238"

#TPC <- TPC[,1:2]
TPC <- unique(TPC)
TPC$cell <- "TPC"

#Nthy <- Nthy[,1:2]
Nthy <- unique(Nthy)
Nthy$cell <- "Nthy"

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


# Overlap three databases -------------------------------------------------
# Three databases overlap
result <- which(adj_df== 1, arr.ind = TRUE) %>%
  as.data.frame() %>%
  mutate(
    protein_A = rownames(adj_df)[row],
    protein_B = colnames(adj_df)[col]
  ) %>%
  filter(protein_A != protein_B) %>%  
  select(protein_A, protein_B)

print(result)

# Overlapped number in each zone
# All count
overlap_counts <- list(
  "TPC_only" = setdiff(pair_list$TPC, union(pair_list$Nthy, pair_list$FTC238)) %>% length(),
  "Nthy_only" = setdiff(pair_list$Nthy, union(pair_list$TPC, pair_list$FTC238)) %>% length(),
  "FTC238_only" = setdiff(pair_list$FTC238, union(pair_list$TPC, pair_list$Nthy)) %>% length(),
  "TPC_Nthy" = length(intersect(pair_list$TPC, pair_list$Nthy)) - length(Reduce(intersect, pair_list)),
  "TPC_FTC238" = length(intersect(pair_list$TPC, pair_list$FTC238)) - length(Reduce(intersect, pair_list)),
  "Nthy_FTC238" = length(intersect(pair_list$Nthy, pair_list$FTC238)) - length(Reduce(intersect, pair_list)),
  "All_three" = length(Reduce(intersect, pair_list)))

overlap_df <- data.frame(
  Region = names(overlap_counts),
  Count = unlist(overlap_counts))

# All pairs
overlap_pairs <- list(
  "TPC_only" = setdiff(pair_list$TPC, union(pair_list$Nthy, pair_list$FTC238)),
  "Nthy_only" = setdiff(pair_list$Nthy, union(pair_list$TPC, pair_list$FTC238)),
  "FTC238_only" = setdiff(pair_list$FTC238, union(pair_list$TPC, pair_list$Nthy)),
  "TPC_Nthy_only" = setdiff(intersect(pair_list$TPC, pair_list$Nthy), pair_list$FTC238),
  "TPC_FTC238_only" = setdiff(intersect(pair_list$TPC, pair_list$FTC238), pair_list$Nthy),
  "Nthy_FTC238_only" = setdiff(intersect(pair_list$Nthy, pair_list$FTC238), pair_list$TPC),
  "All_three" = Reduce(intersect, pair_list))


# Draw overlap ------------------------------------------------------------
# library(ggvenn)
# ggvenn(pair_list,
#        fill_color = c("blue", "red", "green"),
#        stroke_size = 0.5,
#        set_name_size = 4,
#        show_percentage = TRUE,  
#        text_size = 5,            
#        label_sep = "\n") +       
#   labs(title = "Protein Pair Overlap") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# result$pair <- apply(result[, c("protein_A", "protein_B")], 1, 
#                      function(x) paste(sort(x), collapse = "_"))
# venn_regions <- list(
#   "TPC_only" = setdiff(pair_list$TPC, union(pair_list$Nthy, pair_list$FTC238)),
#   "Nthy_only" = setdiff(pair_list$Nthy, union(pair_list$TPC, pair_list$FTC238)),
#   "FTC238_only" = setdiff(pair_list$FTC238, union(pair_list$TPC, pair_list$Nthy)),
#   "TPC_Nthy_only" = setdiff(intersect(pair_list$TPC, pair_list$Nthy), pair_list$FTC238),
#   "TPC_FTC238_only" = setdiff(intersect(pair_list$TPC, pair_list$FTC238), pair_list$Nthy),
#   "Nthy_FTC238_only" = setdiff(intersect(pair_list$Nthy, pair_list$FTC238), pair_list$TPC),
#   "All_three" = Reduce(intersect, pair_list))
# 
# overlap_results <- sapply(venn_regions, function(x) {
#   sum(x %in% result$pair)})
# 
# overlap_df <- data.frame(
#   Region = names(venn_regions),
#   Total_in_Venn = lengths(venn_regions),
#   In_results = overlap_results,
#   Percentage = round(overlap_results / lengths(venn_regions) * 100, 1))
# 
# 
# overlap_df <- overlap_df %>% arrange(desc(In_results))
# 
# print(overlap_df)

# Bar plot
library(ggplot2)
ggplot(overlap_df, aes(x = reorder(Region, -In_results), y = In_results)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(In_results, "\n(", Percentage, "%)")), 
            vjust = -0.5, size = 3) +
  labs(x = "Venn Diagram Region", 
       y = "Number of Pairs in Results",
       title = "Overlap between Venn Regions and Results") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Generate three cell line all interactions with string -------------------
Fshort <- FTC238[,1:2]
Fshort <- unique(Fshort)
Fshort$label <- "FTC238"

Tshort <- TPC[,1:2]
Tshort <- unique(Tshort)
Tshort$label <- "TPC"

Nshort <- Nthy[,1:2]
Nshort <- unique(Nshort)
Nshort$label <- "Nthy"

comb <- rbind(Fshort,Tshort,Nshort)
comb_pairs <- unique(comb[,1:2])
protein <-unique(c(comb_pairs$protein_A,comb_pairs$protein_B)) # 4603
protein_df <- as.data.frame(protein) # 4603

colnames(protein_df) <- "Prot"
rio::export(protein_df, "string_protein.tsv")

# STRING mapping ----------------------------------------------------------
string_g2p <- rio::import("E:\\MenggeLYU\\CHINDA\\STRING\\final\\STRING_protein_list_withlabel_MGL_20240611.tsv")
protein_df2g <- left_join(protein_df,string_g2p, by ="Prot")
match <- protein_df2g[which(is.na(protein_df2g$preferred_name)),]

# unmatched protein
prot2 <- as.data.frame(match$Prot)
rio::export(prot2, "string_protein_nolabel.tsv")

# string results
stringmap <- rio::import("string_mapping_unmatched_proteins.tsv")
second <- stringmap[,c(2,4)]
colnames(second) <- c("Prot", "Gene")
first <- protein_df2g[-which(is.na(protein_df2g$preferred_name)),]
first <- first[,c(1,3)]
colnames(first) <- c("Prot", "Gene")
g2p <- rbind(first,second)
gene <- unique(g2p$Gene) # 3581
rio::export(g2p,"string_gene.tsv")

edge <- rio::import("ensp_id.xlsx")
node <- rio::import("STRING network default node.csv")
node <- node[,c(15,19)]
colnames(node)[2] <- "A"
edge1 <- left_join(edge,node, by = "A")
colnames(edge1)[3]<-'gene_A' 
colnames(node)[2] <- "B"
edge2 <- left_join(edge1,node, by = "B")
colnames(edge2)[4]<-'gene_B' 
colnames(g2p)[2] <-"gene_A" 
edge3 <- left_join(edge2,g2p, by = "gene_A")
colnames(edge3)[5]<-'protein_A' 
colnames(g2p)[2] <-"gene_B" 
edge4 <- left_join(edge3,g2p, by = "gene_B")
colnames(edge4)[6]<-'protein_B' 
STRING_FINAL <- edge4[,5:6]
STRING_FINAL <- unique(STRING_FINAL) # 44658
STRING_FINAL$label <- 'STRING'
# Three all interactions + STRING
COMB_ALL <- rbind(comb,STRING_FINAL)
table(COMB_ALL$label) 
COMB_ALL$pair <- apply(COMB_ALL[, c("protein_A", "protein_B")], 1, function(x) paste(sort(x), collapse = "_"))
# FTC238      N STRING    TPC 
# 1699   2798  44658   2800 
rio::export(COMB_ALL,"network_pair_all.tsv") # Large network - give up

COMB_ALL$color <- 0

pair_labels <- split(COMB_ALL$label, COMB_ALL$pair)
for (pair in names(pair_labels)) {
  labels <- unique(pair_labels[[pair]])
  
  if (all(c("N", "TPC", "FTC238") %in% labels)) {
    color_val <- 7
  } else if (all(c("TPC", "FTC238") %in% labels)) {
    color_val <- 6
  } else if (all(c("N", "FTC238") %in% labels)) {
    color_val <- 5
  } else if (all(c("N", "TPC") %in% labels)) {
    color_val <- 4
  } else if ("FTC238" %in% labels) {
    color_val <- 3
  } else if ("TPC" %in% labels) {
    color_val <- 2
  } else if ("N" %in% labels) {
    color_val <- 1
  }
  
  COMB_ALL$color[COMB_ALL$pair == pair] <- color_val}

table(COMB_ALL$color)
# 1     2     3     4     5     6     7 
# 16264 16349  8883  3761  1811  1466  3421 

# add edge color
COMB_ALL <- COMB_ALL %>%
  mutate(
    resource = case_when(
      label %in% c("FTC238", "N", "TPC") ~ "CF-MS",
      label == "STRING" ~ "STRING",
      TRUE ~ NA_character_))

table(COMB_ALL$resource)

# CF-MS STRING 
# 7297  44658 

all_protein <- as.data.frame(unique(c(COMB_ALL$protein_A,COMB_ALL$protein_B)))
rio::export(COMB_ALL, "thyroid_cell_lineandstring_protein_edge.tsv")
colnames(all_protein)[1] <- "Prot"
all_cano <- as.data.frame(all_protein[-grep("-",all_protein$Prot),])
colnames(all_cano) <- "Prot"
all_iso <- as.data.frame(all_protein[grep("-",all_protein$Prot),])
colnames(all_iso) <- "Prot"
all_iso$Prot_clean <- gsub("-(.*)","",all_iso$Prot)
rio::export(all_cano,"thyroid_cell_lineandstring_protein_edge_canonical_prots.tsv")
rio::export(all_iso,"thyroid_cell_lineandstring_protein_edge_isoform_prots.tsv")
# iso results loading
mapping_iso <- rio::import("thyroid_cell_lineandstring_protein_egde_isoform_prots_mapping results.tsv")
colnames(mapping_iso)[2] <- "Prot_clean"
all_iso <- left_join(all_iso,mapping_iso, by = "Prot_clean")
all_iso <- all_iso[,c(1,5)]
colnames(all_iso)[2] <- "Gene"
all_iso <- unique(all_iso) #612
mapping_cano <- rio::import("thyroid_cell_lineandstring_protein_egde_canonical_prots_mapping results.tsv.tsv")
colnames(mapping_cano)[2] <- "Prot"
all_cano <- left_join(all_cano,mapping_cano, by = "Prot")
all_cano <- all_cano[,c(1,4)]
all_cano<- unique(all_cano) #1870
colnames(all_cano) <- c("Prot","Gene")
is.na(all_cano$Gene)
all <- rbind(all_cano,all_iso)
all2 <- all[-which(is.na(all$Gene)),] #2465
colnames(all2)[1] <- "protein_A"
COMB_ALL1 <- left_join(COMB_ALL,all2, by = 'protein_A')
colnames(COMB_ALL1)[7] <- "gene_A"
colnames(all2)[1] <- "protein_B"
COMB_ALL2 <- left_join(COMB_ALL1,all2, by = 'protein_B')
colnames(COMB_ALL2)[8] <- "gene_B"
COMB_ALL2 <- unique(COMB_ALL2) #51955
rio::export(COMB_ALL2,"thyroid_cell_lineandstring_protein_egde_all_columns.tsv")
gAB <- unique(COMB_ALL2[,7:8]) #50235
dup_rows <- COMB_ALL2[duplicated(COMB_ALL2[, 7:8]) | duplicated(COMB_ALL2[, 7:8], fromLast = TRUE), ]
forstring <- COMB_ALL2[,c(3,5:8)]
forstring_unique <- unique(forstring)
forstring_unique2 <- forstring_unique %>%
  group_by(gene_A, gene_B) %>%
  mutate(
    e_color = case_when(
      all(resource == "CF-MS") ~ 1,
      all(resource == "STRING") ~ 2,
      "CF-MS" %in% resource & "STRING" %in% resource ~ 3,
      TRUE ~ NA_integer_  
    )
  ) %>%
  ungroup()

forstring_unique3 <- forstring_unique2[-which(is.na(forstring_unique2$gene_A)),]
forstring_unique4 <- forstring_unique3[-which(is.na(forstring_unique3$gene_B)),]
forstring_unique5 <- unique(forstring_unique4[,c(4,5,6,2)]) #50203
forstring_unique_test <- unique(forstring_unique5[,1:2]) #50182

dup_rows <- forstring_unique5 [duplicated(forstring_unique5[, 1:2]) | duplicated(forstring_unique5[, 1:2], fromLast = TRUE), ]
rio::export(dup_rows,"hand_correct.xlsx")

dup_rows_after <- rio::import("hand_correct-after.xlsx")
dup_rows_after <- unique(dup_rows_after) # 21 matched

dup_color_mapping <- dup_rows_after %>%
  select(gene_A, gene_B, e_color, color)

# node color
forstring_unique5_updated <- forstring_unique5 %>%
  left_join(dup_color_mapping, by = c("gene_A", "gene_B", "e_color"), suffix = c("_old", "")) %>%
  mutate(color = ifelse(!is.na(color), color, color_old)) %>%  
  select(-color_old)  

forstring_unique6 <- unique(forstring_unique5_updated) # 50182 correct with forstring_unique_test
rio::export(forstring_unique6, "all_with_gene_for_string.tsv")

# string - many egdes - delete string
string <- rio::import("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\Network_fig\\all_with_gene_for_string.tsv")
string <- string[-which(string$e_color==2),]
string <- string %>%
  mutate(
    node_min = pmin(gene_A, gene_B),
    node_max = pmax(gene_A, gene_B),
    edge_id = paste(node_min, node_max, sep = "|"))

duplicate_edges <- string %>%
  group_by(edge_id) %>%
  filter(n() > 1) %>%  
  arrange(edge_id) %>% 
  ungroup() # zero row
rio::export(string, "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\Network_fig\\all_with_gene_for_withoutstring-20250715.tsv")

# generate node anno
all_genes <- unique(c(forstring_unique6$gene_A, forstring_unique6$gene_B))

node_anno<- tibble(gene = all_genes) %>%
  left_join(
    forstring_unique6 %>% 
      pivot_longer(cols = c(gene_A, gene_B), names_to = "type", values_to = "gene") %>%
      group_by(gene) %>%
      summarize(
        color = if (any(color == 7)) "7" else paste(sort(unique(color)), collapse = ",")
      ),
    by = "gene")

rio::export(node_anno, "node_anno_color.tsv") # amended version: node_anno_color_after.txt
node_anno_without7 <- node_anno[-which(node_anno$color==7),]
table(node_anno_without7$color)


##KEGG 
kegg <- rio::import("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\Network_fig\\all_three_enrichment.KEGG.tsv")
sorted_kegg <- kegg %>%
  arrange(desc(strength))


top10 <- sorted_kegg[1:10,]

ggplot(top10, aes(x = strength, y = -log10(`false discovery rate`),
                  size = `observed gene count`,
                  color = `signal`)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = `term description`), 
                  size = 3.5,
                  box.padding = 0.5,
                  max.overlaps = Inf) +  # 确保所有标签都显示
  scale_color_gradient2(
    low = "#226B9F", 
    mid = "#635B9F", 
    high = "#B43D6E", 
    midpoint = median(top10$signal),
    name = "Signal\nStrength"
  ) +
  scale_size_continuous(
    range = c(3, 8),
    name = "Gene\nCount"
  ) +
  labs(
    x = "Pathway Interaction Strength", 
    y = "-log10(Adjusted p-value)", 
    title = "Top 10 Significant KEGG Pathways"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90"))


# dup_rows2 <- dup_rows[which(dup_rows$color>3),]
# # The following from UniProt:
# all[rownames(na_list)[2],2] <- "SNRPGP15"
# all[rownames(na_list)[2],3] <- "GSTT2"
# all[rownames(na_list)[2],4] <- "GPX1"
# all[rownames(na_list)[2],5] <- "HSP90AA5P"
# all[rownames(na_list)[2],6] <- "ZFTRAF1"
# all[rownames(na_list)[2],7] <- "IGKV3-20"
# all[rownames(na_list)[2],8] <- "IGKV3-20"


## keggtop1 --------------------------------------------------------------
top1 <- top10[1,]
unique_genes <- top1 %>%
  pull(`matching proteins in your network (labels)`) %>%  
  str_split(",") %>%        
  unlist() %>%             
  unique() %>%              
  sort() %>%               
  as.data.frame() %>%       
  setNames("gene_symbol")   

gene_list <- unique_genes$gene_symbol  
gene_pattern <- paste0(
  "^((", paste(gene_list, collapse="|"), ")\\|[^|]+$)|",  
  "(^[^|]+\\|(", paste(gene_list, collapse="|"), ")$)")

matched_edges <- string %>%
  filter(str_detect(edge_id, gene_pattern))
matched_edges_safe <- string %>%
  separate(edge_id, into = c("gene1", "gene2"), sep = "\\|", remove = FALSE) %>%
  filter(gene1 %in% gene_list | gene2 %in% gene_list) %>%
  select(-gene1, -gene2)
rio::export(matched_edges_safe,"E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\Network_fig\\kegg_top1.tsv")

## keggtop2 --------------------------------------------------------------
top2 <- top10[2,]
unique_genes <- top2 %>%
  pull(`matching proteins in your network (labels)`) %>%  
  str_split(",") %>%        
  unlist() %>%             
  unique() %>%              
  sort() %>%               
  as.data.frame() %>%       
  setNames("gene_symbol")   

gene_list <- unique_genes$gene_symbol  
gene_pattern <- paste0(
  "^((", paste(gene_list, collapse="|"), ")\\|[^|]+$)|",  
  "(^[^|]+\\|(", paste(gene_list, collapse="|"), ")$)")

matched_edges <- string %>%
  filter(str_detect(edge_id, gene_pattern))
matched_edges_safe <- string %>%
  separate(edge_id, into = c("gene1", "gene2"), sep = "\\|", remove = FALSE) %>%
  filter(gene1 %in% gene_list | gene2 %in% gene_list) %>%
  select(-gene1, -gene2)
rio::export(matched_edges_safe,"E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\Network_fig\\kegg_top2.tsv")

## keggtop3 --------------------------------------------------------------
top3 <- top10[3,]
unique_genes <- top3 %>%
  pull(`matching proteins in your network (labels)`) %>%  
  str_split(",") %>%        
  unlist() %>%             
  unique() %>%              
  sort() %>%               
  as.data.frame() %>%       
  setNames("gene_symbol")   

gene_list <- unique_genes$gene_symbol  
gene_pattern <- paste0(
  "^((", paste(gene_list, collapse="|"), ")\\|[^|]+$)|",  
  "(^[^|]+\\|(", paste(gene_list, collapse="|"), ")$)")

matched_edges <- string %>%
  filter(str_detect(edge_id, gene_pattern))
matched_edges_safe <- string %>%
  separate(edge_id, into = c("gene1", "gene2"), sep = "\\|", remove = FALSE) %>%
  filter(gene1 %in% gene_list | gene2 %in% gene_list) %>%
  select(-gene1, -gene2)
rio::export(matched_edges_safe,"E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\Network_fig\\kegg_top3.tsv")

## keggtop4 --------------------------------------------------------------
top4 <- top10[4,]
unique_genes <- top4 %>%
  pull(`matching proteins in your network (labels)`) %>%  
  str_split(",") %>%        
  unlist() %>%             
  unique() %>%              
  sort() %>%               
  as.data.frame() %>%       
  setNames("gene_symbol")   

gene_list <- unique_genes$gene_symbol  
gene_pattern <- paste0(
  "^((", paste(gene_list, collapse="|"), ")\\|[^|]+$)|",  
  "(^[^|]+\\|(", paste(gene_list, collapse="|"), ")$)")

matched_edges <- string %>%
  filter(str_detect(edge_id, gene_pattern))
matched_edges_safe <- string %>%
  separate(edge_id, into = c("gene1", "gene2"), sep = "\\|", remove = FALSE) %>%
  filter(gene1 %in% gene_list | gene2 %in% gene_list) %>%
  select(-gene1, -gene2)
rio::export(matched_edges_safe,"E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\Network_fig\\kegg_top4.tsv")

## keggtop5 --------------------------------------------------------------
top5 <- top10[5,]
unique_genes <- top5 %>%
  pull(`matching proteins in your network (labels)`) %>%  
  str_split(",") %>%        
  unlist() %>%             
  unique() %>%              
  sort() %>%               
  as.data.frame() %>%       
  setNames("gene_symbol")   

gene_list <- unique_genes$gene_symbol  
gene_pattern <- paste0(
  "^((", paste(gene_list, collapse="|"), ")\\|[^|]+$)|",  
  "(^[^|]+\\|(", paste(gene_list, collapse="|"), ")$)")

matched_edges <- string %>%
  filter(str_detect(edge_id, gene_pattern))
matched_edges_safe <- string %>%
  separate(edge_id, into = c("gene1", "gene2"), sep = "\\|", remove = FALSE) %>%
  filter(gene1 %in% gene_list | gene2 %in% gene_list) %>%
  select(-gene1, -gene2)
rio::export(matched_edges_safe,"E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\Network_fig\\kegg_top5.tsv")

# Three cell lines overlapped ---------------------------------------------
comb$protein_pair <- paste(comb$protein_A, comb$protein_B, sep = "_")
grouped <- comb %>%
  group_by(protein_pair, label) %>%
  summarise(count = n(), .groups = 'drop') # label the detected times
all_labels <- c("FTC238", "Nthy", "TPC")
filtered_pairs <- grouped %>%
  group_by(protein_pair) %>%
  filter(all(all_labels %in% label)) %>%
  select(protein_pair) %>%
  distinct()
threetimes <- comb %>%
  filter(protein_pair %in% filtered_pairs$protein_pair) # 2415 rows
threetimes <- select(threetimes, -protein_pair)
result <- threetimes
result$label <- "CF-MS"
result <- unique(result)
rio::export(result,"network_pair_three.tsv")
# test <- result %>%
#   anti_join(result1, by = names(result)) # old version is not ok, caused by the wrong separator "-" not "_"
three_prot <- unique(c(result$protein_A,result$protein_B))
three_prot2 <- gsub("-(.*)","",three_prot)
three_prot3 <- as.data.frame(unique(three_prot2))
colnames(three_prot3)[1] <- "Prot"
rio::export(three_prot3,"network_pair_three_prot.tsv")

# string mapping ----------------------------------------------------------
string_mapping1 <- rio::import("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\Network_fig\\network_pair_three_prot_mapping_gene.tsv")
string_mapping1 <- string_mapping1[,c(2,4)]

result$protein_A_mapping <- gsub("-(.*)","",result$protein_A)
result$protein_B_mapping <- gsub("-(.*)","",result$protein_B)

colnames(string_mapping1) <- c("protein_A_mapping", "#node1")

result <- left_join(result,string_mapping1, by = 'protein_A_mapping')

colnames(string_mapping1) <- c("protein_B_mapping", "#node2")

result <- left_join(result,string_mapping1, by = 'protein_B_mapping')

string_mapping_edge <- rio::import("string_interactions_edge.tsv")

string_mapping_edge <- string_mapping_edge[,1:2]

string_mapping_edge <- unique(string_mapping_edge)

string_mapping_edge$label <- "STRING"

cfms <- result[,c(6,7,3)]

colnames(cfms) <- c("node1","node2","label")
colnames(string_mapping_edge) <- c("node1","node2","label")
three_common <- rbind(cfms,string_mapping_edge)
three_common <- unique(three_common)

three_common <- three_common %>%
  group_by(node1, node2) %>%
  mutate(new_label = ifelse(n() > 1, "Both", label)) %>%
  ungroup() %>%
  distinct(node1, node2, .keep_all = TRUE) %>%
  select(node1, node2, label = new_label)

three_common <- three_common %>%
  mutate(
    sorted_node1 = pmin(node1, node2),
    sorted_node2 = pmax(node1, node2)
  ) %>%
  group_by(sorted_node1, sorted_node2) %>%
  mutate(
    new_label = ifelse(n_distinct(label) > 1, "Both", label)) %>%
  ungroup() %>%
  distinct(sorted_node1, sorted_node2, .keep_all = TRUE) %>%
  select(node1, node2, label = new_label)

# b <- unique(three_common[,c(1,2)]) # 2431 rows
rio::export(three_common,"three_common_edge.tsv")
table(three_common$label)

three_common <- rio::import("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\Network_fig\\three_common_edge.tsv")
# string <- left_join(string_mapping_edg,string_mapping1)
# colnames(string)[3] <- "protein_A"
# colnames(string_mapping1) <- c("protein", "node2")
# string2 <- left_join(string,string_mapping1, by = "node2")
# string2 <- string2[,3:4]
# colnames(string2)[2] <- "protein_B"
# string2 <- unique(string2)
# string2$label <- "STRING"
# comb_final2 <- rbind(result,string2)
# rio::export(comb_final2,"network_pair_three_WITH STRING.tsv")
# table(comb_final2$label)
# cf_ms_rows <- comb_final2 %>%filter(label == "CF-MS")
# string_rows <- comb_final2 %>%
#   filter(label == "STRING")
# result_anti <- cf_ms_rows %>%
#   anti_join(string_rows, by = c("protein_A", "protein_B"))
# 
# result_string_uni <- string_rows %>%
#   anti_join(cf_ms_rows, by = c("protein_A", "protein_B")) # in string not in cf-ms



# table(comb_final2$label)
# # CF-MS STRING 
# # 747    726 


# comb_final2 <- comb_final2 %>%
#   group_by(protein_A, protein_B) %>%
#   mutate(
#     num = case_when(
#       all(label == "CF-MS") ~ 1,
#       all(label == "STRING") ~ 2,
#       "CF-MS" %in% label & "STRING" %in% label ~ 3,
#       TRUE ~ NA_integer_
#     )
#   ) %>%
#   ungroup()
# 
# head(string_mapping1)
# 
# colnames(string_mapping1)[1] <- "protein_A"
# colnames(string_mapping1)[2] <- "gene"
# 
# comb3 <- left_join(comb_final2,string_mapping1, by = "protein_A")
# colnames(comb3)[5] <- "gene_A"


# 
# colnames(string_mapping1)[1] <- "protein_B"
# colnames(string_mapping1)[2] <- "gene_B"
# 
# comb4 <- left_join(comb3,string_mapping1, by = "protein_B")
# 
# table(comb$label)
# comb$label <- gsub("N", "Nthy", comb$label)
# 
# colnames(comb)[3] <- "cell"
# 
# result <- comb4 %>%
#   inner_join(comb, by = c("protein_A", "protein_B"))
# 
# 
# result <- result[,1:7]
# result <- unique(result)
# 
# 
# rio::export(result,"Thyroid_cell_line_20250630.tsv")
# gene_20250630 <- rio::import("final_mapping_20250630.txt")
# gene_20250630 <- gene_20250630[,c(4:6)]
# gene_20250630 <- unique(gene_20250630)
# rio::export(gene_20250630,"Thyroid_cell_line_gene_20250630.tsv")


