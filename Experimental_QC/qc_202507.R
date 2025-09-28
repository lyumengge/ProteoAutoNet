##################################################
## Project: NC_rebuttal
## Script purpose:JAKA-qc
## Date: 2025-04-14
## Author: Mengge LYU
## Version: 1.0
##################################################

# library loading-----------------------------------------------------------------
library(ggplot2)
library(PrInCE)
library(CCprofiler)
library(dplyr)
library(utils)
library(rio)
library(data.table)
library(corrplot)
library(tidyverse)


# data loading ------------------------------------------------------------
test1 <- read.delim("E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\8_protmattrix\\TPC1_protmatrix_20240918.tsv", row.names=1, quote="")
test2 <- read.delim("E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\8_protmattrix\\TPC2_protmatrix_20240918.tsv", row.names=1, quote="")
test3 <- read.delim("E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\8_protmattrix\\TPC3_protmatrix_20240918.tsv", row.names=1, quote="")

# missing value plot test1------------------------------------------------------
na_counts1 <- data.frame(colSums(is.na(test1)))
rownames(na_counts1) <- gsub("TPC.1_NO1_30minDIA_", "", rownames(na_counts1))
df_sorted1 <- na_counts1 %>%
  tibble::rownames_to_column("row_names") %>%
  arrange(as.numeric(row_names)) %>%
  tibble::column_to_rownames("row_names")
colnames(df_sorted1)[1] <- "na_count"
ggplot(df_sorted1, aes(
  x = as.numeric(rownames(df_sorted1)),  
  y = na_count)) +
  geom_point(size = 3, color = "steelblue") +
  labs(title = "The NA number of each fraction", x = "Frac Num", y = "NA Num") +
  scale_x_continuous(breaks = as.numeric(rownames(df_sorted1)),  
                     labels = rownames(df_sorted1)) +  
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# missing value plot test2------------------------------------------------------
na_counts2 <- data.frame(colSums(is.na(test2)))
rownames(na_counts2) <- gsub("TPC.1_NO2_30minDIA_", "", rownames(na_counts2))
df_sorted2 <- na_counts2 %>%
  tibble::rownames_to_column("row_names") %>%
  arrange(as.numeric(row_names)) %>%
  tibble::column_to_rownames("row_names")
colnames(df_sorted2)[1] <- "na_count"
ggplot(df_sorted2, aes(
  x = as.numeric(rownames(df_sorted2)),  
  y = na_count)) +
  geom_point(size = 3, color = "steelblue") +
  labs(title = "The NA number of each fraction", x = "Frac Num", y = "NA Num") +
  scale_x_continuous(breaks = as.numeric(rownames(df_sorted2)),  
                     labels = rownames(df_sorted2)) +  
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# missing value plot test3------------------------------------------------------
na_counts3 <- data.frame(colSums(is.na(test3)))
rownames(na_counts3) <- gsub("TPC.1_NO3_30minDIA_", "", rownames(na_counts3))
df_sorted3 <- na_counts3 %>%
  tibble::rownames_to_column("row_names") %>%
  arrange(as.numeric(row_names)) %>%
  tibble::column_to_rownames("row_names")
colnames(df_sorted3)[1] <- "na_count"
ggplot(df_sorted3, aes(
  x = as.numeric(rownames(df_sorted3)),  
  y = na_count)) +
  geom_point(size = 3, color = "steelblue") +
  labs(title = "The NA number of each fraction", x = "Frac Num", y = "NA Num") +
  scale_x_continuous(breaks = as.numeric(rownames(df_sorted3)),  
                     labels = rownames(df_sorted3)) +  
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# clean top3 frac ---------------------------------------------------------
df1 <- test1[,-grep("TPC.1_NO3_30minDIA_1", colnames(test1))]
# na_proportion <- rowMeans(is.na(test1))
# df_clean <- test1[na_proportion <= 0.9, ]

# precursor matrix ----------------------------------------------------------
df <- fread("E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\report.tsv")
df_pick <- df[df$Q.Value <= 0.05 & df$Global.Q.Value <= 0.05 & df$PG.Q.Value <= 0.05 & df$Global.PG.Q.Value <= 0.05 & df$Protein.Q.Value <= 0.05,]
df_pick1 <- df_pick[,c("File.Name", "Protein.Group", "Protein.Ids","Stripped.Sequence","PG.Quantity")]
df_pick2 <- unique(df_pick1) # from 11934754 to 10669425 rows
colnames(df_pick2)
df_distinct <- df_pick2 %>% 
  distinct(File.Name,Protein.Group,Protein.Ids,Stripped.Sequence,PG.Quantity, .keep_all = T) 
JAKA_qc <- df_distinct[grep("POOL", df_distinct$File.Name),]
rm(df)
rm(df_pick)
rm(df_pick1)
rm(df_pick2)
# long to wide ------------------------------------------------------------
JAKA_qc2 <- dcast(JAKA_qc, Protein.Group + Stripped.Sequence ~ File.Name, value.var = "PG.Quantity", fun.aggregate = max)
is.na(JAKA_qc2)<-sapply(JAKA_qc2, is.infinite)
colnames(JAKA_qc2) <- gsub("(.*)lyumg_thyroid_complex_","",colnames(JAKA_qc2))
colnames(JAKA_qc2) <- gsub(".raw","",colnames(JAKA_qc2))
# peptide matrix ----------------------------------------------------------
peptide_qc <- JAKA_qc2[,2:14]
peptide_qc2 <- peptide_qc[rowSums(is.na(peptide_qc[, 2:13])) != 12, ]
dim(peptide_qc2) # [1] 30593    13
length(unique(peptide_qc2$Stripped.Sequence))
# [1] 30593
non_na_counts <- sapply(peptide_qc2, function(x) sum(!is.na(x)))
non_na_counts_df <- data.frame(
  Column = names(non_na_counts),
  Count = non_na_counts)
pep_draw <- non_na_counts_df[2:13,]
pdf("Cell_line_replicate_MGL_20250424.pdf", width = 10, height = 6)
ggplot(pep_draw, aes(x = Column, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +  
  xlab("Filename") +
  ylab("Number of Peptide") +
  theme_classic() +  
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),  
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),  
    axis.text.y = element_text(size = 12, color = "black") 
  )
dev.off()
save.image("JAKA_QC_MGL_20250424.Rdata")
# protein matrix ----------------------------------------------------------
###pool sample
thy <- rio::import("E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrix_Thyroid_cellline_all_FDR005_20240906.tsv")
thy_pool <- thy[,grep("POOL", colnames(thy))]
thy_pool2 <- cbind(thy[,1:2], thy_pool)
dim(thy_pool2) # 49588 14
thy_pool2_Dall <- thy_pool2 %>%
  filter(!(rowSums(is.na(thy_pool2[, 3:ncol(thy_pool2)])) == (ncol(thy_pool2) - 2))) # 30593 peptides
thy_pool2_prot <- aggregate(thy_pool2_Dall,
                            by = list(thy_pool2_Dall$protein_id),
                            FUN = max,
                            na.rm = TRUE)
is.na(thy_pool2_prot)<-sapply(thy_pool2_prot, is.infinite) # 4300 prots
thy_pool2_pep <- thy_pool2_Dall[,2:length(colnames(thy_pool2_Dall))]
thy_pool2_pep2 <- thy_pool2_pep[,2:length(colnames(thy_pool2_pep))]
rownames(thy_pool2_pep2) <- thy_pool2_pep[,1]
thy_pep_pear <- cor(thy_pool2_pep2, method = c("pearson"), use = "pairwise.complete.obs")
thy_pep_pear_nor <- round(thy_pep_pear, 2)
res_pear <- rcorr(as.matrix(thy_pool2_pep2), type = "pearson")
res_spear <- rcorr(as.matrix(thy_pool2_pep2), type = "spearman")
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
flattenCorrMatrix(res_pear$r, res_pear$P)
corrplot(thy_pep_pear_nor, add=TRUE, type="lower", method="number",order="AOE", col="black",diag=FALSE,tl.pos="n", cl.pos="n")
thy_pep_spear <- cor(thy_pool2_pep2, method = c("spearman"), use = "pairwise.complete.obs")
thy_pep_spear_nor <- round(thy_pep_spear, 2)
corrplot(thy_pep_spear_nor, add=TRUE, type="lower", method="number",order="AOE", col="black",diag=FALSE,tl.pos="n", cl.pos="n")
# library("PerformanceAnalytics")
# chart.Correlation(thy_pool2_pep2, histogram=TRUE, method = c("pearson"),pch=19)
# chart.Correlation(thy_pool2_pep2, histogram=TRUE, method = c("spearman"),pch=19)

# prot correlation
thy_pool2_prot2 <- thy_pool2_prot[,4:length(colnames(thy_pool2_prot))]
rownames(thy_pool2_prot2) <- thy_pool2_prot[,2]
thy_prot_pear <- cor(thy_pool2_prot2, method = c("pearson"), use = "pairwise.complete.obs")
thy_prot_spear <- cor(thy_pool2_prot2, method = c("spearman"), use = "pairwise.complete.obs")
thy_prot_pear_nor <- round(thy_prot_pear, 2)
thy_prot_spear_nor <- round(thy_prot_spear, 2)
corrplot(thy_prot_pear_nor, add=TRUE, type="lower", method="number",order="AOE", col="black",diag=FALSE,tl.pos="n", cl.pos="n")

# replicate protein matrix ------------------------------------------------
JAKA_all <- df_distinct[grep("all", df_distinct$File.Name),]
JAKA_all2 <- dcast(JAKA_all, Protein.Group  ~ File.Name, value.var = "PG.Quantity", fun.aggregate = max)
is.na(JAKA_all2)<-sapply(JAKA_all2, is.infinite)
colnames(JAKA_all2) <- gsub("(.*)lyumg_thyroid_complex_","",colnames(JAKA_all2))
colnames(JAKA_all2) <- gsub(".raw","",colnames(JAKA_all2))
expr <- JAKA_all2 %>% 
  column_to_rownames(var = names(.)[1]) %>%   
  as.matrix()
expr <- expr[,-grep("FTC133", colnames(expr))]
cor_matrix <- cor(
  expr,                 
  method = "pearson",              
  use = "pairwise.complete.obs")
corrplot(
  cor_matrix,
  method = "circle",                
  type = "upper",                  
  order = "hclust",                
  tl.col = "black",                
  tl.cex = 0.6,                    
  addCoef.col = "white",           
  number.cex = 0.5,                
  diag = FALSE)
# TPC-1 pep-------------------------------------------------------------------
library(VennDiagram)
library(data.table)
tpc_columns <- grep("TPC", colnames(peptide_qc2))
selected_columns <- c(1, tpc_columns)
# selected_columns
# [1] 1 2 3 4
tpc_df <- peptide_qc2[,1:4]
set1 <- tpc_df$Stripped.Sequence[which(!is.na(tpc_df[, 2]))]
set2 <- tpc_df$Stripped.Sequence[which(!is.na(tpc_df[, 3]))]
set3 <- tpc_df$Stripped.Sequence[which(!is.na(tpc_df[, 4]))]

set1 <- as.character(set1)
set2 <- as.character(set2)
set3 <- as.character(set3)

area1 <- length(set1)
area2 <- length(set2)
area3 <- length(set3)
n12 <- length(intersect(set1, set2))
n23 <- length(intersect(set2, set3))
n13 <- length(intersect(set1, set3))
n123 <- length(intersect(intersect(set1, set2), set3))
pdf("TPC_pep_venn_plot_20250424.pdf", width = 10, height = 10)
venn_plot <- draw.triple.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  n12 = n12,
  n23 = n23,
  n13 = n13,
  n123 = n123,
  category = c("TPC-1_NO1_30minDIA_POOL", "TPC-1_NO2_30minDIA_POOL", "TPC-1_NO3_30minDIA_POOL"),
  fill = c("skyblue", "pink", "lightgreen"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.5
)

grid.draw(venn_plot)
dev.off()
save.image("JAKA_QC_MGL_20250424.RData")
# prot matrix - final ----------------------------------------------------------
prot_qc <- aggregate(JAKA_qc2,
                     by = list(JAKA_qc2$Protein.Group),
                     FUN = max,
                     na.rm = TRUE)
prot_qc2 <- prot_qc[,c(2,4:15)]
is.na(prot_qc2)<-sapply(prot_qc2, is.infinite)
dim(prot_qc2) # [1] 4300    13
length(unique(prot_qc2$Protein.Group))
# [1] 4300
non_na_counts <- sapply(prot_qc2, function(x) sum(!is.na(x)))
non_na_counts_df <- data.frame(
  Column = names(non_na_counts),
  Count = non_na_counts
)
prot_draw <- non_na_counts_df[2:13,]
# pdf("Cell_line_prot_replicate_MGL_20250425.pdf", width = 10, height = 6)
# ggplot(prot_draw, aes(x = Column, y = Count)) +
#   geom_bar(stat = "identity", fill = "skyblue", color = "black") +  
#   xlab("Filename") +
#   ylab("Number of Peptide") +
#   theme_classic() +  
#   theme(
#     plot.title = element_text(size = 16, face = "bold", color = "black"),  
#     axis.title.x = element_text(size = 14, face = "bold", color = "black"),  
#     axis.title.y = element_text(size = 14, face = "bold", color = "black"),  
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),  
#     axis.text.y = element_text(size = 12, color = "black") 
#   )
# dev.off()
rownames(prot_draw) <- gsub("_30minDIA_POOL","",rownames(prot_draw))
prot_draw$Column <- gsub("_30minDIA_POOL","",prot_draw$Column)
prot_draw <- prot_draw[-c(7,8,9),]
pdf("Three_cell_line_prot_replicate_MGL_20250425.pdf", width = 12, height = 6)
ggplot(prot_draw, aes(x = Column, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.7) +  
  geom_text(
    aes(label = Count), 
    vjust = -0.5,            
    color = "black", 
    size = 6.5,
    fontface = "bold"
  ) +
  xlab("Filename") +
  ylab("Number of Protein") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme_classic() +  
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),  
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),  
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),  
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1, 
      size = 12, 
      color = "black",
      margin = margin(t = 5)  
    ),  
    axis.text.y = element_text(size = 12, color = "black"),
    panel.grid.major.y = element_line(color = "gray90")  
  )
dev.off()

# TPC prot Venn --------------------------------------------------------------------
grep("TPC", colnames(prot_qc2))

# selected_columns
# [1] 1 2 3 4
tpc_df <- prot_qc2[,1:4]
colnames(tpc_df)[2:4] <- gsub("_30minDIA_POOL","",colnames(tpc_df)[2:4] )
set1 <- tpc_df$Protein.Group[which(!is.na(tpc_df[, 2]))]
set2 <- tpc_df$Protein.Group[which(!is.na(tpc_df[, 3]))]
set3 <- tpc_df$Protein.Group[which(!is.na(tpc_df[, 4]))]

set1 <- as.character(set1)
set2 <- as.character(set2)
set3 <- as.character(set3)

area1 <- length(set1)
area2 <- length(set2)
area3 <- length(set3)
n12 <- length(intersect(set1, set2))
n23 <- length(intersect(set2, set3))
n13 <- length(intersect(set1, set3))
n123 <- length(intersect(intersect(set1, set2), set3))
pdf("TPC_prot_venn_plot_20250425.pdf", width = 10, height = 10)
venn_plot <- draw.triple.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  n12 = n12,
  n23 = n23,
  n13 = n13,
  n123 = n123,
  category = c("TPC-1_NO1", "TPC-1_NO2", "TPC-1_NO3"),
  fill = c("skyblue", "pink", "lightgreen"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.5
)

grid.draw(venn_plot)
dev.off()

# FTC prot Venn --------------------------------------------------------------------
grep("FTC", colnames(prot_qc2))

# selected_columns
# [1] 1 2 3 4
tpc_df <- prot_qc2[,c(1,11:13)]
colnames(tpc_df)[2:4] <- gsub("_30minDIA_POOL","",colnames(tpc_df)[2:4] )
set1 <- tpc_df$Protein.Group[which(!is.na(tpc_df[, 2]))]
set2 <- tpc_df$Protein.Group[which(!is.na(tpc_df[, 3]))]
set3 <- tpc_df$Protein.Group[which(!is.na(tpc_df[, 4]))]

set1 <- as.character(set1)
set2 <- as.character(set2)
set3 <- as.character(set3)

area1 <- length(set1)
area2 <- length(set2)
area3 <- length(set3)
n12 <- length(intersect(set1, set2))
n23 <- length(intersect(set2, set3))
n13 <- length(intersect(set1, set3))
n123 <- length(intersect(intersect(set1, set2), set3))
pdf("FTC238_prot_venn_plot_20250425.pdf", width = 10, height = 10)
venn_plot <- draw.triple.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  n12 = n12,
  n23 = n23,
  n13 = n13,
  n123 = n123,
  category = c("FTC238_NO1", "FTC238_NO2", "FTC238_NO3"),
  fill = c("skyblue", "pink", "lightgreen"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.5
)

grid.draw(venn_plot)
dev.off()
# Nthy prot Venn --------------------------------------------------------------------
grep("Nthy", colnames(prot_qc2))

# selected_columns
# [1] 1 2 3 4
tpc_df <- prot_qc2[,c(1,5:7)]
colnames(tpc_df)[2:4] <- gsub("_30minDIA_POOL","",colnames(tpc_df)[2:4] )
set1 <- tpc_df$Protein.Group[which(!is.na(tpc_df[, 2]))]
set2 <- tpc_df$Protein.Group[which(!is.na(tpc_df[, 3]))]
set3 <- tpc_df$Protein.Group[which(!is.na(tpc_df[, 4]))]

set1 <- as.character(set1)
set2 <- as.character(set2)
set3 <- as.character(set3)

area1 <- length(set1)
area2 <- length(set2)
area3 <- length(set3)
n12 <- length(intersect(set1, set2))
n23 <- length(intersect(set2, set3))
n13 <- length(intersect(set1, set3))
n123 <- length(intersect(intersect(set1, set2), set3))
pdf("Nthy_prot_venn_plot_20250425.pdf", width = 10, height = 10)
venn_plot <- draw.triple.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  n12 = n12,
  n23 = n23,
  n13 = n13,
  n123 = n123,
  category = c("Nthy_NO1", "Nthy_NO2", "Nthy_NO3"),
  fill = c("skyblue", "pink", "lightgreen"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.5
)

grid.draw(venn_plot)
dev.off()


# Nthy Jaccard test-----------------------------------------------------------------
nthy_binary <- tpc_df

nthy_binary[, -1] <- ifelse(is.na(tpc_df[, -1]), 0, 1)  
nthy_binary <- nthy_binary[rowSums(nthy_binary[, -1]) > 0, ]  

jaccard <- function(a, b) {
  intersection <- sum(a & b)
  union <- sum(a | b)
  if (union == 0) return(NA)  
  intersection / union
}


nthy1 <- nthy_binary$Nthy_NO1
nthy2 <- nthy_binary$Nthy_NO2
nthy3 <- nthy_binary$Nthy_NO3


jaccard_matrix <- matrix(nrow = 3, ncol = 3, dimnames = list(
  c("Nthy_NO1", "Nthy_NO2", "Nthy_NO3"),
  c("Nthy_NO1", "Nthy_NO2", "Nthy_NO3")
))

jaccard_matrix["Nthy_NO1", "Nthy_NO2"] <- jaccard(nthy1, nthy2)
jaccard_matrix["Nthy_NO1", "Nthy_NO3"] <- jaccard(nthy1, nthy3)
jaccard_matrix["Nthy_NO2", "Nthy_NO3"] <- jaccard(nthy2, nthy3)


jaccard_matrix[lower.tri(jaccard_matrix)] <- t(jaccard_matrix)[lower.tri(jaccard_matrix)]
diag(jaccard_matrix) <- 1  

print(jaccard_matrix)

tpc_df$mean_abundance <- rowMeans(tpc_df[, -1], na.rm = TRUE)

shared_proteins <- nthy_binary[rowSums(nthy_binary[, -1]) == 3, ]
summary(tpc_df[tpc_df$Protein.Group %in% shared_proteins$Protein.Group, "mean_abundance"])

# TPC Jaccard test-----------------------------------------------------------------
tpc_df <- prot_qc2[,1:4]
tpc_binary <- tpc_df

tpc_binary[, -1] <- ifelse(is.na(tpc_df[, -1]), 0, 1)  
tpc_binary <- tpc_binary[rowSums(tpc_binary[, -1]) > 0, ]  

jaccard <- function(a, b) {
  intersection <- sum(a & b)
  union <- sum(a | b)
  if (union == 0) return(NA) 
  intersection / union
}


tpc1 <- tpc_binary$`TPC-1_NO1_30minDIA_POOL`
tpc2 <- tpc_binary$`TPC-1_NO2_30minDIA_POOL`
tpc3 <- tpc_binary$`TPC-1_NO3_30minDIA_POOL`


jaccard_matrix <- matrix(nrow = 3, ncol = 3, dimnames = list(
  c("TPC_NO1", "TPC_NO2", "TPC_NO3"),
  c("TPC_NO1", "TPC_NO2", "TPC_NO3")
))

jaccard_matrix["TPC_NO1", "TPC_NO2"] <- jaccard(tpc1, tpc2)
jaccard_matrix["TPC_NO1", "TPC_NO3"] <- jaccard(tpc1, tpc3)
jaccard_matrix["TPC_NO2", "TPC_NO3"] <- jaccard(tpc2, tpc3)


jaccard_matrix[lower.tri(jaccard_matrix)] <- t(jaccard_matrix)[lower.tri(jaccard_matrix)]
diag(jaccard_matrix) <- 1  

print(jaccard_matrix)

# FTC238 Jaccard test-----------------------------------------------------------------
tpc_df <- prot_qc2[,c(1,11:13)]
tpc_binary <- tpc_df

tpc_binary[, -1] <- ifelse(is.na(tpc_df[, -1]), 0, 1)  
tpc_binary <- tpc_binary[rowSums(tpc_binary[, -1]) > 0, ]  

jaccard <- function(a, b) {
  intersection <- sum(a & b)
  union <- sum(a | b)
  if (union == 0) return(NA)  
  intersection / union
}


tpc1 <- tpc_binary$FTC238_NO1_30minDIA_POOL
tpc2 <- tpc_binary$FTC238_NO2_30minDIA_POOL
tpc3 <- tpc_binary$FTC238_NO3_30minDIA_POOL


jaccard_matrix <- matrix(nrow = 3, ncol = 3, dimnames = list(
  c("FTC238_NO1", "FTC238_NO2", "FTC238_NO3"),
  c("FTC238_NO1", "FTC238_NO2", "FTC238_NO3")
))

jaccard_matrix["FTC238_NO1", "FTC238_NO2"] <- jaccard(tpc1, tpc2)
jaccard_matrix["FTC238_NO1", "FTC238_NO3"] <- jaccard(tpc1, tpc3)
jaccard_matrix["FTC238_NO2", "FTC238_NO3"] <- jaccard(tpc2, tpc3)


jaccard_matrix[lower.tri(jaccard_matrix)] <- t(jaccard_matrix)[lower.tri(jaccard_matrix)]
diag(jaccard_matrix) <- 1  

# Jaccard plot ------------------------------------------------------------
print(jaccard_matrix)
pdf("Jaccard_plot_MGL_20250427.pdf", width = 10, height = 10)
corrplot(jaccard_matrix, 
         method = "color",
         addCoef.col = "black",
         tl.col = "darkblue",
         title = "Jaccard Similarity Between Nthy Replicates")
dev.off()
save.image("JAKA_QC_MGL_20250424.Rdata")

df <- fread("Z:\\members\\lyumengge\\00_NC_rebuttal\\Automatic_part\\JAKA_qc_16files_MGL_20250422\\report.tsv")
df_pick <- df[df$Q.Value <= 0.05 & df$Global.Q.Value <= 0.05 & df$PG.Q.Value <= 0.05 & df$Global.PG.Q.Value <= 0.05 & df$Protein.Q.Value <= 0.05,]
df_pick1 <- df_pick[,c("File.Name", "Protein.Group", "Protein.Ids","Stripped.Sequence","PG.Quantity")]
df_pick2 <- unique(df_pick1) # from 11934754 to 10669425 rows
colnames(df_pick2)
df_distinct <- df_pick2 %>% 
  distinct(File.Name,Protein.Group,Protein.Ids,Stripped.Sequence,PG.Quantity, .keep_all = T) 
df_wide <- dcast(df_distinct, Protein.Group + Stripped.Sequence ~ File.Name, value.var = "PG.Quantity", fun.aggregate = max)
prot_qc <- aggregate(df_wide,
                     by = list(df_wide$Protein.Group),
                     FUN = max,
                     na.rm = TRUE)
prot_qc2 <- prot_qc[,c(2,4:19)]
is.na(prot_qc2)<-sapply(prot_qc2, is.infinite)
dim(prot_qc2) # [1] 3567    17
length(unique(prot_qc2$Protein.Group))
# [1] 3567
colnames(prot_qc2) <- gsub("(.*)B20250421lyumg_thyroid_complex_TPC","",colnames(prot_qc2))
colnames(prot_qc2) <- gsub("_30minDIA_","",colnames(prot_qc2))
colnames(prot_qc2) <- gsub(".raw","",colnames(prot_qc2))

cv_data <- prot_qc2 %>%
  select(Protein.Group, starts_with("auto"), starts_with("hand")) %>%
  rowwise() %>%
  mutate(
    auto_CV = sd(c_across(starts_with("auto")), na.rm = TRUE) / 
      mean(c_across(starts_with("auto")), na.rm = TRUE) * 100,
    hand_CV = sd(c_across(starts_with("hand")), na.rm = TRUE) / 
      mean(c_across(starts_with("hand")), na.rm = TRUE) * 100
  ) %>%
  select(Protein.Group, auto_CV, hand_CV) %>%
  melt(id.vars = "Protein.Group", 
       variable.name = "Group", 
       value.name = "CV") %>%
  mutate(Group = ifelse(Group == "auto_CV", "Automated", "Manual"))

ggplot(cv_data, aes(x = Group, y = CV, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 18, size = 4, color = "red") +
  labs(title = "Coefficient of Variation (CV) of Protein Abundance",
       x = "Processing Method",
       y = "CV (%)",
       caption = "Red diamond indicates median value") +
  scale_fill_manual(values = c("Automated" = "#4E79A7", "Manual" = "#F28E2B")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

cv_summary <- cv_data %>%
  group_by(Group) %>%
  summarise(
    Median_CV = median(CV, na.rm = TRUE),
    IQR_CV = IQR(CV, na.rm = TRUE),
    .groups = 'drop')


non_na_counts <- sapply(prot_qc2, function(x) sum(!is.na(x)))
non_na_counts_df <- data.frame(
  Column = names(non_na_counts),
  Count = non_na_counts
)
prot_draw <- non_na_counts_df[2:17,]
rownames(prot_draw) <- gsub("(.*)B20250421lyumg_thyroid_complex_TPC","",rownames(prot_draw))
rownames(prot_draw) <- gsub("_30minDIA_","",rownames(prot_draw))
rownames(prot_draw) <- gsub(".raw","",rownames(prot_draw))
prot_draw$Column <- gsub("(.*)B20250421lyumg_thyroid_complex_TPC","",prot_draw$Column)
prot_draw$Column <- gsub("_30minDIA_","",prot_draw$Column)
prot_draw$Column <- gsub(".raw","",prot_draw$Column)

prot_draw$Group <- rep(c("hand", "auto"), each = 8)  
group_summary <- prot_draw %>%
  group_by(Group) %>%
  summarize(Mean = mean(Count),
            SD = sd(Count))

pdf("AutovsHand_prot_replicate_MGL_20250505.pdf", width = 12, height = 6)
ggplot(prot_draw, aes(x = Column, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.7) +  
  geom_text(
    aes(label = Count), 
    vjust = -0.5,            
    color = "black", 
    size = 6.5,
    fontface = "bold"
  ) +
  labs(
    x = "Filename",
    y = "Number of Proteins",  
    title = "Your Title Here"  
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme_classic() +  
  theme(  
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),  
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),  
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),  
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1, 
      size = 12, 
      color = "black",
      margin = margin(t = 5)
    ),  
    axis.text.y = element_text(size = 12, color = "black"),
    panel.grid.major.y = element_line(color = "gray90")
  )
dev.off()

pdf("AutovsHand_cv_MGL_20250505.pdf", width = 12, height = 6)
# ggplot(group_summary, aes(x = Group, y = CV, fill = Group)) +
#   geom_bar(stat = "identity") +  
#   labs(title = "Coefficient of Variation by Group", x = "Group", y = "CV (%)") +  
#   theme_minimal()  
# dev.off()

ggplot(prot_draw, aes(x = Group, y = Count, fill = Group)) +
  geom_violin(trim = FALSE, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.7) +  
  geom_jitter(width = 0.1, size = 2, alpha = 0.5) +       
  scale_fill_manual(values = c("auto" = "#FF6B6B", "hand" = "#4ECDC4")) +  
  labs(title = "Distribution of Protein Counts by Group",
       x = "Preparation Method",
       y = "Protein Count") +
  theme_minimal() +
  theme(legend.position = "none")  

result <- prot_draw %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Count),
    SD = sd(Count),
    CV = (SD / Mean) * 100,  
    Q1 = quantile(Count, 0.25),  
    Median = median(Count), 
    Q3 = quantile(Count, 0.75),  
    .groups = 'drop'
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))  


