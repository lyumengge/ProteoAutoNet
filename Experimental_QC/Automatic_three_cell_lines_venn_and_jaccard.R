##################################################
## Project: NC_rebuttal
## Script purpose:Three cell lines venn and jaccard
## Date: 2025-04-29
## Author: Mengge LYU
## Version: 1.0
##################################################

# library loading-----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(utils)
library(rio)
library(data.table)
library(corrplot)

# data loading ------------------------------------------------------------
setwd("E:\\Thyroid_cellline_all")
df <- fread("report.tsv")
df_pick <- df[df$Q.Value <= 0.05 & df$Global.Q.Value <= 0.05 & df$PG.Q.Value <= 0.05 & df$Global.PG.Q.Value <= 0.05 & df$Protein.Q.Value <= 0.05,]
df_pick1 <- df_pick[,c("File.Name", "Protein.Group", "Protein.Ids","Stripped.Sequence","PG.Quantity")]
df_pick2 <- unique(df_pick1) # from 11934754 to 10669425 rows
colnames(df_pick2)
df_distinct <- df_pick2 %>% 
  distinct(File.Name,Protein.Group,Protein.Ids,Stripped.Sequence,PG.Quantity, .keep_all = T) 
JAKA_qc <- df_distinct[grep("POOL", df_distinct$File.Name),]
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
# protein matrix ----------------------------------------------------------
prot_qc <- aggregate(JAKA_qc2,
                     by = list(JAKA_qc2$Protein.Group),
                     FUN = max,
                     na.rm = TRUE)
prot_qc2 <- prot_qc[,c(2,4:15)]
is.na(prot_qc2)<-sapply(prot_qc2, is.infinite)
dim(prot_qc2) # [1] 4300    13
length(unique(prot_qc2$Protein.Group))
# [1] 4300
# protein identification plot ---------------------------------------------
non_na_counts <- sapply(prot_qc2, function(x) sum(!is.na(x)))
non_na_counts_df <- data.frame(
  Column = names(non_na_counts),
  Count = non_na_counts)
prot_draw <- non_na_counts_df[2:13,]
rownames(prot_draw) <- gsub("_30minDIA_POOL","",rownames(prot_draw))
prot_draw$Column <- gsub("_30minDIA_POOL","",prot_draw$Column)
prot_draw <- prot_draw[-c(7,8,9),]

ggplot(prot_draw, aes(x = Column, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.7) +  
  geom_text(aes(label = Count), 
    vjust = -0.5,            
    color = "black", 
    size = 6.5,
    fontface = "bold") +
  xlab("Filename") +
  ylab("Number of Protein") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme_classic() +  
  theme(plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),  
        axis.title.x = element_text(size = 14, face = "bold", color = "black"),  
        axis.title.y = element_text(size = 14, face = "bold", color = "black"),  
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12,color = "black",margin = margin(t = 5)),  
        axis.text.y = element_text(size = 12, color = "black"),
        panel.grid.major.y = element_line(color = "gray90"))

# Venn plot ---------------------------------------------------------------
## TPC prot Venn --------------------------------------------------------------------
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
venn_plot <- draw.triple.venn(area1 = area1,
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
                              cat.cex = 1.5)
grid.draw(venn_plot)
## FTC prot Venn --------------------------------------------------------------------
grep("FTC", colnames(prot_qc2))
ftc_df <- prot_qc2[,c(1,11:13)]
colnames(ftc_df)[2:4] <- gsub("_30minDIA_POOL","",colnames(ftc_df)[2:4] )
set1 <- ftc_df$Protein.Group[which(!is.na(ftc_df[, 2]))]
set2 <- ftc_df$Protein.Group[which(!is.na(ftc_df[, 3]))]
set3 <- ftc_df$Protein.Group[which(!is.na(ftc_df[, 4]))]
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
venn_plot <- draw.triple.venn(area1 = area1,
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
                              cat.cex = 1.5)
grid.draw(venn_plot)
## Nthy prot Venn --------------------------------------------------------------------
grep("Nthy", colnames(prot_qc2))
nthy_df <- prot_qc2[,c(1,5:7)]
colnames(nthy_df)[2:4] <- gsub("_30minDIA_POOL","",colnames(nthy_df)[2:4] )
set1 <- nthy_df$Protein.Group[which(!is.na(nthy_df[, 2]))]
set2 <- nthy_df$Protein.Group[which(!is.na(nthy_df[, 3]))]
set3 <- nthy_df$Protein.Group[which(!is.na(nthy_df[, 4]))]
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
venn_plot <- draw.triple.venn(area1 = area1,area2 = area2,area3 = area3,
                              n12 = n12,n23 = n23,n13 = n13,n123 = n123,
                              category = c("Nthy_NO1", "Nthy_NO2", "Nthy_NO3"),
                              fill = c("skyblue", "pink", "lightgreen"),
                              lty = "blank",cex = 1.5,cat.cex = 1.5)
grid.draw(venn_plot)
# Jaccard ---------------------------------------------------------------
## TPC Jaccard -------------------------------------------------------------
tpc_df <- prot_qc2[,1:4]
tpc_binary <- tpc_df
tpc_binary[, -1] <- ifelse(is.na(tpc_df[, -1]), 0, 1) 
tpc_binary <- tpc_binary[rowSums(tpc_binary[, -1]) > 0, ]  
jaccard <- function(a, b) {
  intersection <- sum(a & b)
  union <- sum(a | b)
  if (union == 0) return(NA)  
  intersection / union}

tpc1 <- tpc_binary$`TPC-1_NO1_30minDIA_POOL`
tpc2 <- tpc_binary$`TPC-1_NO2_30minDIA_POOL`
tpc3 <- tpc_binary$`TPC-1_NO3_30minDIA_POOL`

jaccard_matrix <- matrix(nrow = 3, ncol = 3, dimnames = list(
  c("TPC_NO1", "TPC_NO2", "TPC_NO3"),
  c("TPC_NO1", "TPC_NO2", "TPC_NO3")))

jaccard_matrix["TPC_NO1", "TPC_NO2"] <- jaccard(tpc1, tpc2)
jaccard_matrix["TPC_NO1", "TPC_NO3"] <- jaccard(tpc1, tpc3)
jaccard_matrix["TPC_NO2", "TPC_NO3"] <- jaccard(tpc2, tpc3)

jaccard_matrix[lower.tri(jaccard_matrix)] <- t(jaccard_matrix)[lower.tri(jaccard_matrix)]
diag(jaccard_matrix) <- 1  
print(jaccard_matrix)
# TPC_NO1   TPC_NO2   TPC_NO3
# TPC_NO1 1.0000000 0.8570210 0.8626561
# TPC_NO2 0.8570210 1.0000000 0.8585774
# TPC_NO3 0.8626561 0.8585774 1.0000000

## Nthy Jaccard -------------------------------------------------------------
nthy_binary <- nthy_df
nthy_binary[, -1] <- ifelse(is.na(nthy_df[, -1]), 0, 1)  
nthy_binary <- nthy_binary[rowSums(nthy_binary[, -1]) > 0, ]  

jaccard <- function(a, b) {
  intersection <- sum(a & b)
  union <- sum(a | b)
  if (union == 0) return(NA) 
  intersection / union}

nthy1 <- nthy_binary$Nthy_NO1
nthy2 <- nthy_binary$Nthy_NO2
nthy3 <- nthy_binary$Nthy_NO3

jaccard_matrix <- matrix(nrow = 3, ncol = 3, dimnames = list(
  c("Nthy_NO1", "Nthy_NO2", "Nthy_NO3"),
  c("Nthy_NO1", "Nthy_NO2", "Nthy_NO3")))

jaccard_matrix["Nthy_NO1", "Nthy_NO2"] <- jaccard(nthy1, nthy2)
jaccard_matrix["Nthy_NO1", "Nthy_NO3"] <- jaccard(nthy1, nthy3)
jaccard_matrix["Nthy_NO2", "Nthy_NO3"] <- jaccard(nthy2, nthy3)

jaccard_matrix[lower.tri(jaccard_matrix)] <- t(jaccard_matrix)[lower.tri(jaccard_matrix)]
diag(jaccard_matrix) <- 1  

print(jaccard_matrix)
# Nthy_NO1  Nthy_NO2  Nthy_NO3
# Nthy_NO1 1.0000000 0.8528260 0.8542083
# Nthy_NO2 0.8528260 1.0000000 0.8563678
# Nthy_NO3 0.8542083 0.8563678 1.0000000
## FTC238 Jaccard -------------------------------------------------------------
ftc_df <- prot_qc2[,c(1,11:13)]
ftc_binary <- ftc_df

ftc_binary[, -1] <- ifelse(is.na(ftc_df[, -1]), 0, 1)  
ftc_binary <- ftc_binary[rowSums(ftc_binary[, -1]) > 0, ]  

jaccard <- function(a, b) {
  intersection <- sum(a & b)
  union <- sum(a | b)
  if (union == 0) return(NA) 
  intersection / union}

ftc1 <- ftc_binary$FTC238_NO1_30minDIA_POOL
ftc2 <- ftc_binary$FTC238_NO2_30minDIA_POOL
ftc3 <- ftc_binary$FTC238_NO3_30minDIA_POOL

jaccard_matrix <- matrix(nrow = 3, ncol = 3, dimnames = list(
  c("FTC238_NO1", "FTC238_NO2", "FTC238_NO3"),
  c("FTC238_NO1", "FTC238_NO2", "FTC238_NO3")))

jaccard_matrix["FTC238_NO1", "FTC238_NO2"] <- jaccard(ftc1, ftc2)
jaccard_matrix["FTC238_NO1", "FTC238_NO3"] <- jaccard(ftc1, ftc3)
jaccard_matrix["FTC238_NO2", "FTC238_NO3"] <- jaccard(ftc2, ftc3)

jaccard_matrix[lower.tri(jaccard_matrix)] <- t(jaccard_matrix)[lower.tri(jaccard_matrix)]
diag(jaccard_matrix) <- 1  

print(jaccard_matrix)
# FTC238_NO1 FTC238_NO2 FTC238_NO3
# FTC238_NO1  1.0000000  0.8630740  0.8416735
# FTC238_NO2  0.8630740  1.0000000  0.8522073
# FTC238_NO3  0.8416735  0.8522073  1.0000000
