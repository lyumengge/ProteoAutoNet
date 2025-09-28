# Date: 2024-09-06
# Project: Thyroid cell line - Nthy protein matrix 
# AUTHOR: Mengge LYU

library(data.table)
library(CCprofiler)
library(dplyr)
library(ggplot2)
library(tidyr)

# FDR005
df <- fread("E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\report.tsv")
df_pick <- df[df$Q.Value <= 0.05 & df$Global.Q.Value <= 0.05 & df$PG.Q.Value <= 0.05 & df$Global.PG.Q.Value <= 0.05 & df$Protein.Q.Value <= 0.05,]
df_pick1 <- df_pick[,c("File.Name", "Protein.Group", "Protein.Ids","Stripped.Sequence","PG.Quantity")]
df_pick2 <- unique(df_pick1) # from 11934754 to 10669425 rows
df_distinct <- df_pick2 %>% 
  distinct(File.Name,Protein.Group,Protein.Ids,Stripped.Sequence,PG.Quantity, .keep_all = T) 
df_wide <- dcast(df_distinct, Protein.Group + Stripped.Sequence ~ File.Name, value.var = "PG.Quantity", fun.aggregate = max)
is.na(df_wide)<-sapply(df_wide, is.infinite)
colnames(df_wide) <- gsub("(.*)lyumg_thyroid_complex_","",colnames(df_wide))
colnames(df_wide) <- gsub(".raw","",colnames(df_wide))
aggdata_max <- aggregate(df_wide,
                         by = list(df_wide$Stripped.Sequence),
                         FUN = max,
                         na.rm = TRUE)
is.na(aggdata_max)<-sapply(aggdata_max, is.infinite)
aggdata_max_filter <- aggdata_max %>%
  filter(!(rowSums(is.na(aggdata_max[, 3:ncol(aggdata_max)])) == (ncol(aggdata_max) - 2)))
aggdata_max_final <- aggdata_max_filter[,2:ncol(aggdata_max_filter)]
colnames(aggdata_max_final)[1] <- "protein_id"
colnames(aggdata_max_final)[2] <- "peptide_id"
aggdata_max_final <- unique(aggdata_max_final) # 49588 rows
rio::export(aggdata_max_final,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrix_Thyroid_cellline_all_FDR005_20240906.tsv")
TPC1 <- aggdata_max_final[,grep("TPC-1_NO1", colnames(aggdata_max_final))]
TPC1 <- cbind(aggdata_max_final[,1:2],TPC1)
TPC1 <- TPC1 %>%
  filter(!(rowSums(is.na(TPC1[, 3:ncol(TPC1)])) == (ncol(TPC1) - 2))) # 46462
TPC2 <- aggdata_max_final[,grep("TPC-1_NO2", colnames(aggdata_max_final))]
TPC2 <- cbind(aggdata_max_final[,1:2],TPC2)
TPC2 <- TPC2 %>%
  filter(!(rowSums(is.na(TPC2[, 3:ncol(TPC2)])) == (ncol(TPC2) - 2))) # 47391
TPC3 <- aggdata_max_final[,grep("TPC-1_NO3", colnames(aggdata_max_final))]
TPC3 <- cbind(aggdata_max_final[,1:2],TPC3)
TPC3 <- TPC3 %>%
  filter(!(rowSums(is.na(TPC3[, 3:ncol(TPC3)])) == (ncol(TPC3) - 2))) # 46821
rio::export(TPC1,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrixwithall_TPC-1NO1_FDR005_20240906.tsv")
rio::export(TPC2,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrixwithall_TPC-1NO2_FDR005_20240906.tsv")
rio::export(TPC3,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrixwithall_TPC-1NO3_FDR005_20240906.tsv")
Nthy1 <- aggdata_max_final[,grep("Nthy_NO1", colnames(aggdata_max_final))]
Nthy1 <- cbind(aggdata_max_final[,1:2],Nthy1)
Nthy1 <- Nthy1 %>%
  filter(!(rowSums(is.na(Nthy1[, 3:ncol(Nthy1)])) == (ncol(Nthy1) - 2))) # 48096
Nthy2 <- aggdata_max_final[,grep("Nthy_NO2", colnames(aggdata_max_final))]
Nthy2 <- cbind(aggdata_max_final[,1:2],Nthy2)
Nthy2 <- Nthy2 %>%
  filter(!(rowSums(is.na(Nthy2[, 3:ncol(Nthy2)])) == (ncol(Nthy2) - 2))) # 47684
Nthy3 <- aggdata_max_final[,grep("Nthy_NO3", colnames(aggdata_max_final))]
Nthy3 <- cbind(aggdata_max_final[,1:2],Nthy3)
Nthy3 <- Nthy3 %>%
  filter(!(rowSums(is.na(Nthy3[, 3:ncol(Nthy3)])) == (ncol(Nthy3) - 2))) # 47511
rio::export(Nthy1,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrixwithall_NthyNO1_FDR005_20240906.tsv")
rio::export(Nthy2,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrixwithall_NthyNO2_FDR005_20240906.tsv")
rio::export(Nthy3,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrixwithall_NthyNO3_FDR005_20240906.tsv")
FTC1331 <- aggdata_max_final[,grep("FTC133_NO1", colnames(aggdata_max_final))]
FTC1331 <- cbind(aggdata_max_final[,1:2],FTC1331)
FTC1331 <- FTC1331 %>%
  filter(!(rowSums(is.na(FTC1331[, 3:ncol(FTC1331)])) == (ncol(FTC1331) - 2))) # 47029
FTC1332 <- aggdata_max_final[,grep("FTC133_NO2", colnames(aggdata_max_final))]
FTC1332 <- cbind(aggdata_max_final[,1:2],FTC1332)
FTC1332 <- FTC1332 %>%
  filter(!(rowSums(is.na(FTC1332[, 3:ncol(FTC1332)])) == (ncol(FTC1332) - 2))) # 42003
FTC1333 <- aggdata_max_final[,grep("FTC133_NO3", colnames(aggdata_max_final))]
FTC1333 <- cbind(aggdata_max_final[,1:2],FTC1333)
FTC1333 <- FTC1333 %>%
  filter(!(rowSums(is.na(FTC1333[, 3:ncol(FTC1333)])) == (ncol(FTC1333) - 2))) # 47421
rio::export(FTC1331,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrixwithall_FTC133NO1_FDR005_20240906.tsv")
rio::export(FTC1332,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrixwithall_FTC133NO2_FDR005_20240906.tsv")
rio::export(FTC1333,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrixwiyhall_FTC133NO3_FDR005_20240906.tsv")
FTC2381 <- aggdata_max_final[,grep("FTC238_NO1", colnames(aggdata_max_final))]
FTC2381 <- cbind(aggdata_max_final[,1:2],FTC2381)
FTC2381 <- FTC2381 %>%
  filter(!(rowSums(is.na(FTC2381[, 3:ncol(FTC2381)])) == (ncol(FTC2381) - 2))) # 46046
FTC2382 <- aggdata_max_final[,grep("FTC238_NO2", colnames(aggdata_max_final))]
FTC2382 <- cbind(aggdata_max_final[,1:2],FTC2382)
FTC2382 <- FTC2382 %>%
  filter(!(rowSums(is.na(FTC2382[, 3:ncol(FTC2382)])) == (ncol(FTC2382) - 2))) # 45896
FTC2383 <- aggdata_max_final[,grep("FTC238_NO3", colnames(aggdata_max_final))]
FTC2383 <- cbind(aggdata_max_final[,1:2],FTC2383)
FTC2383 <- FTC2383 %>%
  filter(!(rowSums(is.na(FTC2383[, 3:ncol(FTC2383)])) == (ncol(FTC2383) - 2))) # 45908
rio::export(FTC2381,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrixwithall_FTC238NO1_FDR005_20240906.tsv")
rio::export(FTC2382,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrixwithall_FTC238NO2_FDR005_20240906.tsv")
rio::export(FTC2383,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\peptidematrixwithall_FTC238NO3_FDR005_20240906.tsv")


TPC1_Dall <- TPC1[,-grep("all", colnames(TPC1))]
TPC1_Dall <- TPC1_Dall[,-grep("POOL", colnames(TPC1))]
TPC1_Dall <- TPC1_Dall %>%
  filter(!(rowSums(is.na(TPC1_Dall[, 3:ncol(TPC1_Dall)])) == (ncol(TPC1_Dall) - 2))) # 45856
rio::export(TPC1_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\3_deletall\\peptidematrixforCCprofiler_TPC1_FDR005_20240906.tsv")

TPC2_Dall <- TPC2[,-grep("all", colnames(TPC2))]
TPC2_Dall <- TPC2_Dall[,-grep("POOL", colnames(TPC2))]
TPC2_Dall <- TPC2_Dall %>%
  filter(!(rowSums(is.na(TPC2_Dall[, 3:ncol(TPC2_Dall)])) == (ncol(TPC2_Dall) - 2))) # 46941
rio::export(TPC2_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\3_deletall\\peptidematrixforCCprofiler_TPC2_FDR005_20240906.tsv")

TPC3_Dall <- TPC3[,-grep("all", colnames(TPC3))]
TPC3_Dall <- TPC3_Dall[,-grep("POOL", colnames(TPC3))]
TPC3_Dall <- TPC3_Dall %>%
  filter(!(rowSums(is.na(TPC3_Dall[, 3:ncol(TPC3_Dall)])) == (ncol(TPC3_Dall) - 2))) # 46172
rio::export(TPC3_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\3_deletall\\peptidematrixforCCprofiler_TPC3_FDR005_20240906.tsv")


Nthy1_Dall <- Nthy1[,-grep("all", colnames(Nthy1))]
Nthy1_Dall <- Nthy1_Dall[,-grep("POOL", colnames(Nthy1))]
Nthy1_Dall <- Nthy1_Dall %>%
  filter(!(rowSums(is.na(Nthy1_Dall[, 3:ncol(Nthy1_Dall)])) == (ncol(Nthy1_Dall) - 2))) # 47728
rio::export(Nthy1_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\3_deletall\\peptidematrixforCCprofiler_Nthy1_FDR005_20240906.tsv")

Nthy2_Dall <- Nthy2[,-grep("all", colnames(Nthy2))]
Nthy2_Dall <- Nthy2_Dall[,-grep("POOL", colnames(Nthy2))]
Nthy2_Dall <- Nthy2_Dall %>%
  filter(!(rowSums(is.na(Nthy2_Dall[, 3:ncol(Nthy2_Dall)])) == (ncol(Nthy2_Dall) - 2))) # 47251
rio::export(Nthy2_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\3_deletall\\peptidematrixforCCprofiler_Nthy2_FDR005_20240906.tsv")

Nthy3_Dall <- Nthy3[,-grep("all", colnames(Nthy3))]
Nthy3_Dall <- Nthy3_Dall[,-grep("POOL", colnames(Nthy3))]
Nthy3_Dall <- Nthy3_Dall %>%
  filter(!(rowSums(is.na(Nthy3_Dall[, 3:ncol(Nthy3_Dall)])) == (ncol(Nthy3_Dall) - 2))) # 47154
rio::export(Nthy3_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\3_deletall\\peptidematrixforCCprofiler_Nthy3_FDR005_20240906.tsv")

FTC1331_Dall <- FTC1331[,-grep("all", colnames(FTC1331))]
FTC1331_Dall <- FTC1331_Dall[,-grep("POOL", colnames(FTC1331))]
FTC1331_Dall <- FTC1331_Dall %>%
  filter(!(rowSums(is.na(FTC1331_Dall[, 3:ncol(FTC1331_Dall)])) == (ncol(FTC1331_Dall) - 2))) # 46890
rio::export(FTC1331_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\3_deletall\\peptidematrixforCCprofiler_FTC133NO1_FDR005_20240906.tsv")

FTC1332_Dall <- FTC1332[,-grep("all", colnames(FTC1332))]
FTC1332_Dall <- FTC1332_Dall[,-grep("POOL", colnames(FTC1332))]
FTC1332_Dall <- FTC1332_Dall %>%
  filter(!(rowSums(is.na(FTC1332_Dall[, 3:ncol(FTC1332_Dall)])) == (ncol(FTC1332_Dall) - 2))) # 39671
rio::export(FTC1332_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\3_deletall\\peptidematrixforCCprofiler_FTC133NO2_FDR005_20240906.tsv")

FTC1333_Dall <- FTC1333[,-grep("all", colnames(FTC1333))]
FTC1333_Dall <- FTC1333_Dall[,-grep("POOL", colnames(FTC1333))]
FTC1333_Dall <- FTC1333_Dall %>%
  filter(!(rowSums(is.na(FTC1333_Dall[, 3:ncol(FTC1333_Dall)])) == (ncol(FTC1333_Dall) - 2))) # 47358
rio::export(FTC1333_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\3_deletall\\peptidematrixforCCprofiler_FTC133NO3_FDR005_20240906.tsv")

FTC2381_Dall <- FTC2381[,-grep("all", colnames(FTC2381))]
FTC2381_Dall <- FTC2381_Dall[,-grep("POOL", colnames(FTC2381))]
FTC2381_Dall <- FTC2381_Dall %>%
  filter(!(rowSums(is.na(FTC2381_Dall[, 3:ncol(FTC2381_Dall)])) == (ncol(FTC2381_Dall) - 2))) # 45522
rio::export(FTC2381_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\3_deletall\\peptidematrixforCCprofiler_FTC238NO1_FDR005_20240906.tsv")

FTC2382_Dall <- FTC2382[,-grep("all", colnames(FTC2382))]
FTC2382_Dall <- FTC2382_Dall[,-grep("POOL", colnames(FTC2382))]
FTC2382_Dall <- FTC2382_Dall %>%
  filter(!(rowSums(is.na(FTC2382_Dall[, 3:ncol(FTC2382_Dall)])) == (ncol(FTC2382_Dall) - 2))) # 45577
rio::export(FTC2382_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\3_deletall\\peptidematrixforCCprofiler_FTC238NO2_FDR005_20240906.tsv")

FTC2383_Dall <- FTC2383[,-grep("all", colnames(FTC2383))]
FTC2383_Dall <- FTC2383_Dall[,-grep("POOL", colnames(FTC2383))]
FTC2383_Dall <- FTC2383_Dall %>%
  filter(!(rowSums(is.na(FTC2383_Dall[, 3:ncol(FTC2383_Dall)])) == (ncol(FTC2383_Dall) - 2))) # 45485
rio::export(FTC2383_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\3_deletall\\peptidematrixforCCprofiler_FTC238NO3_FDR005_20240906.tsv")
#######

TPC1_PROT <- as.data.frame(TPC1_Dall$protein_id)
colnames(TPC1_PROT)[1] <- "Prot"
TPC1_PROT <- unique(TPC1_PROT) # 6176
rio::export(TPC1_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\4_proteinname\\TPC1_protein.tsv")
TPC2_PROT <- as.data.frame(TPC2_Dall$protein_id)
TPC2_PROT <- unique(TPC2_PROT) # 6345
colnames(TPC2_PROT)[1] <- "Prot"
rio::export(TPC2_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\4_proteinname\\TPC2_protein.tsv")
TPC3_PROT <- as.data.frame(TPC3_Dall$protein_id)
TPC3_PROT <- unique(TPC3_PROT) # 6210
colnames(TPC3_PROT)[1] <- "Prot"
rio::export(TPC3_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\4_proteinname\\TPC3_protein.tsv")

Nthy1_PROT <- as.data.frame(Nthy1_Dall$protein_id)
colnames(Nthy1_PROT)[1] <- "Prot"
Nthy1_PROT <- unique(Nthy1_PROT) # 6409
rio::export(Nthy1_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\4_proteinname\\Nthy1_protein.tsv")
Nthy2_PROT <- as.data.frame(Nthy2_Dall$protein_id)
Nthy2_PROT <- unique(Nthy2_PROT) # 6413
colnames(Nthy2_PROT)[1] <- "Prot"
rio::export(Nthy2_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\4_proteinname\\Nthy2_protein.tsv")
Nthy3_PROT <- as.data.frame(Nthy3_Dall$protein_id)
Nthy3_PROT <- unique(Nthy3_PROT) # 6419
colnames(Nthy3_PROT)[1] <- "Prot"
rio::export(Nthy3_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\4_proteinname\\Nthy3_protein.tsv")

FTC1331_PROT <- as.data.frame(FTC1331_Dall$protein_id)
colnames(FTC1331_PROT)[1] <- "Prot"
FTC1331_PROT <- unique(FTC1331_PROT) # 6366
rio::export(FTC1331_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\4_proteinname\\FTC133NO1_protein.tsv")
FTC1332_PROT <- as.data.frame(FTC1332_Dall$protein_id)
FTC1332_PROT <- unique(FTC1332_PROT) # 5834
colnames(FTC1332_PROT)[1] <- "Prot"
rio::export(FTC1332_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\4_proteinname\\FTC133NO2_protein.tsv")
FTC1333_PROT <- as.data.frame(FTC1333_Dall$protein_id)
FTC1333_PROT <- unique(FTC1333_PROT) # 6443
colnames(FTC1333_PROT)[1] <- "Prot"
rio::export(FTC1333_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\4_proteinname\\FTC133NO3_protein.tsv")

FTC2381_PROT <- as.data.frame(FTC2381_Dall$protein_id)
colnames(FTC2381_PROT)[1] <- "Prot"
FTC2381_PROT <- unique(FTC2381_PROT) # 6280
rio::export(FTC2381_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\4_proteinname\\FTC238NO1_protein.tsv")
FTC2382_PROT <- as.data.frame(FTC2382_Dall$protein_id)
FTC2382_PROT <- unique(FTC2382_PROT) # 6285
colnames(FTC2382_PROT)[1] <- "Prot"
rio::export(FTC2382_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\4_proteinname\\FTC238NO2_protein.tsv")
FTC2383_PROT <- as.data.frame(FTC2383_Dall$protein_id)
FTC2383_PROT <- unique(FTC2383_PROT) # 6264
colnames(FTC2383_PROT)[1] <- "Prot"
rio::export(FTC2383_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\4_proteinname\\FTC238NO3_protein.tsv")


TPC1_FRAC <- as.data.frame(colnames(TPC1_Dall)[3:length(colnames(TPC1_Dall))])
colnames(TPC1_FRAC)[1] <- "filename"
TPC1_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", TPC1_FRAC$filename)
rio::export(TPC1_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\5_fraction\\TPC1_fraction.tsv")
TPC2_FRAC <- as.data.frame(colnames(TPC2_Dall)[3:length(colnames(TPC2_Dall))])
colnames(TPC2_FRAC)[1] <- "filename"
TPC2_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", TPC2_FRAC$filename)
rio::export(TPC2_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\5_fraction\\TPC2_fraction.tsv")
TPC3_FRAC <- as.data.frame(colnames(TPC3_Dall)[3:length(colnames(TPC3_Dall))])
colnames(TPC3_FRAC)[1] <- "filename"
TPC3_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", TPC3_FRAC$filename)
rio::export(TPC3_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\5_fraction\\TPC3_fraction.tsv")

Nthy1_FRAC <- as.data.frame(colnames(Nthy1_Dall)[3:length(colnames(Nthy1_Dall))])
colnames(Nthy1_FRAC)[1] <- "filename"
Nthy1_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", Nthy1_FRAC$filename)
rio::export(Nthy1_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\5_fraction\\Nthy1_fraction.tsv")
Nthy2_FRAC <- as.data.frame(colnames(Nthy2_Dall)[3:length(colnames(Nthy2_Dall))])
colnames(Nthy2_FRAC)[1] <- "filename"
Nthy2_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", Nthy2_FRAC$filename)
rio::export(Nthy2_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\5_fraction\\Nthy2_fraction.tsv")
Nthy3_FRAC <- as.data.frame(colnames(Nthy3_Dall)[3:length(colnames(Nthy3_Dall))])
colnames(Nthy3_FRAC)[1] <- "filename"
Nthy3_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", Nthy3_FRAC$filename)
rio::export(Nthy3_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\5_fraction\\Nthy3_fraction.tsv")

FTC1331_FRAC <- as.data.frame(colnames(FTC1331_Dall)[3:length(colnames(FTC1331_Dall))])
colnames(FTC1331_FRAC)[1] <- "filename"
FTC1331_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", FTC1331_FRAC$filename)
rio::export(FTC1331_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\5_fraction\\FTC133NO1_fraction.tsv")
FTC1332_FRAC <- as.data.frame(colnames(FTC1332_Dall)[3:length(colnames(FTC1332_Dall))])
colnames(FTC1332_FRAC)[1] <- "filename"
FTC1332_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", FTC1332_FRAC$filename)
rio::export(FTC1332_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\5_fraction\\FTC133NO2_fraction.tsv")
FTC1333_FRAC <- as.data.frame(colnames(FTC1333_Dall)[3:length(colnames(FTC1333_Dall))])
colnames(FTC1333_FRAC)[1] <- "filename"
FTC1333_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", FTC1333_FRAC$filename)
rio::export(FTC1333_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\5_fraction\\FTC133NO3_fraction.tsv")

FTC2381_FRAC <- as.data.frame(colnames(FTC2381_Dall)[3:length(colnames(FTC2381_Dall))])
colnames(FTC2381_FRAC)[1] <- "filename"
FTC2381_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", FTC2381_FRAC$filename)
rio::export(FTC2381_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\5_fraction\\FTC238NO1_fraction.tsv")
FTC2382_FRAC <- as.data.frame(colnames(FTC2382_Dall)[3:length(colnames(FTC2382_Dall))])
colnames(FTC2382_FRAC)[1] <- "filename"
FTC2382_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", FTC2382_FRAC$filename)
rio::export(FTC2382_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\5_fraction\\FTC238NO2_fraction.tsv")
FTC2383_FRAC <- as.data.frame(colnames(FTC2383_Dall)[3:length(colnames(FTC2383_Dall))])
colnames(FTC2383_FRAC)[1] <- "filename"
FTC2383_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", FTC2383_FRAC$filename)
rio::export(FTC2383_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\5_fraction\\FTC238NO3_fraction.tsv")

save.image("E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR005\\Thyroid_cell_line_preprocessing_20240906.RData")

### FDR001
df <- fread("E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\report.tsv")
df_pick <- df[df$Q.Value <= 0.05 & df$Global.Q.Value <= 0.05 & df$PG.Q.Value <= 0.05 & df$Global.PG.Q.Value <= 0.05 & df$Protein.Q.Value <= 0.05,]
df_pick1 <- df_pick[,c("File.Name", "Protein.Group", "Protein.Ids","Stripped.Sequence","PG.Quantity")]
df_pick2 <- unique(df_pick1) # from 11934754 to 10669425 rows
df_distinct <- df_pick2 %>% 
  distinct(File.Name,Protein.Group,Protein.Ids,Stripped.Sequence,PG.Quantity, .keep_all = T) 
df_wide <- dcast(df_distinct, Protein.Group + Stripped.Sequence ~ File.Name, value.var = "PG.Quantity", fun.aggregate = max)
is.na(df_wide)<-sapply(df_wide, is.infinite)
colnames(df_wide) <- gsub("(.*)lyumg_thyroid_complex_","",colnames(df_wide))
colnames(df_wide) <- gsub(".raw","",colnames(df_wide))
aggdata_max <- aggregate(df_wide,
                         by = list(df_wide$Stripped.Sequence),
                         FUN = max,
                         na.rm = TRUE)
is.na(aggdata_max)<-sapply(aggdata_max, is.infinite)
aggdata_max_filter <- aggdata_max %>%
  filter(!(rowSums(is.na(aggdata_max[, 3:ncol(aggdata_max)])) == (ncol(aggdata_max) - 2)))
aggdata_max_final <- aggdata_max_filter[,2:ncol(aggdata_max_filter)]
colnames(aggdata_max_final)[1] <- "protein_id"
colnames(aggdata_max_final)[2] <- "peptide_id"
aggdata_max_final <- unique(aggdata_max_final) # 49588 rows
rio::export(aggdata_max_final,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrix_Thyroid_cellline_all_FDR001_20240906.tsv")
TPC1 <- aggdata_max_final[,grep("TPC-1_NO1", colnames(aggdata_max_final))]
TPC1 <- cbind(aggdata_max_final[,1:2],TPC1)
TPC1 <- TPC1 %>%
  filter(!(rowSums(is.na(TPC1[, 3:ncol(TPC1)])) == (ncol(TPC1) - 2))) # 46462
TPC2 <- aggdata_max_final[,grep("TPC-1_NO2", colnames(aggdata_max_final))]
TPC2 <- cbind(aggdata_max_final[,1:2],TPC2)
TPC2 <- TPC2 %>%
  filter(!(rowSums(is.na(TPC2[, 3:ncol(TPC2)])) == (ncol(TPC2) - 2))) # 47391
TPC3 <- aggdata_max_final[,grep("TPC-1_NO3", colnames(aggdata_max_final))]
TPC3 <- cbind(aggdata_max_final[,1:2],TPC3)
TPC3 <- TPC3 %>%
  filter(!(rowSums(is.na(TPC3[, 3:ncol(TPC3)])) == (ncol(TPC3) - 2))) # 46821
rio::export(TPC1,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrixwithall_TPC-1NO1_FDR001_20240906.tsv")
rio::export(TPC2,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrixwithall_TPC-1NO2_FDR001_20240906.tsv")
rio::export(TPC3,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrixwithall_TPC-1NO3_FDR001_20240906.tsv")
Nthy1 <- aggdata_max_final[,grep("Nthy_NO1", colnames(aggdata_max_final))]
Nthy1 <- cbind(aggdata_max_final[,1:2],Nthy1)
Nthy1 <- Nthy1 %>%
  filter(!(rowSums(is.na(Nthy1[, 3:ncol(Nthy1)])) == (ncol(Nthy1) - 2))) # 48096
Nthy2 <- aggdata_max_final[,grep("Nthy_NO2", colnames(aggdata_max_final))]
Nthy2 <- cbind(aggdata_max_final[,1:2],Nthy2)
Nthy2 <- Nthy2 %>%
  filter(!(rowSums(is.na(Nthy2[, 3:ncol(Nthy2)])) == (ncol(Nthy2) - 2))) # 47684
Nthy3 <- aggdata_max_final[,grep("Nthy_NO3", colnames(aggdata_max_final))]
Nthy3 <- cbind(aggdata_max_final[,1:2],Nthy3)
Nthy3 <- Nthy3 %>%
  filter(!(rowSums(is.na(Nthy3[, 3:ncol(Nthy3)])) == (ncol(Nthy3) - 2))) # 47511
rio::export(Nthy1,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrixwithall_NthyNO1_FDR001_20240906.tsv")
rio::export(Nthy2,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrixwithall_NthyNO2_FDR001_20240906.tsv")
rio::export(Nthy3,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrixwithall_NthyNO3_FDR001_20240906.tsv")
FTC1331 <- aggdata_max_final[,grep("FTC133_NO1", colnames(aggdata_max_final))]
FTC1331 <- cbind(aggdata_max_final[,1:2],FTC1331)
FTC1331 <- FTC1331 %>%
  filter(!(rowSums(is.na(FTC1331[, 3:ncol(FTC1331)])) == (ncol(FTC1331) - 2))) # 47029
FTC1332 <- aggdata_max_final[,grep("FTC133_NO2", colnames(aggdata_max_final))]
FTC1332 <- cbind(aggdata_max_final[,1:2],FTC1332)
FTC1332 <- FTC1332 %>%
  filter(!(rowSums(is.na(FTC1332[, 3:ncol(FTC1332)])) == (ncol(FTC1332) - 2))) # 42003
FTC1333 <- aggdata_max_final[,grep("FTC133_NO3", colnames(aggdata_max_final))]
FTC1333 <- cbind(aggdata_max_final[,1:2],FTC1333)
FTC1333 <- FTC1333 %>%
  filter(!(rowSums(is.na(FTC1333[, 3:ncol(FTC1333)])) == (ncol(FTC1333) - 2))) # 47421
rio::export(FTC1331,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrixwithall_FTC133NO1_FDR001_20240906.tsv")
rio::export(FTC1332,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrixwithall_FTC133NO2_FDR001_20240906.tsv")
rio::export(FTC1333,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrixwiyhall_FTC133NO3_FDR001_20240906.tsv")
FTC2381 <- aggdata_max_final[,grep("FTC238_NO1", colnames(aggdata_max_final))]
FTC2381 <- cbind(aggdata_max_final[,1:2],FTC2381)
FTC2381 <- FTC2381 %>%
  filter(!(rowSums(is.na(FTC2381[, 3:ncol(FTC2381)])) == (ncol(FTC2381) - 2))) # 46046
FTC2382 <- aggdata_max_final[,grep("FTC238_NO2", colnames(aggdata_max_final))]
FTC2382 <- cbind(aggdata_max_final[,1:2],FTC2382)
FTC2382 <- FTC2382 %>%
  filter(!(rowSums(is.na(FTC2382[, 3:ncol(FTC2382)])) == (ncol(FTC2382) - 2))) # 45896
FTC2383 <- aggdata_max_final[,grep("FTC238_NO3", colnames(aggdata_max_final))]
FTC2383 <- cbind(aggdata_max_final[,1:2],FTC2383)
FTC2383 <- FTC2383 %>%
  filter(!(rowSums(is.na(FTC2383[, 3:ncol(FTC2383)])) == (ncol(FTC2383) - 2))) # 45908
rio::export(FTC2381,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrixwithall_FTC238NO1_FDR001_20240906.tsv")
rio::export(FTC2382,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrixwithall_FTC238NO2_FDR001_20240906.tsv")
rio::export(FTC2383,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\peptidematrixwithall_FTC238NO3_FDR001_20240906.tsv")


TPC1_Dall <- TPC1[,-grep("all", colnames(TPC1))]
TPC1_Dall <- TPC1_Dall[,-grep("POOL", colnames(TPC1))]
TPC1_Dall <- TPC1_Dall %>%
  filter(!(rowSums(is.na(TPC1_Dall[, 3:ncol(TPC1_Dall)])) == (ncol(TPC1_Dall) - 2))) # 45856
rio::export(TPC1_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\3_deletall\\peptidematrixforCCprofiler_TPC1_FDR001_20240906.tsv")

TPC2_Dall <- TPC2[,-grep("all", colnames(TPC2))]
TPC2_Dall <- TPC2_Dall[,-grep("POOL", colnames(TPC2))]
TPC2_Dall <- TPC2_Dall %>%
  filter(!(rowSums(is.na(TPC2_Dall[, 3:ncol(TPC2_Dall)])) == (ncol(TPC2_Dall) - 2))) # 46941
rio::export(TPC2_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\3_deletall\\peptidematrixforCCprofiler_TPC2_FDR001_20240906.tsv")

TPC3_Dall <- TPC3[,-grep("all", colnames(TPC3))]
TPC3_Dall <- TPC3_Dall[,-grep("POOL", colnames(TPC3))]
TPC3_Dall <- TPC3_Dall %>%
  filter(!(rowSums(is.na(TPC3_Dall[, 3:ncol(TPC3_Dall)])) == (ncol(TPC3_Dall) - 2))) # 46172
rio::export(TPC3_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\3_deletall\\peptidematrixforCCprofiler_TPC3_FDR001_20240906.tsv")


Nthy1_Dall <- Nthy1[,-grep("all", colnames(Nthy1))]
Nthy1_Dall <- Nthy1_Dall[,-grep("POOL", colnames(Nthy1))]
Nthy1_Dall <- Nthy1_Dall %>%
  filter(!(rowSums(is.na(Nthy1_Dall[, 3:ncol(Nthy1_Dall)])) == (ncol(Nthy1_Dall) - 2))) # 47728
rio::export(Nthy1_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\3_deletall\\peptidematrixforCCprofiler_Nthy1_FDR001_20240906.tsv")

Nthy2_Dall <- Nthy2[,-grep("all", colnames(Nthy2))]
Nthy2_Dall <- Nthy2_Dall[,-grep("POOL", colnames(Nthy2))]
Nthy2_Dall <- Nthy2_Dall %>%
  filter(!(rowSums(is.na(Nthy2_Dall[, 3:ncol(Nthy2_Dall)])) == (ncol(Nthy2_Dall) - 2))) # 47251
rio::export(Nthy2_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\3_deletall\\peptidematrixforCCprofiler_Nthy2_FDR001_20240906.tsv")

Nthy3_Dall <- Nthy3[,-grep("all", colnames(Nthy3))]
Nthy3_Dall <- Nthy3_Dall[,-grep("POOL", colnames(Nthy3))]
Nthy3_Dall <- Nthy3_Dall %>%
  filter(!(rowSums(is.na(Nthy3_Dall[, 3:ncol(Nthy3_Dall)])) == (ncol(Nthy3_Dall) - 2))) # 47154
rio::export(Nthy3_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\3_deletall\\peptidematrixforCCprofiler_Nthy3_FDR001_20240906.tsv")

FTC1331_Dall <- FTC1331[,-grep("all", colnames(FTC1331))]
FTC1331_Dall <- FTC1331_Dall[,-grep("POOL", colnames(FTC1331))]
FTC1331_Dall <- FTC1331_Dall %>%
  filter(!(rowSums(is.na(FTC1331_Dall[, 3:ncol(FTC1331_Dall)])) == (ncol(FTC1331_Dall) - 2))) # 46890
rio::export(FTC1331_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\3_deletall\\peptidematrixforCCprofiler_FTC133NO1_FDR001_20240906.tsv")

FTC1332_Dall <- FTC1332[,-grep("all", colnames(FTC1332))]
FTC1332_Dall <- FTC1332_Dall[,-grep("POOL", colnames(FTC1332))]
FTC1332_Dall <- FTC1332_Dall %>%
  filter(!(rowSums(is.na(FTC1332_Dall[, 3:ncol(FTC1332_Dall)])) == (ncol(FTC1332_Dall) - 2))) # 39671
rio::export(FTC1332_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\3_deletall\\peptidematrixforCCprofiler_FTC133NO2_FDR001_20240906.tsv")

FTC1333_Dall <- FTC1333[,-grep("all", colnames(FTC1333))]
FTC1333_Dall <- FTC1333_Dall[,-grep("POOL", colnames(FTC1333))]
FTC1333_Dall <- FTC1333_Dall %>%
  filter(!(rowSums(is.na(FTC1333_Dall[, 3:ncol(FTC1333_Dall)])) == (ncol(FTC1333_Dall) - 2))) # 47358
rio::export(FTC1333_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\3_deletall\\peptidematrixforCCprofiler_FTC133NO3_FDR001_20240906.tsv")

FTC2381_Dall <- FTC2381[,-grep("all", colnames(FTC2381))]
FTC2381_Dall <- FTC2381_Dall[,-grep("POOL", colnames(FTC2381))]
FTC2381_Dall <- FTC2381_Dall %>%
  filter(!(rowSums(is.na(FTC2381_Dall[, 3:ncol(FTC2381_Dall)])) == (ncol(FTC2381_Dall) - 2))) # 45522
rio::export(FTC2381_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\3_deletall\\peptidematrixforCCprofiler_FTC238NO1_FDR001_20240906.tsv")

FTC2382_Dall <- FTC2382[,-grep("all", colnames(FTC2382))]
FTC2382_Dall <- FTC2382_Dall[,-grep("POOL", colnames(FTC2382))]
FTC2382_Dall <- FTC2382_Dall %>%
  filter(!(rowSums(is.na(FTC2382_Dall[, 3:ncol(FTC2382_Dall)])) == (ncol(FTC2382_Dall) - 2))) # 45577
rio::export(FTC2382_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\3_deletall\\peptidematrixforCCprofiler_FTC238NO2_FDR001_20240906.tsv")

FTC2383_Dall <- FTC2383[,-grep("all", colnames(FTC2383))]
FTC2383_Dall <- FTC2383_Dall[,-grep("POOL", colnames(FTC2383))]
FTC2383_Dall <- FTC2383_Dall %>%
  filter(!(rowSums(is.na(FTC2383_Dall[, 3:ncol(FTC2383_Dall)])) == (ncol(FTC2383_Dall) - 2))) # 45485
rio::export(FTC2383_Dall,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\3_deletall\\peptidematrixforCCprofiler_FTC238NO3_FDR001_20240906.tsv")


TPC1_PROT <- as.data.frame(TPC1_Dall$protein_id)
colnames(TPC1_PROT)[1] <- "Prot"
TPC1_PROT <- unique(TPC1_PROT) # 6176
rio::export(TPC1_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\4_proteinname\\TPC1_protein.tsv")
TPC2_PROT <- as.data.frame(TPC2_Dall$protein_id)
TPC2_PROT <- unique(TPC2_PROT) # 6345
colnames(TPC2_PROT)[1] <- "Prot"
rio::export(TPC2_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\4_proteinname\\TPC2_protein.tsv")
TPC3_PROT <- as.data.frame(TPC3_Dall$protein_id)
TPC3_PROT <- unique(TPC3_PROT) # 6210
colnames(TPC3_PROT)[1] <- "Prot"
rio::export(TPC3_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\4_proteinname\\TPC3_protein.tsv")

Nthy1_PROT <- as.data.frame(Nthy1_Dall$protein_id)
colnames(Nthy1_PROT)[1] <- "Prot"
Nthy1_PROT <- unique(Nthy1_PROT) # 6409
rio::export(Nthy1_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\4_proteinname\\Nthy1_protein.tsv")
Nthy2_PROT <- as.data.frame(Nthy2_Dall$protein_id)
Nthy2_PROT <- unique(Nthy2_PROT) # 6413
colnames(Nthy2_PROT)[1] <- "Prot"
rio::export(Nthy2_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\4_proteinname\\Nthy2_protein.tsv")
Nthy3_PROT <- as.data.frame(Nthy3_Dall$protein_id)
Nthy3_PROT <- unique(Nthy3_PROT) # 6419
colnames(Nthy3_PROT)[1] <- "Prot"
rio::export(Nthy3_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\4_proteinname\\Nthy3_protein.tsv")

FTC1331_PROT <- as.data.frame(FTC1331_Dall$protein_id)
colnames(FTC1331_PROT)[1] <- "Prot"
FTC1331_PROT <- unique(FTC1331_PROT) # 6366
rio::export(FTC1331_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\4_proteinname\\FTC133NO1_protein.tsv")
FTC1332_PROT <- as.data.frame(FTC1332_Dall$protein_id)
FTC1332_PROT <- unique(FTC1332_PROT) # 5834
colnames(FTC1332_PROT)[1] <- "Prot"
rio::export(FTC1332_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\4_proteinname\\FTC133NO2_protein.tsv")
FTC1333_PROT <- as.data.frame(FTC1333_Dall$protein_id)
FTC1333_PROT <- unique(FTC1333_PROT) # 6443
colnames(FTC1333_PROT)[1] <- "Prot"
rio::export(FTC1333_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\4_proteinname\\FTC133NO3_protein.tsv")

FTC2381_PROT <- as.data.frame(FTC2381_Dall$protein_id)
colnames(FTC2381_PROT)[1] <- "Prot"
FTC2381_PROT <- unique(FTC2381_PROT) # 6280
rio::export(FTC2381_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\4_proteinname\\FTC238NO1_protein.tsv")
FTC2382_PROT <- as.data.frame(FTC2382_Dall$protein_id)
FTC2382_PROT <- unique(FTC2382_PROT) # 6285
colnames(FTC2382_PROT)[1] <- "Prot"
rio::export(FTC2382_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\4_proteinname\\FTC238NO2_protein.tsv")
FTC2383_PROT <- as.data.frame(FTC2383_Dall$protein_id)
FTC2383_PROT <- unique(FTC2383_PROT) # 6264
colnames(FTC2383_PROT)[1] <- "Prot"
rio::export(FTC2383_PROT,"E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\4_proteinname\\FTC238NO3_protein.tsv")


TPC1_FRAC <- as.data.frame(colnames(TPC1_Dall)[3:length(colnames(TPC1_Dall))])
colnames(TPC1_FRAC)[1] <- "filename"
TPC1_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", TPC1_FRAC$filename)
rio::export(TPC1_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\5_fraction\\TPC1_fraction.tsv")
TPC2_FRAC <- as.data.frame(colnames(TPC2_Dall)[3:length(colnames(TPC2_Dall))])
colnames(TPC2_FRAC)[1] <- "filename"
TPC2_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", TPC2_FRAC$filename)
rio::export(TPC2_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\5_fraction\\TPC2_fraction.tsv")
TPC3_FRAC <- as.data.frame(colnames(TPC3_Dall)[3:length(colnames(TPC3_Dall))])
colnames(TPC3_FRAC)[1] <- "filename"
TPC3_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", TPC3_FRAC$filename)
rio::export(TPC3_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\5_fraction\\TPC3_fraction.tsv")

Nthy1_FRAC <- as.data.frame(colnames(Nthy1_Dall)[3:length(colnames(Nthy1_Dall))])
colnames(Nthy1_FRAC)[1] <- "filename"
Nthy1_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", Nthy1_FRAC$filename)
rio::export(Nthy1_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\5_fraction\\Nthy1_fraction.tsv")
Nthy2_FRAC <- as.data.frame(colnames(Nthy2_Dall)[3:length(colnames(Nthy2_Dall))])
colnames(Nthy2_FRAC)[1] <- "filename"
Nthy2_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", Nthy2_FRAC$filename)
rio::export(Nthy2_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\5_fraction\\Nthy2_fraction.tsv")
Nthy3_FRAC <- as.data.frame(colnames(Nthy3_Dall)[3:length(colnames(Nthy3_Dall))])
colnames(Nthy3_FRAC)[1] <- "filename"
Nthy3_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", Nthy3_FRAC$filename)
rio::export(Nthy3_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\5_fraction\\Nthy3_fraction.tsv")

FTC1331_FRAC <- as.data.frame(colnames(FTC1331_Dall)[3:length(colnames(FTC1331_Dall))])
colnames(FTC1331_FRAC)[1] <- "filename"
FTC1331_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", FTC1331_FRAC$filename)
rio::export(FTC1331_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\5_fraction\\FTC133NO1_fraction.tsv")
FTC1332_FRAC <- as.data.frame(colnames(FTC1332_Dall)[3:length(colnames(FTC1332_Dall))])
colnames(FTC1332_FRAC)[1] <- "filename"
FTC1332_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", FTC1332_FRAC$filename)
rio::export(FTC1332_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\5_fraction\\FTC133NO2_fraction.tsv")
FTC1333_FRAC <- as.data.frame(colnames(FTC1333_Dall)[3:length(colnames(FTC1333_Dall))])
colnames(FTC1333_FRAC)[1] <- "filename"
FTC1333_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", FTC1333_FRAC$filename)
rio::export(FTC1333_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\5_fraction\\FTC133NO3_fraction.tsv")

FTC2381_FRAC <- as.data.frame(colnames(FTC2381_Dall)[3:length(colnames(FTC2381_Dall))])
colnames(FTC2381_FRAC)[1] <- "filename"
FTC2381_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", FTC2381_FRAC$filename)
rio::export(FTC2381_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\5_fraction\\FTC238NO1_fraction.tsv")
FTC2382_FRAC <- as.data.frame(colnames(FTC2382_Dall)[3:length(colnames(FTC2382_Dall))])
colnames(FTC2382_FRAC)[1] <- "filename"
FTC2382_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", FTC2382_FRAC$filename)
rio::export(FTC2382_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\5_fraction\\FTC238NO2_fraction.tsv")
FTC2383_FRAC <- as.data.frame(colnames(FTC2383_Dall)[3:length(colnames(FTC2383_Dall))])
colnames(FTC2383_FRAC)[1] <- "filename"
FTC2383_FRAC$fraction_number <- gsub(".*_(\\d+).*", "\\1", FTC2383_FRAC$filename)
rio::export(FTC2383_FRAC, "E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\5_fraction\\FTC238NO3_fraction.tsv")

save.image("E:\\MenggeLYU\\NC\\Thyroid_cellline_all_20240617\\FDR001\\Thyroid_cell_line_preprocessing_20240906.RData")
