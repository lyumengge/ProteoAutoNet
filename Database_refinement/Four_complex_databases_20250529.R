##################################################
## Project: ProteoAutoNet
## Script purpose: Four databases
## Date: 2025-05-29
## Author: MenggeLYU
## Version: 3.0
##################################################


# path and library --------------------------------------------------------
setwd("E://NC//CHINDA")
library(rio)
library(dplyr)
library(tidyr)
library(stringr)
library(corrplot)

# UniProt -----------------------------------------------------------------
uniprot <- rio::import("uniprot_protein_list_for_NC_20240608.tsv")
newcorum_prot_left <- left_join(newcorum_prot_split, uniprot, by = "Prot")
is.na(newcorum_prot_left$Reviewed)
null <- newcorum_prot_left[which(is.na(newcorum_prot_left$Reviewed)),]

# CORUM -------------------------------------------------------------------
newcorum <- rio::import("CORUM download 2022_09_12.xlsx") 
newcorum <- newcorum[which(newcorum$Organism=="Human"),] 
newcorum_prot <- as.data.frame(newcorum$`subunits(UniProt IDs)`)
colnames(newcorum_prot)[1] <- "Prot"
newcorum_prot_split <- newcorum_prot %>% as_tibble() %>% 
  separate_rows(Prot, sep = ";")
newcorum_prot_split <- unique(newcorum_prot_split)

# huMAP -------------------------------------------------------------------
humap <- rio::import("humap2_complexes_20200809.txt")
humap <- humap[,c(1,3)]
colnames(humap)[2] <- "Prot"
humap_split <- humap %>% as_tibble() %>% 
  separate_rows(Prot, sep = " ")
humap_left <- left_join(humap_split, uniprot)
humap_left2 <- humap_left[-grep("^$", humap_left$Prot),]
humap_prot <- unique(humap_left2$Prot)
hunull <- humap_left2[which(is.na(humap_left2$Reviewed)),]
humap_left3 <- humap_left2[-which(is.na(humap_left2$Reviewed)),]
humap_left_prot <- unique(humap_left3$Prot)
humap_left3_twocol <- humap_left3[,1:2]
unique_HuMAP2_IDs <- humap_left3_twocol %>%
  group_by(HuMAP2_ID) %>%
  filter(n_distinct(Prot) == 1) %>%
  pull(HuMAP2_ID)
humap_left3_twocol2 <- humap_left3_twocol[-which(humap_left3_twocol$HuMAP2_ID%in%unique_HuMAP2_IDs),]
length(unique(humap_left3_twocol2$Prot))
humap_left3_twocol3 <- left_join(humap_left3_twocol2,uniprot)
humap_prot <- as.data.frame(unique(humap_left3_twocol3[,2:4]))
table(humap_prot$Reviewed)
table(humap_prot$Sequence)
length(unique(humap_left3_twocol3$HuMAP2_ID))

# Complex Portal ----------------------------------------------------------
cp <- rio::import("complex_portal_9606.tsv")
unique(cp$`Taxonomy identifier`)
cp <- cp[,c(1,19)]
length(unique(cp$Prot))
colnames(cp)[2] <- "Prot"
cp_split <- cp %>% as_tibble() %>% 
  separate_rows(Prot, sep = "\\|")
cp_split$Prot <- gsub("\\(.*?\\)","",cp_split$Prot)
cp_split2 <- cp_split %>% as_tibble() %>% 
  separate_rows(Prot, sep = ",")
cp_split2$Prot <- gsub("\\[","",cp_split2$Prot)
cp_split2$Prot <- gsub("\\]","",cp_split2$Prot)
cp_split2_left <- left_join(cp_split2,uniprot)
cpnull <- cp_split2_left[which(is.na(cp_split2_left$Reviewed)),]
length(unique(cpnull$Prot))
length(unique(cpnull$`#Complex ac`))
cp_left <- cp_split2_left[-which(is.na(cp_split2_left$Reviewed)),]
colnames(cp_left)[1] <- "complex_name"
unique_cp_IDs <- cp_left %>%
  group_by(complex_name) %>%
  filter(n_distinct(Prot) == 1) %>%
  pull(complex_name)
cp_left2 <- cp_left[-which(cp_left$complex_name%in%unique_cp_IDs),]
cp_unique_prot_left2 <- as.data.frame(cp_left2[,2:4])
cp_unique_prot_left2 <- unique(cp_unique_prot_left2)
length(unique(cp_unique_prot_left2$Prot))
table(cp_unique_prot_left2$Sequence)
length(unique(cp_left2$complex_name))

# iRefIndex ---------------------------------------------------------------
# To extract the complexes from iRefIndex which distributed in PSI-MITAB format (Version 2.5), the function named, convert_MITAB_to_complexList() in R package iRef (version 1.13) was used to preprocess the iRefIndex to ‘complex list’ format in 2.13.1 of R before matching UniprotKB identifier from Homo Sapiens. 
iRef <- rio::import("iRefIndex_complex_list_human_uniprot_ID_nolabel_20240608.tsv")
iRef$comma_count <- str_count(iRef$Uniprot_ID, ", ")
iRef_pick <- iRef[,c(1,3)]
iRef_pick$Complex_ID <- paste("iRef", iRef_pick$Complex_ID, sep = "_")
iRef_split <- iRef %>% as_tibble() %>% 
  separate_rows(Uniprot_ID, sep = ", ")
head(iRef_split)
unique_iRef_IDs <- iRef_split %>%
  group_by(Complex_ID) %>%
  filter(n_distinct(Uniprot_ID) == 1) %>%
  pull(Complex_ID)
colnames(iRef_split)[2] <- "Prot"
iRef_left <- iRef_split[-which(iRef_split$Complex_ID%in%unique_iRef_IDs),]
length(unique(iRef_left$Complex_ID)) 
length(unique(iRef_left$Prot)) 
# Random Walk with Restart results - only left node of complexes < 200 proteins
# To select the core proteins in protein interactions, a custom function named random_walk_with_restart has been implemented to carry out the algorithm in Python 3.10.9.
iRef_rwr <- rio::import("output_final.tsv")
iRef_rwr$comma_count <- str_count(iRef_rwr$interacting_proteins, ";")
iRef_rwr$comma_count <- str_count(iRef_rwr$interacting_proteins, ";")
max(iRef_rwr$comma_count)
View(iRef)
iRef$comma_count <- str_count(iRef$Uniprot_ID, ", ")
iRef_split_left <- left_join(iRef_split, uniprot)
iRef_split_left2 <- iRef_split_left[-which(iRef_split_left$Complex_ID%in%unique_iRef_IDs),]
iRef_split_prot <- iRef_split_left2[,2:4]
reviewed_iRef <- iRef_split_left2[which(iRef_split_left2$reviewed=="reviewed"),]
unique_reviewediRef_IDs <- reviewed_iRef %>%
  group_by(Complex_ID) %>%
  filter(n_distinct(Prot) == 1) %>%
  pull(Complex_ID)
iRef_split_prot2 <- unique(iRef_split_prot)
table(iRef_split_prot2$Reviewed)
table(iRef_split_prot2$Sequence)
length(unique(iRef_split$Uniprot_ID))
iRef_protein <- rio::import("iRefIndex_human_uniprot_list_20240608.tsv")