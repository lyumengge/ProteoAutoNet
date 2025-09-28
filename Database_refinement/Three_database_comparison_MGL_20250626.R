##################################################
## Project: ProteoAutoNet
## Script purpose: Three databases comparison
## Date: 2025-06-06 
## Author: MenggeLYU
## Version: 1.0
##################################################



# Path and library --------------------------------------------------------
setwd("E://MenggeLYU//CHINDA")
library(rio)
library(dplyr)
library(tidyr)
library(stringr)
library(corrplot)

# corum -------------------------------------------------------------------
newcorum <- rio::import("E://MenggeLYU//CHINDA//CORUM//CORUM download 2022_09_12.xlsx") #5204
newcorum <- newcorum[which(newcorum$Organism=="Human"),] #3637
newcorum_prot <- as.data.frame(newcorum$`subunits(UniProt IDs)`)
colnames(newcorum_prot)[1] <- "Prot"
newcorum_prot_split <- newcorum_prot %>% as_tibble() %>% 
  separate_rows(Prot, sep = ";")
newcorum_prot_split <- unique(newcorum_prot_split)

uniprot <- rio::import("uniprot_protein_list_for_CHINDA_MGL_20240608.tsv")
newcorum_prot_left <- left_join(newcorum_prot_split, uniprot, by = "Prot")
is.na(newcorum_prot_left$Reviewed)
null <- newcorum_prot_left[which(is.na(newcorum_prot_left$Reviewed)),]
rio::export(null, "CORUM_unmatched_proteins_MGL_20240923.tsv")

# humap -------------------------------------------------------------------
humap <- rio::import("E://MenggeLYU//CHINDA//Hu.map//humap2_complexes_20200809.txt")
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
cp <- rio::import("E:\\MenggeLYU\\CHINDA\\Complex_portal\\complex_portal_9606.tsv")
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


# Venn plot ---------------------------------------------------------------
library(VennDiagram)
library(RColorBrewer)
humap_prot_pickone <- as.data.frame(humap_prot[,1])
colnames(humap_prot_pickone)[1] <- "Prot"
cp_prot_pickone <- as.data.frame(cp_unique_prot_left2[,1])
colnames(cp_prot_pickone)[1] <- "Prot"

veen1 <- venn.diagram(list(CORUM = t(newcorum_prot_split),
                           hu.MAP = t(humap_prot_pickone),
                           Complex_Portal = t(cp_prot_pickone)),
                      fill=c(brewer.pal(7,"Set1")[1:3]),
                      alpha = 0.5,cex = 2,
                      cat.cex = 2,cat.fontface = "bold",lty=2,cat.pos = 0,
                      resolution = 300,filename = NULL)
pdf('Figure1_Veen.pdf', width = 10, height = 10)
grid.draw(veen1)
dev.off()



