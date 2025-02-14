library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggExtra)

##-----------------------------------------
## Fig. 1d: Genetic PCA plot
##-----------------------------------------

### load metadata
load("./output/preprpreprocess/metadata_all.Rdata")

### load Genetic PCs 
genetic_PC <- fread("./data/genetics.bed/filter/teddy_out8.qc.PCA_PC1-PC5.tsv")
genetic_PC <- genetic_PC %>% rename(subject_id = IID)
subject <- fread("./data/TEDDY_metadata/metadata_subject.txt")
all <- genetic_PC %>% left_join(subject)


pdf("./figures/Fig1d_geneticPCA.pdf",7,6)
# Scatter plot
p <- ggplot(all, aes(x = PC1, y = PC2, color = Country)) +
  geom_point(alpha = 1, size = 2) +  
  theme_minimal() +
  labs(x = "Genetic PC1", y = "Genetic PC2") +
  theme(legend.position = "right") +
  scale_color_manual(values = c("#1e6091","#e63946","#fcbf49","#6d597a")) +
  theme_cowplot()

# Add density plots 
ggMarginal(p, type = "density", margins = "both", groupColour = TRUE, groupFill = TRUE)
dev.off()

##-----------------------------------------
## Fig. 1e: Genetic PC1 and Micorbiome PCo1
##-----------------------------------------

### load microbiome PCos
load("./output/microbiome_PC.Rdata")
metadata_select <- metadata %>% select(subject_id,sample_id,mgx_age)
microbiomePCo <- microbiomePCo %>% inner_join(metadata_select)
##Combined Microbiome PCoA and Genetic PC Data
all <- all %>% left_join(microbiomePCo)

pdf("./figures/Fig1e_microbiome_gentic_pc.pdf",5,4)
ggplot(all,aes(x= PC1,y = Microbiome_PCo1,color = mgx_age)) +geom_point(size= 1) + 
  theme_cowplot()+xlab("Genetic PC1") +ylab("micorbiome PCo1")+ 
  scale_color_viridis(option="magma")+labs(color = "Days") 
dev.off()


##-----------------------------------------
## Fig. 1f: PERMANOVA results
##-----------------------------------------

Permanova <- fread("./output//Permanova.txt")
##Select univariate results of PERMANOVA
uni <- Permanova %>% filter(type == "Uni")
uni <- uni %>% arrange(R2)
uni$Covariate <- factor(uni$Covariate,levels = uni$Covariate)

pdf("./figures/Fig1f_permanova.pdf",6,8)
ggplot(uni,aes(x = R2*100, y = Covariate,fill = type)) +
  geom_bar(stat = "identity",position = "dodge") +
  #geom_text(data = b,aes(x = R2*100, y = Covariate,label = sig), vjust = 1.2)+
  geom_text(data = uni,aes(x = R2*100, y = Covariate,label = sig), vjust = 0.4) +
  theme_cowplot() + scale_fill_manual(values =c( "#457b9d","#bc4749"))+
  xlab("% variance explained")
dev.off()


