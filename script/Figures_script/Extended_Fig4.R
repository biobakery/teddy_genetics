library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringi)
library(stringr)
setwd("~/ddong/TEDDY_old//")
load("./output_2024/0129/Cluster800.Rdata")
bin <- cut(test$Days,breaks = seq(100,800,100),include.lowest = T)
test$bins <- bin
binnames <- names(table(bin))

setwd("~/ddong/TEDDY_old/output_2024/maaslin_species_0301_new_separate/")
load('./all_species.Rdata')
all_species.df <- all_species.df %>% filter(Type != "Cluster1 vs Cluster2")
sig <- all_species.df %>% filter(qval < 0.25)
alldata_select_species <- all_species.df %>% filter(feature %in% sig$feature)
alldata_select_species <-alldata_select_species %>% mutate(sig = case_when(qval < 0.001 ~ "***",
                                                                           qval < 0.05 ~ "**",
                                                                           qval < 0.25 ~ "*"))
alldata_select_species$bin <- factor(alldata_select_species$bin,levels = binnames)
alldata_select_species <- alldata_select_species %>% mutate(Species = gsub("s__","",feature))
pdf("~/ddong/TEDDY_project//output/Extended_Fig4",23,20)
ggplot(alldata_select_species,aes(x = bin,y = Species,fill = coef))+ geom_tile()+
  scale_fill_gradient2(low = "darkgreen",mid = "white",high = "#7209b7") +   theme_cowplot()  +
  geom_text(aes(label = sig),size = 10) +facet_wrap(~ Type) + theme(axis.text.y = element_text(size =20))
dev.off()






