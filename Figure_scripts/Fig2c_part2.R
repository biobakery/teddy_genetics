library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringi)
library(stringr)
library(cowplot)
setwd("~/ddong/TEDDY_project/")

load('./output/traj_cluster.Rdata')
bin <- cut(test$Days,breaks = seq(100,800,100),include.lowest = T)
test$bins <- bin
binnames <- names(table(bin))

setwd("~/ddong/TEDDY/output/maaslin_species/")

##### load all MaAslin results
data <- list()
all_species.df <- data.frame()
for(i in 1:length(binnames)){
  b_n <- binnames[i]
  data[[i]] <- fread(paste("./", paste("Cluster_together",b_n,sep = "_"),"/all_results.tsv",sep = "") )
  if(i <= 3){
    data[[i]]  <- data[[i]]  %>% filter(metadata == "Cluster_together")  %>%
      mutate(Type = paste(value," vs Cluster1",sep =""),
             bin = binnames[i])
  }else{
    data[[i]]  <- data[[i]]  %>% filter(metadata == "Cluster_together")  %>%
      mutate(Type = paste(value," vs Cluster1&2",sep =""),
             bin = binnames[i])
  }
  all_species.df <- rbind(all_species.df,data[[i]])
}
data2 <-list()
table(all_species.df$Type)

sig <- all_species.df %>% filter(qval < 0.25)
alldata_select_species <- all_species.df %>% filter(feature %in% sig$feature)
alldata_select_species <-alldata_select_species %>% mutate(sig = case_when(qval < 0.001 ~ "***",
                                                                           qval < 0.05 ~ "**",
                                                                           qval < 0.25 ~ "*"))
alldata_select_species$bin <- factor(alldata_select_species$bin,levels = binnames)
alldata_select_species <- alldata_select_species %>% mutate(Species = gsub("s__","",feature))


species_curated <- c("Dorea_longicatena","Bifidobacterium_breve","Bifidobacterium_bifidum",
                     "Lactobacillus_rhamnosus","Oscillibacter_unclassified","Coprobacillus_unclassified",
                     "Clostridium_symbiosum","Clostridium_nexile","Clostridium_hathewayi",
                     "Clostridium_bolteae","Bacteroides_fragilis","Ruminococcus_gnavus",
                     "Bacteroides_thetaiotaomicron","Bacteroides_ovatus","Bifidobacterium_dentium",
                     "Bifidobacterium_catenulatum")



sig2 <- all_species.df %>% filter(qval < 0.25) %>% filter(abs(coef)> 0.5)
first <- sig2 %>% filter(Type == "Cluster2_3 vs Cluster1")
sp1 <- names(sort(-table(first$feature)))[1:10]
second <- sig2 %>% filter(Type == "Cluster3 vs Cluster1&2")
sp2<- table(second$feature)

sp <- table(sig2$feature)
sp <- names(sp[sp > 2])
sp2 <- names(sp2[sp2 > 1])
sp <- unique(c(sp,sp2))
sp <- c(sp,"s__Bifidobacterium_longum")
sp <- gsub("s__","",sp)


alldata_select_species2 <- alldata_select_species %>% filter(Species %in% sp)
alldata_select_species2$Type2 <-factor(alldata_select_species2$Type,levels = c("Cluster3 vs Cluster1&2","Cluster2_3 vs Cluster1"))
x <- alldata_select_species2 %>% arrange((Type2),desc(bin),-coef) 
alldata_select_species2$Species <- factor(alldata_select_species2$Species,levels = rev(unique(x$Species)))
alldata_select_species2$rank2 <- as.numeric(alldata_select_species2$Species)

alldata_select_species2 <- alldata_select_species2 %>% group_by(bin,Type) %>%
  mutate(rank = rank(-coef))


color = c("#2c6e49","#023047","#219ebc",
          "#ee964b","#f4d35e","#b6ad90","#5e548e","#a53860","#da627d")
split_species <- data.frame(str_split_fixed(alldata_select_species2$Species,"_",2))
split_species <- split_species %>% mutate(genus = str_sub(X1,1,1))
alldalldata_select_species2$Species <- factor(alldata_select_species2$Species,
                                          levels = names(table(as.character(alldata_select_species2$Species))))ata_select_species2$Species_anno <- paste(split_species$genus,split_species$X2,sep = ".")
alldata_select_species2$Species <- factor(alldata_select_species2$Species,
                                          levels = names(table(as.character(alldata_select_species2$Species))))
pdf("~/ddong/TEDDY_project//figures//Fig2b_part2",14,6)
ggplot(alldata_select_species2,aes(x = coef,y = as.factor(rank),fill = Species))+geom_bar(stat = "identity")  +
  facet_wrap(Type~bin,nrow = 2,scales = "free_y") + theme_bw()  +
  geom_text(aes(label = (sig)),size =6,color ="darkred")+
  scale_fill_manual(values = color)+theme_cowplot()+ 
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank())+ ylab("Species")
dev.off()
