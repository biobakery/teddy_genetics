setwd("~/ddong/TEDDY_project/")

library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringi)
library(stringr)
library(Maaslin2)
library(cowplot )

##---------------------------------------------------
## Analysis: MaAsLin
##---------------------------------------------------

### Data loading: Microbiome
microbiome <- fread('./output/metaphlan2_afterQC.tsv',header = T)
microbiome <- microbiome %>% column_to_rownames("Species")
microbiome <- microbiome %>% t(.) %>% as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  #setNames(paste0("microbiome_", names(.))) %>%
  rename(sample_id=1)
microbiome[,-1] <-  sweep(microbiome[,-1], 1, STATS = rowSums(microbiome[,-1]), FUN = "/") *100
species <- colnames(microbiome)[-1]

### Data loading: Clustering results.
load("./output/traj_cluster_results.Rdata")
bin <- cut(metadata_cluster$Days,breaks = seq(100,800,100),include.lowest = T)
metadata_cluster$bins <- bin
binnames <- names(table(bin))


### MaAsLin association analyses 
setwd("./output/maaslin_species")
for ( i in 1:length(binnames)){
  b_n <- binnames[i]
  samples <- metadata_cluster %>% filter(bins == b_n)
  sample_select <- microbiome %>% filter(sample_id %in% samples$sample_id)
  
  meta <- samples[match(sample_select$sample_id,samples$sample_id),]
  
  ### Compare catch-up and incomplete vs. healthy in the first 400 days, 
  ### and compare incomplete vs. catch-up and healthy in the later 400 days
  if(i <= 3){
    meta  <- meta %>% mutate(Cluster_together = case_when(Cluster == 1 ~ "Cluster_base",
                                                          Cluster == 2 ~ "Cluster2_3",
                                                          Cluster == 3 ~ "Cluster2_3" ))
  } else{
    meta  <- meta %>% mutate(Cluster_together = case_when(Cluster == 1 ~ "Cluster_base",
                                                          Cluster == 2 ~ "Cluster_base",
                                                          Cluster == 3 ~ "Cluster3" ))
  }
  sample_select <- column_to_rownames(sample_select,"sample_id")
  meta <- data.frame(meta) 
  rownames(meta) <- meta$sample_id
  cluster_masslin <- Maaslin2(sample_select,
                              meta,
                              paste("Cluster_together",b_n,sep = "_"),
                              min_abundance=0.0001,
                              min_prevalence=0.1,
                              fixed_effects = c("Cluster_together","Clinical_Center","Sex",
                                                "delivery_class", "abx_this_month",
                                                "solid_food_now","current_brst_fed",),
                              reference      = c("Country,Sweden",
                                                 "Cluster_together,Cluster_base"),
                              random_effects = c("subject_id"))
}

### Combine all MaAsLin results 
data <- list()
all_results.df <- data.frame()
for(i in 1:length(binnames)){
  b_n <- binnames[i]
  data[[i]] <- fread(paste("./", paste("Cluster_together",b_n,sep = "_"),"/all_results.tsv",sep = "") )
  data[[i]]  <- data[[i]]  %>% filter(metadata == "Cluster_together")  %>% 
    mutate(bin = binnames[i])
  all_results.df <- rbind(all_results.df,data[[i]])
}

all_results.df$Type <- all_results.df$value
save.image(file = "./output/MaAslin_species.Rdata")


##---------------------------------------------------
## Fig. 2c: Bar plot of MaAsLin results
##---------------------------------------------------

load("./output/MaAslin_species.Rdata")
setwd("~/ddong/TEDDY_project/")
### Filter Out Significant Results
sig <- all_results.df %>% filter(qval < 0.25)
alldata_select_species <- all_results.df %>% filter(feature %in% sig$feature)

select_species <- names(sort(-table(sig$feature)))[1:8] #add s__Dorea_longicatena
select_species <- gsub("s__","",select_species)
alldata_select_species$Species <- gsub("s__","",alldata_select_species$feature)
all_species.df <- alldata_select_species %>% filter(Species %in% select_species)
all_species.df$bin <- factor(all_species.df$bin,levels = binnames)


all_species.df <- all_species.df %>%
  group_by(bin) %>%
  mutate(rank = rank(-coef)) 

color = (c("#2c6e49","#023047","#219ebc",
           "#ee964b","#f4d35e","#b6ad90",
           "#5e548e","#a53860","#da627d"))

pdf("~/ddong/TEDDY_project//figures//Fig2b_Masslin.pdf",20,3)
ggplot(all_species.df, aes(y = as.factor(rank), x= coef, fill = Species)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ bin,scale = "free_y",nrow = 1) +  # Each facet has its own x order
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() +
  scale_fill_manual(values = color)+
  theme(
    panel.grid = element_blank(),          # Remove all grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Keep border
    axis.title.y = element_blank(),        # Remove y-axis title
    axis.text.y = element_blank(),         # Remove y-axis labels
    axis.ticks.y = element_blank()         # Remove y-axis ticks
  )
dev.off()




##-----------------------------------------------------
## Fig. 2c: Bar Plot of Species Relative Abundance
##-----------------------------------------------------

### Data Loading: Cluster Results, Microbiome, and Metadata
load("./output/traj_cluster_results.Rdata") 
microbiome <- fread('./output/metaphlan2_afterQC.tsv',header = T) 
load('./output/preprpreprocess/metadata_all.Rdata')

microbiome <- microbiome %>% column_to_rownames("ID")
microbiome <- microbiome %>% t(.) %>% as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  rename(sample_id=1)
microbiome[,-1] <- microbiome[,-1]*100


### Adding Time Information for Samples
map.df <- metadata[,c("subject_id","sample_id","mgx_age")]
map.df$sample_id <- as.character(map.df$sample_id)
microbiome <- microbiome %>% inner_join(map.df)
microbiome$sample_id <- as.character(microbiome$sample_id )
microbiome$subject_id <- as.character(microbiome$subject_id )
microbiome$Days <- as.character(microbiome$mgx_age )
microbiome$mgx_age <- NULL
microbiome.df <- melt(microbiome)
microbiome.df$Days <- as.numeric(microbiome.df$Days)
microbiome.df<- microbiome.df %>% filter(Days < 800)

### Combine with clustering results
Cluster.df <- Cluster_results %>% select(subject_id,Cluster)
Cluster.df$subject_id <- as.character(Cluster.df$subject_id)
microbiome.df <- microbiome.df %>% inner_join(Cluster.df)

### Stratified by time into bins
bin <-  cut(as.numeric(microbiome.df$Days),breaks = c(seq(100,800,100)),include.lowest = T)
microbiome.df$bins <- bin
microbiome.df <- microbiome.df %>% filter(!is.na(bins))

microbiome_summary.df <- microbiome.df %>% group_by(bins, variable,Cluster) %>% 
  summarise(mean = mean(value),median = median(value))

### Select Four Highlighted Species
select_species <- c("Bifidobacterium_bifidum","Ruminococcus_gnavus",    
                    "Clostridium_hathewayi", "Dorea_longicatena" )
select_species <- gsub("_"," ",select_species)
color = c("#023047","#a53860","#ee964b","#5e548e")
names(color) <- select_species

microbiome_summary.df$species <- gsub("s__","",microbiome_summary.df$variable)
microbiome_summary.df$species <- gsub("_"," ",microbiome_summary.df$species)
microbiome_summary_selected.df <- microbiome_summary.df %>% 
  filter(species %in% select_species)

pdf("./figures/Fig2b_barplot.pdf",12,15)
ggplot(microbiome_summary_selected.df,aes(x = bins,y = mean*100,fill = color))+
  geom_bar(stat = "identity") +
  geom_smooth(data = microbiome_summary_selected.df,
              aes(x = as.numeric(bins),y = mean*100),method = "loess",se = FALSE)+
  facet_grid(color~Cluster,scales = "free")+scale_fill_manual(values = color)+
  xlab("Time (Days)")+ theme_cowplot()+theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1))+
  ylab("Mean Relative aboundance") + geom_point() 
dev.off()

save.image("./output/Fig2c.Rdata")
