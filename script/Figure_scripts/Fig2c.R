source("~/script/packages_load.R")
setwd("~/ddong/TEDDY/")

##loading data（Clusteing results, microbiome, metadata）
load("./output/output/Cluster800.Rdata") 
major <- fread('~/ddong/TEDDY/data/data_derived/metaphlan2_major.tsv',header = T) 

metadata <- fread("./output/metadata_sample.txt") 

major <- major %>% column_to_rownames("ID")
major <- major %>% t(.) %>% as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  rename(sample_id=1)
major_old <- major
major[,-1] <-  sweep(major[,-1], 1, STATS = rowSums(major[,-1]), FUN = "/") *100
map.df <- metadata[,c("subject_id","sample_id","Days")]
map.df$sample_id <- as.character(map.df$sample_id)

major <- major %>% inner_join(map.df)
major$sample_id <- as.character(major$sample_id )
major$subject_id <- as.character(major$subject_id )
major$Days <- as.character(major$Days )
major.df <- melt(major)
major.df$Days <- as.numeric(major.df$Days)
major.df<- major.df %>% filter(Days < 800)


Cluster.df <- Cluster800 %>% select(subject_id,Cluster)
Cluster.df$subject_id <- as.character(Cluster.df$subject_id)
major.df <- major.df %>% left_join(Cluster.df)
major.df <- major.df %>% filter(!is.na(Cluster))



bc <- fread("./output/bc_subject_sample.txt")
bc$sample_id <- as.character(bc$sample_id )
bc$subject_id <- as.character(bc$subject_id )
major.df <- major.df %>% inner_join(Cluster.df)
bin <-  cut(as.numeric(major.df$Days),breaks = c(seq(100,800,100)),include.lowest = T)
major.df$bins <- bin
major.df <- major.df %>% filter(!is.na(bins))

major_summary.df <- major.df %>% group_by(bins, variable,Cluster) %>% summarise(mean = mean(value),
                                                                                median = median(value))



select_species <- c("Bifidobacterium bifidum","Bifidobacterium breve","Bifidobacterium longum",
                    "Ruminococcus gnavus","Ruminococcus torques", "Ruminococcus bromii","Escherichia coli",
                    "Faecalibacterium prausnitzii" ,"Klebsiella oxytoca","Klebsiella unclassified" ,
                    "Coprobacillus unclassified")
select_species <- c("Bacteroides_thetaiotaomicron", "Bifidobacterium_bifidum","Bifidobacterium_breve"  ,     
                    "Clostridium_hathewayi","Clostridium_symbiosum","Coprobacillus_unclassified"  ,
                    "Ruminococcus_gnavus","Ruminococcus_torques"    ,     "Dorea_longicatena" ,       
                    "Bifidobacterium_longum")
select_species <- gsub("_"," ",select_species)
select_species <-names(table(select_species))
major_summary.df$species <- gsub("s__","",major_summary.df$variable)
major_summary.df$species <- gsub("_"," ",major_summary.df$species)
major_summary.df <- major_summary.df %>% mutate(color =  case_when(species %in% select_species ~ species,
                                                                   !species %in% select_species ~ "others"))
major_summary.df$color <- factor(major_summary.df$color,levels = (c(select_species,"others")))
color = (c("#2c6e49","#023047","#219ebc","#8ecae6",
           "#ee964b","#f4d35e","#b6ad90","#5e548e","#a53860","#da627d","#ffa5ab",
           "grey"))

pdf("./output/Fig2B_barplot.pdf",10,23)
major_summary_noother.df <- major_summary.df %>% filter(color != "others")
ggplot(major_summary_noother.df,aes(x = bin_numeric,y = mean,fill = color))+
  geom_bar(stat = "identity") +
  geom_smooth(data = major_select.df,aes(x = Day,y = value,fill= color),method = "loess",se = FALSE)+
  facet_grid(color~Cluster,scales = "free")+scale_fill_manual(values = color)+theme_bw()+
  
  #scale_x_discrete(labels = c(2:8)) +
  xlab("100 Days bin")+ theme_cowplot()+theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )+
  ylab("Mean Relative aboundance with Cluster") + geom_point() 
dev.off()



