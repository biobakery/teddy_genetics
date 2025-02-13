source("~/script/packages_load.R")
load("~/ddong/TEDDY_project/output/microbiome_PC.Rdata")
load("~/ddong/TEDDY_project/output/tax_phylum_sample.Rdata")
PC <- fread("~/ddong/TEDDY_project/genetics.bed/filter/teddy_out8.qc.PCA_PC1-PC5.tsv")
PC <- PC %>% rename(subject_id = IID)
subject <- fread("~/ddong/TEDDY_project/TEDDY_metadata/metadata_subject.txt")
all <- PC %>% left_join(subject)
metadata <- fread("~/ddong/TEDDY_project_old/output_2024/metadata_sample.txt")
metadata_select <- metadata# %>% select(subject_id,sample_id,Days)
mds.data$sample_id <- as.numeric(mds.data$sample_id )
all <- metadata_select %>% inner_join(mds.data)
all <- all %>% inner_join(PC) %>% inner_join(subject)
all <- all %>% inner_join(result)
all <- all %>% inner_join(tax_phylum.df)
permanova <- fread("~/ddong/TEDDY_project/output/file.txt")

pdf("/figures/Fig1c",5,3.5)
ggplot(all,aes(x= PC1,y = PC2,color = Country)) +geom_point() + 
  theme_cowplot()+
  #theme_bw()+
  geom_ysidedensity(aes(fill = Country),alpha = 0.5) +
  geom_xsidedensity(aes(fill = Country),alpha = 0.5) +
  scale_fill_manual(values = c("#1e6091","#e63946","#fcbf49","#6d597a"))+
  scale_color_manual(values = c("#1e6091","#e63946","#fcbf49","#6d597a"))
dev.off()

pdf("./figures/Fig1e",5,4)
ggplot(all,aes(x= PC1,y = Microbiome_PCo1,color = Days)) +geom_point(size= 1) + 
  theme_cowplot()+xlab("Genetic PC1") +ylab("micorbiome PCo1")+ scale_color_viridis(option="magma")
dev.off()





#####Figure 1f
all <- fread("~/ddong/TEDDY/output/Permanova.txt")
all <- all %>% filter(Covariate != "delivery")
b <- all %>% filter(type == "multi")
a <- all %>% filter(type == "Uni")
a <- a %>% arrange(R2)
all$Covariate <- factor(all$Covariate,levels = a$Covariate)
plot1 <- ggplot(all,aes(x = R2*100, y = Covariate,fill = type)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_text(data = b,aes(x = R2*100, y = Covariate,label = sig), vjust = 1.2)+
  geom_text(data = a,aes(x = R2*100, y = Covariate,label = sig), vjust = 0.4) +
  theme_cowplot() + scale_fill_manual(values =c( "#457b9d","#bc4749"))

all2 <- all 
all2$R2[all2$R2 > 0.02] <- 0.02
b <- all2 %>% filter(type == "multi")
a <- all2 %>% filter(type == "Uni")
plot2 <- ggplot(all2,aes(x = R2*100, y = Covariate,fill = type)) +
  geom_bar(stat = "identity",position = "dodge") +
  # geom_text(data = b,aes(x = R2*100, y = Covariate,label = sig), vjust = 1.3)+
  #geom_text(data = a,aes(x = R2*100, y = Covariate,label = sig), vjust = 0.3) +
  theme_cowplot() + scale_fill_manual(values =c( "#457b9d","#bc4749"))

pdf("./figures/Fig1f.pdf",12,8)
plot_grid(plot1,plot2,nrow= 1)
dev.off()
