source("~/script/packages_load.R")
load("~/ddong/TEDDY/output/microbiome_PC.Rdata")
load("~/ddong/TEDDY/output/tax_phylum_sample.Rdata")
PC <- fread("~/ddong/TEDDY/genetics.bed/filter/teddy_out8.qc.PCA_PC1-PC5.tsv")
PC <- PC %>% rename(subject_id = IID)
subject <- fread("~/ddong/TEDDY/TEDDY_metadata/metadata_subject.txt")
all <- PC %>% left_join(subject)
metadata <- fread("~/ddong/TEDDY_old/output_2024/metadata_sample.txt")
metadata_select <- metadata# %>% select(subject_id,sample_id,Days)
mds.data$sample_id <- as.numeric(mds.data$sample_id )
all <- metadata_select %>% inner_join(mds.data)
all <- all %>% inner_join(PC) %>% inner_join(subject)
all <- all %>% inner_join(result)
all <- all %>% inner_join(tax_phylum.df)
permanova <- fread("~/ddong/TEDDY/output/file.txt")

pdf("./output/microbiome_PCoA.pdf",5,5)
ggplot(data=all, aes(x=Microbiome_PCo1, y=Microbiome_PCo2)) +
  geom_point(size=2.5, aes(fill=Days), shape=21) +
  # scale_fill_manual(values = c("#1e6091","#e63946","#fcbf49","#6d597a"))+
  scale_fill_viridis_c(option="magma")+ 
  #scale_fill_viridis_c(direction=-1, limits=c(0,60), breaks=c(0,20,40,60)) +
  xlab(paste("Microbiome PCo1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("Microbiome PCo2 - ", mds.var.per[2], "%", sep="")) +
  #labs(title=spe_lab[i]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 14, face = "italic", hjust = 0.5),
        legend.position="bottom")

ggplot(data=all, aes(x=Microbiome_PCo1, y=Microbiome_PCo2)) +
  geom_point(size=2.5, aes(fill=Country), shape=21) +
  scale_fill_manual(values = c("#1e6091","#e63946","#fcbf49","#6d597a"))+
  #scale_fill_viridis_c(direction=-1, limits=c(0,60), breaks=c(0,20,40,60)) +
  xlab(paste("Microbiome PCo1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("Microbiome PCo2 - ", mds.var.per[2], "%", sep="")) +
  #labs(title=spe_lab[i]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 14, face = "italic", hjust = 0.5),
        legend.position="bottom")

ggplot(data=all, aes(x=Microbiome_PCo1, y=Microbiome_PCo2)) +
  geom_point(size=2.5, aes(fill=Bacteroidetes), shape=21,alpha = 0.8) +
  # scale_fill_manual(values = c("#1e6091","#e63946","#fcbf49","#6d597a"))+
  #scale_fill_viridis_c(option="magma")+ 
  #scale_fill_viridis_c(direction=-1, limits=c(0,60), breaks=c(0,20,40,60)) +
  xlab(paste("Microbiome PCo1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("Microbiome PCo2 - ", mds.var.per[2], "%", sep="")) +
  #labs(title=spe_lab[i]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 14, face = "italic", hjust = 0.5),
        legend.position="bottom")

ggplot(data=all, aes(x=Microbiome_PCo1, y=Microbiome_PCo2)) +
  geom_point(size=2.5, aes(fill=Firmicutes), shape=21,alpha = 0.8) +
  # scale_fill_manual(values = c("#1e6091","#e63946","#fcbf49","#6d597a"))+
  #scale_fill_viridis_c(option="magma")+ 
  #scale_fill_viridis_c(direction=-1, limits=c(0,60), breaks=c(0,20,40,60)) +
  xlab(paste("Microbiome PCo1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("Microbiome PCo2 - ", mds.var.per[2], "%", sep="")) +
  #labs(title=spe_lab[i]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 14, face = "italic", hjust = 0.5),
        legend.position="bottom")

ggplot(data=all, aes(x=Microbiome_PCo1, y=Microbiome_PCo2)) +
  geom_point(size=2.5, aes(fill=Actinobacteria), shape=21,alpha = 0.8) +
  # scale_fill_manual(values = c("#1e6091","#e63946","#fcbf49","#6d597a"))+
  #scale_fill_viridis_c(option="magma")+ 
  #scale_fill_viridis_c(direction=-1, limits=c(0,60), breaks=c(0,20,40,60)) +
  xlab(paste("Microbiome PCo1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("Microbiome PCo2 - ", mds.var.per[2], "%", sep="")) +
  #labs(title=spe_lab[i]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 14, face = "italic", hjust = 0.5),
        legend.position="bottom")

dev.off()
