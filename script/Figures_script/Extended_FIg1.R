source("~/script/packages_load.R")
setwd("~/ddong/TEDDY/")
load("./output/perm.Rdata")

source("~/script/packages_load.R")
load("~/ddong/TEDDY/output/microbiome_PC.Rdata")
PC <- fread("~/ddong/TEDDY/genetics.bed/filter/teddy_out8.qc.PCA_PC1-PC5.tsv")
PC <- PC %>% rename(subject_id = IID)
subject <- fread("~/ddong/TEDDY/TEDDY_metadata/metadata_subject.txt")
all <- PC %>% left_join(subject)
metadata <- fread("~/ddong/TEDDY_old/output_2024/metadata_sample.txt")
sex <- metadata %>% select(subject_id,Sex,time_to_brstfed_stop) %>% unique()
metadata_select <- metadata %>% select(subject_id,sample_id,Days)
mds.data$sample_id <- as.numeric(mds.data$sample_id )
all <- metadata_select %>% inner_join(mds.data)
all <- all %>% inner_join(PC) %>% inner_join(subject)
all$Country <- as.factor(all$Country)


species_rank <- data.frame(aveRA = colMeans(df_sp))
species_rank$species <- rownames(species_rank)
species_rank <- species_rank %>% arrange(desc(aveRA))


top10 <- df_sp %>% select(one_of(species_rank$species[1:10]))
top10$sample_id <- rownames(top10)
top10.df <- melt(top10) %>% rename(species = variable,RA=value)
top10.df$sample_id <- as.numeric(top10.df$sample_id)
all_select <- all %>% select(Days,Country,sample_id)
top10.df$RA[top10.df$RA == 0] <- 6e-07
top10.df <- top10.df %>% inner_join(all_select)
top10.df <- top10.df %>% left_join(b1)
top10.df <- top10.df %>% arrange(Country,Days)
top10.df$species <- factor(top10.df$species,levels = rev(species_rank$species[1:10]))
top10.df$sample_id <- factor(top10.df$sample_id,levels = unique(top10.df$sample_id))


plot_species <- ggplot(top10.df, aes(sample_id,species)) +
  geom_tile(aes(fill = log10(RA*100))) +
  geom_xsidetile(aes(y = "Country", xfill = Country))+
  scale_xfill_manual(values =c("#1e6091","#e63946","#fcbf49","#6d597a"))+
  #  geom_xsidetile(aes(y = "Days", xfill = Days))+
  scale_fill_viridis_c(name ="Relative abundance\nlog10 scale", na.value = "white",option = "rocket") +
  theme(axis.text.x = element_blank(),  # 隐藏 x 轴文字
        axis.ticks.x = element_blank()) 


EC <- fread("~/ddong/TEDDY/output_2024/ec_all_qc_1e6.txt")
ec_name <- fread("~/ddong/TEDDY/output_2024/ec_name.txt")
EC_rank <- data.frame(aveRA = colMeans(EC[,-1]))
EC_rank$EC_feature <- rownames(EC_rank)                
ec_name <- ec_name %>% mutate(EC_feature = paste("X",EC,sep = ""))
EC_rank <- EC_rank %>% left_join(ec_name)
EC_rank <- EC_rank %>% arrange(desc(aveRA))
EC_rank <- EC_rank %>% filter(!EC_feature %in% c("UNGROUPED","UNMAPPED"))


top10EC <- EC %>% select(sample_id,one_of(EC_rank$EC_feature[1:10]))
top10EC$sample_id <- as.character(top10EC$sample_id)
top10EC.df <- melt(top10EC) %>% rename(EC_feature = variable,ECRA=value)
top10EC.df$sample_id <- as.numeric(top10EC.df$sample_id)
all_select <- all %>% select(Days,Country,sample_id)
top10EC.df <- top10EC.df %>% inner_join(all_select)
top10EC.df <- top10EC.df %>% filter(sample_id %in% top10.df$sample_id)
top10EC.df <- top10EC.df %>% inner_join(ec_name)
top10EC.df$sample_id <- factor(top10EC.df$sample_id,levels = unique(top10.df$sample_id))
top10EC.df$EC_name <- factor(top10EC.df$EC_name,levels = unique(top10EC.df$EC_name))
EC_rank <- EC_rank %>% mutate(ID =paste(EC,EC_name,sep = "_"))
top10EC.df <- top10EC.df %>% mutate(ID =paste(EC,EC_name,sep = "_"))
level <- rev(EC_rank$ID[1:10])
top10EC.df$ID <- factor(top10EC.df$ID,levels = rev(EC_rank$ID[1:10]))

plotEC <- ggplot(top10EC.df, aes(sample_id,ID)) +
  geom_tile(aes(fill = log10(ECRA*100))) +
  geom_xsidetile(aes(y = "Days", xfill = Days))+
  scale_fill_viridis_c(name ="Relative abundance\nlog10 scale", na.value = "white", option = "mako")+
  scale_xfill_gradient(low ="#dedbd2", high = "#4a5759") +
  theme(axis.text.x = element_blank(),  # 隐藏 x 轴文字
        axis.ticks.x = element_blank()) 


pdf("./output/top10_distirbution.pdf",14,6)
plot_grid(plot_species, plotEC, ncol = 1,align = 1)
dev.off()
