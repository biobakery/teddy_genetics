library(dplyr)
library(data.table)
library(stringr)
library(ggside)
setwd("~/ddong/TEDDY_old//")

####load data
load("./output/infant_species_all_nobrst.Rdata")
load("./output/infant_species_all_brst.Rdata")
load("./output/toddler_species_all.Rdata")


df_sig_infant_brst$Class <- "Infant Breastfeeding"
df_sig_infant_nobrst$Class <- "Infant Post-breastfeeding"
df_sig_Toddler <- df_sig_Toddler %>% group_by(Type) %>% 
  mutate(P_adj_bon = p.adjust(P,method = "bonferroni"),
         P_adj  = p.adjust(P,method = "fdr"))

##combine infant results
col <- intersect(colnames(df_sig_Toddler),colnames(df_sig_infant_brst))
all <- rbind(df_sig_infant_brst[,col],df_sig_infant_nobrst[,col],df_sig_Toddler[,col])
all <- all %>% mutate(star = case_when(P_adj_bon < 0.001 ~ "***",
                                       P_adj_bon < 0.01 ~ "**",
                                       P_adj_bon < 0.05 ~ "*"))
all <- all %>% mutate(species = str_extract(Model,"(?<=lmer\\().*(?=\\ \\~)"))
all$species <- gsub("`","",all$species)

##load tax pologentic information
tax <- fread("./output/species_major_tree.txt")
all <- all %>% left_join(tax,by = c("species" = "species"))

###keep the significant results
all <- all %>% mutate(paste= paste(Predictor,Class))
sig <- all %>%filter(P_adj_bon< 0.05) 
Species <- all %>%group_by(Type) %>%  filter(paste %in% sig$paste &
                                               MicroFeat %in% sig$MicroFeat)
Species$species <- gsub("s__","",Species$species)
Species$species <- gsub("_"," ",Species$species)
Species <- Species %>% arrange(phylum)
Species$classtype <- factor(Species$Class ,levels =  c("Infant Breastfeeding",
                                                       "Infant Post-breastfeeding",
                                                       "toddler" ))
Species$species <- factor(Species$species,levels = unique(Species$species))

pdf("./figures/Fig5A.pdf",width = 4.5,6)
ggplot(data=Species, aes(y=species, x=Predictor, fill=value)) +
  geom_tile(colour="black") +
  facet_grid(~classtype, space="free", scales="free") +
  scale_fill_gradient2(low = "#023047", 
                       mid = "white", 
                       high = "#9e2a2b") +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_discrete(expand=c(0, 0)) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.4),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.4),
        axis.text.y = element_text(size = 8,color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8, colour = "black"),
        axis.ticks.length=unit(0.04,"inch"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6),
        legend.position = "bottom",
        legend.text = element_text(colour = "black", size = 8),
        legend.title = element_text(colour = "black", size = 8)) +
  guides(fill = guide_colourbar(title.position="top",
                                title.hjust = 0.5,
                                barheight = 1.5,
                                barwidth = 10))+
  geom_text(aes(label=star), color="black", size=3, nudge_y = -0.2)  +
  geom_ysidetile(aes(x = "phylum", yfill = phylum)) +
  labs(fill="sign(coef) * -log10(adj_pvalue)", x="", y="")
dev.off()
fwrite(all,file ="./output/linearmix_infant_toddler.txt",sep = '\t',quote = F,col.names = T)
