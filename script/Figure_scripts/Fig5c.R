source("~/script/packages_load.R")
library(ggrepel)
setwd("~/ddong/TEDDY_project//")

###load SNP * Time Species level results
all <-fread("./output_new/SNP_linear_Species_toddler.txt")
all  <- all %>% mutate(rsid = case_when( grepl(":",Predictor) ~ str_extract(pattern = ".*(?=\\:)",Predictor),
                                         !grepl(":",Predictor) ~ Predictor))
position <- read.csv("./genetics.bed/filter/teddy-out8_map.map",header=F,sep="\t") %>%
  rename(chr=1, snp=2, morgans=3, bp=4)
sig <- all %>% filter(P <  5e-08)
select_species <- all %>% filter(Species %in% sig$Species)

####map to gene
map <- fread("./map.txt")
map <- map %>% select(V4,V8,V9)
map <- map %>% rename(snp = V4, symbol = V8, gene_id = V9)
map <- map %>% filter(!duplicated(snp))


####
fix <- select_species %>% filter(grepl(":",Predictor))
fix <- fix %>% rename(snp = rsid)

fix <- fix %>% mutate(log10P = -log10(P))
sig <- fix %>% filter(P < 5e-8)
species <- names(table(sig$Species))
fix_sepcies <- fix %>% filter(Species %in% species)
man <- read.csv("./genetics.bed/filter/teddy-out8_map.map",
                header=F,
                sep="\t") %>%
  rename(chr=1, snp=2, morgans=3, bp=4) %>%
  arrange(chr, bp) %>%
  mutate(bp = cumsum(as.numeric(bp))) %>%
  merge(., fix_sepcies, by="snp")

tax <- fread("./output_1028/species_major_tree.txt")
tax_select <- tax %>% filter(species %in% man$Species)
tax_select <- tax_select %>% select(phylum,genus,species)
man <- man %>% left_join(tax_select,by = c("Species" = "species"))

mostsig <- man %>% group_by(Species) %>% mutate(minP=min(P,na.rm = T) )
mostsig <- mostsig %>% filter(P == minP)

sig_man <- mostsig %>% left_join(map)
sig_man <- sig_man %>% mutate(anno = paste(Species,symbol,sep = "/"))
sig_man$anno <- gsub("s__","",sig_man$anno)
sig_man$anno <- gsub("_"," ",sig_man$anno)
sig_man <- sig_man %>% filter(!duplicated(Species))
man <- man %>% filter(chr != 25)
axis.set <- man %>% 
  group_by(chr) %>% 
  summarise(center = (max(bp) + min(bp)) / 2, .groups="keep") 

axis.breaks <- man %>% 
  group_by(chr) %>% 
  summarise(min=min(bp), max=max(bp), .groups="keep") 

nCHR <- nlevels(as.factor(man$chr))
man  <- man %>% mutate(Color = case_when(chr %in% c(seq(1,22,2),25) ~ "a",
                                         chr %in% c(seq(2,22,2)) ~ "b")) 
motsig_pos <- mostsig %>% filter(Coefficient >= 0)
motsig_neg <- mostsig %>% filter(Coefficient <0)
g <- names(table(man$genus))
g <- gsub("g__","",g)
man <- man %>% mutate(Color_combine = case_when(genus  == "g__Bacteroides" & P < 5e-8 ~ "c",
                                                genus  == "g__Barnesiella" & P < 5e-8 ~ "d",
                                                genus  == "g__Bifidobacterium" & P < 5e-8 ~ "e",
                                                genus  == "g__Collinsella" & P < 5e-8 ~ "f",
                                                genus  == "g__Coprococcus" & P < 5e-8 ~ "g",
                                                genus  == "g__Dorea" & P < 5e-8 ~ "h",
                                                genus  == "g__Eubacterium" & P < 5e-8 ~ "i",
                                                genus  == "g__Faecalibacterium" & P < 5e-8 ~ "j",
                                                genus  == "g__Parabacteroides" & P < 5e-8 ~ "k",
                                                genus  == "g__Prevotella" & P < 5e-8 ~ "l",
                                                genus  == "g__Ruminococcus" & P < 5e-8 ~ "m",
                                                genus  == "g__Streptococcus" & P < 5e-8 ~ "n",
                                                genus  == "g__Subdoligranulum" & P < 5e-8 ~ "o"
))

man <- man %>% mutate(Color_all = case_when(is.na(Color_combine) ~ Color,
                                            !is.na(Color_combine) ~ Color_combine))


colors <- rainbow(13)
colors <- c(colors,"#00C4FF")

manhplot <- ggplot(man, aes(x=bp, y=log10P)) +
  geom_point(aes(colour=Color_all, fill=Color_all),alpha = 0.8,size = 1.2) +
  geom_hline(yintercept = -log10(5e-8), color = "grey20", linetype = "dashed") +
  #facet_wrap(~Species, scales="fixed", nrow=length(species)) +
  scale_color_manual(values = c("#ccc5b9","#8d99ae",colors),
                     breaks = c(names(table(man$Color_all))),
                     labels = c("a","b",g))+
  #guides(fill = "none",color = "none")+
  #scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred") +
  #scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred") +
  theme_bw() +ggtitle("Toddler_Fixed")+
  theme( 
    legend.title = element_text(size=12.5, face="bold", vjust=0.75),
    legend.text = element_text(size=8.75),
    legend.position = "top",
    panel.grid.minor.x = element_line(colour="black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
    axis.text.y = element_text(size=10),
    axis.title = element_text(size=12.5, face="bold"),
    strip.text = element_text(size=12.5, face="bold")) +
  scale_x_continuous(
    expand=c(0, 0),
    label = axis.set$chr, 
    breaks = axis.set$center,
    guide = guide_axis(check.overlap = TRUE),
    minor_breaks = axis.breaks$min)+
  labs(x="Chromosome", y="-log10Pvalue")#+

ggsave(plot=manhplot, paste0("./figures/fig5c.png"), dpi=300, height=7, width=13)

