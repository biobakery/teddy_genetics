setwd("~/ddong/TEDDY_project/")
library(tidyverse)
library(reshape2)
library(ggrepel)
library(scales)
library(dplyr)
library(data.table)


##--------------------------------------------------------------------------------
## Fig. 4c: Manhattan plot of the top 5 SNP loading factors on chromosomes 3, 6, and 8.
##--------------------------------------------------------------------------------

###Genetic PC loadings
loadings <- fread("./data/genetics.bed/filter/teddy.qc.excl_outliers.eigenvec.var", header=T) %>%
  select(2, 5:9) %>%
  rename(snp=VAR)

### SNP-to-gene mapping
map <- fread("./output/map.txt")
map <- map %>% select(V4,V8,V9)
map <- map %>% select(snp,hgnc_symbol,gene_id)

ml <- loadings %>%
  melt(id="snp", variable="pc", value.name="value")

ml_ann <- loadings %>%
  melt(id="snp", variable.name="pc", value.name="value") %>%
  mutate(v2 = value^2) %>%
  group_by(pc) %>%
  top_n(v2, n=10) %>%
  ungroup() %>%
  merge(., map, by="snp") %>%
  dplyr::select(pc, snp, hgnc_symbol) %>%
  merge(., ml, by=c("pc", "snp"))


# SNPs positions
allpos <- fread("./data/genetics.bed/filter/teddy_qc.map",
                header=F,sep="\t")
allpos <- allpos %>% filter(V1 %in% c(3,6,8))

man <- allpos %>%
  rename(chr=1, snp=2, morgans=3, bp=4) %>%
  arrange(chr, bp) %>%
  mutate(bp = cumsum(as.numeric(bp))) %>%
  merge(., ml, by="snp")

# set x-axis labels to the center of each chromosome
axis.set <- man %>% 
  group_by(chr) %>% 
  summarise(center = (max(bp) + min(bp)) / 2, .groups="keep")  %>% rename(chr = chr)

axis.breaks <- man %>% 
  group_by(chr) %>% 
  summarise(min=min(bp), max=max(bp), .groups="keep") %>% rename(chr = chr)

# the number of chromosomes
nCHR <- nlevels(as.factor(man$chr))

# annotation Top 10 SNPs
man_ann <- loadings %>%
  melt(id="snp", variable.name="pc", value.name="value") %>%
  mutate(v2 = value^2) %>%
  group_by(pc) %>%
  top_n(v2, n=10) %>%
  ungroup() %>%
  merge(., map, by="snp") %>%
  dplyr::select(pc, snp, hgnc_symbol) %>%
  merge(., man, by=c("pc", "snp")) %>% mutate(id = paste(pc,hgnc_symbol)) %>%
  filter(!duplicated(id)) %>%
  as.data.frame() 

##################
# Manhattan plot #
##################

pdf("./figures/Fig3c.pdf",5,4.5)
ggplot(man, aes(x=bp, y=value)) +
  facet_grid(pc~., scales="fixed") +
  geom_point(aes(colour=as.factor(chr), fill=as.factor(chr)),size =0.7,) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  scale_colour_manual(values = c("#457b9d","#bc4749","#e9c46a"))+
  theme_bw() +
  theme( 
    legend.title = element_text(size=12.5, face="bold", vjust=0.75),
    legend.text = element_text(size=8.75),
    legend.position = "top",
    panel.grid.minor.x = element_line(colour="black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5),
    axis.text.y = element_text(size=10),
    axis.title = element_text(size=12.5, face="bold"),
    strip.text = element_text(size=12.5, face="bold")) +
  scale_x_continuous(
    expand=c(0, 0),
    label = axis.set$chr, 
    breaks = axis.set$center,
    guide = guide_axis(check.overlap = TRUE),
    minor_breaks = axis.breaks$min) +
  labs(x="Chromosome", y="Loading", colour="Loading factors", fill="Loading factors") +
  geom_label_repel(data=man_ann, aes(label=hgnc_symbol),max.overlaps = 20,size = 3)+
  scale_fill_brewer(palette = "Set3")

dev.off()
