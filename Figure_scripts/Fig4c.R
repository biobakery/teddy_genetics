source("~/ddong/packages_load.R")
setwd("~/ddong/TEDDY_project")
opts <- list()
opts$i <- "./genetics.bed/filter/teddy_out8.qc.excl_outliers.eigenvec.var"
opts$m <- "./genetics.bed/filter/teddy-out8_map.map"
opts$a <- "./map.txt"
library(tidyverse)
library(reshape2)
library(ggrepel)
library(scales)

loadings <- read.table(opts$i, header=T) %>%
  select(2, 5:9) %>%
  rename(snp=VAR)

map <- fread(opts$a)
map <- map %>% select(V4,hgnc_symbol,V9)
map <- map %>% rename(snp = V4,  gene_id = V9)

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

##################
# Manhattan plot #
##################

# positions
allpos <- fread(opts$m,
                header=F,
                sep="\t")
allpos <- allpos %>% filter(V1 %in% c(3,6,8))

man <- allpos %>%
  rename(chr=1, snp=2, morgans=3, bp=4) %>%
  arrange(chr, bp) %>%
  mutate(bp = cumsum(as.numeric(bp))) %>%
  merge(., ml, by="snp")

# set x-axis labels to the centre of each chromosome

axis.set <- man %>% 
  group_by(chr) %>% 
  summarise(center = (max(bp) + min(bp)) / 2, .groups="keep")  %>% rename(chr = chr)

axis.breaks <- man %>% 
  group_by(chr) %>% 
  summarise(min=min(bp), max=max(bp), .groups="keep") %>% rename(chr = chr)

# the number of chromosomes

nCHR <- nlevels(as.factor(man$chr))

# annotation file

man_ann <- loadings %>%
  melt(id="snp", variable.name="pc", value.name="value") %>%
  mutate(v2 = value^2) %>%
  group_by(pc) %>%
  top_n(v2, n=15) %>%
  ungroup() %>%
  merge(., map, by="snp") %>%
  dplyr::select(pc, snp, hgnc_symbol) %>%
  merge(., man, by=c("pc", "snp")) %>% mutate(id = paste(pc,hgnc_symbol)) %>%
  filter(!duplicated(id)) %>%
  as.data.frame() 


source("~/ddong//packages_load.R")
setwd("~/ddong/TEDDY/TEDDY/")
library(fgsea)
library(clusterProfiler)
#a <- fread("~/ddong/TEDDY/TEDDY/output_2024/GO_fugassem/my_output_20_loading_hgnc.txt",header = F,fill = T)
a <- fread("~/ddong/TEDDY_2024/GO/hg19_geneticPCsnp20_GO.txt",header = F,fill = T)
a <- melt(a,"V1")
gene_sets <- a %>% dplyr::select(V1,value) %>% dplyr::rename(Term = V1,Gene = value)
gene_sets <- gene_sets %>% filter(Gene != "")
grouped_list <- split(gene_sets$Gene,gene_sets$Term)

opts <- list()
opts$i <- "./genetics.bed/filter/teddy_out8.qc.excl_outliers.eigenvec.var"
opts$m <- "./genetics.bed/filter/teddy-out8_map.map"

loadings <- read.table(opts$i, header=T) %>%
  dplyr::select(2, 5:9) %>%
  dplyr::rename(snp=VAR)
loadings <- loadings %>% mutate(`PC3-PC2` = PC3-PC2,
                                `PC3-PC4` = PC3- PC4, 
                                `PC3-PC5` = PC3 - PC5)


map_select <- map %>% filter(snp %in% loadings$snp)
map_select <- map_select %>% select(snp,hgnc_symbol)
gene_sets_snp <- gene_sets %>% left_join(map_select,by = c("Gene" = "hgnc_symbol"))
gene_sets_snp <- gene_sets_snp %>% filter(!is.na(snp))
grouped_list_snp <- split(gene_sets_snp$snp,gene_sets_snp$Term)

load(file = "~/ddong/TEDDY_2024/output/snpfgsea_allmap.Rdata")
all_sig <- all_snpgsea %>% filter(padj< 0.1) 
man_ann <- man_ann %>% filter(chr %in% c(6,8) )  %>% filter(pc != "PC5")
man_ann <- man_ann %>% left_join(gene_sets_snp,by = c("hgnc_symbol" = "Gene",
                                                      "snp"))
man_ann <- man_ann %>% filter(!is.na(Term)) 
allsig2 <- all_sig %>% select(pathway,PC) %>% mutate(sig = T)
man_ann <- man_ann %>% left_join(allsig2,c("pc" = "PC","Term" = "pathway"))
man_ann <- man_ann %>% arrange(sig) %>% mutate(gene_pc = paste(pc,hgnc_symbol))
man_ann <- man_ann %>% filter(!duplicated(gene_pc))



man_ann_anti <- man_ann %>% filter(Term == "GO:0019730: [BP] antimicrobial humoral response")

anno <- fread("~/ddong/TEDDY_2024/anno.csv",header = F)
colnames(anno) <- colnames(man_ann)
man_ann <- man_ann 


pdf("~/ddong/TEDDY_project/figures/Fig3c",5,4.5)
ggplot(man, aes(x=bp, y=value)) +
  facet_grid(pc~., scales="fixed") +
  geom_point(aes(colour=as.factor(chr), fill=as.factor(chr)),size =0.7,) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  #scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred") +
  #scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred") +
  #scale_fill_manual(values = c("#184e77","#34a0a4","#d9ed92"))+
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
  labs(x="Chromosome", y="Loading", colour="Loading", fill="Loading") +
  geom_label_repel(data=anno2, aes(label=hgnc_symbol),max.overlaps = 20,size = 3)+
  scale_fill_brewer(palette = "Set3")

dev.off()
