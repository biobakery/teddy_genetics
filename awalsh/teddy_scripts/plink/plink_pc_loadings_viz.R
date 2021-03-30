#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plink_pc_loadings_viz.R [-i <input_loadings> -a <rsid_annotations> -m <rsid_map> -o <output_prefix>]

Options:
   -i loadings (e.g. plink.eigenvec.var)
   -a RsID annotations (.rsid_to_gene.map.txt)
   -m RsID map (e.g. plink.map)
   -o prefix for outputs
   
' -> doc

opts <- docopt(doc)

###

library(tidyverse)
library(reshape2)
library(ggrepel)
library(scales)

loadings <- read.table(opts$i, header=T) %>%
	select(2, 5:9) %>%
	rename(snp=VAR)

map <- read.csv(opts$a, sep="\t", header=T)

ml <- loadings %>%
	melt(id="snp", variable="pc", value.name="value")

ml_ann <- loadings %>%
	melt(id="snp", variable.name="pc", value.name="value") %>%
	mutate(v2 = value^2) %>%
	group_by(pc) %>%
	top_n(v2, n=10) %>%
	ungroup() %>%
	merge(., map, by="snp") %>%
	select(pc, snp, symbol, name) %>%
	merge(., ml, by=c("pc", "snp"))

##################
# Manhattan plot #
##################

# positions

man <- read.csv(opts$m,
		header=F,
		sep="\t") %>%
	rename(chr=1, snp=2, morgans=3, bp=4) %>%
	arrange(chr, bp) %>%
	mutate(bp = cumsum(as.numeric(bp))) %>%
	merge(., ml, by="snp")

# set x-axis labels to the centre of each chromosome

axis.set <- man %>% 
  group_by(chr) %>% 
  summarise(center = (max(bp) + min(bp)) / 2, .groups="keep") 

axis.breaks <- man %>% 
  group_by(chr) %>% 
  summarise(min=min(bp), max=max(bp), .groups="keep") 

# the number of chromosomes

nCHR <- nlevels(as.factor(man$chr))

# annotation file

man_ann <- loadings %>%
	melt(id="snp", variable.name="pc", value.name="value") %>%
	mutate(v2 = value^2) %>%
	group_by(pc) %>%
	top_n(v2, n=10) %>%
	ungroup() %>%
	merge(., map, by="snp") %>%
	select(pc, snp, symbol, name) %>%
	merge(., man, by=c("pc", "snp")) %>%
	as.data.frame()

# plot

manhplot <- ggplot(man, aes(x=bp, y=value)) +
	facet_wrap(~pc, scales="fixed", nrow=5) +
	geom_point(aes(colour=value, fill=value)) +
	geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
	scale_colour_gradient2(low = "steelblue", mid = "white", high = "darksalmon") +
	scale_fill_gradient2(low = "steelblue", mid = "white", high = "darksalmon") +
	theme_bw() +
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
		minor_breaks = axis.breaks$min) +
	labs(x="Chromosome", y="Loading", colour="Loading", fill="Loading") +
	geom_label_repel(data=man_ann, aes(label=symbol))

ggsave(plot=manhplot, paste0(opts$o, ".pc_loadings_manhattan.png"), dpi=300, height=2/3*15, width=2/3*20)

###
