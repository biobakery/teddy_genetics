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
	rename(snp=VAR) %>%
	mutate(index=row_number())

map <- read.csv(opts$a, sep="\t", header=T)

################
# scatter plot #
################

ml <- loadings %>%
	melt(id=c("snp", "index"), variable="pc", value.name="value")

ml_ann <- loadings %>%
	melt(id="snp", variable.name="pc", value.name="value") %>%
	mutate(v2 = value^2) %>%
	group_by(pc) %>%
	top_n(v2, n=10) %>%
	ungroup() %>%
	merge(., map, by="snp") %>%
	select(pc, snp, symbol, name) %>%
	merge(., ml, by=c("pc", "snp"))

p <- ggplot(ml, aes(x=index, y=value)) +
	facet_wrap(~pc, nrow=nlevels(as.factor(ml$pc)), scales="free_y") +
	geom_point(aes(colour=value, fill=value)) +
	scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred") +
	scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred") +
	theme_bw() +
	theme(
		axis.title = element_text(size=12.5, face="bold"),
		axis.text = element_text(size=12.5),
		legend.title = element_text(size=12.5, face="bold", vjust=0.75),
		legend.text = element_text(size=8.75),
		legend.position = "top",
		strip.text = element_text(size=12.5, face="bold", angle=0),
		panel.background=element_rect(fill="white"),
		panel.grid=element_blank()) +
	geom_hline(yintercept=0, linetype="dashed", colour="grey") +
	geom_label_repel(data=ml_ann, aes(x=index, y=value, label=symbol), segment.alpha=0.5) +
	scale_x_continuous(labels=scientific) +
	scale_y_continuous(expand=c(0.1, 0.1)) +
	labs(x="Index", y="Loading", colour="Loading", fill="Loading")

ggsave(plot=p, paste0(opts$o, "_pc_loadings_scatter.png"), dpi=300, height=20, width=15)

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
	merge(., man, by=c("pc", "snp"))

# plot

manhplot <- ggplot(man, aes(x=bp, y=value)) +
	facet_wrap(~pc, scales="fixed", nrow=5) +
	geom_point(aes(colour=value, fill=value)) +
	geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
	scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred") +
	scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred") +
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

ggsave(plot=manhplot, paste0(opts$o, "_pc_loadings_man_pruned.png"), dpi=300, height=15, width=20)

###