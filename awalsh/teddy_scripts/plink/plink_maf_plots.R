#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plink_maf_plots.R [-i <input> -t <threshold>]

Options:
   -i input (.frq)
   -t threshold for maf

' -> doc 

opts <- docopt(doc)

#

library(tidyverse)
library(reshape2)
library(scales)



maf <- read.table(opts$i, header=T, as.is=T)

opts$t <- as.numeric(opts$t)

# hist

maf_hist <- ggplot(maf, aes(x=MAF)) +
	geom_histogram(fill="goldenrod", colour="black") +
	theme_bw() +
	theme(panel.grid=element_blank(),
		plot.title=element_text(size=15, face="bold", hjust=0.5),
		axis.title=element_text(size=12.5, face="bold"),
		axis.text=element_text(size=12.5),
		strip.text=element_text(size=12.5, face="bold"),
		legend.title=element_text(size=12.5, face="bold"),
		legend.text=element_text(size=12.5),
		aspect.ratio=0.5) +
	labs(title="Minor allele frequency distribution", 
		x="Minor allele frequency", 
		y="Number of SNPs") +
	scale_y_continuous(labels=scientific) +
	geom_vline(xintercept=opts$t, linetype="dashed", colour="firebrick")

# bar

maf_summary <- maf %>%
	mutate(Status = case_when(
		MAF < opts$t ~ "False", MAF >= opts$t ~ "True")) %>%
	group_by(Status) %>%
	tally() %>%
	ungroup() %>%
	as.data.frame()

maf_bar <- ggplot(maf_summary, aes(x=Status, y=n)) +
	geom_bar(stat="identity", width=0.75, color="black", fill="goldenrod") +
	geom_text(aes(label=scientific(n)), size=5, vjust=-0.25) +
	theme_bw() +
	theme(panel.grid=element_blank(),
		plot.title=element_text(size=15, face="bold", hjust=0.5),
		axis.title=element_text(size=12.5, face="bold"),
		axis.text=element_text(size=12.5),
		strip.text=element_text(size=12.5, face="bold"),
		legend.title=element_text(size=12.5, face="bold"),
		legend.text=element_text(size=12.5),
		aspect.ratio=1) +
	scale_y_continuous(expand = expansion(mult = c(0, .1)),
		labels=scientific) +
	ggtitle("Summary of MAF") +
	labs(x=paste0("MAF > ", opts$t), y="Number of pairs")
maf_bar

# save

prefix <- str_remove(opts$i, ".frq")

ggsave(maf_hist, file=paste0(prefix, "_plink_maf_hist.png"), height=5, dpi=300)

ggsave(maf_bar, file=paste0(prefix, "_plink_maf_bar.png"), height=5, dpi=300)
