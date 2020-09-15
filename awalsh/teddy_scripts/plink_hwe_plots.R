#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plink_hwe_plot.R [-i <individuals> -t <SNPs>]

Options:
   -i input (.hwe)
   -t threshold for hwe [default: 0.001]

' -> doc 

opts <- docopt(doc)

#

library(tidyverse)
library(scales)

hardy <- read.table(opts$i, header=T, as.is=T) %>%
	mutate(P_log = -log10(P))

opts$t <- as.numeric(opts$t)

hwe_hist <- ggplot(hardy, aes(x=P_log)) +
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
	labs(title="The p-value of the exact test for HWE", 
		x=expression(bold(-log[10]~p-value)), 
		y="Number of SNPs") +
	scale_y_continuous(trans="sqrt", label=scientific) +
	geom_vline(xintercept=-log10(as.numeric(opts$t)), linetype="dashed", colour="firebrick")

#

hardy_summary <- hardy %>%
	mutate(Status = case_when(
		P < opts$t ~ "False", P >= opts$t ~ "True")) %>%
	group_by(Status) %>%
	tally() %>%
	ungroup() %>%
	as.data.frame()

hardy_bar <- ggplot(hardy_summary, aes(x=Status, y=n)) +
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
		trans="sqrt",
		labels=scientific) +
	ggtitle("Summary of HWE") +
	labs(x=paste0("p-value > ", opts$t), y="Number of SNPs")

prefix <- str_remove(opts$i, ".hwe")

ggsave(hwe_hist, file=paste0(prefix, "_plink_hwe_hist.png"), height=5, dpi=300)

ggsave(hardy_bar, file=paste0(prefix, "_plink_hwe_bar.png"), height=5, dpi=300)

###