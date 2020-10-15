#!/usr/bin/env Rscript

require(docopt)
'Usage:
   filter_ab_prev_list.R [-i <microbiome> -a <abundance> -p <prevalence> -t <transformation>]

Options:
   -i input from MetaPhlAn2 (features as rows, samples as columns)
   -a abundance [default: 0]
   -p prevalence [default: 0.1]
   -t transformation (abs, log, log10, sqrt) [default: log]

' -> doc 

opts <- docopt(doc)

#

library(dplyr)
library(reshape2)

MicroFeat <- read.csv(opts$i, sep="\t", header=T, check.names=F) %>%
	rename(MicroFeat=1)

ab <- as.numeric(opts$a)
prev <- as.numeric(opts$p)
N <- prev * ncol(MicroFeat[-1])

prefix <- gsub(".*/", "", stringr::str_remove(opts$i, ".metaphlan2.tsv"))

MicroFeat_list <- MicroFeat %>%
	filter(grepl("\\|s__", MicroFeat) & !grepl("\\|t__", MicroFeat)) %>%
	mutate(MicroFeat = gsub(".*\\|", "", MicroFeat)) %>%
	melt() %>%
	filter(value > ab) %>%
	group_by(MicroFeat) %>%
	tally() %>%
	ungroup() %>%
	filter(n > N) %>%
	select(MicroFeat) %>%
	as.data.frame()

write.table(MicroFeat_list, paste0(prefix, ".major_bugs_list.tsv"), quote=F, row.names=F)
	
#