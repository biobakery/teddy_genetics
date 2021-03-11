#!/usr/bin/env Rscript

require(docopt)
'Usage:
   filter_ab_prev_list.R [-i <microbiome> -m <metadata> -f <feature> -a <abundance> -p <prevalence> -o <out_dir>]

Options:
   -i input from MetaPhlAn2/HUMAnN2 (features as rows, samples as columns)
   -m metadata
   -a abundance [default: 1]
   -p prevalence [default: 0.25]
   -o output directory

' -> doc 

opts <- docopt(doc)

#

library(dplyr)
library(reshape2)

# microbiome data

MicroFeat <- read.csv(opts$i, sep="\t", header=T, check.names=F) %>%
	rename(MicroFeat=1)

#

MicroFeat_melt <- MicroFeat %>%
	melt(id="MicroFeat", variable.name="sample_id", value="value") %>%
	select(sample_id)

# metadata

metadata <- read.csv(opts$m, header=T, sep="\t") %>%
	select(subject_id, sample_id)

# count the number subjects

DF <- metadata %>%
	filter(sample_id %in% MicroFeat_melt$sample_id) %>%
	select(subject_id)

nsubs <- nlevels(as.factor(DF$subject_id))

# abundance + prevalence

ab <- as.numeric(opts$a)
prev <- as.numeric(opts$p)
N <- prev * nsubs

# prefix for output

prefix <- stringr::str_remove(opts$i, ".tsv") %>%
	gsub(".*\\/", "", .)

paste0(opts$o, "/",  prefix, ".major_list.tsv")

# list the features that meet the criteria

MicroFeat_list <- MicroFeat %>%
	filter(grepl("\\|s__", MicroFeat) & !grepl("\\|t__", MicroFeat)) %>%
	mutate(MicroFeat = gsub(".*\\|", "", MicroFeat)) %>%
	melt(id="MicroFeat", variable.name="sample_id", value="value") %>%
	merge(., read.csv(opts$m, header=T, sep="\t"), by="sample_id") %>%
	filter(value > ab) %>%
	select(MicroFeat, subject_id) %>%
	distinct() %>%
	group_by(MicroFeat) %>%
	tally() %>%
	ungroup() %>%
	filter(n > N) %>%
	select(MicroFeat) %>%
	as.data.frame()

write.table(MicroFeat_list, paste0(opts$o, "/", prefix, ".major_list.tsv"), quote=F, row.names=F)
	
#
