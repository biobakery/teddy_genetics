#!/usr/bin/env Rscript

require(docopt)

'Usage:
   samples_0-3_y.R [--metadata <metadata> --microbiome <microbiome>]

Options:
   --metadata sample metadata
   --microbiome microbiome data

' -> doc

opts <- docopt(doc)

###

library(tidyverse)

opts$microbiome <- ("~/Desktop/teddy2/data_derived/metaphlan2.tsv")
opts$metadata <- ("~/Desktop/teddy2/test/metadata.tsv")

# metadata

metadata_for_samples <- read.table(opts$metadata, header=T)

# samples < 3.0 y/o

samples <- metadata_for_samples %>%
	filter(subject_age_days <= 365 * 3) %>%
	select(sample_id)

# raw microbiome data

micro_raw <- read.table(opts$microbiome, header=T, row.names=1, check.names=F) %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("sample_id") %>%
	mutate(sample_id = gsub("_Abundance-RPKs", "", sample_id))

head(micro_raw)

# filtered microbiome data

micro_filtered <- merge(micro_raw, samples, by="sample_id") %>%
	column_to_rownames("sample_id") %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("MicroFeat")

micro_filtered[1:3, 1:3]

# filter by abundance + prevalence

ab <- as.numeric(opts$a)
prev <- as.numeric(opts$p)
N <- prev * ncol(MicroFeat[-1])

# list the features that meet the criteria

if(grepl("Species", opts$f) == TRUE){

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

} else {

MicroFeat_list <- MicroFeat %>%
	melt() %>%
	filter(value > ab) %>%
	group_by(MicroFeat) %>%
	tally() %>%
	ungroup() %>%
	filter(n > N) %>%
	filter(!grepl("UN", MicroFeat))  %>%
	select(MicroFeat) %>%
	as.data.frame()
}

prefix <- stringr::str_remove(opts$i, ".tsv")

write.table(MicroFeat_list, paste0(prefix, ".major_species_list.tsv"))

#