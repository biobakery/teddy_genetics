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

# filtered microbiome data

micro_filtered <- merge(micro_raw, samples, by="sample_id") %>%
	column_to_rownames("sample_id") %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("MicroFeat")

# write table

ext <- gsub(".*/", "", opts$microbiome)

write.table(micro_filtered, paste0("days_0-1095.", ext), sep="\t", quote=F, row.names=F)

#