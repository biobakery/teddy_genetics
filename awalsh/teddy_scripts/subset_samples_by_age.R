#!/usr/bin/env Rscript

require(docopt)

'Usage:
   subset_samples_by_age.R [--metadata <metadata> --min-a <min_1> --max-a <max_1> --min-b <min_2> --max-b <max_2> --microbiome <microbiome>]

Options:
   --metadata sample metadata
   --min-a minimum age in days for infant subset
   --max-a maximum age in days for infant subset
   --min-b minimum age in days for toddler subset
   --max-b maximum age in days for toddler subset
   --microbiome microbiome data

' -> doc

opts <- docopt(doc)

###

library(tidyverse)

opts$min_a <- as.numeric(opts$min_a)
opts$max_a <- as.numeric(opts$max_a)
opts$min_b <- as.numeric(opts$min_b)
opts$max_b <- as.numeric(opts$max_b)

#########
# samples

metadata_for_samples <- read.table(opts$metadata, header=T)

# infant samples

samples_1 <- metadata_for_samples %>%
	filter(subject_age_days >= opts$min_a & subject_age_days <= opts$max_a) %>%
	select(sample_id)

# toddler samples

samples_2 <- metadata_for_samples %>%
	filter(subject_age_days >= opts$min_b & subject_age_days <= opts$max_b) %>%
	select(sample_id)

############
# microbiome

suffix <- stringr::str_remove(opts$microbiome, ".tsv") %>%
	gsub(".*/", "", .)

microbiome <- read.table(opts$microbiome, header=T, row.names=1, check.names=F) %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("sample_id") %>%
	mutate(sample_id = gsub("_Abundance-RPKs", "", sample_id))

micro_subset_1 <- merge(microbiome, samples_1, by="sample_id") %>%
	column_to_rownames("sample_id") %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("MicoFeat")

write.table(micro_subset_1, 
	paste0("days_", opts$min_a, "-", opts$max_a, ".", suffix, ".tsv"),
	sep="\t", quote=F, row.names=F)

micro_subset_2 <- merge(microbiome, samples_2, by="sample_id") %>%
	column_to_rownames("sample_id") %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("MicoFeat")

write.table(micro_subset_2, 
	paste0("days_", opts$min_b, "-", opts$max_b, ".", suffix, ".tsv"),
	sep="\t", quote=F, row.names=F)

###