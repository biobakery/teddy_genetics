#!/usr/bin/env Rscript

require(docopt)

'Usage:
   subset_samples_by_age.R [--sam <subject_metadata> --minI <min_age_1> --maxI <max_age_1> --minT <min_age_2> --maxT <max_age_2> --micro <microbiome>]

Options:
   --sam sample metadata
   --minI minimum age in days for infant subset
   --maxI maximum age in days for infant subset
   --minT minimum age in days for toddler subset
   --maxT maximum age in days for toddler subset
   --micro microbiome data

' -> doc

opts <- docopt(doc)

###

library(tidyverse)

opts$minI <- as.numeric(opts$minI)
opts$maxI <- as.numeric(opts$maxI)
opts$minT <- as.numeric(opts$minT)
opts$maxT <- as.numeric(opts$maxT)

#########
# samples

metadata_for_samples <- read.table(opts$sam, header=T)
head(metadata_for_samples)

# infant samples

samples_1 <- metadata_for_samples %>%
	filter(subject_age_days >= opts$minI & subject_age_days <= opts$maxI) %>%
	select(sample_id)

# toddler samples

samples_2 <- metadata_for_samples %>%
	filter(subject_age_days >= opts$minT & subject_age_days <= opts$maxT) %>%
	select(sample_id)

############
# microbiome

suffix <- str_remove(opts$micro, ".tsv")

microbiome <- read.table(opts$micro, header=T, row.names=1, check.names=F) %>%
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
	paste0("days_", opts$minI, "-", opts$maxI, "_", suffix, ".tsv"),
	sep="\t", quote=F, row.names=F)

micro_subset_2 <- merge(microbiome, samples_2, by="sample_id") %>%
	column_to_rownames("sample_id") %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("MicoFeat")

write.table(micro_subset_2, 
	paste0("days_", opts$minT, "-", opts$maxT, "_", suffix, ".tsv"),
	sep="\t", quote=F, row.names=F)

###