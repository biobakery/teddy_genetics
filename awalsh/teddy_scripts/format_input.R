#!/usr/bin/env Rscript

require(docopt)
'Usage:
   format_input.R [-m <microbiome> -g <genetics> -d <metadata> -o <output>]

Options:
   -m microbiome features
   -g genetic features
   -d metadata
   -o output

' -> doc 

opts <- docopt(doc)

#

library(tidyverse)
library(reshape2)

# microbiome data

microbiome <- read.csv(opts$m, sep="\t", header=T, check.names=F) %>%
	setNames(paste0("microbiome:", names(.))) %>%
	rename(sample_id=1)

# genetics data

snps <- read.csv(opts$g, sep="\t", header=T, check.names=F) %>%
	setNames(paste0("genetics:", names(.))) %>%
	rename(Subject=1)

# metadata

metadata <- read.csv(opts$d, "\t", header=T) %>%
	rename(Subject=subject_id, Days=subject_age_days, Country=country)

# merge data and metadata

data <- merge(metadata, snps, by="Subject") %>%
	merge(., microbiome, by="sample_id")

write.table(data, opts$o, sep="\t", quote=F, row.names=F)

###