#!/usr/bin/env Rscript

require(docopt)

'Usage:
   plink_make_genetics_table.R [-i <input_prefix> -o <output.tsv>]

Options:
   -i input prefix
   -o output (.tsv)

' -> doc

opts <- docopt(doc)

#####

library(tidyverse)
library(data.table)
library(reshape2)

# make a map to recode values

homozygous <- as.data.frame(expand.grid(c("A", "T", "C", "G"), c("A", "T", "C", "G"))) %>%
	filter(Var1 == Var2) %>%
	mutate(Combination = paste0(Var1, Var2)) %>%
	mutate(Genotype = Combination) %>%
	select(Genotype, Combination)

heterozygous <- as.data.frame(t(combn((unique(c(c('A', 'G', 'T', 'C'), c('A', 'G', 'T', 'C')))), 2))) %>%
	mutate(V3 = V2) %>%
	mutate(V4 = V1) %>%
	mutate(Forward = paste0(V1, V2)) %>%
	mutate(Reverse = paste0(V3, V4)) %>%
	mutate(Genotype = Forward) %>%
	select(Genotype, Forward, Reverse) %>%
	reshape2::melt(id="Genotype", value.name="Combination") %>%
	select(-variable)

indel_homozygous <- as.data.frame(expand.grid(c("I", "D"), c("I", "D"))) %>%
	filter(Var1 == Var2) %>%
	mutate(Combination = paste0(Var1, Var2)) %>%
	mutate(Genotype = Combination) %>%
	select(Genotype, Combination)
	
indel_heterozygous <- as.data.frame(t(combn((unique(c(c('I', 'D'), c('I', 'D')))), 2))) %>%
	mutate(V3 = V2) %>%
	mutate(V4 = V1) %>%
	mutate(Forward = paste0(V1, V2)) %>%
	mutate(Reverse = paste0(V3, V4)) %>%
	mutate(Genotype = Forward) %>%
	select(Genotype, Forward, Reverse) %>%
	reshape2::melt(id="Genotype", value.name="Combination") %>%
	select(-variable)

combinations <- rbind(homozygous, heterozygous, indel_homozygous, indel_heterozygous)

# define the inputs

in_map <- paste0(opts$i, ".map")

in_ped <- paste0(opts$i, ".ped")

# plink.map file

map <- fread(in_map)

# plink.ped file

ped <- fread(in_ped, header=F) %>%
	select(-c(1, 3, 4, 5, 6)) %>%
	column_to_rownames("V2") %>%
	setNames(map$V2) %>%
	rownames_to_column("subject_id") %>%
	reshape2::melt(id="subject_id", variable.name="Variant", value.name="Combination") %>%
	filter(!grepl("0", Combination)) %>%
	mutate(Combination = gsub(" ", "", Combination)) %>%
	merge(., combinations, by="Combination") %>%
	select(-Combination) %>%
	reshape2::dcast(Variant ~ subject_id, value.var="Genotype")

# output

write.table(ped, opts$o, sep="\t", row.names=F, quote=F)

#
