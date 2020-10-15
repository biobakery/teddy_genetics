#!/usr/bin/env Rscript

require(docopt)

'Usage:
   plink_make_genetics_table.R [-i <input.lgen> -r <lgen-ref> -o <output.tsv>]

Options:
   -i input (.lgen)
   -r lgen-ref (.ref)
   -o output (.tsv)

' -> doc

opts <- docopt(doc)

#####

library(tidyverse)
library(reshape2)

data <- data.table::fread(opts$i, sep=" ", header=F) %>%
	mutate(V1 = gsub("0\t", "", V1))

df <- data %>%
	mutate(V5 = paste0(V3, V4)) %>%
	filter(!grepl("00", V5)) %>%
	select(V1, V2, V5) %>%
	rename(Subject=V1, Variant=V2, Genotype=V5)

map <- data.table::fread(opts$r, header=F) %>%
	rename(Variant=1, Major=2, Minor=3)

gen <- df %>%
	select(Variant, Genotype) %>%
	distinct() %>%
	inner_join(., map, by="Variant") %>%
	group_by(Variant, Genotype) %>%
	mutate(Allele =	case_when(grepl(Major, Genotype) ~ "Major", !grepl(Major, Genotype) ~ "Minor")) %>%
	ungroup() %>%
	select(Variant, Genotype, Allele)

DF <- inner_join(gen, df, by=c("Variant", "Genotype")) %>%
	select(Subject, Variant, Allele) %>%
	dcast(Variant ~ Subject, value.var="Allele")

write.table(DF, opts$o, sep="\t", row.names=F, quote=F)

#####