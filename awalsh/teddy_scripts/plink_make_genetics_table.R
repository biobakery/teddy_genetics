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

data <- read.csv(opts$i, header=F, sep="\t")

df <- data %>%
	select(2,3,5) %>%
	filter(!grepl("00", V5)) %>%
	mutate(V5 = gsub(" ", "", V5)) %>%
	rename(Subject=1, Variant=2, Genotype=3)

map <- read.csv(opts$r, header=F, sep="\t") %>%
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