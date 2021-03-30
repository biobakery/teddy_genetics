#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plink_rsid_to_gene.R [-i <input> -o <output>]

Options:
   -i input
   -o output

' -> doc

opts <- docopt(doc)

#

library(tidyverse)

map <- read.csv("tmp/hgnc_complete_set.txt", header=T, sep="\t")

df <- read.csv(opts$i, header=F, sep="\t") %>%
	filter(grepl("HGNC", V13)) %>%
	mutate(V13 = gsub(".*,HGNC\\:", "", V13)) %>%
	mutate(V13 = gsub("\\;.*", "", V13)) %>%
	mutate(V13 = gsub("\\,.*", "", V13)) %>%
	rename(hgnc_id=V13) %>%
	merge(map, ., by="hgnc_id") %>%
	select(7, 1:3, 16) %>%
	rename(snp=1, distance=5) %>%
	distinct()

write.table(df, file=paste0(opts$o), sep="\t", row.names=F, quote=F)

#
