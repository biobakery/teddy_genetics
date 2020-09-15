#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plink_map.R [-i <input_map> -o <output_bed>]

Options:
   -i input (.map)
   -o output (.bed)

' -> doc

opts <- docopt(doc)

#

library(tidyverse)

map <- read.csv("GCF_000001405.39_GRCh38.p13_assembly_report.txt", sep="\t", header=T) %>%
	dplyr::select(7, 10) %>%
	dplyr::rename(RefSeq=1, UCSC=2) %>%
	distinct()
head(map)

genetics <- read.csv(opts$i, sep="\t", header=F) %>%
	rename(UCSC=V1, start=V4, rs=2) %>%
	select(UCSC, start, rs) %>%
	mutate(end = start+1) %>%
	select(UCSC, start, end, rs) %>%
	mutate(UCSC = paste0("chr", UCSC)) %>%
	mutate(UCSC = gsub("25", "X", UCSC)) %>%
	merge(., map, by="UCSC") %>%
	select(RefSeq, start, end, rs)
head(genetics)

paste(opts$o)

write.table(genetics, file=paste0(opts$o), col.names=F, row.names=F, sep="\t", quote=F)

###