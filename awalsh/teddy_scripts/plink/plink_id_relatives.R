#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plink_id_relatives.R [-g <genome> -t <pihat> -m <missingness> -o <output>]

Options:
   -g genome (.genome)
   -t threshold (pihat)
   -m missingness (.imiss)
   -o output

' -> doc 

opts <- docopt(doc)

#

library(tidyverse)
library(reshape2)

threshold <- as.numeric(opts$t)

subjects_to_remove <- read.table(opts$g, header=T, as.is=T) %>%
	filter(PI_HAT > threshold) %>%
	mutate(Pair = paste0("Pair", row_number())) %>%
	select(Pair, IID1, IID2) %>%
	melt(id="Pair", variable.name="variable", value.name="IID") %>%
	select(Pair, IID) %>%
	merge(., read.table(opts$m, header=T, as.is=T), by="IID") %>%
	group_by(Pair) %>%
	slice(which.max(F_MISS)) %>%
	ungroup() %>%
	mutate(FID=c(0)) %>%
	select(FID, IID) %>%
	as.data.frame()

write.table(subjects_to_remove, 
	opts$o,
	col.names=F,
	row.names=F,
	quote=F)

#
