#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plink_het.R [-i <input> -o <output>]

Options:
   -i input (.het)
   -o output
   
' -> doc 

opts <- docopt(doc)

###

library(tidyverse)

het <- read.table(opts$i, header=T) %>%
	mutate(HET_RATE = (N.NM. - O.HOM.)/N.NM.) %>%
	filter(
		HET_RATE < mean(HET_RATE) - 3 * sd(HET_RATE) | 
		HET_RATE > mean(HET_RATE) + 3 * sd(HET_RATE)
		) %>%
	select(FID, IID)

write.table(het,
	opts$o,
	col.names=FALSE,
	row.names=FALSE,
	quote=F)

###
