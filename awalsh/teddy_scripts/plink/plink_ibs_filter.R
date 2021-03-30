#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plink_ibs_filter.R [-i <input> -o <output>]

Options:
   -i input (.genome)
   -o output
   
' -> doc 

opts <- docopt(doc)

###

library(tidyverse)
library(cowplot)

genome <- data.table::fread(opts$i)

# calculate the median distance of each subject to every other subject
dst <- genome %>%
	group_by(IID1) %>%
	summarise(DST_median = median(DST)) %>%
	ungroup()

# copy the input to dst_in
dst_in <- dst %>%
	mutate(Stage = "Pre-QC")
	
# define outliers
outliers <- dst %>%
		filter(DST_median < quantile(DST_median, probs=c(.25, .75), na.rm = FALSE)[1] - 1.5 * IQR(DST_median))

# remove outliers iteratively
while (nrow(outliers) > 0) {

	dst <- dst %>%
		filter(DST_median > quantile(DST_median, probs=c(.25, .75), na.rm = FALSE)[1] - 1.5 * IQR(DST_median))

	outliers <- dst %>%
		filter(DST_median < quantile(DST_median, probs=c(.25, .75), na.rm = FALSE)[1] - 1.5 * IQR(DST_median))
	
	}

subjects <- dst %>%
	select(1) %>%
	mutate(fam = 0) %>%
	select(2, 1)

write.table(subjects, opts$o, col.names=F, row.names=F, quote=F)

###
