#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plink_remove_outliers.R [-i <input> -s <subjects> -o <output>]

Options:
   -i input
   -o output

' -> doc 

opts <- docopt(doc)

#

library(tidyverse)
library(reshape2)

# PCA

plink <- read.csv(opts$i, sep="\t", header=T)

genetics_pca <- plink[,2:7] %>%
	rename(subject_id=1)

# remove PCA outliers

df <- genetics_pca %>%
	mutate(subject_id = factor(subject_id)) %>%
	melt() %>%
	select(subject_id, variable, value) %>%
	group_by(variable) %>%
	filter(value > quantile(value, probs=c(.25, .75), na.rm = FALSE)[1] - 1.5*IQR(value) &
	       value < quantile(value, probs=c(.25, .75), na.rm = FALSE)[2] + 1.5*IQR(value)) %>%
	ungroup() %>%
	dcast(subject_id ~ variable, value.var="value") %>%
	drop_na() %>%
	select(subject_id) %>%
	mutate(col1 = c(0)) %>%
	select(col1, subject_id)

write.table(df, opts$o, col.names=F, row.names=F, quote=F, sep="\t")

#
