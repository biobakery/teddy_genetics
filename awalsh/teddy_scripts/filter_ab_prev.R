#!/usr/bin/env Rscript

require(docopt)
'Usage:
   filter_ab_prev.R [-i <microbiome> -f <feature> -a <abundance> -p <prevalence> -t <transformation>]

Options:
   -i input from MetaPhlAn2/HUMAnN2 (features as rows, samples as columns)
   -f type of feature (EC, Pfam, Species)
   -a abundance [default: 0]
   -p prevalence [default: 0.1]
   -t transformation (abs, log, log10, sqrt) [default: log]

' -> doc 

opts <- docopt(doc)

#

library(dplyr)
library(reshape2)

MicroFeat <- read.csv(opts$i, sep="\t", header=T, check.names=F) %>%
	rename(MicroFeat=1)

ab <- as.numeric(opts$a)
prev <- as.numeric(opts$p)
N <- prev * ncol(MicroFeat[-1])

prefix <- stringr::str_remove(opts$i, ".tsv")

# list the features that meet the criteria

if(grepl("Species", opts$f) == TRUE){

MicroFeat_list <- MicroFeat %>%
	filter(grepl("\\|s__", MicroFeat) & !grepl("\\|t__", MicroFeat)) %>%
	mutate(MicroFeat = gsub(".*\\|", "", MicroFeat)) %>%
	melt() %>%
	filter(value > ab) %>%
	group_by(MicroFeat) %>%
	tally() %>%
	ungroup() %>%
	filter(n > N) %>%
	select(MicroFeat) %>%
	as.data.frame()

write.table(MicroFeat_list, paste0(prefix, ".major_list.tsv"), quote=F, row.names=F)

} else {

MicroFeat_list <- MicroFeat %>%
	melt() %>%
	filter(value > ab) %>%
	group_by(MicroFeat) %>%
	tally() %>%
	ungroup() %>%
	filter(n > N) %>%
	filter(!grepl("UN", MicroFeat))  %>%
	select(MicroFeat) %>%
	as.data.frame()
}

# remove features that do not meet the criteria

if(grepl("Species", opts$f) == TRUE){

MicroFeat_df <- MicroFeat %>%
	filter(grepl("\\|s__", MicroFeat) & !grepl("\\|t__", MicroFeat)) %>%
	mutate(MicroFeat = gsub(".*\\|", "", MicroFeat)) %>%
	merge(., MicroFeat_list, by="MicroFeat") 

} else {
	
MicroFeat_df <- merge(MicroFeat_list, MicroFeat, by="MicroFeat") %>%
	mutate(MicroFeat = paste0(opts$f, ".", MicroFeat))

}

# calculate eps (the smallest non-zero value)

MicroFeat_melt <-  MicroFeat_df %>%
	melt() %>%
	dplyr::filter(value > 0)

eps <- min(MicroFeat_melt$value)

# transformation method

tran_meth <- opts$t

tran_func <- function(x) {
	sapply(x, function(x){
		eval(parse(text=paste0(tran_meth, "(", x, ")")))
})}

# transform the data

MicroTrans <- MicroFeat_df %>%
	melt() %>%
	mutate(value = value + eps) %>%
	mutate(value = tran_func(value)) %>%
	dcast(variable ~ MicroFeat, value.var="value") %>%
	rename(sample_id=1)

write.table(MicroTrans, paste0(prefix, ".major.tsv"), sep="\t", quote=F, row.names=F)
	
#