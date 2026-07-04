#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plink_king_tidy.R [-p <prefix_output>]

Options:
   -p prefix for outputs
   
' -> doc 

opts <- docopt(doc)

###

library(tidyverse)
library(reshape2)
library(scales)

prefix <- opts$p

king_raw <- read.csv("tmp/temp_autoQC_smry.txt", header=F, sep="\t") %>%
	mutate(V4 = case_when(
		grepl("SNP",V1) ~ "SNPs",
		!grepl("SNP",V1) & !grepl("data", V1) ~ "Subjects",
		grepl("data", V1) ~ "Summary")) %>%
	mutate(V1 = gsub(" \\(removed\\)", "", V1)) %>%
	mutate(V2 = gsub("\\(|\\)", "", V2))

###########
# summary #
###########

king_summary <- king_raw %>%
	filter(grepl("Raw data counts|Final QC'ed data", V1)) %>%
	mutate(V1 = gsub("Raw data counts", "Raw", V1)) %>%
	mutate(V1 = gsub("Final QC'ed data", "QC'ed", V1)) %>%
	mutate(V1 = factor(V1, levels=c("Raw", "QC'ed"))) %>%
	mutate(V2 = as.numeric(V2)) %>%
	rename(Data=V1, Subjects=V2, SNPs=V3) %>%
	select(1, 2, 3) %>%
	melt(id="Data", variable.name="Type", value.name="Number")

write.table(king_summary, file=paste0(prefix, "_king_autoQC_Summary.txt"),
	sep="\t",
	quote=F,
	row.names=F)

#############
# breakdown #
#############

king_snp <- king_raw %>%
	filter(grepl("SNPs", V4) & V2 > 0) %>%
	select(-V3) %>%
	rename(Step=V1, Number=V2)

king_ind <- king_raw %>%
	filter(grepl("Subjects", V4) & V2 > 0) %>%
	select(-V3) %>%
	rename(Step=V1, Number=V2)
	
king_snp_ind <- rbind(king_snp, king_ind) %>%
	select(3, 1, 2) %>%
	rename(Type=V4)
	
write.table(king_snp_ind, file=paste0(prefix, "_king_autoQC_Breakdown.txt"),
	sep="\t",
	quote=F,
	row.names=F)

###