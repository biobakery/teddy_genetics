#!/usr/bin/env Rscript

require(docopt)
'Usage:
   run_models_snps.R [-m <microbiome> -f <feature> -g <genetics> -d <metadata> -c <cpus> -o <out_dir>]

Options:
   -m microbiome data
   -f microbiome feature
   -g genetics data
   -d metadata
   -c number of cpus [default: 1]
   -o output directory

' -> doc

opts <- docopt(doc)

#

opts

library(tidyverse)
library(reshape2)
library(doParallel)

# microbiome data

if (opts$f == "Bray-Curtis") {
	
	microbiome <- data.table::fread(opts$m, sep="\t", header=T, check.names=F) %>%
		rename(sample_id=1) %>%
		column_to_rownames("sample_id") %>%
		t(.) %>%
		as.data.frame() %>%
		rownames_to_column("sample_id") %>%
		setNames(paste0("microbiome_", names(.))) %>%
		rename(sample_id=1)
		
	} else {
	
	microbiome <- data.table::fread(opts$m, sep="\t", header=T, check.names=F) %>%
		setNames(paste0("microbiome_", names(.))) %>%
		rename(sample_id=1)
	
	}

# genetics data

snps <- data.table::fread(opts$g, sep="\t", header=T, check.names=F) %>%
	setNames(paste0("genetics_", names(.))) %>%
	rename(Subject=1)

# metadata

metadata <- read.csv(opts$d, "\t", header=T) %>%
	rename(Subject=subject_id, Days=subject_age_days, Country=country)

# merge data and metadata

df <- merge(metadata, snps, by="Subject") %>%
	merge(., microbiome, by="sample_id")

#############################
# make a list of the models #
#############################

micro_list <- colnames(select(microbiome, contains("microbiome_"))) %>%
	gsub("microbiome_", "", .) %>%
	gsub("-", "_", .) %>%
	gsub("\\.", "_", .)

genet_list <- colnames(select(snps, contains("genetics_"))) %>%
	gsub("genetics_", "", .) %>%
	gsub("-", "_", .) %>%
	gsub("\\.", "_", .)

data <- df %>%
	setNames(gsub("microbiome_", "", names(.))) %>%
	setNames(gsub("genetics_", "", names(.))) %>%
	setNames(gsub("\\.|-", "_", names(.)))

formulae <- list()

for (i in micro_list) {
	for (j in genet_list) {
		formula <- data.frame(formula = paste0("lmer(", i, " ~ ", j, " * ", "Days + Country + ( 1 + Days | Subject ), data=data, na.action=na.omit)"))	
		formulae <- rbind(formulae, formula)
		}
	}

formulae <- levels(as.factor(formulae$formula))

#####################################
# loop through the list in parallel #
#####################################

registerDoParallel(cores=opts$c)

model_anova_total <- foreach(i=formulae, 
			      .combine=rbind, 
			      .packages=c("lmerTest", "tidyverse", "performance"),
			      .errorhandling = 'remove'
			      )	%dopar% {

	model <- eval(parse( text= i ))
	
	model_coef <- as.data.frame(coef(summary(model))[ , "Estimate"]) %>%
		rename(Coefficient=1) %>%
		rownames_to_column("Predictor") %>%
		filter(!grepl("Intercept|Country", Predictor))
	
	model_anova <- anova(model) %>%
		rownames_to_column("Predictor") %>%
		select(Predictor, `Pr(>F)`) %>%
		rename(P = `Pr(>F)`) %>%
		mutate(Model = i) %>%
		mutate(Coefficient = c(.subset2(model_coef, 2)[1], .subset2(model_coef, 2)[2], NA, .subset2(model_coef, 2)[3])) %>%
		as.data.frame()

	}

# output

dir.create(opts$o)

last_char <- str_sub(opts$o, -1, -1)

if (last_char == "/") {
	out_dir <- str_sub(opts$o, end=-2)
} else {
	out_dir <- opts$o
	}

part_a <- opts$m %>%
	gsub(".*\\/", "", .) %>%
	gsub("\\..*", "", .)

part_b <- opts$g %>%
	gsub(".*\\/", "", .) %>%
	gsub(".txt", ".tsv", .)

output <- paste0(out_dir, "/", part_a, ".", part_b, ".gz")

write.table(model_anova_total, gzfile(output), sep="\t", quote=F, row.names=F)

#
