library(tidyverse)
library(data.table)
library(doParallel)

#setwd("~/Desktop/teddy2/test/")

metadata_updated <- fread("permanova_metadata.tsv")

plink <- fread("genetics.PCA_PC1-PC20.tsv") %>%
	rename(subject_id=1) %>%
	select(1:6)

DF <- merge(metadata_updated, plink, by="subject_id") %>%
	arrange(subject_id, sample_id) %>%
	distinct()

mpa <- read.csv("metaphlan2_major.tsv", sep="\t", header=T, check.names=F, row.names=1) %>% 
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("sample_id") %>%
	merge(., DF, by="sample_id") %>%
	column_to_rownames("sample_id") %>%
	drop_na()

df_md <- mpa %>%
	select(colnames(DF[,-2])) %>%
	mutate_if(is.character, as.factor) %>%
	rename(age=mgx_age)

md <- colnames(df_md)

df_sp <- mpa %>%
	select(-colnames(DF[,-2]))

# generate distance matrix

vd <- vegan::vegdist(df_sp, method="bray")

# run loop

ncores <- detectCores()

if (ncores < length(md)) { n = ncores / 2 }

if (ncores >= length(md)) { n = length(md)}

registerDoParallel(cores=n)

adonis_combined <- foreach(i=md, 
			      .combine=rbind, 
			      .packages=c("tidyverse", "vegan"),
			      .errorhandling="remove"
			      )	%dopar% {

formula <- paste0("vegan::adonis2(vd ~ age + ", i, ", data=df_md, permutations=perm, method='bray')")

perm <- permute::how(nperm = 999)
permute::setBlocks(perm) <- with(df_md, subject_id)

adonis_res <- data.frame(eval(parse( text = formula ))) %>%
	tibble::rownames_to_column("Covariate") %>%
	mutate(Formula = paste0("age ~ ", i))

}

write.table(adonis_combined, "adonis_univariate.tsv", sep="\t", quote=F, row.names=F)

#
