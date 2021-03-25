library(tidyverse)
library(data.table)

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

df_sp <- mpa %>%
	select(-colnames(DF[,-2]))
	
vd <- vegan::vegdist(df_sp, method="bray")

perm <- permute::how(nperm = 999)
permute::setBlocks(perm) <- with(df_md, subject_id)

adonis_res <- data.frame(vegan::adonis2(vd ~ age + country + t1d + persist_conf_ab + ever_brstfed + abx_consumed + delivery + current_t1d_diag + current_persist_conf_ab + current_brst_fed + abx_this_month + PC1 + PC2 + PC3 + PC4 + PC5, 
		data=df_md, 
		permutations=perm,
		by="margin",
		parallel=1)) %>%
	tibble::rownames_to_column("Covariate")

write.table(adonis_res, "adonis_multivariate.tsv", sep="\t", quote=F, row.names=F)

#