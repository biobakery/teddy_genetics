library(tidyverse)
library(viridis)
library(scales)
library(wesanderson)

setwd("~/Desktop/teddy2/test/")

# metadata

metadata <- data.table::fread("~/Desktop/teddy2/teddy_scripts/metadata.tsv")

########
# PCoA #
########

mpa <- read.csv("metaphlan2_major.tsv", sep="\t", header=T, check.names=F, row.names=1) %>% 
	t(.) %>%
	as.data.frame()

mpa.norm <- vegan::wisconsin(mpa)

mpa.norm.dist <- vegan::vegdist(mpa.norm, method="bray")

mpa.norm.pcoa <- ape::pcoa(mpa.norm.dist)

pcoa_vectors <- as.data.frame(mpa.norm.pcoa$vectors)[1:2] %>%
	rownames_to_column("sample_id")

write.table(pcoa_vectors, "metaphlan2_pcoa_vectors.tsv", sep="\t", quote=F, row.names=F)

pcoa_values <- data.frame(Relative_eig = mpa.norm.pcoa$values$Relative_eig) %>%
	mutate(Axis = paste0("Axis.", row_number())) %>%
	mutate(Relative_eig = 100 * Relative_eig) %>%
	select(Axis, Relative_eig)

write.table(pcoa_values, "metaphlan2_pcoa_values.tsv", sep="\t", quote=F, row.names=F)

#############
# plink pca #
#############

plink <- data.table::fread("~/Desktop/teddy2/test/genetics.PCA_PC1-PC20.tsv") %>%
	select(1:3) %>%
	rename(subject_id=IID) %>%
	merge(., metadata, by="subject_id") %>%
	select(-sample_id, -subject_age_days) %>%
	distinct() %>%
	mutate(country = recode(country, FIN="Finland", GER="Germany", SWE="Sweden", USA="USA"))
plink

p1 <- ggplot(plink, aes(x=X1, y=Y1)) +
	geom_point(colour="black", fill="firebrick", pch=21) +
	theme_bw() +
	labs(x="PC1", y="PC2")
p1

p1 <- ggplot(plink, aes(x=PC1, y=PC2)) +
	geom_point(aes(fill=country, colour=country), alpha=0.75, pch=21) +
	theme_bw() +
	theme(
	#	aspect.ratio=0.75,
		legend.position="bottom",
		legend.title = element_text(face="bold")
		) +
	scale_fill_viridis(discrete=T) +
	scale_colour_viridis(discrete=T)
p1

country_genetics <- ggExtra::ggMarginal(p1, groupColour = TRUE, groupFill = TRUE)

##################
# metaphlan2 pca #
##################

pcoa_plot_in <- read.csv("metaphlan2_pcoa_vectors.tsv", header=T, sep="\t") %>%
	merge(., metadata, by="sample_id") %>%
	filter(subject_age_days <= 1095) %>%
	mutate(country = recode(country, FIN="Finland", GER="Germany", SWE="Sweden", USA="USA"))

p_age <- ggplot(pcoa_plot_in, aes(x=Axis.1, y=Axis.2)) +
	geom_point(aes(colour=subject_age_days, fill=subject_age_days), alpha=0.75, pch=21) +
	scale_colour_distiller(palette="RdBu", labels=comma, breaks=c(365, 2*365)) +
	scale_fill_distiller(palette="RdBu", labels=comma, breaks=c(365, 2*365)) +
	theme_bw() +
	theme(
	#	aspect.ratio=0.75,
		legend.position="bottom",
		legend.title = element_text(face="bold")
		) +
	labs(x="PC1", y="PC2", colour="Age (days)", fill="Age (days)")

age_pcoa <- ggExtra::ggMarginal(p_age, colour=NA)

#

pcoa_country <- ggplot(pcoa_plot_in, aes(x=Axis.1, y=Axis.2)) +
	geom_point(aes(colour=country, fill=country), pch=21, alpha=0.75) +
	theme_bw() +
	theme(
	#	aspect.ratio=0.75,
		legend.position="bottom",
		legend.title = element_text(face="bold")
		) +
	labs(x="PC1", y="PC2", colour="Country", fill="Country") +
	scale_fill_viridis(discrete=T) +
	scale_colour_viridis(discrete=T)

country_bugs <- ggExtra::ggMarginal(pcoa_country, groupColour = TRUE, groupFill = TRUE)

pc_gen_bug <- cowplot::plot_grid(age_pcoa, country_bugs, country_genetics,
#	align="h", 
#	axis="tblr",
	nrow=2,
	scale=0.95,
	labels="AUTO")

ggsave(plot=pc_gen_bug, file="pca_gen_pcoa_bug_spectral.png", dpi=300, height=6.25, width=15)

###############
# Mantel test #
###############

###

library(tidyverse)
library(wesanderson)

############
# metadata #
############

md <- read.table("~/Desktop/teddy2/October_Models/metadata.tsv", header=T) %>%
	rename(Days=subject_age_days) %>%
	filter(Days <= 1095) %>%
	mutate(Bin = cut(Days, breaks = seq(min(Days), max(Days), by = (max(Days) - min(Days))/12 ))) %>%
	drop_na()

##############
# metaphlan2 #
##############

mpa <- read.csv("~/Desktop/teddy2/data_derived/metaphlan2_major.tsv", 
		sep="\t", header=T, check.names=F) %>% 
	rename(species=ID) %>%
	reshape2::melt(variable.name="sample_id", value.name="abundance") %>%
	merge(., md, by="sample_id")
	
mpa.mean <- mpa %>%
	group_by(subject_id, Bin, species) %>%
	summarise(abundance = mean(abundance))

mpa.median <- mpa %>%
	group_by(subject_id, Bin, species) %>%
	summarise(abundance = median(abundance))

levels(as.factor(mpa.mean$Bin)) == levels(as.factor(mpa.median$Bin)) 

#########
# PLINK #
#########

plink_dist_id <- read.table("~/Desktop/teddy2/test/plink.dist.id", header=F)

plink_dist <- read.table("~/Desktop/teddy2/test/plink.dist", header=F)
rownames(plink_dist) <- levels(as.factor(plink_dist_id$V2))
colnames(plink_dist) <- levels(as.factor(plink_dist_id$V2))

########
# loop #
########

mantel_res <- list()

for (i in levels(as.factor(mpa$Bin))) {
	
	DF.mean <- mpa.mean %>%
		filter(Bin == i) %>%
		reshape2::dcast(subject_id ~ species, value.var="abundance") %>%
		arrange(subject_id)

	DF.median <- mpa.median %>%
		filter(Bin == i) %>%
		reshape2::dcast(subject_id ~ species, value.var="abundance") %>%
		arrange(subject_id)
	
	pd.mean <- plink_dist %>%
		rownames_to_column("subject_id") %>%
		reshape2::melt() %>%
		filter( subject_id %in% DF.mean$subject_id ) %>%
		filter( variable %in% DF.mean$subject_id ) %>%
		reshape2::dcast(subject_id ~ variable, value.var="value") %>%
		arrange(subject_id) %>%
		column_to_rownames("subject_id")
	
	pd.median <- plink_dist %>%
		rownames_to_column("subject_id") %>%
		reshape2::melt() %>%
		filter( subject_id %in% DF.median$subject_id ) %>%
		filter( variable %in% DF.median$subject_id ) %>%
		reshape2::dcast(subject_id ~ variable, value.var="value") %>%
		arrange(subject_id) %>%
		column_to_rownames("subject_id")

	DF.mean <- DF.mean %>%
		filter( subject_id %in% rownames(pd.mean) ) %>%
		column_to_rownames("subject_id")

	DF.median <- DF.median %>%
		filter( subject_id %in% rownames(pd.median) ) %>%
		column_to_rownames("subject_id")
	
	PD.mean <- pd.mean %>% as.dist()
	PD.median <- pd.median %>% as.dist()
	
	mpa.mean.norm <- vegan::wisconsin(DF.mean)
	mpa.median.norm <- vegan::wisconsin(DF.median)

	mpa.mean.norm.dist <- vegan::vegdist(mpa.mean.norm, method="bray")
	mpa.median.norm.dist <- vegan::vegdist(mpa.median.norm, method="bray")
		
	MT.mean <- ade4::mantel.rtest(PD.mean, mpa.mean.norm.dist)
	MT.median <- ade4::mantel.rtest(PD.median, mpa.median.norm.dist)
	
	r.mean <- MT.mean$obs
	r.median <- MT.median$obs
	
	p.mean <- MT.mean$pvalue
	p.median <- MT.median$pvalue
	
	mt_df.mean <- data.frame(Bin=i, r=r.mean, p=p.mean, n=nrow(DF.mean), stat="mean / subject")
	mt_df.median <- data.frame(Bin=i, r=r.median, p=p.median, n=nrow(DF.median), stat="median / subject")
	
	mr <- rbind(mt_df.mean, mt_df.median)
	
	mantel_res <- rbind(mantel_res, mr)
	
	}

my_breaks = c(
	1/max(mantel_res$p),
	(1 /max(mantel_res$p) + 1/min(mantel_res$p))/2,
	1/min(mantel_res$p)
	)

my_labels = c(
	round(max(mantel_res$p), 2),
	round((max(mantel_res$p) + min(mantel_res$p))/2, 2),
	round(min(mantel_res$p), 2)
	)

head(mantel_res)

mantel_p <- ggplot(mantel_res, aes(x=Bin, y=100*r, colour=stat, fill=stat)) +
	geom_hline(yintercept=0, linetype="dashed", colour="grey") +
	geom_line(aes(group=stat)) +
	geom_point(aes(size=n), pch=21, fill="white") +
	geom_point(aes(size=n, alpha=1/p, fill=stat), pch=21, colour="black") +
	theme_bw() +
	theme(
		axis.text.x=element_text(angle=45, hjust=1)) +
	labs(x="Age (days)", y="r", alpha="p-value", colour="Statistic", fill="Statistic") +
	scale_alpha_continuous(breaks=my_breaks, labels=my_labels) +
	scale_colour_manual(values=c("steelblue", "goldenrod")) +
	scale_fill_manual(values=c("steelblue", "goldenrod"))
mantel_p

ggsave(plot=mantel_p, file="~/Desktop/teddy2/test/mantel_test_results.png", dpi=300, height=5, width=5/0.75)

pc_gene_bug <- cowplot::plot_grid(age_pcoa, country_bugs, country_genetics, mantel_p,
	align="hv", 
	axis="tblr",
	scale=0.95,
	labels="AUTO", nrow=2)

ggsave(plot=pc_gene_bug, file="~/Desktop/teddy2/test/fig_1_rough.png", dpi=300, height=10, width=12.5)

###