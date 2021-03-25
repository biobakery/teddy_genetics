#!/usr/bin/env Rscript

require(docopt)
'Usage:
   diagnostics_snps.R [-i <input> -p <p_value> -m <microbiome>  -c <cores> -f <feature> -g <genetics> -d <metadata> -a <age_group> -n <names> -o <out_dir>]

Options:
   -i model stats
   -p p-value [default: 5e-08]
   -m microbiome data
   -c cores for modeling [default: 1]
   -f microbiome feature (Bray-Curtis / Species / EC / Pfam)
   -g genetic data
   -d metadata
   -a age group (infants or toddlers)
   -n SNP annotations
   -o output directory [default: diagnostic_plots]
   
' -> doc

opts <- docopt(doc)

###

# setwd("~/Dropbox (Harvard University)/teddy_scripts/diagnostic_plots/test_data/")

# opts$i <- "species_pvalues_toddlers.tsv"
# opts$m <- "days_549-1095.metaphlan2.major.tsv"
# opts$f <- "Species"
# opts$d <- "metadata.tsv"
# opts$g <- "genetics.SNPs.tsv"
# opts$a <- "toddlers"
# opts$p <- "5e-7"
# opts$c <- 12
# opts$n <- "../../../tmp/genetics.maf-0.05.hwe-0.001.ibd-0.2.qc.rsid_to_gene.map.txt"

opts$c <- as.numeric(opts$c)
opts$p <- as.numeric(opts$p)

# output directories

dir.create(opts$o)
dir.create(paste0(opts$o, "/", opts$a))
dir.create(paste0(opts$o, "/", opts$a, "/fixed_effects"))
dir.create(paste0(opts$o, "/", opts$a, "/interaction_effects"))	

fe_dir <- paste0(opts$o, "/", opts$a, "/fixed_effects")

ie_dir <- paste0(opts$o, "/", opts$a, "/interaction_effects")

# colour palette

my_pal <- c("Homozygous major"="cornflowerblue", "Heterozygous"="#876192", "Homozygous minor"="firebrick")

# libraries

library(data.table)
library(tidyverse)
library(reshape2)
library(doParallel)
library(scales)
library(purrr)
library(broom)

# SNP annotations

snp_ann <- fread(opts$n) %>%
	rename(SNP = 1) %>%
	mutate(SNP = gsub("-", "_", SNP))

############
# p-values #
############

model_stats <- data.table::fread(opts$i, sep="\t", header=T)

if (opts$f == "EC" | opts$f == "Pfam") {
	
	path_map <- model_stats %>%
		mutate(Pathway = paste0(PathID, ": ", PathName)) %>%
		select(MicroFeat, Pathway) %>%
		distinct()
	
	}

# significant results

sig_micro <- model_stats %>%
	filter(P < opts$p) %>%
	select(MicroFeat, SNP, Type, P) %>%
	mutate(Model = paste0(MicroFeat, " ~ ", SNP))

sig_bugs <- c(sig_micro$MicroFeat) %>%
	as.data.frame() %>%
	distinct() %>%
	rename(MicroFeat = 1)

sig_snps <- c(sig_micro$SNP) %>%
	as.data.frame() %>%
	distinct() %>%
	rename(SNP = 1)

####################
# input for models #
####################

# genetics

genetics <- read.csv(opts$g, header=T, sep="\t", check.names=F) %>%
	rename(SNP=1) %>%
	mutate(SNP = gsub("-", "_", SNP)) %>%
	filter(SNP %in% levels(as.factor(sig_snps$SNP))) %>%
	column_to_rownames("SNP") %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("Subject")

# microbiome

micro <- read.csv(opts$m, header=T, sep="\t", row.names=1, check.names=F) %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("MicroFeat") %>%
	mutate(MicroFeat = gsub("\\.", "_", MicroFeat)) %>%
	filter(MicroFeat %in% levels(as.factor(sig_bugs$MicroFeat))) %>%
	column_to_rownames("MicroFeat") %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("sample_id")
	
# metadata

metadata <- read.csv(opts$d, header=T, sep="\t") %>%
	rename(Subject=subject_id, Days=subject_age_days, Country=country)

md_list <- c("Subject", "sample_id", "Country", "Days")

# merge

data <- merge(genetics, metadata, by="Subject") %>%
	merge(., micro, by="sample_id")
	
# list features

genet_list <- colnames(genetics[,-1])
micro_list <- colnames(micro[,-1])

####################
# rerun the models #
####################

# calculate the allele frequencies for each SNP

allele_map <- genetics %>%
	reshape2::melt(id="Subject", variable.name="SNP", value.name="Allele") %>%
	filter(Subject %in% data$Subject) %>%
	group_by(SNP, Allele) %>%
	tally() %>%
	ungroup() %>%
	drop_na() %>%
	mutate(A1 = str_sub(Allele, 1, 1)) %>%
	mutate(A2 = str_sub(Allele, 2, 2)) %>%
	mutate(Zygosity = case_when(A1 == A2 ~ "Homozygous", A1 != A2 ~ "Heterozygous")) %>%
	group_by(SNP) %>%
	mutate(temp = case_when(
		n == min(n) & Zygosity != "Heterozygous" ~ "minor",
		n != min(n) & Zygosity != "Heterozygous" ~ "major",
		Zygosity == "Heterozygous" ~ "")) %>%
	mutate(AlleleType = str_trim(paste0(Zygosity, " ", temp))) %>%
	mutate(AlleleType = gsub("Heterozygous .*", "Heterozygous", AlleleType)) %>%
	select(SNP, Allele, AlleleType, n)

# specify the minimum frequency of the minor homozygous allele

minor_hom <- allele_map %>%
	group_by(SNP) %>%
	filter(n == min(n)) %>%
	filter(n > 10) %>%
	.$SNP

# recode the alleles

allele_map_i <- allele_map %>%
	filter(AlleleType != "Homozygous minor") %>%
	mutate(value = Allele) %>%
	select(SNP, Allele, value)
	
allele_map_ii <- allele_map %>%
	filter(AlleleType == "Homozygous minor") %>%
	mutate(value = case_when(SNP %in% minor_hom ~ Allele)) %>%
	select(SNP, Allele, value)

allele_map_updated <- rbind(allele_map_i, allele_map_ii) %>%
	distinct()

# update the genetics table

genetics_updated <- genetics %>%
	filter(Subject %in% data$Subject) %>%
	reshape2::melt(id="Subject", variable.name="SNP", value.name="Allele") %>%
	merge(., allele_map_updated, by=c("SNP", "Allele")) %>%
	distinct() %>%
	reshape2::dcast(Subject ~ SNP, value.var="value")

# update the data

data <- merge(genetics_updated, metadata, by="Subject") %>%
	merge(., micro, by="sample_id")

# list the models

models <- sig_micro %>%
	mutate(Model = paste0("lmer(", MicroFeat, " ~ ", SNP, " * Days + Country + ( 1 + Days | Subject ), data=data, na.action=na.omit)")) %>%
	select(Model) %>%
	distinct() %>%
	.$Model

# run the models in parallel

registerDoParallel(cores=opts$c)

#models[1]

model_anova_total <- foreach(i=models, 
			      .combine=rbind, 
			      .packages=c("lmerTest", "tidyverse", "performance", "DHARMa"),
			      .errorhandling = 'remove'
			      )	%dopar% {

	model <- eval(parse( text= i ))

	model_coef <- as.data.frame(coef(summary(model))[ , "Estimate"]) %>%
		rename(Coefficient=1) %>%
		rownames_to_column("Predictor") %>%
		filter(!grepl("Intercept|Country", Predictor))
	
	R2_conditional <- (data.frame(r2(model)$R2_conditional))[1,1]
	
	simulationOutput <- simulateResiduals(fittedModel = model, refit = F)
	p_uni <- testUniformity(simulationOutput, plot=F)$p.value
	p_out <- testOutliers(simulationOutput, plot=F, type="binomial")$p.value
	p_dis <- testDispersion(simulationOutput, plot=F)$p.value
	
	model_anova <- anova(model) %>%
		rownames_to_column("Predictor") %>%
		select(Predictor, `Pr(>F)`) %>%
		rename(P = `Pr(>F)`) %>%
		mutate(Model = i) %>%
		mutate(Coefficient = c(.subset2(model_coef, 2)[1], .subset2(model_coef, 2)[2], NA, .subset2(model_coef, 2)[3])) %>%
		mutate(R2_conditional=R2_conditional, Uniformity=p_uni, Outliers=p_out, Dispersion=p_dis) %>%
		as.data.frame()

	}
	
# updated results

sig_micro <- model_anova_total %>%
	filter(Predictor != "Days" & Predictor != "Country") %>%
	filter(P < 5e-07) %>%
	mutate(MicroFeat = gsub(".*lmer\\(| .*", "", Model)) %>%
	mutate(SNP = gsub(".* ~ | .*", "", Model)) %>%
	mutate(Type = case_when(
		grepl(":", Predictor) ~ "Interaction", 
		!grepl(":", Model) ~ "Fixed")) %>%
	select(-Model, -Predictor) %>%
	mutate(Model = paste0(MicroFeat, " ~ ", SNP))

# add pathway names

if (opts$f == "EC" | opts$f == "Pfam") {
	
	sig_micro <- merge(sig_micro, path_map, by="MicroFeat")
	
	}
	
sm_updated <- merge(sig_micro, snp_ann, by="SNP")

write.table(sm_updated, paste0(str_to_lower(opts$f), "_p_values_", opts$a, ".QC.tsv"), sep="\t", row.names=F, quote=F)

# sig snps

genet_list_updated <- levels(as.factor(sig_micro$SNP))

allele_map <- genetics_updated %>%
	select(Subject, all_of(genet_list_updated)) %>%
	reshape2::melt(id="Subject", variable.name="SNP", value.name="Allele") %>%
	group_by(SNP, Allele) %>%
	tally() %>%
	ungroup() %>%
	drop_na() %>%
	mutate(A1 = str_sub(Allele, 1, 1)) %>%
	mutate(A2 = str_sub(Allele, 2, 2)) %>%
	mutate(Zygosity = case_when(A1 == A2 ~ "Homozygous", A1 != A2 ~ "Heterozygous")) %>%
	group_by(SNP) %>%
	mutate(temp = case_when(
		n == min(n) & Zygosity != "Heterozygous" ~ "minor",
		n != min(n) & Zygosity != "Heterozygous" ~ "major",
		Zygosity == "Heterozygous" ~ "")) %>%
	mutate(AlleleType = str_trim(paste0(Zygosity, " ", temp))) %>%
	mutate(AlleleType = gsub("Heterozygous .*", "Heterozygous", AlleleType)) %>%
	select(SNP, Allele, AlleleType, n) %>%
	mutate(Order = case_when(
		AlleleType == "Homozygous major" ~ 1,
		AlleleType == "Heterozygous" ~ 2,
		AlleleType == "Homozygous minor" ~ 3)) %>%
	ungroup() %>%
	arrange(SNP, Order) %>%
	mutate(X = paste0(SNP, "_", Order)) %>%
	mutate(Allele = as.character(unname(Allele))) %>%
	droplevels() %>%
	merge(., sig_micro, by="SNP") %>%
	mutate(n = paste0("n=", n))

if (opts$f == "EC" | opts$f == "Pfam") {
		
	allele_map <- allele_map %>%
		mutate(Model = paste0(gsub("\\:.*", "", Pathway), " ~ ", SNP))

	}

# updated results for fixed effects

sig_fe <- sig_micro %>%
	filter(Type == "Fixed")

# updated results for interaction effects
	
sig_int <- sig_micro %>%
	filter(Type == "Interaction")

#################
# fixed effects #
#################

species <- levels(as.factor(sig_fe$MicroFeat))

for (i in species){
	
	dir.create(paste0(fe_dir, "/", i))	
	
	}

allele_map_fe <- allele_map %>%
	filter(Type == "Fixed")

models_fe <- levels(as.factor(sig_fe$Model))

foreach(i=models_fe,  
	.packages=c("tidyverse"),
	.errorhandling = 'remove'
	)	%dopar% {
	
	microfeat <- gsub(" .*", "", i)
	
	snp <- gsub(".* ", "", i)
	
	# subset the data
	
	mdf <- data %>%
		select(all_of(md_list), all_of(microfeat), all_of(snp)) %>%
		rename(Abundance=5, Allele=6) %>%
		mutate(MicroFeat = microfeat) %>%
		mutate(SNP = snp) %>%
		group_by(MicroFeat, SNP, Allele, Subject) %>%
		summarise(Abundance = mean(Abundance)) %>%
		ungroup() %>%
		merge(., allele_map_fe, by=c("MicroFeat", "SNP", "Allele")) %>%
		mutate(AlleleType = factor(AlleleType, levels=c("Homozygous major", "Heterozygous", "Homozygous minor")))
	
	# calculate the maximum abundance of the feature
	
	max_ab <- mdf %>%
		mutate(Delta = max(Abundance) - min(Abundance)) %>%
		filter(Abundance == max(Abundance)) %>%
		mutate(Abundance = (0.1 * Delta) + Abundance) %>%
		select(Model, Abundance, Delta)
	
	# the p-value for the model
		
	pval_map <- sig_fe %>%
		filter(Model == i) %>%
		filter(P == min(P)) %>%
		select(Model, P)
	
	# annotation for the allele frequencies
	
	ann_n <- merge(allele_map_fe, max_ab, by="Model")
	
	# annotation for the p-values
		
	ann_p <- ann_n %>%
		mutate(P = paste0("p=", scientific(P))) %>%
		mutate(Abundance = (0.2 * Delta) + (Abundance)) %>%
		filter(Order == 2) %>%
		distinct() %>%
		filter(Type == "Fixed")
		
	if (nlevels(as.factor(mdf$X)) == 2) { nudge = -0.5 } 
	if (nlevels(as.factor(mdf$X)) == 3) { nudge = 0}
	
	am <- allele_map_fe %>% 
		filter(MicroFeat == microfeat & SNP == snp)
	
	custom_labels <- am$Allele
	names(custom_labels) <- am$X
	
	p <- ggplot(data=mdf, aes(x=X, y=Abundance)) +
		scale_x_discrete(labels = custom_labels) +
		facet_wrap(~Model, scales="free", ncol=1) +
		geom_jitter(aes(colour=AlleleType, fill=AlleleType), pch=21, alpha=0.5, width=0.25) +
		geom_boxplot(colour="black", fill=NA, outlier.shape=NA) +
		theme_bw() +
		theme(
			aspect.ratio = 1,
			plot.title = element_text(size=12.5, face="bold", hjust=0.5),
			axis.title = element_text(size=12.5, face="bold"),
			axis.text = element_text(size=12.5),
			legend.title = element_blank(),
			legend.text = element_text(size=12.5),
			strip.text = element_text(size=10, face="bold")) +
		scale_colour_manual(values = my_pal) +
		scale_fill_manual(values = my_pal) +
		geom_text(data=ann_n, aes(x=X, y=Abundance, label=n, colour=AlleleType), size=5, show.legend = FALSE) +
		geom_text(data=ann_p, aes(x=X, y=Abundance, label=P), size=5, nudge_x=nudge) +
		labs(x="Allele", y="log(abundance [%]) per subject")  +
			guides(color = guide_legend(override.aes=list(size = 3)))
	
	ggsave(plot=p, paste0(fe_dir, "/", microfeat, "/", microfeat, "_", snp, "_fe.png"), dpi=300, height=5, width=5/0.75)
		
	}

#####################
# top fixed effects #
#####################

if (opts$f == "EC" | opts$f == "Pfam") { 
	top_p <- 5e-08
} else {
	top_p <- opts$p
}

models_fe_top <- sig_fe %>%
	filter(P < top_p) %>%
	group_by(MicroFeat) %>%
	filter(P == min(P)) %>%
	ungroup() %>%
	.$Model

plot_list = list()

for (i in models_fe_top) {

	microfeat <- gsub(" .*", "", i)
	
	snp <- gsub(".* ", "", i)
	
	# subset the data
	
	mdf <- data %>%
		select(all_of(md_list), all_of(microfeat), all_of(snp)) %>%
		rename(Abundance=5, Allele=6) %>%
		mutate(MicroFeat = microfeat) %>%
		mutate(SNP = snp) %>%
		group_by(MicroFeat, SNP, Allele, Subject) %>%
		summarise(Abundance = mean(Abundance)) %>%
		ungroup() %>%
		merge(., allele_map_fe, by=c("MicroFeat", "SNP", "Allele")) %>%
		mutate(AlleleType = factor(AlleleType, levels=c("Homozygous major", "Heterozygous", "Homozygous minor")))
	
	# calculate the maximum abundance of the feature
	
	max_ab <- mdf %>%
		mutate(Delta = max(Abundance) - min(Abundance)) %>%
		filter(Abundance == max(Abundance)) %>%
		mutate(Abundance = (0.1 * Delta) + Abundance) %>%
		select(Model, Abundance, Delta)
	
	# the p-value for the model
		
	pval_map <- sig_fe %>%
		filter(Model == i) %>%
		filter(P == min(P)) %>%
		select(Model, P)
	
	# annotation for the allele frequencies
	
	ann_n <- merge(allele_map_fe, max_ab, by="Model")
	
	# annotation for the p-values
		
	ann_p <- ann_n %>%
		mutate(P = paste0("p=", scientific(P))) %>%
		mutate(Abundance = (0.2 * Delta) + (Abundance)) %>%
		filter(Order == 2) %>%
		distinct() %>%
		filter(Type == "Fixed")
		
	if (nlevels(as.factor(mdf$X)) == 2) { nudge = -0.5 } 
	if (nlevels(as.factor(mdf$X)) == 3) { nudge = 0 }
	
	am <- allele_map_fe %>% 
		filter(MicroFeat == microfeat & SNP == snp)
	
	custom_labels <- am$Allele
	names(custom_labels) <- am$X
	
	p <- ggplot(data=mdf, aes(x=X, y=Abundance)) +
		scale_x_discrete(labels = custom_labels) +
		facet_wrap(~Model, scales="free", ncol=1) +
		geom_jitter(aes(colour=AlleleType, fill=AlleleType), pch=21, alpha=0.5, width=0.25) +
		geom_boxplot(colour="black", fill=NA, outlier.shape=NA) +
		theme_bw() +
		theme(
			aspect.ratio = 1,
			plot.title = element_text(size=12.5, face="bold", hjust=0.5),
			axis.title = element_text(size=12.5, face="bold"),
			axis.text = element_text(size=12.5),
			legend.title = element_blank(),
			legend.text = element_text(size=12.5),
			strip.text = element_text(size=10, face="bold")) +
		scale_colour_manual(values = my_pal) +
		scale_fill_manual(values = my_pal) +
		geom_text(data=ann_n, aes(x=X, y=Abundance, label=n, colour=AlleleType), size=5, show.legend = FALSE) +
		geom_text(data=ann_p, aes(x=X, y=Abundance, label=P), size=5, nudge_x=nudge) +
		labs(x="Allele", y="log(abundance [%]) per subject")  +
		guides(color = guide_legend(override.aes=list(size = 3)))
	
	plot_list[[i]] = p
		
	}
	
pdf(paste0(opts$o, "/", opts$a, "/", str_to_lower(opts$f), "_diagnostic_plots_fixed_effects_", opts$a, ".pdf"), height=5, width=5/0.75)
plot_list
dev.off()

################
# interactions #
################

species <- levels(as.factor(sig_int$MicroFeat))

for (i in species){
	
	dir.create(paste0(ie_dir, "/", i))
	
	}

allele_map_int <- allele_map %>%
	filter(Type == "Interaction")

models_int <- levels(as.factor(sig_int$Model))

foreach(i=models_int,  
	.packages=c("tidyverse", "purrr", "broom", "scales"),
	.errorhandling = 'remove'
	)	%dopar% {

	microfeat <- gsub(" .*", "", i)
	
	snp <- gsub(".* ", "", i)

	mdi <- data %>%
		select(all_of(md_list), all_of(microfeat), all_of(snp)) %>%
		rename(Abundance=5, Allele=6) %>%
		mutate(MicroFeat = microfeat) %>%
		mutate(SNP = snp)

	# slope

	slope_df <- mdi %>%
		split(., list(.$Subject, .$MicroFeat), drop = TRUE) %>%
		map(~lm(Abundance ~ Days, data = .x)) %>% 
		map_df(tidy, .id = "Subject.MicroFeat") %>%
		filter(term == "Days") %>%
		mutate(Subject = gsub("\\..*", "", Subject.MicroFeat)) %>%
		mutate(MicroFeat = gsub(".*\\.", "", Subject.MicroFeat)) %>%
		select(-1, -2) %>%
		merge(., mdi, by=c("Subject", "MicroFeat")) %>%
		select(Subject, MicroFeat, SNP, Allele, estimate) %>%
		replace(is.na(.), 0) %>%
		distinct() %>%
		merge(., allele_map_int, by=c("MicroFeat", "SNP", "Allele")) %>%
		mutate(AlleleType = factor(AlleleType, levels=c("Homozygous major", "Heterozygous", "Homozygous minor")))
		
	# calculate the maximum slope for the feature

	max_slope <- slope_df %>%
		mutate(Delta = max(estimate) - min(estimate)) %>%
		filter(estimate == max(estimate)) %>%
		mutate(estimate = (0.1 * Delta) + estimate) %>%
		select(Model, estimate, Delta)

	# the p-value for the model
		
	pval_map <- sig_int %>%
		filter(Model == i) %>%
		filter(P == min(P)) %>%
		select(Model, P)
	
	# annotation for the allele frequencies
	
	ann_n <- merge(allele_map_int, max_slope, by="Model")
	
	# annotation for the p-values
		
	ann_p <- ann_n %>%
		mutate(P = paste0("p=", scientific(P))) %>%
		mutate(estimate = (0.2 * Delta) + (estimate)) %>%
		filter(Order == 2)
	
	# specify nudge_x for geom_text()
	
	if (nlevels(as.factor(slope_df$X)) == 2) { nudge = -0.5 } 
	if (nlevels(as.factor(slope_df$X)) == 3) { nudge = 0}
	
	# custom labels
	
	am <- allele_map_int %>% 
		filter(MicroFeat == microfeat & SNP == snp)
	
	custom_labels <- am$Allele
	names(custom_labels) <- am$X
	
	p <- ggplot(data=slope_df, aes(x=X, y=estimate)) +
		scale_x_discrete(labels = custom_labels) +
		facet_wrap(~Model, scales="free", ncol=1) +
		geom_jitter(aes(colour=AlleleType, fill=AlleleType), pch=21, alpha=0.5, width=0.25) +
		geom_boxplot(colour="black", fill=NA, outlier.shape=NA) +
		theme_bw() +
		theme(
			aspect.ratio = 1,
			plot.title = element_text(size=12.5, face="bold", hjust=0.5),
			axis.title = element_text(size=12.5, face="bold"),
			axis.text = element_text(size=12.5),
			legend.title = element_blank(),
			legend.text = element_text(size=12.5),
			strip.text = element_text(size=10, face="bold")) +
		scale_colour_manual(values = my_pal) +
		scale_fill_manual(values = my_pal) +
		geom_text(data=ann_n, aes(x=X, y=estimate, label=n, colour=AlleleType), size=5, show.legend = FALSE) +
		geom_text(data=ann_p, aes(x=X, y=estimate, label=P), size=5, nudge_x=nudge) +
		labs(x="Allele", y="Slope per subject")  +
			guides(color = guide_legend(override.aes=list(size = 3)))
	
	ggsave(plot=p, paste0(ie_dir, "/", microfeat, "/", microfeat, "_", snp, "_ie.png"), dpi=300, height=5, width=5/0.75)
	
	}

####################
# top interactions #
####################

models_int_top <- sig_int %>%
	filter(P < top_p) %>%
	group_by(MicroFeat) %>%
	filter(P == min(P)) %>%
	ungroup() %>%
	.$Model

plot_list = list()

for (i in models_int_top) {

	microfeat <- gsub(" .*", "", i)
	
	snp <- gsub(".* ", "", i)

	mdi <- data %>%
		select(all_of(md_list), all_of(microfeat), all_of(snp)) %>%
		rename(Abundance=5, Allele=6) %>%
		mutate(MicroFeat = microfeat) %>%
		mutate(SNP = snp)

	# slope

	slope_df <- mdi %>%
		split(., list(.$Subject, .$MicroFeat), drop = TRUE) %>%
		map(~lm(Abundance ~ Days, data = .x)) %>% 
		map_df(tidy, .id = "Subject.MicroFeat") %>%
		filter(term == "Days") %>%
		mutate(Subject = gsub("\\..*", "", Subject.MicroFeat)) %>%
		mutate(MicroFeat = gsub(".*\\.", "", Subject.MicroFeat)) %>%
		select(-1, -2) %>%
		merge(., mdi, by=c("Subject", "MicroFeat")) %>%
		select(Subject, MicroFeat, SNP, Allele, estimate) %>%
		replace(is.na(.), 0) %>%
		distinct() %>%
		merge(., allele_map_int, by=c("MicroFeat", "SNP", "Allele")) %>%
		mutate(AlleleType = factor(AlleleType, levels=c("Homozygous major", "Heterozygous", "Homozygous minor")))
		
	# calculate the maximum slope for the feature

	max_slope <- slope_df %>%
		mutate(Delta = max(estimate) - min(estimate)) %>%
		filter(estimate == max(estimate)) %>%
		mutate(estimate = (0.1 * Delta) + estimate) %>%
		select(Model, estimate, Delta)

	# the p-value for the model
		
	pval_map <- sig_int %>%
		filter(Model == i) %>%
		filter(P == min(P)) %>%
		select(Model, P)
	
	# annotation for the allele frequencies
	
	ann_n <- merge(allele_map_int, max_slope, by="Model")
	
	# annotation for the p-values
		
	ann_p <- ann_n %>%
		mutate(P = paste0("p=", scientific(P))) %>%
		mutate(estimate = (0.2 * Delta) + (estimate)) %>%
		filter(Order == 2)
	
	# specify nudge_x for geom_text()
	
	if (nlevels(as.factor(slope_df$X)) == 2) { nudge = -0.5 } 
	if (nlevels(as.factor(slope_df$X)) == 3) { nudge = 0}
	
	# custom labels
	
	am <- allele_map_int %>% 
		filter(MicroFeat == microfeat & SNP == snp)
	
	custom_labels <- am$Allele
	names(custom_labels) <- am$X
	
	p <- ggplot(data=slope_df, aes(x=X, y=estimate)) +
		scale_x_discrete(labels = custom_labels) +
		facet_wrap(~Model, scales="free", ncol=1) +
		geom_jitter(aes(colour=AlleleType, fill=AlleleType), pch=21, alpha=0.5, width=0.25) +
		geom_boxplot(colour="black", fill=NA, outlier.shape=NA) +
		theme_bw() +
		theme(
			aspect.ratio = 1,
			plot.title = element_text(size=12.5, face="bold", hjust=0.5),
			axis.title = element_text(size=12.5, face="bold"),
			axis.text = element_text(size=12.5),
			legend.title = element_blank(),
			legend.text = element_text(size=12.5),
			strip.text = element_text(size=10, face="bold")) +
		scale_colour_manual(values = my_pal) +
		scale_fill_manual(values = my_pal) +
		geom_text(data=ann_n, aes(x=X, y=estimate, label=n, colour=AlleleType), size=5, show.legend = FALSE) +
		geom_text(data=ann_p, aes(x=X, y=estimate, label=P), size=5, nudge_x=nudge) +
		labs(x="Allele", y="Slope per subject")  +
			guides(color = guide_legend(override.aes=list(size = 3)))
		
	plot_list[[i]] = p

	}

pdf(paste0(opts$o, "/", opts$a, "/", str_to_lower(opts$f), "_diagnostic_plots_interaction_", opts$a, ".pdf"), height=5, width=5/0.75)
plot_list
dev.off()

#
