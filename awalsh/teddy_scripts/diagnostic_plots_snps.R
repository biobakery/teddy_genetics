#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plot_results_snps.R [-i <input> -p <p_value> -m <microbiome> -g <genetics> -d <metadata> -a <age_group> -o <out_dir>]

Options:
   -i model stats
   -p p-value [default: 5e-08]
   -m microbiome data
   -g genetic data
   -d metadata
   -a age group (infants or toddlers)
   -o output directory [default: diagnostic_plots]
   
' -> doc

opts <- docopt(doc)

###

opts$i <- "species_pvalues_infants.tsv"
opts$m <- "days_0-548.metaphlan2.major.tsv"
opts$d <- "metadata.tsv"
opts$g <- "../../../plink_output/genetics.SNPs.tsv"
opts$a <- "toddler"
opts$p <- 5e-07

opts$p <- as.numeric(opts$p)

dir.create(opts$o)

library(tidyverse)
library(reshape2)
library(ggrepel)
library(lmerTest)
library(scales)
library(cowplot)

############
# p-values #
############

model_stats <- data.table::fread(opts$i, sep="\t", header=T)

# significant results

sig_micro <- model_stats %>%
	filter(P < opts$p) %>%
	select(MicroFeat, SNP, Type, P)

my_pal <- c("Homozygous major"="cornflowerblue", "Heterozygous"="#876192", "Homozygous minor"="firebrick")

if (nrow(sig_micro) > 0) {
	
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

	# metaphlan2

	mphlan <- read.csv(opts$m, header=T, sep="\t", row.names=1) %>%
		t(.) %>%
		as.data.frame() %>%
		rownames_to_column("Species") %>%
		mutate(Species = gsub("microbiome\\.", "", Species)) %>%
		filter(Species %in% levels(as.factor(sig_bugs$MicroFeat))) %>%
		column_to_rownames("Species") %>%
		t(.) %>%
		as.data.frame() %>%
		rownames_to_column("sample_id")

	# metadata

	metadata <- read.csv(opts$d, header=T, sep="\t") %>%
		rename(Subject=subject_id, Days=subject_age_days, Country=country)

	md_list <- c("Subject", "sample_id", "Country", "Days")

	# merge

	data <- merge(genetics, metadata, by="Subject") %>%
		merge(., mphlan, by="sample_id")

	# list features

	genet_list <- colnames(genetics[,-1])

	micro_list <- colnames(mphlan[,-1])

	#################
	# fixed effects #
	#################
	
	# fixed effect

	sig_fe <- sig_micro %>%
		filter(Type == "Fixed") %>%
		mutate(Model = paste0("lmer(", MicroFeat, " ~ ", SNP, " * ", "Day + Country + ( 1 + Day | Subject ), data=data, na.action=na.omit)")) %>%
		mutate(combo = paste0(MicroFeat, " ~ ", SNP)) %>%
		as.data.frame()
	
	DF <- list()
	
	for (i in levels(as.factor(sig_fe$combo))){
		
		species <- gsub(" .*", "", i)
		
		snp <- gsub(".* ", "", i)
		
		my_df <- data %>%
			select(all_of(md_list), all_of(species), all_of(snp)) %>%
			rename(Abundance=5, Allele=6) %>%
			mutate(MicroFeat=species) %>%
			mutate(SNP=snp) %>%
			select(Subject, MicroFeat, SNP, Allele) %>%
			distinct() %>%
			drop_na() %>%
			group_by(MicroFeat, SNP, Allele) %>%
			tally() %>%
			filter(n == min(n)) %>%
			filter(n >= 5) %>%
			ungroup() %>%
			group_by(MicroFeat) %>%
			filter(n == max(n))
		
		DF <- rbind(DF, my_df)
		
	}
	
	DF2 <- merge(sig_fe, DF, by=c("MicroFeat", "SNP")) %>%
		group_by(MicroFeat) %>%
		filter(P == min(P)) %>%
		filter(n == max(n))
	
	if (nrow(DF2) > 0) {

		fnc <- nlevels(as.factor(DF2$combo))

		fh <- 5 * ceiling(fnc/5)

		if (fnc < 5) {fw <- 5 * fnc} else {fw <- 25}
		
		allele_map <- data %>%
			reshape2::melt(id=c(md_list, micro_list), variable.name="SNP", value.name="Allele") %>%
			reshape2::melt(id=c(md_list, "SNP", "Allele"), variable.name="MicroFeat", value.name="Abundance") %>%
			mutate(combo = paste0(MicroFeat, " ~ ", SNP)) %>%
			filter(combo %in% DF2$combo) %>%
			drop_na() %>%
			group_by(combo, SNP, Allele) %>%
			tally() %>%
			ungroup() %>%
			mutate(A1 = str_sub(Allele, 1, 1)) %>%
			mutate(A2 = str_sub(Allele, 2, 2)) %>%
			mutate(Zygosity = case_when(A1 == A2 ~ "Homozygous", A1 != A2 ~ "Heterozygous")) %>%
			group_by(SNP) %>%
			mutate(temp = case_when(n == min(n) ~ "minor", n == max(n) ~ "major", n != min(n) & n != max(n) ~ "")) %>%
			mutate(AlleleType = str_trim(paste0(Zygosity, " ", temp))) %>%
			select(combo, SNP, Allele, AlleleType) %>%
			mutate(Order = case_when(
				AlleleType == "Homozygous major" ~ 1,
				AlleleType == "Heterozygous" ~ 2,
				AlleleType == "Homozygous minor" ~ 3)) %>%
			ungroup() %>%
			arrange(combo, SNP, Order) %>%
			mutate(X = paste0(SNP, "_", Order))

		mdf <- data %>%
			reshape2::melt(id=c(md_list, micro_list), variable.name="SNP", value.name="Allele") %>%
			reshape2::melt(id=c(md_list, "SNP", "Allele"), variable.name="MicroFeat", value.name="Abundance") %>%
			mutate(combo = paste0(MicroFeat, " ~ ", SNP)) %>%
			filter(combo %in% DF2$combo) %>%
			drop_na() %>%
			merge(., allele_map, by=c("combo", "SNP", "Allele"))

		# boxplot

		fixed_plot <- ggplot(data=mdf, aes(x=X, y=Abundance)) +
			facet_wrap(~combo, scales="free", ncol=5) +
			geom_jitter(aes(colour=AlleleType, fill=AlleleType), pch=21, alpha=0.5, width=0.1) +
			geom_boxplot(colour="black", fill=NA, outlier.shape=NA) +
			theme_bw() +
			theme(
				legend.position = "top",
				aspect.ratio = 0.75,
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=12.5),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=10, face="bold")) +
			scale_colour_manual(values = my_pal) +
			scale_fill_manual(values = my_pal) +
			scale_x_discrete(labels = unique(allele_map$Allele))

		ggsave(plot=fixed_plot, paste0(opts$o, "/diagnostic_fe_box_", opts$a, ".png"), dpi=300, height=fh, width=fw)
		
		# boxplot (average per subject)

		
		mdf_2 <- mdf %>%
			group_by(Subject, MicroFeat, SNP, X, Allele, AlleleType, combo) %>%
			summarise(Abundance = mean(Abundance)) %>%
			group_by(combo, SNP, Allele, AlleleType) %>%
			mutate(AlleleType = factor(AlleleType, levels=c("Homozygous major", "Heterozygous", "Homozygous minor")))
		
		custom_labels <- allele_map$Allele
		names(custom_labels) <- allele_map$X
		
		fixed_plot_sub <- ggplot(data=mdf_2, aes(x=X, y=Abundance)) +
			scale_x_discrete(labels = custom_labels) +
			facet_wrap(~combo, scales="free", ncol=5) +
			geom_jitter(aes(colour=AlleleType, fill=AlleleType), pch=21, alpha=0.5, width=0.25) +
			geom_boxplot(colour="black", fill=NA, outlier.shape=NA) +
			theme_bw() +
			theme(
				legend.position = "top",
				aspect.ratio = 0.75,
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=12.5),
				legend.title = element_blank(),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=10, face="bold")) +
			scale_colour_manual(values = my_pal) +
			scale_fill_manual(values = my_pal)

		ggsave(plot=fixed_plot_sub, paste0(opts$o, "/diagnostic_fe_box_avg_per_sub_", opts$a, ".png"), dpi=300, height=fh, width=fw)

		# line plot

		fixed_plot_line <- ggplot(data=mdf, aes(x=Days, y=Abundance)) +
			facet_wrap(~combo, scales="free", ncol=5) +
			geom_smooth(aes(colour=AlleleType, fill=AlleleType), method="lm") +
			theme_bw() +
			theme(
				legend.position = "top",
				aspect.ratio = 0.75,
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=12.5),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=10, face="bold"))

		ggsave(plot=fixed_plot_line, paste0(opts$o, "/diagnostic_fe_line_", opts$a, ".png"), dpi=300, height=fh, width=fw)

		# zero v non-zero

		fixed_zeros <- mdf %>%
			mutate(Zero = case_when(
				Abundance == min(Abundance) ~ "zero",
				Abundance > min(Abundance) ~ "non_zero")) %>%
			group_by(MicroFeat, SNP, Allele, Zero) %>%
			tally() %>%
			ungroup() %>%
			reshape2::dcast(MicroFeat + SNP + Allele ~ Zero, value.var="n") %>%
			mutate(Percent = 100 * (non_zero / (non_zero + zero)))
	
		fixed_plot_zeros <- ggplot(data=fixed_zeros, aes(x=AlleleType, y=Percent)) +
			scale_x_discrete(labels = unique(allele_map$Allele)) +
			facet_wrap(MicroFeat ~ SNP, scales="free", ncol=5) +
			geom_bar(stat="identity", aes(fill=AlleleType), colour="black") +
			theme_bw() +
			theme(
				legend.position = "top",
				aspect.ratio = 0.75,
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=12.5),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=10, face="bold")) +
			scale_colour_manual(values = my_pal) +
			scale_fill_manual(values = my_pal)

		ggsave(plot=fixed_plot_zeros, paste0(opts$o, "/diagnostic_fe_zeros_bar_", opts$a, ".png"), dpi=300, height=fh, width=fw)
	
		}
	
	################
	# interactions #
	################
	
	# interaction effect

	sig_int <- sig_micro %>%
		filter(Type == "Interaction") %>%
		group_by(MicroFeat) %>%
		filter(P == min(P)) %>%
		ungroup() %>%
		distinct() %>%
		mutate(Model = paste0("lmer(", MicroFeat, " ~ ", SNP, " * ", "Day + Country + ( 1 + Day | Subject ), data=data, na.action=na.omit)")) %>%
		as.data.frame()

	sig_int_map <- sig_int %>%
		mutate(combo = paste0(MicroFeat, " ~ ", SNP)) %>%
		select(combo)
	
	if (nrow(sig_int_map) > 0) {

		inc <- nlevels(as.factor(sig_int_map$combo))

		ih <- 5 * ceiling(fnc/5)

		if (inc < 5) {iw <- 5 * inc} else {iw <- 25}

		mdi <- data %>%
			reshape2::melt(id=c(md_list, micro_list), variable.name="SNP", value.name="Allele") %>%
			reshape2::melt(id=c(md_list, "SNP", "Allele"), variable.name="MicroFeat", value.name="Abundance") %>%
			mutate(combo = paste0(MicroFeat, " ~ ", SNP)) %>%
			filter(combo %in% sig_int_map$combo) %>%
			drop_na()
		
		# interaction line
		
		interaction_plot <- ggplot(data=mdi, aes(x=Days, y=Abundance)) +
			facet_wrap(MicroFeat ~ SNP, scales="free", ncol=5) +
			geom_smooth(aes(colour=Allele, fill=Allele), method="lm") +
			theme_bw() +
			theme(
				legend.position = "top",
				aspect.ratio = 0.75,
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=12.5),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=12.5, face="bold"))

		ggsave(plot=interaction_plot, paste0(opts$o, "/diagnostic_ie_line_", opts$a, ".png"), dpi=300, height=ih, width=iw)

		# slope

		library(purrr)
		library(broom)		

		slope_df <- mdi %>%
			split(., list(.$Subject, .$MicroFeat), drop = TRUE) %>%
			map(~lm(Abundance ~ Days, data = .x)) %>% 
			map_df(tidy, .id = "Subject.MicroFeat") %>%
			filter(term == "Days") %>%
			mutate(Subject = gsub("\\..*", "", Subject.MicroFeat)) %>%
			mutate(MicroFeat = gsub(".*\\.", "", Subject.MicroFeat)) %>%
			select(-1, -2) %>%
			merge(., mdi, by=c("Subject", "MicroFeat")) %>%
			drop_na() %>%
			select(-sample_id, -Days, -Abundance) %>%
			distinct()
	
		bugs <- slope_df %>%
			group_by(MicroFeat, SNP, Allele) %>%
			tally() %>%
			filter(n == min(n)) %>%
			filter(n > 10) %>%
			ungroup() %>%
			select(MicroFeat)
		
		sdf <- slope_df %>%
			filter(MicroFeat %in% bugs$MicroFeat)
		
		if (nrow(sdf) > 0) {
			
			snc <- nlevels(as.factor(sdf$combo))
		
			sh <- 5 * ceiling(snc/5)

			if (snc < 5) {sw <- 5 * snc} else {sw <- 25}

			# slope density plot

			slope_plot <- ggplot(sdf, aes(x=estimate)) +
				facet_wrap(~combo, scales="free", ncol=5) +
				geom_density(aes(colour=Allele), size=1) +
				theme_bw() +
				theme(
					legend.position = "top",
					aspect.ratio = 0.75,
					plot.title = element_text(size=12.5, face="bold", hjust=0.5),
					axis.title = element_text(size=12.5, face="bold"),
					axis.text = element_text(size=12.5),
					legend.title = element_text(size=12.5, face="bold", vjust=0.75),
					legend.text = element_text(size=12.5),
					strip.text = element_text(size=10, face="bold")) +
				labs(x="Slope", y="Density")

			ggsave(plot=slope_plot, paste0(opts$o, "/diagnostic_ie_slope_density_", opts$a, ".png"), dpi=300, height=sh, width=sw)

			# slope boxplot

			slope_violin <- ggplot(data=sdf, aes(x=Allele, y=estimate)) +
				facet_wrap(~combo, scales="free", ncol=5) +
				geom_violin(aes(fill=Allele), colour="black", outlier.shape=NA) +
				theme_bw() +
				theme(
					legend.position = "top",
					aspect.ratio = 0.75,
					plot.title = element_text(size=12.5, face="bold", hjust=0.5),
					axis.title = element_text(size=12.5, face="bold"),
					axis.text = element_text(size=12.5),
					legend.title = element_text(size=12.5, face="bold", vjust=0.75),
					legend.text = element_text(size=12.5),
					strip.text = element_text(size=10, face="bold")) +
				labs(x="Allele", y="Slope")

			ggsave(plot=slope_violin, paste0(opts$o, "/diagnostic_ie_slope_violin_", opts$a, ".png"), dpi=300, height=sh, width=sw)
		
			}
		
		}
		
	}

#
