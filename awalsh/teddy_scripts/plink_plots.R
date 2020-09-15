#!/usr/bin/env Rscript

require(docopt)

'Usage:
   plink_plots.R [--ki <txt> --mi <frq> --mt <maf_min> --ei <hwe> --et <hwe_min> --gi <genome> --gt <pihat_max> --hi <het> --pc <pca>]

Options:
   --ki output from KING autoQC (_king_autoQC_Summary.txt)
   --mi minor allele frequencies (.frq)
   --mt threshold for MAF
   --ei Hardy-Weinberg equilibrium values (.hwe)
   --et threshold for Hardy-Weinberg equilibrium
   --gi PIHAT values (.genome)
   --gt threshold for PI_HAT
   --hi heterozygosity values (.het)
   --pc prefix of pca files

 ' -> doc

opts <- docopt(doc)

#############
# libraries #
#############

library(tidyverse)
library(reshape2)
library(scales)
library(cowplot)

###############
# KING autoQC #
###############

king_summary <- read.csv(opts$ki, sep="\t", header=T)

# KING SNPs summary

king_summary_snps <- king_summary %>%
	filter(grepl("SNPs", Type)) %>%
	mutate(Data = gsub("QC'ed", "Pass", Data)) %>%
	select(Data, Number) %>%
	column_to_rownames("Data") %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("Stage") %>%
	mutate(Fail = Raw - Pass) %>%
	mutate(Stage = gsub("Number", "1: KING", Stage)) %>%
	select(-Raw)

# KING Subjects summary

king_summary_subjects <- king_summary %>%
	filter(grepl("Subjects", Type)) %>%
	mutate(Data = gsub("QC'ed", "Pass", Data)) %>%
	select(Data, Number) %>%
	column_to_rownames("Data") %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("Stage") %>%
	mutate(Fail = Raw - Pass) %>%
	mutate(Stage = gsub("Number", "1: KING", Stage)) %>%
	select(-Raw)

############################
# minor allele frequencies #
############################

maf <- read.table(opts$mi, header=T, as.is=T)

opts$mt <- as.numeric(opts$mt)

# maf hist

maf_hist <- ggplot(maf, aes(x=MAF)) +
	geom_histogram(fill="goldenrod", colour="black") +
	theme_bw() +
	theme(panel.grid=element_blank(),
		plot.title=element_text(size=15, face="bold", hjust=0.5),
		axis.title=element_text(size=12.5, face="bold"),
		axis.text=element_text(size=12.5),
		strip.text=element_text(size=12.5, face="bold"),
		legend.title=element_text(size=12.5, face="bold"),
		legend.text=element_text(size=12.5)) +
	labs(title="Minor allele frequency distribution", 
		x="Minor allele frequency", 
		y="Number of SNPs") +
	scale_y_continuous(expand = expansion(mult = c(0, .1)),	labels=comma) +
	geom_vline(xintercept=opts$mt, linetype="dashed", colour="firebrick")

# maf summary

maf_summary <- maf %>%
	mutate(Status = case_when(
		MAF < opts$mt ~ "Fail", MAF >= opts$mt ~ "Pass")) %>%
	group_by(Status) %>%
	tally() %>%
	ungroup() %>%
	column_to_rownames("Status") %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("Stage") %>%
	mutate(Stage = gsub("n", "2: MAF", Stage)) %>%
	select(Stage, Pass, Fail)

##############################
# Hardy-Weinberg equilibrium #
##############################

hardy <- read.table(opts$ei, header=T, as.is=T) %>%
	mutate(P_log = -log10(P))

opts$et <- as.numeric(opts$et)

# hwe hist

hwe_hist <- ggplot(hardy, aes(x=P_log)) +
	geom_histogram(fill="goldenrod", colour="black") +
	theme_bw() +
	theme(panel.grid=element_blank(),
		plot.title=element_text(size=15, face="bold", hjust=0.5),
		axis.title=element_text(size=12.5, face="bold"),
		axis.text=element_text(size=12.5),
		strip.text=element_text(size=12.5, face="bold"),
		legend.title=element_text(size=12.5, face="bold"),
		legend.text=element_text(size=12.5)) +
	labs(title="The p-value of the exact test for HWE", 
		x=expression(bold(-log[10]~p-value)), 
		y="Number of SNPs") +
	scale_y_continuous(expand = expansion(mult = c(0, .1)), label=comma) +
	geom_vline(xintercept=-log10(as.numeric(opts$et)), linetype="dashed", colour="firebrick")

# hwe summary

hardy_summary <- hardy %>%
	mutate(Status = case_when(
		P < opts$et ~ "Fail", P >= opts$et ~ "Pass")) %>%
	group_by(Status) %>%
	tally() %>%
	ungroup() %>%
	column_to_rownames("Status") %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("Stage") %>%
	mutate(Stage = gsub("n", "3: HWE", Stage)) %>%
	select(Stage, Pass, Fail)

##################
# heterozygosity #
##################

het <- read.table(opts$hi, header=T) %>%
	mutate(HET_RATE = (N.NM. - O.HOM.)/N.NM.)

mint <- mean(het$HET_RATE) - 3 * sd(het$HET_RATE)
maxt <- mean(het$HET_RATE) + 3 * sd(het$HET_RATE)

# heterozygosity hist

het_hist <- ggplot(het, aes(x=HET_RATE)) +
	geom_histogram(colour="black", fill="steelblue") +
	theme_bw() +
	theme(panel.grid=element_blank(),
		plot.title=element_text(size=15, face="bold", hjust=0.5),
		axis.title=element_text(size=12.5, face="bold"),
		axis.text=element_text(size=12.5),
		strip.text=element_text(size=12.5, face="bold"),
		legend.title=element_text(size=12.5, face="bold"),
		legend.text=element_text(size=12.5)) +
	labs(title="The rate of heterozygosity for samples", 
		x="Rate of heterozygosity", 
		y="Number of subjects") +
	scale_y_continuous(expand = expansion(mult = c(0, .1)),	labels=comma) +
	geom_vline(xintercept=mint, linetype="dashed", colour="firebrick") +
	geom_vline(xintercept=maxt, linetype="dashed", colour="firebrick")

# heterozygosity summary

fail <- het %>%
	filter(HET_RATE < mint | HET_RATE > maxt)

pass <- het %>%
	filter(HET_RATE > mint & HET_RATE < maxt)

het_summary <- data.frame(Stage="2: Heterozygosity", Pass=nrow(pass), Fail=nrow(fail))

#######################
# identity by descent #
#######################

genome <- read.table(opts$gi, header=T, as.is=T)

opts$gt <- as.numeric(opts$gt)

# ibd hist

ibd_hist <- ggplot(genome, aes(x=PI_HAT)) +
	geom_histogram(colour="black", fill="steelblue") +
	theme_bw() +
	theme(panel.grid=element_blank(),
		plot.title=element_text(size=15, face="bold", hjust=0.5),
		axis.title=element_text(size=12.5, face="bold"),
		axis.text=element_text(size=12.5),
		strip.text=element_text(size=12.5, face="bold"),
		legend.title=element_text(size=12.5, face="bold"),
		legend.text=element_text(size=12.5)) +
	labs(title="IBD for all sample pairs", 
		x="Estimate pairwise IBD (PIHAT)", 
		y="Number of pairs") +
	scale_y_continuous(expand = expansion(mult = c(0, .1)), labels=comma) +
	geom_vline(xintercept=opts$gt, linetype="dashed", colour="firebrick")

# ibd summary

ibd_fail <- nrow(filter(genome, PI_HAT > opts$gt))

ibd_summary <- data.frame(Stage="3: IBD", Pass=(het_summary[1,2] - ibd_fail), Fail=ibd_fail)

################################
# Principal Component Analysis #
################################

# eigenval

eigenval_eo <- read.table(paste0(opts$pc, ".excl_outliers.eigenval"), header=F) %>%
	mutate(PC=row_number())  %>%
	mutate(Outliers="Excluded")

scree_plot <- ggplot(eigenval_eo, aes(x=as.factor(PC), y=V1)) +
	geom_bar(stat="identity", colour="black", fill="skyblue") +
	geom_line(aes(x=as.numeric(PC), y=V1), colour="firebrick") +
	geom_point(pch=21, colour="black", fill="firebrick") +
	theme_bw() +
	theme(plot.title=element_text(face="bold", size=15, hjust=0.5),
         	axis.title = element_text(face="bold", size=12.5),
         	axis.text = element_text(size=12.5),
         	panel.grid = element_blank()) +
	labs(title="PLINK eigenval", x="PC", y="%") +
	scale_y_continuous(expand = expansion(mult = c(0, .1)))

# eigenvec

eigenvec_io <- read.table(paste0(opts$pc, ".incl_outliers.eigenvec"), header=T)

eigenvec_eo <- read.table(paste0(opts$pc, ".excl_outliers.eigenvec"), header=T) %>%
	select(-FID) %>%
	mutate(Outliers="Excl. outliers")
	
pca_summary <- data.frame(Stage = "PCA",
	Pass = nrow(eigenvec_eo),
	Fail= nrow(eigenvec_io) - nrow(eigenvec_eo))

# scatterplot

pca_plot <- ggplot(eigenvec_eo, aes(x=PC1, y=PC2)) +
	geom_point(pch=21, colour="black", fill="firebrick") +
	theme_bw() +
	theme(plot.title=element_text(face="bold", size=15, hjust=0.5),
		axis.title = element_text(face="bold", size=12.5),
		axis.text = element_text(size=12.5)) +
	labs(title="PLINK eigenvec", 
		x=paste0("PC1 ", round(eigenval_eo[1,1], 2), "%"),
		y=paste0("PC1 ", round(eigenval_eo[2,1], 2), "%"))

# combined pca figure

pca_fig <- plot_grid(scree_plot, pca_plot, rel_widths=c(1.25,1), scale=0.9)

ggsave("PLINK_PCA_fig.png",
	height=5, width=12.5,
	dpi=300)

# pca summary

eigenvec_io <- read.table(paste0(opts$pc, ".incl_outliers.eigenvec"), header=T)

pca_summary <- data.frame(Stage = "4: PCA",
	Pass = nrow(eigenvec_eo),
	Fail= nrow(eigenvec_io) - nrow(eigenvec_eo))

#########################
# combine the summaries #
#########################

snps_summary <- rbind(king_summary_snps, maf_summary, hardy_summary) %>%
	mutate(Type=c("SNPs"))

subs_summary <- rbind(king_summary_subjects, het_summary, ibd_summary, pca_summary) %>%
	mutate(Type=c("Subjects"))

combined_summary <- rbind(snps_summary, subs_summary) %>%
	melt(id=c("Stage", "Type"), variable.name="QC", value.name="Number") %>%
	mutate(QC = factor(QC, levels=c("Pass", "Fail")))

summary_bar <- ggplot(combined_summary, 
		aes(x=Stage, y=Number, fill=QC, label=comma(Number, accuracy=1))) +
	geom_bar(stat="identity",
		position=position_dodge(width=0.9),
		color="black") +
	facet_wrap(~Type, scales="free") +
	geom_text(position=position_dodge(width=0.9), vjust=-0.25, size=5) +
	theme_bw() +
	theme(panel.grid=element_blank(),
		plot.title=element_text(size=15, face="bold", hjust=0.5),
		axis.title=element_text(size=12.5, face="bold"),
		axis.text=element_text(size=12.5),
		strip.text=element_text(size=12.5, face="bold"),
		legend.title=element_text(size=12.5, face="bold"),
		legend.text=element_text(size=12.5)) +
	scale_y_continuous(expand = expansion(mult = c(0, .1)),  labels=comma) +
	scale_fill_manual(values=c("Pass"="forestgreen", "Fail"="grey")) +
	ggtitle("QC summary") +
	labs(x="QC step", y="Number")
summary_bar

ggsave("QC_bar_summary.png",
	height=5, width=12.5,
	dpi=300)

#######################################
# combine the histograms with cowplot #
#######################################

# histograms

combined_hist <- plot_grid(
	maf_hist, hwe_hist,
	ibd_hist, het_hist,
	nrow=2,
	scale=0.9,
	align="hv", axis="tblr")
combined_hist

ggsave("QC_histograms.png",
	height=10, width=15,
	dpi=300)

#############################
# PCA plot (before v after) #
#############################

