#!/usr/bin/env Rscript

require(docopt)

'Usage:
   plink_plots.R [--ki <txt> --mi <frq> --mt <maf_min> --ei <hwe> --et <hwe_min> --gi <genome> --gt <pihat_max> --hi <het> --pc <pca> --md <metadata> --out <output>]

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
   --md metadata for samples (.tsv)
   --out output directory [default: plink_output/]

 ' -> doc

opts <- docopt(doc)

#############
# libraries #
#############

library(tidyverse)
library(reshape2)
library(scales)
library(data.table)
library(ggridges)
library(wesanderson)
library(cowplot)

prefix <- gsub(".*\\/", "", opts$pc)

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

opts$mt <- as.numeric(opts$mt)

maf <- read.table(opts$mi, header=T, as.is=T)

# maf hist

maf_hist <- ggplot(maf, aes(x=MAF)) +
	geom_histogram(fill="goldenrod", colour="black") +
	theme_bw() +
	theme(panel.grid=element_blank(),
		plot.title=element_text(size=15, face="bold", hjust=0.5),
		axis.title=element_text(size=12.5, face="bold"),
		axis.text=element_text(size=12.5),
		legend.position="none") +
	labs(title="Minor allele frequency distribution", 
		x="Minor allele frequency", 
		y="Number of SNPs") +
	scale_y_continuous(expand = expansion(mult = c(0, .1)),	labels=comma) +
	geom_vline(xintercept=opts$mt, linetype="dashed", colour="firebrick", size=0.75) +
	annotate(geom="label", label=paste0("threshold \U2265 ", opts$mt), 
		x=Inf, y=Inf, size=5, hjust=1, vjust=1, label.size = 0)

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

opts$et <- as.numeric(opts$et)

hardy <- read.table(opts$ei, header=T, as.is=T) %>%
	mutate(P_log = -log10(P))

# hwe hist

hwe_hist <- ggplot(hardy, aes(x=P_log)) +
	geom_histogram(fill="goldenrod", colour="black") +
	theme_bw() +
	theme(panel.grid=element_blank(),
		plot.title=element_text(size=15, face="bold", hjust=0.5),
		axis.title=element_text(size=12.5, face="bold"),
		axis.text=element_text(size=12.5),
		legend.position="none") +
	labs(title="The p-value of the exact test for HWE", 
		x=expression(bold(-log[10]~p-value)), 
		y="Number of SNPs") +
	scale_y_continuous(expand = expansion(mult = c(0, .1)), label=comma) +
	geom_vline(xintercept=-log10(as.numeric(opts$et)), linetype="dashed", colour="firebrick", size=0.75) +
	annotate(geom="label", label=paste0("threshold \U2264 \U2212log(", opts$et, ")"),
		x=Inf, y=Inf, size=5, hjust=1, vjust=1, label.size = 0)

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
	geom_vline(xintercept=mint, linetype="dashed", colour="firebrick", size=0.75) +
	geom_vline(xintercept=maxt, linetype="dashed", colour="firebrick", size=0.75) +
	annotate(geom="label", label=paste0("threshold = mean ", "\U00B1", " 3 * sd"), 
		x=Inf, y=Inf, size=5, hjust=1, vjust=1, label.size = 0)
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
	geom_vline(xintercept=opts$gt, linetype="dashed", colour="firebrick", size=0.75) +
	annotate(geom="label", label=paste0("threshold \U2264 ", opts$gt),
		x=Inf, y=Inf, size=5, hjust=1, vjust=1, label.size = 0)

# ibd summary

ibd_fail <- nrow(filter(genome, PI_HAT > opts$gt))

ibd_summary <- data.frame(Stage="3: IBD", Pass=(het_summary[1,2] - ibd_fail), Fail=ibd_fail)

#####################
# identity by state #
#####################

# calculate the median distance of each subject to every other subject

dst <- genome %>%
	group_by(IID1) %>%
	summarise(DST_median = median(DST)) %>%
	ungroup()

# define outliers

outliers <- dst %>%
		filter(DST_median < quantile(DST_median, probs=c(.25, .75), na.rm = FALSE)[1] - 1.5 * IQR(DST_median))

#copy input

dst_in <- dst

# remove outliers iteratively

while (nrow(outliers) > 0) {

	dst <- dst %>%
		filter(DST_median > quantile(DST_median, probs=c(.25, .75), na.rm = FALSE)[1] - 1.5 * IQR(DST_median))

	outliers <- dst %>%
		filter(DST_median < quantile(DST_median, probs=c(.25, .75), na.rm = FALSE)[1] - 1.5 * IQR(DST_median))
	
	}

subjects <- dst %>%
	select(1)

# summary

ibs_summary <- data.frame( Stage = "4: IBS", Pass = nrow(subjects), Fail = ibd_summary[1,2] - nrow(subjects))

# plot

ibs_vline <- quantile(dst$DST_median, probs=c(.25, .75), na.rm = FALSE)[1] - 1.5 * IQR(dst$DST_median)

ibs_hist <- ggplot(dst_in, aes(x=DST_median)) +
 	geom_histogram(fill="steelblue", colour="black") +
 	scale_y_continuous(expand = expansion(mult = c(0, .1))) +
 	theme_bw() +
 	theme(panel.grid=element_blank(),
 		plot.title=element_text(size=15, face="bold", hjust=0.5),
 		axis.title=element_text(size=12.5, face="bold"),
 		axis.text=element_text(size=12.5),
 		strip.text=element_text(size=12.5, face="bold"),
 		legend.title=element_text(size=12.5, face="bold"),
 		legend.text=element_text(size=12.5),
 		legend.position="none") +
 	ggtitle("IBS for all samples") +
 	labs(x="Median identity of individual to others", y="Number of subjects") +
 	geom_vline(xintercept=ibs_vline, linetype="dashed", colour="firebrick", size=0.75) +
 	annotate(geom="label", label=paste0("threshold \U2265 Q1 \U2212 1.5 * IQR"),
		x=Inf, y=Inf, size=5, hjust=1, vjust=1, label.size = 0)

################################
# Principal Component Analysis #
################################

# eigenvec

eigenvec <- fread(paste0(opts$pc, ".qc.excl_outliers.eigenvec")) %>%
	rename(subject_id=IID) %>%
	select(2:7)

# eigenval

eigenval <- fread(paste0(opts$pc, ".qc.excl_outliers.eigenval")) %>%
	mutate(pc = row_number()) %>%
	rename(eigenval=1) %>%
	mutate(pc_eigenval = paste0("PC", pc, " [", round(eigenval, 2), "%]"))

eigenval_2 <- eigenval %>%
	mutate(pc = paste0("PC", pc))

# metadata

metadata <- fread(opts$md) %>%
	select(subject_id, country) %>%
	distinct()

# combine eigenvec + metadata + eigenval

df <- merge(eigenvec, metadata, by="subject_id") %>%
	melt(id=c("subject_id", "country"), variable.name="pc", value.name="value") %>%
	merge(., eigenval_2, by="pc")

# scree plot

scree_plot <- ggplot(eigenval, aes(x=as.factor(pc), y=eigenval)) +
	geom_bar(stat="identity", colour="black", fill="skyblue") +
	geom_line(aes(x=as.numeric(pc), y=eigenval), colour="firebrick") +
	geom_point(pch=21, colour="black", fill="firebrick") +
	theme_bw() +
	theme(plot.title=element_text(face="bold", size=15, hjust=0.5),
         	axis.title = element_text(face="bold", size=12.5),
         	axis.text = element_text(size=12.5),
         	panel.grid = element_blank()) +
	labs(title="PLINK eigenval", x="PC", y="%") +
	scale_y_continuous(expand = expansion(mult = c(0, .1)))

# scatter plot

pca_plot <- ggplot(eigenvec, aes(x=PC1, y=PC2)) +
	geom_point(pch=21, colour="black", fill="firebrick") +
	theme_bw() +
	theme(plot.title=element_text(face="bold", size=15, hjust=0.5),
		axis.title = element_text(face="bold", size=12.5),
		axis.text = element_text(size=12.5)) +
	labs(title="PLINK eigenvec", 
		x=paste0("PC1 ", "[", round(eigenval[1,1], 2), "%]"),
		y=paste0("PC2 ", "[", round(eigenval[2,1], 2), "%]"))

# ridge plot

pc_ridges <- ggplot(df, aes(x=value, y=pc_eigenval)) +
	geom_density_ridges2(aes(fill=country), colour="black") +
	theme_bw() +
	scale_fill_manual(values=rev(wes_palette("GrandBudapest2"))) +
	theme(legend.position = "top") +
	labs(x="Eigenvector", y="", fill="Country") +
	scale_x_continuous(expand = expansion(mult = c(0, 0))) +
	scale_y_discrete(expand = expansion(mult = c(0, .1))) +
	theme(plot.title=element_text(size=15, face="bold", hjust=0.5),
		axis.title=element_text(size=12.5, face="bold"),
		axis.text=element_text(size=12.5),
		strip.text=element_text(size=12.5, face="bold"),
		legend.title=element_text(size=12.5, face="bold"),
		legend.text=element_text(size=12.5))

# combine PCA plots

pca_fig_ab <- plot_grid(scree_plot, pca_plot, nrow=2, align="v", scale=0.95)

pca_fig_abc <- plot_grid(pca_fig_ab, pc_ridges, nrow=1, align="h", axis="tblr", scale=0.95)

ggsave(plot=pca_fig_abc,
	file=paste0(opts$out, prefix, ".PLINK_PCA_fig.png"),
	height=10, width=15,
	dpi=300)

#########################
# combine the summaries #
#########################

snps_summary <- rbind(king_summary_snps, maf_summary, hardy_summary) %>%
	mutate(Type=c("SNPs"))

subs_summary <- rbind(king_summary_subjects, het_summary, ibd_summary, ibs_summary) %>%
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
	geom_text(position=position_dodge(width=0.9), vjust=-0.25, size=4) +
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

gp <- ggplotGrob(summary_bar)

facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]

x.var <- sapply(ggplot_build(summary_bar)$layout$panel_scales_x,
                function(l) length(l$range$range))

gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var

ggsave(plot=gp, paste0(opts$out, prefix, ".QC_bar_summary.png"),
	height=5, width=12.5,
	dpi=300)

#######################################
# combine the histograms with cowplot #
#######################################

top_row <- plot_grid(NULL, maf_hist, hwe_hist, NULL, ncol=4, rel_widths=c(1/6, 1/3, 1/3, 1/6))

bottom_row <- plot_grid(het_hist, ibd_hist, ibs_hist, ncol=3)

combined_hist <- plot_grid(top_row, bottom_row, nrow=2, scale=0.95)

ggsave(plot=combined_hist, paste0(opts$out, prefix, ".QC_histograms.png"),
	height=10, width=20,
	dpi=300)

#
