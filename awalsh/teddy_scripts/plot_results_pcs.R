#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plot_results_pcs.R [-s <stats> -d <data> -o <out_dir> -p <padj_method> -g <group_by> -t <trans> -f <feature> --humann <humann>]

Options:
   -s stats (output from run_models.R)
   -d data
   -o output directory for plots [default: out_dir]
   -p method for p-value adjustment [default: bonferroni]
   -g for p-value adjustment, group by Type or Predictor [default: Type] 
   -t transformation method
   -f microbiome feature (Bray-Curtis, EC, Pfam, Species)
   --humann HUMAnN2 pathway names

' -> doc

opts <- docopt(doc)

#

library(tidyverse)
library(reshape2)
library(ggridges)
library(viridis)
library(scales)
library(lmerTest)
library(cowplot)

#################
# load the data #
#################

# data (input for models)

data <- read.csv(opts$d, sep="\t", header=T) %>%
	setNames(gsub("microbiome\\.|genetics_", "", names(.))) %>%
	setNames(gsub("\\.", "_", names(.)))

colnames(data) <- gsub("microbiome\\.|genetics_", "", colnames(data))

############################
# significant associations #
############################

# input

model_stats <- read.csv(opts$s, sep="\t", header=T, check.names=F)

# output directory

dir.create(opts$o)

######################
# p-value adjustment #
######################

df_sig <- model_stats %>%
	filter(grepl("PC", Predictor)) %>%
	mutate(Type = case_when(
		grepl(":", Predictor) ~ "PC:Days",
		!grepl(":", Predictor) ~ "PC")) %>%
	group_by(!!sym(opts$g)) %>%
	mutate(P_adj = p.adjust(P, method=opts$p)) %>%
	ungroup() %>%
	filter(P_adj < 0.05) %>%
	mutate(Coefficient = case_when((Coefficient > 0) ~ 1, (Coefficient < 0) ~ -1)) %>%
	mutate(value = Coefficient * -log(P_adj)) %>%
	mutate(MicroFeat = gsub("lmer\\(", "", Model)) %>%
	mutate(MicroFeat = gsub(" .*", "", MicroFeat)) %>%
	mutate(MicroFeat = gsub("_noname_unclassified", "_?_?", MicroFeat)) %>%
	mutate(MicroFeat = gsub("_unclassified", "_?", MicroFeat))  %>%
	mutate(Model = gsub("genetics_", "", Model)) %>%
	mutate(PCN = gsub("genetics_|\\:.*", "", Predictor)) %>%
	mutate(Predictor = gsub("genetics_", "", Predictor)) %>%
	droplevels() %>%
	as.data.frame()

# label

if (opts$f == "Bray-Curtis"){

	my_label <- "Bray-Curtis Distance"

	}
	
if (opts$f == "EC" | opts$f == "Pfam"){

	my_label <- "log(Abundance [CPM])"

	}

if (opts$f == "Species"){

	my_label <- "log(Abundance [%])"

	}
	
#
	
if (nrow(df_sig) > 0) {

	#############
	# add names #
	#############

	# EC names

	if(grepl("EC_", df_sig$MicroFeat[1]) == TRUE){

		pathways <- read.csv(opts$humann, sep="\t", header=F) %>%
		rename(PathID=1, Name=2) %>%
		mutate(MicroFeat = paste0("EC_", gsub("\\.", "_", PathID)))

		df_sig <- merge(df_sig, pathways, by="MicroFeat") %>%
			select(-MicroFeat) %>%
			mutate(MicroFeat = paste0("EC ", PathID, ": ", Name)) %>%
			droplevels()

	}

	# Pfam names

	if(grepl("Pfam_", df_sig$MicroFeat[1]) == TRUE){

		pathways <- read.csv(opts$humann, sep="\t", header=F) %>%
			rename(PathID=1, Name=2) %>%
			mutate(MicroFeat = paste0("Pfam_", gsub("\\.", "_", PathID)))

		df_sig <- merge(df_sig, pathways, by="MicroFeat") %>%
			select(-MicroFeat) %>%
			mutate(MicroFeat = paste0(PathID, ": ", Name)) %>%
			droplevels()
	
	}

	write.table(df_sig, file=paste0("./", opts$o, "/model_stats.tsv"), sep="\t", quote=F, row.names=F)

	########################################
	# visualization of significant results #
	########################################

	###########
	# heatmap #
	###########

	if (opts$f != "Bray-Curtis"){

	heatmap <- ggplot(data=df_sig, aes(x=Predictor, y=MicroFeat, fill=value)) +
		geom_tile(colour="black", size=0.125) +
		facet_grid(. ~ Type, space="free", scales="free") +
		scale_fill_gradient2(low = "steelblue", 
			mid = "white",
			high = "darksalmon") +
		scale_x_discrete(expand=c(0, 0)) +
		scale_y_discrete(expand=c(0, 0), position = "right") +
		theme_bw() +
		theme(
			axis.title = element_text(size=10, face="bold"),
			axis.text.x = element_text(size=10),
			axis.text.y = element_text(size=8.75),
			legend.title = element_text(size=10, face="bold", vjust=0.75),
			legend.text = element_text(size=8.75),
			legend.position = "top",
			strip.text = element_text(size=10, face="bold", angle=0),
			panel.background=element_rect(fill="white"),
			panel.grid=element_blank()) +
		labs(fill="sign(coef) * -log(q)", x="", y="") +
		guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth=7.5))

	hh <- 2.5 + ceiling(nlevels(as.factor(df_sig$MicroFeat)) / 7.5)

	wh <- 1.5 * (2.5 + ceiling(nlevels(as.factor(df_sig$Predictor)) * 0.5))

	ggsave(plot=heatmap, file=paste0("./", opts$o, "/heatmap.png"), dpi=300, height=hh, width=wh)

	}

	###############
	# line graphs #
	###############

	# subset the dataframe

	# fixed effects

	df_pc <- df_sig %>%
		filter(grepl("PC", Type)) %>%
		filter(!grepl("\\:", Type)) %>%
		top_n(P_adj, n=-25)

	if (ceiling(nrow(df_pc)/5) > 1){
		
		fe_height <- 3 * ceiling(nrow(df_pc)/5)
	
	} else {
	
		fe_height <- 5
	
	}

	# interaction effects

	df_pcxday <- df_sig %>%
		filter(grepl("PC", Type)) %>%
		filter(grepl("\\:", Type)) %>%
		top_n(P_adj, n=-25)

	if (ceiling(nrow(df_pcxday)/5) > 1){
	
		ie_height <- 3 * ceiling(nrow(df_pcxday)/5)
	
	} else {
	
		ie_height <- 5
	
	}

	##########################
	# plot the fixed effects #
	##########################

	if (nrow(df_pc) > 0){

		combined <- list()

		for (i in df_pc$Model){

			model <- eval(parse(text=i))

			aug <- broom.mixed::augment(model) %>%
				rename(Actual=1) %>%
				rename(Predicted=.fitted) %>%
				rename(PC=2) %>%
				mutate(Model = i) %>%
				mutate(MicroFeat = gsub("lmer\\(| .*", "", Model)) %>%
				mutate(MicroFeat = gsub("_noname_unclassified", "_?_?", MicroFeat)) %>%
				mutate(MicroFeat = gsub("_unclassified", "_?", MicroFeat)) %>%
				mutate(PCN = gsub("*.\\~  | \\*.*", "", Model)) %>%
				mutate(PCN = gsub(".*P", "P", PCN)) %>%
				group_by(PCN) %>%
				mutate(PC_Section = case_when(
					PC <= 0.25 * max(PC) + min(PC) ~ "0.00-0.25", 
					PC  > 0.25 * max(PC) + min(PC) & PC <= 0.50 * max(PC) + min(PC) ~ "0.25-0.50",
					PC  > 0.50 * max(PC) + min(PC) & PC <= 0.75 * max(PC) + min(PC) ~ "0.50-0.75",
					PC  > 0.75 * max(PC) + min(PC) ~ "0.75-1.00"
					)) %>%
				mutate(PC_Quartile = ntile(PC, 4)) %>%
				ungroup() %>%
				group_by(PCN, MicroFeat) %>%
				mutate(Abundance_Quartile = as.factor(ntile(Actual, 4))) %>%
				ungroup() %>%
				select(Subject, MicroFeat, Actual, Abundance_Quartile, Predicted, PC, PCN, PC_Section, PC_Quartile)
	
			combined <- rbind(combined, aug)
	
		}

		if (grepl("EC|Pfam", opts$f) == TRUE){
	
			combined <- combined %>%
				merge(., pathways, by="MicroFeat") %>%
				mutate(FullName = paste0("EC ", PathID, ": ", Name)) %>%
				select(-MicroFeat) %>%
				rename(MicroFeat=PathID)
	
		}

		# map to add q value to plot

		if (grepl("Species|Bray-Curtis", opts$f) == TRUE){
	
			ann_sig <- df_pc %>%
				mutate(P_adj = paste0("q=", scientific(P_adj)))
	
			}

		if (grepl("EC|Pfam", opts$f) == TRUE){
	
			ann_sig <- df_pc %>%
				mutate(P_adj = paste0("q=", scientific(P_adj))) %>%
				rename(FullName=MicroFeat, MicroFeat=PathID)
	
			}
	
		# number of columns in facetted plot

		if (nrow(df_pc) < 5){
		
			ncol = nrow(df_pc)
	
		} else {
	
			ncol = 5	

		}

		######################
		# fixed effects line #
		######################

		combo <- combined %>%
			select(Subject, MicroFeat, PCN, PC, Actual, Predicted) %>%
			melt(id=c("Subject", "PCN", "PC", "MicroFeat"), variable.name="Data", value.name="Value")

		fe_line <- ggplot(combo, aes(x=PC, y=Value)) +
			facet_wrap(PCN ~ MicroFeat, scales="free", ncol=ncol, labeller = label_wrap_gen(width=25)) +
			geom_smooth(aes(colour=Data, fill=Data), method="lm", size=0.5, alpha=0.5) +
			theme_bw() +
			theme(
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=8.75),
				strip.text = element_text(size=8.75, face="bold"),
				legend.title = element_text(size=12.5),
				legend.text = element_text(size=12.5),
				legend.position="top",
				aspect.ratio=0.75
				) +
			labs(title="Fixed effect(s)", y=my_label) +
			geom_text(data=ann_sig, aes(x=Inf, y=Inf, label=P_adj), 
				hjust=1.5,
				vjust=1.5) +
			scale_colour_manual(values=c("Actual"="orangered", "Predicted"="midnightblue")) +
			scale_fill_manual(values=c("Actual"="orangered", "Predicted"="midnightblue"))

		ggsave(plot=fe_line, file=paste0("./", opts$o, "/fixed_effects.png"), height=fe_height, width=15, dpi=300)

		#####################
		# fixed effects box #
		#####################

		# jittered boxplot

		fe_box_1 <- ggplot(combined, aes(x=as.factor(PC_Quartile), y=Actual)) +
			facet_wrap(PCN ~ MicroFeat, scales="free_y", ncol=ncol, labeller = label_wrap_gen(width=25)) +
			geom_jitter(aes(colour=as.factor(PC_Quartile), fill=as.factor(PC_Quartile)), pch=21, alpha=0.125) +
			geom_boxplot(aes(fill=as.factor(PC_Quartile)), colour="darkgrey", outlier.shape=NA, alpha=0.5) +
			theme_bw() +
			theme(
				panel.grid=element_blank(),
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=8.75),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=8.75, face="bold"),
				legend.position = "top",
				aspect.ratio=0.75
				) +
			scale_colour_viridis(discrete=T) +
			scale_fill_viridis(discrete=T) +
			labs(x="PC Quartile", y=my_label, colour="PC Quartile", fill="PC Quartile") +
			geom_text(data=ann_sig, aes(x=Inf, y=Inf, label=P_adj),
				colour="black",
				hjust=1.5,
				vjust=1.5)

		ggsave(plot=fe_box_1, file=paste0("./", opts$o, "/fixed_effects_box_jitter.png"), height=fe_height, width=15, dpi=300)


		# classic boxplot

		fe_box_2 <- ggplot(combined, aes(x=as.factor(PC_Quartile), y=Actual)) +
			facet_wrap(PCN ~ MicroFeat, scales="free_y", ncol=ncol, labeller = label_wrap_gen(width=25)) +
			geom_boxplot(aes(fill=as.factor(PC_Quartile)), colour="darkgrey", outlier.shape=NA) +
			theme_bw() +
			theme(
				panel.grid=element_blank(),
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=8.75),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=8.75, face="bold"),
				legend.position = "top",
				aspect.ratio=0.75
				) +
			scale_colour_viridis(discrete=T) +
			scale_fill_viridis(discrete=T) +
			labs(x="PC Quartile", y=my_label, colour="PC Quartile", fill="PC Quartile") +
			geom_text(data=ann_sig, aes(x=Inf, y=Inf, label=P_adj),
				colour="black",
				hjust=1.5,
				vjust=1.5)

		ggsave(plot=fe_box_2, file=paste0("./", opts$o, "/fixed_effects_box_classic.png"), height=fe_height, width=15, dpi=300)

		#############################
		# fixed effects ribbon plot #
		#############################

		ribbon_in <- combined %>%
			group_by(PCN) %>%
			mutate(PC_BIN = cut(PC, breaks = seq(min(PC), max(PC), by = (max(PC) - min(PC))/100))) %>%
			ungroup() %>%
			group_by(PCN, PC_BIN, MicroFeat) %>%
			mutate(PC_POS = gsub("\\(|\\]", "", PC_BIN)) %>%
			mutate(PC_POS = as.numeric(gsub(".*\\,", "", PC_POS)) - ((as.numeric(gsub(".*\\,", "", PC_POS)) - as.numeric(gsub("\\,.*", "", PC_POS))) / 2)) %>%
			mutate(BIN_MEAN = mean(Actual)) %>%
			mutate(BIN_MIN = min(Actual)) %>%
			mutate(BIN_MAX = max(Actual)) %>%
			ungroup()  %>%
			select(MicroFeat, PCN, PC, PC_BIN, PC_POS, BIN_MEAN, BIN_MIN, BIN_MAX) %>%
			distinct() %>%
			as.data.frame()
	
		fe_ribbon <- ggplot(ribbon_in, aes(x=PC_POS, y=BIN_MEAN)) +
			facet_wrap(PCN ~ MicroFeat, scales="free", ncol=5, labeller = label_wrap_gen(width=25)) +
			geom_ribbon(aes(ymin=BIN_MIN, ymax=BIN_MAX, x=PC, fill = "band"), fill="grey")  +
			geom_line(aes(y=BIN_MEAN), colour="black") +
			theme_bw() +
			theme(
				panel.grid=element_blank(),
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=8.75),
				strip.text = element_text(size=8.75, face="bold"),
				aspect.ratio=0.75
				) +
			geom_text(data=ann_sig, aes(x=Inf, y=Inf, label=P_adj),
				hjust=1.5,
				vjust=1.5) +
			labs(x="PC", y=my_label)

		ggsave(plot=fe_ribbon, file=paste0("./", opts$o, "/fixed_effects_ribbon.png"),  height=fe_height, width=15, dpi=300)

		###################################
		# fixed effects line + error bars #
		###################################

		combined_smry <- combined %>%
			group_by(MicroFeat, PCN, PC_Section) %>%
			summarise_each(funs(mean, sd, se=sd(.)/sqrt(n())), Actual)

		fe_err <- ggplot(combined_smry, aes(x=PC_Section, y=mean, group=1)) +
			facet_wrap(PCN ~ MicroFeat, scales="free_y", ncol=5, labeller = label_wrap_gen(width=25)) +
			geom_line(colour="#899DA4") +
			geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=0.125) +
			theme_bw() +
			theme(
				panel.grid=element_blank(),
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text.x = element_text(angle=45, hjust=1, size=8.75),
				axis.text.y = element_text(size=8.75),
				strip.text = element_text(size=8.75, face="bold"),
				aspect.ratio=0.75
				) +
			labs(x="PC Quarter", y=my_label) +
			geom_label(data=ann_sig, aes(x=Inf, y=Inf, label=P_adj), 
				hjust=2.5,
				vjust=1.5,
				colour="firebrick",
				fill="white")

		ggsave(plot=fe_err, file=paste0("./", opts$o, "/fixed_effects_line_error.png"),  height=fe_height, width=15, dpi=300)

		############################
		# fixed effects ridge plot #
		############################

		ann_sig_2 <- combined %>%
			group_by(PCN) %>%
			summarise( PC = ( max(PC) + min(PC) ) / 2 ) %>%
			ungroup() %>%
			mutate(Abundance_Quartile=1) %>% 
			merge(., ann_sig, by="PCN")

		fe_ridge <- ggplot(combined, aes(x=PC, y=as.factor(Abundance_Quartile))) +
			facet_wrap(PCN ~ MicroFeat, scales="free", ncol=5, labeller = label_wrap_gen(width=25)) +
			geom_density_ridges2(aes(fill=as.factor(Abundance_Quartile)), colour="darkgrey") +
			theme_bw() +
			theme(
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=8.75),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=8.75, face="bold"),
				legend.position = "top",
				aspect.ratio=0.75) +
			scale_colour_viridis(option="magma", discrete=T) +
			scale_fill_viridis(option="magma", discrete=T) +
			labs(x="PC", y="Abundance Quartile", colour="Abundance Quartile", fill="Abundance Quartile") +
			geom_text(data=ann_sig_2, aes(x=PC, y=as.factor(Abundance_Quartile), label=P_adj),
				hjust=0.5,
				vjust=1.5) +
			scale_x_continuous(expand = c(0, 0))
	
		ggsave(plot=fe_ridge, file=paste0("./", opts$o, "/fixed_effects_density_ridges.png"), height=fe_height, width=15, dpi=300)

		#################################
		# fixed effects stacked density #
		#################################

		fe_den_stack_2 <- ggplot(combined, aes(x=PC)) +
			facet_wrap(PCN ~ MicroFeat, scales="free_x", ncol=5, labeller = label_wrap_gen(width=25)) +
			geom_density(aes(fill=as.factor(Abundance_Quartile)), colour="darkgrey", position="fill") +
			theme_bw() +
			theme(
				panel.grid=element_blank(),
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=8.75),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=8.75, face="bold"),
				panel.spacing = unit(2, "lines"),
				legend.position = "top",
				aspect.ratio=0.75) +
			scale_x_continuous(expand = c(0, 0)) +
			scale_y_continuous(expand = expansion(mult = c(0, 0))) +
			scale_colour_viridis(option="magma", discrete=T) +
			scale_fill_viridis(option="magma", discrete=T) +
			labs(colour="Abundance Quartile", fill="Abundance Quartile", x="PC", y="Scaled")

		ggsave(plot=fe_den_stack_2, file=paste0("./", opts$o, "/fixed_effects_density_stack.png"),  height=fe_height, width=15, dpi=300)

		#########################
		# fixed effects density #
		#########################

		ann_sig_d <- combined %>%
			group_by(PCN, MicroFeat) %>%
			summarise( Actual = ( max(Actual) + min(Actual) ) / 2 ) %>%
			ungroup() %>%
			merge(., ann_sig, by=c("PCN", "MicroFeat"))

		fe_den_col <- ggplot(combined, aes(x=Actual)) +
			facet_wrap(PCN ~ MicroFeat, scales="free_y", labeller = label_wrap_gen(width=25), ncol=ncol) +
			geom_density(aes(colour=as.factor(PC_Quartile)), size=0.875) +
			theme_bw() +
			theme(
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=8.75),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text.x = element_text(size=8.75, face="bold"),
				strip.text.y = element_text(angle=0, size=8.75, face="bold"),
				legend.position = "top",
				aspect.ratio=0.75) +
			scale_colour_viridis(discrete=T) +
			scale_fill_viridis(discrete=T) +
			scale_x_continuous(expand = c(0, 0)) +
			scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
			labs(x="log(Abundance [%])", y="Density", colour="PC Quartile", fill="PC Quartile") +
			geom_text(data=ann_sig_d, aes(x=Actual, y=Inf, label=P_adj), 
				alpha=0.75, hjust=0.5, vjust=1.5)

		ggsave(plot=fe_den_col, file=paste0("./", opts$o, "/fixed_effects_density.png"),  height=fe_height, width=15, dpi=300)

		#######################
		# zeros/non-zeros bar #
		#######################

		combined_smry <- combined %>%
			mutate(Zeros = case_when(Actual > min(Actual) ~ "non_zero", Actual == min(Actual) ~ "zero" )) %>%
			group_by(PC_Quartile, MicroFeat, PCN, Zeros) %>%
			tally() %>%
			ungroup() %>%
			dcast(PC_Quartile + MicroFeat + PCN ~ Zeros, value.var="n") %>%
			mutate(perc_non_zero = 100 * non_zero/(non_zero + zero) ) %>%
			arrange(MicroFeat) %>%
			as.data.frame()

		fe_zero <- ggplot(combined_smry, aes(x=as.factor(PC_Quartile), y=perc_non_zero)) +
			facet_wrap(PCN ~ MicroFeat, ncol="5") +
			geom_bar(aes(fill=as.factor(PC_Quartile)), stat="identity", colour="black") +
			theme_bw() +
			theme(
				panel.grid=element_blank(),
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=8.75),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=7.5, face="bold"),
				legend.position = "top",
				aspect.ratio=0.75
				) +
			scale_fill_viridis(discrete=T) +
			labs(x="PC Quartile", y="Percentage of non-zeroes (%)", fill="PC Quartile") +
			scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=c(0, 25, 50, 75, 100))

		ggsave(plot=fe_zero, file=paste0("./", opts$o, "/fixed_effects_non_zeros_bar.png"),  height=fe_height, width=15, dpi=300)

		}

	########################
	# plot the interaction #
	########################

	if (nrow(df_pcxday) > 0) {

		combined_int <- list()

		for (i in df_pcxday$Model) {

			model <- eval(parse(text=i))

			slope <- as.data.frame(coef(model)$Subject) %>%
				tibble::rownames_to_column("Subject") %>%
				mutate(Model = i)

			aug <- broom.mixed::augment(model) %>%
				rename(RelAb=1, PC=2) %>%
				mutate(Model = i) %>%
				mutate(MicroFeat = gsub("lmer\\(| .*", "", Model)) %>%
				mutate(MicroFeat = gsub("_noname_unclassified", "_?_?", MicroFeat)) %>%
				mutate(MicroFeat = gsub("_unclassified", "_?", MicroFeat)) %>%
				mutate(pcn = gsub("*.\\~  | \\*.*", "", Model)) %>%
				mutate(pcn = gsub(".*P", "P", pcn)) %>%
				group_by(pcn) %>%
				mutate(Section = case_when(
					PC <= 0.25 * max(PC) + min(PC) ~ "0.00-0.25", 
					PC  > 0.25 * max(PC) + min(PC) & PC <= 0.50 * max(PC) + min(PC) ~ "0.25-0.50",
					PC  > 0.50 * max(PC) + min(PC) & PC <= 0.75 * max(PC) + min(PC) ~ "0.50-0.75",
					PC  > 0.75 * max(PC) + min(PC) ~ "0.75-1.00"
					)) %>%
				mutate(Quartile = ntile(PC, 4)) %>%
				ungroup() %>%
				select(Subject, MicroFeat, pcn, Section, Quartile, Days, .fitted, RelAb, PC)

			combined_int <- rbind(combined_int, aug)	

			}

		head(combined_int)

		if (grepl("EC|Pfam", opts$f) == TRUE){
	
				combined_int <- combined_int %>%
					merge(., pathways, by="MicroFeat") %>%
					mutate(FullName = paste0("EC ", PathID, ": ", Name)) %>%
					select(-MicroFeat) %>%
					rename(MicroFeat=PathID)
	
			}

		if (grepl("Species|Bray-Curtis", opts$f) == TRUE){
	
			ann_sig_int <- df_pcxday %>%
				mutate(P_adj = paste0("q=", scientific(P_adj)))
	
			}

		if (grepl("EC|Pfam", opts$f) == TRUE){
	
			ann_sig_int <- df_pcxday %>%
				mutate(P_adj = paste0("q=", scientific(P_adj))) %>%
				rename(FullName=MicroFeat, MicroFeat=PathID)
	
			}

		interaction <- ggplot(combined_int, aes(x=Days, y=RelAb)) +
			facet_wrap(pcn ~ MicroFeat, scales="free", ncol=5) +
			geom_smooth(aes(colour=factor(Quartile), fill=factor(Quartile)), method="lm") +
			scale_colour_viridis(discrete=TRUE) +
			scale_fill_viridis(discrete=TRUE) +
			theme_bw() +
			theme(
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=8.75),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=7.5, face="bold"),
				legend.position = "top",
				aspect.ratio=0.75) +
			labs(title="Interaction effect(s)", y=my_label,
				colour="Quartile",
				fill="Quartile") +
			geom_text(data=ann_sig_int, aes(x=Inf, y=Inf, label=P_adj), 
				hjust=1.5,
				vjust=1.5)

		ggsave(plot=interaction, file=paste0("./", opts$o, "/interaction.png"), height=ie_height, width=15, dpi=300)

		########################
		# slope per individual #
		########################

		quartile_map <- combined_int %>% 
			select(Subject, Quartile) %>%
			distinct()

		library(purrr)
		library(broom)

		slope_df <- combined_int %>%
		  split(., list(.$Subject, .$MicroFeat, .$pcn), drop = TRUE) %>%
		  map(~lm(RelAb ~ Days, data = .x)) %>% 
		  map_df(tidy, .id = "Subject.MicroFeat.pcn") %>%
		  filter(term == "Days") %>%
		  mutate(PC = gsub(".*\\.", "", Subject.MicroFeat.pcn)) %>%
		  mutate(Subject.MicroFeat = gsub(".PC.*", "", Subject.MicroFeat.pcn)) %>%
		  mutate(Subject = gsub("\\..*", "", Subject.MicroFeat)) %>%
		  as.data.frame() %>%
		  drop_na()


		if (grepl("Species|Bray-Curtis", opts$f) == TRUE){
	
			slope_df <- slope_df %>%
				mutate(MicroFeat = gsub(".*\\.", "", Subject.MicroFeat))
	
			}
 
		if (grepl("EC|Pfam", opts$f) == TRUE){
	
			slope_df <- slope_df %>%
				mutate(MicroFeat = str_sub(Subject.MicroFeat, 8) )
	
			}

		DF <- merge(slope_df, quartile_map, by="Subject") %>%
			mutate(Quartile = factor(Quartile))

		ann_sig_sb <- ann_sig_int %>%
			rename(PC=PCN)

		slope_box <- ggplot(data=DF, aes(x=Quartile, y=estimate)) +
			facet_wrap(PC~MicroFeat, scales="free", ncol=5) +
			geom_jitter(aes(colour=Quartile, fill=Quartile), pch=21, alpha=0.125) +
			geom_boxplot(aes(fill=Quartile), colour="darkgrey", outlier.shape=NA) +
			theme_bw() +
			theme(
				panel.grid=element_blank(),
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=8.75),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=7.5, face="bold"),
				legend.position = "top",
				aspect.ratio=0.75
				) +
			scale_colour_viridis(discrete=T) +
			scale_fill_viridis(discrete=T) +
			labs(title="Slope per subject", x="Quartile", y="Slope", colour="Quartile", fill="Quartile") +
			scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
			geom_text(data=ann_sig_sb, aes(x=Inf, y=Inf, label=P_adj), 
				hjust=1.5,
				vjust=1.5)

		ggsave(plot=slope_box, file=paste0("./", opts$o, "/slope_per_subject_box.png"), height=ie_height, width=15, dpi=300)

		# slope ridge plot

		ann_sig_sb <- ann_sig_int %>%
			rename(PC=PCN)

		ann_sig_2b <- slope_df %>%
			group_by(PC, MicroFeat) %>%
			summarise( estimate = ( max(estimate) + min(estimate) ) / 2 ) %>%
			ungroup() %>%
			mutate(Quartile=1) %>% 
			merge(., ann_sig_sb, by=c("PC", "MicroFeat"))
	
		slope_ridge <- ggplot(DF, aes(x=estimate, y=as.factor(Quartile))) +
			facet_wrap(PC ~ MicroFeat, scales="free", ncol=5, labeller = label_wrap_gen(width=25)) +
			geom_density_ridges2(aes(fill=as.factor(Quartile)), colour="darkgrey") +
			theme_bw() +
			theme(
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=8.75),
				legend.title = element_text(size=12.5, face="bold", vjust=0.75),
				legend.text = element_text(size=12.5),
				strip.text = element_text(size=7.5, face="bold"),
				legend.position = "top",
				aspect.ratio=0.75) +
			scale_colour_viridis(discrete=T) +
			scale_fill_viridis(discrete=T) +
			geom_vline(xintercept=0, linetype="dashed", colour="darkgrey") +
			labs(x="Slope", y="PC Quartile", fill="PC Quartile") +
			geom_text(data=ann_sig_2b, aes(x=estimate, y=as.factor(Quartile), label=P_adj),
				hjust=0.5,
				vjust=1.25)

		ggsave(plot=slope_ridge, file=paste0("./", opts$o, "/slope_per_subject_ridges.png"), height=ie_height, width=15, dpi=300)

		############################
		# combined interaction fig #
		############################

		int_combo <- cowplot::plot_grid(interaction, slope_ridge, slope_box, axis="tblr", align="hv", nrow=3, scale=0.95)

		ggsave(plot=int_combo, file=paste0("./", opts$o, "/interaction_figure.png"), height=3*ie_height, width=15, dpi=300)

		######################
		# random slopes line #
		######################

		slope_sub <- slope_df %>%
			select(Subject, MicroFeat, PC, estimate) %>%
			rename(pcn=PC)

		rand_subs <- combined_int %>%
			group_by(Quartile, Subject) %>%
			mutate(max_day = max(Days)) %>%
			ungroup() %>%
			group_by(Quartile) %>%
			sample_n(10) %>%
			ungroup() %>%
			droplevels() %>%
			select(Subject) %>%
			merge(., combined_int, by="Subject") %>%
			merge(., slope_sub, by=c("Subject", "MicroFeat", "pcn")) %>%
			arrange(Subject) %>%
			group_by(Quartile) %>%
			mutate(rank=min_rank(Subject)) %>%
			ungroup() %>%
			arrange(Quartile, rank) %>%
			as.data.frame()

		ind_ids <- rand_subs %>% 
			select(Quartile, Subject, rank) %>% 
			distinct() %>% 
			group_by(Quartile) %>%
			mutate(ind_id = min_rank(rank)) %>%
			filter(ind_id < 11) %>%
			select(-rank) %>%
			data.frame()

		RandSubs <- merge(rand_subs, ind_ids, by=c("Subject", "Quartile")) %>%
			mutate(ind_id = factor(ind_id)) %>%
			mutate(Quartile = paste0("Q", Quartile))

		rand_subs_fit <- RandSubs %>%
			select(-RelAb) %>%
			rename(RelAb = .fitted)

		rand_subs_p <- ggplot(RandSubs, aes(x=Days, y=RelAb)) +
			facet_grid(MicroFeat ~ Quartile, scales="free_x") +
			geom_line(data=rand_subs_fit, aes(x=Days, y=RelAb, colour=ind_id)) +
			theme_bw() +
			theme(
				panel.grid = element_blank(),
				plot.title = element_text(size=12.5, face="bold", hjust=0.5),
				axis.title = element_text(size=12.5, face="bold"),
				axis.text = element_text(size=8.75),
				strip.text.x = element_text(size=8.75, face="bold"),
				strip.text.y = element_text(size=8.75, face="bold", angle=0),
				legend.position="none") +
			labs(title="Random slopes (10 individuals from each PC Quartile)") +
			scale_colour_brewer(palette="Set3")

		ggsave(plot=rand_subs_p, file=paste0("./", opts$o, "/slope_per_subject_line.png"), height=10, width=10, dpi=300)

		}
	}

####################
# plot time effect #
####################

if (!grepl("EC|Pfam", opts$f) == TRUE){
	
	df_time <- model_stats %>%
		filter(!grepl("PC", Predictor)) %>%
		mutate(PC = gsub("*.\\~  | \\*.*", "", Model)) %>%
		mutate(PC = gsub(".*P", "P", PC)) %>%
		filter(!grepl("PCX", PC)) %>%
		group_by(PC) %>%
		mutate(P_adj = p.adjust(P, method="bonferroni")) %>%
		ungroup() %>%
		filter(P_adj < 0.05) %>%
		mutate(Species = gsub("lmer\\(", "", Model) ) %>%
		mutate(Species = gsub(" .*", "", Species) ) %>%
		mutate(Species = gsub("_noname_unclassified", "_?_?", Species)) %>%
		mutate(Species = gsub("_unclassified", "_?", Species)) %>%
		mutate(Coefficient = case_when((Coefficient > 0)~1, (Coefficient < 0)~-1)) %>%
		mutate(value = Coefficient * -log(P_adj)) %>%
		filter(!grepl("Country", Predictor)) %>%
		arrange(desc(value)) %>%
		as.data.frame()

	width_days <- nlevels(as.factor(df_time$Species)) / 5

	hmap_time <- ggplot(data=df_time, aes(x=reorder(Species, value), y=PC, fill=value)) +
		geom_tile(colour="black") +
		scale_fill_gradient2(low = "steelblue", 
			mid = "white", 
			high = "darksalmon", 
			na.value = NA, breaks=c(round(0.75 * min(df_time$value)), 0, round(0.75 * max(df_time$value)))) +
		scale_x_discrete(expand=c(0, 0)) +
		scale_y_discrete(expand=c(0, 0)) +
		theme_bw() +
		theme(
			plot.title = element_text(size=12.5, face="bold", hjust=0.5),
			plot.subtitle = element_text(size=12.5, face="bold", hjust=0.5),
			axis.title = element_text(size=12.5, face="bold"),
			axis.text.x = element_text(size=7.5, angle=45, hjust=1),
			axis.text.y = element_text(size=7.5),
			legend.title = element_text(size=12.5, face="bold", vjust=0.75),
			legend.text = element_text(size=8.75),
			legend.position = "top",
			strip.text.y = element_text(size=12.5, face="bold", angle=0),
			panel.background=element_rect(fill="white"),
			panel.grid=element_blank(),
			plot.margin=unit(c(0.5, 0.5, 0.5, 1), "cm")) +
		ggtitle("Age effect") +
		labs(x="Species", y="PC")

	ggsave(plot=hmap_time, file=paste0("./", opts$o, "/heatmap_day_effect.png"),
		height=6.25, width=width_days, dpi=300)

	}

#
