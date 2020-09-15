#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plot_results_pcs.R [-s <stats> -d <data> -o <out_dir> -p <padj_method>]

Options:
   -s stats (output from run_models.R)
   -d data
   -o output directory for plots [default: out_dir]
   -p method for p-value adjustment [default: bonferroni]

' -> doc

opts <- docopt(doc)

#################
# load the data #
#################

library(tidyverse)
library(reshape2)

# data (input for models)

data <- read.csv(opts$d, sep="\t", header=T) %>%
	setNames(gsub("microbiome\\.|genetics\\.", "", names(.))) %>%
	setNames(gsub("\\.", "_", names(.)))

colnames(data) <- gsub("microbiome\\.|genetics\\.", "", colnames(data))

############################
# significant associations #
############################

library(dplyr)
library(reshape2)
library(ggplot2)
library(viridis)
library(scales)
library(lmerTest)

# input

model_stats <- read.csv(opts$s, sep="\t", header=T)

# output directory

dir.create(opts$o)

####################
# p-value adjustment

df_sig <- model_stats %>%
	group_by(Type) %>%
	mutate(P_adj = p.adjust(P, method=opts$p)) %>%
	mutate(P_fdr = p.adjust(P, method="fdr")) %>%
	mutate(P_bon = p.adjust(P, method="bonferroni")) %>%
	ungroup() %>%
	filter(P_adj < 0.05) %>%
	filter(grepl("PC", Predictor)) %>%
	mutate(Coefficient = case_when((Coefficient > 0) ~ 1, (Coefficient < 0) ~ -1)) %>%
	mutate(value = Coefficient * -log(P_adj)) %>%
	mutate(MicroFeat = gsub("lmer\\(", "", Model)) %>%
	mutate(MicroFeat = gsub(" .*", "", MicroFeat)) %>%
	mutate(MicroFeat = gsub("_noname_unclassified", "_?_?", MicroFeat)) %>%
	mutate(MicroFeat = gsub("_unclassified", "_?", MicroFeat))  %>%
	mutate(pcn = gsub("\\:.*", "", Predictor))

###########
# add names

# EC names

if(grepl("EC_", df_sig$MicroFeat[1]) == TRUE){

pathways <- read.csv("humann2_names.txt", sep="\t", header=T)

df_sig <- merge(df_sig, pathways, by="MicroFeat")

}

# Pfam names

if(grepl("Pfam_", df_sig$MicroFeat[1]) == TRUE){

pathways <- read.csv("map_pfam_names.tsv", sep="\t", header=T)

df_sig <- merge(df_sig, pathways, by="MicroFeat")

}

#

write.table(df_sig, file=paste0("./", opts$o, "/model_stats.tsv"), sep="\t", quote=F, row.names=F)

########################################
# visualization of significant results #
########################################

#########
# heatmap

if(!grepl("EC_|Pfam_", df_sig$MicroFeat[1]) == TRUE){

heatmap <- ggplot(data=df_sig, aes(x=MicroFeat, y=Predictor, fill=value)) +
	geom_tile(colour="black") +
	facet_grid(Type ~ ., space="free", scales="free") +
	scale_fill_gradient2(low = "steelblue", 
		mid = "white", 
		high = "darksalmon") +
	scale_x_discrete(expand=c(0, 0)) +
	scale_y_discrete(expand=c(0, 0)) +
	theme_bw() +
	theme(
		axis.title = element_text(size=12.5, face="bold"),
		axis.text.x = element_text(size=10, angle=45, hjust=1),
		axis.text.y = element_text(size=10),
		legend.title = element_text(size=12.5, face="bold", vjust=0.75),
		legend.text = element_text(size=8.75),
		legend.position = "top",
		strip.text.y = element_text(size=12.5, face="bold", angle=0),
		panel.background=element_rect(fill="white"),
		panel.grid=element_blank(),
		aspect.ratio=1) +
	labs(fill="sign(coef) * -log(q)", x="", y="")

} else {

heatmap <- ggplot(data=df_sig, aes(x=Name, y=Predictor, fill=value)) +
	geom_tile(colour="black") +
	facet_grid(Type ~ ., space="free", scales="free") +
	scale_fill_gradient2(low = "steelblue", 
		mid = "white", 
		high = "darksalmon") +
	scale_x_discrete(expand=c(0, 0)) +
	scale_y_discrete(expand=c(0, 0)) +
	theme_bw() +
	theme(
		axis.title = element_text(size=12.5, face="bold"),
		axis.text.x = element_text(size=10, angle=45, hjust=1),
		axis.text.y = element_text(size=10),
		legend.title = element_text(size=12.5, face="bold", vjust=0.75),
		legend.text = element_text(size=8.75),
		legend.position = "top",
		strip.text.y = element_text(size=12.5, face="bold", angle=0),
		panel.background=element_rect(fill="white"),
		panel.grid=element_blank(),
		aspect.ratio=1) +
	labs(fill="sign(coef) * -log(q)", x="", y="")

}

ggsave(plot=heatmap, file=paste0("./", opts$o, "/heatmap.png"), dpi=300, width=10)

##############################################
# line graphs for each significant association

# subset the dataframe

df_pc <- df_sig %>%
	filter(grepl("PC", Type)) %>%
	filter(!grepl("\\:", Type)) %>%
	droplevels()

df_pcxday <- df_sig %>%
	filter(grepl("PC", Type)) %>%
	filter(grepl("\\:", Type))

########################
# plot the fixed effects

combined <- list()

for (i in df_pc$Model) {

model <- eval(parse(text=i))

aug <- broom.mixed::augment(model) %>%
	rename(PC=2) %>%
	mutate(Model = i) %>%
	mutate(MicroFeat = gsub("lmer\\(| .*", "", Model)) %>%
	mutate(MicroFeat = gsub("_noname_unclassified", "_?_?", MicroFeat)) %>%
	mutate(MicroFeat = gsub("_unclassified", "_?", MicroFeat)) %>%
	mutate(PCN = gsub("*.\\~  | \\*.*", "", Model)) %>%
	mutate(PCN = gsub(".*P", "P", PCN)) %>%
	group_by(PCN) %>%
	mutate(Section = case_when(
		PC <= 0.25 * max(PC) + min(PC) ~ "0.00-0.25", 
		PC  > 0.25 * max(PC) + min(PC) & PC <= 0.50 * max(PC) + min(PC) ~ "0.25-0.50",
		PC  > 0.50 * max(PC) + min(PC) & PC <= 0.75 * max(PC) + min(PC) ~ "0.50-0.75",
		PC  > 0.75 * max(PC) + min(PC) ~ "0.75-1.00"
		)) %>%
	ungroup() %>%
	rename(Abundance=1) %>%
	select(Subject, MicroFeat, Abundance, PCN, .fitted, PC, Section)

combined <- rbind(combined, aug)
	
}

# map to add q value to plot

ann_sig <- df_sig %>%
	 filter(!grepl(":", Predictor)) %>%
	 mutate(PCN = Predictor) %>%
	 mutate(P_adj = paste0("q=", scientific(P_adj))) %>%
	 select(MicroFeat, PCN, P_adj, Coefficient)

# fixed effects line

fe_line <- ggplot(combined, aes(x=PC, y=.fitted)) +
	facet_wrap(~ MicroFeat + PCN, scales="free", ncol=5) +
#	geom_point(colour="grey", alpha=0.25) +
	geom_smooth(method="lm", colour="black", fill="firebrick", size=0.5) +
	theme_bw() +
	theme(
		plot.title = element_text(size=12.5, face="bold", hjust=0.5),
		axis.title = element_text(size=12.5, face="bold"),
		axis.text = element_text(size=8.75),
		legend.title = element_text(size=12.5, face="bold", vjust=0.75),
		legend.text = element_text(size=12.5),
		strip.text = element_text(size=8.75, face="bold"),
		legend.position = c(7/8, 1/6),
		aspect.ratio=0.75) +
	labs(title="Fixed effect(s)") +
	geom_text(data=ann_sig, aes(x=Inf, y=Inf, label=P_adj), 
		hjust=1.5,
		vjust=1.5)

ggsave(plot=fe_line, file=paste0("./", opts$o, "/fixed_effects.png"), width=15, dpi=300)

# fixed effects box

fe_box <- ggplot(combined, aes(x=Section, y=.fitted)) +
	facet_wrap(~ MicroFeat + PCN, scales="free_y", ncol=5) +
	geom_jitter(aes(colour=Section, fill=Section), pch=21, alpha=0.125) +
	geom_boxplot(aes(colour=Section, fill=Section), colour="black", outlier.shape=NA) +
	theme_bw() +
	theme(
		panel.grid=element_blank(),
		plot.title = element_text(size=12.5, face="bold", hjust=0.5),
		axis.title = element_text(size=12.5, face="bold"),
		axis.text.x = element_text(angle=45, hjust=1, size=8.75),
		axis.text.y = element_text(size=8.75),
		legend.title = element_text(size=12.5, face="bold", vjust=0.75),
		legend.text = element_text(size=12.5),
		strip.text = element_text(size=8.75, face="bold"),
		legend.position = "none",
		aspect.ratio=0.75
		) +
	scale_colour_viridis(discrete=T) +
	scale_fill_viridis(discrete=T) +
	labs(x="Relative PC", y=".fitted", colour="Relative PC", fill="Relative PC")

ggsave(plot=fe_box, file=paste0("./", opts$o, "/fixed_effects_box.png"), width=15, dpi=300)

# line error bars

combined_smry <- combined %>%
	group_by(MicroFeat, PCN, Section) %>%
	summarise_each(funs(mean, sd, se=sd(.)/sqrt(n())), .fitted)

fe_err <- ggplot(combined_smry, aes(x=Section, y=mean, group=1)) +
	facet_wrap(MicroFeat~PCN, scales="free_y", ncol=5) +
	geom_line(colour="#C93312") +
	geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.125) +
	geom_point(size=1.75, pch=21, colour="black", fill="#899DA4") +
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
	labs(x="Relative PC", y=".fitted")

ggsave(plot=fe_err, file=paste0("./", opts$o, "/fixed_effects_line_error.png"), width=15, dpi=300)

######################
# plot the interaction

combined <- list()

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
	mutate(Quartile = ntile(PC, 4)) %>%
	select(Subject, MicroFeat, pcn, Quartile, Day, .fitted, RelAb, PC)

combined <- rbind(combined, aug)	

}

ann_sig <- df_sig %>%
	filter(grepl(":", Predictor)) %>%
	mutate(pcn = Predictor) %>%
	mutate(pcn = gsub("\\:.*", "", pcn)) %>%
	mutate(P_adj = paste0("q=", scientific(P_adj))) %>%
	select(MicroFeat, pcn, P_adj, Coefficient)

interaction <- ggplot(combined, aes(x=Day, y=.fitted)) +
	facet_wrap(~ MicroFeat + pcn, scales="free", ncol=5) +
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
		strip.text = element_text(size=8.75, face="bold"),
		legend.position = "top",
		aspect.ratio=0.75) +
	labs(title="Interaction effect(s)", 
		colour="PC Quartile",
		fill="PC Quartile") +
	geom_text(data=ann_sig, aes(x=Inf, y=Inf, label=P_adj), 
		hjust=1.5,
		vjust=1.5)

ggsave(plot=interaction, file=paste0("./", opts$o, "/interaction.png"), width=15, dpi=300)

######################
# slope per individual

quartile_map <- combined %>% 
	select(Subject, Quartile) %>%
	distinct()

library(purrr)
library(broom)

slope_df <- combined %>%
  split(., list(.$Subject, .$MicroFeat, .$pcn), drop = TRUE) %>%
  map(~lm(.fitted ~ Day, data = .x)) %>% 
  map_df(tidy, .id = "Subject.MicroFeat.pcn") %>%
  filter(term == "Day") %>%
  mutate(PC = gsub(".*\\.", "", Subject.MicroFeat.pcn)) %>%
  mutate(Subject.MicroFeat = gsub(".PC.", "", Subject.MicroFeat.pcn)) %>%
  mutate(Subject = gsub("\\..*", "", Subject.MicroFeat)) %>%
  mutate(MicroFeat = gsub(".*\\.", "", Subject.MicroFeat)) %>%
  as.data.frame()

DF <- merge(slope_df, quartile_map, by="Subject") %>%
	mutate(Quartile = factor(Quartile))

slope_box <- ggplot(data=DF, aes(x=Quartile, y=estimate)) +
	facet_wrap(MicroFeat~PC, scales="free", ncol=5) +
	geom_jitter(aes(colour=Quartile, fill=Quartile), pch=21, alpha=0.125) +
	geom_boxplot(aes(fill=Quartile), colour="black", outlier.shape=NA) +
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
	labs(title="Slope per subject",
		y="Slope")

ggsave(plot=slope_box, file=paste0("./", opts$o, "/slope_per_subject_box.png"), width=15, dpi=300)

# combined interaction fig

int_combo <- plot_grid(interaction, slope_box, scale=0.95, align="v", nrow=2)
int_combo

ggsave(plot=int_combo, file=paste0("./", opts$o, "/interaction_figure.png"), width=15, dpi=300)

####################
# random slopes line

slope_sub <- slope_df %>%
	select(Subject, MicroFeat, PC, estimate) %>%
	rename(pcn=PC)

rand_subs <- combined %>%
	group_by(Quartile, Subject) %>%
	mutate(max_day = max(Day)) %>%
	ungroup() %>%
	#filter(max_day > 500 & max_day < 730) %>%
	group_by(Quartile) %>%
	sample_n(10) %>%
	ungroup() %>%
	droplevels() %>%
	select(Subject) %>%
	merge(., combined, by="Subject") %>%
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

rand_subs_p <- ggplot(RandSubs, aes(x=Day, y=RelAb)) +
	facet_grid(MicroFeat ~ Quartile, scales="free_x") +
	geom_line(data=rand_subs_fit, aes(x=Day, y=RelAb, colour=ind_id)) +
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
	scale_colour_viridis(option="magma", discrete=T)

ggsave(plot=rand_subs_p, file=paste0("./", opts$o, "/slope_per_subject_line.png"), height=10, width=10, dpi=300)

##################
# plot time effect

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

hmap_time <- ggplot(data=df_time, aes(x=reorder(Species, value), y=PC, fill=value)) +
	geom_tile(colour="black") +
	scale_fill_gradient2(low = "steelblue", 
		mid = "white", 
		high = "darksalmon", 
		na.value = NA) +
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
		plot.margin=unit(c(0.5, 0.5, 0.5, 1),"cm"),
		aspect.ratio=0.1) +
	ggtitle("Species ~ PC * Day + ( 1 + Day | Subject )") +
	labs(subtitle="Day effect", fill="sign(coef) * -log(q)", x="Species")

ggsave(plot=hmap_time, file=paste0("./", opts$o, "/heatmap_day_effect.png"), height=7.5, width=12.5, dpi=300)