#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plot_path_bars.R [--inp <input> --met <metadata> --names <pathway_names> --min <min_age> --max <max_age> --feat <feature>]

Options:
   --inp input
   --met metadata
   --names pathway names
   --min minimum age in subset
   --max maximum age in subset
   --feat feature (EC, Pfam)
   
' -> doc

opts <- docopt(doc)

#####

library(tidyverse)
library(reshape2)
library(cowplot)

############
# metadata #
############

metadata_for_samples <- read.csv(opts$met, header=T, sep="\t")

########
# ages #
########

min_age <- as.numeric(opts$min)
max_age <- as.numeric(opts$max)

##################
# filter samples #
##################

df_samples <- metadata_for_samples %>%
	filter(subject_age_days > min_age & subject_age_days <= max_age) %>%
	select(sample_id) %>%
	distinct()

#################
# pathway names #
#################

ec_map <- read.csv(opts$names, header=F, sep="\t") %>%
	rename(EC=1, EC_Name=2) %>%
	mutate(FullName = paste0(EC, ": ", EC_Name))

######################
# pathway abundances #
######################

strat <- read.csv(opts$inp, header=T, sep="\t", check.names=F) %>%
	rename(EC=1) %>%
	melt(variable.name="sample_id", value.name="cpm") %>%
	mutate(sample_id = gsub("_.*", "", sample_id)) %>%
	merge(., df_samples, by="sample_id") %>%
	filter(cpm > 0) %>%
	mutate(Species = gsub(".*\\.s__", "s__", EC)) %>%
	mutate(Species = gsub(".*\\|unclassified", "unclassified", Species)) %>%
	mutate(EC = gsub("\\|.*", "", EC)) %>%
	group_by(EC, sample_id) %>%
	mutate(relab = cpm / sum(cpm)) %>%
	ungroup() %>%
	select(sample_id, EC, Species, relab) %>%
	merge(., ec_map, by="EC") %>%
	as.data.frame()

#############################
# list species with pathway #
#############################

species_list <- strat %>% 
	select(Species) %>%
	distinct()
	
write.table(species_list,
	paste0("days_", min_age, "-", max_age, ".", opts$feat, ".species_bar.tsv"), 
	sep="\t", row.names=F, quote=F)

###################################
# assign a colour to each species #
###################################

my_pal <- c(
	's__Bacteroides_fragilis'='#E6194B',
	's__Bacteroides_ovatus'='#3CB44B',
	's__Bacteroides_uniformis'='#FFE119',
	's__Bacteroides_vulgatus'='#4363D8',
	's__Bifidobacterium_bifidum'='#F58231',
	's__Bifidobacterium_breve'='#911EB4',
	's__Bifidobacterium_longum'='#42D4F4',
	's__Bifidobacterium_pseudocatenulatum'='#F032E6',
	's__Clostridium_nexile'='#BFEF45',
	's__Eubacterium_rectale'='#FABED4',
	's__Faecalibacterium_prausnitzii'='#469990',
	's__Ruminococcus_bromii'='#DCBEFF',
	's__Ruminococcus_gnavus'='#9A6324',
	's__Ruminococcus_sp_5_1_39BFAA'='#FFFAC8',
	's__Ruminococcus_torques'='#800000',
	's__Collinsella_aerofaciens'='#AAFFC3',
	's__Escherichia_coli'='#808000',
	's__Clostridium_ramosum'='#FFD8B1',
	's__Veillonella_parvula'='#000075',
	's__Akkermansia_muciniphila'='#A9A9A9'
	)

library(tidyverse)

head(legend_in) %>%
	mutate(Total = paste0("n=", n()))

##########################
# run loop to make plots #
##########################

dir.create(paste0("./days_", min_age, "-", max_age, "_", opts$feat, "_stacked_bar_charts"), showWarnings = FALSE)

ec_list <- levels(factor(strat$EC))

plot_list = list()

for (i in ec_list) {

DF <- strat %>%
	filter(grepl(i, EC)) %>%
	select(EC, FullName, Species) %>%
	distinct() %>%
	merge(., strat, by=c("EC", "FullName", "Species")) %>%
	filter(relab > 0) %>%
	droplevels() %>%
	data.frame()

df_dcast <- DF %>% 
	dcast(Species ~ sample_id, value.var="relab") %>%
	replace(is.na(.), 0) %>%
	tibble::column_to_rownames("Species")

if (nlevels(factor(DF$Species)) > 1){
	hc <- hclust(dist(df_dcast))
	dd <- as.dendrogram(hc)
	order_bug <- rev(c(labels(dd)))
	DF$Species = factor(DF$Species, levels=order_bug)
	
	hc2 <- hclust(dist(t(df_dcast)))
	dd2 <- as.dendrogram(hc2)
	order_samples <- rev(c(labels(dd2)))
	DF$sample_id = factor(DF$sample_id, levels=order_samples)

}

p <- ggplot(DF, aes(x=factor(sample_id), y=relab, colour=Species, fill=Species)) +
	facet_grid(~FullName, space="free", scales="free") +
	geom_bar(stat="identity", position="stack", width=1) +
	scale_colour_manual(values=my_pal) +
	scale_fill_manual(values=my_pal) +
	theme_bw() +
	theme(
		axis.text.x=element_blank(),
		axis.title.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.text.y=element_text(size=10),
		axis.title.y=element_text(size=10, face="bold"),
		strip.text=element_text(size=10, face="bold"),
		legend.position="none") +
	scale_x_discrete(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	labs(y="Relative abundance")

plot_list[[i]] = p

ggsave(plot=p, 
	file=paste0("./days_", min_age, "-", max_age, "_", opts$feat, "_stacked_bar_charts/plot_", i,".png"), 
	height=5, width=10, dpi=300)
	
}

###################################
# make the legend from dummy data #
###################################

legend_in <- strat %>% 
	select(Species) %>%
	mutate(Group = "Group") %>%
	mutate(Value = 1) %>%
	distinct()

legend_in <- data.frame(Group="Group", Value=1, Species=species_list$Species)

legend <- ggplot(legend_in, aes(x=Group, y=Value)) +
		geom_bar(aes(fill=Species), colour="black", size=0.5, stat="identity", position="stack", width=1) +
		theme(
			legend.title=element_text(size=10, face="bold"),
			legend.text=element_text(size=10)) +
		scale_fill_manual(values=my_pal)

###################################
# set the direction of the legend #
###################################

if ((length(ec_list)/5) %% 1 == 0) {
	if (nlevels(as.factor(species_list$Species)) < 5){
			legend <- legend + guides(fill = guide_legend(direction="horizontal", nrow=1)) 
		} else {
			legend <- legend + 	guides(fill = guide_legend(direction="horizontal", ncol=5))
	}
} else {
	if (nlevels(as.factor(legend_in$Species)) < 10){
		legend <- legend + guides(fill = guide_legend(ncol=1))
	} else {
		legend <- legend + guides(fill = guide_legend(ncol=2, hjust=0.5))
	}
}

################################
# merge the plots with cowplot #
################################

LEGEND <- cowplot::get_legend(legend)

if ((length(ec_list)/5) %% 1 == 0) {

	p_combo <- cowplot::plot_grid(plotlist=plot_list, align="hv", axis="tblr", scale=0.9, ncol=5, labels="AUTO")

	p_combo <- cowplot::plot_grid(p_combo, LEGEND, nrow=2, rel_heights=c(1, 0.125))

} else {
	
	p_combo <- cowplot::plot_grid(plotlist=c(plot_list, LEGEND), align="hv", axis="tblr", scale=0.9, ncol=5, labels="AUTO")
	
}

height <- 3 * ceiling(length(ec_list) / 5)

ggsave(plot=p_combo, 
	file=paste0("./days_", min_age, "-", max_age, "_", opts$feat, "_stacked_bar_charts/cowplot.png"), 
	height=height, width=30, dpi=300)

###
