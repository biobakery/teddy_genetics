#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plot_path_bars.R [--inp <input> --met <metadata> --names <pathway_names> --min <min_age> --max <max_age>]

Options:
   --inp input
   --met metadata
   --names pathway names
   --min minimum age in subset
   --max maximum age in subset
   
' -> doc

opts <- docopt(doc)

#####

library(tidyverse)
library(reshape2)

##########
# metadata

metadata_for_samples <- read.csv(opts$met, header=T, sep="\t")

head(metadata_for_samples)

min_age <- as.numeric(opts$min)
max_age <- as.numeric(opts$max)

df_samples <- metadata_for_samples %>%
	filter(Day > min_age & Day <= max_age) %>%
	select(sample_id) %>%
	distinct()

head(df_samples)

#####
# ECs

ec_map <- read.csv(opts$names, header=F, sep="\t") %>%
	rename(EC=1, EC_Name=2) %>%
	mutate(FullName = paste0(EC, ": ", EC_Name))

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

#######
# plots

dir.create(paste0("./days_", min_age, "-", max_age, "_stacked_bar_charts"), showWarnings = FALSE)

ec_list <- levels(factor(strat$EC))

plot_list = list()

for (i in ec_list) {

DF <- strat %>%
	filter(grepl(i, EC)) %>%
	select(EC, FullName, Species) %>%
	distinct() %>%
	merge(., strat, by=c("EC", "FullName", "Species")) %>%
	data.frame()

df_dcast <- DF %>% 
	dcast(Species ~ sample_id, value.var="relab") %>%
	replace(is.na(.), 0) %>%
	tibble::column_to_rownames("Species")

hc <- hclust(dist(df_dcast))
dd <- as.dendrogram(hc)
order_bug <- rev(c(labels(dd)))

hc2 <- hclust(dist(t(df_dcast)))
dd2 <- as.dendrogram(hc2)
order_samples <- rev(c(labels(dd2)))

DF$Species = factor(DF$Species, levels=order_bug)
DF$sample_id = factor(DF$sample_id, levels=order_samples)

p <- ggplot(DF, aes(x=factor(sample_id), y=relab, fill=Species)) +
	facet_grid(~FullName, space="free", scales="free") +
	geom_bar(stat="identity", position="stack", width=1) +
	scale_fill_brewer(palette="Set3") +
	theme_bw() +
	theme(panel.grid=element_blank(),
		plot.title=element_text(face="bold", hjust=0.5),
		axis.text.x=element_blank(),
		axis.title.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.title.y=element_text(face="bold"),
		strip.text=element_text(face="bold"),
		legend.title=element_text(face="bold")) +
	scale_y_continuous(expand = c(0,0)) +
	labs(y="Relative abundance")

plot_list[[i]] = p

ggsave(p, 
	file=paste0("./days_", min_age, "-", max_age, "_stacked_bar_charts/plot_", i,".png"), 
	height=5, width=10, dpi=300)
	
}

pdf("EC_stacked_bar_charts.pdf", height=5, width=10)
for (i in 1:length(ec_list)) {
    print(plot_list[[i]])
}
dev.off()

###