#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plot_path_smry.R [--input <input> --feat <feature> --names <pathway_names> --strat <path_strat> --padj <padj> --group <group_by> --mpa <mpa_species> --out <out_dir>]

Options:
   --input input
   --feat feature (EC or Pfam)
   --names pathway names
   --strat list of stratified pathways
   --padj p-value adjustment method
   --group for p-value adjustment, group by Type or Predictor [default: Type] 
   --mpa metaphlan2 major
   --out output directory
   
' -> doc

opts <- docopt(doc)

###

library(tidyverse)
library(reshape2)

dir.create(opts$o, showWarnings=FALSE)

model_stats <- read.csv(opts$input, sep="\t", header=T)

paths_map <- read.csv(opts$names, header=F, sep="\t") %>%
	rename(Path=1, Path_Name=2) %>%
	mutate(FullName = paste0(Path, ": ", Path_Name)) %>%
	mutate(Path = paste0(opts$feat, "_", Path)) %>%
	mutate(Path = gsub("\\.", "_", Path))

df_a <- model_stats %>%
	filter(Predictor != "Country" & Predictor != "Days") %>%
	mutate(Predictor = gsub("genetics_", "", Predictor)) %>%
	mutate(Type = case_when(grepl(":", Predictor) ~ "Interaction", !grepl(":", Predictor) ~ "Fixed")) %>%
	filter(grepl("PC", Predictor)) %>%
	mutate(Type = gsub(".*\\:", "PC\\:", Type)) %>%
	mutate(Type = case_when( grepl(":", Type) ~ "PC:Days" , !grepl(":", Type) ~ "PC" )) %>%
	group_by(!!sym(opts$g)) %>%
	mutate(P_adj = p.adjust(P, method=opts$p)) %>%
	ungroup() %>%
	filter(P_adj < 0.05) %>%
	mutate(Path = gsub("lmer\\(", "", Model)) %>%
	mutate(Path = gsub(" .*", "", Path)) %>%
	droplevels()

head(df_a)

if (nrow(df_a) == 0){
	
	print("No significant features.")
	
	write.table(df_a, paste0(gsub("\\..*", "", opts$i), ".humann2_stratified_names.sig_", opts$feat, ".tsv"), 
		quote=F,
		col.names=F,
		row.names=F)
	
	}

if (nrow(df_a) > 0) {
	
	paths_map_strat <- read.csv(opts$strat, header=T, sep="\t") %>%
		rename(Stratified=1) %>%
		mutate(Path = gsub("\\|.*", "", Stratified)) %>%
		mutate(Path = paste0(opts$feat, "_", Path)) %>%
		mutate(Path = gsub("\\.", "_", Path)) %>%
		mutate(Species = gsub(".*\\|", "", Stratified))

	paths_sig <- df_a %>% 
		select(Path) %>% 
		distinct() %>%
		merge(., paths_map_strat, by="Path") %>%
		mutate(variable = gsub(".*\\.s__", "s__", Species))

	metaphlan2_major <- read.csv(opts$mpa, sep="\t", header=T) %>%
		rename(variable=1)

	paths_bugs_num <- paths_sig %>%
		merge(., metaphlan2_major, by="variable") %>%
		group_by(Path) %>%
		tally() %>%
		arrange(n) %>%
		merge(., paths_map, by="Path") %>%
		as.data.frame()

	paths_bar <- ggplot(paths_bugs_num, aes(x=reorder(FullName, n), y=n, label=n)) +
		geom_bar(stat="identity", colour="black", fill="firebrick") +
		theme_classic() +
		ggtitle("Number of species with EC") +
		labs(x="", y="Number of species") +
		theme(plot.title = element_text(size=15, face="bold", hjust=0.5),
			axis.title = element_text(size=12.5, face="bold"),
			axis.text.x = element_text(size=12.5),
			axis.text.y = element_text(size=8.75)) +
		scale_y_continuous(expand=expansion(mult=c(0, 0.075))) +
		geom_text(hjust=-0.125) +
		coord_flip()

	x <- nrow(paths_bugs_num)/4.5
	
	if (x < 5) {
		bar_height = 5		
	} else {
		bar_height = x
	}

	ggsave(plot=paths_bar, file=paste0("./", opts$out, "/contributions_bar.png"), dpi=300, width=10, height=5)

	paths_bugs_list <- paths_sig %>%
		merge(., metaphlan2_major, by="variable")
	
	paths_rare <- paths_sig %>%
		merge(., metaphlan2_major, by="variable") %>%
		group_by(Path) %>%
		tally() %>%
		filter(n <= 12) %>%
		as.data.frame() %>%
		merge(., paths_bugs_list, by="Path") %>%
		select(Stratified)

	write.table(paths_rare, paste0(gsub("\\..*", "", opts$i), ".humann2_stratified_names.sig_", opts$feat, ".tsv"), 
		quote=F,
		col.names=F,
		row.names=F)

	# cluster the dataframe

	df <- paths_sig %>%
		select(1, 4) %>%
		merge(., metaphlan2_major, by="variable") %>%
		merge(., paths_map, by="Path") %>%
		rename(Species=variable) %>%
		mutate(Presence=c(1)) %>%
		select(FullName, Species, Presence) %>%
		group_by(FullName) %>%
		mutate(n = sum(Presence)) %>%
		ungroup()

	df_dcast <- df %>%
		dcast(FullName~Species, value.var="Presence") %>%
		replace(is.na(.), 0) %>%
		column_to_rownames("FullName")

	if (nlevels(factor(df$FullName)) > 1){
		# order_pathways
		hc <- hclust(dist(df_dcast))
		dd <- as.dendrogram(hc)
		order_paths <- rev(c(labels(dd)))
		df$FullName = factor(df$FullName, levels=order_paths)
		}

	if (nlevels(factor(df$Species)) > 1){
		# order species
		hc2 <- hclust(dist(t(df_dcast)))
		dd2 <- as.dendrogram(hc2)
		order_bug <- rev(c(labels(dd2)))
		df$Species = factor(df$Species, levels=order_bug)
		}

	hmap_cont <- ggplot(data=df, aes(y=reorder(FullName, n), x=Species, fill=Presence)) +
		geom_tile(colour="black") +
		scale_fill_gradientn(colours=c("firebrick")) +
		scale_x_discrete(expand=c(0, 0)) +
		scale_y_discrete(expand=c(0, 0), position="right") +
		theme_bw() +
		theme(
			plot.title = element_text(size=12.5, face="bold", hjust=0.5),
			axis.title = element_text(size=12.5, face="bold"),
			axis.text.x = element_text(size=10, angle=45, hjust=1),
			axis.text.y = element_text(size=8.75),
			legend.position = "none",
			strip.text = element_text(size=12.5, face="bold", angle=0),
			panel.background=element_rect(fill="white"),
			panel.grid=element_blank(),
			plot.margin=unit(c(0.5, 0.5, 0.5, 2.5),"cm")) +
		labs(x=opts$feat, y="Species")
	
	y <- nlevels(factor(df$FullName))/2.5
	
	if (y < 5){
		h_height = 5
	} else {
		h_height = y	
	}
	
	h_width <- nlevels(factor(df$Species))/2
	h_width

	ggsave(file=paste0("./", opts$out, "/contributions_heatmap.png"), dpi=300, height=h_height, width=h_width)

	}

#
