#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plot_path_smry.R [--input <input> --feat <feature> --names <pathway_names> --strat <path_strat> --padj <padj> --group <group_by> --mpa <mpa_species> --out <prefix_out>]

Options:
   --input input
   --feat feature (EC or Pfam)
   --names pathway names
   --strat list of stratified pathways
   --padj p-value adjustment method
   --group for p-value adjustment, group by Type or Predictor [default: Type] 
   --mpa metaphlan2 major
   --out output prefix
   
' -> doc

opts <- docopt(doc)

###

library(tidyverse)
library(reshape2)

model_anova_total <- read.csv(opts$input, sep="\t", header=T)

paths_map <- read.csv(opts$names, header=F, sep="\t") %>%
	rename(Path=1, Path_Name=2) %>%
	mutate(FullName = paste0(Path, ": ", Path_Name)) %>%
	mutate(Path = paste0(opts$feat, "_", Path))

df_a <- model_anova_total %>%
	filter(!grepl(":", Predictor) & grepl("PC", Predictor)) %>%
	group_by(!!sym(opts$g)) %>%
	mutate(P_adj = p.adjust(P, method=opts$padj)) %>%
	ungroup() %>%
	filter(P_adj < 0.05) %>%
	mutate(Coefficient = case_when((Coefficient > 0)~1, (Coefficient < 0)~-1)) %>%
	mutate(value = Coefficient * -log10(P_adj)) %>%
	mutate(Path = gsub("lmer\\(", "", Model) ) %>%
	mutate(Path = gsub(" .*", "", Path) ) %>%
	merge(., paths_map, by="Path") %>%
	as.data.frame()

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

bar_height <- nrow(paths_num_bug)/4.5

ggsave(file=paste0(opts$out, "_contributions_bar.png"), dpi=300, width=10, height=bar_height)

###

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

write.table(paths_rare, paste0(opts$out, "_humann2_stratified_names_sig.tsv"), 
	quote=F,
	col.names=F,
	row.names=F)

# cluster the dataframe

df <- paths_sig %>%
	select(1, 4) %>%
	merge(., metaphlan2_major, by="variable") %>%
	mutate(Path = gsub("EC_|Pfam_", "", Path)) %>%
	rename(Species=variable) %>%
	mutate(Presence=c(1)) %>%
	select(Path, Species, Presence) %>%
	group_by(Path) %>%
	mutate(n = sum(Presence)) %>%
	ungroup()

df_dcast <- df %>%
	dcast(Path~Species, value.var="Presence") %>%
	replace(is.na(.), 0) %>%
	column_to_rownames("Path")

# order_pathways
hc <- hclust(dist(df_dcast))
dd <- as.dendrogram(hc)
order_paths <- rev(c(labels(dd)))

# order species
hc2 <- hclust(dist(t(df_dcast)))
dd2 <- as.dendrogram(hc2)
order_bug <- rev(c(labels(dd2)))

df$Path = factor(df$Path, levels=order_paths)
df$Species = factor(df$Species, levels=order_bug)

hmap_cont <- ggplot(data=df, aes(x=reorder(Path, -n), y=Species, fill=Presence)) +
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
		plot.margin=unit(c(0.5, 0.5, 0.5, 1),"cm")) +
	labs(x=opts$feat, y="Species") +
	coord_fixed(ratio=0.5)

h_height <- nlevels(factor(df$Path))/5
h_width <- nlevels(factor(df$Species))/5

ggsave(file=paste0(opts$out, "_contributions_heatmap.png"), dpi=300, height=h_height, width=h_width)

###