############################################################
## Figure 2c (combined)
## (1) Relative abundance of selected species across 100-day bins
##     stratified by microbiome maturation clusters
## (2) MaAsLin effect sizes (coefficients) across bins and contrasts
############################################################
## Set project root
setwd("~/ddong/TEDDY_project/")

###########################
## 0. Load packages and set project root
###########################
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(cowplot)

###########################
## 1. Load shared inputs
###########################
## Trajectory clustering results (provides Cluster_results and metadata_cluster)
load("./output/traj_cluster_results.Rdata")

###########################
## 2. Panel c: Relative abundance across time bins by cluster
###########################
## Input: post-QC MetaPhlAn species table + sample metadata (mgx_age)
microbiome <- fread("./data/processed/metaphlan2_afterQC.tsv", header = TRUE)
load("./data/metadata/metadata_all.Rdata")

## Reshape to sample-by-species table (percent scale)
microbiome <- microbiome %>% column_to_rownames("ID")
microbiome <- microbiome %>%
  t(.) %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  rename(sample_id = 1)

## Convert to percentage scale
microbiome[, -1] <- microbiome[, -1] * 100

map.df <- metadata[, c("subject_id", "sample_id", "mgx_age")]
map.df$sample_id <- as.character(map.df$sample_id)

microbiome <- microbiome %>% inner_join(map.df)

## Standardize variable names and types
microbiome$sample_id  <- as.character(microbiome$sample_id)
microbiome$subject_id <- as.character(microbiome$subject_id)

## Use mgx_age as the analysis time variable (Days)
microbiome$Days <- as.character(microbiome$mgx_age)
microbiome$mgx_age <- NULL

microbiome.df <- melt(microbiome)
microbiome.df$Days <- as.numeric(microbiome.df$Days)

## Restrict to < 800 days
microbiome.df <- microbiome.df %>% filter(Days < 800)


## Add cluster labels (subject-level) to each sample
Cluster.df <- Cluster_results %>% select(subject_id, Cluster)
Cluster.df$subject_id <- as.character(Cluster.df$subject_id)
microbiome.df <- microbiome.df %>% inner_join(Cluster.df)

## Bin time into 100-day windows (100–800 days)
microbiome.df$bins <- cut(
  microbiome.df$Days,
  breaks = c(seq(100, 800, 100)),
  include.lowest = TRUE
)
microbiome.df <- microbiome.df %>% filter(!is.na(bins))

## Summarize abundance within each bin and cluster
microbiome_summary.df <- microbiome.df %>%
  group_by(bins, variable, Cluster) %>%
  summarise(
    mean   = mean(value),
    median = median(value),
    .groups = "drop"
  )

## Select highlighted species and define a color palette
select_species <- c(
  "Bifidobacterium_bifidum",
  "Ruminococcus_gnavus",
  "Clostridium_hathewayi",
  "Dorea_longicatena"
)
select_species <- gsub("_", " ", select_species)

color <- c("#023047", "#a53860", "#ee964b", "#5e548e")
names(color) <- select_species

microbiome_summary.df$species <- gsub("s__", "", microbiome_summary.df$variable)
microbiome_summary.df$species <- gsub("_", " ", microbiome_summary.df$species)

microbiome_summary_selected.df <- microbiome_summary.df %>%
  filter(species %in% select_species)

## Panel c plot object
p2c_part1 <- ggplot(microbiome_summary_selected.df,
                    aes(x = bins, y = mean, fill = species)) +
  geom_bar(stat = "identity") +
  geom_smooth(
    data = microbiome_summary_selected.df,
    aes(x = as.numeric(bins), y = mean),
    method = "loess",
    se = FALSE
  ) +
  facet_grid(species ~ Cluster, scales = "free") +
  scale_fill_manual(values = color) +
  xlab("Time (Days)") +
  ylab("Mean relative abundance (%)") +
  theme_cowplot() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  geom_point()

###########################
## 3. Panel d: MaAsLin coefficients across bins and contrasts
###########################
## Define bins from metadata_cluster (consistent with MaAsLin binning)
metadata_cluster$bins <- cut(
  metadata_cluster$Days,
  breaks = seq(100, 800, 100),
  include.lowest = TRUE
)
binnames <- names(table(metadata_cluster$bins))

## Root directory containing per-bin MaAsLin outputs
maaslin_root <- "./output/maaslin_species/"

## Load and stack MaAsLin results across bins
all_species.df <- data.frame()

for (i in seq_along(binnames)) {
  b_n <- binnames[i]
  fpath <- file.path(
    maaslin_root,
    paste("Cluster_together", b_n, sep = "_"),
    "all_results.tsv"
  )
  
  tmp <- fread(fpath)
  
  ## Harmonize contrast labels used downstream
  if (i <= 3) {
    tmp <- tmp %>%
      filter(metadata == "Cluster_together") %>%
      mutate(
        Type = paste(value, "vs Cluster1"),
        bin  = b_n
      )
  } else {
    tmp <- tmp %>%
      filter(metadata == "Cluster_together") %>%
      mutate(
        Type = paste(value, "vs Cluster1&2"),
        bin  = b_n
      )
  }
  
  all_species.df <- rbind(all_species.df, tmp)
}

## Select significant species (q < 0.25) and create significance labels
sig <- all_species.df %>% filter(qval < 0.25)

alldata_select_species <- all_species.df %>%
  filter(feature %in% sig$feature) %>%
  mutate(
    sig = case_when(
      qval < 0.001 ~ "***",
      qval < 0.05  ~ "**",
      qval < 0.25  ~ "*",
      TRUE         ~ ""
    ),
    bin     = factor(bin, levels = binnames),
    Species = gsub("^s__", "", feature)
  )

## Additional filtering: strong effects (|coef| > 0.5) among significant features
sig_strong_effect <- all_species.df %>%
  filter(qval < 0.25) %>%
  filter(abs(coef) > 0.5)

## Contrast-specific subsets
sig_cluster23_vs_1 <- sig_strong_effect %>%
  filter(Type == "Cluster2_3 vs Cluster1")

sig_cluster3_vs_12 <- sig_strong_effect %>%
  filter(Type == "Cluster3 vs Cluster1&2")

## Species repeatedly detected across bins (overall)
recurrent_species_all <- table(sig_strong_effect$feature)
recurrent_species_all <- names(recurrent_species_all[recurrent_species_all > 2])

## Species repeatedly detected in Cluster3 vs Cluster1&2 contrast
recurrent_species_cluster3 <- table(sig_cluster3_vs_12$feature)
recurrent_species_cluster3 <- names(recurrent_species_cluster3[recurrent_species_cluster3 > 1])

## Final species set used for plotting
selected_species_recurrent <- unique(c(
  recurrent_species_all,
  recurrent_species_cluster3,
  "s__Bifidobacterium_longum"   # manually retained key species
))

## Clean species names (remove MetaPhlAn prefix)
selected_species_recurrent <- gsub("^s__", "", selected_species_recurrent)
alldata_select_species2 <- alldata_select_species %>%
  filter(Species %in% selected_species_recurrent)

## Order species and create within-panel ranks
alldata_select_species2$Type2 <- factor(
  alldata_select_species2$Type,
  levels = c("Cluster3 vs Cluster1&2", "Cluster2_3 vs Cluster1")
)

x <- alldata_select_species2 %>% arrange(Type2, desc(bin), -coef)
alldata_select_species2$Species <- factor(
  alldata_select_species2$Species,
  levels = sort(unique(x$Species))
)

alldata_select_species2 <- alldata_select_species2 %>%
  group_by(bin, Type) %>%
  mutate(rank = rank(-coef)) %>%
  ungroup()

## Color palette for species bars (ensure sufficient colors for plotted species)
col_part2 <- c(
  "#2c6e49", "#023047", "#219ebc",
  "#ee964b", "#f4d35e", "#b6ad90",
  "#5e548e", "#a53860", "#da627d"
)

## Panel d plot object
p2c_part2 <- ggplot(alldata_select_species2,
                    aes(x = coef, y = as.factor(rank), fill = Species)) +
  geom_bar(stat = "identity") +
  facet_wrap(Type ~ bin, nrow = 2, scales = "free_y") +
  geom_text(aes(label = sig), size = 5, color = "darkred") +
  scale_fill_manual(values = col_part2) +
  theme_bw() +
  theme_cowplot() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.6)
  ) +
  xlab("MaAsLin coefficient") +
  ylab("Species (ranked)")

###########################
## 4. Combine panels into a single Figure 2c
###########################
p2c_combined <- plot_grid(
  p2c_part1,
  p2c_part2,
  ncol = 1,
  rel_heights = c(1.8, 1),
  labels = c("c", "d"),
  label_size = 14
)

###########################
## 5. Save combined figure
###########################

pdf("./figures/main/Fig2c_combined.pdf", width = 14, height = 18)
print(p2c_combined)
dev.off()
