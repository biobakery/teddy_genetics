############################################################
## Extended Figure 4
##
## Significant microbial species associated with microbiome
## maturation trajectories across age bins
############################################################

###########################
## 0. Load libraries
###########################
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringi)
library(stringr)
library(cowplot)

###########################
## 1. Load trajectory clustering results
###########################
setwd("~/ddong/TEDDY_project/")
load("./output/traj_cluster_results.Rdata")

###########################
## 2. Define age bins
###########################
bin <- cut(
  metadata_cluster$Days,
  breaks = seq(100, 800, 100),
  include.lowest = TRUE
)

metadata_cluster$bins <- bin
binnames <- names(table(bin))

###########################
## 3. Load MaAsLin species association results
###########################
maaslin_root <- "./output/maaslin_species/"
setwd(maaslin_root)
load("./all_species.Rdata")

###########################
## 4. Select significant species associations
###########################
all_species.df <- all_species.df %>%
  filter(Type != "Cluster1 vs Cluster2")

sig <- all_species.df %>%
  filter(qval < 0.25)

alldata_select_species <- all_species.df %>%
  filter(feature %in% sig$feature)

alldata_select_species <- alldata_select_species %>%
  mutate(
    sig = case_when(
      qval < 0.001 ~ "***",
      qval < 0.05 ~ "**",
      qval < 0.25 ~ "*"
    )
  )

alldata_select_species$bin <- factor(
  alldata_select_species$bin,
  levels = binnames
)

alldata_select_species <- alldata_select_species %>%
  mutate(
    Species = gsub("s__", "", feature)
  )

###########################
## 5. Generate heatmap
###########################
pdf(
  "./figures/Extended/Extended_Fig4_species.pdf",
  width = 23,
  height = 20
)

ggplot(
  alldata_select_species,
  aes(
    x = bin,
    y = Species,
    fill = coef
  )
) +
  geom_tile() +
  geom_text(
    aes(label = sig),
    size = 10
  ) +
  facet_wrap(~Type) +
  scale_fill_gradient2(
    low = "darkgreen",
    mid = "white",
    high = "#7209b7"
  ) +
  theme_cowplot() +
  theme(
    axis.text.y = element_text(size = 20)
  )

dev.off()