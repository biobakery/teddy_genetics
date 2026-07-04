############################################################
## SupFigure 9: Genetics main and interaction effects (heatmaps)
############################################################

## Set project root
setwd("~/ddong/TEDDY_project/")

###########################
## 0. Load packages
###########################
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggside)
library(cowplot)
library(grid)   


###########################
## 1. Load data objects from .Rdata
###########################
load("./data/SupFig12/infant_species_all_nobrst.Rdata")
load("./data/SupFig12/infant_species_all_brst.Rdata")
load("./data/SupFig12/toddler_species_all.Rdata")

###########################
## 2. Harmonize labels 
###########################
df_sig_infant_brst$Class   <- "Infant Breastfeeding"
df_sig_infant_nobrst$Class <- "Infant Post-breastfeeding"

###########################
## 3. Combine infant + toddler results (column-consistent)
###########################
shared_cols <- intersect(colnames(df_sig_Toddler), colnames(df_sig_infant_brst))
results_all <- rbind(
  df_sig_infant_brst[, shared_cols],
  df_sig_infant_nobrst[, shared_cols],
  df_sig_Toddler[, shared_cols]
)

## Significance stars (Bonferroni)
results_all <- results_all %>%
  mutate(
    star = case_when(
      P_adj_bon < 0.001 ~ "***",
      P_adj_bon < 0.01  ~ "**",
      P_adj_bon < 0.05  ~ "*",
      TRUE              ~ ""
    )
  )

## Extract species name from the model string (same logic as your original)
results_all <- results_all %>%
  mutate(
    species = str_extract(Model, "(?<=lmer\\().*(?=\\ \\~)"),
    species = gsub("`", "", species)
  )

###########################
## 4. Load taxonomy/phylogeny annotation and merge
###########################
tax_tbl <- fread("./output/species_major_tree.txt")
results_all <- results_all %>%
  left_join(tax_tbl, by = c("species" = "species"))

###########################
## 5. Keep significant combinations and format species labels
###########################
results_all <- results_all %>%
  mutate(group_key = paste(Predictor, Class))

sig_combos <- results_all %>%
  filter(P_adj_bon < 0.05)

results_sig <- results_all %>%
  group_by(Type) %>%
  filter(
    group_key %in% sig_combos$group_key,
    MicroFeat %in% sig_combos$MicroFeat
  ) %>%
  ungroup()

## Clean species names for display
results_sig <- results_sig %>%
  mutate(
    species = gsub("^s__", "", species),
    species = gsub("_", " ", species)
  ) %>%
  arrange(phylum)

## Class ordering
results_sig$classtype <- factor(
  results_sig$Class,
  levels = c("Infant Breastfeeding", "Infant Post-breastfeeding", "toddler")
)

## Species ordering (y-axis)
results_sig$species <- factor(results_sig$species, levels = unique(results_sig$species))

###########################
## 6. Plotting function (reduces duplicated code)
###########################
make_heatmap_plot <- function(df, panel_title) {
  ggplot(df, aes(y = species, x = Predictor, fill = Coe)) +
    geom_tile(colour = "black") +
    facet_grid(~classtype, space = "free", scales = "free") +
    scale_fill_gradient2(
      low = "#023047",
      mid = "white",
      high = "#9e2a2b"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_line(colour = "black", size = 0.4),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_line(colour = "black", size = 0.4),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8, colour = "black"),
      axis.ticks.length = unit(0.04, "inch"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),
      legend.position = "bottom",
      legend.text = element_text(colour = "black", size = 8),
      legend.title = element_text(colour = "black", size = 8),
      plot.title = element_text(size = 11, face = "bold", hjust = 0)
    ) +
    guides(
      fill = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = 1.5,
        barwidth = 10
      )
    ) +
    geom_text(aes(label = star), color = "black", size = 3, nudge_y = -0.2) +
    geom_ysidetile(aes(x = "phylum", yfill = phylum)) +
    labs(fill = "Coefficient", x = "", y = "", title = panel_title)
}

###########################
## 7. Build panel plots
###########################
df_main_effect <- results_sig %>% filter(Type == "genetics")
df_interaction <- results_sig %>% filter(Type == "genetics:time")

p_main_effect        <- make_heatmap_plot(df_main_effect,  "Main genetic effects")
p_interaction_effect <- make_heatmap_plot(df_interaction, "Genetic × time interaction effects")

###########################
## 8. Combine and save SupFigure9
###########################
SupFig12 <- plot_grid(
  p_main_effect,
  p_interaction_effect,
  ncol = 2,
  rel_heights = c(1, 1),
  labels = c("a","b")
)


pdf("./figures/supplement/SupFigure9.pdf", width = 9, height = 7)
print(SupFig12)
dev.off()
