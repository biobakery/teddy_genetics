############################################################
## Supplementary Figure 2
## Association between host genetics and gut microbiome
## structure across developmental time points.
############################################################

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)

###########################
## 1. Set input/output paths
###########################
set.seed(1)
setwd("~/ddong/TEDDY_project/")

metadata_file  <- "./output/preprocess/metadata_microbiome.Rdata"
microbiome_file <- "./output/microbiome_PC.Rdata"
genetics_file  <- "./data/genetics.bed/filter/teddy.qc.PCA_PC1-PC5.tsv"
output_dir     <- "./figures/supplement/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

###########################
## 2. Load and merge inputs
###########################

load(metadata_file)     # -> metadata
load(microbiome_file)   # -> microbiomePCo (sample_id, Microbiome_PCo1)

genetic_PC <- fread(genetics_file) %>%
  rename(subject_id = 1) %>%
  select(subject_id, PC1)

meta_sel <- metadata %>%
  select(subject_id, sample_id, mgx_age) %>%
  mutate(sample_id = as.numeric(sample_id))

microbiomePCo <- microbiomePCo %>%
  mutate(sample_id = as.numeric(sample_id))

merged <- meta_sel %>%
  inner_join(microbiomePCo, by = "sample_id") %>%
  inner_join(genetic_PC,   by = "subject_id") %>%
  filter(!is.na(mgx_age), !is.na(Microbiome_PCo1), !is.na(PC1))

###########################
## 3. Assign developmental time windows
###########################

breaks <- c(0, 365, 730, 1095, 1460,1825,2190)
labels <- c("0-1y", "1-2y", "2-3y",
            "3-4y", "4-5y", "5-6y")

merged <- merged %>%
  mutate(window = cut(mgx_age, breaks = breaks,
                      labels = labels, include.lowest = TRUE)) %>%
  filter(!is.na(window))

###########################
## 4. One sample per participant per window
##    (earliest sample within each window)
###########################

merged <- merged %>%
  group_by(window, subject_id) %>%
  slice_min(mgx_age, n = 1, with_ties = FALSE) %>%
  ungroup()

###########################
## 5. Generate faceted figure
###########################

n_label <- merged %>%
  group_by(window) %>%
  summarise(
    n = sum(complete.cases(PC1, Microbiome_PCo1)),
    .groups = "drop"
  )

p <- ggplot(merged, aes(x = PC1, y = Microbiome_PCo1)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_smooth(method = "lm", alpha = 0.15, linewidth = 0.7) +
  stat_cor(
    method = "spearman",
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 3
  ) +
  geom_text(
    data = n_label,
    aes(
      x = -Inf,
      y = Inf,
      label = paste0("N = ", n)
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 3.2,
    size = 3
  ) +
  facet_wrap(~ window, scales = "free") +
  labs(x = "Genetic PC1", y = "Microbiome PCo1") +
  theme_cowplot() +
  theme(strip.background = element_rect(fill = "grey92"))

###########################
## 6. Per-window correlation table
###########################

cor.df <- merged %>%
  group_by(window) %>%
  summarise(
    n = n(),
    r = cor(PC1, Microbiome_PCo1, method = "pearson"),
    p = cor.test(PC1, Microbiome_PCo1)$p.value,
    .groups = "drop"
  )

###########################
## 7. Save outputs
###########################
ggsave("./figures/supplement/SupFig2.pdf", p, width = 10, height = 6.5)
