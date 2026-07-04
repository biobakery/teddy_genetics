############################################################
## Extended Data  Figures1: Top 10 species and EC feature distribution heatmaps
## a) Top 10  species relative abundance across samples
## 5) Top 10 EC relative abundance across matched samples
############################################################

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggside)
library(viridis)
library(cowplot)

###########################
## 1. Set project root and load input data
###########################
setwd("~/ddong/TEDDY_project/")
load("./output/Permanova.Rdata")
load("./output/microbiome_PC.Rdata")

genetic_pc <- fread("./data/genetics.bed/filter/teddy.qc.PCA_PC1-PC5.tsv") %>%
  rename(subject_id = IID)

load("./output/preprocess/metadata_all.Rdata")

sex_metadata <- metadata %>%
  select(subject_id, Sex, time_to_brstfed_stop) %>%
  unique()

sample_metadata_select <- metadata %>%
  select(subject_id, sample_id, Days)

mds.data$sample_id <- as.numeric(mds.data$sample_id)

sample_metadata <- sample_metadata_select %>%
  inner_join(mds.data) %>%
  inner_join(genetic_pc) %>%
  inner_join(metadata_subject)

sample_metadata$Country <- as.factor(sample_metadata$Country)

###########################
## 2. Rank species by mean relative abundance
###########################
species_rank_df <- data.frame(
  aveRA = colMeans(df_sp)
)

species_rank_df$species <- rownames(species_rank_df)

species_rank_df <- species_rank_df %>%
  arrange(desc(aveRA))

###########################
## 3. Prepare top 10 species heatmap data
###########################
top10_species <- df_sp %>%
  select(one_of(species_rank_df$species[1:10]))

top10_species$sample_id <- rownames(top10_species)

top10_species_long <- melt(top10_species) %>%
  rename(
    species = variable,
    RA = value
  )

top10_species_long$sample_id <- as.numeric(top10_species_long$sample_id)

sample_info <- sample_metadata %>%
  select(Days, Country, sample_id)

top10_species_long$RA[top10_species_long$RA == 0] <- 6e-07

top10_species_long <- top10_species_long %>%
  inner_join(sample_info) %>%
  left_join(b1) %>%
  arrange(Country, Days)

top10_species_long$species <- factor(
  top10_species_long$species,
  levels = rev(species_rank_df$species[1:10])
)

top10_species_long$sample_id <- factor(
  top10_species_long$sample_id,
  levels = unique(top10_species_long$sample_id)
)

###########################
## 4. Plot top 10 species heatmap
###########################
p_species <- ggplot(top10_species_long, aes(sample_id, species)) +
  geom_tile(aes(fill = log10(RA * 100))) +
  geom_xsidetile(aes(y = "Country", xfill = Country)) +
  scale_xfill_manual(
    values = c("#1e6091", "#e63946", "#fcbf49", "#6d597a")
  ) +
  scale_fill_viridis_c(
    name = "Relative abundance\nlog10 scale",
    na.value = "white",
    option = "rocket"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

###########################
## 5. Load and rank EC features by mean relative abundance
###########################
ec_abundance <- fread("~/ddong/TEDDY/output_2024/ec_all_qc_1e6.txt")
ec_annotation <- fread("~/ddong/TEDDY/output_2024/ec_name.txt")

ec_rank_df <- data.frame(
  aveRA = colMeans(ec_abundance[, -1])
)

ec_rank_df$EC_feature <- rownames(ec_rank_df)

ec_annotation <- ec_annotation %>%
  mutate(EC_feature = paste("X", EC, sep = ""))

ec_rank_df <- ec_rank_df %>%
  left_join(ec_annotation) %>%
  arrange(desc(aveRA)) %>%
  filter(!EC_feature %in% c("UNGROUPED", "UNMAPPED"))

###########################
## 6. Prepare top 10 EC heatmap data
###########################
top10_ec <- ec_abundance %>%
  select(sample_id, one_of(ec_rank_df$EC_feature[1:10]))

top10_ec$sample_id <- as.character(top10_ec$sample_id)

top10_ec_long <- melt(top10_ec) %>%
  rename(
    EC_feature = variable,
    ECRA = value
  )

top10_ec_long$sample_id <- as.numeric(top10_ec_long$sample_id)

top10_ec_long <- top10_ec_long %>%
  inner_join(sample_info) %>%
  filter(sample_id %in% top10_species_long$sample_id) %>%
  inner_join(ec_annotation)

top10_ec_long$sample_id <- factor(
  top10_ec_long$sample_id,
  levels = unique(top10_species_long$sample_id)
)

top10_ec_long$EC_name <- factor(
  top10_ec_long$EC_name,
  levels = unique(top10_ec_long$EC_name)
)

ec_rank_df <- ec_rank_df %>%
  mutate(ID = paste(EC, EC_name, sep = "_"))

top10_ec_long <- top10_ec_long %>%
  mutate(ID = paste(EC, EC_name, sep = "_"))

top10_ec_long$ID <- factor(
  top10_ec_long$ID,
  levels = rev(ec_rank_df$ID[1:10])
)

###########################
## 7. Plot top 10 EC heatmap
###########################
p_ec <- ggplot(top10_ec_long, aes(sample_id, ID)) +
  geom_tile(aes(fill = log10(ECRA * 100))) +
  geom_xsidetile(aes(y = "Days", xfill = Days)) +
  scale_fill_viridis_c(
    name = "Relative abundance\nlog10 scale",
    na.value = "white",
    option = "mako"
  ) +
  scale_xfill_gradient(
    low = "#dedbd2",
    high = "#4a5759"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

###########################
## 8. Save combined heatmap
###########################
pdf("./figures/Extended/Extended_Fig1_top10_distribution.pdf", 14, 6)

plot_grid(
  p_species,
  p_ec,
  ncol = 1,
  align = 1
)

dev.off()