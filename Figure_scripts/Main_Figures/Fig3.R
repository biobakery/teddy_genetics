############################################################
## Figure 3
## Functional groups and EC-level patterns across time and clusters
##
## Panels:
##  (a-1) Functional groups total abundance trajectories across days (LOESS)
##  (a-2) Species composition within Functional groups (stacked fraction bars)
##  (b) Heatmap of scaled EC abundance across time bins and maturation clusters, stratified by functional group
############################################################

## Set project root
setwd("~/ddong/TEDDY_project/")

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(scales)
library(ggnewscale)
library(ggside)
library(viridis)

###########################
## 1. Load trajectory / metadata objects
###########################
load("./output/traj_cluster_results.Rdata")  # expects metadata_cluster at least

metadata_cluster <- metadata_cluster %>%
  mutate(
    sample_id = as.character(sample_id),
    Cluster   = as.character(Cluster),
    Days      = as.numeric(Days)
  )

###########################
## 2. Load KO annotation and define functional groups
###########################
KO <- fread("./data/external/KO_03202022.csv") %>%
  mutate(
    EC_number = gsub("^EC", "", ec),
    feature   = paste0("X", EC_number)
  )

select_pathway <- c(
  "00260 Glycine, serine and threonine metabolism [PATH:ko00260]",
  "00270 Cysteine and methionine metabolism [PATH:ko00270]",
  "00280 Valine, leucine and isoleucine degradation [PATH:ko00280]",
  "00290 Valine, leucine and isoleucine biosynthesis [PATH:ko00290]",
  "00300 Lysine biosynthesis [PATH:ko00300]",
  "00310 Lysine degradation [PATH:ko00310]",
  "00340 Histidine metabolism [PATH:ko00340]",
  "00360 Phenylalanine metabolism [PATH:ko00360]",
  "00380 Tryptophan metabolism [PATH:ko00380]",
  "00400 Phenylalanine, tyrosine and tryptophan biosynthesis [PATH:ko00400]",
  "00670 One carbon pool by folate [PATH:ko00670]",
  "00730 Thiamine metabolism [PATH:ko00730]",
  "00740 Riboflavin metabolism [PATH:ko00740]",
  "00750 Vitamin B6 metabolism [PATH:ko00750]",
  "00760 Nicotinate and nicotinamide metabolism [PATH:ko00760]",
  "00770 Pantothenate and CoA biosynthesis [PATH:ko00770]",
  "00780 Biotin metabolism [PATH:ko00780]",
  "00790 Folate biosynthesis [PATH:ko00790]",
  "00052 Galactose metabolism [PATH:ko00052]"
)

KO_select <- KO %>%
  filter(cat3 %in% select_pathway) %>%
  mutate(
    EC_top = case_when(
      cat3 %in% "00290 Valine, leucine and isoleucine biosynthesis [PATH:ko00290]" ~ "BCAA biosynthesis",
      cat3 %in% "00400 Phenylalanine, tyrosine and tryptophan biosynthesis [PATH:ko00400]" ~ "AAA biosynthesis",
      cat3 %in% select_pathway[1:10] ~ "Essential amino acid metabolism",
      cat3 %in% select_pathway[11:18] ~ "VitaminB and co-factor",
      cat3 %in% select_pathway[19] ~ "Galactose metabolism",
      TRUE ~ NA_character_
    )
  )

ec_top_levels <- c(
  "BCAA biosynthesis",
  "Galactose metabolism",
  "VitaminB and co-factor",
  "AAA biosynthesis"
)

###########################
## 3. Prepare time bins
###########################
metadata_cluster <- metadata_cluster %>%
  mutate(
    bins = cut(Days, breaks = seq(100, 800, 100), include.lowest = TRUE)
  )

binnames <- names(table(metadata_cluster$bins))
metadata_cluster$bins <- factor(metadata_cluster$bins, levels = binnames)

###########################
## 4. Load QC-filtered EC table (this is your "ec_all_qc_1e6")
###########################
ec_df <- fread("./data/processed/EC_afterQC.txt")
stopifnot("sample_id" %in% colnames(ec_df))

ec_df$sample_id <- as.character(ec_df$sample_id)

###########################
## 5. Build EC long table (reused by Fig3a panel a and Fig3b)
###########################
ec_abundance_long <- melt(ec_df, id.vars = "sample_id")
colnames(ec_abundance_long) <- c("sample_id", "feature", "abundance")

ec_abundance_long <- ec_abundance_long %>%
  mutate(
    sample_id = as.character(sample_id),
    feature   = as.character(feature)
  ) %>%
  inner_join(
    metadata_cluster %>% select(sample_id, Days, bins, Cluster),
    by = "sample_id"
  )

###########################
## 6. Figure 3a - Panel a: EC_top total abundance trajectories (LOESS)
###########################
ec_top_by_feature <- KO_select %>%
  select(feature, EC_top) %>%
  distinct() %>%
  filter(!is.na(EC_top)) %>%
  filter(EC_top %in% ec_top_levels)

ec_abundance_top <- ec_abundance_long %>%
  inner_join(ec_top_by_feature, by = "feature") %>%
  group_by(sample_id, Days, bins, Cluster, EC_top) %>%
  summarise(sum_group = sum(abundance, na.rm = TRUE), .groups = "drop")

ec_abundance_top$EC_top <- factor(ec_abundance_top$EC_top, levels = ec_top_levels)

p3a_traj <- ggplot(ec_abundance_top, aes(x = Days, y = sum_group, color = Cluster)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, alpha = 0.25) +
  facet_wrap(~ EC_top, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("1" = "#4c956c", "2" = "#457b9d", "3" = "#e26d5c")) +
  labs(x = "Time (Days)", y = "Total abundance (QC-based)") +
  theme_cowplot() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  )

###########################
## 7. Figure 3a - Panel b: Species contribution within EC_top (stacked fractions)
###########################
ec_strat <- fread("./data/raw/humann2_ecs_stratified.tsv")

## Keep only samples included in metadata_cluster
sample_keep <- intersect(
  gsub("_Abundance-RPKs$", "", colnames(ec_strat)),
  as.character(metadata_cluster$sample_id)
)

## Parse stratified rows: "EC|g__Genus.s__Species"
ec_row <- ec_strat[, 1, drop = FALSE]
colnames(ec_row) <- "Gene_family"

ec_tax <- ec_row %>%
  separate(Gene_family, into = c("EC", "Taxon"), sep = "\\|", fill = "right", remove = FALSE) %>%
  separate(Taxon, into = c("Genus", "Species"), sep = "\\.", fill = "right") %>%
  mutate(
    Genus   = ifelse(is.na(Genus) | Genus == "", "unclassified", Genus),
    Species = ifelse(is.na(Species) | Species == "", "unclassified", Species)
  )

ec_mat <- ec_strat %>%
  select(1, one_of(paste0(sample_keep, "_Abundance-RPKs")))

ec_tax_ab <- bind_cols(ec_tax, ec_mat[, -1, drop = FALSE])

## Use KO_select as the EC set for Fig3a (no separate sig_EC filter)
KO_select_sig <- KO_select %>%
  filter(!is.na(EC_top)) %>%
  filter(EC_top %in% ec_top_levels)

ec_tax_ab <- ec_tax_ab %>%
  filter(EC %in% KO_select_sig$EC_number)

ec_long <- melt(ec_tax_ab, id.vars = c("Gene_family", "EC", "Genus", "Species")) %>%
  mutate(
    sample_id = gsub("_Abundance-RPKs$", "", variable),
    RA = value
  ) %>%
  select(-variable, -value) %>%
  inner_join(metadata_cluster %>% select(sample_id, bins, Cluster), by = "sample_id")

ec_mean <- ec_long %>%
  group_by(EC, Genus, Species, Cluster, bins) %>%
  summarise(meanRA = mean(RA, na.rm = TRUE), .groups = "drop")

ec_mean <- ec_mean %>%
  inner_join(
    KO_select_sig %>% select(EC_number, EC_top) %>% distinct(),
    by = c("EC" = "EC_number")
  ) %>%
  filter(!is.na(EC_top)) %>%
  mutate(EC_top = factor(EC_top, levels = ec_top_levels))

## Compute within-group fractions
ec_group <- ec_mean %>%
  group_by(Genus, Species, Cluster, bins, EC_top) %>%
  summarise(groupRA = sum(meanRA), .groups = "drop") %>%
  group_by(Cluster, bins, EC_top) %>%
  mutate(
    rank_species = ifelse(Species == "unclassified", Inf, rank(-groupRA)),
    select_species = ifelse(rank_species <= 9, Species, NA_character_),
    sum_groupRA = sum(groupRA),
    fraction = ifelse(sum_groupRA > 0, groupRA / sum_groupRA, 0)
  ) %>%
  ungroup()

## All unselected categories (including unclassified) collapsed into “Others”
ec_group <- ec_group %>%
  mutate(select_species = ifelse(is.na(select_species) | select_species == "unclassified", "Others", select_species))

ec_group$bins <- factor(ec_group$bins, levels = binnames)

## Reference line per functional group (Cluster 1, first bin, excluding Others)
ref_lines <- ec_group %>%
  filter(
    Cluster == "1",
    bins == binnames[1],
    select_species != "Others"
  ) %>%
  group_by(EC_top) %>%
  summarise(ref_y = sum(fraction, na.rm = TRUE), .groups = "drop")

## Species color palette
color.df <- fread("./data/external/species_color.txt")
species_colors <- color.df$color5
names(species_colors) <- color.df$species

## Ensure “Others” exists and is gray; keep Others last in legend
species_colors <- c(species_colors, Others = "grey70")

p3a_contrib <- ggplot(ec_group, aes(x = bins, y = fraction, fill = select_species)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  facet_wrap(EC_top ~ Cluster, ncol = 3) +
  geom_hline(
    data = ref_lines,
    aes(yintercept = ref_y),
    linetype = "dashed",
    color = "black"
  ) +
  scale_fill_manual(
    values = species_colors,
    breaks = legend_levels,
    labels = function(x) gsub("_", " ", gsub("^s__", "", x))
  ) +
  scale_x_discrete(labels = c(seq(200, 800, 100))) +
  xlab("Time (Days)") +
  ylab("Species contribution (%)") +
  guides(fill = guide_legend(reverse = FALSE)) +
  theme_cowplot()

###########################
## 8. Save Figure 3a
###########################
fig3a <- plot_grid(
  p3a_traj,
  p3a_contrib,
  ncol = 2,
  rel_widths = c(0.25, 1),
  labels = c("a", "b"),
  label_size = 14
)

pdf("./figures/main/Fig3a.pdf", width = 18, height = 11)
print(fig3a)
dev.off()

###########################
## 9. Figure 3b: EC heatmap across bins and clusters
###########################
## EC name table
ec_name <- fread("./data/processed/ec_name.txt")  # expects columns EC, EC_name

###########################
## 10. Scale abundance within each EC (z-score per feature)
###########################
ec_abundance_long2 <- ec_abundance_long %>%
  inner_join(KO_select %>% select(feature, cat3, EC_top), by = "feature") %>%
  filter(EC_top != "Essential amino acid metabolism") %>%
  group_by(feature) %>%
  mutate(scale_abund = as.numeric(scale(abundance))) %>%
  ungroup()

###########################
## 11. Summarise mean scaled abundance per EC × bin × cluster
###########################
## Map feature -> EC_number using KO table 
ec_abundance_mean <- ec_abundance_long2 %>%
  mutate(EC_number = gsub("^X", "", feature)) %>%
  group_by(feature, EC_top, EC_number, bins, Cluster, cat3) %>%
  summarise(mean = mean(scale_abund, na.rm = TRUE), .groups = "drop")

ec_abundance_mean <- ec_abundance_mean %>%
  left_join(ec_name, by = c("EC_number" = "EC"))

###########################
## 12. Load MaAsLin results (significant EC set)
###########################
load("./data/processed/allEC_MaAsLin_results.Rdata")
sig_maaslin <- all_EC.df %>% filter(qval < 0.25)
sig_features <- unique(sig_maaslin$feature)

## Keep only significant features
ec_abundance_mean <- ec_abundance_mean %>% filter(feature %in% sig_features)

###########################
## 13. Load curated EC ordering (pathway-based)
###########################
order_list_based_on_pathway <- fread("./data/external/EC-order-pathway.csv", header = FALSE)
colnames(order_list_based_on_pathway) <- c("EC_id", "EC_group", "KEGG_Pathway")
order_list_based_on_pathway <- order_list_based_on_pathway %>% filter(!duplicated(EC_id))

ec_abundance_mean2 <- ec_abundance_mean %>%
  mutate(EC_id = paste0("EC", EC_number)) %>%
  left_join(order_list_based_on_pathway, by = "EC_id") %>%
  filter(!is.na(KEGG_Pathway))

## Ordering
ec_abundance_mean2$EC_id <- factor(
  ec_abundance_mean2$EC_id,
  levels = rev(unique(order_list_based_on_pathway$EC_id))
)

###########################
## 14. Add MaAsLin significance stars (bin × contrast)
###########################
sig_select <- sig_maaslin %>%
  select(feature, value, bin, qval) %>%
  rename(Cluster = value, bins = bin)

ec_abundance_mean2 <- ec_abundance_mean2 %>%
  left_join(sig_select, by = c("feature", "Cluster", "bins")) %>%
  mutate(sig = case_when(
    qval < 0.001 ~ "***",
    qval < 0.05  ~ "**",
    qval < 0.25  ~ "*",
    TRUE ~ NA_character_
  ))

ec_abundance_mean2$bins <- factor(ec_abundance_mean2$bins, levels = binnames)

###########################
## 15. Build y-axis label
###########################
ec_abundance_mean2 <- ec_abundance_mean2 %>%
  mutate(
    EC_name2 = ifelse(is.na(EC_name) | EC_name == "", "(no name)", EC_name),
    EC_ID_name = paste0("EC ", EC_number, ": ", EC_name2)
  )

ec_abundance_mean3 <- ec_abundance_mean2 %>% arrange(EC_id)
ec_abundance_mean2$EC_ID_name <- factor(
  ec_abundance_mean2$EC_ID_name,
  levels = unique(ec_abundance_mean3$EC_ID_name)
)

###########################
## 16. Side-tile colors
###########################
col_EC_top <- c(
  "BCAA biosynthesis"        = "#3a62a0",
  "Galactose metabolism"     = "#b759a2",
  "Vitamin B and co-factor"   = "#b13a25",
  "AAA biosynthesis"         = "#e57f25"
)

col_EC_cat3 <- c(
  "00052 Galactose metabolism [PATH:ko00052]" = "#D9A3BA",
  "00290 Valine, leucine and isoleucine biosynthesis [PATH:ko00290]" = "#7A7FC2",
  "00400 Phenylalanine, tyrosine and tryptophan biosynthesis [PATH:ko00400]" = "#EFC96B",
  "00670 One carbon pool by folate [PATH:ko00670]" = "#9E7C1F",
  "00730 Thiamine metabolism [PATH:ko00730]" = "#8F4F2A",
  "00740 Riboflavin metabolism [PATH:ko00740]" = "#A85F33",
  "00750 Vitamin B6 metabolism [PATH:ko00750]" = "#B06A3B",
  "00760 Nicotinate and nicotinamide metabolism [PATH:ko00760]" = "#C47E4A",
  "00770 Pantothenate and CoA biosynthesis [PATH:ko00770]" = "#D49461",
  "00780 Biotin metabolism [PATH:ko00780]" = "#E0AC7F",
  "00790 Folate biosynthesis [PATH:ko00790]" = "#C89A2A"
)

colorall <- c(col_EC_top, col_EC_cat3)


###########################
## 17. Figure 3b: EC heatmap
###########################
manual_keep_path <- "./data/external/EC_order.txt"

ec_keep <- fread(manual_keep_path) %>% filter(KEEP == 1)

ec_keep_manual <- ec_abundance_mean2 %>%
  filter(EC_number %in% as.character(ec_keep$EC))

ec_keep_manual$EC_group <- factor(
  ec_keep_manual$EC_group,
  levels = c("BCAA biosynthesis", "Galactose metabolism", "Vitamin B and co-factor", "AAA biosynthesis")
)

pdf("./figures/main/Fig3b.pdf", width = 11, height = 10)

ggplot(ec_keep_manual, aes(x = bins, y = EC_ID_name, fill = mean)) +
  geom_tile() +
  facet_grid(EC_group ~ Cluster, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#d9ed92", mid = "#52b69a", high = "#184e77") +
  # geom_text(aes(label = sig))+
  ggnewscale::new_scale_fill() +
  geom_ysidetile(aes(x = "EC_group", yfill = EC_group)) +
  geom_ysidetile(aes(x = "EC_subgroup", yfill = KEGG_Pathway)) +
  scale_yfill_manual(values = colorall, name = "KEGG pathway") +
  scale_x_discrete(labels = c(seq(200, 800, 100))) +
  xlab("Bin (Days)") +
  theme_cowplot() +
  theme(
    legend.position = "top",
    panel.spacing.y = unit(0.01, "cm")
  )

dev.off()

