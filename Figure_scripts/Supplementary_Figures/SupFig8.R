############################################################
## Supplementary Figure 8
## Time-stratified, pathway-level Shannon diversity of EC
## composition across microbiome maturation clusters.
############################################################

###########################
## 0. Load packages
###########################

library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(vegan)          
library(ggplot2)
library(ggpubr)         
library(cowplot)

set.seed(1)

###########################
## 1. Set input/output paths
###########################

setwd("~/ddong/TEDDY_project")

output_dir <- "./figures/supplement/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

###########################
## 2. Load metadata + define bins
###########################

load("./output/traj_cluster_results.Rdata")   # -> metadata_cluster

metadata_cluster <- metadata_cluster %>%
  mutate(
    sample_id = as.character(sample_id),
    Cluster   = factor(Cluster,
                       levels = c("Early Matured", "Late Matured", "Early Plateaued")),
    bins      = cut(Days, breaks = seq(0, 800, by = 100),
                    right = FALSE, include.lowest = TRUE),
    bin_start = as.numeric(gsub("\\[|\\)|\\,", "",
                                sub(",.*", "", as.character(bins)))),
    bin_mid   = bin_start + 50
  )

sample_bin <- metadata_cluster %>%
  select(sample_id, bins, bin_start, bin_mid, Cluster3) %>%
  filter(!is.na(bins), !is.na(Cluster3)) %>%
  distinct()

###########################
## 3. Define EC functional categories (EC_top), as in Fig3a
###########################

KO <- fread("./data/external/KO_03202022.csv") %>%
  mutate(EC_number = gsub("^EC", "", ec))

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
      cat3 %in% select_pathway[11:18] ~ "VitaminB and co-factor",
      cat3 %in% select_pathway[19]    ~ "Galactose metabolism",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(EC_top))

ec_top_levels <- c("BCAA biosynthesis",
                   "Galactose metabolism",
                   "VitaminB and co-factor",
                   "AAA biosynthesis")

###########################
## 4. Load species-stratified EC abundance table
##    (rows: "EC|g__Genus.s__Species", cols: <sample>_Abundance-RPKs)
###########################

ec_strat <- fread("./data/raw/humann2_ecs_stratified.tsv")
colnames(ec_strat)[1] <- "Gene_family"

ec_renorm <- ec_strat %>%
  separate(Gene_family, into = c("EC", "Taxon"),
           sep = "\\|", fill = "right", remove = TRUE) %>%
  separate(Taxon, into = c("Genus", "Species"),
           sep = "\\.", fill = "right") %>%
  mutate(Species = ifelse(is.na(Species) | Species == "", "unclassified", Species))

## Clean sample column names to match sample_id
sample_cols <- setdiff(colnames(ec_renorm), c("EC", "Genus", "Species"))
setnames(ec_renorm, sample_cols, gsub("_Abundance-RPKs$", "", sample_cols))
sample_cols <- gsub("_Abundance-RPKs$", "", sample_cols)

###########################
## 5. Per-category Shannon diversity of EC composition
##    (sum abundance by Species within a category, then diversity())
###########################

alpha_list <- list()

for (grp in ec_top_levels) {
  ecs <- KO_select %>% filter(EC_top == grp) %>% pull(EC_number) %>% unique()
  
  sub <- ec_renorm %>%
    filter(EC %in% ecs) %>%
    select(Species, all_of(sample_cols)) %>%
    group_by(Species) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  
  if (nrow(sub) < 2) next
  
  mat <- t(as.matrix(sub[, -1]))          # rows = samples, cols = species
  shannon_vec <- diversity(mat, index = "shannon")
  
  alpha_list[[grp]] <- data.frame(
    sample_id = sample_cols,
    Group     = grp,
    shannon   = as.numeric(shannon_vec),
    stringsAsFactors = FALSE
  )
}

alpha_all <- rbind(alpha_list$`BCAA biosynthesis`,alpha_list$`Galactose metabolism`,
                   alpha_list$`VitaminB and co-factor`,alpha_list$`AAA biosynthesis`)
alpha_all$Group <- factor(alpha_all$Group, levels = ec_top_levels)
alpha.df <- alpha_all%>%
  inner_join(sample_bin, by = "sample_id")  

alpha.df <- alpha.df %>%
  mutate(
    Cluster = factor(
      Cluster,
      levels = c(1, 2, 3),
      labels = c(
        "Early Matured",
        "Late Matured",
        "Early Plateaued"
      )
    )
  )

###########################
## 6. Panel B: per-bin boxplots across clusters (pairwise Wilcoxon)
###########################

pairs_list <- list(
  c("Early Matured", "Late Matured"),
  c("Early Matured", "Early Plateaued"),
  c("Late Matured",  "Early Plateaued")
)

cols_cluster <- c("Early Matured"   = "#4c956c",
                  "Late Matured"    = "#457b9d",
                  "Early Plateaued" = "#e26d5c")

p_box <- ggplot(alpha.df, aes(x = Cluster, y = shannon, color = Cluster)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
  stat_compare_means(comparisons = pairs_list, method = "wilcox.test",
                     size = 2.6, hide.ns = TRUE) +
  facet_grid(Group ~ bins, scales = "free_y") +
  scale_color_manual(values = cols_cluster) +
  labs(x = "Maturation pattern", y = "Shannon diversity of EC composition") +
  theme_cowplot() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

###########################
## 7. Assemble and save
###########################

pdf("./figures/supplement/SupFig8.pdf", width = 14, height = 12)
p_box
dev.off()
