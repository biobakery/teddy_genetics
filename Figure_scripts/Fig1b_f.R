############################################################
## Figure 1 (panels b–f):
## b) MGX sample age distribution
## c) ImmunoChip SNP distribution (5 Mb bins)
## d) Genetic PCA (PC1 vs PC2) with marginal densities
## e) Genetic PC1 vs microbiome PCo1 
## f) PERMANOVA % variance explained (univariate)
############################################################
## Set project root
setwd("~/ddong/TEDDY_project")

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggExtra)

###########################
## 1. Shared inputs
###########################
load("./output/preprpreprocess/metadata_all.Rdata")

###########################
## 2. Panel b: MGX sample age histogram
###########################
p_b <- ggplot(metadata, aes(x = mgx_age, fill = Country)) +
  geom_histogram(bins = 40) +
  theme_cowplot() +
  ylab("Number of samples") +
  xlab("Time (Days)") +
  scale_fill_manual(values = c("#556fa1", "#c04d4d", "#e6ad5b", "#a26162"))

###########################
## 3. Panel c: ImmunoChip SNP distribution (5 Mb bins)
###########################
genetic <- fread("./data/genetics.bed/teddy-final.bim")
colnames(genetic) <- c("Chrom", "rsID", "cM", "Position", "REF", "ALT")
genetic <- genetic %>% filter(Chrom %in% 1:22)

breaks <- seq(0, 250000000, 5000000)
genetic$bin <- findInterval(genetic$Position, breaks)

genetic_summary <- genetic %>%
  group_by(Chrom, bin) %>%
  count(name = "n")

genetic_summary$Chrom <- factor(genetic_summary$Chrom, levels = rev(1:22))

maxbin <- genetic %>%
  group_by(Chrom) %>%
  summarise(maxbin = max(bin), .groups = "drop")

allbin <- data.frame(
  Chrom = rep(maxbin$Chrom, maxbin$maxbin),
  bin = unlist(lapply(maxbin$maxbin, function(x) seq_len(x)))
)

allbin$Chrom <- factor(allbin$Chrom, levels = rev(1:22))

allbin_df <- allbin %>%
  left_join(genetic_summary, by = c("Chrom", "bin")) %>%
  mutate(
    n = ifelse(is.na(n), 0, n),
    log_n = log10(n)
  )

p_c <- ggplot(allbin_df, aes(x = bin, y = Chrom, fill = log_n)) +
  geom_tile(alpha = 0.9, color = "grey") +
  theme_cowplot() +
  ylab("Chromosome") +
  xlab("Position on the chromosome (5e+06)") +
  scale_fill_viridis(option = "magma", name = "# of variants")

###########################
## 4. Panel d: Genetic PCA (PC1 vs PC2) with marginal densities
###########################
genetic_PC <- fread("./data/genetics.bed/filter/teddy_out8.qc.PCA_PC1-PC5.tsv")
genetic_PC <- genetic_PC %>% rename(subject_id = IID)

subject <- fread("./data/TEDDY_metadata/metadata_subject.txt")
all_pc <- genetic_PC %>% left_join(subject, by = "subject_id")

p_d_scatter <- ggplot(all_pc, aes(x = PC1, y = PC2, color = Country)) +
  geom_point(alpha = 1, size = 2) +
  theme_cowplot() +
  labs(x = "Genetic PC1", y = "Genetic PC2") +
  theme(legend.position = "right") +
  scale_color_manual(values = c("#1e6091", "#e63946", "#fcbf49", "#6d597a"))

p_d <- ggMarginal(
  p_d_scatter,
  type = "density",
  margins = "both",
  groupColour = TRUE,
  groupFill = TRUE
)

###########################
## 5. Panel e: Genetic PC1 vs Microbiome PCo1 (colored by days)
###########################
load("./output/microbiome_PC.Rdata")

metadata_select <- metadata %>%
  select(subject_id, sample_id, mgx_age) %>%
  mutate(sample_id = as.numeric(sample_id))

microbiomePCo <- microbiomePCo %>% inner_join(metadata_select, by = "sample_id")
all_pc_mgx <- all_pc %>% left_join(microbiomePCo, by = "subject_id")

p_e <- ggplot(all_pc_mgx, aes(x = PC1, y = Microbiome_PCo1, color = mgx_age)) +
  geom_point(size = 1) +
  theme_cowplot() +
  xlab("Genetic PC1") +
  ylab("Micorbiome PCo1") +
  scale_color_viridis(option = "magma") +
  labs(color = "Days")

###########################
## 6. Panel f: PERMANOVA (univariate % variance explained)
###########################
Permanova <- fread("./output/Permanova.txt")

uni <- Permanova %>%
  filter(type == "Uni") %>%
  arrange(R2)

uni$Covariate <- factor(uni$Covariate, levels = uni$Covariate)

p_f <- ggplot(uni, aes(x = R2 * 100, y = Covariate, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sig), vjust = 0.4) +
  theme_cowplot() +
  scale_fill_manual(values = c("#bc4749")) +
  xlab("% variance explained")

p_f_labeled <- ggdraw(p_f) +
  draw_label("f", x = 0, y = 1, hjust = 0, vjust = 1, fontface = "bold", size = 18)

###########################
## 7. Assemble multi-panel figure (b–f)
###########################
p_left <- plot_grid(
  p_b, p_d,
  ncol = 1,
  labels = c("b", "d"),
  label_fontface = "bold",
  label_size = 18,
  label_x = 0.02,
  label_y = 0.98
)

p_mid <- plot_grid(
  p_c, p_e,
  ncol = 1,
  labels = c("c", "e"),
  label_fontface = "bold",
  label_size = 18,
  label_x = 0.02,
  label_y = 0.98
)

fig1_bf <- plot_grid(
  p_left, p_mid, p_f_labeled,
  nrow = 1,
  rel_widths = c(1, 1, 1.15)
)

###########################
## 8. Save
###########################
ggsave(
  filename = "./figures/main/Fig1_b-f.pdf",
  plot = fig1_bf,
  width = 14,
  height = 6.5
)
