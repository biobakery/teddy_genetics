############################################################
## SupFigure 7. Genetic PC2/PC3: chr6 zoom + rank scatter + top30 loadings
## Panels:
##   (a) Chr6 zoom-in: PC2- vs PC3-specific top SNPs (non-overlapping) with gene labels
##   (b) SNP rank scatter: PC2 vs PC3 (top 50 per PC; label top 10 per PC)
##   (c) Point plots: top 30 annotated SNPs for PC2 and PC3 (show loadings on both PCs)
############################################################

## Set project root
setwd("~/ddong/TEDDY_project/")

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

###########################
## 1. Read inputs (once)
###########################

## 1.1 Genetic PC loadings
loadings <- fread(
  "data/genetics.bed/filter/teddy.qc.excl_outliers.eigenvec.var",
  header = TRUE
) %>%
  select(2, 5:9) %>%           # col2: VAR (rsID); cols 5–9: PC1–PC5
  rename(snp = VAR)

## 1.2 SNP-to-gene mapping
map <- fread("output/map.txt") %>%
  rename(snp = rs) %>%
  select(snp, hgnc_symbol, hgnc_id) %>%
  distinct(snp, .keep_all = TRUE)

## 1.3 SNP positions (.map) and add two anchors for the HLA region (chr6)
allpos <- fread(
  "data/genetics.bed/filter/teddy_map.map",
  header = FALSE, sep = "\t"
)

hla_anchors <- data.frame(
  V1 = c(6, 6),
  V2 = c("HLA1", "HLA2"),
  V3 = c(0, 0),
  V4 = c(25900000, 33400000)  # HLA interval (bp)
)

allpos <- rbind(allpos, hla_anchors)

## Keep only chromosomes used in the original figure (chr3/6/8)
allpos <- allpos %>% filter(V1 %in% c(3, 6, 8))

man_pos <- allpos %>%
  rename(chr = V1, snp = V2, morgans = V3, bp = V4) %>%
  arrange(chr, bp)

## HLA interval in raw chr6 coordinates
hla_raw <- man_pos %>% filter(snp %in% c("HLA1", "HLA2"))
HLA_min_raw <- min(hla_raw$bp_raw)
HLA_max_raw <- max(hla_raw$bp_raw)

###########################
## 2. Merge loadings + annotations + positions
###########################
dat_all <- loadings %>%
  left_join(map, by = "snp") %>%
  left_join(man_pos %>% select(snp, chr, bp_raw, bp_genome), by = "snp")

###########################
## 3. Panel (a): chr6 zoom-in with PC-specific top SNPs (non-overlapping)
###########################
topN_perPC <- 20  # number of SNPs to highlight/label per PC (tunable)

dat_chr6 <- dat_all %>%
  filter(chr == 6) %>%
  mutate(
    diff_abs = abs(PC2) - abs(PC3),
    abs_diff = abs(diff_abs)
  )

## Long format for plotting
chr6_long <- dat_chr6 %>%
  select(snp, chr, bp_raw, hgnc_symbol, PC2, PC3) %>%
  pivot_longer(
    cols = c(PC2, PC3),
    names_to  = "PC",
    values_to = "loading"
  )

## Identify top SNPs per PC on chr6 among annotated SNPs only
chr6_top_PC2 <- dat_chr6 %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  slice_max(order_by = abs(PC2), n = topN_perPC) %>%
  pull(snp)

chr6_top_PC3 <- dat_chr6 %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  slice_max(order_by = abs(PC3), n = topN_perPC) %>%
  pull(snp)

## Keep PC-specific (non-overlapping) sets
chr6_top_PC2_only <- setdiff(chr6_top_PC2, chr6_top_PC3)
chr6_top_PC3_only <- setdiff(chr6_top_PC3, chr6_top_PC2)
chr6_top_union_only <- union(chr6_top_PC2_only, chr6_top_PC3_only)

## Background points: only SNPs with gene annotation (to avoid unlabeled rsIDs)
chr6_long_annot <- chr6_long %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "")

## Points to highlight (PC-specific top SNPs)
chr6_highlight <- chr6_long_annot %>%
  filter(
    (PC == "PC2" & snp %in% chr6_top_PC2_only) |
      (PC == "PC3" & snp %in% chr6_top_PC3_only)
  )

## Labels: gene symbols for highlighted SNPs
chr6_labels <- chr6_highlight %>%
  mutate(label = hgnc_symbol)

p7a <- ggplot(chr6_long_annot, aes(x = bp_raw, y = loading)) +
  annotate(
    "rect",
    xmin = HLA_min_raw, xmax = HLA_max_raw,
    ymin = -Inf, ymax = Inf,
    fill = "grey90", alpha = 0.6
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(color = "grey80", size = 0.6) +
  geom_point(
    data = chr6_highlight,
    aes(color = PC),
    size = 1.7
  ) +
  geom_label_repel(
    data = chr6_labels,
    aes(label = label, color = PC),
    size = 3,
    max.overlaps = 50,
    show.legend = FALSE
  ) +
  facet_wrap(~PC, ncol = 1, scales = "fixed") +
  scale_color_manual(values = c(PC2 = "#1f77b4", PC3 = "#d62728")) +
  labs(
    title = "Chromosome 6: PC2- vs PC3-specific top SNPs (|loading|)",
    x = "Chromosome 6 position (bp)",
    y = "Loading"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 10, face = "bold")
  )

###########################
## 4. Panel (b): PC2 vs PC3 SNP rank scatter (top 50 per PC)
###########################
rank_dat <- dat_all %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  mutate(
    rank_PC2 = rank(-abs(PC2), ties.method = "first"),
    rank_PC3 = rank(-abs(PC3), ties.method = "first"),
    label_combo = paste(snp, hgnc_symbol, sep = "_")
  )

PC2_top50 <- rank_dat %>% slice_min(rank_PC2, n = 50) %>% pull(snp)
PC3_top50 <- rank_dat %>% slice_min(rank_PC3, n = 50) %>% pull(snp)

sel_dat <- rank_dat %>%
  filter(snp %in% union(PC2_top50, PC3_top50))

PC2_top10 <- sel_dat %>% slice_min(rank_PC2, n = 10) %>% pull(snp)
PC3_top10 <- sel_dat %>% slice_min(rank_PC3, n = 10) %>% pull(snp)
label_snps <- union(PC2_top10, PC3_top10)

p7b <- ggplot(sel_dat, aes(x = log10(rank_PC2), y = log10(rank_PC3))) +
  geom_point(color = "grey70", alpha = 0.5, size = 2) +
  geom_point(
    data = sel_dat %>% filter(snp %in% PC2_top50),
    color = "#1f77b4", size = 2, alpha = 0.9
  ) +
  geom_point(
    data = sel_dat %>% filter(snp %in% PC3_top50),
    color = "#d62728", size = 2, alpha = 0.9
  ) +
  geom_text_repel(
    data = sel_dat %>% filter(snp %in% label_snps),
    aes(label = label_combo),
    size = 3,
    max.overlaps = 100
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "PC2 vs PC3 SNP rank scatter (top 50 per PC; |loading|-based)",
    x = "log10(PC2 rank) (1 = highest |loading|)",
    y = "log10(PC3 rank) (1 = highest |loading|)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 11, face = "bold")
  )

###########################
## 5. Panel (c): top 30 annotated SNPs for PC2 and PC3 (show both loadings)
###########################
pc2pc3_annot <- dat_all %>%
  select(snp, PC2, PC3, hgnc_symbol) %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  distinct(snp, .keep_all = TRUE) %>%
  mutate(label_combo = paste0(snp, "_", hgnc_symbol))

## PC2 top 30 (among annotated SNPs)
pc2_top30 <- pc2pc3_annot %>%
  slice_max(order_by = abs(PC2), n = 30) %>%
  pull(snp)

pc2_plot_dat <- pc2pc3_annot %>%
  filter(snp %in% pc2_top30) %>%
  pivot_longer(cols = c(PC2, PC3), names_to = "PC", values_to = "loading")

pc2_plot_dat$label_combo <- factor(
  pc2_plot_dat$label_combo,
  levels = pc2_plot_dat %>%
    filter(PC == "PC2") %>%
    arrange(loading) %>%
    pull(label_combo)
)

gg_pc2_top30 <- ggplot(pc2_plot_dat, aes(x = label_combo, y = loading, color = PC)) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c(PC2 = "#1f77b4", PC3 = "#d62728")) +
  labs(
    title = "PC2 top 30 loci (with gene annotation)",
    x = "SNP (rsID_gene)",
    y = "Loading",
    color = "PC"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    plot.title  = element_text(size = 11, face = "bold"),
    legend.position = "top"
  )

## PC3 top 30 (among annotated SNPs)
pc3_top30 <- pc2pc3_annot %>%
  slice_max(order_by = abs(PC3), n = 30) %>%
  pull(snp)

pc3_plot_dat <- pc2pc3_annot %>%
  filter(snp %in% pc3_top30) %>%
  pivot_longer(cols = c(PC2, PC3), names_to = "PC", values_to = "loading")

pc3_plot_dat$label_combo <- factor(
  pc3_plot_dat$label_combo,
  levels = pc3_plot_dat %>%
    filter(PC == "PC3") %>%
    arrange(loading) %>%
    pull(label_combo)
)

gg_pc3_top30 <- ggplot(pc3_plot_dat, aes(x = label_combo, y = loading, color = PC)) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c(PC2 = "#1f77b4", PC3 = "#d62728")) +
  labs(
    title = "PC3 top 30 loci (with gene annotation)",
    x = "SNP (rsID_gene)",
    y = "Loading",
    color = "PC"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    plot.title  = element_text(size = 11, face = "bold"),
    legend.position = "top"
  )

p7c <- gg_pc2_top30 + gg_pc3_top30 +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Top 30 annotated SNPs for PC2 and PC3: loadings on both PCs",
    theme = theme(plot.title = element_text(size = 11, face = "bold"))
  )

###########################
## 6. Combine panels and save as Sup_Figure10
###########################
supfig7 <- (p7a | p7b) / p7c +
  plot_layout(heights = c(1, 1))  

pdf("figures/supplement/SupFig10.pdf", width = 10, height = 14)
print(supfig7)
dev.off()

