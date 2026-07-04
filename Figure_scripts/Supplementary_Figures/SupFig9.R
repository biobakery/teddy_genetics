############################################################
## Supplementary Figure 9
## Genome-wide distribution of SNP loadings for principal components
## (Manhattan-style visualization across autosomes)
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

###########################
## 1. Read PC loadings and SNP annotations
###########################
## Genetic PC loadings (PC1–PC5)
loadings <- fread(
  "./data/genetics.bed/filter/teddy.qc.excl_outliers.eigenvec.var",
  header = TRUE
) %>%
  select(2, 5:9) %>%              # VAR = rsID; PC1–PC5
  rename(snp = VAR)

## SNP-to-gene mapping
map <- fread("./output/map.txt") %>%
  rename(snp = rs) %>%
  select(snp, hgnc_symbol, hgnc_id)

## Convert loadings to long format
ml <- loadings %>%
  melt(id = "snp", variable.name = "pc", value.name = "value")

###########################
## 2. SNP positions and genome-wide coordinates
###########################
## Read PLINK .map file
allpos <- fread(
  "./data/genetics.bed/filter/teddy_map.map",
  header = FALSE,
  sep = "\t"
)

## Add dummy anchor SNPs to indicate the HLA region on chr6
add <- data.frame(
  V1 = c(6, 6),
  V2 = c("HLA1", "HLA2"),
  V3 = c(0, 0),
  V4 = c(25900000, 33400000)
)

allpos <- rbind(allpos, add)

## Rename columns and sort by chromosome and position
man_pos <- allpos %>%
  rename(chr = V1, snp = V2, morgans = V3, bp = V4) %>%
  arrange(chr, bp)

## Compute chromosome offsets for a concatenated (Manhattan-style) x-axis
chr_offset <- man_pos %>%
  group_by(chr) %>%
  summarise(chr_len = max(bp), .groups = "drop") %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

man_pos <- man_pos %>%
  left_join(chr_offset, by = "chr") %>%
  mutate(
    bp_raw = bp,          # within-chromosome position
    bp     = bp + offset  # genome-wide position
  )

## Define the HLA interval in genome-wide coordinates
hla_range <- man_pos %>% filter(snp %in% c("HLA1", "HLA2"))
HLA_start <- min(hla_range$bp)
HLA_end   <- max(hla_range$bp)

## Merge SNP positions with PC loadings
man <- man_pos %>%
  inner_join(ml, by = "snp")

###########################
## 3. X-axis settings (chromosome centers and boundaries)
###########################
## Chromosome centers for x-axis labels
axis.set <- man %>%
  group_by(chr) %>%
  summarise(center = (max(bp) + min(bp)) / 2, .groups = "keep")

## Chromosome boundaries for minor breaks
axis.breaks <- man %>%
  group_by(chr) %>%
  summarise(min = min(bp), max = max(bp), .groups = "keep")

###########################
## 4. Select top SNPs per PC for annotation
###########################
## For each PC, select the top 10 SNPs based on squared loadings
man_ann <- loadings %>%
  melt(id = "snp", variable.name = "pc", value.name = "value") %>%
  mutate(v2 = value^2) %>%
  group_by(pc) %>%
  top_n(v2, n = 10) %>%
  ungroup() %>%
  inner_join(map, by = "snp") %>%
  select(pc, snp, hgnc_symbol) %>%
  inner_join(man, by = c("pc", "snp")) %>%
  mutate(id = paste(pc, hgnc_symbol)) %>%
  filter(!duplicated(id)) %>%
  as.data.frame()

###########################
## 5. Manhattan-style plot of PC loadings
###########################
y_min <- min(man$value, na.rm = TRUE)
y_max <- max(man$value, na.rm = TRUE)

## Color chromosomes by parity
man <- man %>%
  mutate(
    color = case_when(
      chr %% 2 == 1 ~ "Odd",
      chr %% 2 == 0 ~ "Even"
    )
  ) %>%
  filter(chr <= 22)

pdf("./figures/supplement/SupFig9.pdf", width = 15, height = 7)

ggplot(man, aes(x = bp, y = value)) +
  facet_grid(pc ~ ., scales = "fixed") +
  geom_point(aes(colour = color, fill = color), size = 0.7) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  scale_colour_manual(values = c("#274c77", "#a3cef1")) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.minor.x = element_line(colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.title  = element_text(size = 12.5, face = "bold"),
    strip.text  = element_text(size = 12.5, face = "bold")
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = axis.set$center,
    labels = axis.set$chr,
    minor_breaks = axis.breaks$min,
    guide = guide_axis(check.overlap = TRUE)
  ) +
  labs(
    x = "Chromosome",
    y = "Loading",
    colour = "Chromosome",
    fill   = "Chromosome"
  ) +
  geom_label_repel(
    data = man_ann,
    aes(label = hgnc_symbol),
    size = 3,
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  scale_fill_brewer(palette = "Set3")+
  guides(color = "none",fill = "none")

dev.off()
