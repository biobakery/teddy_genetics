############################################################
## Figure 4c: 
## Genome-wide SNP loadings of genetic PCs with highlighted HLA region
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
## 1. Read in PC loadings & gene annotation
###########################
## Genetic PC loadings
loadings <- fread("./data/genetics.bed/filter/teddy.qc.excl_outliers.eigenvec.var",
                  header = TRUE) %>%
  select(2, 5:9) %>%             # Column 2 = VAR (SNP); columns 5:9 = PC1–PC5
  rename(snp = VAR)

## SNP-to-gene mapping
map <- fread("./output/map.txt") %>%
  rename(snp = rs) %>%
  select(snp, hgnc_symbol, hgnc_id)

## Convert loadings to long format
ml <- loadings %>%
  melt(id = "snp", variable.name = "pc", value.name = "value")

###########################
## 2. SNP positions & chromosome offsets
###########################
## Original .map file
allpos <- fread("./data/genetics.bed/filter/teddy_map.map",
                header = FALSE, sep = "\t")

## Add two dummy anchor points for the HLA region
add <- data.frame(
  V1 = c(6, 6),
  V2 = c("HLA1", "HLA2"),
  V3 = c(0, 0),
  V4 = c(25900000, 33400000)  # HLA start/end (bp)
)

allpos <- rbind(allpos, add)

## Keep only chromosomes of interest (here: 3, 6, 8)
allpos <- allpos %>% filter(V1 %in% c(3, 6, 8))

## Rename & sort
man_pos <- allpos %>%
  rename(chr = V1, snp = V2, morgans = V3, bp = V4) %>%
  arrange(chr, bp)

## Compute per-chromosome length and cumulative offset (for genome-wide concatenation)
chr_offset <- man_pos %>%
  group_by(chr) %>%
  summarise(chr_len = max(bp), .groups = "drop") %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

## Add offsets; keep raw within-chromosome coordinate
man_pos <- man_pos %>%
  left_join(chr_offset, by = "chr") %>%
  mutate(
    bp_raw = bp,          # Raw bp (within chromosome)
    bp     = bp + offset  # Genome-wide bp (for Manhattan-style x-axis)
  )

## HLA region start/end in genome-wide coordinates
hla_range <- man_pos %>% filter(snp %in% c("HLA1", "HLA2"))
HLA_start <- min(hla_range$bp)
HLA_end   <- max(hla_range$bp)

## Merge loadings (long) with SNP positions
## Note: HLA1/HLA2 are used only for coordinate anchoring and typically have no loadings
man <- man_pos %>%
  inner_join(ml, by = "snp")

###########################
## 3. X-axis ticks (chromosome centers & boundaries)
###########################
axis.set <- man %>%
  group_by(chr) %>%
  summarise(center = (max(bp) + min(bp)) / 2,
            .groups = "keep") %>%
  rename(chr = chr)

axis.breaks <- man %>%
  group_by(chr) %>%
  summarise(min = min(bp), max = max(bp),
            .groups = "keep") %>%
  rename(chr = chr)

###########################
## 4. Annotation: top SNPs (largest squared loading) per PC
###########################
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
## 5. Manhattan-style PC loading plot with HLA block
###########################
y_min <- min(man$value, na.rm = TRUE)
y_max <- max(man$value, na.rm = TRUE)

pdf("./figures/main/Fig3c.pdf", 12, 5)

ggplot(man, aes(x = bp, y = value)) +
  facet_grid(pc ~ ., scales = "fixed") +
  geom_point(aes(colour = as.factor(chr),
                 fill   = as.factor(chr)),
             size = 0.7) +
  geom_hline(yintercept = 0,
             color = "grey",
             linetype = "dashed") +
  
  scale_colour_manual(values = c("#457b9d", "#bc4749", "#e9c46a")) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 11, face = "bold", vjust = 0.75),
    legend.text  = element_text(size = 8.75),
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
    label  = axis.set$chr,
    breaks = axis.set$center,
    guide  = guide_axis(check.overlap = TRUE),
    minor_breaks = axis.breaks$min
  ) +
  labs(x = "Chromosome",
       y = "Loading",
       colour = "Chromosome",
       fill   = "Chromosome") +
  geom_label_repel(data = man_ann,
                   aes(label = hgnc_symbol),
                   max.overlaps = 20,
                   size = 3,
                   show.legend = FALSE) +
  scale_fill_brewer(palette = "Set3")+
  ## HLA block at the bottom
  annotate("rect",
           xmin = HLA_start, xmax = HLA_end,
           ymin = y_min,
           ymax = y_min + 0.08 * (y_max - y_min),
           fill = "grey", alpha = 0.8) +
  annotate("text",
           x = (HLA_start + HLA_end) / 2,
           y = y_min + 0.04 * (y_max - y_min),
           label = "HLA region",
           size = 3) 

dev.off()
