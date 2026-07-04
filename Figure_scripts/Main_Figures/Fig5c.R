############################################################
## Fig 5c:
## Distribution of microbiome-associated signals across ImmunoChip loci
############################################################

## Set project root
setwd("~/ddong/TEDDY_project/")

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)

###########################
## 1. Read results and parse rsid
###########################
all <- fread("./output/SNP_linear_Species_toddler.txt")

all <- all %>%
  mutate(
    rsid = case_when(
      grepl(":", Predictor) ~ str_extract(Predictor, ".*(?=\\:)"),
      TRUE ~ Predictor
    )
  )

###########################
## 2. Keep SNP:time rows and compute -log10(P)
###########################
fix <- all %>%
  filter(grepl(":", Predictor)) %>%
  rename(snp = rsid) %>%
  mutate(log10P = -log10(P))

###########################
## 3. Keep only species with at least one genome-wide significant hit
###########################
sig_hits <- fix %>% filter(P < 5e-9)
keep_species <- unique(sig_hits$Species)

fix_species <- fix %>% filter(Species %in% keep_species)

###########################
## 4. SNP positions (chr, bp) and cumulative genomic coordinate
###########################
pos <- fread("./data/genetics.bed/filter/teddy_map.map", header = FALSE, sep = "\t") %>%
  rename(chr = V1, snp = V2, morgans = V3, bp = V4) %>%
  mutate(
    chr = as.integer(chr),
    bp  = as.numeric(bp)
  ) %>%
  filter(!is.na(chr), !is.na(bp)) %>%
  arrange(chr, bp)

pos <- pos %>% filter(chr != 25)

chr_info <- pos %>%
  group_by(chr) %>%
  summarise(chr_len = max(bp, na.rm = TRUE), .groups = "drop") %>%
  arrange(chr) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

pos_cum <- pos %>%
  left_join(chr_info, by = "chr") %>%
  mutate(bp_cum = bp + offset)

axis_set <- pos_cum %>%
  group_by(chr) %>%
  summarise(center = (min(bp_cum) + max(bp_cum)) / 2, .groups = "drop")

axis_breaks <- pos_cum %>%
  group_by(chr) %>%
  summarise(min = min(bp_cum), max = max(bp_cum), .groups = "drop")

###########################
## 5. Merge positions with association results
###########################
man <- pos_cum %>%
  inner_join(fix_species, by = "snp")

###########################
## 6. Taxonomy join
###########################
tax <- fread("./output/species_major_tree.txt")
tax_select <- tax %>%
  filter(species %in% man$Species) %>%
  select(phylum, genus, species)

man <- man %>% left_join(tax_select, by = c("Species" = "species"))

###########################
## 7. Color logic
###########################
man <- man %>%
  mutate(
    Color = case_when(
      chr %in% c(seq(1, 22, 2)) ~ "a",
      chr %in% c(seq(2, 22, 2)) ~ "b",
      TRUE ~ "a"
    ),
    Color_combine = case_when(
      genus == "g__Bacteroides"      & P < 5e-9 ~ "c",
      genus == "g__Barnesiella"      & P < 5e-9 ~ "d",
      genus == "g__Bifidobacterium"  & P < 5e-9 ~ "e",
      genus == "g__Collinsella"      & P < 5e-9 ~ "f",
      genus == "g__Coprococcus"      & P < 5e-9 ~ "g",
      genus == "g__Dorea"            & P < 5e-9 ~ "h",
      genus == "g__Eubacterium"      & P < 5e-9 ~ "i",
      genus == "g__Faecalibacterium" & P < 5e-9 ~ "j",
      genus == "g__Parabacteroides"  & P < 5e-9 ~ "k",
      genus == "g__Prevotella"       & P < 5e-9 ~ "l",
      genus == "g__Ruminococcus"     & P < 5e-9 ~ "m",
      genus == "g__Streptococcus"    & P < 5e-9 ~ "n",
      genus == "g__Subdoligranulum"  & P < 5e-9 ~ "o",
      TRUE ~ NA_character_
    ),
    Color_all = if_else(is.na(Color_combine), Color, Color_combine)
  )

###########################
## 8. Plot: Manhattan-style on cumulative coordinate
###########################
colors <- rainbow(13)
colors <- c(colors, "#00C4FF")

manhplot <- ggplot(man, aes(x = bp_cum, y = log10P)) +
  geom_point(aes(colour = Color_all), alpha = 0.8, size = 1.2) +
  geom_hline(yintercept = -log10(5e-9), color = "grey20", linetype = "dashed") +
  scale_color_manual(values = c("#ccc5b9", "#8d99ae", colors)) +
  theme_bw() +
  ggtitle("Toddler_Fixed") +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12.5, face = "bold")
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = axis_set$center,
    labels = axis_set$chr,
    minor_breaks = axis_breaks$min,
    guide = guide_axis(check.overlap = TRUE)
  ) +
  labs(x = "Chromosome", y = "-log10(P)") +
  geom_vline(
    xintercept = axis_breaks$min,
    colour = "black",
    linewidth = 0.3
  )

ggsave(
  plot = manhplot,
  filename = "./figures/main/fig5c.png",
  dpi = 300,
  height = 7,
  width = 13
)

###########################
## 9. Legend-only export
###########################
legend_gtable <- cowplot::get_legend(manhplot)

ggsave(
  filename = "./figures/main/fig5c_legend_only.pdf",
  plot = legend_gtable,
  width = 5,
  height = 1
)
