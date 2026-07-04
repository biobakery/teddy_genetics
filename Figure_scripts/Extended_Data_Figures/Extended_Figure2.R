############################################################
## Extended Data Figure 2: Microbiome PCoA visualization
############################################################

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

###########################
## 1. Load input data
###########################
setwd("~/ddong/TEDDY_project/")
load("./output/microbiome_PC.Rdata")
load("./output/tax_phylum_sample.Rdata")

genetic_pc <- fread("./data/genetics.bed/filter/teddy.qc.PCA_PC1-PC5.tsv") %>%
  rename(subject_id = IID)

# metadata_subject <- fread("~/ddong/TEDDY/TEDDY_metadata/metadata_subject.txt")
# 
# metadata <- fread("~/ddong/TEDDY_old/output_2024/metadata_sample.txt")
load("../output/preprocess/metadata_all.Rdata")
###########################
## 2. Merge metadata and microbiome PCoA data
###########################
mds.data$sample_id <- as.numeric(mds.data$sample_id)

merged_metadata <- metadata %>%
  inner_join(mds.data) %>%
  inner_join(genetic_pc) %>%
  inner_join(metadata_subject) %>%
  inner_join(result) %>%
  inner_join(tax_phylum.df)

###########################
## 3. Define shared PCoA plot settings
###########################
pcoa_theme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 14, colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(size = 14, face = "italic", hjust = 0.5),
    legend.position = "bottom"
  )

x_lab <- paste("Microbiome PCo1 - ", mds.var.per[1], "%", sep = "")
y_lab <- paste("Microbiome PCo2 - ", mds.var.per[2], "%", sep = "")

###########################
## 4. Plot PCoA colored by sample age
###########################
p_days <- ggplot(
  data = merged_metadata,
  aes(x = Microbiome_PCo1, y = Microbiome_PCo2)
) +
  geom_point(size = 2.5, aes(fill = Days), shape = 21) +
  scale_fill_viridis_c(option = "magma") +
  xlab(x_lab) +
  ylab(y_lab) +
  pcoa_theme

###########################
## 5. Plot PCoA colored by country
###########################
p_country <- ggplot(
  data = merged_metadata,
  aes(x = Microbiome_PCo1, y = Microbiome_PCo2)
) +
  geom_point(size = 2.5, aes(fill = Country), shape = 21) +
  scale_fill_manual(
    values = c("#1e6091", "#e63946", "#fcbf49", "#6d597a")
  ) +
  xlab(x_lab) +
  ylab(y_lab) +
  pcoa_theme

###########################
## 6. Plot PCoA colored by Bacteroidetes abundance
###########################
p_bacteroidetes <- ggplot(
  data = merged_metadata,
  aes(x = Microbiome_PCo1, y = Microbiome_PCo2)
) +
  geom_point(size = 2.5, aes(fill = Bacteroidetes), shape = 21, alpha = 0.8) +
  xlab(x_lab) +
  ylab(y_lab) +
  pcoa_theme

###########################
## 7. Plot PCoA colored by Firmicutes abundance
###########################
p_firmicutes <- ggplot(
  data = merged_metadata,
  aes(x = Microbiome_PCo1, y = Microbiome_PCo2)
) +
  geom_point(size = 2.5, aes(fill = Firmicutes), shape = 21, alpha = 0.8) +
  xlab(x_lab) +
  ylab(y_lab) +
  pcoa_theme

###########################
## 8. Assemble multi-panel figure
###########################

panel_a <- plot_grid(
  p_firmicutes,
  p_bacteroidetes,     
  labels = c("a", ""),
  nrow = 1,
  label_fontface = "bold",
  label_size = 16,
  label_x = 0.02,
  label_y = 0.98
)

panel_b <- plot_grid(
  p_days,
  p_country,
  labels = c("b", ""),
  nrow = 1,
  label_fontface = "bold",
  label_size = 16,
  label_x = 0.02,
  label_y = 0.98
)

extended_fig2 <- plot_grid(
  panel_a,
  panel_b,
  ncol = 1
)

###########################
## 9. Save figure
###########################

ggsave(
  filename = "./figures/Extended/Extended_Fig2_microbiome_PCoA.pdf",
  plot = extended_fig2,
  width = 10,
  height = 10
)