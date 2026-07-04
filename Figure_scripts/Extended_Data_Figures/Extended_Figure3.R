############################################################
## Extended Figure 3
##
## a) Shannon diversity across microbiome maturation clusters
##    before 400 days of age
## b) Shannon diversity across microbiome maturation clusters
##    after 400 days of age
############################################################

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggpubr)
library(cowplot)

###########################
## 1. Set project root and load input data
###########################
setwd("~/ddong/TEDDY_project/")

microbiome <- fread(
  "~/ddong/TEDDY_project/output/metaphlan2_afterQC.tsv",
  header = TRUE
)

load("./output/traj_cluster_results.Rdata")

###########################
## 2. Calculate Shannon diversity
###########################
shannon <- diversity(
  t(microbiome[, -1]),
  index = "shannon"
)

shannon.df <- data.frame(
  sample_id = colnames(microbiome)[-1],
  shannon = shannon
)

###########################
## 3. Merge Shannon diversity with cluster metadata
###########################
metadata_cluster$sample_id <- as.character(metadata_cluster$sample_id)

shannon.df <- shannon.df %>%
  inner_join(metadata_cluster)

shannon.df$Cluster <- factor(
  shannon.df$Cluster,
  levels = c(1, 2, 3),
  labels = c(
    "Early Matured",
    "Late Matured",
    "Early Plateaued"
  )
)

###########################
## 4. Split samples by age
###########################
shannon.df_before400 <- shannon.df %>%
  filter(Days <= 400)

shannon.df_after400 <- shannon.df %>%
  filter(Days > 400)

###########################
## 5. Define cluster comparisons and colors
###########################
cluster_comparisons <- list(
  c("Early Matured", "Late Matured"),
  c("Early Matured", "Early Plateaued"),
  c("Late Matured", "Early Plateaued")
)

cluster_colors <- c(
  "Early Matured" = "#4c956c",
  "Late Matured" = "#457b9d",
  "Early Plateaued" = "#e26d5c"
)

###########################
## 6. Generate Shannon diversity boxplots
###########################
plot <- list()

plot[[1]] <- ggplot(
  shannon.df_before400,
  aes(x = Cluster, y = shannon, color = Cluster)
) +
  geom_boxplot() +
  stat_compare_means(
    method = "t.test",
    comparisons = cluster_comparisons
  ) +
  scale_color_manual(values = cluster_colors) +
  ggtitle("Before 400 Days") +
  ylab("Shannon Diversity") +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      vjust = 1
    )
  )

plot[[2]] <- ggplot(
  shannon.df_after400,
  aes(x = Cluster, y = shannon, color = Cluster)
) +
  geom_boxplot() +
  stat_compare_means(
    method = "t.test",
    comparisons = cluster_comparisons
  ) +
  scale_color_manual(values = cluster_colors) +
  ggtitle("After 400 Days") +
  ylab("Shannon Diversity") +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      vjust = 1
    )
  )

###########################
## 7. Assemble Extended Figure 3
###########################
extended_fig3 <- plot_grid(
  plotlist = plot,
  ncol = 2,
  labels = c("a", "b"),
  label_fontface = "bold",
  label_size = 18
)

extended_fig3

###########################
## 8. Save Extended Figure 3
###########################
ggsave(
  filename = "./figures/Extended/Extended_Fig3_shannon_diversity.pdf",
  plot = extended_fig3,
  width = 8,
  height = 4
)