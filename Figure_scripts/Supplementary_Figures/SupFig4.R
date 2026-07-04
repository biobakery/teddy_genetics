############################################################
## Supplementary Figure 4
## Trajectory features ('measures of change') identified and
## selected for clustering of microbiome maturational patterns.
############################################################

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(traj)
library(ggplot2)
library(cowplot)

###########################
## 1. Set input/output paths
###########################
set.seed(1)
setwd("~/ddong/TEDDY_project")

genetics_file <- "./data/genetics.bed/filter/teddy_out8.qc.PCA_PC1-PC5.tsv"
metadata_file <- "./output/traj_cluster_results.Rdata"
output_dir    <- "./figures/supplement/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

###########################
## 2. Load and preprocess data (identical to SupFig3 / run_trajectory)
###########################

load(metadata_file)          # -> metadata_cluster
metadata <- metadata_cluster

geneticpc <- fread(genetics_file) %>%
  rename(subject_id = 1)

allphenotype <- metadata %>%
  inner_join(geneticpc, by = "subject_id") %>%
  filter(
    Days <= 800,
    (is.na(T1D) | T1D == FALSE),
    (is.na(AB)  | AB  == FALSE),
    !is.na(bc_dists_init)
  ) %>%
  group_by(subject_id) %>%
  filter(n() > 4) %>%
  ungroup()

allphenotype.df <- allphenotype %>%
  dcast(subject_id ~ Days, value.var = "bc_dists_init") %>%
  rename(ID = subject_id)

###########################
## 3. Trajectory feature computation and selection
###########################

step1 <- Step1Measures(
  allphenotype.df,
  Time = as.numeric(colnames(allphenotype.df)[-1]),
  ID   = TRUE
)

step2 <- Step2Selection(step1)
select_feature <- step2$selection


###########################
## 4. Attach final cluster labels
###########################
step3 <- Step3Clusters(step2,nclusters = 3)
cluster.df <- step3$partition 

select_feature.df <- cluster.df %>% 
  left_join(select_feature)

plot.df <- select_feature.df %>%
  rename(subject_id = ID) %>%
  melt(id.vars = c("subject_id", "Cluster"),
       variable.name = "Measure", value.name = "Value") %>%
  mutate(Cluster = factor(Cluster))

## z-score each measure so panels share a comparable scale
plot.df <- plot.df %>%
  group_by(Measure) %>%
  mutate(Value = as.numeric(scale(Value))) %>%
  ungroup()

###########################
## 5. Generate figure
###########################

pal <- c("#4c956c", "#457b9d", "#e26d5c")

p <- ggplot(plot.df, aes(x = Cluster, y = Value, fill = Cluster)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85) +
  facet_wrap(~ Measure, scales = "free") +
  scale_fill_manual(values = pal[seq_len(nlevels(plot.df$Cluster))]) +
  labs(x = "Maturational pattern",y="") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "grey92"))+
  theme_cowplot()

###########################
## 6. Save outputs
###########################
out_pdf <- file.path(output_dir, "SupFig4.pdf")
ggsave(out_pdf, p, width = 10, height = 6.5)

