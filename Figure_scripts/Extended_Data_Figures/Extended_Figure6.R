############################################################
## Extended Figure 6
##
## Significant microbial metabolic pathways associated with
## microbiome maturation trajectories across age bins
############################################################

###########################
## 0. Load libraries
###########################
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ggside)
library(scales)

###########################
## 1. Load pathway association results
###########################
setwd("~/ddong/TEDDY_project/")
pathway <- fread("./output/maaslin_pathway005_sig.txt")
pathway_name <- fread("./data/pathway_map.txt")

###########################
## 2. Select curated significant pathways
###########################
sig005 <- pathway %>%
  filter(select == 1)

pathway_group <- sig005 %>%
  select(feature, Group) %>%
  unique()

alldata_select_pathway <- all_EC.df %>%
  filter(feature %in% sig005$feature)

alldata_select_pathway <- alldata_select_pathway %>%
  left_join(pathway_group) %>%
  left_join(pathway_name)

alldata_select_pathway_all <- alldata_select_pathway

###########################
## 3. Add significance labels
###########################
alldata_select_pathway <- alldata_select_pathway %>%
  mutate(
    sig = case_when(
      qval < 0.001 ~ "***",
      qval < 0.05  ~ "**",
      qval < 0.25  ~ "*"
    )
  )

alldata_select_pathway$bin <- factor(
  alldata_select_pathway$bin,
  levels = binnames
)

###########################
## 4. Cluster pathways for visualization
###########################
alldata_keep_pathway <- alldata_select_pathway %>%
  mutate(value = coef)

pathway_cluster.df <- alldata_keep_pathway %>%
  filter(Type == "Cluster1_Cluster3")

pathway_cluster_dcast.df <- dcast(
  pathway_cluster.df,
  "feature ~ bin"
)

pathway_cluster_dcast.df <- pathway_cluster_dcast.df[
  complete.cases(pathway_cluster_dcast.df),
]

hc.cols <- hclust(
  dist(pathway_cluster_dcast.df[, -1])
)

pathway_order <- as.character(
  pathway_cluster_dcast.df$feature[hc.cols$order]
)

###########################
## 5. Prepare pathway heatmap data
###########################
alldata_select_pathway$feature <- factor(
  alldata_select_pathway$feature,
  levels = pathway_order
)

alldata_select_pathway <- alldata_select_pathway %>%
  arrange(Group, feature)

alldata_select_pathway$pathway <- factor(
  alldata_select_pathway$pathway,
  levels = unique(alldata_select_pathway$pathway)
)

alldata_select_pathway <- alldata_select_pathway %>%
  filter(!is.na(feature))

alldata_select_pathway2 <- alldata_select_pathway

alldata_select_pathway2$coef[
  alldata_select_pathway$coef >= 1.1
] <- 1.1

alldata_select_pathway2$coef[
  alldata_select_pathway$coef <= -2.1
] <- -2.1

###########################
## 6. Generate pathway heatmap
###########################
p_pathway <- ggplot(
  alldata_select_pathway2,
  aes(x = bin, y = pathway, fill = coef)
) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#0053a8", "white", "#a80022"),
    values = rescale(c(-2.1, 0, 1.1)),
    breaks = c(-2, -1, 0, 1),
    limits = c(-2.1, 1.1)
  ) +
  scale_x_discrete(
    labels = c(seq(200, 800, 100))
  ) +
  xlab("Time(Days)") +
  geom_text(
    aes(label = sig),
    size = 6
  ) +
  facet_wrap(~ Type) +
  geom_ysidetile(
    aes(x = "Pathway group", yfill = Group)
  ) +
  theme_cowplot()

p_pathway

###########################
## 7. Save Extended Figure 6
###########################
ggsave(
  filename = "./figures/Extended/Extended_Fig6_pathway.pdf",
  plot = p_pathway,
  width = 16,
  height = 7
)