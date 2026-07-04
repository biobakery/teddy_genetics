############################################################
## Extended Figure 5
##
## Significant microbial EC features associated with microbiome
## maturation trajectories across age bins
############################################################

###########################
## 0. Load packages
###########################
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringi)
library(stringr)
library(cowplot)
library(scales)

###########################
## 1. Load trajectory clustering results
###########################
setwd("~/ddong/TEDDY_project/")
load("./output/traj_cluster_results.Rdata")

###########################
## 2. Load MaAsLin EC association results
###########################
maaslin_root <- "./output/maaslin_EC/"
setwd(maaslin_root)
load("./allEC_result.Rdata")
ec_name <-  fread("./data/processed/ec_name.txt") 

###########################
## 3. Define age bins
###########################
bin <- cut(
  metadata_cluster$Days,
  breaks = seq(100, 800, 100),
  include.lowest = TRUE
)

metadata_cluster$bins <- bin
binnames <- names(table(bin))

###########################
## 4. Annotate EC features
###########################
table(all_EC.df$Type)

ec_name <- ec_name %>%
  mutate(feature = paste("X", EC, sep = ""))

all_EC.df <- all_EC.df %>%
  left_join(ec_name)

all_EC.df <- all_EC.df %>%
  filter(Type != "Cluster2_Cluster1")

###########################
## 5. Select significant EC features
###########################
sig <- all_EC.df %>%
  filter(qval < 0.25)

sig_summary <- sig %>%
  mutate(
    bin_group = case_when(
      bin %in% c(binnames[1:3]) ~ "Early",
      bin %in% c(binnames[4:7]) ~ "Late"
    )
  )

sig_summary <- sig_summary %>%
  group_by(EC, EC_name, bin_group) %>%
  count()

sig_summary.df <- dcast(
  sig_summary,
  "EC +EC_name ~ bin_group"
)

sig_summary_select.df <- sig_summary.df %>%
  filter(Early > 3 | !is.na(Late))

alldata_select_EC <- all_EC.df %>%
  filter(EC %in% sig$EC)

alldata_select_EC <- all_EC.df %>%
  filter(EC %in% sig_summary_select.df$EC)

###########################
## 6. Remove EC features dominated by one genus
###########################
fraction <- fread("./data/ec_genus_fraction.txt")

fraction_dominant <- fraction %>%
  filter(fraction > 50)

alldata_select_EC <- alldata_select_EC %>%
  filter(!EC %in% fraction_dominant$ec)

###########################
## 7. Add significance labels
###########################
alldata_select_EC <- alldata_select_EC %>%
  mutate(
    sig = case_when(
      qval < 0.001 ~ "***",
      qval < 0.05 ~ "**",
      qval < 0.25 ~ "*"
    )
  )

alldata_select_EC$bin <- factor(
  alldata_select_EC$bin,
  levels = binnames
)

###########################
## 8. Load curated EC feature order
###########################
ec_keep <- fread("./data/external/EC_order.txt")

ec_keep <- ec_keep %>%
  filter(KEEP == 1)

alldata_keep_EC <- all_EC.df %>%
  filter(EC %in% ec_keep$EC)

alldata_keep_EC <- alldata_keep_EC %>%
  mutate(
    sig = case_when(
      qval < 0.001 ~ "***",
      qval < 0.05 ~ "**",
      qval < 0.25 ~ "*"
    )
  )

alldata_keep_EC_all <- alldata_keep_EC

alldata_keep_EC <- alldata_keep_EC %>%
  mutate(value = coef)

alldata_keep_EC <- alldata_keep_EC %>%
  left_join(ec_keep)

###########################
## 9. Cluster EC features within functional groups
###########################
vitamin_ec.df <- alldata_keep_EC %>%
  filter(EC_top == "VitaminB and co-factor") %>%
  filter(Type == "Cluster1_Cluster3")

vitamin_ec_dcast.df <- dcast(
  vitamin_ec.df,
  "EC~bin"
)

hc.cols <- hclust(
  dist(vitamin_ec_dcast.df[, -1])
)

cluster_vt <- as.character(
  vitamin_ec_dcast.df$EC[hc.cols$order]
)

Galactose_ec.df <- alldata_keep_EC %>%
  filter(EC_top == "Galactose") %>%
  filter(Type == "Cluster1_Cluster3")

Galactose_ec_dcast.df <- dcast(
  Galactose_ec.df,
  "EC~bin"
)

hc.cols <- hclust(
  dist(Galactose_ec_dcast.df[, -1])
)

cluster_gl <- as.character(
  Galactose_ec_dcast.df$EC[hc.cols$order]
)

BACC_ec.df <- alldata_keep_EC %>%
  filter(EC_top == "BACC biosynthes") %>%
  filter(Type == "Cluster1_Cluster3")

BACC_ec_dcast.df <- dcast(
  BACC_ec.df,
  "EC~bin"
)

hc.cols <- hclust(
  dist(BACC_ec_dcast.df[, -1])
)

cluster_bacc <- as.character(
  BACC_ec_dcast.df$EC[hc.cols$order]
)

AAA_ec.df <- alldata_keep_EC %>%
  filter(EC_top == "AAA biosynthes") %>%
  filter(Type == "Cluster1_Cluster3")

AAA_ec_dcast.df <- dcast(
  AAA_ec.df,
  "EC~bin"
)

hc.cols <- hclust(
  dist(AAA_ec_dcast.df[, -1])
)

cluster_AAA <- as.character(
  AAA_ec_dcast.df$EC[hc.cols$order]
)

level <- c(
  cluster_bacc,
  cluster_AAA,
  cluster_gl,
  rev(cluster_vt)
)

###########################
## 10. Prepare curated EC heatmap data
###########################
alldata_keep_EC <- alldata_keep_EC %>%
  left_join(ec_keep)

alldata_keep_EC <- alldata_keep_EC %>%
  mutate(
    EC_ID_name = paste(EC, EC_name, sep = "_")
  )

alldata_keep_EC$EC <- factor(
  alldata_keep_EC$EC,
  levels = ec_keep$EC
)

alldata_keep_EC <- alldata_keep_EC %>%
  arrange(EC)

alldata_keep_EC$EC_ID_name <- factor(
  alldata_keep_EC$EC_ID_name,
  levels = rev(unique(alldata_keep_EC$EC_ID_name))
)

alldata_keep_EC$bin <- factor(
  alldata_keep_EC$bin,
  levels = binnames
)

alldata_keep_EC$coef[alldata_keep_EC$coef > 0.7] <- 0.7
alldata_keep_EC$coef[alldata_keep_EC$coef < (-1.7)] <- (-1.7)

###########################
## 11. Save curated EC heatmap using manual order
###########################
pdf(
  "./figures/Extended/Extended_Fig5_EC.pdf",
  width = 17,
  height = 14
)

ggplot(
  alldata_keep_EC,
  aes(x = bin, y = EC_ID_name, fill = coef)
) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#0053a8", "white", "#a80022"),
    values = rescale(c(-1.7, 0, 0.7)),
    breaks = c(-1.5, -1, -0.5, 0, 0.5),
    limits = c(-1.7, 0.7)
  ) +
  scale_x_discrete(labels = c(seq(200, 800, 100))) +
  xlab("Time(Days)") +
  geom_ysidetile(
    aes(x = "EC_group", yfill = EC_top)
  ) +
  geom_text(
    aes(label = sig),
    size = 6,
    vjust = -0.2
  ) +
  facet_wrap(~ Type) +
  theme_cowplot()

dev.off()
# 
# ###########################
# ## 12. Save curated EC heatmap using clustered order
# ###########################
# alldata_keep_EC$EC <- factor(
#   alldata_keep_EC$EC,
#   levels = level
# )
# 
# alldata_keep_EC <- alldata_keep_EC %>%
#   arrange(EC)
# 
# alldata_keep_EC$EC_ID_name <- factor(
#   alldata_keep_EC$EC_ID_name,
#   levels = rev(unique(alldata_keep_EC$EC_ID_name))
# )
# 
# pdf(
#   "./figures/Extended/Extended_Fig5_EC_clustered.pdf",
#   width = 17,
#   height = 14
# )
# 
# ggplot(
#   alldata_keep_EC,
#   aes(x = bin, y = EC_ID_name, fill = coef)
# ) +
#   geom_tile() +
#   scale_fill_gradientn(
#     colors = c("#0053a8", "white", "#a80022"),
#     values = rescale(c(-1.7, 0, 0.7)),
#     breaks = c(-1.5, -1, -0.5, 0, 0.5),
#     limits = c(-1.7, 0.7)
#   ) +
#   scale_x_discrete(labels = c(seq(200, 800, 100))) +
#   xlab("Time(Days)") +
#   geom_ysidetile(
#     aes(x = "EC_group", yfill = EC_top)
#   ) +
#   geom_text(
#     aes(label = sig),
#     size = 6,
#     vjust = -0.2
#   ) +
#   facet_wrap(~ Type) +
#   theme_cowplot()
# 
# dev.off()