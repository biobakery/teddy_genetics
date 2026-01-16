############################################################
## Main Figures 2a and 2b
## Microbiome trajectory and cumulative incidence curves
############################################################

## Set project root
setwd("~/ddong/TEDDY_project/")

###########################
## 0. Load packages and set seed
###########################
library(dplyr)
library(data.table)
library(nnet)
library(survival)
library(survminer)
library(patchwork)
library(cowplot)

set.seed(2025)

###########################
## 1. Load trajectory clustering results
###########################
load("~/ddong/TEDDY_project/TEDDY_project/output/traj_cluster_results.Rdata")

###########################
## 2. Estimate propensity scores using multinomial logistic regression
###########################
ps_model <- multinom(
  Cluster ~ Clinical_Center + Sex + FDR +
    delivery + PC1 + PC2 + PC3 + PC4 + PC5 + first_day +
    Antibiotics + Probiotic + time_to_brstfed_stop + age_startsolid,
  data = Cluster_results
)

## Extract propensity scores for each cluster
ps <- predict(ps_model, type = "probs")
colnames(ps) <- c("ps_0", "ps_1", "ps_2")

## Append propensity scores to the analysis dataset
data <- cbind(Cluster_results, ps)

###########################
## 3. Weighted Cox proportional hazards model
###########################
##Weighted Cox model stratified by cluster, adjusting for propensity scores
cox_model <- coxph(Surv(IA_T1D_time,AB_T1D) ~ strata(Cluster) +
                     ps_0+ps_1+ps_2, weights = weight,data = data)

###########################
## 4. Figure 2a: Smoothed microbiome trajectory plot
###########################
pdf("./figures/main/Fig2a.pdf", width = 5, height = 3.5)

ggplot(metadata_cluster,
       aes(x = Days, y = bc_dists_init, color = Cluster)) +
  geom_smooth(method = "loess") +
  theme_cowplot() +
  scale_color_manual(values = c("#4c956c", "#457b9d", "#e26d5c"))

dev.off()

###########################
## 5. Figure 2b: Cumulative incidence curves (propensity score–adjusted)
###########################
pdf("./figures/main/Fig2b.pdf", width = 4.5, height = 4.5)

ggsurvplot(
  survfit(cox_model),
  data = data,
  fun = "event",
  palette = c("#4c956c", "#457b9d", "#e26d5c"),
  conf.int = TRUE,
  xlab = "Time",
  ylab = "Cumulative Event Rate",
  xlim = c(0, 1825),
  ylim = c(0, 1),
  ggtheme = theme_cowplot()
)

dev.off()
