library(dplyr)
library(ggplot2)
library(cowplot)
library(survival)
library(survminer)

setwd("~/ddong/TEDDY_project/")
##---------------------------------------------------------------
## Fig. 4b: Cumulative event rate plot (Stratified by PC3 level)
##---------------------------------------------------------------
load("./output/Fig2a_b.Rdata")
# Calculate the quantiles of the PCs and divide the data into two parts (50% quantile)
quantiles <- quantile(data$PC1, probs = c(0, 1/2,1))
data$PC1_half <- cut(data$PC1, breaks = quantiles, labels = c("Q1", "Q2" ),include.lowest = T)
quantiles <- quantile(data$PC2, probs = c(0, 1/2,1))
data$PC2_half <- cut(data$PC2, breaks = quantiles, labels = c("Q1", "Q2" ),include.lowest = T)
quantiles <- quantile(data$PC3, probs = c(0, 1/2,1))
data$PC3_half <- cut(data$PC3, breaks = quantiles, labels = c("Q1", "Q2" ),include.lowest = T)
quantiles <- quantile(data$PC4, probs = c(0, 1/2,1))
data$PC4_half <- cut(data$PC4, breaks = quantiles, labels = c("Q1", "Q2" ),include.lowest = T)
quantiles <- quantile(data$PC5, probs = c(0, 1/2,1))
data$PC5_half <- cut(data$PC5, breaks = quantiles, labels = c("Q1", "Q2" ),include.lowest = T)

Cluster800PC3Q1 <- data %>% filter(PC3_half == "Q1")
Cluster800PC3Q2 <- data %>% filter(PC3_half == "Q2")
cox_modelQ1 <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ strata(Cluster) +
                       ps_0+ps_1+ps_2, data = Cluster800PC3Q1)
cox_modelQ2 <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~strata(Cluster) +
                       ps_0+ps_1+ps_2, data = Cluster800PC3Q2)

pdf("./figures/Fig4b",width = 4.5,height = 4.5)
ggsurvplot(survfit(cox_modelQ1), 
           data = Cluster800PC3Q1, 
           fun = "event", 
           palette = c("#4c956c","#457b9d","#e26d5c"),
           # conf.int.alpha = 0.1,
           conf.int = T, # 添加置信区间
           xlab = "Time", 
           ylab = "Cumulative Event Rate",
           ggtheme = theme_cowplot()) # 使用简洁的主题

ggsurvplot(survfit(cox_modelQ2), 
           data = Cluster800PC3Q2, 
           fun = "event", 
           palette = c("#4c956c","#457b9d","#e26d5c"),
           # conf.int.alpha = 0.1,
           conf.int = T, # 添加置信区间
           xlab = "Time", 
           ylab = "Cumulative Event Rate",
           ggtheme = theme_cowplot()) # 使用简洁的主题

dev.off()



